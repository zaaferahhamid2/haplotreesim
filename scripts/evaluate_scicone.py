"""
Evaluate SCICoNE output against HaploTreeSim ground truth.
Uses SCICoNE's own cell_node_labels, segmented_regions, and tree edges.
"""

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from dataset_io import load_dataset


def parse_scicone_tree_str(path: Path):
    """Parse SCICoNE's tree_str.txt into parent links and per-node CN deltas."""
    parents = {}
    events = {}
    pattern = re.compile(r"node\s+(\d+):\s+p_id:([^,]+),\[(.*)\]$")

    with open(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            match = pattern.match(line)
            if match is None:
                continue

            node_id = int(match.group(1))
            parent_text = match.group(2)
            parent_id = None if parent_text == "NULL" else int(parent_text)
            node_events = {}
            events_text = match.group(3).strip()
            if events_text:
                for item in events_text.split(","):
                    region, delta = item.split(":")
                    node_events[int(region)] = int(delta)

            parents[node_id] = parent_id
            events[node_id] = node_events

    if not parents:
        raise ValueError(f"No SCICoNE node lines found in {path}")

    edges = [(parent, child) for child, parent in parents.items() if parent is not None]
    return parents, events, edges


def cumulative_scicone_profiles(parents, events, n_regions: int):
    """Reconstruct each SCICoNE tree node's cumulative total-CN profile by region."""
    profiles = {}

    def profile_for(node_id):
        if node_id in profiles:
            return profiles[node_id]
        parent_id = parents[node_id]
        if parent_id is None:
            profile = np.full(n_regions, 2, dtype=int)
        else:
            profile = profile_for(parent_id).copy()
        for region, delta in events[node_id].items():
            if 0 <= region < n_regions:
                profile[region] += delta
        profiles[node_id] = profile
        return profile

    for node_id in parents:
        profile_for(node_id)
    return profiles


def bins_to_scicone_regions(bin_profile: np.ndarray, region_sizes: np.ndarray) -> np.ndarray:
    """Aggregate a bin-level SCICoNE CN profile to SCICoNE segmented regions."""
    if len(bin_profile) == len(region_sizes):
        return bin_profile.astype(int)
    if int(region_sizes.sum()) != len(bin_profile):
        raise ValueError(
            "Cannot align inferred_cnvs.csv with segmented_region_sizes.txt: "
            f"{len(bin_profile)} bins vs {int(region_sizes.sum())} region-sized bins"
        )

    region_profile = []
    start = 0
    for size in region_sizes:
        end = start + int(size)
        region_profile.append(int(round(float(np.mean(bin_profile[start:end])))))
        start = end
    return np.array(region_profile, dtype=int)


def modal_label_profiles(node_labels, tcn_pred_bins: np.ndarray, region_sizes: np.ndarray):
    """Return one representative region-level CN profile for each SCICoNE label."""
    labels = sorted(set(node_labels))
    profiles = {}
    labels_array = np.array(node_labels)

    for label in labels:
        label_rows = tcn_pred_bins[labels_array == label]
        region_rows = np.vstack([
            bins_to_scicone_regions(row, region_sizes) for row in label_rows
        ])
        unique_rows, counts = np.unique(region_rows, axis=0, return_counts=True)
        profiles[label] = unique_rows[int(np.argmax(counts))]

    return profiles


def match_scicone_labels_to_tree_nodes(label_profiles, tree_profiles):
    """Map SCICoNE labels such as A/B/C to numeric tree node IDs by CN profile."""
    labels = sorted(label_profiles)
    tree_nodes = sorted(tree_profiles)
    costs = np.array([
        [np.sum(np.abs(label_profiles[label] - tree_profiles[node])) for node in tree_nodes]
        for label in labels
    ])

    label_to_tree_node = {}
    if len(tree_nodes) >= len(labels):
        row_ind, col_ind = linear_sum_assignment(costs)
        for row, col in zip(row_ind, col_ind):
            label_to_tree_node[labels[row]] = tree_nodes[col]
    else:
        for row, label in enumerate(labels):
            label_to_tree_node[label] = tree_nodes[int(np.argmin(costs[row]))]

    total_l1 = int(sum(
        np.sum(np.abs(label_profiles[label] - tree_profiles[node]))
        for label, node in label_to_tree_node.items()
    ))
    return label_to_tree_node, total_l1


def positive_hungarian_mapping(true_labels, pred_labels):
    """One-to-one pred-label to true-clone mapping, excluding zero-overlap matches."""
    true_unique = sorted(np.unique(true_labels).tolist())
    pred_unique = sorted(set(pred_labels))
    contingency = np.zeros((len(true_unique), len(pred_unique)), dtype=int)
    true_to_idx = {label: i for i, label in enumerate(true_unique)}
    pred_to_idx = {label: i for i, label in enumerate(pred_unique)}

    for true_label, pred_label in zip(true_labels, pred_labels):
        contingency[true_to_idx[true_label], pred_to_idx[pred_label]] += 1

    row_ind, col_ind = linear_sum_assignment(-contingency)
    mapping = {}
    total_overlap = 0
    for row, col in zip(row_ind, col_ind):
        overlap = int(contingency[row, col])
        if overlap <= 0:
            continue
        mapping[pred_unique[col]] = int(true_unique[row])
        total_overlap += overlap

    accuracy = total_overlap / len(true_labels) if len(true_labels) else 0.0
    return mapping, accuracy


def attached_leaf_bipartitions(edges, attachments):
    """
    Build rooted bipartitions after attaching one comparable leaf per clone.

    This lets us compare cell-bearing internal nodes, including the normal/root
    clone, as labeled leaves without pretending SCICoNE and HaploTreeSim share
    raw node IDs.
    """
    children = {}
    all_nodes = set()
    child_nodes = set()
    for parent, child in edges:
        children.setdefault(parent, []).append(child)
        children.setdefault(child, [])
        all_nodes.update([parent, child])
        child_nodes.add(child)

    for clone_id, node_id in attachments.items():
        leaf_id = ("clone", int(clone_id))
        children.setdefault(node_id, []).append(leaf_id)
        children[leaf_id] = []
        all_nodes.update([node_id, leaf_id])
        child_nodes.add(leaf_id)

    roots = sorted([node for node in all_nodes if node not in child_nodes], key=str)
    if not roots:
        return set()
    leaf_labels = {int(clone_id) for clone_id in attachments}

    def descendants(node_id):
        if isinstance(node_id, tuple) and node_id[0] == "clone":
            return {int(node_id[1])}
        desc = set()
        for child in children.get(node_id, []):
            desc.update(descendants(child))
        return desc

    bipartitions = set()
    for node_id in all_nodes:
        if isinstance(node_id, tuple) and node_id[0] == "clone":
            continue
        desc = descendants(node_id)
        if 0 < len(desc) < len(leaf_labels):
            bipartitions.add(frozenset(desc))
    return bipartitions


def ancestor_relation_accuracy(true_edges, pred_edges, true_attachments, pred_attachments):
    """Compare rooted ancestor/descendant relationships among matched clones."""
    def parent_map(edges):
        return {child: parent for parent, child in edges}

    def is_ancestor(ancestor, descendant, parents):
        current = descendant
        while current in parents:
            current = parents[current]
            if current == ancestor:
                return True
        return False

    true_parents = parent_map(true_edges)
    pred_parents = parent_map(pred_edges)
    clones = sorted(set(true_attachments) & set(pred_attachments))
    total = 0
    correct = 0

    for i, clone_a in enumerate(clones):
        for clone_b in clones[i + 1:]:
            true_a = true_attachments[clone_a]
            true_b = true_attachments[clone_b]
            pred_a = pred_attachments[clone_a]
            pred_b = pred_attachments[clone_b]

            true_rel = (
                is_ancestor(true_a, true_b, true_parents),
                is_ancestor(true_b, true_a, true_parents),
            )
            pred_rel = (
                is_ancestor(pred_a, pred_b, pred_parents),
                is_ancestor(pred_b, pred_a, pred_parents),
            )
            total += 1
            correct += int(true_rel == pred_rel)

    return float(correct / total) if total else float("nan")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, required=True)
    parser.add_argument("--scicone-output-dir", type=Path, required=True)
    parser.add_argument("--metrics-out", type=Path, default=None)
    parser.add_argument("--tolerance", type=int, default=2)
    return parser.parse_args()


def main():
    args = parse_args()
    metrics_out = args.metrics_out or args.scicone_output_dir / "metrics.json"

    print(f"Loading ground truth from {args.dataset_dir}")
    dataset = load_dataset(args.dataset_dir)
    truth = dataset["truth_cell_hscn_segments"]
    cells = dataset["cells"]
    segments = dataset["segments"]
    bins = dataset["bins"]

    n_cells = len(cells)
    n_segs = len(segments)
    clone_labels = cells["clone_assignment"].values
    cell_names = cells["cell"].tolist()

    # True total CN per cell per segment
    tcn_true = np.zeros((n_cells, n_segs), dtype=int)
    for _, row in truth.iterrows():
        ci = int(row["cell_index"])
        si = int(row["segment"])
        tcn_true[ci, si] = int(row["cn_A"]) + int(row["cn_B"])

    # True breakpoints
    true_bps = sorted([int(seg["start_bin"]) for _, seg in segments.iterrows() if seg["segment"] > 0])

    # ── 1. Load SCICoNE cell_node_labels ──────────────────────────────────────
    print("Loading SCICoNE output...")
    labels_file = args.scicone_output_dir / "cell_node_labels.csv"
    if not labels_file.exists():
        raise FileNotFoundError("cell_node_labels.csv not found — run run_scicone.py first")

    node_labels = pd.read_csv(labels_file, header=None).iloc[:, 0].tolist()
    print(f"  Cell node labels: {len(node_labels)} cells, nodes: {sorted(set(node_labels))}")

    # Convert string labels to integers
    unique_nodes = sorted(set(node_labels))
    node_to_int = {n: i for i, n in enumerate(unique_nodes)}
    pred_labels = np.array([node_to_int[l] for l in node_labels])

    # ── 2. Clone ARI/NMI using SCICoNE's own labels ───────────────────────────
    ari = float(adjusted_rand_score(clone_labels, pred_labels))
    nmi = float(normalized_mutual_info_score(clone_labels, pred_labels))

    # ── 3. TCN MSE from inferred CNVs ─────────────────────────────────────────
    cnv_file = args.scicone_output_dir / "inferred_cnvs.csv"
    tcn_pred_bins = np.loadtxt(cnv_file, delimiter=",", dtype=int)
    print(f"  CNV matrix: {tcn_pred_bins.shape}")

    n_bins_pred = tcn_pred_bins.shape[1]
    tcn_pred = np.zeros((n_cells, n_segs), dtype=int)
    for si, seg_row in segments.iterrows():
        start = int(seg_row["start_bin"])
        end = int(seg_row["end_bin"])
        if start < n_bins_pred and end < n_bins_pred:
            tcn_pred[:, si] = np.round(
                np.mean(tcn_pred_bins[:, start:end+1], axis=1)).astype(int)
        else:
            tcn_pred[:, si] = 2

    tcn_mse = float(np.mean((tcn_pred - tcn_true) ** 2))

    # ── 4. Breakpoints from SCICoNE's segmented_regions ──────────────────────
    seg_regions_file = args.scicone_output_dir / "segmented_regions.txt"
    if seg_regions_file.exists():
        seg_regions = np.loadtxt(seg_regions_file, dtype=int)
        pred_bps = sorted(seg_regions.tolist()) if seg_regions.ndim > 0 else []
    else:
        # Fallback: detect from CN changes
        pred_bps = set()
        for ci in range(tcn_pred_bins.shape[0]):
            for b in range(1, tcn_pred_bins.shape[1]):
                if tcn_pred_bins[ci, b] != tcn_pred_bins[ci, b-1]:
                    pred_bps.add(b)
        pred_bps = sorted(pred_bps)

    # Breakpoint F1
    tol = args.tolerance
    matched_true = set()
    matched_pred = set()
    for pb in pred_bps:
        for tb in true_bps:
            if abs(pb - tb) <= tol and tb not in matched_true:
                matched_true.add(tb)
                matched_pred.add(pb)
                break

    tp = len(matched_true)
    precision = tp / len(pred_bps) if pred_bps else 0.0
    recall = tp / len(true_bps) if true_bps else 0.0
    f1 = 2*precision*recall/(precision+recall) if (precision+recall) > 0 else 0.0

    # ── 5. Tree metrics with reconciliation ───────────────────────────────────
    tree_metrics = {}
    tree_str_file = args.scicone_output_dir / "tree_str.txt"
    tree_struct_file = args.dataset_dir / "tree_structure.json"

    if tree_str_file.exists() and tree_struct_file.exists():
        with open(tree_struct_file) as handle:
            tree_struct = json.load(handle)
        true_edges = [(n["parent_id"], n["node_id"])
                      for n in tree_struct["nodes"] if n["parent_id"] != -1]
        true_node_ids = {n["node_id"] for n in tree_struct["nodes"]}

        pred_parents, pred_events, pred_edges = parse_scicone_tree_str(tree_str_file)
        region_sizes_file = args.scicone_output_dir / "segmented_region_sizes.txt"
        if not region_sizes_file.exists():
            raise FileNotFoundError(
                "segmented_region_sizes.txt is required to align SCICoNE labels "
                "with tree_str.txt nodes"
            )
        region_sizes = np.atleast_1d(np.loadtxt(region_sizes_file, dtype=int))

        tree_profiles = cumulative_scicone_profiles(
            pred_parents,
            pred_events,
            n_regions=len(region_sizes),
        )
        label_profiles = modal_label_profiles(node_labels, tcn_pred_bins, region_sizes)
        label_to_tree_node, label_tree_l1 = match_scicone_labels_to_tree_nodes(
            label_profiles,
            tree_profiles,
        )

        label_to_true_clone, match_acc = positive_hungarian_mapping(
            clone_labels,
            node_labels,
        )

        # K = truth clone labels observed in cells, including normal/root cells.
        K = len(set(clone_labels))
        K_hat = len(unique_nodes)

        true_attachments = {}
        pred_attachments = {}
        for label, true_clone in label_to_true_clone.items():
            if true_clone not in true_node_ids or label not in label_to_tree_node:
                continue
            true_attachments[int(true_clone)] = int(true_clone)
            pred_attachments[int(true_clone)] = int(label_to_tree_node[label])

        matched_true_clones = set(true_attachments)
        tree_coverage = len(matched_true_clones) / K if K else 0.0

        true_bip = attached_leaf_bipartitions(true_edges, true_attachments)
        pred_bip = attached_leaf_bipartitions(pred_edges, pred_attachments)
        rf_raw = len(true_bip.symmetric_difference(pred_bip))
        rf_max = max(1, len(true_bip) + len(pred_bip))
        nrf = rf_raw / rf_max
        ancestor_acc = ancestor_relation_accuracy(
            true_edges,
            pred_edges,
            true_attachments,
            pred_attachments,
        )

        tree_metrics = {
            "nrf_distance": nrf,
            "rf_raw": rf_raw,
            "rf_max": rf_max,
            "tree_coverage": tree_coverage,
            "cell_node_match_accuracy": match_acc,
            "ancestor_relation_accuracy": ancestor_acc,
            "n_pred_nodes": K_hat,
            "n_true_clones": K,
            "n_tree_matched_clones": len(matched_true_clones),
            "scicone_label_to_tree_node": {
                str(label): int(node) for label, node in label_to_tree_node.items()
            },
            "scicone_label_to_true_clone": {
                str(label): int(clone) for label, clone in label_to_true_clone.items()
            },
            "scicone_label_tree_profile_l1": label_tree_l1,
        }

    print(f"\n{'='*50}")
    print("SCICONE EVALUATION RESULTS")
    print(f"{'='*50}")
    print(f"  TCN MSE:          {tcn_mse:.4f}  (0=perfect)")
    print(f"  Clone ARI:        {ari:.4f}  (1=perfect, from cell_node_labels)")
    print(f"  Clone NMI:        {nmi:.4f}  (1=perfect)")
    print(f"  Breakpoint F1:    {f1:.4f}  (P={precision:.3f} R={recall:.3f})")
    print(f"  True breakpoints: {len(true_bps)}")
    print(f"  Pred breakpoints: {len(pred_bps)}")
    if tree_metrics:
        print(f"  nRF Distance:     {tree_metrics['nrf_distance']:.4f}  (0=perfect)")
        print(f"  TreeCoverage:     {tree_metrics['tree_coverage']:.4f}  (1=perfect)")
        print(f"  Node match acc:   {tree_metrics['cell_node_match_accuracy']:.4f}")
        print(f"  Ancestor acc:     {tree_metrics['ancestor_relation_accuracy']:.4f}")
        print(
            f"  Pred/True nodes:  {tree_metrics['n_pred_nodes']}/"
            f"{tree_metrics['n_true_clones']} "
            f"({tree_metrics['n_tree_matched_clones']} matched)"
        )
    print(f"{'='*50}")

    metrics = {
        "tcn_mse": tcn_mse,
        "clone_ari": ari,
        "clone_nmi": nmi,
        "breakpoint_precision": precision,
        "breakpoint_recall": recall,
        "breakpoint_f1": f1,
        "n_true_breakpoints": len(true_bps),
        "n_pred_breakpoints": len(pred_bps),
        **tree_metrics
    }
    with open(metrics_out, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nWrote metrics to {metrics_out}")


if __name__ == "__main__":
    main()
