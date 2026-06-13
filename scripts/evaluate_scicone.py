"""
Evaluate SCICoNE output against HaploTreeSim ground truth.
Uses SCICoNE's own cell_node_labels, segmented_regions, and tree edges.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from dataset_io import load_dataset
from haplotreesim.metrics.tree_metrics import compute_robinson_foulds_distance
from haplotreesim.metrics.clone_metrics import hungarian_matching


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
    tree_edges_file = args.scicone_output_dir / "tree_edges.txt"
    tree_struct_file = args.dataset_dir / "tree_structure.json"

    if tree_edges_file.exists() and tree_struct_file.exists():
        import json as json_mod
        tree_struct = json_mod.load(open(tree_struct_file))
        true_edges = [(n["parent_id"], n["node_id"])
                      for n in tree_struct["nodes"] if n["parent_id"] != -1]
        true_leaves = [n["node_id"] for n in tree_struct["nodes"] if n["is_leaf"]]

        pred_edges_raw = []
        with open(tree_edges_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) == 2:
                    pred_edges_raw.append((int(parts[0]), int(parts[1])))

        # Reconcile: map SCICoNE nodes (A,B,C...) to true clone IDs
        # Use hungarian matching on cell assignments
        mapping, match_acc = hungarian_matching(clone_labels, pred_labels)
        # mapping: pred_int -> true_clone_id

        # Map pred node IDs to true clone leaf IDs
        pred_node_ids = sorted(set([p for p,c in pred_edges_raw] + [c for p,c in pred_edges_raw]))
        pred_id_map = {nid: i for i, nid in enumerate(pred_node_ids)}
        seq_pred_edges = [(pred_id_map[p], pred_id_map[c]) for p,c in pred_edges_raw]

        rf = compute_robinson_foulds_distance(true_edges, seq_pred_edges)
        tree_metrics = {
            "rf_distance": rf,
            "cell_node_match_accuracy": match_acc,
            "n_pred_nodes": len(unique_nodes),
            "n_true_nodes": len(tree_struct["nodes"])
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
        print(f"  RF Distance:      {tree_metrics['rf_distance']:.4f}  (0=perfect)")
        print(f"  Node match acc:   {tree_metrics['cell_node_match_accuracy']:.4f}")
        print(f"  Pred/True nodes:  {tree_metrics['n_pred_nodes']}/{tree_metrics['n_true_nodes']}")
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
