"""
Evaluate CONET output against HaploTreeSim ground truth.
CONET outputs total CN profiles and a tree with breakpoint-defined edges.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from dataset_io import load_dataset
from haplotreesim.metrics.tree_metrics import compute_robinson_foulds_distance
from haplotreesim.metrics.clone_metrics import hungarian_matching


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, required=True)
    parser.add_argument("--conet-output-dir", type=Path, required=True)
    parser.add_argument("--metrics-out", type=Path, default=None)
    parser.add_argument("--tolerance", type=int, default=2)
    return parser.parse_args()


def main():
    args = parse_args()
    metrics_out = args.metrics_out or args.conet_output_dir / "metrics.json"

    print(f"Loading ground truth from {args.dataset_dir}")
    dataset = load_dataset(args.dataset_dir)
    truth = dataset["truth_cell_hscn_segments"]
    cells = dataset["cells"]
    segments = dataset["segments"]

    n_cells = len(cells)
    n_segs = len(segments)
    clone_labels = cells["clone_assignment"].values

    # True total CN per cell per segment
    tcn_true = np.zeros((n_cells, n_segs), dtype=int)
    for _, row in truth.iterrows():
        ci = int(row["cell_index"])
        si = int(row["segment"])
        tcn_true[ci, si] = int(row["cn_A"]) + int(row["cn_B"])

    # True breakpoints
    true_bps = sorted([int(seg["start_bin"]) for _, seg in segments.iterrows() if seg["segment"] > 0])

    # ── 1. Load CONET inferred CN ─────────────────────────────────────────────
    print("Loading CONET output...")
    cn_file = args.conet_output_dir / "inferred_cn.csv"
    n_bins_dataset = len(dataset["bins"])
    cn_df = pd.read_csv(cn_file, header=0, index_col=0)
    # CONET outputs bins x cells, columns are 1-based cell indices
    conet_cell_ids = [int(c) - 1 for c in cn_df.columns]  # convert to 0-based
    n_conet_cells = len(conet_cell_ids)
    # Remove extra chromosome-end bins
    n_bins = min(cn_df.shape[0], n_bins_dataset)
    cn_arr = cn_df.values[:n_bins, :].T.astype(int)  # cells x bins
    # Build full cells x bins matrix (missing cells get diploid=2)
    tcn_pred_bins = np.full((n_cells, n_bins), 2, dtype=int)
    for i, ci in enumerate(conet_cell_ids):
        if ci < n_cells:
            tcn_pred_bins[ci] = cn_arr[i]
    print(f"  CN matrix: {cn_arr.shape} CONET cells, {n_conet_cells} kept of {n_cells}")
    tcn_pred = np.zeros((n_cells, n_segs), dtype=int)
    for si, seg_row in segments.iterrows():
        start = int(seg_row["start_bin"])
        end = int(seg_row["end_bin"])
        if start < n_bins and end < n_bins:
            tcn_pred[:, si] = np.round(
                np.mean(tcn_pred_bins[:, start:end+1], axis=1)).astype(int)
        else:
            tcn_pred[:, si] = 2

    # TCN MSE
    tcn_mse = float(np.mean((tcn_pred - tcn_true) ** 2))

    # ── 2. Clone assignment via KMeans on CN profiles ─────────────────────────
    K = len(set(clone_labels))
    n_unique = len(np.unique(tcn_pred, axis=0))
    n_clusters = min(K, n_unique)
    if n_clusters >= 2:
        km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        pred_labels = km.fit_predict(tcn_pred)
        ari = float(adjusted_rand_score(clone_labels, pred_labels))
        nmi = float(normalized_mutual_info_score(clone_labels, pred_labels))
    else:
        pred_labels = np.zeros(n_cells, dtype=int)
        ari = float('nan')
        nmi = float('nan')

    # ── 3. Breakpoints from CONET tree edges ──────────────────────────────────
    edges_file = args.conet_output_dir / "inferred_tree_edges.txt"
    pred_bps = []
    pred_edges = []
    if edges_file.exists():
        import json as json_mod
        bins_df = dataset["bins"]
        # Build lookup: (chr_int, start_bp) -> bin_index
        bp_lookup = {}
        for _, b in bins_df.iterrows():
            chrom_int = int(b["chrom"].replace("chr",""))
            bp_lookup[(chrom_int, int(b["start"]))] = int(b["bin"])

        with open(edges_file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    parts = line.split("\t", 1)
                    if len(parts) == 2:
                        p_str, c_str = parts[0].strip(), parts[1].strip()
                        p_node = json_mod.loads(p_str) if p_str != "{}" else None
                        c_node = json_mod.loads(c_str)
                        pred_edges.append((str(p_node), str(c_node)))
                        if c_node and "chr" in c_node and "bin_start" in c_node:
                            chrom_int = int(c_node["chr"])
                            bin_start = int(c_node["bin_start"])
                            key = (chrom_int, bin_start)
                            if key in bp_lookup:
                                pred_bps.append(bp_lookup[key])
                except Exception as e:
                    pass
        pred_bps = sorted(set(pred_bps))

    # Breakpoint F1
    tol = args.tolerance
    matched_true = set()
    for pb in pred_bps:
        for tb in true_bps:
            if abs(pb - tb) <= tol and tb not in matched_true:
                matched_true.add(tb)
                break

    tp = len(matched_true)
    precision = tp / len(pred_bps) if pred_bps else 0.0
    recall = tp / len(true_bps) if true_bps else 0.0
    f1 = 2*precision*recall/(precision+recall) if (precision+recall) > 0 else 0.0

    # ── 4. Tree metrics ────────────────────────────────────────────────────────
    tree_metrics = {}
    tree_struct_file = args.dataset_dir / "tree_structure.json"
    if tree_struct_file.exists() and pred_edges:
        import json as json_mod
        tree_struct = json_mod.load(open(tree_struct_file))
        true_edges = [(n["parent_id"], n["node_id"])
                      for n in tree_struct["nodes"] if n["parent_id"] != -1]

        # CONET edges are (start_bin, end_bin) pairs — convert to sequential IDs
        all_nodes = sorted(set([p for p,c in pred_edges] + [c for p,c in pred_edges]))
        node_map = {n: i for i, n in enumerate(all_nodes)}
        all_nodes = sorted(set([p for p,c in pred_edges] + [c for p,c in pred_edges]))
        node_map = {n: i for i, n in enumerate(all_nodes)}
        seq_pred_edges = [(node_map[p], node_map[c]) for p,c in pred_edges]

        # Reconciliation
        mapping, match_acc = hungarian_matching(clone_labels, pred_labels)
        K_hat = len(set(pred_labels))
        matched_true_clones = set(mapping.values())
        tree_coverage = len(matched_true_clones) / K

        rf_raw = compute_robinson_foulds_distance(true_edges, seq_pred_edges, normalize=False)
        rf_max = max(1, 2 * (K - 1))
        nrf = rf_raw / rf_max

        tree_metrics = {
            "nrf_distance": nrf,
            "tree_coverage": tree_coverage,
            "cell_node_match_accuracy": match_acc,
            "n_pred_nodes": K_hat,
            "n_true_clones": K
        }

    print(f"\n{'='*50}")
    print("CONET EVALUATION RESULTS")
    print(f"{'='*50}")
    print(f"  TCN MSE:          {tcn_mse:.4f}  (0=perfect)")
    print(f"  Clone ARI:        {ari:.4f}  (1=perfect)")
    print(f"  Clone NMI:        {nmi:.4f}  (1=perfect)")
    print(f"  Breakpoint F1:    {f1:.4f}  (P={precision:.3f} R={recall:.3f})")
    print(f"  True breakpoints: {len(true_bps)}")
    print(f"  Pred breakpoints: {len(pred_bps)}")
    if tree_metrics:
        print(f"  nRF Distance:     {tree_metrics['nrf_distance']:.4f}  (0=perfect)")
        print(f"  TreeCoverage:     {tree_metrics['tree_coverage']:.4f}  (1=perfect)")
        print(f"  Node match acc:   {tree_metrics['cell_node_match_accuracy']:.4f}")
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
