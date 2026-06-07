"""
Evaluate SCICoNE output against HaploTreeSim ground truth.
SCICoNE outputs total CN (not haplotype-specific), so we evaluate TCN-MSE,
clone ARI, and breakpoint accuracy.
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
from dataset_io import load_dataset


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, required=True)
    parser.add_argument("--scicone-output-dir", type=Path, required=True)
    parser.add_argument("--metrics-out", type=Path, default=None)
    parser.add_argument("--tolerance", type=int, default=2,
                        help="Breakpoint tolerance in bins")
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

    # True breakpoints (bin indices)
    true_bps = set()
    for _, seg in segments.iterrows():
        if seg["segment"] > 0:
            true_bps.add(int(seg["start_bin"]))
    true_bps = sorted(true_bps)

    # Load SCICoNE CNV output
    print("Loading SCICoNE output...")
    cnv_file = args.scicone_output_dir / "scicone_inferred_cnvs.csv"
    if not cnv_file.exists():
        raise FileNotFoundError(f"scicone_inferred_cnvs.csv not found in {args.scicone_output_dir}")

    cnvs = pd.read_csv(cnv_file, header=None)
    print(f"  CNV matrix: {cnvs.shape[0]} cells x {cnvs.shape[1]} regions")

    # Load region sizes to map regions back to bins
    reg_sizes_file = args.scicone_output_dir / "scicone_segmented_region_sizes.txt"
    reg_sizes = []
    with open(reg_sizes_file) as f:
        for line in f:
            line = line.strip()
            if line:
                reg_sizes.append(int(line))

    # SCICoNE outputs bin-level CNVs directly
    n_bins = cnvs.shape[1]
    tcn_pred_bins = cnvs.values.astype(int)

    # Aggregate to segments
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

    # Clone ARI using predicted CN profiles
    n_clusters = min(len(np.unique(clone_labels)),
                     len(np.unique(tcn_pred, axis=0)))
    if n_clusters >= 2:
        km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        pred_labels = km.fit_predict(tcn_pred)
        ari = float(adjusted_rand_score(clone_labels, pred_labels))
        nmi = float(normalized_mutual_info_score(clone_labels, pred_labels))
    else:
        ari = float('nan')
        nmi = float('nan')

    # Predicted breakpoints: where CN changes across bins
    pred_bps = set()
    for cell_idx in range(tcn_pred_bins.shape[0]):
        row = tcn_pred_bins[cell_idx]
        for b in range(1, len(row)):
            if row[b] != row[b-1]:
                pred_bps.add(b)
    pred_bps = sorted(pred_bps)

    # Breakpoint F1 with tolerance
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
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    print(f"\n{'='*50}")
    print("SCICONE EVALUATION RESULTS")
    print(f"{'='*50}")
    print(f"  TCN MSE:          {tcn_mse:.4f}  (0=perfect)")
    print(f"  Clone ARI:        {ari:.4f}  (1=perfect)")
    print(f"  Clone NMI:        {nmi:.4f}  (1=perfect)")
    print(f"  Breakpoint F1:    {f1:.4f}  (P={precision:.3f} R={recall:.3f})")
    print(f"  True breakpoints: {len(true_bps)}")
    print(f"  Pred breakpoints: {len(pred_bps)}")
    print(f"{'='*50}")
    print(f"\nTrue TCN sample (cell 0): {tcn_true[0,:10].tolist()}")
    print(f"Pred TCN sample (cell 0): {tcn_pred[0,:10].tolist()}")

    metrics = {
        "tcn_mse": tcn_mse,
        "clone_ari": ari,
        "clone_nmi": nmi,
        "breakpoint_precision": precision,
        "breakpoint_recall": recall,
        "breakpoint_f1": f1,
        "n_true_breakpoints": len(true_bps),
        "n_pred_breakpoints": len(pred_bps),
        "n_cells": n_cells,
        "n_segments": n_segs
    }
    with open(metrics_out, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nWrote metrics to {metrics_out}")


if __name__ == "__main__":
    main()
