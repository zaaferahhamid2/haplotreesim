"""
Evaluate CNRein (CNNaive) output against HaploTreeSim ground truth.
CNRein outputs haplotype-specific CN per cell per region.
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
from haplotreesim.metrics.clone_metrics import hungarian_matching


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, required=True)
    parser.add_argument("--cnrein-output-dir", type=Path, required=True)
    parser.add_argument("--metrics-out", type=Path, default=None)
    parser.add_argument("--tolerance", type=int, default=2)
    return parser.parse_args()


def cnrein_region_to_genomic_interval(region, bins_df):
    n_bins = len(bins_df)
    start_bin = int(region[1])
    end_bin = int(region[2])

    if start_bin < 0 or start_bin >= n_bins:
        raise ValueError(f"CNRein region start bin {start_bin} is outside 0..{n_bins - 1}")
    if end_bin <= start_bin:
        raise ValueError(f"CNRein region has non-positive length: start={start_bin}, end={end_bin}")

    last_bin = min(end_bin, n_bins) - 1
    start_row = bins_df.iloc[start_bin]
    end_row = bins_df.iloc[last_bin]
    if start_row["chrom"] != end_row["chrom"]:
        raise ValueError(
            "CNRein region crosses chromosome boundary: "
            f"{start_row['chrom']}:{start_bin} to {end_row['chrom']}:{last_bin}"
        )

    chrom = str(start_row["chrom"]).replace("chr", "")
    return [chrom, int(start_row["start"]), int(end_row["end"])]


def load_cnrein_predictions(cnrein_output_dir, bins_df, cells):
    initial_cna_path = cnrein_output_dir / "binScale" / "initialCNA.npz"
    regions_path = cnrein_output_dir / "binScale" / "regions.npz"
    if initial_cna_path.exists() and regions_path.exists():
        pred = np.load(initial_cna_path)["arr_0"]
        regions = np.load(regions_path)["arr_0"]
        if pred.ndim != 3 or pred.shape[2] < 2:
            raise ValueError(f"Unexpected CNRein prediction shape in {initial_cna_path}: {pred.shape}")
        if regions.ndim != 2 or regions.shape[1] < 3:
            raise ValueError(f"Unexpected CNRein regions shape in {regions_path}: {regions.shape}")
        if pred.shape[1] != regions.shape[0]:
            raise ValueError(
                f"Prediction has {pred.shape[1]} regions, but regions.npz has {regions.shape[0]}"
            )
        cell_list = cells["cell"].tolist()
        if pred.shape[0] != len(cell_list):
            raise ValueError(
                f"Prediction has {pred.shape[0]} cells, but dataset has {len(cell_list)}"
            )
        region_list = np.array(
            [cnrein_region_to_genomic_interval(region, bins_df) for region in regions],
            dtype=object,
        )
        return cell_list, region_list, pred[:, :, 0], pred[:, :, 1], str(initial_cna_path)

    pred_file = cnrein_output_dir / "finalPrediction/CNNaivePrediction.csv"
    pred_df = pd.read_csv(pred_file)
    print(f"  Predictions: {pred_df.shape}")

    cell_list = pred_df["Cell barcode"].unique().tolist()
    region_list = pred_df[pred_df["Cell barcode"] == cell_list[0]][["Chromosome","Start","End"]].values
    n_regions = len(region_list)
    n_pred_cells = len(cell_list)

    cn_A_pred_reg = np.zeros((n_pred_cells, n_regions), dtype=int)
    cn_B_pred_reg = np.zeros((n_pred_cells, n_regions), dtype=int)

    for ci, cell in enumerate(cell_list):
        cell_df = pred_df[pred_df["Cell barcode"] == cell]
        cn_A_pred_reg[ci] = cell_df["Haplotype 1"].values
        cn_B_pred_reg[ci] = cell_df["Haplotype 2"].values

    return cell_list, region_list, cn_A_pred_reg, cn_B_pred_reg, str(pred_file)


def main():
    args = parse_args()
    metrics_out = args.metrics_out or args.cnrein_output_dir / "metrics.json"

    print(f"Loading ground truth from {args.dataset_dir}")
    dataset = load_dataset(args.dataset_dir)
    truth = dataset["truth_cell_hscn_segments"]
    cells = dataset["cells"]
    segments = dataset["segments"]
    bins_df = dataset["bins"]

    n_cells = len(cells)
    n_segs = len(segments)
    clone_labels = cells["clone_assignment"].values
    n_bins = len(bins_df)
    bin_chroms = bins_df["chrom"].astype(str).to_numpy()

    def is_internal_breakpoint(bp):
        return bp > 0 and bp < n_bins and bin_chroms[bp - 1] == bin_chroms[bp]

    # True HSCN per cell per segment
    cn_A_true = np.zeros((n_cells, n_segs), dtype=int)
    cn_B_true = np.zeros((n_cells, n_segs), dtype=int)
    for _, row in truth.iterrows():
        ci, si = int(row["cell_index"]), int(row["segment"])
        cn_A_true[ci, si] = int(row["cn_A"])
        cn_B_true[ci, si] = int(row["cn_B"])
    tcn_true = cn_A_true + cn_B_true

    true_bps = sorted(
        int(seg["start_bin"])
        for _, seg in segments.iterrows()
        if seg["segment"] > 0 and is_internal_breakpoint(int(seg["start_bin"]))
    )

    # Load CNRein predictions
    print("Loading CNRein output...")
    cell_list, region_list, cn_A_pred_reg, cn_B_pred_reg, prediction_source = load_cnrein_predictions(
        args.cnrein_output_dir,
        bins_df,
        cells,
    )
    n_regions = len(region_list)
    n_pred_cells = len(cell_list)
    print(f"  Predictions: {n_pred_cells} cells x {n_regions} regions from {prediction_source}")

    # Map regions to segments
    # Build lookup: chrom, start -> bin index
    chrom_start_to_bin = {}
    for _, b in bins_df.iterrows():
        chrom_start_to_bin[(b["chrom"].replace("chr",""), int(b["start"]))] = int(b["bin"])

    cn_A_pred = np.full((n_cells, n_segs), 1, dtype=int)
    cn_B_pred = np.full((n_cells, n_segs), 1, dtype=int)

    bin_to_segment = np.full(n_bins, -1, dtype=int)
    for si, seg_row in segments.reset_index(drop=True).iterrows():
        start_bin = int(seg_row["start_bin"])
        end_bin = min(int(seg_row["end_bin"]), n_bins - 1)
        bin_to_segment[start_bin:end_bin + 1] = si

    cell_to_truth_idx = {
        str(cell): idx
        for idx, cell in enumerate(cells["cell"].tolist())
    }
    pred_cell_indices = np.array(
        [cell_to_truth_idx.get(str(cell), -1) for cell in cell_list],
        dtype=int,
    )
    valid_pred_cells = pred_cell_indices >= 0

    for ri, (chrom, start, end) in enumerate(region_list):
        key = (str(int(chrom)), int(start))
        bin_idx = chrom_start_to_bin.get(key, None)
        if bin_idx is None:
            continue
        si = bin_to_segment[bin_idx]
        if si < 0:
            continue
        cn_A_pred[pred_cell_indices[valid_pred_cells], si] = cn_A_pred_reg[valid_pred_cells, ri]
        cn_B_pred[pred_cell_indices[valid_pred_cells], si] = cn_B_pred_reg[valid_pred_cells, ri]

    tcn_pred = cn_A_pred + cn_B_pred

    # HSCN swap-invariant error
    err1 = np.mean(np.abs(cn_A_pred - cn_A_true) + np.abs(cn_B_pred - cn_B_true))
    err2 = np.mean(np.abs(cn_A_pred - cn_B_true) + np.abs(cn_B_pred - cn_A_true))
    hscn_error = float(min(err1, err2))

    # LOH F1
    loh_true = (cn_A_true == 0) | (cn_B_true == 0)
    loh_pred = (cn_A_pred == 0) | (cn_B_pred == 0)
    tp = np.sum(loh_true & loh_pred)
    fp = np.sum(~loh_true & loh_pred)
    fn = np.sum(loh_true & ~loh_pred)
    loh_prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    loh_rec = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    loh_f1 = 2*loh_prec*loh_rec/(loh_prec+loh_rec) if (loh_prec+loh_rec) > 0 else 0.0

    # TCN MSE
    tcn_mse = float(np.mean((tcn_pred - tcn_true) ** 2))

    # Clone ARI
    K = len(set(clone_labels))
    n_unique = len(np.unique(tcn_pred, axis=0))
    n_clusters = min(K, n_unique)
    if n_clusters >= 2:
        km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        pred_labels = km.fit_predict(tcn_pred)
        ari = float(adjusted_rand_score(clone_labels, pred_labels))
        nmi = float(normalized_mutual_info_score(clone_labels, pred_labels))
    else:
        ari = nmi = float('nan')

    # Breakpoints from region boundaries
    pred_bps = []
    for chrom, start, _ in region_list:
        key = (str(int(chrom)), int(start))
        bin_idx = chrom_start_to_bin.get(key, None)
        if bin_idx is not None and is_internal_breakpoint(bin_idx):
            pred_bps.append(bin_idx)
    pred_bps = sorted(set(pred_bps))

    tol = args.tolerance
    matched_true = set()
    for pb in pred_bps:
        for tb in true_bps:
            if abs(pb - tb) <= tol and tb not in matched_true:
                matched_true.add(tb)
                break
    tp_bp = len(matched_true)
    bp_prec = tp_bp / len(pred_bps) if pred_bps else 0.0
    bp_rec = tp_bp / len(true_bps) if true_bps else 0.0
    bp_f1 = 2*bp_prec*bp_rec/(bp_prec+bp_rec) if (bp_prec+bp_rec) > 0 else 0.0

    print(f"\n{'='*50}")
    print("CNREIN (CNNaive) EVALUATION RESULTS")
    print(f"{'='*50}")
    print(f"  HSCN Error:       {hscn_error:.4f}  (0=perfect)")
    print(f"  LOH F1:           {loh_f1:.4f}")
    print(f"  TCN MSE:          {tcn_mse:.4f}  (0=perfect)")
    print(f"  Clone ARI:        {ari:.4f}  (1=perfect)")
    print(f"  Clone NMI:        {nmi:.4f}  (1=perfect)")
    print(f"  Breakpoint F1:    {bp_f1:.4f}  (P={bp_prec:.3f} R={bp_rec:.3f})")
    print(f"  True breakpoints: {len(true_bps)}")
    print(f"  Pred breakpoints: {len(pred_bps)}")
    print(f"{'='*50}")

    metrics = {
        "hscn_error": hscn_error,
        "loh_f1": loh_f1,
        "tcn_mse": tcn_mse,
        "clone_ari": ari,
        "clone_nmi": nmi,
        "breakpoint_f1": bp_f1,
        "breakpoint_precision": bp_prec,
        "breakpoint_recall": bp_rec,
        "n_true_breakpoints": len(true_bps),
        "n_pred_breakpoints": len(pred_bps)
    }
    with open(metrics_out, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nWrote metrics to {metrics_out}")


if __name__ == "__main__":
    main()
