"""
Evaluate SEACON output against a saved general HaploTreeSim dataset.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from haplotreesim.metrics import (
    compute_all_hscn_metrics,
    compute_breakpoint_sensitivity_curve,
    compute_recurrent_breakpoint_metrics,
    compute_clone_assignment_metrics,
)
from dataset_io import load_dataset


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-dir",
        type=Path,
        required=True,
        help="General HaploTreeSim dataset directory.",
    )
    parser.add_argument(
        "--seacon-output-dir",
        type=Path,
        required=True,
        help="Directory containing SEACON calls.tsv/calls_flat.tsv.",
    )
    parser.add_argument("--breakpoint-tolerance", type=int, default=2)
    parser.add_argument("--breakpoint-threshold", type=float, default=0.05)
    parser.add_argument(
        "--breakpoint-thresholds",
        type=float,
        nargs="+",
        default=[0.01, 0.05, 0.10, 0.20],
        help="Recurrence thresholds for breakpoint sensitivity reporting.",
    )
    parser.add_argument("--clusters", type=int, default=None, help="KMeans clusters; default uses truth clone count.")
    parser.add_argument("--metrics-out", type=Path, default=None, help="Optional JSON metrics output path.")
    return parser.parse_args()


def parse_cn(value):
    if pd.isna(value):
        return None
    cn = str(value).strip().strip("()[]")
    if cn == "" or cn.upper() == "NA":
        return None
    parts = cn.split(",")
    if len(parts) != 2:
        return None
    return int(float(parts[0])), int(float(parts[1]))


def json_metric(value):
    if value is None:
        return None
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    return float(value)


def format_metric(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return "NA"
    return f"{float(value):.4f}"


def parse_flat_bin_index(column_name, fallback_idx: int, bins: pd.DataFrame) -> int:
    try:
        return int(column_name)
    except (TypeError, ValueError):
        pass

    text = str(column_name)
    if ":" in text:
        coord = text.split(":", 1)[1].split("-", 1)[0].replace(",", "")
        try:
            start = int(coord)
        except ValueError:
            return fallback_idx
        matches = bins.index[bins["start"] == start].tolist()
        if matches:
            return int(matches[0])
    return fallback_idx


def interval_to_bins(start: int, end: int, bins: pd.DataFrame, chrom: str | None = None) -> tuple[int, int]:
    starts = bins["start"].to_numpy(dtype=int)
    ends = bins["end"].to_numpy(dtype=int)
    overlaps = (starts <= end) & (ends >= start)
    if chrom is not None and "chrom" in bins.columns:
        overlaps &= bins["chrom"].astype(str).to_numpy() == str(chrom)
    overlapping = np.where(overlaps)[0]
    if len(overlapping) == 0:
        return -1, -2
    return int(overlapping[0]), int(overlapping[-1])


def chromosome_boundary_mask(bins: pd.DataFrame) -> np.ndarray:
    """Return True for adjacent-bin boundaries inside a chromosome."""
    chroms = bins["chrom"].astype(str).to_numpy()
    if len(chroms) < 2:
        return np.array([], dtype=bool)
    return chroms[:-1] == chroms[1:]


def true_cna_breakpoints_from_segments(segments: pd.DataFrame, bins: pd.DataFrame) -> np.ndarray:
    """Return segment starts excluding chromosome-start boundaries."""
    candidate_breakpoints = segments["start_bin"].to_numpy(dtype=int)[1:]
    boundary_mask = chromosome_boundary_mask(bins)
    return np.array([
        int(bp)
        for bp in candidate_breakpoints
        if bp > 0 and bp - 1 < len(boundary_mask) and boundary_mask[bp - 1]
    ], dtype=int)


def load_seacon_bin_calls(output_dir: Path, cell_names: list[str], bins: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, int]:
    n_cells = len(cell_names)
    n_bins = len(bins)
    cell_to_idx = {cell: idx for idx, cell in enumerate(cell_names)}

    cn_A = np.full((n_cells, n_bins), np.nan)
    cn_B = np.full((n_cells, n_bins), np.nan)
    parsed_rows = 0

    calls_flat_path = output_dir / "calls_flat.tsv"
    calls_path = output_dir / "calls.tsv"

    use_calls_flat = False
    calls_flat = None
    if calls_flat_path.exists():
        calls_flat = pd.read_csv(calls_flat_path, sep="\t", index_col=0)
        use_calls_flat = calls_flat.shape[1] == n_bins

    if use_calls_flat:
        for row_label, row in calls_flat.iterrows():
            cell = str(row_label)
            if cell not in cell_to_idx:
                continue
            cell_idx = cell_to_idx[cell]
            for fallback_bin_idx, (column_name, value) in enumerate(row.items()):
                bin_idx = parse_flat_bin_index(column_name, fallback_bin_idx, bins)
                if bin_idx < 0 or bin_idx >= n_bins:
                    continue
                parsed = parse_cn(value)
                if parsed is None:
                    continue
                cn_A[cell_idx, bin_idx] = parsed[0]
                cn_B[cell_idx, bin_idx] = parsed[1]
                parsed_rows += 1
    elif calls_path.exists():
        calls = pd.read_csv(calls_path, sep="\t")
        for _, row in calls.iterrows():
            cell = str(row["cell"])
            if cell not in cell_to_idx:
                continue
            parsed = parse_cn(row["CN"])
            if parsed is None:
                continue

            chrom = row["chrom"] if "chrom" in row.index else None
            start_bin, end_bin = interval_to_bins(int(row["start"]), int(row["end"]), bins, chrom=chrom)
            if start_bin > end_bin:
                continue
            cell_idx = cell_to_idx[cell]
            cn_A[cell_idx, start_bin:end_bin + 1] = parsed[0]
            cn_B[cell_idx, start_bin:end_bin + 1] = parsed[1]
            parsed_rows += end_bin - start_bin + 1
    else:
        raise FileNotFoundError(f"Expected calls_flat.tsv or calls.tsv in {output_dir}")

    return cn_A, cn_B, parsed_rows


def truth_segment_matrices(dataset: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    truth = dataset["truth_cell_hscn_segments"]
    cells = dataset["cells"]["cell"].tolist()
    segments = dataset["segments"]["segment"].astype(int).tolist()

    cell_to_idx = {cell: idx for idx, cell in enumerate(cells)}
    segment_to_idx = {segment: idx for idx, segment in enumerate(segments)}
    cn_A = np.zeros((len(cells), len(segments)), dtype=int)
    cn_B = np.zeros((len(cells), len(segments)), dtype=int)

    for _, row in truth.iterrows():
        cell_idx = cell_to_idx[row["cell"]]
        seg_idx = segment_to_idx[int(row["segment"])]
        cn_A[cell_idx, seg_idx] = int(row["cn_A"])
        cn_B[cell_idx, seg_idx] = int(row["cn_B"])

    clone_labels = dataset["cells"]["clone_assignment"].to_numpy(dtype=int)
    return cn_A, cn_B, clone_labels


def aggregate_prediction_to_segments(
    pred_A_bins: np.ndarray,
    pred_B_bins: np.ndarray,
    segments: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Map SEACON bin-level calls onto true simulator segments by majority state.

    This isolates HSCN/LOH state accuracy from breakpoint oversegmentation:
    for each cell and true segment, the predicted segment state is the most
    frequent finite (A, B) bin state inside that true segment.
    """
    n_cells = pred_A_bins.shape[0]
    n_segs = len(segments)
    pred_A = np.zeros((n_cells, n_segs), dtype=int)
    pred_B = np.zeros((n_cells, n_segs), dtype=int)

    for seg_idx, segment in enumerate(segments.itertuples(index=False)):
        start_bin = int(segment.start_bin)
        end_bin = int(segment.end_bin)
        pred_A_slice = pred_A_bins[:, start_bin:end_bin + 1]
        pred_B_slice = pred_B_bins[:, start_bin:end_bin + 1]

        for cell_idx in range(n_cells):
            states = []
            for a, b in zip(pred_A_slice[cell_idx], pred_B_slice[cell_idx]):
                if np.isfinite(a) and np.isfinite(b):
                    states.append((int(a), int(b)))

            if not states:
                continue

            counts = {}
            first_seen = {}
            for order, state in enumerate(states):
                counts[state] = counts.get(state, 0) + 1
                first_seen.setdefault(state, order)

            majority_state = max(
                counts,
                key=lambda state: (counts[state], -first_seen[state])
            )
            pred_A[cell_idx, seg_idx], pred_B[cell_idx, seg_idx] = majority_state

    return pred_A, pred_B


def main() -> None:
    args = parse_args()
    dataset = load_dataset(args.dataset_dir)
    cells = dataset["cells"]["cell"].tolist()
    segments = dataset["segments"]
    bins = dataset["bins"]

    print("Loading saved ground truth...")
    cn_A_true, cn_B_true, clone_labels = truth_segment_matrices(dataset)
    boundary_mask = chromosome_boundary_mask(bins)
    true_breakpoints = true_cna_breakpoints_from_segments(segments, bins)
    print(f"  {len(cells)} cells, {len(segments)} segments, {len(true_breakpoints)} true breakpoints")

    print("\nParsing SEACON output...")
    pred_A_bins, pred_B_bins, parsed_rows = load_seacon_bin_calls(args.seacon_output_dir, cells, bins)
    cn_A_pred, cn_B_pred = aggregate_prediction_to_segments(pred_A_bins, pred_B_bins, segments)
    print(f"  Parsed {parsed_rows} cell-bin CN states from SEACON output")

    print("\nComputing metrics...")
    hscn = compute_all_hscn_metrics(cn_A_true, cn_B_true, cn_A_pred, cn_B_pred, clone_labels)
    bp = compute_recurrent_breakpoint_metrics(
        true_breakpoints,
        pred_cn_A=pred_A_bins,
        pred_cn_B=pred_B_bins,
        threshold=args.breakpoint_threshold,
        tolerance=args.breakpoint_tolerance,
        boundary_mask=boundary_mask,
    )
    bp_curve = compute_breakpoint_sensitivity_curve(
        true_breakpoints,
        pred_cn_A=pred_A_bins,
        pred_cn_B=pred_B_bins,
        thresholds=args.breakpoint_thresholds,
        tolerance=args.breakpoint_tolerance,
        boundary_mask=boundary_mask,
    )
    pred_breakpoints = bp["pred_breakpoints"]
    print(
        f"  Detected {len(pred_breakpoints)} predicted breakpoints "
        f"at recurrence threshold tau={args.breakpoint_threshold:g}"
    )

    from sklearn.cluster import KMeans

    n_clusters = args.clusters or len(np.unique(clone_labels))
    tcn_profiles = cn_A_pred + cn_B_pred
    pred_labels = KMeans(n_clusters=n_clusters, random_state=42, n_init=10).fit_predict(tcn_profiles)
    clone_metrics = compute_clone_assignment_metrics(clone_labels, pred_labels)

    metrics = {
        "hscn_error": float(hscn["hscn_error"]),
        "loh_precision": json_metric(hscn["loh_precision"]),
        "loh_recall": json_metric(hscn["loh_recall"]),
        "loh_f1": json_metric(hscn["loh_f1"]),
        "loh_n_true": int(hscn["loh_n_true"]),
        "loh_n_pred": int(hscn["loh_n_pred"]),
        "msr": json_metric(hscn["msr"]),
        "n_mirrored_triples": int(hscn["n_mirrored_triples"]),
        "breakpoint_precision": float(bp["precision"]),
        "breakpoint_recall": float(bp["recall"]),
        "breakpoint_f1": float(bp["f1"]),
        "breakpoint_threshold": float(args.breakpoint_threshold),
        "breakpoint_n_true": int(bp["n_true"]),
        "breakpoint_n_pred": int(bp["n_pred"]),
        "breakpoint_sensitivity": {
            str(threshold): {
                "precision": float(curve_metrics["precision"]),
                "recall": float(curve_metrics["recall"]),
                "f1": float(curve_metrics["f1"]),
                "tp": int(curve_metrics["tp"]),
                "fp": int(curve_metrics["fp"]),
                "fn": int(curve_metrics["fn"]),
                "n_pred": int(curve_metrics["n_pred"]),
            }
            for threshold, curve_metrics in bp_curve.items()
        },
        "clone_ari": float(clone_metrics["ari"]),
        "clone_nmi": float(clone_metrics["nmi"]),
        "true_breakpoints": true_breakpoints.tolist(),
        "pred_breakpoints": pred_breakpoints.tolist(),
    }

    metrics_out = args.metrics_out or args.seacon_output_dir / "evaluation_metrics.json"
    with open(metrics_out, "w") as handle:
        json.dump(metrics, handle, indent=2)
        handle.write("\n")

    print("\n" + "=" * 50)
    print("SEACON EVALUATION RESULTS")
    print("=" * 50)
    print(f"  HSCN Error:       {metrics['hscn_error']:.4f}  (0=perfect)")
    print(f"  LOH F1:           {format_metric(metrics['loh_f1'])}  (NA if no true LOH)")
    print("n_mirrored_triples:", metrics["n_mirrored_triples"])
    print(f"  MSR:              {format_metric(metrics['msr'])}  (NA if no mirrored triples)")
    print(
        f"  Breakpoint F1:    {metrics['breakpoint_f1']:.4f}  "
        f"(tau={metrics['breakpoint_threshold']:.3f}, "
        f"P={metrics['breakpoint_precision']:.3f} R={metrics['breakpoint_recall']:.3f})"
    )
    print(f"  Clone ARI:        {metrics['clone_ari']:.4f}  (1=perfect)")
    print(f"  Clone NMI:        {metrics['clone_nmi']:.4f}  (1=perfect)")
    print("=" * 50)
    print(f"\nWrote metrics to {metrics_out}")
    print("\nTrue breakpoints: ", true_breakpoints)
    print("Pred breakpoints: ", pred_breakpoints[:20], "..." if len(pred_breakpoints) > 20 else "")
    print("\nTrue CN sample (cell 0):", list(zip(cn_A_true[0], cn_B_true[0])))
    print("Pred CN sample (cell 0):", list(zip(cn_A_pred[0], cn_B_pred[0])))


if __name__ == "__main__":
    main()
