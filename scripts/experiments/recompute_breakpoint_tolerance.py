#!/usr/bin/env python3
"""
Recompute breakpoint precision/recall/F1 across bin tolerances from existing outputs.

This does not rerun any tool and does not rerun full evaluation. It reads the
saved HaploTreeSim experiment files plus existing tool outputs under results/.

Default paths match the benchmark HPC layout used by the rerun scripts.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_BASE_DIR = Path("/blue/lzhang.uwf/zh34.uwf/haplotreesim")
DEFAULT_TOOLS = ["SEACON", "CNRein", "CONET", "SCICoNE"]
DEFAULT_TOLERANCES = [1, 3, 5, 10]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help=f"Project base directory. Default: {DEFAULT_BASE_DIR}",
    )
    parser.add_argument(
        "--experiments-dir",
        type=Path,
        default=None,
        help="Experiment directory. Default: <base-dir>/experiments.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help="Results directory. Default: <base-dir>/results.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV. Default: <results-dir>/breakpoint_tolerance_summary.csv.",
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=None,
        help="Datasets to process. Default: all clone_4_rep* experiment directories.",
    )
    parser.add_argument(
        "--tools",
        nargs="+",
        default=DEFAULT_TOOLS,
        choices=["SEACON", "CNRein", "CONET", "SCICoNE"],
        help=f"Tools to process. Default: {' '.join(DEFAULT_TOOLS)}.",
    )
    parser.add_argument(
        "--tolerances",
        type=int,
        nargs="+",
        default=DEFAULT_TOLERANCES,
        help=f"Breakpoint tolerances in bins. Default: {DEFAULT_TOLERANCES}.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on missing/invalid tool outputs instead of skipping them.",
    )
    return parser.parse_args()


def dataset_rep(dataset: str) -> int | None:
    match = re.search(r"_rep(\d+)$", dataset)
    return int(match.group(1)) if match else None


def internal_boundary_mask(bins: pd.DataFrame) -> np.ndarray:
    chroms = bins["chrom"].astype(str).to_numpy()
    if len(chroms) < 2:
        return np.array([], dtype=bool)
    return chroms[:-1] == chroms[1:]


def is_internal_boundary(bp: int, boundary_mask: np.ndarray) -> bool:
    return bp > 0 and bp - 1 < len(boundary_mask) and bool(boundary_mask[bp - 1])


def filter_internal_breakpoints(breakpoints, boundary_mask: np.ndarray) -> list[int]:
    filtered = {
        int(bp)
        for bp in breakpoints
        if pd.notna(bp) and is_internal_boundary(int(bp), boundary_mask)
    }
    return sorted(filtered)


def true_breakpoints(dataset_dir: Path, boundary_mask: np.ndarray) -> list[int]:
    segments = pd.read_csv(dataset_dir / "segments.tsv", sep="\t")
    return filter_internal_breakpoints(
        segments["start_bin"].to_numpy(dtype=int)[1:],
        boundary_mask,
    )


def bin_start_lookup(bins: pd.DataFrame) -> dict[tuple[str, int], int]:
    lookup = {}
    for row in bins.itertuples(index=False):
        chrom = str(row.chrom).replace("chr", "")
        lookup[(chrom, int(row.start))] = int(row.bin)
    return lookup


def predicted_seacon(result_dir: Path) -> tuple[list[int], str]:
    metrics_path = result_dir / "metrics.json"
    with open(metrics_path, encoding="utf-8") as handle:
        metrics = json.load(handle)
    if "pred_breakpoints" not in metrics:
        raise KeyError(f"{metrics_path} does not contain pred_breakpoints")
    return [int(bp) for bp in metrics["pred_breakpoints"]], str(metrics_path)


def predicted_cnrein(result_dir: Path, bins: pd.DataFrame) -> tuple[list[int], str]:
    pred_path = result_dir / "finalPrediction" / "CNNaivePrediction.csv"
    pred_df = pd.read_csv(pred_path)
    if pred_df.empty:
        return [], str(pred_path)

    lookup = bin_start_lookup(bins)
    first_cell = pred_df["Cell barcode"].iloc[0]
    regions = pred_df.loc[
        pred_df["Cell barcode"] == first_cell,
        ["Chromosome", "Start", "End"],
    ].to_numpy()

    breakpoints = []
    for idx in range(1, len(regions)):
        chrom_prev = str(int(regions[idx - 1][0]))
        chrom_curr = str(int(regions[idx][0]))
        if chrom_prev != chrom_curr:
            continue
        bin_idx = lookup.get((chrom_curr, int(regions[idx][1])))
        if bin_idx is not None:
            breakpoints.append(bin_idx)
    return breakpoints, str(pred_path)


def predicted_conet(result_dir: Path, bins: pd.DataFrame) -> tuple[list[int], str]:
    edges_path = result_dir / "inferred_tree_edges.txt"
    lookup = bin_start_lookup(bins)
    breakpoints = []

    with open(edges_path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            child = json.loads(parts[1])
            if child and "chr" in child and "bin_start" in child:
                key = (str(int(child["chr"])), int(child["bin_start"]))
                bin_idx = lookup.get(key)
                if bin_idx is not None:
                    breakpoints.append(bin_idx)
    return breakpoints, str(edges_path)


def predicted_scicone(result_dir: Path) -> tuple[list[int], str]:
    regions_path = result_dir / "segmented_regions.txt"
    values = np.atleast_1d(np.loadtxt(regions_path, dtype=int))
    if values.size == 0:
        return [], str(regions_path)
    return [int(bp) for bp in values.tolist()], str(regions_path)


def predicted_breakpoints(tool: str, result_dir: Path, bins: pd.DataFrame) -> tuple[list[int], str]:
    if tool == "SEACON":
        return predicted_seacon(result_dir)
    if tool == "CNRein":
        return predicted_cnrein(result_dir, bins)
    if tool == "CONET":
        return predicted_conet(result_dir, bins)
    if tool == "SCICoNE":
        return predicted_scicone(result_dir)
    raise ValueError(f"Unsupported tool: {tool}")


def breakpoint_metrics(true_bps: list[int], pred_bps: list[int], tolerance: int) -> dict[str, int | float]:
    true_sorted = sorted(int(bp) for bp in true_bps)
    pred_sorted = sorted(int(bp) for bp in pred_bps)

    if not true_sorted and not pred_sorted:
        return {
            "precision": 1.0,
            "recall": 1.0,
            "f1": 1.0,
            "tp": 0,
            "fp": 0,
            "fn": 0,
            "n_true": 0,
            "n_pred": 0,
        }

    candidate_pairs = []
    for pred_idx, pred_bp in enumerate(pred_sorted):
        for true_idx, true_bp in enumerate(true_sorted):
            distance = abs(pred_bp - true_bp)
            if distance <= tolerance:
                candidate_pairs.append((distance, pred_idx, true_idx))

    matched_pred = set()
    matched_true = set()
    for _, pred_idx, true_idx in sorted(candidate_pairs):
        if pred_idx in matched_pred or true_idx in matched_true:
            continue
        matched_pred.add(pred_idx)
        matched_true.add(true_idx)

    tp = len(matched_pred)
    fp = len(pred_sorted) - tp
    fn = len(true_sorted) - len(matched_true)
    precision = tp / (tp + fp) if (tp + fp) else 0.0
    recall = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    return {
        "precision": float(precision),
        "recall": float(recall),
        "f1": float(f1),
        "tp": int(tp),
        "fp": int(fp),
        "fn": int(fn),
        "n_true": len(true_sorted),
        "n_pred": len(pred_sorted),
    }


def discover_datasets(experiments_dir: Path) -> list[str]:
    datasets = [
        path.name
        for path in experiments_dir.glob("clone_4_rep*")
        if path.is_dir() and (path / "segments.tsv").exists() and (path / "bins.tsv").exists()
    ]
    return sorted(datasets)


def handle_error(message: str, strict: bool) -> None:
    if strict:
        raise RuntimeError(message)
    print(f"SKIP: {message}", file=sys.stderr)


def main() -> None:
    args = parse_args()
    base_dir = args.base_dir
    experiments_dir = args.experiments_dir or base_dir / "experiments"
    results_dir = args.results_dir or base_dir / "results"
    output_path = args.output or results_dir / "breakpoint_tolerance_summary.csv"
    datasets = args.datasets or discover_datasets(experiments_dir)

    if not datasets:
        raise SystemExit(f"No clone_4_rep* datasets found under {experiments_dir}")

    rows = []
    for dataset in datasets:
        dataset_dir = experiments_dir / dataset
        if not dataset_dir.exists():
            handle_error(f"{dataset}: missing experiment directory {dataset_dir}", args.strict)
            continue

        bins = pd.read_csv(dataset_dir / "bins.tsv", sep="\t")
        boundary_mask = internal_boundary_mask(bins)
        true_bps = true_breakpoints(dataset_dir, boundary_mask)

        for tool in args.tools:
            result_dir = results_dir / tool / dataset
            if not result_dir.exists():
                handle_error(f"{tool} {dataset}: missing result directory {result_dir}", args.strict)
                continue

            try:
                pred_raw, source_file = predicted_breakpoints(tool, result_dir, bins)
            except Exception as exc:
                handle_error(f"{tool} {dataset}: cannot extract breakpoints ({exc})", args.strict)
                continue

            pred_bps = filter_internal_breakpoints(pred_raw, boundary_mask)
            for tolerance in args.tolerances:
                metrics = breakpoint_metrics(true_bps, pred_bps, tolerance=tolerance)
                rows.append(
                    {
                        "method": tool,
                        "dataset": dataset,
                        "rep": dataset_rep(dataset),
                        "tolerance": tolerance,
                        "precision": metrics["precision"],
                        "recall": metrics["recall"],
                        "f1": metrics["f1"],
                        "tp": metrics["tp"],
                        "fp": metrics["fp"],
                        "fn": metrics["fn"],
                        "n_true": metrics["n_true"],
                        "n_pred": metrics["n_pred"],
                        "n_pred_raw": len(set(int(bp) for bp in pred_raw if pd.notna(bp))),
                        "source_file": source_file,
                    }
                )

    if not rows:
        raise SystemExit("No breakpoint rows generated.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "method",
        "dataset",
        "rep",
        "tolerance",
        "precision",
        "recall",
        "f1",
        "tp",
        "fp",
        "fn",
        "n_true",
        "n_pred",
        "n_pred_raw",
        "source_file",
    ]
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {output_path}")


if __name__ == "__main__":
    main()
