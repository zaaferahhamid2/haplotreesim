#!/usr/bin/env python3
"""
Aggregate per-tool benchmark metrics into summary CSV files.

This reads result files produced by the per-tool evaluation scripts:

    results/<Tool>/<dataset>/metrics.json
    results/<Tool>/<dataset>/metrics.json.runtime

and writes:

    results/master_summary.csv
    results/metrics_summary_full.csv
    results/metrics_summary.csv
    results/runtime_summary.csv

Dataset parameter values are inferred from names such as clone_4_rep42,
coverage_0_005_rep42, and wgd_1_rep42 using the same defaults as
scripts/experiments/generate_experiments.sh.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


TOOL_ORDER = ["SEACON", "Alleloscope", "SCICoNE", "CONET", "CNRein"]

DEFAULT_PARAMS = {
    "clones": 4,
    "beta": 0.5,
    "lambdaE": 40.0,
    "N": 200,
    "pNormal": 0.1,
    "lambdaDelta": 2.0,
    "focal": 0.7,
    "mean_library_size": 100.0,
    "mean_allelic_coverage": 0.2,
    "phase_switch": 0.01,
    "wgd": 0.0,
    "doublet": 0.0,
}

COVERAGE_VALUES = {
    "0_005": (16.67, 0.005),
    "0_02": (66.67, 0.02),
    "0_05": (166.67, 0.05),
    "0_1": (333.33, 0.1),
}

PARAM_COLUMNS = [
    "rep",
    "clones",
    "beta",
    "lambdaE",
    "N",
    "pNormal",
    "lambdaDelta",
    "focal",
    "mean_library_size",
    "mean_allelic_coverage",
    "phase_switch",
    "wgd",
    "doublet",
]

MASTER_METRIC_COLUMNS = [
    "hscn_error",
    "tcn_mse",
    "clone_ari",
    "clone_nmi",
    "loh_f1",
    "loh_precision",
    "loh_recall",
    "msr",
    "breakpoint_f1",
    "breakpoint_precision",
    "breakpoint_recall",
    "nrf_distance",
    "tree_coverage",
    "cell_node_match_accuracy",
    "retention_rate",
]

METRICS_FULL_COLUMNS = [
    "breakpoint_f1",
    "breakpoint_precision",
    "breakpoint_recall",
    "cell_node_match_accuracy",
    "clone_ari",
    "clone_nmi",
    "hscn_error",
    "loh_f1",
    "loh_precision",
    "loh_recall",
    "msr",
    "nrf_distance",
    "retention_rate",
    "tcn_mse",
    "tree_coverage",
]

METRICS_COMPACT_COLUMNS = [
    "breakpoint_f1",
    "clone_ari",
    "hscn_error",
    "loh_f1",
    "msr",
    "nrf_distance",
    "tcn_mse",
    "tree_coverage",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("results"),
        help="Directory containing per-tool result folders.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory for summary CSVs. Defaults to --results-dir.",
    )
    parser.add_argument(
        "--include-extra-metrics",
        action="store_true",
        help=(
            "Append scalar metrics not in the current master schema, such as "
            "ancestor_relation_accuracy."
        ),
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail if a metrics.json file is present but cannot be parsed.",
    )
    return parser.parse_args()


def parse_tag_number(text: str) -> float:
    return float(text.replace("_", "."))


def parse_dataset_name(dataset: str) -> dict[str, Any]:
    match = re.search(r"_rep(\d+)$", dataset)
    if match is None:
        raise ValueError(f"Dataset name does not end with _rep<seed>: {dataset}")

    params = dict(DEFAULT_PARAMS)
    params["rep"] = int(match.group(1))
    stem = dataset[: match.start()]

    if stem.startswith("clone_"):
        params["clones"] = int(stem.removeprefix("clone_"))
    elif stem.startswith("beta_"):
        params["beta"] = parse_tag_number(stem.removeprefix("beta_"))
    elif stem.startswith("events_"):
        params["lambdaE"] = float(stem.removeprefix("events_"))
    elif stem.startswith("cells_"):
        params["N"] = int(stem.removeprefix("cells_"))
    elif stem.startswith("normal_"):
        params["pNormal"] = parse_tag_number(stem.removeprefix("normal_"))
    elif stem.startswith("amplitude_"):
        params["lambdaDelta"] = parse_tag_number(stem.removeprefix("amplitude_"))
    elif stem.startswith("focal_"):
        params["focal"] = parse_tag_number(stem.removeprefix("focal_"))
    elif stem.startswith("coverage_"):
        tag = stem.removeprefix("coverage_")
        if tag in COVERAGE_VALUES:
            params["mean_library_size"], params["mean_allelic_coverage"] = COVERAGE_VALUES[tag]
        else:
            coverage = parse_tag_number(tag)
            params["mean_allelic_coverage"] = coverage
            params["mean_library_size"] = round(coverage * 3333.33, 2)
    elif stem.startswith("phaseswitch_"):
        params["phase_switch"] = parse_tag_number(stem.removeprefix("phaseswitch_"))
    elif stem.startswith("wgd_"):
        params["wgd"] = parse_tag_number(stem.removeprefix("wgd_"))
    elif stem.startswith("doublet_"):
        params["doublet"] = parse_tag_number(stem.removeprefix("doublet_"))
    else:
        raise ValueError(f"Unrecognized dataset prefix: {dataset}")

    return params


def load_metrics(path: Path, strict: bool) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        message = f"Skipping invalid metrics file {path}: {exc}"
        if strict:
            raise RuntimeError(message) from exc
        print(message, file=sys.stderr)
        return None


def scalar_metrics(metrics: dict[str, Any]) -> dict[str, Any]:
    scalars = {}
    for key, value in metrics.items():
        if isinstance(value, (dict, list, tuple)):
            continue
        scalars[key] = value
    return scalars


def parse_runtime(path: Path) -> float | None:
    if not path.exists():
        return None

    runtime = None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("RUNTIME_SECONDS="):
                runtime = float(line.split("=", 1)[1])
    return runtime


def format_value(value: Any) -> Any:
    if value is None:
        return "NA"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "NA"
    return value


def tool_sort_key(tool: str) -> tuple[int, str]:
    try:
        return (TOOL_ORDER.index(tool), tool)
    except ValueError:
        return (len(TOOL_ORDER), tool)


def collect_rows(results_dir: Path, strict: bool):
    rows = []
    runtime_rows = []
    skipped = []

    for tool_dir in sorted([p for p in results_dir.iterdir() if p.is_dir()], key=lambda p: tool_sort_key(p.name)):
        tool = tool_dir.name
        for dataset_dir in sorted([p for p in tool_dir.iterdir() if p.is_dir()], key=lambda p: p.name):
            dataset = dataset_dir.name
            metrics_path = dataset_dir / "metrics.json"
            metrics = load_metrics(metrics_path, strict=strict)
            if metrics is None:
                skipped.append(str(metrics_path))
                continue

            try:
                params = parse_dataset_name(dataset)
            except ValueError as exc:
                if strict:
                    raise
                print(f"Skipping {metrics_path}: {exc}", file=sys.stderr)
                skipped.append(str(metrics_path))
                continue

            row = {
                "method": tool,
                "tool": tool,
                "dataset": dataset,
                **params,
                **scalar_metrics(metrics),
            }
            rows.append(row)

            runtime = parse_runtime(dataset_dir / "metrics.json.runtime")
            if runtime is not None:
                runtime_rows.append({
                    "tool": tool,
                    "dataset": dataset,
                    "runtime_seconds": runtime,
                    "runtime_minutes": runtime / 60.0,
                })

    rows.sort(key=lambda r: (tool_sort_key(str(r["method"])), str(r["dataset"])))
    runtime_rows.sort(key=lambda r: (tool_sort_key(str(r["tool"])), str(r["dataset"])))
    return rows, runtime_rows, skipped


def extra_metric_columns(rows: list[dict[str, Any]]) -> list[str]:
    excluded = {"method", "tool", "dataset", *PARAM_COLUMNS}
    known = set(MASTER_METRIC_COLUMNS)
    extras = set()
    for row in rows:
        for key, value in row.items():
            if key in excluded or key in known:
                continue
            if isinstance(value, (dict, list, tuple)):
                continue
            extras.add(key)
    return sorted(extras)


def write_csv(path: Path, columns: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: format_value(row.get(column)) for column in columns})


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir or args.results_dir

    rows, runtime_rows, skipped = collect_rows(args.results_dir, strict=args.strict)
    extras = extra_metric_columns(rows) if args.include_extra_metrics else []

    master_columns = ["method", "dataset", *PARAM_COLUMNS, *MASTER_METRIC_COLUMNS, *extras]
    full_columns = ["tool", "dataset", *METRICS_FULL_COLUMNS, *extras]
    compact_columns = ["tool", "dataset", *METRICS_COMPACT_COLUMNS]

    full_rows = [{**row, "tool": row["method"]} for row in rows]
    compact_rows = full_rows

    write_csv(output_dir / "master_summary.csv", master_columns, rows)
    write_csv(output_dir / "metrics_summary_full.csv", full_columns, full_rows)
    write_csv(output_dir / "metrics_summary.csv", compact_columns, compact_rows)
    write_csv(
        output_dir / "runtime_summary.csv",
        ["tool", "dataset", "runtime_seconds", "runtime_minutes"],
        runtime_rows,
    )

    print(f"Wrote {len(rows)} metric rows to {output_dir / 'master_summary.csv'}")
    print(f"Wrote {len(runtime_rows)} runtime rows to {output_dir / 'runtime_summary.csv'}")
    if skipped:
        print(f"Skipped {len(skipped)} missing/invalid metrics files", file=sys.stderr)


if __name__ == "__main__":
    main()
