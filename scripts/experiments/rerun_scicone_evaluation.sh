#!/usr/bin/env bash
set -euo pipefail

# Rerun only the SCICoNE evaluation step.
# Defaults match scripts/experiments/generate_slurm_jobs.py.
#
# Usage:
#   bash scripts/experiments/rerun_scicone_evaluation.sh
#   bash scripts/experiments/rerun_scicone_evaluation.sh clone_4_rep42
#   BASE_DIR=$PWD bash scripts/experiments/rerun_scicone_evaluation.sh clone_4_rep42

BASE_DIR="${BASE_DIR:-/blue/lzhang.uwf/zh34.uwf/haplotreesim}"
EXP_DIR="${EXP_DIR:-$BASE_DIR/experiments}"
SCICONE_OUT_DIR="${SCICONE_OUT_DIR:-$BASE_DIR/results/SCICoNE}"
EVAL_SCRIPT="${EVAL_SCRIPT:-$BASE_DIR/scripts/evaluate_scicone.py}"
PYTHON_CMD="${PYTHON_CMD:-python3}"
DATASET_FILTER="${1:-all}"

if [[ ! -f "$EVAL_SCRIPT" ]]; then
    echo "Missing evaluation script: $EVAL_SCRIPT" >&2
    exit 1
fi

if [[ "$DATASET_FILTER" == "all" ]]; then
    mapfile -t DATASETS < <(
        find "$SCICONE_OUT_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; | sort
    )
else
    DATASETS=("$DATASET_FILTER")
fi

if [[ "${#DATASETS[@]}" -eq 0 ]]; then
    echo "No SCICoNE result directories found under $SCICONE_OUT_DIR" >&2
    exit 1
fi

cd "$BASE_DIR"

failed=0
for dataset in "${DATASETS[@]}"; do
    dataset_dir="$EXP_DIR/$dataset"
    output_dir="$SCICONE_OUT_DIR/$dataset"
    metrics_out="$output_dir/metrics.json"

    if [[ ! -d "$dataset_dir" ]]; then
        echo "SKIP $dataset: missing dataset dir $dataset_dir" >&2
        failed=$((failed + 1))
        continue
    fi

    if [[ ! -f "$output_dir/cell_node_labels.csv" ]]; then
        echo "SKIP $dataset: missing SCICoNE output $output_dir/cell_node_labels.csv" >&2
        failed=$((failed + 1))
        continue
    fi

    echo "=== Evaluating SCICoNE: $dataset ==="
    $PYTHON_CMD "$EVAL_SCRIPT" \
        --dataset-dir "$dataset_dir" \
        --scicone-output-dir "$output_dir" \
        --metrics-out "$metrics_out"
done

if [[ "$failed" -gt 0 ]]; then
    echo "Completed with $failed skipped/failed dataset(s)." >&2
    exit 1
fi

echo "SCICoNE evaluation rerun complete."
