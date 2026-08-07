#!/usr/bin/env bash
set -euo pipefail

# Rerun only the CNRein evaluation step from existing CNRein outputs.
# Defaults match scripts/experiments/generate_slurm_jobs.py.
#
# Usage:
#   bash scripts/experiments/rerun_cnrein_evaluation.sh
#   bash scripts/experiments/rerun_cnrein_evaluation.sh clone_4_rep42
#   BASE_DIR=$PWD PYTHON_CMD="conda run -n haplotreesim python" \
#     bash scripts/experiments/rerun_cnrein_evaluation.sh clone_4_rep42

BASE_DIR="${BASE_DIR:-/blue/lzhang.uwf/zh34.uwf/haplotreesim}"
EXP_DIR="${EXP_DIR:-$BASE_DIR/experiments}"
CNREIN_OUT_DIR="${CNREIN_OUT_DIR:-$BASE_DIR/results/CNRein}"
EVAL_SCRIPT="${EVAL_SCRIPT:-$BASE_DIR/scripts/evaluate_cnrein.py}"
PYTHON_CMD="${PYTHON_CMD:-python3}"
ALLOW_SKIPS="${ALLOW_SKIPS:-1}"
DATASET_FILTER="${1:-all}"

if [[ ! -f "$EVAL_SCRIPT" ]]; then
    echo "Missing evaluation script: $EVAL_SCRIPT" >&2
    exit 1
fi

if [[ "$DATASET_FILTER" == "all" ]]; then
    DATASETS=()
    while IFS= read -r dataset; do
        DATASETS+=("$dataset")
    done < <(
        find "$CNREIN_OUT_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; |
            grep -E '^(clone|beta|events|cells|normal|amplitude|focal|coverage|phaseswitch|wgd|doublet)_' |
            sort
    )
else
    DATASETS=("$DATASET_FILTER")
fi

if [[ "${#DATASETS[@]}" -eq 0 ]]; then
    echo "No CNRein result directories found under $CNREIN_OUT_DIR" >&2
    exit 1
fi

cd "$BASE_DIR"

skipped=0
echo "Found ${#DATASETS[@]} CNRein dataset(s) to evaluate."
for dataset in "${DATASETS[@]}"; do
    dataset_dir="$EXP_DIR/$dataset"
    output_dir="$CNREIN_OUT_DIR/$dataset"
    metrics_out="$output_dir/metrics.json"
    initial_cna="$output_dir/binScale/initialCNA.npz"
    regions="$output_dir/binScale/regions.npz"
    pred_csv="$output_dir/finalPrediction/CNNaivePrediction.csv"

    if [[ ! -d "$dataset_dir" ]]; then
        echo "SKIP $dataset: missing dataset dir $dataset_dir" >&2
        skipped=$((skipped + 1))
        continue
    fi

    if [[ ! -f "$initial_cna" || ! -f "$regions" ]] && [[ ! -f "$pred_csv" ]]; then
        echo "SKIP $dataset: missing CNRein output. Expected $initial_cna + $regions, or $pred_csv" >&2
        skipped=$((skipped + 1))
        continue
    fi

    mkdir -p "$output_dir"
    echo "=== Evaluating CNRein: $dataset ==="
    $PYTHON_CMD "$EVAL_SCRIPT" \
        --dataset-dir "$dataset_dir" \
        --cnrein-output-dir "$output_dir" \
        --metrics-out "$metrics_out"
done

if [[ "$skipped" -gt 0 ]]; then
    echo "Completed with $skipped skipped evaluation(s)." >&2
fi

if [[ "$skipped" -gt 0 && "$ALLOW_SKIPS" == "0" ]]; then
    exit 1
fi

echo "CNRein evaluation rerun complete."
