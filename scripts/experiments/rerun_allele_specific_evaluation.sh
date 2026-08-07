#!/usr/bin/env bash
set -euo pipefail

# Rerun only the allele-specific tool evaluation steps.
# Defaults match scripts/experiments/generate_slurm_jobs.py.
#
# Usage:
#   bash scripts/experiments/rerun_allele_specific_evaluation.sh
#   bash scripts/experiments/rerun_allele_specific_evaluation.sh clone_4_rep42
#   TOOLS="SEACON" bash scripts/experiments/rerun_allele_specific_evaluation.sh
#   BASE_DIR=$PWD PYTHON_CMD="conda run -n haplotreesim python" \
#     bash scripts/experiments/rerun_allele_specific_evaluation.sh clone_4_rep42

BASE_DIR="${BASE_DIR:-/blue/lzhang.uwf/zh34.uwf/haplotreesim}"
EXP_DIR="${EXP_DIR:-$BASE_DIR/experiments}"
RESULTS_DIR="${RESULTS_DIR:-$BASE_DIR/results}"
SEACON_OUTPUT_ROOT="${SEACON_OUTPUT_ROOT:-$BASE_DIR/SEACON}"
ALLELOSCOPE_OUTPUT_ROOT="${ALLELOSCOPE_OUTPUT_ROOT:-$RESULTS_DIR/Alleloscope}"
PYTHON_CMD="${PYTHON_CMD:-python3}"
RSCRIPT_CMD="${RSCRIPT_CMD:-Rscript}"
TOOLS="${TOOLS:-SEACON Alleloscope}"
ALLOW_SKIPS="${ALLOW_SKIPS:-1}"
DATASET_FILTER="${1:-all}"

if [[ "$DATASET_FILTER" == "all" ]]; then
    DATASETS=()
    while IFS= read -r dataset; do
        DATASETS+=("$dataset")
    done < <(
        find "$EXP_DIR" -mindepth 1 -maxdepth 1 -type d -exec basename {} \; |
            grep -E '^(clone|beta|events|cells|normal|amplitude|focal|coverage|phaseswitch|wgd|doublet)_' |
            sort
    )
else
    DATASETS=("$DATASET_FILTER")
fi

if [[ "${#DATASETS[@]}" -eq 0 ]]; then
    echo "No experiment dataset directories found under $EXP_DIR" >&2
    exit 1
fi

cd "$BASE_DIR"

skipped=0
echo "Found ${#DATASETS[@]} dataset(s) to evaluate."
for dataset in "${DATASETS[@]}"; do
    dataset_dir="$EXP_DIR/$dataset"

    if [[ ! -d "$dataset_dir" ]]; then
        echo "SKIP $dataset: missing dataset dir $dataset_dir" >&2
        skipped=$((skipped + 1))
        continue
    fi

    for tool in $TOOLS; do
        case "$tool" in
            SEACON)
                seacon_output="$SEACON_OUTPUT_ROOT/$dataset"
                metrics_dir="$RESULTS_DIR/SEACON/$dataset"
                metrics_out="$metrics_dir/metrics.json"
                if [[ ! -f "$seacon_output/calls.tsv" && ! -f "$seacon_output/calls_flat.tsv" ]]; then
                    echo "SKIP SEACON $dataset: missing calls.tsv/calls_flat.tsv in $seacon_output" >&2
                    skipped=$((skipped + 1))
                    continue
                fi
                mkdir -p "$metrics_dir"
                echo "=== Evaluating SEACON: $dataset ==="
                $PYTHON_CMD "$BASE_DIR/scripts/evaluate_seacon.py" \
                    --dataset-dir "$dataset_dir" \
                    --seacon-output-dir "$seacon_output" \
                    --metrics-out "$metrics_out"
                ;;
            Alleloscope)
                alleloscope_output="$ALLELOSCOPE_OUTPUT_ROOT/$dataset"
                metrics_out="$alleloscope_output/metrics.json"
                if [[ ! -f "$alleloscope_output/genotype_values.tsv" ]]; then
                    echo "SKIP Alleloscope $dataset: missing genotype_values.tsv in $alleloscope_output" >&2
                    skipped=$((skipped + 1))
                    continue
                fi
                echo "=== Evaluating Alleloscope: $dataset ==="
                $RSCRIPT_CMD "$BASE_DIR/scripts/evaluate_alleloscope.R" \
                    --dataset-dir "$dataset_dir" \
                    --output-dir "$alleloscope_output" \
                    --metrics-out "$metrics_out"
                ;;
            *)
                echo "Unknown tool in TOOLS: $tool" >&2
                exit 1
                ;;
        esac
    done
done

if [[ "$skipped" -gt 0 ]]; then
    echo "Completed with $skipped skipped evaluation(s)." >&2
fi

if [[ "$skipped" -gt 0 && "$ALLOW_SKIPS" == "0" ]]; then
    exit 1
fi

echo "Allele-specific evaluation rerun complete."
