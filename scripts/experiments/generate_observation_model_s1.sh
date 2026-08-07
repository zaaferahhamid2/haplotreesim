#!/usr/bin/env bash
set -euo pipefail

# Generate Supplementary Figure S1 from simulated experiment observation files.
# Defaults match the benchmark HPC layout used by the rerun scripts.
#
# Usage:
#   bash scripts/experiments/generate_observation_model_s1.sh
#   BASE_DIR=$PWD PYTHON_CMD="conda run -n haplotreesim python" \
#     bash scripts/experiments/generate_observation_model_s1.sh

BASE_DIR="${BASE_DIR:-/blue/lzhang.uwf/zh34.uwf/haplotreesim}"
PYTHON_CMD="${PYTHON_CMD:-python3}"
MPLCONFIGDIR="${MPLCONFIGDIR:-$BASE_DIR/.matplotlib-cache}"

mkdir -p "$MPLCONFIGDIR"
cd "$BASE_DIR"

echo "Generating Supplementary Figure S1 from $BASE_DIR/experiments"
MPLCONFIGDIR="$MPLCONFIGDIR" $PYTHON_CMD "$BASE_DIR/plots/generate_figures.py" --figure s1
