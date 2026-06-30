#!/bin/bash
# generate_experiments.sh
# Creates all 16 simulation datasets for the benchmark grid.
# Run from the haplotreesim root directory: bash scripts/experiments/generate_experiments.sh

set -e
BASE_DIR="experiments"
mkdir -p "$BASE_DIR"

SIM="python3 scripts/simulate_dataset.py"

# Shared defaults
DEFAULT_CLONES=4
DEFAULT_CELLS=200
DEFAULT_EVENTS=40
DEFAULT_AMP=2.0
DEFAULT_NORMAL=0.1
DEFAULT_FOCAL=0.7
DEFAULT_ALPHA=0.5
DEFAULT_BETA=0.5
SEED=42

echo "=== Default dataset ==="
$SIM --output-dir $BASE_DIR/default \
    --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
    --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
    --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
    --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
    --random-seed $SEED --overwrite

echo "=== 1. Clone number ==="
for K in 2 6 8; do
    $SIM --output-dir $BASE_DIR/clones_${K} \
        --whole-genome --num-clones $K --num-cells $DEFAULT_CELLS \
        --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
        --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
        --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
        --random-seed $SEED --overwrite
done

echo "=== 2. Tree balance ==="
$SIM --output-dir $BASE_DIR/balance_moderate \
    --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
    --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
    --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
    --alpha-tree 0.5 --beta-tree 0.3 \
    --random-seed $SEED --overwrite

$SIM --output-dir $BASE_DIR/balance_caterpillar \
    --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
    --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
    --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
    --alpha-tree 0.5 --beta-tree 0.1 \
    --random-seed $SEED --overwrite

echo "=== 3. Event rate ==="
for E in 10 80; do
    $SIM --output-dir $BASE_DIR/events_${E} \
        --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
        --lambda-events $E --lambda-amplitude $DEFAULT_AMP \
        --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
        --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
        --random-seed $SEED --overwrite
done

echo "=== 4. Cell count ==="
for N in 50 500; do
    $SIM --output-dir $BASE_DIR/cells_${N} \
        --whole-genome --num-clones $DEFAULT_CLONES --num-cells $N \
        --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
        --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
        --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
        --random-seed $SEED --overwrite
done

echo "=== 5. Normal cell fraction ==="
for NRM in 0.0 0.3; do
    NRM_TAG=$(echo $NRM | tr '.' '_')
    $SIM --output-dir $BASE_DIR/normal_${NRM_TAG} \
        --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
        --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
        --prob-normal $NRM --prob-focal $DEFAULT_FOCAL \
        --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
        --random-seed $SEED --overwrite
done

echo "=== 6. Event amplitude ==="
for AMP in 1.0 4.0; do
    AMP_TAG=$(echo $AMP | tr '.' '_')
    $SIM --output-dir $BASE_DIR/amplitude_${AMP_TAG} \
        --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
        --lambda-events $DEFAULT_EVENTS --lambda-amplitude $AMP \
        --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
        --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
        --random-seed $SEED --overwrite
done

echo "=== 7. Focal vs broad events ==="
for F in 0.3 1.0; do
    F_TAG=$(echo $F | tr '.' '_')
    $SIM --output-dir $BASE_DIR/focal_${F_TAG} \
        --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
        --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
        --prob-normal $DEFAULT_NORMAL --prob-focal $F \
        --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
        --random-seed $SEED --overwrite
done

echo ""
echo "Done. Created $(ls $BASE_DIR | wc -l) experiment directories in $BASE_DIR/"
ls $BASE_DIR
