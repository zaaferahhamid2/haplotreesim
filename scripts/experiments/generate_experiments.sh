#!/bin/bash
# generate_experiments.sh
# Run from haplotreesim root directory.

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

SEEDS=(42 43 44 45 46)

echo "=== 1. Clone number (generates reference tree structures) ==="
for SEED in "${SEEDS[@]}"; do
    for K in 2 4 8; do
        DIR="$BASE_DIR/clone_${K}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $K --num-cells $DEFAULT_CELLS \
            --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
            --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
            --random-seed $SEED --overwrite
        echo "  Created $DIR (tree: $DIR/tree_structure.json)"
    done
done

echo "=== 2. Tree balance ==="
for SEED in "${SEEDS[@]}"; do
    for BETA in 0.1 0.3 0.5; do
        BETA_TAG=$(echo $BETA | tr '.' '_')
        DIR="$BASE_DIR/beta_${BETA_TAG}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
            --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
            --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $BETA \
            --random-seed $SEED --overwrite
    done
done

echo "=== 3. Event rate (reuse tree from clones_4_rep{seed}) ==="
for SEED in "${SEEDS[@]}"; do
    TREE="$BASE_DIR/clone_4_rep${SEED}/tree_structure.json"
    for E in 10 80; do
        DIR="$BASE_DIR/events_${E}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
            --lambda-events $E --lambda-amplitude $DEFAULT_AMP \
            --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
            --tree-structure $TREE \
            --random-seed $SEED --overwrite
    done
done

echo "=== 4. Cell count (reuse tree from clones_4_rep{seed}) ==="
for SEED in "${SEEDS[@]}"; do
    TREE="$BASE_DIR/clone_4_rep${SEED}/tree_structure.json"
    for N in 50 500; do
        DIR="$BASE_DIR/cells_${N}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $DEFAULT_CLONES --num-cells $N \
            --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
            --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
            --tree-structure $TREE \
            --random-seed $SEED --overwrite
    done
done

echo "=== 5. Normal cell fraction (reuse tree from clones_4_rep{seed}) ==="
for SEED in "${SEEDS[@]}"; do
    TREE="$BASE_DIR/clone_4_rep${SEED}/tree_structure.json"
    for NRM in 0.0 0.3; do
        NRM_TAG=$(echo $NRM | tr '.' '_')
        DIR="$BASE_DIR/normal_${NRM_TAG}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
            --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
            --prob-normal $NRM --prob-focal $DEFAULT_FOCAL \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
            --tree-structure $TREE \
            --random-seed $SEED --overwrite
    done
done

echo "=== 6. Event amplitude (reuse tree from clones_4_rep{seed}) ==="
for SEED in "${SEEDS[@]}"; do
    TREE="$BASE_DIR/clone_4_rep${SEED}/tree_structure.json"
    for AMP in 1.0 4.0; do
        AMP_TAG=$(echo $AMP | tr '.' '_')
        DIR="$BASE_DIR/amplitude_${AMP_TAG}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
            --lambda-events $DEFAULT_EVENTS --lambda-amplitude $AMP \
            --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
            --tree-structure $TREE \
            --random-seed $SEED --overwrite
    done
done

echo "=== 7. Focal vs broad events (reuse tree from clones_4_rep{seed}) ==="
for SEED in "${SEEDS[@]}"; do
    TREE="$BASE_DIR/clone_4_rep${SEED}/tree_structure.json"
    for F in 0.3 1.0; do
        F_TAG=$(echo $F | tr '.' '_')
        DIR="$BASE_DIR/focal_${F_TAG}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
            --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
            --prob-normal $DEFAULT_NORMAL --prob-focal $F \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
            --tree-structure $TREE \
            --random-seed $SEED --overwrite
    done
done

echo "=== 8. Coverage (read-depth; reuse tree from clones_4_rep{seed}) ==="
for SEED in "${SEEDS[@]}"; do
    TREE="$BASE_DIR/clone_4_rep${SEED}/tree_structure.json"
    for COV in 25 200; do
        DIR="$BASE_DIR/coverage_${COV}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
            --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
            --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
            --mean-library-size $COV \
            --tree-structure $TREE \
            --random-seed $SEED --overwrite
    done
done

echo "=== 9. Phase switch error rate (reuse tree from clones_4_rep{seed}) ==="
for SEED in "${SEEDS[@]}"; do
    TREE="$BASE_DIR/clone_4_rep${SEED}/tree_structure.json"
    for PS in 0.001 0.05; do
        PS_TAG=$(echo $PS | tr '.' '_')
        DIR="$BASE_DIR/phaseswitch_${PS_TAG}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
            --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
            --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
            --phase-switch-prob $PS \
            --tree-structure $TREE \
            --random-seed $SEED --overwrite
    done
done

echo "=== 10. WGD probability (reuse tree from clones_4_rep{seed}) ==="
for SEED in "${SEEDS[@]}"; do
    TREE="$BASE_DIR/clone_4_rep${SEED}/tree_structure.json"
    for WGD in 0.0 0.4; do
        WGD_TAG=$(echo $WGD | tr '.' '_')
        DIR="$BASE_DIR/wgd_${WGD_TAG}_rep${SEED}"
        $SIM --output-dir $DIR \
            --whole-genome --num-clones $DEFAULT_CLONES --num-cells $DEFAULT_CELLS \
            --lambda-events $DEFAULT_EVENTS --lambda-amplitude $DEFAULT_AMP \
            --prob-normal $DEFAULT_NORMAL --prob-focal $DEFAULT_FOCAL \
            --alpha-tree $DEFAULT_ALPHA --beta-tree $DEFAULT_BETA \
            --prob-wgd $WGD \
            --tree-structure $TREE \
            --random-seed $SEED --overwrite
    done
done

echo ""
echo "Done. Created $(ls $BASE_DIR | wc -l) experiment directories."
ls $BASE_DIR
