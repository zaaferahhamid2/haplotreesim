#!/bin/bash
# convert_all.sh
# Converts all experiment datasets to each tool's input format.
# Run from haplotreesim root directory.
# Usage: bash scripts/experiments/convert_all.sh

set -e
EXP_DIR="experiments"
CONV_DIR="."

module load conda 2>/dev/null || true
module load R/4.1 2>/dev/null || true

# Get all experiment directories
DATASETS=$(ls $EXP_DIR | grep -E "^(clone|beta|events|cells|normal|amplitude|focal|coverage|phaseswitch|wgd|doublet)_")

echo "Found $(echo $DATASETS | wc -w) datasets"

# ── SEACON ───────────────────────────────────────────────────────────────────
echo "=== Converting for SEACON ==="
source /apps/conda/25.7.0/etc/profile.d/conda.sh && conda activate haplotreesim
for DS in $DATASETS; do
    IN="$EXP_DIR/$DS"
    OUT="$CONV_DIR/SEACON/$DS"
    if [ -f "$OUT/BAF.tsv" ]; then
        echo "  Skipping $DS (already converted)"
        continue
    fi
    echo "  Converting $DS..."
    python3 scripts/convert_to_seacon.py \
        --input-dir $IN \
        --output-dir $OUT \
        --overwrite

    # Write BAF.tsv
    python3 - << EOF
import pandas as pd, numpy as np
sd = '$OUT'
rc = pd.read_csv(f'{sd}/readcounts.tsv', sep='\t', index_col=0)
baf_in = pd.read_csv(f'{sd}/precomputed_baf.tsv', sep='\t')
regions = pd.read_csv(f'{sd}/filtered_regions.bed', sep='\t', header=None, names=['chrom','start','end'])
regions['bin_idx'] = range(len(regions))
lookup = {}
for _, row in regions.iterrows():
    lookup[(row['chrom'], int(row['start']))] = int(row['bin_idx'])
    lookup[(row['chrom'], int(row['start'])-1)] = int(row['bin_idx'])
baf_in['bin_idx'] = baf_in.apply(lambda r: lookup.get((r['chrom'], int(r['start'])), lookup.get((r['chrom'], int(r['start'])+1), None)), axis=1)
baf_in = baf_in.dropna(subset=['bin_idx'])
baf_in['bin_idx'] = baf_in['bin_idx'].astype(int)
baf_mat = baf_in.pivot(index='cell', columns='bin_idx', values='BAF')
baf_mat = baf_mat.reindex(index=rc.index.tolist(), columns=range(rc.shape[1]), fill_value=0.5)
baf_mat.to_csv(f'{sd}/BAF.tsv', sep='\t')
print(f'  BAF.tsv written: {baf_mat.shape}')
EOF
done

# ── Alleloscope ──────────────────────────────────────────────────────────────
echo "=== Converting for Alleloscope ==="
for DS in $DATASETS; do
    IN="$EXP_DIR/$DS"
    OUT="$CONV_DIR/Alleloscope/$DS"
    if [ -d "$OUT" ] && [ -f "$OUT/var_all.vcf" ]; then
        echo "  Skipping $DS (already converted)"
        continue
    fi
    echo "  Converting $DS..."
    mkdir -p "$OUT"
    Rscript scripts/convert_to_alleloscope.R \
        --input-dir $IN \
        --output-dir $OUT 2>&1 | tail -3
done

# ── SCICoNE ──────────────────────────────────────────────────────────────────
echo "=== Converting for SCICoNE ==="
source /apps/conda/25.7.0/etc/profile.d/conda.sh && conda activate haplotreesim
for DS in $DATASETS; do
    IN="$EXP_DIR/$DS"
    OUT="$CONV_DIR/SCICoNE/$DS"
    if [ -f "$OUT/readcounts.csv" ]; then
        echo "  Skipping $DS (already converted)"
        continue
    fi
    echo "  Converting $DS..."
    python3 scripts/convert_to_scicone.py \
        --input-dir $IN \
        --output-dir $OUT \
        --overwrite
done

# ── CONET ────────────────────────────────────────────────────────────────────
echo "=== Converting for CONET ==="
source /apps/conda/25.7.0/etc/profile.d/conda.sh && conda activate CONET
for DS in $DATASETS; do
    IN="$EXP_DIR/$DS"
    OUT="$CONV_DIR/CONET/$DS"
    if [ -f "$OUT/corrected_counts.csv" ]; then
        echo "  Skipping $DS (already converted)"
        continue
    fi
    echo "  Converting $DS..."
    python3 scripts/convert_to_conet.py \
        --input-dir $IN \
        --output-dir $OUT \
        --overwrite
done

# ── CNRein ───────────────────────────────────────────────────────────────────
echo "=== Converting for CNRein ==="
source /apps/conda/25.7.0/etc/profile.d/conda.sh && conda activate haplotreesim
for DS in $DATASETS; do
    IN="$EXP_DIR/$DS"
    OUT="$CONV_DIR/CNRein/$DS"
    if [ -f "$OUT/initial/RDR_1M.npz" ]; then
        echo "  Skipping $DS (already converted)"
        continue
    fi
    echo "  Converting $DS..."
    python3 scripts/convert_to_cnrein.py \
        --input-dir $IN \
        --output-dir $OUT \
        --overwrite
done

echo ""
echo "Conversion complete."
