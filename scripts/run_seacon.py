"""
Convert HaploTreeSim output to SEACON input and run SEACON.
Skips prep_readcount (needs BAM) by writing files directly.
Uses --precomputed-baf to skip CHISEL BAF estimation.
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from haplotreesim import SimulationConfig, HaploTreeSimulator

# ── 1. Run simulator ──────────────────────────────────────────────────────────
print("Running simulator...")
config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=100,
    lambda_events=1.5,
    lambda_amplitude=1.0,
    prob_wgd=0.0,
    random_seed=42
)

sim = HaploTreeSimulator(config)
read_counts, allele_counts = sim.run()
alt_counts, ref_counts, total_counts = allele_counts
ground_truth = sim.get_ground_truth()
segments = ground_truth['segments']

print(f"  {read_counts.shape[0]} cells, {read_counts.shape[1]} bins, {len(segments)} segments")

# ── 2. Set up paths ───────────────────────────────────────────────────────────
seacon_dir = Path(os.environ.get('SEACON_DIR', Path.home() / 'Documents/SEACON'))
out_dir = Path('seacon_output')
out_dir.mkdir(exist_ok=True)

n_cells, n_bins = read_counts.shape
cell_names = [f'cell{i}' for i in range(n_cells)]
bin_size = config.bin_width

# Genomic coordinates for chr1 bins (1-based, matching SEACON convention)
bin_coords = []
for b in range(n_bins):
    start = b * bin_size + 1
    end = (b + 1) * bin_size
    bin_coords.append(('chr1', start, end))

# ── 3. Write cells.txt ────────────────────────────────────────────────────────
cells_file = out_dir / 'cells.txt'
with open(cells_file, 'w') as f:
    for c in cell_names:
        f.write(c + '\n')
print(f"  Wrote {cells_file}")

# ── 4. Write filtered_regions.bed ────────────────────────────────────────────
regions_file = out_dir / 'filtered_regions.bed'
with open(regions_file, 'w') as f:
    for chrom, start, end in bin_coords:
        f.write(f'{chrom}\t{start}\t{end}\n')
print(f"  Wrote {regions_file}")

# ── 5. Write readcounts.tsv (cells x bins, integer counts) ───────────────────
# SEACON flat format: index=cell, columns=0..J-1 (integer bin indices)
rc_df = pd.DataFrame(read_counts, index=cell_names,
                     columns=list(range(n_bins)))
rc_df.index.name = None
rc_file = out_dir / 'readcounts.tsv'
rc_df.to_csv(rc_file, sep='\t')
print(f"  Wrote {rc_file}")

# ── 6. Write RDR.tsv (read depth ratio normalized per cell) ──────────────────
cell_avgs = rc_df.mean(axis=1).replace(0, 1)
rdr_df = rc_df.div(cell_avgs, axis=0).round(5)
rdr_file = out_dir / 'RDR.tsv'
rdr_df.to_csv(rdr_file, sep='\t')
print(f"  Wrote {rdr_file}")

# ── 7. Write precomputed BAF file ─────────────────────────────────────────────
# Format: chrom  start  end  cell  BAF  (one row per cell per segment)
# BAF = alt / (alt + ref), mirrored to [0, 0.5]
baf_rows = []
for seg_idx, seg in enumerate(segments):
    seg_start_bp = seg.start_bin * bin_size + 1
    seg_end_bp   = (seg.end_bin + 1) * bin_size

    for cell_idx in range(n_cells):
        alt = alt_counts[cell_idx, seg_idx]
        ref = ref_counts[cell_idx, seg_idx]
        tot = alt + ref
        if tot == 0:
            baf = 0.5
        else:
            raw_baf = alt / tot
            baf = min(raw_baf, 1 - raw_baf)   # mirror to [0, 0.5]
        baf_rows.append({
            'chrom': 'chr1',
            'start': seg_start_bp,
            'end':   seg_end_bp,
            'cell':  cell_names[cell_idx],
            'BAF':   round(baf, 5)
        })

baf_df = pd.DataFrame(baf_rows, columns=['chrom', 'start', 'end', 'cell', 'BAF'])
baf_file = out_dir / 'precomputed_baf.tsv'
baf_df.to_csv(baf_file, sep='\t', index=False)
print(f"  Wrote {baf_file} ({len(baf_rows)} rows)")

# ── 8. Run seacon prep_baf then seacon call ───────────────────────────────────
print("\nRunning seacon prep_baf...")
ret = os.system(
    f'seacon prep_baf -o {out_dir} --precomputed-baf {baf_file} --no-normal 2>&1'
)
if ret != 0:
    print("  WARNING: prep_baf returned non-zero. Checking if BAF.tsv was created...")

baf_out = out_dir / 'BAF.tsv'
if not baf_out.exists():
    print("  BAF.tsv not found - copying precomputed BAF directly")
    # Pivot to SEACON's internal format: index=cell, columns=bin_idx
    # Map segment BAFs back to bins
    baf_bin_rows = {}
    for seg_idx, seg in enumerate(segments):
        for b in range(seg.start_bin, seg.end_bin + 1):
            for cell_idx in range(n_cells):
                cell = cell_names[cell_idx]
                if cell not in baf_bin_rows:
                    baf_bin_rows[cell] = {}
                baf_bin_rows[cell][b] = baf_rows[seg_idx * n_cells + cell_idx]['BAF']
    baf_internal = pd.DataFrame(baf_bin_rows).T
    baf_internal.columns = list(range(n_bins))
    baf_internal.to_csv(baf_out, sep='\t')
    print(f"  Wrote {baf_out}")

print("\nRunning seacon call...")
ret = os.system(
    f'seacon call -o {out_dir} --no-normal --upper-filter 5 --tolerance 0.15 --max-wgd 1 2>&1'
)

# ── 9. Check output ───────────────────────────────────────────────────────────
calls_file = out_dir / 'calls.tsv'
if calls_file.exists():
    calls = pd.read_csv(calls_file, sep='\t')
    print(f"\n✓ SEACON complete. calls.tsv has {len(calls)} rows.")
    print(calls.head(10).to_string())
else:
    print("\n✗ calls.tsv not found. Check seacon_output/ for logs.")
    for f in sorted(out_dir.iterdir()):
        print(f"  {f.name}")
