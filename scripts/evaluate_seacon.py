"""
Parse SEACON calls.tsv and evaluate against ground truth.
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))

from haplotreesim import SimulationConfig, HaploTreeSimulator
from haplotreesim.metrics import (
    compute_all_hscn_metrics,
    compute_breakpoint_metrics,
    compute_clone_assignment_metrics
)

# ── 1. Re-run simulator with same seed to get ground truth ────────────────────
print("Regenerating ground truth...")
config = SimulationConfig(
    chromosome='chr1', bin_width=500000, num_clones=5, num_cells=100,
    lambda_events=1.5, lambda_amplitude=1.0, prob_wgd=0.0, random_seed=42
)
sim = HaploTreeSimulator(config)
sim.run()
ground_truth = sim.get_ground_truth()
segments = ground_truth['segments']
n_cells = 100
n_segs = len(segments)
bin_size = config.bin_width

# Build ground truth CN matrices (cells x segments)
clone_labels = np.array(ground_truth['clone_assignments'])
cn_A_true = np.zeros((n_cells, n_segs), dtype=int)
cn_B_true = np.zeros((n_cells, n_segs), dtype=int)
for cell_idx, clone_idx in enumerate(clone_labels):
    clone = sim.clones[clone_idx]
    for seg_idx, seg in enumerate(segments):
        cn_A_true[cell_idx, seg_idx] = int(np.mean(clone.cn_profile_A[seg.start_bin:seg.end_bin+1]))
        cn_B_true[cell_idx, seg_idx] = int(np.mean(clone.cn_profile_B[seg.start_bin:seg.end_bin+1]))

true_breakpoints = np.array([seg.start_bin for seg in segments[1:]])
print(f"  {n_cells} cells, {n_segs} segments, {len(true_breakpoints)} true breakpoints")

# ── 2. Parse SEACON calls.tsv ─────────────────────────────────────────────────
print("\nParsing SEACON output...")
calls = pd.read_csv('seacon_output/calls.tsv', sep='\t')
cell_names = [f'cell{i}' for i in range(n_cells)]

# Build per-cell bin-level CN from SEACON (bin index from start position)
cn_A_pred_bins = np.zeros((n_cells, 498), dtype=int)
cn_B_pred_bins = np.zeros((n_cells, 498), dtype=int)

for _, row in calls.iterrows():
    cell_idx = int(row['cell'].replace('cell', ''))
    bin_idx = int((row['start'] - 1) // bin_size)
    if bin_idx >= 498:
        continue
    a, b = str(row['CN']).split(',')
    cn_A_pred_bins[cell_idx, bin_idx] = int(a)
    cn_B_pred_bins[cell_idx, bin_idx] = int(b)

# Aggregate to segment level (mean per segment)
cn_A_pred = np.zeros((n_cells, n_segs), dtype=int)
cn_B_pred = np.zeros((n_cells, n_segs), dtype=int)
for seg_idx, seg in enumerate(segments):
    cn_A_pred[:, seg_idx] = np.round(
        np.mean(cn_A_pred_bins[:, seg.start_bin:seg.end_bin+1], axis=1)).astype(int)
    cn_B_pred[:, seg_idx] = np.round(
        np.mean(cn_B_pred_bins[:, seg.start_bin:seg.end_bin+1], axis=1)).astype(int)

print(f"  Parsed {len(calls)} rows from calls.tsv")

# ── 3. Extract predicted breakpoints from SEACON ─────────────────────────────
pred_breakpoints = set()
for cell_name, cell_calls in calls.groupby('cell'):
    cell_calls = cell_calls.sort_values('start')
    prev_cn = None
    for _, row in cell_calls.iterrows():
        cn = row['CN']
        bin_idx = int((row['start'] - 1) // bin_size)
        if prev_cn is not None and cn != prev_cn:
            pred_breakpoints.add(bin_idx)
        prev_cn = cn
pred_breakpoints = np.array(sorted(pred_breakpoints))
print(f"  Detected {len(pred_breakpoints)} predicted breakpoints")

# ── 4. Clone assignment from SEACON (majority CN state per cell) ──────────────
# Use total CN profile per cell as a simple fingerprint for clustering
tcn_profiles = cn_A_pred + cn_B_pred
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
km = KMeans(n_clusters=5, random_state=42, n_init=10)
pred_labels = km.fit_predict(tcn_profiles)

# ── 5. Compute metrics ────────────────────────────────────────────────────────
print("\nComputing metrics...")

hscn = compute_all_hscn_metrics(cn_A_true, cn_B_true, cn_A_pred, cn_B_pred, clone_labels)
bp   = compute_breakpoint_metrics(true_breakpoints, pred_breakpoints, tolerance=2)
clone_metrics = compute_clone_assignment_metrics(clone_labels, pred_labels)

print("\n" + "="*50)
print("SEACON EVALUATION RESULTS")
print("="*50)
print(f"  HSCN Error:       {hscn['hscn_error']:.4f}  (0=perfect)")
print(f"  LOH F1:           {hscn['loh_f1']:.4f}  (1=perfect)")
print(f"  MSR:              {hscn['msr']:.4f}  (1=perfect)")
print(f"  Breakpoint F1:    {bp['f1']:.4f}  (P={bp['precision']:.3f} R={bp['recall']:.3f})")
print(f"  Clone ARI:        {clone_metrics['ari']:.4f}  (1=perfect)")
print(f"  Clone NMI:        {clone_metrics['nmi']:.4f}  (1=perfect)")
print("="*50)
print("\nTrue breakpoints: ", true_breakpoints)
print("Pred breakpoints: ", pred_breakpoints[:20], "..." if len(pred_breakpoints) > 20 else "")
print("\nTrue CN sample (cell 0):", list(zip(cn_A_true[0], cn_B_true[0])))
print("Pred CN sample (cell 0):", list(zip(cn_A_pred[0], cn_B_pred[0])))
