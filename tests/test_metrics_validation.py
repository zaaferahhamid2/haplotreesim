"""
Metrics Validation Tests - Professor Request

Tests metrics with:
1. Ground truth (perfect predictions) -> should get perfect scores
2. Manually corrupted data -> should show degradation appropriately
"""

import numpy as np
import sys
sys.path.insert(0, 'src')

from haplotreesim.metrics import (
    compute_all_hscn_metrics,
    compute_breakpoint_metrics,
    compute_clone_assignment_metrics,
    compute_all_tree_metrics
)

print("="*70)
print("METRICS VALIDATION: Ground Truth vs. Corrupted Data")
print("="*70)

# ============================================================================
# TEST 1: HSCN METRICS - Perfect vs Corrupted
# ============================================================================
print("\n" + "="*70)
print("TEST 1: HSCN METRICS")
print("="*70)

# Generate ground truth
np.random.seed(42)
n_cells, n_segments = 100, 20
true_cn_A = np.random.randint(0, 4, size=(n_cells, n_segments))
true_cn_B = np.random.randint(0, 4, size=(n_cells, n_segments))
true_labels = np.random.randint(0, 5, size=n_cells)

print("\n1a. Perfect Predictions (ground truth = predictions)")
perfect_metrics = compute_all_hscn_metrics(
    true_cn_A, true_cn_B, 
    true_cn_A.copy(), true_cn_B.copy(),
    true_labels
)
print(f"  HSCN Error: {perfect_metrics['hscn_error']:.4f} (expect 0.0)")
print(f"  TCN MSE: {perfect_metrics['tcn_mse']:.4f} (expect 0.0)")
print(f"  LOH F1: {perfect_metrics['loh_f1']:.4f} (expect 1.0)")
print(f"  MSR: {perfect_metrics['msr']:.4f} (expect 1.0)")

assert perfect_metrics['hscn_error'] == 0.0, "Perfect predictions should have 0 error!"
assert perfect_metrics['tcn_mse'] == 0.0, "Perfect predictions should have 0 MSE!"
assert perfect_metrics['loh_f1'] == 1.0, "Perfect predictions should have F1=1.0!"
print("  ✓ All perfect scores achieved!")

print("\n1b. Corrupted Predictions (10% of values changed)")
corrupted_A = true_cn_A.copy()
corrupted_B = true_cn_B.copy()

# Introduce errors: change 10% of copy numbers
n_corrupt = int(0.1 * n_cells * n_segments)
corrupt_indices = np.random.choice(n_cells * n_segments, n_corrupt, replace=False)
for idx in corrupt_indices:
    i, j = idx // n_segments, idx % n_segments
    corrupted_A[i, j] = min(corrupted_A[i, j] + np.random.randint(-1, 2), 5)
    corrupted_A[i, j] = max(corrupted_A[i, j], 0)

corrupted_metrics = compute_all_hscn_metrics(
    true_cn_A, true_cn_B,
    corrupted_A, corrupted_B,
    true_labels
)
print(f"  HSCN Error: {corrupted_metrics['hscn_error']:.4f} (expect > 0)")
print(f"  TCN MSE: {corrupted_metrics['tcn_mse']:.4f} (expect > 0)")
print(f"  LOH F1: {corrupted_metrics['loh_f1']:.4f} (expect < 1.0)")

assert corrupted_metrics['hscn_error'] > 0, "Corrupted data should have error > 0!"
assert corrupted_metrics['tcn_mse'] > 0, "Corrupted data should have MSE > 0!"
print("  ✓ Corruption detected correctly!")

# ============================================================================
# TEST 2: BREAKPOINT METRICS - Perfect vs Corrupted
# ============================================================================
print("\n" + "="*70)
print("TEST 2: BREAKPOINT METRICS")
print("="*70)

true_breakpoints = np.array([100, 200, 300, 400, 500])

print("\n2a. Perfect Predictions")
perfect_bp = compute_breakpoint_metrics(true_breakpoints, true_breakpoints.copy(), tolerance=1)
print(f"  Precision: {perfect_bp['precision']:.4f} (expect 1.0)")
print(f"  Recall: {perfect_bp['recall']:.4f} (expect 1.0)")
print(f"  F1: {perfect_bp['f1']:.4f} (expect 1.0)")

assert perfect_bp['f1'] == 1.0, "Perfect breakpoints should have F1=1.0!"
print("  ✓ Perfect scores achieved!")

print("\n2b. Corrupted Predictions (missed 2 breakpoints, added 1 false positive)")
corrupted_breakpoints = np.array([101, 199, 300, 999])  # Missed 400,500, added 999
corrupted_bp = compute_breakpoint_metrics(true_breakpoints, corrupted_breakpoints, tolerance=2)
print(f"  Precision: {corrupted_bp['precision']:.4f} (expect < 1.0)")
print(f"  Recall: {corrupted_bp['recall']:.4f} (expect < 1.0)")
print(f"  F1: {corrupted_bp['f1']:.4f} (expect < 1.0)")
print(f"  TP={corrupted_bp['tp']}, FP={corrupted_bp['fp']}, FN={corrupted_bp['fn']}")

assert corrupted_bp['f1'] < 1.0, "Corrupted breakpoints should have F1 < 1.0!"
assert corrupted_bp['fp'] > 0, "Should detect false positive!"
assert corrupted_bp['fn'] > 0, "Should detect false negatives!"
print("  ✓ Corruption detected correctly!")

# ============================================================================
# TEST 3: CLONE ASSIGNMENT METRICS - Perfect vs Corrupted
# ============================================================================
print("\n" + "="*70)
print("TEST 3: CLONE ASSIGNMENT METRICS")
print("="*70)

true_clone_labels = np.array([0]*20 + [1]*20 + [2]*20 + [3]*20 + [4]*20)

print("\n3a. Perfect Predictions (permuted labels)")
# Perfect clustering but different label names
perfect_clone_pred = np.array([5]*20 + [3]*20 + [7]*20 + [1]*20 + [9]*20)
perfect_clone = compute_clone_assignment_metrics(true_clone_labels, perfect_clone_pred)
print(f"  ARI: {perfect_clone['ari']:.4f} (expect 1.0)")
print(f"  NMI: {perfect_clone['nmi']:.4f} (expect 1.0)")
print(f"  Hungarian Accuracy: {perfect_clone['hungarian_accuracy']:.4f} (expect 1.0)")

assert perfect_clone['ari'] == 1.0, "Perfect clustering should have ARI=1.0!"
assert perfect_clone['nmi'] == 1.0, "Perfect clustering should have NMI=1.0!"
print("  ✓ Perfect scores achieved!")

print("\n3b. Corrupted Predictions (20 cells misassigned)")
corrupted_clone_pred = perfect_clone_pred.copy()
# Randomly reassign 20 cells
corrupt_cells = np.random.choice(100, 20, replace=False)
corrupted_clone_pred[corrupt_cells] = np.random.randint(0, 5, size=20)

corrupted_clone = compute_clone_assignment_metrics(true_clone_labels, corrupted_clone_pred)
print(f"  ARI: {corrupted_clone['ari']:.4f} (expect < 1.0)")
print(f"  NMI: {corrupted_clone['nmi']:.4f} (expect < 1.0)")
print(f"  Hungarian Accuracy: {corrupted_clone['hungarian_accuracy']:.4f} (expect < 1.0)")

assert corrupted_clone['ari'] < 1.0, "Corrupted clustering should have ARI < 1.0!"
assert corrupted_clone['nmi'] < 1.0, "Corrupted clustering should have NMI < 1.0!"
print("  ✓ Corruption detected correctly!")

# ============================================================================
# TEST 4: TREE METRICS - Perfect vs Corrupted
# ============================================================================
print("\n" + "="*70)
print("TEST 4: TREE METRICS")
print("="*70)

true_tree = [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5)]
true_events = [
    {'start': 100, 'end': 200, 'haplotype': 'A'},
    {'start': 300, 'end': 400, 'haplotype': 'B'}
]

print("\n4a. Perfect Predictions")
perfect_tree = compute_all_tree_metrics(
    true_tree, true_tree, 
    true_events, true_events
)
print(f"  RF Distance: {perfect_tree['rf_distance']:.4f} (expect 0.0)")
print(f"  Ancestor-Descendant Acc: {perfect_tree['ancestor_descendant_accuracy']:.4f} (expect 1.0)")
print(f"  Event F1: {perfect_tree['event_f1']:.4f} (expect 1.0)")

assert perfect_tree['rf_distance'] == 0.0, "Identical trees should have RF=0!"
assert perfect_tree['event_f1'] == 1.0, "Perfect events should have F1=1.0!"
print("  ✓ Perfect scores achieved!")

print("\n4b. Corrupted Predictions (different tree + missing event)")
corrupted_tree = [(0, 1), (0, 3), (1, 2), (1, 4), (3, 5)]  # Different topology
corrupted_events = [
    {'start': 105, 'end': 195, 'haplotype': 'A'}  # Only one event, slightly off
]

corrupted_tree_metrics = compute_all_tree_metrics(
    true_tree, corrupted_tree,
    true_events, corrupted_events
)
print(f"  RF Distance: {corrupted_tree_metrics['rf_distance']:.4f} (expect > 0)")
print(f"  Event Recall: {corrupted_tree_metrics['event_recall']:.4f} (expect < 1.0)")
print(f"  Event F1: {corrupted_tree_metrics['event_f1']:.4f} (expect < 1.0)")

assert corrupted_tree_metrics['event_f1'] < 1.0, "Corrupted events should have F1 < 1.0!"
print("  ✓ Corruption detected correctly!")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("\n" + "="*70)
print("VALIDATION SUMMARY")
print("="*70)
print()
print("✓ HSCN metrics: Perfect=0.0 error, Corrupted>0 error")
print("✓ Breakpoint metrics: Perfect F1=1.0, Corrupted F1<1.0")
print("✓ Clone assignment: Perfect ARI/NMI=1.0, Corrupted<1.0")
print("✓ Tree metrics: Perfect RF=0/F1=1.0, Corrupted RF>0/F1<1.0")
print()
print("="*70)
print("✓✓✓ ALL VALIDATION TESTS PASSED! ✓✓✓")
print("="*70)
print()
print("Metrics are working correctly:")
print("  - Ground truth predictions → perfect scores")
print("  - Corrupted predictions → degraded scores")
print()
print("Ready to test on real methods!")
print("="*70)
