"""
Tests for HSCN metrics (Week 11)
"""

import numpy as np
import sys
sys.path.insert(0, 'src')

from haplotreesim.metrics.hscn_metrics import (
    compute_hscn_error,
    compute_tcn_mse,
    compute_loh_metrics,
    compute_mirrored_subclone_resolution,
    compute_all_hscn_metrics
)

def test_hscn_error_perfect():
    """Test HSCN error with perfect predictions"""
    true_A = np.array([[2, 1, 0]])
    true_B = np.array([[0, 1, 2]])
    pred_A = np.array([[2, 1, 0]])
    pred_B = np.array([[0, 1, 2]])
    
    error = compute_hscn_error(true_A, true_B, pred_A, pred_B)
    assert error == 0.0, f"Expected 0, got {error}"
    print("✓ Test 1: Perfect prediction (no swap)")

def test_hscn_error_swapped():
    """Test HSCN error with swapped haplotypes"""
    true_A = np.array([[2, 1, 0]])
    true_B = np.array([[0, 1, 2]])
    pred_A = np.array([[0, 1, 2]])  # Swapped!
    pred_B = np.array([[2, 1, 0]])
    
    error = compute_hscn_error(true_A, true_B, pred_A, pred_B)
    assert error == 0.0, f"Expected 0 (swap-invariant), got {error}"
    print("✓ Test 2: Perfect prediction (swapped)")

def test_hscn_error_off_by_one():
    """Test HSCN error with off-by-one errors"""
    true_A = np.array([[2, 1, 0]])
    true_B = np.array([[0, 1, 2]])
    pred_A = np.array([[3, 1, 0]])  # Off by 1 in first position
    pred_B = np.array([[0, 1, 2]])
    
    error = compute_hscn_error(true_A, true_B, pred_A, pred_B)
    expected = 1.0 / 3  # Error of 1 in 1 out of 3 positions
    assert abs(error - expected) < 1e-6, f"Expected {expected}, got {error}"
    print("✓ Test 3: Off-by-one error")

def test_tcn_mse_perfect():
    """Test TCN MSE with perfect predictions"""
    true_tcn = np.array([[2, 3, 4]])
    pred_tcn = np.array([[2, 3, 4]])
    
    mse = compute_tcn_mse(true_tcn, pred_tcn)
    assert mse == 0.0, f"Expected 0, got {mse}"
    print("✓ Test 4: TCN MSE perfect")

def test_tcn_mse_errors():
    """Test TCN MSE with errors"""
    true_tcn = np.array([[2, 3, 4]])
    pred_tcn = np.array([[2, 4, 3]])  # Errors of 0, 1, 1
    
    mse = compute_tcn_mse(true_tcn, pred_tcn)
    expected = (0**2 + 1**2 + 1**2) / 3
    assert abs(mse - expected) < 1e-6, f"Expected {expected}, got {mse}"
    print("✓ Test 5: TCN MSE with errors")

def test_loh_perfect():
    """Test LOH detection with perfect predictions"""
    true_A = np.array([[2, 0, 1, 0]])
    true_B = np.array([[0, 2, 1, 0]])  # LOH at positions 0, 1
    pred_A = np.array([[2, 0, 1, 0]])
    pred_B = np.array([[0, 2, 1, 0]])
    
    metrics = compute_loh_metrics(true_A, true_B, pred_A, pred_B)
    assert metrics['precision'] == 1.0, f"Expected precision 1.0, got {metrics['precision']}"
    assert metrics['recall'] == 1.0, f"Expected recall 1.0, got {metrics['recall']}"
    assert metrics['f1'] == 1.0, f"Expected F1 1.0, got {metrics['f1']}"
    print("✓ Test 6: LOH detection perfect")

def test_loh_partial():
    """Test LOH detection with partial correctness"""
    # True: (2,0)=LOH, (0,2)=LOH, (1,1)=not LOH
    true_A = np.array([[2, 0, 1]])
    true_B = np.array([[0, 2, 1]])
    # Pred: (2,0)=LOH, (1,1)=not LOH, (1,1)=not LOH
    # Position 1 should be LOH but predicted as (1,1)
    pred_A = np.array([[2, 1, 1]])
    pred_B = np.array([[0, 1, 1]])
    
    metrics = compute_loh_metrics(true_A, true_B, pred_A, pred_B)
    # TP=1 (position 0 correct), FN=1 (position 1 missed)
    assert metrics['tp'] == 1, f"Expected TP=1, got {metrics['tp']}"
    assert metrics['fn'] == 1, f"Expected FN=1, got {metrics['fn']}"
    assert metrics['precision'] == 1.0, f"Expected precision=1.0, got {metrics['precision']}"
    assert metrics['recall'] == 0.5, f"Expected recall=0.5, got {metrics['recall']}"
    print("✓ Test 7: LOH detection partial")

def test_msr_no_mirrored():
    """Test MSR when no mirrored clones exist"""
    true_A = np.array([[2, 1], [2, 1]])  # Both clones identical
    true_B = np.array([[0, 1], [0, 1]])
    pred_A = np.array([[2, 1], [2, 1]])
    pred_B = np.array([[0, 1], [0, 1]])
    labels = np.array([0, 1])
    
    msr = compute_mirrored_subclone_resolution(true_A, true_B, pred_A, pred_B, labels)
    assert msr == 1.0, f"Expected 1.0 (no mirrored pairs), got {msr}"
    print("✓ Test 8: MSR with no mirrored clones")

def test_msr_mirrored_resolved():
    """Test MSR with mirrored clones correctly resolved"""
    # Clone 0: (2,0), Clone 1: (0,2) - mirrored!
    true_A = np.array([[2], [0]])
    true_B = np.array([[0], [2]])
    pred_A = np.array([[2], [0]])  # Correctly distinguished
    pred_B = np.array([[0], [2]])
    labels = np.array([0, 1])
    
    msr = compute_mirrored_subclone_resolution(true_A, true_B, pred_A, pred_B, labels)
    assert msr == 1.0, f"Expected 1.0 (resolved), got {msr}"
    print("✓ Test 9: MSR with mirrored clones resolved")

def test_msr_mirrored_unresolved():
    """Test MSR with mirrored clones not resolved"""
    # Clone 0: (2,0), Clone 1: (0,2) - mirrored!
    true_A = np.array([[2], [0]])
    true_B = np.array([[0], [2]])
    pred_A = np.array([[2], [2]])  # Failed to distinguish
    pred_B = np.array([[0], [0]])
    labels = np.array([0, 1])
    
    msr = compute_mirrored_subclone_resolution(true_A, true_B, pred_A, pred_B, labels)
    assert msr == 0.0, f"Expected 0.0 (unresolved), got {msr}"
    print("✓ Test 10: MSR with mirrored clones unresolved")

def test_all_metrics():
    """Test computing all metrics at once"""
    true_A = np.array([[2, 0, 1]])
    true_B = np.array([[0, 2, 1]])
    pred_A = np.array([[2, 0, 1]])
    pred_B = np.array([[0, 2, 1]])
    labels = np.array([0])
    
    metrics = compute_all_hscn_metrics(true_A, true_B, pred_A, pred_B, labels)
    
    assert 'hscn_error' in metrics
    assert 'tcn_mse' in metrics
    assert 'loh_f1' in metrics
    assert 'msr' in metrics
    print("✓ Test 11: All metrics computed")

if __name__ == '__main__':
    print("="*70)
    print("HSCN METRICS TESTS")
    print("="*70)
    
    test_hscn_error_perfect()
    test_hscn_error_swapped()
    test_hscn_error_off_by_one()
    test_tcn_mse_perfect()
    test_tcn_mse_errors()
    test_loh_perfect()
    test_loh_partial()
    test_msr_no_mirrored()
    test_msr_mirrored_resolved()
    test_msr_mirrored_unresolved()
    test_all_metrics()
    
    print("="*70)
    print("✓✓✓ ALL TESTS PASSED! ✓✓✓")
    print("="*70)
