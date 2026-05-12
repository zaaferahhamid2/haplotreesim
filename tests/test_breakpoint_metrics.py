"""
Tests for breakpoint metrics
"""

import numpy as np
import sys
sys.path.insert(0, 'src')

from haplotreesim.metrics.breakpoint_metrics import compute_breakpoint_metrics

def test_perfect_match():
    """Test with perfect breakpoint detection"""
    true_bp = np.array([100, 200, 300])
    pred_bp = np.array([100, 200, 300])
    
    metrics = compute_breakpoint_metrics(true_bp, pred_bp, tolerance=1)
    assert metrics['precision'] == 1.0
    assert metrics['recall'] == 1.0
    assert metrics['f1'] == 1.0
    print("✓ Test 1: Perfect match")

def test_within_tolerance():
    """Test breakpoints within tolerance"""
    true_bp = np.array([100, 200, 300])
    pred_bp = np.array([101, 199, 301])  # All within ±2
    
    metrics = compute_breakpoint_metrics(true_bp, pred_bp, tolerance=2)
    assert metrics['precision'] == 1.0
    assert metrics['recall'] == 1.0
    print("✓ Test 2: Within tolerance")

def test_missed_breakpoint():
    """Test with missed breakpoint (FN)"""
    true_bp = np.array([100, 200, 300])
    pred_bp = np.array([100, 200])  # Missed 300
    
    metrics = compute_breakpoint_metrics(true_bp, pred_bp, tolerance=1)
    assert metrics['tp'] == 2
    assert metrics['fn'] == 1
    assert metrics['recall'] == 2/3
    print("✓ Test 3: Missed breakpoint")

def test_false_positive():
    """Test with false positive breakpoint"""
    true_bp = np.array([100, 200])
    pred_bp = np.array([100, 200, 999])  # Extra FP
    
    metrics = compute_breakpoint_metrics(true_bp, pred_bp, tolerance=1)
    assert metrics['tp'] == 2
    assert metrics['fp'] == 1
    assert metrics['precision'] == 2/3
    print("✓ Test 4: False positive")

def test_empty():
    """Test with no breakpoints"""
    true_bp = np.array([])
    pred_bp = np.array([])
    
    metrics = compute_breakpoint_metrics(true_bp, pred_bp, tolerance=1)
    assert metrics['precision'] == 1.0
    assert metrics['recall'] == 1.0
    print("✓ Test 5: No breakpoints")

if __name__ == '__main__':
    print("="*70)
    print("BREAKPOINT METRICS TESTS")
    print("="*70)
    
    test_perfect_match()
    test_within_tolerance()
    test_missed_breakpoint()
    test_false_positive()
    test_empty()
    
    print("="*70)
    print("✓✓✓ ALL TESTS PASSED! ✓✓✓")
    print("="*70)
