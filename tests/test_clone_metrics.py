"""
Tests for clone assignment metrics (Week 12)
"""

import numpy as np
import sys
sys.path.insert(0, 'src')

from haplotreesim.metrics.clone_metrics import (
    compute_contingency_matrix,
    hungarian_matching,
    compute_ari,
    compute_nmi,
    handle_unequal_K,
    compute_clone_assignment_metrics
)

def test_contingency_matrix():
    """Test contingency matrix computation"""
    true = np.array([0, 0, 1, 1, 2])
    pred = np.array([0, 0, 1, 2, 2])
    
    C = compute_contingency_matrix(true, pred)
    
    expected = np.array([
        [2, 0, 0],  # True cluster 0: 2 in pred 0
        [0, 1, 1],  # True cluster 1: 1 in pred 1, 1 in pred 2
        [0, 0, 1]   # True cluster 2: 1 in pred 2
    ])
    
    assert np.array_equal(C, expected), f"Expected:\n{expected}\nGot:\n{C}"
    print("✓ Test 1: Contingency matrix")

def test_hungarian_perfect():
    """Test Hungarian matching with perfect clustering"""
    true = np.array([0, 0, 1, 1, 2, 2])
    pred = np.array([5, 5, 3, 3, 7, 7])  # Different labels, same clustering
    
    mapping, accuracy = hungarian_matching(true, pred)
    
    assert accuracy == 1.0, f"Expected accuracy 1.0, got {accuracy}"
    assert len(mapping) == 3, f"Expected 3 mappings, got {len(mapping)}"
    print("✓ Test 2: Hungarian perfect match")

def test_hungarian_imperfect():
    """Test Hungarian with some errors"""
    true = np.array([0, 0, 0, 1, 1, 1])
    pred = np.array([0, 0, 1, 1, 1, 1])  # One cell misassigned
    
    mapping, accuracy = hungarian_matching(true, pred)
    
    expected_accuracy = 5.0 / 6.0  # 5 correct out of 6
    assert abs(accuracy - expected_accuracy) < 1e-6, \
        f"Expected accuracy {expected_accuracy}, got {accuracy}"
    print("✓ Test 3: Hungarian with errors")

def test_hungarian_unequal_K_under():
    """Test Hungarian when K_pred < K_true (undersegmentation)"""
    true = np.array([0, 0, 1, 1, 2, 2])
    pred = np.array([0, 0, 0, 0, 1, 1])  # Only 2 clusters (merged 0,1)
    
    mapping, accuracy = hungarian_matching(true, pred)
    
    assert len(mapping) == 2, f"Expected 2 mappings, got {len(mapping)}"
    assert accuracy < 1.0, "Accuracy should be < 1.0 with merged clusters"
    print("✓ Test 4: Hungarian undersegmentation")

def test_hungarian_unequal_K_over():
    """Test Hungarian when K_pred > K_true (oversegmentation)"""
    true = np.array([0, 0, 0, 0, 1, 1])
    pred = np.array([0, 0, 1, 1, 2, 2])  # Split cluster 0 into 0,1
    
    mapping, accuracy = hungarian_matching(true, pred)
    
    assert len(mapping) <= 3, f"Expected ≤3 mappings, got {len(mapping)}"
    assert accuracy < 1.0, "Accuracy should be < 1.0 with split clusters"
    print("✓ Test 5: Hungarian oversegmentation")

def test_ari_perfect():
    """Test ARI with perfect clustering"""
    true = np.array([0, 0, 1, 1, 2, 2])
    pred = np.array([5, 5, 3, 3, 7, 7])
    
    ari = compute_ari(true, pred)
    assert ari == 1.0, f"Expected ARI 1.0, got {ari}"
    print("✓ Test 6: ARI perfect")

def test_ari_random():
    """Test ARI with random clustering"""
    true = np.array([0, 0, 0, 1, 1, 1])
    pred = np.array([0, 0, 1, 0, 1, 1])  # Somewhat random
    
    ari = compute_ari(true, pred)
    # ARI around 0 for random
    assert -0.5 < ari < 0.5, f"Expected ARI near 0, got {ari}"
    print("✓ Test 7: ARI random")

def test_nmi_perfect():
    """Test NMI with perfect clustering"""
    true = np.array([0, 0, 1, 1, 2, 2])
    pred = np.array([5, 5, 3, 3, 7, 7])
    
    nmi = compute_nmi(true, pred)
    assert nmi == 1.0, f"Expected NMI 1.0, got {nmi}"
    print("✓ Test 8: NMI perfect")

def test_handle_unequal_K_merged():
    """Test merge detection"""
    true = np.array([0, 0, 1, 1, 2, 2])
    pred = np.array([0, 0, 0, 0, 1, 1])  # Clusters 0,1 merged
    
    info = handle_unequal_K(true, pred)
    
    assert info['K_true'] == 3
    assert info['K_pred'] == 2
    assert len(info['merged_clones']) > 0, "Should detect merged clusters"
    print("✓ Test 9: Detect merged clusters")

def test_handle_unequal_K_split():
    """Test split detection"""
    true = np.array([0, 0, 0, 0, 1, 1])
    pred = np.array([0, 0, 1, 1, 2, 2])  # Cluster 0 split
    
    info = handle_unequal_K(true, pred)
    
    assert info['K_true'] == 2
    assert info['K_pred'] == 3
    assert len(info['split_clones']) > 0, "Should detect split clusters"
    print("✓ Test 10: Detect split clusters")

def test_all_metrics():
    """Test computing all metrics at once"""
    true = np.array([0, 0, 1, 1, 2, 2])
    pred = np.array([5, 5, 3, 3, 7, 7])
    
    metrics = compute_clone_assignment_metrics(true, pred)
    
    assert 'ari' in metrics
    assert 'nmi' in metrics
    assert 'hungarian_accuracy' in metrics
    assert 'K_true' in metrics
    assert 'K_pred' in metrics
    assert 'mapping' in metrics
    
    assert metrics['ari'] == 1.0
    assert metrics['nmi'] == 1.0
    assert metrics['hungarian_accuracy'] == 1.0
    print("✓ Test 11: All metrics computed")

if __name__ == '__main__':
    print("="*70)
    print("CLONE ASSIGNMENT METRICS TESTS (WEEK 12)")
    print("="*70)
    
    test_contingency_matrix()
    test_hungarian_perfect()
    test_hungarian_imperfect()
    test_hungarian_unequal_K_under()
    test_hungarian_unequal_K_over()
    test_ari_perfect()
    test_ari_random()
    test_nmi_perfect()
    test_handle_unequal_K_merged()
    test_handle_unequal_K_split()
    test_all_metrics()
    
    print("="*70)
    print("✓✓✓ ALL TESTS PASSED! ✓✓✓")
    print("="*70)
