"""
Haplotype-Specific Copy Number (HSCN) Metrics

Implements metrics from Section 4.3 of the paper:
- Swap-invariant HSCN error (Equation 57)
- Total copy number MSE (Equation 58)
- LOH detection precision/recall/F1
- Mirrored-subclone resolution (Equation 60)
"""

import numpy as np
from typing import Tuple, Dict, Optional


def compute_hscn_error(
    true_cn_A: np.ndarray,
    true_cn_B: np.ndarray,
    pred_cn_A: np.ndarray,
    pred_cn_B: np.ndarray,
    mask: Optional[np.ndarray] = None
) -> float:
    """
    Compute swap-invariant haplotype-specific copy number error (Equation 57).
    
    The error is the minimum of:
    1. Direct error: |true_A - pred_A| + |true_B - pred_B|
    2. Swapped error: |true_A - pred_B| + |true_B - pred_A|
    
    Args:
        true_cn_A: True copy number on haplotype A, shape (n_cells, n_segments)
        true_cn_B: True copy number on haplotype B, shape (n_cells, n_segments)
        pred_cn_A: Predicted copy number on haplotype A, shape (n_cells, n_segments)
        pred_cn_B: Predicted copy number on haplotype B, shape (n_cells, n_segments)
        mask: Optional boolean mask for valid positions
    
    Returns:
        Mean swap-invariant error across all cell-segment pairs
    """
    assert true_cn_A.shape == true_cn_B.shape == pred_cn_A.shape == pred_cn_B.shape
    
    # Compute direct error (no swap)
    error_direct = np.abs(true_cn_A - pred_cn_A) + np.abs(true_cn_B - pred_cn_B)
    
    # Compute swapped error
    error_swapped = np.abs(true_cn_A - pred_cn_B) + np.abs(true_cn_B - pred_cn_A)
    
    # Take minimum (swap-invariant)
    error = np.minimum(error_direct, error_swapped)
    
    # Apply mask if provided
    if mask is not None:
        error = error[mask]
    
    return float(np.mean(error))


def compute_tcn_mse(
    true_tcn: np.ndarray,
    pred_tcn: np.ndarray,
    mask: Optional[np.ndarray] = None
) -> float:
    """
    Compute total copy number (TCN) mean squared error (Equation 58).
    
    Args:
        true_tcn: True total copy number, shape (n_cells, n_segments)
        pred_tcn: Predicted total copy number, shape (n_cells, n_segments)
        mask: Optional boolean mask for valid positions
    
    Returns:
        Mean squared error
    """
    assert true_tcn.shape == pred_tcn.shape
    
    squared_error = (true_tcn - pred_tcn) ** 2
    
    if mask is not None:
        squared_error = squared_error[mask]
    
    return float(np.mean(squared_error))


def compute_loh_metrics(
    true_cn_A: np.ndarray,
    true_cn_B: np.ndarray,
    pred_cn_A: np.ndarray,
    pred_cn_B: np.ndarray,
    mask: Optional[np.ndarray] = None
) -> Dict[str, float]:
    """
    Compute LOH (Loss of Heterozygosity) detection metrics.
    
    LOH occurs when one haplotype has 0 copies but total CN > 0.
    
    Returns:
        Dictionary with keys: 'precision', 'recall', 'f1', 'tp', 'fp', 'fn'
    """
    assert true_cn_A.shape == true_cn_B.shape == pred_cn_A.shape == pred_cn_B.shape
    
    # Define LOH: one haplotype is 0, total CN > 0
    true_loh = ((true_cn_A == 0) | (true_cn_B == 0)) & ((true_cn_A + true_cn_B) > 0)
    pred_loh = ((pred_cn_A == 0) | (pred_cn_B == 0)) & ((pred_cn_A + pred_cn_B) > 0)
    
    if mask is not None:
        true_loh = true_loh[mask]
        pred_loh = pred_loh[mask]
    
    true_loh = true_loh.flatten()
    pred_loh = pred_loh.flatten()
    
    # Compute TP, FP, FN
    true_positive = np.sum(true_loh & pred_loh)
    false_positive = np.sum(~true_loh & pred_loh)
    false_negative = np.sum(true_loh & ~pred_loh)
    
    precision = true_positive / (true_positive + false_positive) if (true_positive + false_positive) > 0 else 0.0
    recall = true_positive / (true_positive + false_negative) if (true_positive + false_negative) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    
    return {
        'precision': float(precision),
        'recall': float(recall),
        'f1': float(f1),
        'tp': int(true_positive),
        'fp': int(false_positive),
        'fn': int(false_negative)
    }


def compute_mirrored_subclone_resolution(
    true_cn_A: np.ndarray,
    true_cn_B: np.ndarray,
    pred_cn_A: np.ndarray,
    pred_cn_B: np.ndarray,
    true_clone_labels: np.ndarray,
    mask: Optional[np.ndarray] = None
) -> float:
    """
    Compute mirrored-subclone resolution rate (Equation 60).
    
    A mirrored subclone pair has identical total CN but different haplotype states.
    Example: Clone 1 has (2,0) and Clone 2 has (0,2) in a segment.
    
    The key insight: if predictions are (2,0) vs (0,2), they ARE distinguished,
    even though one could be "swapped" - they're still different cell profiles.
    
    Returns:
        MSR rate: fraction of mirrored clone pairs correctly distinguished
    """
    assert true_cn_A.shape == true_cn_B.shape == pred_cn_A.shape == pred_cn_B.shape
    assert true_clone_labels.shape[0] == true_cn_A.shape[0]
    
    # Compute total CN
    true_tcn = true_cn_A + true_cn_B
    
    unique_clones = np.unique(true_clone_labels)
    
    if len(unique_clones) < 2:
        return 1.0
    
    # Find mirrored pairs
    mirrored_pairs = []
    
    for i, k1 in enumerate(unique_clones):
        for k2 in unique_clones[i+1:]:
            cells_k1 = np.where(true_clone_labels == k1)[0]
            cells_k2 = np.where(true_clone_labels == k2)[0]
            
            if len(cells_k1) == 0 or len(cells_k2) == 0:
                continue
            
            cell1 = cells_k1[0]
            cell2 = cells_k2[0]
            
            # Check if total CN is identical
            tcn_identical = np.all(true_tcn[cell1] == true_tcn[cell2])
            
            if not tcn_identical:
                continue
            
            # Check if haplotype states differ
            haplotype_differs = np.any(
                (true_cn_A[cell1] != true_cn_A[cell2]) | 
                (true_cn_B[cell1] != true_cn_B[cell2])
            )
            
            if haplotype_differs:
                mirrored_pairs.append((k1, k2, cell1, cell2))
    
    if len(mirrored_pairs) == 0:
        return 1.0
    
    # Check if each mirrored pair is correctly distinguished
    correctly_resolved = 0
    
    for k1, k2, cell1, cell2 in mirrored_pairs:
        # Check if predictions differ between the two cells
        # A method resolves the mirror if it predicts different CN profiles
        
        pred_diff = (
            np.any(pred_cn_A[cell1] != pred_cn_A[cell2]) or
            np.any(pred_cn_B[cell1] != pred_cn_B[cell2])
        )
        
        # Resolved if predictions differ
        if pred_diff:
            correctly_resolved += 1
    
    return correctly_resolved / len(mirrored_pairs)


def compute_all_hscn_metrics(
    true_cn_A: np.ndarray,
    true_cn_B: np.ndarray,
    pred_cn_A: np.ndarray,
    pred_cn_B: np.ndarray,
    true_clone_labels: Optional[np.ndarray] = None,
    mask: Optional[np.ndarray] = None
) -> Dict[str, float]:
    """
    Compute all HSCN-related metrics at once.
    
    Returns:
        Dictionary with all metric values
    """
    results = {}
    
    results['hscn_error'] = compute_hscn_error(true_cn_A, true_cn_B, pred_cn_A, pred_cn_B, mask)
    
    true_tcn = true_cn_A + true_cn_B
    pred_tcn = pred_cn_A + pred_cn_B
    results['tcn_mse'] = compute_tcn_mse(true_tcn, pred_tcn, mask)
    
    loh_metrics = compute_loh_metrics(true_cn_A, true_cn_B, pred_cn_A, pred_cn_B, mask)
    results.update({f'loh_{k}': v for k, v in loh_metrics.items()})
    
    if true_clone_labels is not None:
        results['msr'] = compute_mirrored_subclone_resolution(
            true_cn_A, true_cn_B, pred_cn_A, pred_cn_B, true_clone_labels, mask
        )
    
    return results
