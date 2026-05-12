"""
Breakpoint Detection Metrics

Measures accuracy of detecting segment boundaries (CNA breakpoints).
"""

import numpy as np
from typing import List, Tuple, Dict


def compute_breakpoint_metrics(
    true_breakpoints: np.ndarray,
    pred_breakpoints: np.ndarray,
    tolerance: int = 1,
    chromosome_length: int = None
) -> Dict[str, float]:
    """
    Compute breakpoint detection precision, recall, and F1.
    
    A predicted breakpoint matches a true breakpoint if they are within
    'tolerance' bins of each other.
    
    Args:
        true_breakpoints: Array of true breakpoint positions (bin indices)
        pred_breakpoints: Array of predicted breakpoint positions (bin indices)
        tolerance: Maximum distance (in bins) for a match
        chromosome_length: Total chromosome length in bins (for FP normalization)
    
    Returns:
        Dictionary with 'precision', 'recall', 'f1', 'tp', 'fp', 'fn'
    
    Example:
        >>> true_bp = np.array([100, 200, 300])
        >>> pred_bp = np.array([101, 199, 450])  # First two match, third is FP
        >>> metrics = compute_breakpoint_metrics(true_bp, pred_bp, tolerance=2)
        >>> print(f"Precision: {metrics['precision']:.2f}")
        >>> print(f"Recall: {metrics['recall']:.2f}")
    """
    true_breakpoints = np.sort(true_breakpoints)
    pred_breakpoints = np.sort(pred_breakpoints)
    
    if len(true_breakpoints) == 0 and len(pred_breakpoints) == 0:
        return {
            'precision': 1.0,
            'recall': 1.0,
            'f1': 1.0,
            'tp': 0,
            'fp': 0,
            'fn': 0
        }
    
    # Match predicted to true breakpoints
    matched_true = set()
    matched_pred = set()
    
    for i, pred_bp in enumerate(pred_breakpoints):
        # Find closest true breakpoint
        if len(true_breakpoints) > 0:
            distances = np.abs(true_breakpoints - pred_bp)
            closest_idx = np.argmin(distances)
            
            if distances[closest_idx] <= tolerance and closest_idx not in matched_true:
                matched_true.add(closest_idx)
                matched_pred.add(i)
    
    # Compute TP, FP, FN
    tp = len(matched_pred)
    fp = len(pred_breakpoints) - tp
    fn = len(true_breakpoints) - len(matched_true)
    
    # Compute precision, recall, F1
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    
    return {
        'precision': float(precision),
        'recall': float(recall),
        'f1': float(f1),
        'tp': int(tp),
        'fp': int(fp),
        'fn': int(fn)
    }
