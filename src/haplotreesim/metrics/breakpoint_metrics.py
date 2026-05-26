"""
Breakpoint Detection Metrics

Measures accuracy of detecting segment boundaries (CNA breakpoints).
"""

import numpy as np
from typing import List, Tuple, Dict


def _validate_boundary_mask(boundary_mask, n_boundaries: int) -> np.ndarray:
    """Validate an optional mask over adjacent-bin boundaries."""
    if boundary_mask is None:
        return np.ones(n_boundaries, dtype=bool)

    mask = np.asarray(boundary_mask, dtype=bool)
    if mask.shape != (n_boundaries,):
        raise ValueError("boundary_mask must have length n_bins - 1")
    return mask


def compute_tcn_breakpoint_frequencies(pred_tcn: np.ndarray, boundary_mask=None) -> np.ndarray:
    """
    Compute recurrence frequency of predicted total-CN changes per bin boundary.

    Boundary coordinate b denotes the boundary between bins b-1 and b in
    0-based Python arrays, matching segment start-bin breakpoint coordinates.
    Returns an array of length n_bins - 1, where entry j corresponds to
    breakpoint coordinate j + 1.
    """
    pred_tcn = np.asarray(pred_tcn, dtype=float)

    if pred_tcn.ndim != 2:
        raise ValueError("pred_tcn must have shape (n_cells, n_bins)")

    if pred_tcn.shape[1] < 2:
        return np.array([], dtype=float)

    boundary_mask = _validate_boundary_mask(boundary_mask, pred_tcn.shape[1] - 1)
    valid = np.isfinite(pred_tcn[:, :-1]) & np.isfinite(pred_tcn[:, 1:])
    valid &= boundary_mask[None, :]

    changes = valid & (pred_tcn[:, :-1] != pred_tcn[:, 1:])

    numerator = changes.sum(axis=0)
    denominator = valid.sum(axis=0)

    frequencies = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=float),
        where=denominator > 0,
    )

    return frequencies

def compute_hscn_breakpoint_frequencies(
    pred_cn_A: np.ndarray,
    pred_cn_B: np.ndarray,
    boundary_mask=None,
) -> np.ndarray:
    """
    Compute recurrence frequency of predicted HSCN changes per bin boundary.

    A cell contributes a predicted breakpoint at boundary b when either
    haplotype-specific copy-number track changes across that boundary.

    Boundary coordinate b denotes the boundary between bins b-1 and b in
    0-based Python arrays, matching segment start-bin breakpoint coordinates.
    Returns an array of length n_bins - 1, where entry j corresponds to
    breakpoint coordinate j + 1.
    """
    pred_cn_A = np.asarray(pred_cn_A, dtype=float)
    pred_cn_B = np.asarray(pred_cn_B, dtype=float)

    if pred_cn_A.shape != pred_cn_B.shape:
        raise ValueError("pred_cn_A and pred_cn_B must have the same shape")

    if pred_cn_A.ndim != 2:
        raise ValueError("pred_cn_A and pred_cn_B must have shape (n_cells, n_bins)")

    if pred_cn_A.shape[1] < 2:
        return np.array([], dtype=float)

    boundary_mask = _validate_boundary_mask(boundary_mask, pred_cn_A.shape[1] - 1)
    valid = (
        np.isfinite(pred_cn_A[:, :-1])
        & np.isfinite(pred_cn_A[:, 1:])
        & np.isfinite(pred_cn_B[:, :-1])
        & np.isfinite(pred_cn_B[:, 1:])
    )
    valid &= boundary_mask[None, :]

    changes = valid & (
        (pred_cn_A[:, :-1] != pred_cn_A[:, 1:])
        | (pred_cn_B[:, :-1] != pred_cn_B[:, 1:])
    )

    numerator = changes.sum(axis=0)
    denominator = valid.sum(axis=0)

    frequencies = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=float),
        where=denominator > 0,
    )

    return frequencies

def breakpoints_from_frequencies(frequencies: np.ndarray, threshold: float = 0.05) -> np.ndarray:
    """
    Threshold boundary recurrence frequencies into breakpoint start-bin indices.
    """
    if threshold < 0 or threshold > 1:
        raise ValueError("threshold must be in [0, 1]")
    frequencies = np.asarray(frequencies, dtype=float)
    return np.flatnonzero(frequencies >= threshold) + 1


def predicted_breakpoints_from_tcn(
    pred_tcn: np.ndarray,
    threshold: float = 0.05,
    boundary_mask=None,
) -> np.ndarray:
    """Compute thresholded recurrent breakpoints from cell-level total CN profiles."""
    return breakpoints_from_frequencies(
        compute_tcn_breakpoint_frequencies(pred_tcn, boundary_mask=boundary_mask),
        threshold=threshold
    )


def predicted_breakpoints_from_hscn(
    pred_cn_A: np.ndarray,
    pred_cn_B: np.ndarray,
    threshold: float = 0.05,
    boundary_mask=None,
) -> np.ndarray:
    """Compute thresholded recurrent breakpoints from cell-level HSCN profiles."""
    return breakpoints_from_frequencies(
        compute_hscn_breakpoint_frequencies(pred_cn_A, pred_cn_B, boundary_mask=boundary_mask),
        threshold=threshold
    )


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
            'fn': 0,
            'n_true': 0,
            'n_pred': 0,
        }
    
    # Match predicted to true breakpoints one-to-one by minimum distance.
    matched_true = set()
    matched_pred = set()

    candidate_pairs = []
    for pred_idx, pred_bp in enumerate(pred_breakpoints):
        for true_idx, true_bp in enumerate(true_breakpoints):
            distance = abs(int(pred_bp) - int(true_bp))
            if distance <= tolerance:
                candidate_pairs.append((distance, pred_idx, true_idx))

    for _, pred_idx, true_idx in sorted(candidate_pairs):
        if pred_idx in matched_pred or true_idx in matched_true:
            continue
        matched_pred.add(pred_idx)
        matched_true.add(true_idx)
    
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
        'fn': int(fn),
        'n_true': int(len(true_breakpoints)),
        'n_pred': int(len(pred_breakpoints)),
    }


def compute_recurrent_breakpoint_metrics(
    true_breakpoints: np.ndarray,
    *,
    pred_tcn: np.ndarray = None,
    pred_cn_A: np.ndarray = None,
    pred_cn_B: np.ndarray = None,
    threshold: float = 0.05,
    tolerance: int = 2,
    boundary_mask=None,
) -> Dict[str, float]:
    """
    Compute breakpoint metrics from cell-level bin CN profiles using recurrence.

    Provide either ``pred_tcn`` for total-CN outputs or both ``pred_cn_A`` and
    ``pred_cn_B`` for HSCN outputs.
    """
    if pred_tcn is not None:
        frequencies = compute_tcn_breakpoint_frequencies(pred_tcn, boundary_mask=boundary_mask)
    elif pred_cn_A is not None and pred_cn_B is not None:
        frequencies = compute_hscn_breakpoint_frequencies(
            pred_cn_A,
            pred_cn_B,
            boundary_mask=boundary_mask,
        )
    else:
        raise ValueError("Provide pred_tcn or both pred_cn_A and pred_cn_B")

    pred_breakpoints = breakpoints_from_frequencies(frequencies, threshold=threshold)
    metrics = compute_breakpoint_metrics(
        true_breakpoints=true_breakpoints,
        pred_breakpoints=pred_breakpoints,
        tolerance=tolerance,
    )
    metrics.update({
        'threshold': float(threshold),
        'n_true': int(len(true_breakpoints)),
        'n_pred': int(len(pred_breakpoints)),
        'pred_breakpoints': pred_breakpoints,
        'breakpoint_frequencies': frequencies,
    })
    return metrics


def compute_breakpoint_sensitivity_curve(
    true_breakpoints: np.ndarray,
    *,
    pred_tcn: np.ndarray = None,
    pred_cn_A: np.ndarray = None,
    pred_cn_B: np.ndarray = None,
    thresholds=(0.01, 0.05, 0.10, 0.20),
    tolerance: int = 2,
    boundary_mask=None,
) -> Dict[float, Dict[str, float]]:
    """Compute recurrent breakpoint metrics over a threshold grid."""
    return {
        float(threshold): compute_recurrent_breakpoint_metrics(
            true_breakpoints=true_breakpoints,
            pred_tcn=pred_tcn,
            pred_cn_A=pred_cn_A,
            pred_cn_B=pred_cn_B,
            threshold=float(threshold),
            tolerance=tolerance,
            boundary_mask=boundary_mask,
        )
        for threshold in thresholds
    }
