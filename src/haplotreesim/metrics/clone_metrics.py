"""
Clone Assignment Metrics (Week 12)

Implements:
- Hungarian algorithm for clone matching when K̂ ≠ K
- Adjusted Rand Index (ARI)
- Normalized Mutual Information (NMI)
- Merge/split reconciliation
"""

import numpy as np
from typing import Dict, Tuple, Optional
from scipy.optimize import linear_sum_assignment
from scipy.special import comb
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def compute_contingency_matrix(
    true_labels: np.ndarray,
    pred_labels: np.ndarray
) -> np.ndarray:
    """
    Compute contingency matrix between true and predicted cluster assignments.
    
    Entry [i,j] = number of cells assigned to true cluster i AND predicted cluster j
    
    Args:
        true_labels: True cluster assignments, shape (n_cells,)
        pred_labels: Predicted cluster assignments, shape (n_cells,)
    
    Returns:
        Contingency matrix, shape (K_true, K_pred)
    
    Example:
        >>> true = np.array([0, 0, 1, 1, 2])
        >>> pred = np.array([0, 0, 1, 2, 2])
        >>> C = compute_contingency_matrix(true, pred)
        >>> print(C)
        [[2 0 0]
         [0 1 1]
         [0 0 1]]
    """
    true_unique = np.unique(true_labels)
    pred_unique = np.unique(pred_labels)
    
    K_true = len(true_unique)
    K_pred = len(pred_unique)
    
    # Map labels to indices
    true_to_idx = {label: i for i, label in enumerate(true_unique)}
    pred_to_idx = {label: i for i, label in enumerate(pred_unique)}
    
    contingency = np.zeros((K_true, K_pred), dtype=int)
    
    for true_label, pred_label in zip(true_labels, pred_labels):
        i = true_to_idx[true_label]
        j = pred_to_idx[pred_label]
        contingency[i, j] += 1
    
    return contingency


def hungarian_matching(
    true_labels: np.ndarray,
    pred_labels: np.ndarray
) -> Tuple[Dict[int, int], float]:
    """
    Find optimal matching between predicted and true clusters using Hungarian algorithm.
    
    Maximizes overlap (number of cells correctly assigned after matching).
    Handles K̂ ≠ K by allowing unmatched clusters.
    
    Args:
        true_labels: True cluster assignments, shape (n_cells,)
        pred_labels: Predicted cluster assignments, shape (n_cells,)
    
    Returns:
        - mapping: Dict mapping pred_label -> true_label
        - accuracy: Fraction of cells correctly assigned after matching
    
    Example:
        >>> true = np.array([0, 0, 1, 1, 2, 2])
        >>> pred = np.array([5, 5, 3, 3, 7, 7])  # Different labels, same clustering
        >>> mapping, acc = hungarian_matching(true, pred)
        >>> print(f"Accuracy: {acc:.2f}")  # Should be 1.0
    """
    contingency = compute_contingency_matrix(true_labels, pred_labels)
    
    # Hungarian algorithm minimizes cost, so negate to maximize overlap
    # Pad to square matrix if needed
    K_true, K_pred = contingency.shape
    K_max = max(K_true, K_pred)
    
    cost_matrix = np.zeros((K_max, K_max))
    cost_matrix[:K_true, :K_pred] = -contingency  # Negate for maximization
    
    # Solve assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    # Build mapping (only for valid assignments)
    true_unique = sorted(np.unique(true_labels))
    pred_unique = sorted(np.unique(pred_labels))
    
    mapping = {}
    total_overlap = 0
    
    for i, j in zip(row_ind, col_ind):
        if i < K_true and j < K_pred:
            pred_label = pred_unique[j]
            true_label = true_unique[i]
            mapping[pred_label] = true_label
            total_overlap += contingency[i, j]
    
    accuracy = total_overlap / len(true_labels)
    
    return mapping, accuracy


def compute_ari(
    true_labels: np.ndarray,
    pred_labels: np.ndarray
) -> float:
    """
    Compute Adjusted Rand Index (ARI).
    
    ARI measures clustering similarity adjusted for chance.
    Range: [-1, 1], where 1 = perfect match, 0 = random, <0 = worse than random
    
    Uses sklearn implementation which is efficient and well-tested.
    
    Args:
        true_labels: True cluster assignments
        pred_labels: Predicted cluster assignments
    
    Returns:
        ARI score
    
    Example:
        >>> true = np.array([0, 0, 1, 1, 2, 2])
        >>> pred = np.array([0, 0, 1, 1, 2, 2])
        >>> ari = compute_ari(true, pred)
        >>> print(f"ARI: {ari:.3f}")  # Should be 1.0
    """
    return adjusted_rand_score(true_labels, pred_labels)


def compute_nmi(
    true_labels: np.ndarray,
    pred_labels: np.ndarray,
    average_method: str = 'arithmetic'
) -> float:
    """
    Compute Normalized Mutual Information (NMI).
    
    NMI measures clustering similarity based on mutual information.
    Range: [0, 1], where 1 = perfect match, 0 = independent
    
    Args:
        true_labels: True cluster assignments
        pred_labels: Predicted cluster assignments
        average_method: How to average MI ('arithmetic', 'geometric', 'min', 'max')
    
    Returns:
        NMI score
    
    Example:
        >>> true = np.array([0, 0, 1, 1, 2, 2])
        >>> pred = np.array([0, 0, 1, 1, 2, 2])
        >>> nmi = compute_nmi(true, pred)
        >>> print(f"NMI: {nmi:.3f}")  # Should be 1.0
    """
    return normalized_mutual_info_score(
        true_labels, 
        pred_labels,
        average_method=average_method
    )


def handle_unequal_K(
    true_labels: np.ndarray,
    pred_labels: np.ndarray,
    contingency: Optional[np.ndarray] = None
) -> Dict[str, any]:
    """
    Analyze merge/split behavior when K̂ ≠ K.
    
    When K̂ < K (undersegmentation): Some true clusters are merged
    When K̂ > K (oversegmentation): Some true clusters are split
    
    Args:
        true_labels: True cluster assignments
        pred_labels: Predicted cluster assignments
        contingency: Pre-computed contingency matrix (optional)
    
    Returns:
        Dictionary with:
        - 'K_true': Number of true clusters
        - 'K_pred': Number of predicted clusters
        - 'merged_clones': List of true clusters merged together
        - 'split_clones': List of true clusters split into multiple predicted
        - 'pure_matches': List of (true, pred) pairs with 1-to-1 mapping
    
    Example:
        >>> # True: 3 clusters, Pred: 2 clusters (merged)
        >>> true = np.array([0, 0, 1, 1, 2, 2])
        >>> pred = np.array([0, 0, 0, 0, 1, 1])  # Clusters 0,1 merged
        >>> info = handle_unequal_K(true, pred)
        >>> print(f"K_true: {info['K_true']}, K_pred: {info['K_pred']}")
    """
    if contingency is None:
        contingency = compute_contingency_matrix(true_labels, pred_labels)
    
    K_true, K_pred = contingency.shape
    
    true_unique = sorted(np.unique(true_labels))
    pred_unique = sorted(np.unique(pred_labels))
    
    # Analyze each true cluster
    merged_clones = []
    split_clones = []
    pure_matches = []
    
    for i, true_label in enumerate(true_unique):
        # How many predicted clusters does this true cluster map to?
        pred_matches = np.where(contingency[i, :] > 0)[0]
        
        if len(pred_matches) == 1:
            # This true cluster maps to exactly one predicted cluster
            j = pred_matches[0]
            # Check if it's a 1-to-1 match
            true_matches_for_pred = np.where(contingency[:, j] > 0)[0]
            if len(true_matches_for_pred) == 1:
                pure_matches.append((true_label, pred_unique[j]))
            # If multiple true clusters map to this pred cluster, it's merged
            elif len(true_matches_for_pred) > 1:
                merged_true = [true_unique[k] for k in true_matches_for_pred]
                if merged_true not in merged_clones:
                    merged_clones.append(merged_true)
        
        elif len(pred_matches) > 1:
            # This true cluster is split across multiple predicted clusters
            split_clones.append({
                'true_label': true_label,
                'pred_labels': [pred_unique[j] for j in pred_matches],
                'cell_counts': [contingency[i, j] for j in pred_matches]
            })
    
    return {
        'K_true': K_true,
        'K_pred': K_pred,
        'merged_clones': merged_clones,
        'split_clones': split_clones,
        'pure_matches': pure_matches
    }


def compute_clone_assignment_metrics(
    true_labels: np.ndarray,
    pred_labels: np.ndarray
) -> Dict[str, float]:
    """
    Compute all clone assignment metrics.
    
    Args:
        true_labels: True clone assignments, shape (n_cells,)
        pred_labels: Predicted clone assignments, shape (n_cells,)
    
    Returns:
        Dictionary with:
        - 'ari': Adjusted Rand Index
        - 'nmi': Normalized Mutual Information
        - 'hungarian_accuracy': Accuracy after optimal matching
        - 'K_true': Number of true clusters
        - 'K_pred': Number of predicted clusters
        - 'mapping': Dict mapping pred -> true labels (via Hungarian)
    
    Example:
        >>> true = np.array([0, 0, 1, 1, 2, 2])
        >>> pred = np.array([5, 5, 3, 3, 7, 7])
        >>> metrics = compute_clone_assignment_metrics(true, pred)
        >>> print(f"ARI: {metrics['ari']:.3f}")
        >>> print(f"NMI: {metrics['nmi']:.3f}")
    """
    # Basic counts
    K_true = len(np.unique(true_labels))
    K_pred = len(np.unique(pred_labels))
    
    # Compute metrics
    ari = compute_ari(true_labels, pred_labels)
    nmi = compute_nmi(true_labels, pred_labels)
    mapping, accuracy = hungarian_matching(true_labels, pred_labels)
    
    # Analyze K differences
    unequal_info = handle_unequal_K(true_labels, pred_labels)
    
    return {
        'ari': float(ari),
        'nmi': float(nmi),
        'hungarian_accuracy': float(accuracy),
        'K_true': int(K_true),
        'K_pred': int(K_pred),
        'mapping': mapping,
        'n_pure_matches': len(unequal_info['pure_matches']),
        'n_merged': len(unequal_info['merged_clones']),
        'n_split': len(unequal_info['split_clones'])
    }
