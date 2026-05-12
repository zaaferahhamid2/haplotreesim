"""
Clone Assignment Metrics (Week 12)

Will implement:
- Hungarian matching for predicted-to-true clone mapping
- Adjusted Rand Index (ARI)
- Normalized Mutual Information (NMI)
"""

import numpy as np
from typing import Dict


def compute_clone_assignment_metrics(
    true_labels: np.ndarray,
    pred_labels: np.ndarray
) -> Dict[str, float]:
    """
    Compute clone assignment accuracy metrics.
    
    TODO: Week 12 implementation
    - Hungarian matching
    - ARI, NMI
    - Handle K_hat != K
    
    Args:
        true_labels: True clone assignments, shape (n_cells,)
        pred_labels: Predicted clone assignments, shape (n_cells,)
    
    Returns:
        Dictionary with 'ari', 'nmi', and other metrics
    """
    raise NotImplementedError("Week 12: Clone mapping not yet implemented")
