"""
Evaluation metrics for HaploTreeSim benchmark.

Week 11: HSCN metrics + breakpoint metrics
Week 12: Clone assignment metrics
Week 13: Tree metrics (upcoming)
"""

from .hscn_metrics import (
    compute_hscn_error,
    compute_tcn_mse,
    compute_loh_metrics,
    compute_mirrored_subclone_resolution,
    compute_all_hscn_metrics
)

from .breakpoint_metrics import (
    compute_breakpoint_metrics
)

from .clone_metrics import (
    compute_contingency_matrix,
    hungarian_matching,
    compute_ari,
    compute_nmi,
    handle_unequal_K,
    compute_clone_assignment_metrics
)

__all__ = [
    # HSCN metrics (Week 11)
    'compute_hscn_error',
    'compute_tcn_mse',
    'compute_loh_metrics',
    'compute_mirrored_subclone_resolution',
    'compute_all_hscn_metrics',
    
    # Breakpoint metrics (Week 11)
    'compute_breakpoint_metrics',
    
    # Clone assignment metrics (Week 12)
    'compute_contingency_matrix',
    'hungarian_matching',
    'compute_ari',
    'compute_nmi',
    'handle_unequal_K',
    'compute_clone_assignment_metrics',
]
