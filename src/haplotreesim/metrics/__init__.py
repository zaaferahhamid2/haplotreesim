"""
Evaluation metrics for HaploTreeSim benchmark.
"""

from .hscn_metrics import (
    compute_hscn_error,
    compute_tcn_mse,
    compute_loh_metrics,
    find_mirrored_subclone_triples,
    compute_mirrored_subclone_resolution,
    compute_all_hscn_metrics
)

from .breakpoint_metrics import (
    breakpoints_from_frequencies,
    compute_breakpoint_metrics,
    compute_breakpoint_sensitivity_curve,
    compute_hscn_breakpoint_frequencies,
    compute_recurrent_breakpoint_metrics,
    compute_tcn_breakpoint_frequencies,
    predicted_breakpoints_from_hscn,
    predicted_breakpoints_from_tcn,
)

from .clone_metrics import (
    compute_contingency_matrix,
    hungarian_matching,
    compute_ari,
    compute_nmi,
    handle_unequal_K,
    compute_clone_assignment_metrics
)

from .tree_metrics import (
    compute_robinson_foulds_distance,
    compute_ancestor_descendant_accuracy,
    compute_event_matching,
    compute_all_tree_metrics
)

__all__ = [
    'compute_hscn_error',
    'compute_tcn_mse',
    'compute_loh_metrics',
    'find_mirrored_subclone_triples',
    'compute_mirrored_subclone_resolution',
    'compute_all_hscn_metrics',
    'breakpoints_from_frequencies',
    'compute_breakpoint_metrics',
    'compute_breakpoint_sensitivity_curve',
    'compute_hscn_breakpoint_frequencies',
    'compute_recurrent_breakpoint_metrics',
    'compute_tcn_breakpoint_frequencies',
    'predicted_breakpoints_from_hscn',
    'predicted_breakpoints_from_tcn',
    'compute_contingency_matrix',
    'hungarian_matching',
    'compute_ari',
    'compute_nmi',
    'handle_unequal_K',
    'compute_clone_assignment_metrics',
    'compute_robinson_foulds_distance',
    'compute_ancestor_descendant_accuracy',
    'compute_event_matching',
    'compute_all_tree_metrics',
]
