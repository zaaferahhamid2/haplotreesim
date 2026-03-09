"""
HaploTreeSim: A Controlled, End-to-End Benchmark for Low-Pass scDNA-seq
"""

__version__ = "0.4.0"

from .data_models import (
    Bin,
    Segment,
    HaplotypeBlock,
    CNAEvent,
    Clone,
    Cell,
    SimulationConfig,
    Haplotype,
)

from .simulator import HaploTreeSimulator
from .event_generator import EventGenerator
from .event_applier import EventApplier
from .tree_builder import TreeBuilder
from .beta_tree_builder import BetaSplittingTreeBuilder
from .segment_detector import SegmentDetector

__all__ = [
    "Bin",
    "Segment", 
    "HaplotypeBlock",
    "CNAEvent",
    "Clone",
    "Cell",
    "SimulationConfig",
    "Haplotype",
    "HaploTreeSimulator",
    "EventGenerator",
    "EventApplier",
    "TreeBuilder",
    "BetaSplittingTreeBuilder",
    "SegmentDetector",
]
