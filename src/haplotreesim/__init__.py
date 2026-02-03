"""
HaploTreeSim: A Controlled, End-to-End Benchmark for Low-Pass scDNA-seq

A simulator and benchmark suite that generates synthetic scDNA-seq datasets 
under explicit haplotype-specific CNA evolution on a ground-truth clone tree.
"""

__version__ = "0.1.0"

from .data_models import (
    Bin,
    Segment,
    HaplotypeBlock,
    CNAEvent,
    Clone,
    Cell,
    SimulationConfig
)

from .simulator import HaploTreeSimulator

__all__ = [
    "Bin",
    "Segment", 
    "HaplotypeBlock",
    "CNAEvent",
    "Clone",
    "Cell",
    "SimulationConfig",
    "HaploTreeSimulator",
]
