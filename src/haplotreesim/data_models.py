"""
Core data models for HaploTreeSim.

This module defines the fundamental data structures used throughout the simulator:
- Bins: Fixed-width genomic regions
- Segments: Contiguous bin intervals with uniform copy number
- Haplotype Blocks: Groups of segments with shared phasing
- CNA Events: Copy number alteration events on tree edges
- Clones: Nodes in the clone tree with haplotype-specific copy numbers
- Cells: Individual simulated cells with clone assignments
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
from enum import Enum
import numpy as np


class Haplotype(Enum):
    """Enum for haplotype labels."""
    A = "A"
    B = "B"
    WGD = "WGD"  # Whole genome duplication (affects both haplotypes)


@dataclass
class Bin:
    """
    Represents a fixed-width genomic bin.
    
    Attributes:
        index: Bin index (0-based, corresponds to b ∈ {1, ..., B} in paper as b-1)
        chromosome: Chromosome name/number
        start: Start position in base pairs
        end: End position in base pairs
        length: Length in base pairs (ℓ_b in paper)
        gc_content: GC content fraction (optional, for bias modeling)
        mappability: Mappability score (optional, for bias modeling)
    """
    index: int
    chromosome: str
    start: int
    end: int
    length: int
    gc_content: Optional[float] = None
    mappability: Optional[float] = None
    
    def __post_init__(self):
        """Validate bin properties."""
        assert self.start < self.end, "Bin start must be < end"
        assert self.length == self.end - self.start, "Length must match end - start"


@dataclass
class Segment:
    """
    Represents a contiguous interval of bins with uniform copy number.
    
    Segments are defined by the union of CNA changepoints across all clones.
    
    Attributes:
        index: Segment index (corresponds to s ∈ {1, ..., S} in paper as s-1)
        bin_indices: Set of bin indices in this segment (B(s) in paper)
        start_bin: First bin index in segment
        end_bin: Last bin index in segment (inclusive)
        length: Total length in base pairs (L_s in paper)
        haplotype_block: Index of haplotype block containing this segment
    """
    index: int
    bin_indices: Set[int]
    start_bin: int
    end_bin: int
    length: int
    haplotype_block: Optional[int] = None
    
    def __post_init__(self):
        """Validate segment properties."""
        assert self.start_bin <= self.end_bin, "Segment start_bin must be <= end_bin"
        assert self.start_bin in self.bin_indices, "start_bin must be in bin_indices"
        assert self.end_bin in self.bin_indices, "end_bin must be in bin_indices"


@dataclass
class HaplotypeBlock:
    """
    Represents a contiguous group of segments with shared phasing.
    
    Attributes:
        index: Block index (corresponds to h ∈ {1, ..., H} in paper as h-1)
        segment_indices: List of segment indices in this block
        orientation: Phase orientation (+1 or -1, η_h in paper)
        alternate_haplotype: Which haplotype carries the alternate allele (ϕ in paper)
    """
    index: int
    segment_indices: List[int]
    orientation: int = 1  # +1 or -1
    alternate_haplotype: Haplotype = Haplotype.A
    
    def __post_init__(self):
        """Validate haplotype block properties."""
        assert self.orientation in [-1, 1], "Orientation must be +1 or -1"
        assert self.alternate_haplotype in [Haplotype.A, Haplotype.B], \
            "Alternate haplotype must be A or B (not WGD)"


@dataclass
class CNAEvent:
    """
    Represents a copy number alteration event.
    
    Events occur on edges of the clone tree and modify haplotype-specific copy numbers.
    
    Attributes:
        start_bin: First bin affected by event (b_min in paper)
        end_bin: Last bin affected by event (b_max in paper, inclusive)
        haplotype: Which haplotype is affected (γ ∈ {A, B, WGD} in paper)
        amplitude: Magnitude of gain/loss (Δ in paper, ignored for WGD)
        event_id: Unique identifier for this event
    """
    start_bin: int
    end_bin: int
    haplotype: Haplotype
    amplitude: int = 0  # Ignored if haplotype == WGD
    event_id: Optional[str] = None
    
    def __post_init__(self):
        """Validate event properties."""
        assert self.start_bin <= self.end_bin, "Event start_bin must be <= end_bin"
        if self.haplotype != Haplotype.WGD:
            assert self.amplitude != 0, "Non-WGD events must have non-zero amplitude"
        
        # Generate event_id if not provided
        if self.event_id is None:
            hap_str = self.haplotype.value
            self.event_id = f"{hap_str}_{self.start_bin}_{self.end_bin}_{self.amplitude}"
    
    @property
    def length(self) -> int:
        """Return the length of the event in bins (ℓ_e in paper)."""
        return self.end_bin - self.start_bin + 1
    
    def is_gain(self) -> bool:
        """Check if this is a gain event."""
        return self.haplotype != Haplotype.WGD and self.amplitude > 0
    
    def is_loss(self) -> bool:
        """Check if this is a loss event."""
        return self.haplotype != Haplotype.WGD and self.amplitude < 0
    
    def is_wgd(self) -> bool:
        """Check if this is a whole genome duplication."""
        return self.haplotype == Haplotype.WGD


@dataclass
class Clone:
    """
    Represents a clone (node in the clone tree).
    
    Each clone has haplotype-specific copy numbers at each bin and a set of
    events that occurred on the edge from its parent.
    
    Attributes:
        index: Clone index (k ∈ {1, ..., K} in paper as k-1)
        parent_index: Index of parent clone (π(k) in paper, None for root)
        cn_profile_A: Copy number profile for haplotype A (c^(A)_{k,b} in paper)
        cn_profile_B: Copy number profile for haplotype B (c^(B)_{k,b} in paper)
        events: List of CNA events on edge from parent to this clone (E_k in paper)
        is_root: Whether this is the root clone
    """
    index: int
    parent_index: Optional[int]
    cn_profile_A: np.ndarray  # Shape: (B,) where B is number of bins
    cn_profile_B: np.ndarray  # Shape: (B,)
    events: List[CNAEvent] = field(default_factory=list)
    is_root: bool = False
    
    def __post_init__(self):
        """Validate clone properties."""
        assert self.cn_profile_A.shape == self.cn_profile_B.shape, \
            "Haplotype A and B profiles must have same shape"
        if self.is_root:
            assert self.parent_index is None, "Root clone cannot have a parent"
        else:
            assert self.parent_index is not None, "Non-root clone must have a parent"
    
    @property
    def num_bins(self) -> int:
        """Return the number of bins in the genome."""
        return len(self.cn_profile_A)
    
    def total_cn(self) -> np.ndarray:
        """
        Compute total copy number profile (TCN_{k,b} in paper).
        
        Returns:
            Array of shape (B,) with total CN at each bin
        """
        return self.cn_profile_A + self.cn_profile_B
    
    def get_segment_cn(self, segment: Segment) -> Tuple[float, float, float]:
        """
        Compute average haplotype-specific and total CN over a segment.
        
        Args:
            segment: Segment to compute CN for
            
        Returns:
            Tuple of (avg_cn_A, avg_cn_B, avg_total_cn) for the segment
        """
        bin_list = sorted(segment.bin_indices)
        avg_cn_A = np.mean(self.cn_profile_A[bin_list])
        avg_cn_B = np.mean(self.cn_profile_B[bin_list])
        avg_total = avg_cn_A + avg_cn_B
        return avg_cn_A, avg_cn_B, avg_total
    
    def is_loh_at_segment(self, segment: Segment) -> bool:
        """
        Check if segment shows loss of heterozygosity (LOH).
        
        LOH is defined as min(c^(A), c^(B)) = 0 and total CN >= 1.
        
        Args:
            segment: Segment to check
            
        Returns:
            True if segment shows LOH
        """
        cn_A, cn_B, cn_total = self.get_segment_cn(segment)
        return min(cn_A, cn_B) < 0.5 and cn_total >= 1.0


@dataclass
class Cell:
    """
    Represents a simulated single cell.
    
    Attributes:
        index: Cell index (n ∈ {1, ..., N} in paper as n-1)
        clone_assignment: Index of assigned clone (z_n in paper)
        library_size: Library size factor (α_n in paper)
        allelic_coverage: Allelic coverage factor (β_n in paper)
        read_counts: Read depth counts per bin (x_{n,b} in paper), shape (B,)
        allele_counts: Tuple of (alternate, reference, total) counts per segment
                      (a_{n,s}, r_{n,s}, t_{n,s} in paper), each shape (S,)
        is_normal: Whether this is a normal diploid cell
        is_doublet: Whether this is a doublet (fusion of two cells)
        doublet_clones: If doublet, the two clone indices involved
    """
    index: int
    clone_assignment: int
    library_size: float
    allelic_coverage: float
    read_counts: Optional[np.ndarray] = None  # Shape: (B,)
    allele_counts: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None  # (a, r, t), each (S,)
    is_normal: bool = False
    is_doublet: bool = False
    doublet_clones: Optional[Tuple[int, int]] = None
    
    def __post_init__(self):
        """Validate cell properties."""
        assert self.library_size > 0, "Library size must be positive"
        assert self.allelic_coverage > 0, "Allelic coverage must be positive"
        if self.is_doublet:
            assert self.doublet_clones is not None, \
                "Doublet cells must have doublet_clones specified"
            assert len(self.doublet_clones) == 2, \
                "doublet_clones must contain exactly 2 clone indices"


@dataclass
class SimulationConfig:
    """
    Configuration parameters for HaploTreeSim simulation.
    
    This class encapsulates all parameters needed to run a simulation,
    organized by the sections in the paper.
    """
    
    # Genome representation (Section 3.1)
    num_bins: int = 10000  # B in paper
    bin_length: int = 50000  # Length of each bin in bp (50kb default)
    genome_length: int = 500_000_000  # Total genome length (e.g., 500Mb for testing)
    
    # Clone tree and states (Section 3.2)
    num_clones: int = 5  # K in paper
    max_copy_number: int = 8  # C_max in paper
    root_type: str = "diploid"  # "diploid" or "wgd" - root initialization
    
    # CNA event model (Section 3.3)
    lambda_events: float = 1.5  # λ_E: mean number of events per edge
    lambda_amplitude: float = 1.0  # λ_Δ: parameter for amplitude distribution
    prob_wgd: float = 0.0  # p_WGD: probability of WGD event
    prob_mirror: float = 0.0  # p_mirror: probability of mirrored subclones
    gain_prob: float = 0.5  # Probability that an event is a gain (vs loss)
    
    # Event size distribution (focal/arm/chromosomal)
    focal_prob: float = 0.7  # Probability of focal event
    arm_prob: float = 0.2  # Probability of arm-level event
    chrom_prob: float = 0.1  # Probability of chromosomal event
    focal_size_mean: int = 50  # Mean size of focal events in bins
    arm_size_mean: int = 500  # Mean size of arm events in bins
    
    # Cell sampling (Section 3.4)
    num_cells: int = 200  # N in paper
    prob_normal: float = 0.0  # p_normal: fraction of normal cells
    prob_doublet: float = 0.0  # p_doublet: fraction of doublets
    
    # Read-depth observation model (Section 3.5)
    mean_library_size: float = 100.0  # Mean value for α_n
    library_size_cv: float = 0.3  # Coefficient of variation for library sizes
    theta_x: float = 10.0  # Overdispersion parameter for read counts
    use_gc_bias: bool = False  # Whether to model GC bias (κ_b)
    
    # Allelic observation model (Section 3.6)
    snp_density: float = 0.001  # ρ_SNP: heterozygous SNPs per bp
    mean_allelic_coverage: float = 50.0  # Mean value for β_n
    allelic_coverage_cv: float = 0.3  # Coefficient of variation for allelic coverage
    nu_a: float = 20.0  # Concentration parameter for Beta-Binomial
    prob_phase_switch: float = 0.01  # p_switch: phase switch probability
    
    # Random seed
    random_seed: Optional[int] = None
    
    def __post_init__(self):
        """Validate configuration parameters."""
        assert self.num_bins > 0, "num_bins must be positive"
        assert self.num_clones > 0, "num_clones must be positive"
        assert self.num_cells > 0, "num_cells must be positive"
        assert self.root_type in ["diploid", "wgd"], "root_type must be 'diploid' or 'wgd'"
        assert 0 <= self.prob_normal <= 1, "prob_normal must be in [0, 1]"
        assert 0 <= self.prob_doublet <= 1, "prob_doublet must be in [0, 1]"
        assert 0 <= self.prob_wgd <= 1, "prob_wgd must be in [0, 1]"
        assert 0 <= self.prob_mirror <= 1, "prob_mirror must be in [0, 1]"
        assert abs(self.focal_prob + self.arm_prob + self.chrom_prob - 1.0) < 1e-6, \
            "Event size probabilities must sum to 1"
    
    def get_root_cn(self) -> Tuple[int, int]:
        """
        Get the root copy number for both haplotypes.
        
        Returns:
            Tuple of (cn_A, cn_B) for root clone
        """
        if self.root_type == "diploid":
            return (1, 1)
        else:  # wgd
            return (2, 2)
