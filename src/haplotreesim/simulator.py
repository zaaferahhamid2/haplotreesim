"""
HaploTreeSim simulator core.

This module implements the main simulator class that generates synthetic
scDNA-seq data with haplotype-specific CNAs and clone tree ground truth.
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
from .data_models import (
    Bin, Segment, HaplotypeBlock, CNAEvent, Clone, Cell,
    SimulationConfig, Haplotype
)


class HaploTreeSimulator:
    """
    Main simulator class for HaploTreeSim.
    
    This class orchestrates the simulation process:
    1. Initialize genome (bins, segments, haplotype blocks)
    2. Generate clone tree and haplotype-specific CN profiles
    3. Sample cells from clones
    4. Generate observations (read counts and allele counts)
    
    Attributes:
        config: Simulation configuration parameters
        rng: NumPy random number generator
        bins: List of genomic bins
        segments: List of segments
        haplotype_blocks: List of haplotype blocks
        clones: List of clones (tree nodes)
        cells: List of simulated cells
    """
    
    def __init__(self, config: SimulationConfig):
        """
        Initialize the simulator.
        
        Args:
            config: Simulation configuration
        """
        self.config = config
        self.rng = np.random.default_rng(config.random_seed)
        
        # Data structures (to be populated)
        self.bins: List[Bin] = []
        self.segments: List[Segment] = []
        self.haplotype_blocks: List[HaplotypeBlock] = []
        self.clones: List[Clone] = []
        self.cells: List[Cell] = []
        
        # Mappings
        self.bin_to_segment: Dict[int, int] = {}  # bin index -> segment index
        self.segment_to_block: Dict[int, int] = {}  # segment index -> block index
    
    def run(self) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        Run the full simulation pipeline.
        
        Returns:
            Tuple of:
                - read_counts: Array of shape (N, B) with read counts
                - allele_counts: Tuple of (alternate, reference, total) arrays,
                                each of shape (N, S)
        """
        print("Initializing genome...")
        self._initialize_genome()
        
        print(f"Generating clone tree with {self.config.num_clones} clones...")
        self._generate_clone_tree()
        
        print(f"Sampling {self.config.num_cells} cells...")
        self._sample_cells()
        
        print("Generating read-depth observations...")
        read_counts = self._generate_read_counts()
        
        print("Generating allelic observations...")
        allele_counts = self._generate_allele_counts()
        
        print("Simulation complete!")
        return read_counts, allele_counts
    
    def _initialize_genome(self):
        """
        Initialize genome representation (bins, segments, haplotype blocks).
        
        For Week 5 deliverable: creates uniform bins with a single segment
        covering the entire genome.
        """
        # Create bins
        self.bins = []
        for i in range(self.config.num_bins):
            start = i * self.config.bin_length
            end = start + self.config.bin_length
            bin_obj = Bin(
                index=i,
                chromosome="chr1",  # Simplified: single chromosome
                start=start,
                end=end,
                length=self.config.bin_length
            )
            self.bins.append(bin_obj)
        
        # Create a single segment covering all bins (for now)
        # In future weeks, segments will be defined by CNA breakpoints
        segment = Segment(
            index=0,
            bin_indices=set(range(self.config.num_bins)),
            start_bin=0,
            end_bin=self.config.num_bins - 1,
            length=self.config.num_bins * self.config.bin_length,
            haplotype_block=0
        )
        self.segments = [segment]
        
        # Create a single haplotype block
        hap_block = HaplotypeBlock(
            index=0,
            segment_indices=[0],
            orientation=1,
            alternate_haplotype=Haplotype.A
        )
        self.haplotype_blocks = [hap_block]
        
        # Build mappings
        for bin_idx in range(self.config.num_bins):
            self.bin_to_segment[bin_idx] = 0
        self.segment_to_block[0] = 0
        
        print(f"  Created {len(self.bins)} bins, {len(self.segments)} segments, "
              f"{len(self.haplotype_blocks)} haplotype blocks")
    
    def _generate_clone_tree(self):
        """
        Generate clone tree with haplotype-specific copy number profiles.
        
        For Week 5 deliverable: creates a single root clone with diploid CN
        (no CNAs applied yet).
        """
        # Get root copy number
        cn_A_init, cn_B_init = self.config.get_root_cn()
        
        # Create root clone (diploid across all bins)
        root_cn_A = np.full(self.config.num_bins, cn_A_init, dtype=int)
        root_cn_B = np.full(self.config.num_bins, cn_B_init, dtype=int)
        
        root_clone = Clone(
            index=0,
            parent_index=None,
            cn_profile_A=root_cn_A,
            cn_profile_B=root_cn_B,
            events=[],
            is_root=True
        )
        self.clones = [root_clone]
        
        # TODO (future weeks): Generate tree structure and apply CNA events
        # For now, just create additional clones as copies of root
        for k in range(1, self.config.num_clones):
            clone = Clone(
                index=k,
                parent_index=0,  # All children of root for now
                cn_profile_A=root_cn_A.copy(),
                cn_profile_B=root_cn_B.copy(),
                events=[],
                is_root=False
            )
            self.clones.append(clone)
        
        print(f"  Created {len(self.clones)} clones (all diploid for now)")
    
    def _sample_cells(self):
        """
        Sample cells from clones with optional normal/doublet contamination.
        
        For Week 5: uniform sampling from clones, no contamination.
        """
        # Compute clone frequencies (uniform for now)
        clone_frequencies = np.ones(self.config.num_clones) / self.config.num_clones
        
        # Sample library sizes and allelic coverage factors
        library_sizes = self._sample_library_sizes(self.config.num_cells)
        allelic_coverages = self._sample_allelic_coverages(self.config.num_cells)
        
        # Assign cells to clones
        clone_assignments = self.rng.choice(
            self.config.num_clones,
            size=self.config.num_cells,
            p=clone_frequencies
        )
        
        # Create cell objects
        self.cells = []
        for n in range(self.config.num_cells):
            cell = Cell(
                index=n,
                clone_assignment=clone_assignments[n],
                library_size=library_sizes[n],
                allelic_coverage=allelic_coverages[n],
                is_normal=False,
                is_doublet=False
            )
            self.cells.append(cell)
        
        print(f"  Sampled {len(self.cells)} cells")
        print(f"  Clone distribution: {np.bincount(clone_assignments)}")
    
    def _sample_library_sizes(self, num_cells: int) -> np.ndarray:
        """
        Sample library size factors (α_n) for cells.
        
        Uses log-normal distribution to ensure positive values.
        """
        mean = self.config.mean_library_size
        cv = self.config.library_size_cv
        
        # Log-normal parameters
        sigma = np.sqrt(np.log(1 + cv**2))
        mu = np.log(mean) - 0.5 * sigma**2
        
        return self.rng.lognormal(mu, sigma, size=num_cells)
    
    def _sample_allelic_coverages(self, num_cells: int) -> np.ndarray:
        """
        Sample allelic coverage factors (β_n) for cells.
        
        Uses log-normal distribution to ensure positive values.
        """
        mean = self.config.mean_allelic_coverage
        cv = self.config.allelic_coverage_cv
        
        # Log-normal parameters
        sigma = np.sqrt(np.log(1 + cv**2))
        mu = np.log(mean) - 0.5 * sigma**2
        
        return self.rng.lognormal(mu, sigma, size=num_cells)
    
    def _generate_read_counts(self) -> np.ndarray:
        """
        Generate read-depth counts (x_{n,b}) via negative binomial model.
        
        Returns:
            Array of shape (N, B) with read counts
        """
        N = self.config.num_cells
        B = self.config.num_bins
        read_counts = np.zeros((N, B), dtype=int)
        
        for n, cell in enumerate(self.cells):
            clone = self.clones[cell.clone_assignment]
            
            for b in range(B):
                # Total copy number at this bin
                tcn = clone.cn_profile_A[b] + clone.cn_profile_B[b]
                
                # Bin bias (κ_b, set to 1 for now)
                kappa_b = 1.0
                
                # Mean read count: μ_{n,b} = α_n * κ_b * (TCN / 2)
                mu = cell.library_size * kappa_b * (tcn / 2.0)
                
                # Negative binomial: variance = μ + μ²/θ
                # We use the (n, p) parameterization where:
                #   n = θ, p = θ/(θ + μ)
                if mu > 0:
                    n_param = self.config.theta_x
                    p_param = n_param / (n_param + mu)
                    read_counts[n, b] = self.rng.negative_binomial(n_param, p_param)
                else:
                    read_counts[n, b] = 0
            
            # Store in cell object
            cell.read_counts = read_counts[n, :]
        
        return read_counts
    
    def _generate_allele_counts(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate allele counts (a_{n,s}, r_{n,s}) via Beta-Binomial model.
        
        Returns:
            Tuple of (alternate, reference, total) arrays, each of shape (N, S)
        """
        N = self.config.num_cells
        S = len(self.segments)
        
        alternate_counts = np.zeros((N, S), dtype=int)
        reference_counts = np.zeros((N, S), dtype=int)
        total_counts = np.zeros((N, S), dtype=int)
        
        for n, cell in enumerate(self.cells):
            clone = self.clones[cell.clone_assignment]
            
            for s, segment in enumerate(self.segments):
                # Sample number of heterozygous SNPs in this segment
                M_s = self.rng.poisson(self.config.snp_density * segment.length)
                
                if M_s == 0:
                    continue
                
                # Get average copy numbers over segment
                cn_A, cn_B, tcn = clone.get_segment_cn(segment)
                
                # Total allelic depth: Poisson(β_n * M_s * TCN / 2)
                if tcn > 0:
                    t_ns = self.rng.poisson(cell.allelic_coverage * M_s * (tcn / 2.0))
                else:
                    t_ns = 0
                
                if t_ns == 0:
                    continue
                
                # Expected alternate fraction (assuming phase orientation = 1 for now)
                hap_block = self.haplotype_blocks[segment.haplotype_block]
                if hap_block.alternate_haplotype == Haplotype.A:
                    p_alt = cn_A / (tcn + 1e-10)
                else:
                    p_alt = cn_B / (tcn + 1e-10)
                
                # Apply phase orientation
                if hap_block.orientation == -1:
                    p_alt = 1.0 - p_alt
                
                # Beta-Binomial: sample q ~ Beta(p*ν, (1-p)*ν), then a ~ Binomial(t, q)
                alpha_beta = p_alt * self.config.nu_a
                beta_beta = (1.0 - p_alt) * self.config.nu_a
                
                q = self.rng.beta(alpha_beta, beta_beta)
                a_ns = self.rng.binomial(t_ns, q)
                r_ns = t_ns - a_ns
                
                alternate_counts[n, s] = a_ns
                reference_counts[n, s] = r_ns
                total_counts[n, s] = t_ns
            
            # Store in cell object
            cell.allele_counts = (
                alternate_counts[n, :],
                reference_counts[n, :],
                total_counts[n, :]
            )
        
        return alternate_counts, reference_counts, total_counts
    
    def get_ground_truth(self) -> Dict:
        """
        Extract ground truth data for evaluation.
        
        Returns:
            Dictionary containing:
                - clone_tree: Adjacency list representation of tree
                - clone_assignments: Array of shape (N,) with cell->clone mapping
                - cn_profiles_A: Array of shape (K, B) with haplotype A CNs
                - cn_profiles_B: Array of shape (K, B) with haplotype B CNs
                - segments: List of segment objects
                - events: Dict mapping clone_idx -> list of events
        """
        # Clone tree as adjacency list
        clone_tree = {}
        for clone in self.clones:
            if not clone.is_root:
                parent = clone.parent_index
                if parent not in clone_tree:
                    clone_tree[parent] = []
                clone_tree[parent].append(clone.index)
        
        # Clone assignments
        clone_assignments = np.array([cell.clone_assignment for cell in self.cells])
        
        # CN profiles
        cn_profiles_A = np.vstack([clone.cn_profile_A for clone in self.clones])
        cn_profiles_B = np.vstack([clone.cn_profile_B for clone in self.clones])
        
        # Events
        events = {clone.index: clone.events for clone in self.clones}
        
        return {
            "clone_tree": clone_tree,
            "clone_assignments": clone_assignments,
            "cn_profiles_A": cn_profiles_A,
            "cn_profiles_B": cn_profiles_B,
            "segments": self.segments,
            "events": events,
        }
