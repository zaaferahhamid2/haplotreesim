"""
Event Generator for HaploTreeSim.

Implements Section 3.4: Event attribute priors including:
- Scale class sampling (focal/arm/chromosome)
- Focal event length from truncated log-uniform
- Amplitude from 1 + Poisson(λ_Δ)
- Haplotype selection
"""

import numpy as np
from typing import List, Tuple
from .data_models import CNAEvent, Haplotype


class EventGenerator:
    """
    Generates CNA events with updated event attribute priors.
    
    Implements Equations 20-28 from Section 3.4.
    """
    
    def __init__(
        self,
        rng: np.random.Generator,
        num_bins: int,
        bin_length: int = 500000,  # 500 kb bins default
        chromosome: str = 'chr1',  # Which chromosome
        lambda_events: float = 2.0,
        lambda_amplitude: float = 1.0,
        prob_wgd: float = 0.0,
        gain_prob: float = 0.5,
        prob_focal: float = 0.7,
        prob_arm_given_broad: float = 0.75,
        focal_length_min: float = 0.5e6,
        focal_length_max_fraction: float = 0.5,
        prob_haplotype_A: float = 0.5
    ):
        """
        Initialize event generator with updated priors.
        
        Args:
            rng: Random number generator
            num_bins: Total number of genomic bins
            bin_length: Length of each bin in base pairs
            lambda_events: Event rate parameter (λ_E)
            lambda_amplitude: Amplitude parameter (λ_Δ)
            prob_wgd: Probability of WGD
            gain_prob: Probability of gain vs loss (p_gain)
            prob_focal: Focal event probability (p_focal, Eq 21)
            prob_arm_given_broad: Arm prob among broad events (q_arm, Eq 22)
            focal_length_min: Minimum focal event length in bp
            focal_length_max_fraction: Max focal as fraction of arm length
            prob_haplotype_A: Probability of haplotype A (p_A, Eq 28)
        """
        self.rng = rng
        self.num_bins = num_bins
        self.bin_length = bin_length
        self.lambda_events = lambda_events
        self.lambda_amplitude = lambda_amplitude
        self.prob_wgd = prob_wgd
        self.gain_prob = gain_prob
        self.prob_focal = prob_focal
        self.prob_arm_given_broad = prob_arm_given_broad
        self.focal_length_min = max(bin_length, focal_length_min)
        self.focal_length_max_fraction = focal_length_max_fraction
        self.prob_haplotype_A = prob_haplotype_A
        
        # Compute scale class weights (Equations 21-23)
        self.w_focal = prob_focal
        self.w_arm = (1 - prob_focal) * prob_arm_given_broad
        self.w_chr = (1 - prob_focal) * (1 - prob_arm_given_broad)
        
        # Get real chromosome arm lengths
        from .chromosome_data import get_chromosome_length, get_arm_boundaries
        
        self.chromosome = chromosome
        self.genome_length = get_chromosome_length(chromosome)
        
        # Get centromere position
        p_end, q_start = get_arm_boundaries(chromosome)
        self.p_arm_length = p_end
        self.q_arm_length = self.genome_length - q_start
        self.arm_length = max(self.p_arm_length, self.q_arm_length)
    
    def generate_events(self, num_events: int) -> List[CNAEvent]:
        """
        Generate a list of CNA events.
        
        Args:
            num_events: Number of events to generate
            
        Returns:
            List of CNAEvent objects
        """
        events = []
        for _ in range(num_events):
            event = self._generate_single_event()
            events.append(event)
        return events
    
    def _generate_single_event(self) -> CNAEvent:
        """
        Generate a single CNA event following updated priors.
        
        Returns:
            CNAEvent object
        """
        # Sample scale class (Equation 20)
        scale_class = self.rng.choice(
            ['focal', 'arm', 'chr'],
            p=[self.w_focal, self.w_arm, self.w_chr]
        )
        
        # Generate event based on scale class
        if scale_class == 'focal':
            start_bin, end_bin = self._sample_focal_event()
        elif scale_class == 'arm':
            start_bin, end_bin = self._sample_arm_event()
        else:  # chr
            start_bin, end_bin = self._sample_chromosome_event()
        
        # Sample amplitude (Equations 26-27)
        amplitude_magnitude = 1 + self.rng.poisson(self.lambda_amplitude)
        
        # Sample sign
        is_gain = self.rng.random() < self.gain_prob
        amplitude = amplitude_magnitude if is_gain else -amplitude_magnitude
        
        # Sample haplotype (Equation 28)
        is_haplotype_A = self.rng.random() < self.prob_haplotype_A
        haplotype = Haplotype.A if is_haplotype_A else Haplotype.B
        
        return CNAEvent(
            start_bin=start_bin,
            end_bin=end_bin,
            haplotype=haplotype,
            amplitude=amplitude
        )
    
    def _sample_focal_event(self) -> Tuple[int, int]:
        """
        Sample focal event using truncated log-uniform distribution.
        
        Implements Equations 24-25.
        
        Returns:
            (start_bin, end_bin) tuple
        """
        # L_max = 0.5 * arm_length (Equation 24)
        L_max = self.focal_length_max_fraction * self.arm_length
        L_min = self.focal_length_min
        
        # Sample from truncated log-uniform (Equation 24)
        # f(L) = 1/(L * log(L_max/L_min)) for L in [L_min, L_max]
        # Use inverse transform: L = L_min * (L_max/L_min)^U where U ~ Uniform(0,1)
        u = self.rng.random()
        L = L_min * (L_max / L_min) ** u
        
        # Convert to bin span (Equation 25)
        m = max(1, int(np.ceil(L / self.bin_length)))
        
        # Ensure m doesn't exceed genome
        m = min(m, self.num_bins)
        
        # Sample start bin uniformly
        # "choose a valid start bin uniformly among placements"
        max_start = max(0, self.num_bins - m)
        start_bin = self.rng.integers(0, max_start + 1)
        end_bin = min(start_bin + m - 1, self.num_bins - 1)
        
        return start_bin, end_bin
    
    def _sample_arm_event(self) -> Tuple[int, int]:
        """
        Sample arm-level event using real chromosome arms.
        
        Returns:
            (start_bin, end_bin) tuple
        """
        from .chromosome_data import get_arm_boundaries
        
        p_end, q_start = get_arm_boundaries(self.chromosome)
        
        # Convert to bin indices
        p_arm_end_bin = p_end // self.bin_length
        q_arm_start_bin = q_start // self.bin_length
        
        # Choose p or q arm with probability proportional to arm length
        p_prob = self.p_arm_length / (self.p_arm_length + self.q_arm_length)
        
        if self.rng.random() < p_prob:
            # P arm
            start_bin = 0
            end_bin = min(p_arm_end_bin, self.num_bins - 1)
        else:
            # Q arm
            start_bin = max(q_arm_start_bin, 0)
            end_bin = self.num_bins - 1
        
        return start_bin, end_bin
    
    def _sample_chromosome_event(self) -> Tuple[int, int]:
        """
        Sample chromosome-level event.
        
        In simplified model, covers entire genome.
        
        Returns:
            (start_bin, end_bin) tuple
        """
        return 0, self.num_bins - 1
