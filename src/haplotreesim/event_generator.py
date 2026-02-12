"""
CNA Event Generation for HaploTreeSim.

This module handles the generation of copy number alteration events
according to the event model in Section 3.3 of the paper.
"""

import numpy as np
from typing import List, Tuple, Optional
from .data_models import CNAEvent, Haplotype


class EventGenerator:
    """
    Generates CNA events on tree edges.
    
    Attributes:
        rng: NumPy random number generator
        num_bins: Total number of genomic bins
        lambda_events: Mean number of events per edge (λ_E in paper)
        lambda_amplitude: Parameter for amplitude distribution (λ_Δ in paper)
        prob_wgd: Probability of WGD event (p_WGD in paper)
        focal_prob: Probability of focal event
        arm_prob: Probability of arm-level event
        chrom_prob: Probability of chromosomal event
        focal_size_mean: Mean size of focal events in bins
        arm_size_mean: Mean size of arm events in bins
        gain_prob: Probability that an event is a gain (vs loss)
    """
    
    def __init__(
        self,
        rng: np.random.Generator,
        num_bins: int,
        lambda_events: float = 1.5,
        lambda_amplitude: float = 1.0,
        prob_wgd: float = 0.0,
        focal_prob: float = 0.7,
        arm_prob: float = 0.2,
        chrom_prob: float = 0.1,
        focal_size_mean: int = 50,
        arm_size_mean: int = 500,
        gain_prob: float = 0.5,
    ):
        """Initialize the event generator."""
        self.rng = rng
        self.num_bins = num_bins
        self.lambda_events = lambda_events
        self.lambda_amplitude = lambda_amplitude
        self.prob_wgd = prob_wgd
        self.focal_prob = focal_prob
        self.arm_prob = arm_prob
        self.chrom_prob = chrom_prob
        self.focal_size_mean = focal_size_mean
        self.arm_size_mean = arm_size_mean
        self.gain_prob = gain_prob
        
        # Validate probabilities sum to 1
        total_prob = focal_prob + arm_prob + chrom_prob
        assert abs(total_prob - 1.0) < 1e-6, \
            f"Event size probs must sum to 1, got {total_prob}"
    
    def generate_events_for_edge(
        self,
        allow_wgd: bool = False
    ) -> List[CNAEvent]:
        """
        Generate CNA events for a single tree edge.
        
        Args:
            allow_wgd: Whether to allow WGD events on this edge
            
        Returns:
            List of CNAEvent objects
        """
        events = []
        
        # Sample number of events: |E_k| ~ Poisson(λ_E)
        num_events = self.rng.poisson(self.lambda_events)
        
        # Check for WGD (replaces all other events if it occurs)
        if allow_wgd and self.rng.random() < self.prob_wgd:
            wgd_event = CNAEvent(
                start_bin=0,
                end_bin=self.num_bins - 1,
                haplotype=Haplotype.WGD,
                amplitude=0  # Ignored for WGD
            )
            return [wgd_event]
        
        # Generate regular gain/loss events
        for i in range(num_events):
            event = self._generate_single_event()
            events.append(event)
        
        return events
    
    def _generate_single_event(self) -> CNAEvent:
        """Generate a single gain or loss event."""
        # Sample event size category
        size_category = self.rng.choice(
            ['focal', 'arm', 'chromosomal'],
            p=[self.focal_prob, self.arm_prob, self.chrom_prob]
        )
        
        # Sample event length based on category
        if size_category == 'focal':
            length = max(1, int(self.rng.poisson(self.focal_size_mean)))
        elif size_category == 'arm':
            length = max(1, int(self.rng.poisson(self.arm_size_mean)))
        else:  # chromosomal
            length = self.num_bins  # Entire genome
        
        # Ensure length doesn't exceed genome
        length = min(length, self.num_bins)
        
        # Sample start position
        max_start = self.num_bins - length
        start_bin = self.rng.integers(0, max(1, max_start + 1))
        end_bin = min(start_bin + length - 1, self.num_bins - 1)
        
        # Sample haplotype (A or B)
        haplotype = self.rng.choice([Haplotype.A, Haplotype.B])
        
        # Sample amplitude: |Δ| ~ 1 + Poisson(λ_Δ)
        amplitude_magnitude = 1 + self.rng.poisson(self.lambda_amplitude)
        
        # Determine sign (gain vs loss)
        is_gain = self.rng.random() < self.gain_prob
        amplitude = amplitude_magnitude if is_gain else -amplitude_magnitude
        
        return CNAEvent(
            start_bin=start_bin,
            end_bin=end_bin,
            haplotype=haplotype,
            amplitude=amplitude
        )
    
    def generate_mirrored_events(
        self,
        base_event: CNAEvent
    ) -> Tuple[CNAEvent, CNAEvent]:
        """
        Generate a pair of mirrored subclone events.
        
        Mirrored events affect the same genomic region with the same
        amplitude, but on different haplotypes.
        
        Args:
            base_event: Template event
            
        Returns:
            Tuple of (event_A, event_B) for the two sibling clones
        """
        event_A = CNAEvent(
            start_bin=base_event.start_bin,
            end_bin=base_event.end_bin,
            haplotype=Haplotype.A,
            amplitude=base_event.amplitude
        )
        
        event_B = CNAEvent(
            start_bin=base_event.start_bin,
            end_bin=base_event.end_bin,
            haplotype=Haplotype.B,
            amplitude=base_event.amplitude
        )
        
        return event_A, event_B
