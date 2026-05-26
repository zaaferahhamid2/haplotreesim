"""
Event Generator for HaploTreeSim.

Implements Section 3.4: Event attribute priors including:
- Scale class sampling (focal/arm/chromosome)
- Focal event length from truncated log-uniform
- Amplitude from 1 + Poisson(λ_Δ)
- Haplotype selection
"""

import numpy as np
from typing import List, Optional, Tuple
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
        chromosomes: Optional[List[str]] = None,
        chromosome_intervals: Optional[List[Tuple[str, int, int]]] = None,
        arm_intervals: Optional[List[Tuple[str, str, int, int]]] = None,
        lambda_events: float = 2.0,
        lambda_amplitude: float = 1.0,
        prob_wgd: float = 0.0,
        gain_prob: float = 0.5,
        prob_focal: float = 0.7,
        prob_arm_given_broad: float = 0.75,
        focal_length_min: float = 0.5e6,
        focal_length_max_fraction: float = 0.5,
        prob_haplotype_A: float = 0.5,
        max_copy_number: int = 8,
        max_event_attempts: int = 50
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
            max_copy_number: Maximum allowed copy number (C_max)
            max_event_attempts: Max attempts to draw one valid proposal
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
        self.max_copy_number = max_copy_number
        self.max_event_attempts = max_event_attempts
        
        # Compute scale class weights (Equations 21-23)
        self.w_focal = prob_focal
        self.w_arm = (1 - prob_focal) * prob_arm_given_broad
        self.w_chr = (1 - prob_focal) * (1 - prob_arm_given_broad)
        
        # Genomic intervals are global bin coordinates. Multi-chromosome
        # callers provide these explicitly; single-chromosome callers fall back
        # to intervals derived from hg38 coordinates and the configured bin span.
        from .chromosome_data import get_chromosome_length, get_arm_boundaries, normalize_chromosomes

        self.chromosomes = normalize_chromosomes(chromosomes if chromosomes is not None else [chromosome])
        self.chromosome = self.chromosomes[0]

        if chromosome_intervals is None:
            self.chromosome_intervals = [(self.chromosome, 0, self.num_bins - 1)]
        else:
            self.chromosome_intervals = [
                (chrom, int(start), int(end))
                for chrom, start, end in chromosome_intervals
                if int(start) <= int(end)
            ]

        self._configured_arm_intervals = None
        if arm_intervals is not None:
            self._configured_arm_intervals = [
                (chrom, arm, int(start), int(end))
                for chrom, arm, start, end in arm_intervals
                if int(start) <= int(end)
            ]

        self.genome_length = sum(get_chromosome_length(chrom) for chrom in self.chromosomes)

        p_end, q_start = get_arm_boundaries(self.chromosome)
        self.p_arm_length = p_end
        self.q_arm_length = get_chromosome_length(self.chromosome) - q_start
        self.arm_length = max(self.p_arm_length, self.q_arm_length)
    
    def generate_events(self, num_events: int) -> List[CNAEvent]:
        """
        Reject legacy unordered CNA event generation.

        Non-WGD events require current copy-number state so gains cannot
        amplify absent material and losses can be bounded as they are ordered.
        """
        raise RuntimeError(
            "Unordered event generation is disabled; use generate_valid_events "
            "with current copy-number state."
        )

    def generate_valid_events(
        self,
        num_proposals: int,
        edge_length: float,
        cn_A: np.ndarray,
        cn_B: np.ndarray
    ) -> List[CNAEvent]:
        """
        Generate ordered state-dependent events for one tree edge.

        The input copy-number arrays should already include inherited parent
        state and any WGD applied on this edge. Proposals that would amplify
        absent material or exceed copy-number bounds are rejected; losses are
        capped by the minimum current copy number in the affected interval.
        """
        if num_proposals <= 0:
            return []

        current_A = np.asarray(cn_A, dtype=int).copy()
        current_B = np.asarray(cn_B, dtype=int).copy()
        if len(current_A) != self.num_bins or len(current_B) != self.num_bins:
            raise ValueError("State arrays must have length num_bins")

        event_times = np.sort(self.rng.uniform(0.0, max(edge_length, 0.0), size=num_proposals))
        events: List[CNAEvent] = []

        for event_time in event_times:
            event = self._sample_valid_event(
                current_A=current_A,
                current_B=current_B,
                event_time=float(event_time)
            )
            if event is None:
                continue

            self._apply_event_to_state(current_A, current_B, event)
            events.append(event)

        return events

    def generate_events_for_edge(
        self,
        edge_length: float = 1.0,
        allow_wgd: bool = False
    ) -> List[CNAEvent]:
        """
        Reject legacy unordered edge event generation.

        The simulator owns WGD ordering and supplies inherited copy-number
        state to ``generate_valid_events`` for non-WGD edge events.
        """
        raise RuntimeError(
            "Unordered edge event generation is disabled; use "
            "generate_valid_events with current copy-number state."
        )
    
    def _generate_single_event(self) -> CNAEvent:
        """
        Generate a single CNA event following updated priors.
        
        Returns:
            CNAEvent object
        """
        scale_class, start_bin, end_bin = self._sample_event_interval()
        
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
            amplitude=amplitude,
            scale_class=scale_class
        )

    def _sample_valid_event(
        self,
        current_A: np.ndarray,
        current_B: np.ndarray,
        event_time: float
    ) -> Optional[CNAEvent]:
        """
        Draw one valid ordered event against the current copy-number state.
        """
        for _ in range(self.max_event_attempts):
            scale_class, start_bin, end_bin = self._sample_event_interval()

            is_haplotype_A = self.rng.random() < self.prob_haplotype_A
            haplotype = Haplotype.A if is_haplotype_A else Haplotype.B
            current = current_A if haplotype == Haplotype.A else current_B
            region = current[start_bin:end_bin + 1]
            if region.size == 0:
                continue

            min_current = int(np.min(region))
            if min_current < 1:
                continue

            is_gain = self.rng.random() < self.gain_prob
            proposed_magnitude = 1 + int(self.rng.poisson(self.lambda_amplitude))

            if is_gain:
                if int(np.max(region)) + proposed_magnitude > self.max_copy_number:
                    continue
                amplitude = proposed_magnitude
            else:
                amplitude = -min(proposed_magnitude, min_current)

            return CNAEvent(
                start_bin=start_bin,
                end_bin=end_bin,
                haplotype=haplotype,
                amplitude=amplitude,
                event_time=event_time,
                scale_class=scale_class
            )

        return None

    def _sample_event_interval(self) -> Tuple[str, int, int]:
        """Sample event scale class and genomic interval."""
        scale_class = self.rng.choice(
            ['focal', 'arm', 'chr'],
            p=[self.w_focal, self.w_arm, self.w_chr]
        )

        if scale_class == 'focal':
            start_bin, end_bin = self._sample_focal_event()
        elif scale_class == 'arm':
            start_bin, end_bin = self._sample_arm_event()
        else:
            start_bin, end_bin = self._sample_chromosome_event()

        return scale_class, start_bin, end_bin

    def _apply_event_to_state(
        self,
        current_A: np.ndarray,
        current_B: np.ndarray,
        event: CNAEvent
    ) -> None:
        """Apply an already-validated non-WGD event to mutable state arrays."""
        current = current_A if event.haplotype == Haplotype.A else current_B
        updated = current[event.start_bin:event.end_bin + 1] + event.amplitude

        if np.any(updated < 0) or np.any(updated > self.max_copy_number):
            raise RuntimeError("Validated CNA event would leave copy-number bounds")

        current[event.start_bin:event.end_bin + 1] = updated
    
    def _sample_focal_event(self) -> Tuple[int, int]:
        """
        Sample focal event using truncated log-uniform distribution.
        
        Implements Equations 24-25.
        
        Returns:
            (start_bin, end_bin) tuple
        """
        arm_start, arm_end = self._sample_arm_interval()
        arm_bin_count = arm_end - arm_start + 1

        arm_length_bp = arm_bin_count * self.bin_length
        L_max = min(
            arm_length_bp,
            max(self.bin_length, self.focal_length_max_fraction * arm_length_bp),
        )
        L_min = min(self.focal_length_min, L_max)
        
        # Sample from truncated log-uniform (Equation 24)
        # f(L) = 1/(L * log(L_max/L_min)) for L in [L_min, L_max]
        # Use inverse transform: L = L_min * (L_max/L_min)^U where U ~ Uniform(0,1)
        if L_max <= L_min:
            L = L_min
        else:
            u = self.rng.random()
            L = L_min * (L_max / L_min) ** u
        
        # Convert to bin span (Equation 25)
        m = max(1, int(np.ceil(L / self.bin_length)))
        
        # Restrict each focal placement to the selected chromosome arm.
        m = min(m, arm_bin_count)
        
        # Sample start bin uniformly
        # "choose a valid start bin uniformly among placements"
        max_start = arm_end - m + 1
        start_bin = self.rng.integers(arm_start, max_start + 1)
        end_bin = start_bin + m - 1
        
        return start_bin, end_bin
    
    def _sample_arm_event(self) -> Tuple[int, int]:
        """
        Sample arm-level event using real chromosome arms.
        
        Returns:
            (start_bin, end_bin) tuple
        """
        return self._sample_arm_interval()

    def _sample_arm_interval(self) -> Tuple[int, int]:
        """Choose one simulated chromosome-arm interval weighted by bin count."""
        arm_intervals = self._get_arm_intervals()
        if not arm_intervals:
            return self._sample_chromosome_event()

        weights = np.array([end - start + 1 for start, end in arm_intervals], dtype=float)
        selected_idx = self.rng.choice(len(arm_intervals), p=weights / weights.sum())
        return arm_intervals[selected_idx]

    def _get_arm_intervals(self) -> List[Tuple[int, int]]:
        """Return chromosome-arm intervals that intersect the simulated bins."""
        if self._configured_arm_intervals is not None:
            return [(start, end) for _, _, start, end in self._configured_arm_intervals]

        from .chromosome_data import get_arm_boundaries
        
        p_end, q_start = get_arm_boundaries(self.chromosome)
        
        # Convert to bin indices and keep only arms that intersect this simulation.
        p_arm_end_bin = min(p_end // self.bin_length, self.num_bins - 1)
        q_arm_start_bin = q_start // self.bin_length

        arm_intervals = []
        if p_arm_end_bin >= 0:
            arm_intervals.append((0, p_arm_end_bin))
        if q_arm_start_bin <= self.num_bins - 1:
            arm_intervals.append((q_arm_start_bin, self.num_bins - 1))

        return arm_intervals
    
    def _sample_chromosome_event(self) -> Tuple[int, int]:
        """
        Sample chromosome-level event.
        
        Covers one chromosome in global bin coordinates.
        
        Returns:
            (start_bin, end_bin) tuple
        """
        weights = np.array(
            [end - start + 1 for _, start, end in self.chromosome_intervals],
            dtype=float
        )
        selected_idx = self.rng.choice(len(self.chromosome_intervals), p=weights / weights.sum())
        _, start, end = self.chromosome_intervals[selected_idx]
        return start, end
