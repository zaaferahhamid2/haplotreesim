"""
Event Application for HaploTreeSim.

This module applies CNA events to clone copy number profiles.
"""

import numpy as np
from typing import List
from .data_models import CNAEvent, Clone, Haplotype


class EventApplier:
    """
    Applies CNA events to create child clone profiles from parent profiles.
    
    Implements the state transition rules from Section 3.3 of the paper.
    """
    
    def __init__(self, max_copy_number: int = 8):
        """
        Initialize the event applier.
        
        Args:
            max_copy_number: Maximum allowed copy number (C_max in paper)
        """
        self.max_copy_number = max_copy_number
    
    def apply_events(
        self,
        parent_clone: Clone,
        events: List[CNAEvent]
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Apply a list of events to a parent clone to create child CN profiles.
        
        Args:
            parent_clone: Parent clone
            events: List of events to apply
            
        Returns:
            Tuple of (cn_profile_A, cn_profile_B) for the child clone
        """
        # Start with parent's copy numbers
        cn_A = parent_clone.cn_profile_A.copy()
        cn_B = parent_clone.cn_profile_B.copy()
        
        # Apply each event in order
        for event in events:
            cn_A, cn_B = self._apply_single_event(cn_A, cn_B, event)
        
        return cn_A, cn_B
    
    def _apply_single_event(
        self,
        cn_A: np.ndarray,
        cn_B: np.ndarray,
        event: CNAEvent
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Apply a single event to copy number profiles.
        
        Implements the state transition rules from Section 3.3.
        
        Args:
            cn_A: Haplotype A copy numbers
            cn_B: Haplotype B copy numbers
            event: Event to apply
            
        Returns:
            Updated (cn_A, cn_B)
        """
        if event.is_wgd():
            # WGD: double both haplotypes
            cn_A = self._clip_cn(cn_A * 2)
            cn_B = self._clip_cn(cn_B * 2)
        elif event.haplotype == Haplotype.A:
            # Gain/loss on haplotype A
            cn_A = self._apply_to_region(cn_A, event)
        elif event.haplotype == Haplotype.B:
            # Gain/loss on haplotype B
            cn_B = self._apply_to_region(cn_B, event)
        
        return cn_A, cn_B
    
    def _apply_to_region(
        self,
        cn_profile: np.ndarray,
        event: CNAEvent
    ) -> np.ndarray:
        """
        Apply gain/loss to a specific genomic region.
        
        Implements: c^(γ)_{k,b} = clip(c^(γ)_{π(k),b} + Δ·I[b_min ≤ b ≤ b_max], 0, C_max)
        
        Args:
            cn_profile: Copy number profile for one haplotype
            event: Event specifying region and amplitude
            
        Returns:
            Updated copy number profile
        """
        cn_profile = cn_profile.copy()
        
        # Apply amplitude to affected region
        cn_profile[event.start_bin:event.end_bin + 1] += event.amplitude
        
        # Clip to valid range [0, C_max]
        cn_profile = self._clip_cn(cn_profile)
        
        return cn_profile
    
    def _clip_cn(self, cn_profile: np.ndarray) -> np.ndarray:
        """
        Clip copy numbers to valid range [0, C_max].
        
        Args:
            cn_profile: Copy number profile
            
        Returns:
            Clipped copy number profile
        """
        return np.clip(cn_profile, 0, self.max_copy_number)
