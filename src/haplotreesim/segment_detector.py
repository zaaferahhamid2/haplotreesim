"""
Segment Boundary Detection for HaploTreeSim.

This module detects segment boundaries from CNA breakpoints across clones.
Implements Section 3.1: "Segments are defined as contiguous bin intervals 
induced by the union of CNA changepoints across clones."
"""

import numpy as np
from typing import List, Set, Tuple
from .data_models import Clone, Segment


class SegmentDetector:
    """
    Detects segment boundaries from clone CNA profiles.
    
    Segments are maximal contiguous regions where no clone has a 
    copy number changepoint.
    """
    
    def __init__(
        self,
        num_bins: int,
        bin_width: int = 1,
        mandatory_breakpoints=None,
        bin_lengths=None,
    ):
        """
        Initialize the segment detector.
        
        Args:
            num_bins: Total number of genomic bins
            bin_width: Width of each bin in base pairs
        """
        self.num_bins = num_bins
        self.bin_width = bin_width
        self.mandatory_breakpoints = set(mandatory_breakpoints or [])
        self.bin_lengths = None if bin_lengths is None else np.asarray(bin_lengths, dtype=int)
        if self.bin_lengths is not None and self.bin_lengths.shape != (self.num_bins,):
            raise ValueError("bin_lengths must have length num_bins")
    
    def detect_segments(self, clones: List[Clone]) -> List[Segment]:
        """
        Detect segment boundaries from clone CN profiles.
        
        Args:
            clones: List of clones with CN profiles
            
        Returns:
            List of Segment objects
        """
        # Find all breakpoints across all clones
        breakpoints = self._find_all_breakpoints(clones)
        
        # Convert breakpoints to segments
        segments = self._breakpoints_to_segments(breakpoints)
        
        print(f"  Detected {len(segments)} segments from {len(breakpoints)} breakpoints")
        
        return segments
    
    def _find_all_breakpoints(self, clones: List[Clone]) -> Set[int]:
        """
        Find all copy number changepoints across all clones.
        
        A breakpoint occurs at position b if any clone has different
        CN at bins b-1 and b (for either haplotype).
        
        Args:
            clones: List of clones
            
        Returns:
            Set of bin indices where breakpoints occur
        """
        breakpoints = {0, self.num_bins}  # Always include start and end
        breakpoints.update(
            breakpoint
            for breakpoint in self.mandatory_breakpoints
            if 0 <= breakpoint <= self.num_bins
        )
        
        for clone in clones:
            # Find breakpoints in haplotype A
            for b in range(1, self.num_bins):
                if clone.cn_profile_A[b] != clone.cn_profile_A[b-1]:
                    breakpoints.add(b)
            
            # Find breakpoints in haplotype B
            for b in range(1, self.num_bins):
                if clone.cn_profile_B[b] != clone.cn_profile_B[b-1]:
                    breakpoints.add(b)
        
        return breakpoints
    
    def _breakpoints_to_segments(self, breakpoints: Set[int]) -> List[Segment]:
        """
        Convert breakpoints to segment objects.
        
        Args:
            breakpoints: Set of bin indices
            
        Returns:
            List of Segment objects
        """
        # Sort breakpoints
        sorted_breakpoints = sorted(breakpoints)
        
        segments = []
        
        # Create segments between consecutive breakpoints
        for i in range(len(sorted_breakpoints) - 1):
            start_bin = sorted_breakpoints[i]
            end_bin = sorted_breakpoints[i + 1] - 1  # End is inclusive
            
            # Skip if this would create an empty segment
            if start_bin > end_bin:
                continue
            
            # Create bin set for this segment
            bin_indices = set(range(start_bin, end_bin + 1))
            
            # Segment length is stored in base pairs for the allele-count model.
            if self.bin_lengths is None:
                segment_length = len(bin_indices) * self.bin_width
            else:
                segment_length = int(self.bin_lengths[start_bin:end_bin + 1].sum())
            
            segment = Segment(
                index=len(segments),
                bin_indices=bin_indices,
                start_bin=start_bin,
                end_bin=end_bin,
                length=segment_length,
                haplotype_block=0  # Will be assigned later
            )
            
            segments.append(segment)
        
        return segments
    
    def get_segment_statistics(self, segments: List[Segment]) -> dict:
        """
        Calculate statistics about detected segments.
        
        Args:
            segments: List of segments
            
        Returns:
            Dictionary of statistics
        """
        if not segments:
            return {
                "num_segments": 0,
                "mean_length": 0,
                "median_length": 0,
                "min_length": 0,
                "max_length": 0
            }
        
        lengths = [len(seg.bin_indices) for seg in segments]
        
        return {
            "num_segments": len(segments),
            "mean_length": np.mean(lengths),
            "median_length": np.median(lengths),
            "min_length": min(lengths),
            "max_length": max(lengths)
        }
