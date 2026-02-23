"""
Week 7 Tests: Segment Boundary Detection
"""

import sys
sys.path.insert(0, '/Users/zaaferahhamid/Documents/haplotreesim/src')

import numpy as np
from haplotreesim import (
    SimulationConfig,
    HaploTreeSimulator,
    SegmentDetector,
    Clone
)


def test_segment_detection():
    """Test that segments are detected from CNA breakpoints."""
    print("Testing segment detection...")
    
    config = SimulationConfig(
        num_bins=100,
        num_clones=3,
        num_cells=50,
        lambda_events=2.0,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    read_counts, alleles = sim.run()
    
    # Should have more than 1 segment if events occurred
    assert len(sim.segments) >= 1, "Should have at least 1 segment"
    
    # Check all bins are covered
    total_bins = sum(len(s.bin_indices) for s in sim.segments)
    assert total_bins == config.num_bins, f"All bins should be covered: {total_bins} != {config.num_bins}"
    
    # Check segments don't overlap
    all_bins = set()
    for seg in sim.segments:
        assert len(all_bins & seg.bin_indices) == 0, "Segments should not overlap"
        all_bins.update(seg.bin_indices)
    
    print(f"  ✓ Detected {len(sim.segments)} segments")
    print(f"  ✓ All {config.num_bins} bins covered")
    print(f"  ✓ No overlapping segments")


def test_segment_boundaries():
    """Test that segment boundaries align with CNA breakpoints."""
    print("\nTesting segment boundaries...")
    
    # Create simple test case with known breakpoints
    cn_A_clone0 = np.ones(100, dtype=int)
    cn_B_clone0 = np.ones(100, dtype=int)
    
    clone0 = Clone(
        index=0,
        parent_index=None,
        cn_profile_A=cn_A_clone0,
        cn_profile_B=cn_B_clone0,
        events=[],
        is_root=True
    )
    
    # Clone 1 has a gain at bins 20-40
    cn_A_clone1 = np.ones(100, dtype=int)
    cn_A_clone1[20:41] = 2  # Gain on haplotype A
    cn_B_clone1 = np.ones(100, dtype=int)
    
    clone1 = Clone(
        index=1,
        parent_index=0,
        cn_profile_A=cn_A_clone1,
        cn_profile_B=cn_B_clone1,
        events=[],
        is_root=False
    )
    
    # Detect segments
    detector = SegmentDetector(num_bins=100)
    segments = detector.detect_segments([clone0, clone1])
    
    # Should have 3 segments: [0-19], [20-40], [41-99]
    assert len(segments) == 3, f"Should have 3 segments, got {len(segments)}"
    
    # Check boundaries
    assert segments[0].start_bin == 0 and segments[0].end_bin == 19
    assert segments[1].start_bin == 20 and segments[1].end_bin == 40
    assert segments[2].start_bin == 41 and segments[2].end_bin == 99
    
    print(f"  ✓ Segment boundaries correct")
    print(f"  ✓ Breakpoints at bins 20 and 41 detected")


def test_multiple_clones_segments():
    """Test segment detection with multiple clones."""
    print("\nTesting multiple clone segments...")
    
    config = SimulationConfig(
        num_bins=500,
        num_clones=5,
        num_cells=100,
        lambda_events=2.0,
        random_seed=123
    )
    
    sim = HaploTreeSimulator(config)
    read_counts, alleles = sim.run()
    
    # With 5 clones and lambda_events=2.0, expect multiple segments
    assert len(sim.segments) > 1, "Should have multiple segments with CNAs"
    
    # Check segments have reasonable sizes
    segment_sizes = [len(s.bin_indices) for s in sim.segments]
    assert min(segment_sizes) >= 1, "Segments should have at least 1 bin"
    assert max(segment_sizes) <= config.num_bins, "Segments can't exceed genome size"
    
    print(f"  ✓ Detected {len(sim.segments)} segments")
    print(f"  ✓ Segment sizes range: {min(segment_sizes)}-{max(segment_sizes)} bins")


def test_allele_counts_with_segments():
    """Test that allele counts work correctly with multiple segments."""
    print("\nTesting allele counts with segments...")
    
    config = SimulationConfig(
        num_bins=200,
        num_clones=4,
        num_cells=50,
        lambda_events=1.5,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    read_counts, (alt, ref, total) = sim.run()
    
    # Check shapes
    assert alt.shape == (config.num_cells, len(sim.segments)), \
        f"Alt counts shape mismatch: {alt.shape} vs ({config.num_cells}, {len(sim.segments)})"
    
    # Check that some segments have allele counts
    segments_with_counts = (total > 0).sum(axis=1).mean()
    assert segments_with_counts > 0, "Should have allele counts in some segments"
    
    print(f"  ✓ Allele counts shape: {alt.shape}")
    print(f"  ✓ Average segments with counts per cell: {segments_with_counts:.1f}")


def main():
    """Run all Week 7 tests."""
    print("="*70)
    print("Week 7 Tests: Segment Boundary Detection")
    print("="*70)
    
    test_segment_detection()
    test_segment_boundaries()
    test_multiple_clones_segments()
    test_allele_counts_with_segments()
    
    print("\n" + "="*70)
    print("ALL WEEK 7 TESTS PASSED! ✓")
    print("="*70)


if __name__ == "__main__":
    main()
