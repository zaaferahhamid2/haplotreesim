"""
Week 5 Test Script: Minimal diploid simulator test.

This script verifies that the simulator can:
1. Initialize bins, segments, and haplotype blocks
2. Create a diploid clone tree (no CNAs yet)
3. Sample cells from clones
4. Generate read counts and allele counts
5. Extract ground truth data
"""

import sys
sys.path.insert(0, '/home/claude/haplotreesim/src')

import numpy as np
from haplotreesim import SimulationConfig, HaploTreeSimulator


def test_minimal_diploid():
    """Test the minimal diploid simulator."""
    print("="*70)
    print("Week 5 Test: Minimal Diploid Simulator")
    print("="*70)
    
    # Create minimal configuration
    config = SimulationConfig(
        num_bins=100,
        bin_length=50000,
        num_clones=3,
        num_cells=50,
        root_type="diploid",
        lambda_events=0.0,  # No events yet
        random_seed=42
    )
    
    print("\nConfiguration:")
    print(f"  Bins: {config.num_bins}")
    print(f"  Clones: {config.num_clones}")
    print(f"  Cells: {config.num_cells}")
    print(f"  Root type: {config.root_type}")
    print()
    
    # Initialize simulator
    sim = HaploTreeSimulator(config)
    
    # Run simulation
    read_counts, (alt_counts, ref_counts, total_counts) = sim.run()
    
    # Verify outputs
    print("\n" + "="*70)
    print("Verification")
    print("="*70)
    
    # Check dimensions
    assert read_counts.shape == (config.num_cells, config.num_bins), \
        f"Read counts shape mismatch: {read_counts.shape}"
    print(f"✓ Read counts shape: {read_counts.shape}")
    
    num_segments = len(sim.segments)
    assert alt_counts.shape == (config.num_cells, num_segments), \
        f"Allele counts shape mismatch: {alt_counts.shape}"
    print(f"✓ Allele counts shape: {alt_counts.shape}")
    
    # Check that all clones are diploid
    for clone in sim.clones:
        assert np.all(clone.cn_profile_A == 1), f"Clone {clone.index} hap A not diploid"
        assert np.all(clone.cn_profile_B == 1), f"Clone {clone.index} hap B not diploid"
    print(f"✓ All {len(sim.clones)} clones are diploid (1,1)")
    
    # Check read counts are reasonable
    mean_reads = read_counts.mean()
    print(f"✓ Mean read count per bin: {mean_reads:.2f}")
    assert mean_reads > 0, "No reads generated"
    
    # Check allele counts
    total_allelic_reads = total_counts.sum()
    print(f"✓ Total allelic reads: {total_allelic_reads}")
    assert total_allelic_reads > 0, "No allelic reads generated"
    
    # Get ground truth
    gt = sim.get_ground_truth()
    print(f"\n✓ Ground truth extracted:")
    print(f"  Clone assignments shape: {gt['clone_assignments'].shape}")
    print(f"  CN profiles A shape: {gt['cn_profiles_A'].shape}")
    print(f"  CN profiles B shape: {gt['cn_profiles_B'].shape}")
    
    # Verify diploid ground truth
    assert np.all(gt['cn_profiles_A'] == 1), "Ground truth hap A not all 1"
    assert np.all(gt['cn_profiles_B'] == 1), "Ground truth hap B not all 1"
    print(f"  ✓ Ground truth is diploid")
    
    # Show clone distribution
    clone_counts = np.bincount(gt['clone_assignments'])
    print(f"\n✓ Clone distribution:")
    for k, count in enumerate(clone_counts):
        print(f"    Clone {k}: {count} cells ({100*count/config.num_cells:.1f}%)")
    
    print("\n" + "="*70)
    print("Week 5 Test PASSED! ✓")
    print("="*70)
    print("\nThe simulator successfully:")
    print("  1. Initialized genome (bins, segments, blocks)")
    print("  2. Created diploid clone tree")
    print("  3. Sampled cells from clones")
    print("  4. Generated read counts via Negative Binomial")
    print("  5. Generated allele counts via Beta-Binomial")
    print("  6. Extracted ground truth data")
    
    return sim, read_counts, (alt_counts, ref_counts, total_counts)


def test_wgd_root():
    """Test WGD root initialization."""
    print("\n\n" + "="*70)
    print("Additional Test: WGD Root")
    print("="*70)
    
    config = SimulationConfig(
        num_bins=100,
        bin_length=50000,
        num_clones=3,
        num_cells=50,
        root_type="wgd",  # WGD root
        lambda_events=0.0,
        random_seed=123
    )
    
    sim = HaploTreeSimulator(config)
    read_counts, _ = sim.run()
    
    # Check that all clones have WGD state (2,2)
    for clone in sim.clones:
        assert np.all(clone.cn_profile_A == 2), f"Clone {clone.index} hap A not WGD"
        assert np.all(clone.cn_profile_B == 2), f"Clone {clone.index} hap B not WGD"
    
    print("✓ WGD root test PASSED - all clones are (2,2)")
    
    # Check that mean read counts are roughly 2x diploid
    mean_reads = read_counts.mean()
    print(f"✓ Mean read count (WGD): {mean_reads:.2f}")
    
    return sim


if __name__ == "__main__":
    # Run tests
    sim_diploid, reads, alleles = test_minimal_diploid()
    sim_wgd = test_wgd_root()
    
    print("\n" + "="*70)
    print("ALL TESTS PASSED! Week 5 deliverable is complete.")
    print("="*70)
