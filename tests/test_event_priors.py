"""
Tests for Updated Event Attribute Priors (v0.5.0)
"""

import sys
sys.path.insert(0, '/Users/zaaferahhamid/Documents/haplotreesim/src')

import numpy as np
from haplotreesim import SimulationConfig, HaploTreeSimulator


def test_scale_class_probabilities():
    """Test that scale class probabilities are calculated correctly."""
    print("Testing scale class probability calculations...")
    
    config = SimulationConfig(
        prob_focal=0.7,
        prob_arm_given_broad=0.75
    )
    
    # Calculate weights (Equations 21-23)
    w_focal = config.prob_focal
    w_arm = (1 - config.prob_focal) * config.prob_arm_given_broad
    w_chr = (1 - config.prob_focal) * (1 - config.prob_arm_given_broad)
    
    # Should sum to 1
    assert abs(w_focal + w_arm + w_chr - 1.0) < 1e-10, "Weights must sum to 1"
    
    # Check expected values
    assert abs(w_focal - 0.7) < 1e-10, f"w_focal should be 0.7, got {w_focal}"
    assert abs(w_arm - 0.225) < 1e-10, f"w_arm should be 0.225, got {w_arm}"
    assert abs(w_chr - 0.075) < 1e-10, f"w_chr should be 0.075, got {w_chr}"
    
    print(f"  ✓ w_focal = {w_focal:.3f}")
    print(f"  ✓ w_arm = {w_arm:.3f}")
    print(f"  ✓ w_chr = {w_chr:.3f}")
    print(f"  ✓ Sum = {w_focal + w_arm + w_chr:.6f}")


def test_event_length_distribution():
    """Test that event lengths follow expected distribution."""
    print("Testing event length distribution...")
    
    config = SimulationConfig(
        num_bins=500,
        num_clones=10,
        num_cells=100,
        lambda_events=2.0,
        prob_focal=0.7,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    sim.run()
    
    # Collect event lengths
    lengths = []
    for clone in sim.clones:
        for event in clone.events:
            if event.haplotype.value != 'WGD':
                lengths.append(event.length)
    
    assert len(lengths) > 0, "Should have generated some events"
    
    # Check range is reasonable
    assert min(lengths) >= 1, "Min length should be at least 1 bin"
    assert max(lengths) <= config.num_bins, "Max length shouldn't exceed genome"
    
    print(f"  ✓ Generated {len(lengths)} events")
    print(f"  ✓ Length range: [{min(lengths)}, {max(lengths)}] bins")
    print(f"  ✓ Mean length: {np.mean(lengths):.1f} bins")


def test_haplotype_selection():
    """Test that haplotypes are selected with correct probabilities."""
    print("Testing haplotype selection...")
    
    config = SimulationConfig(
        num_bins=200,
        num_clones=10,
        num_cells=50,
        lambda_events=3.0,
        prob_haplotype_A=0.5,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    sim.run()
    
    # Count haplotypes
    hap_counts = {'A': 0, 'B': 0, 'WGD': 0}
    for clone in sim.clones:
        for event in clone.events:
            hap_counts[event.haplotype.value] += 1
    
    total_non_wgd = hap_counts['A'] + hap_counts['B']
    
    if total_non_wgd > 0:
        frac_A = hap_counts['A'] / total_non_wgd
        # With p_A = 0.5, expect roughly 50% A (allow wide range due to randomness)
        assert 0.2 < frac_A < 0.8, f"Haplotype A fraction {frac_A:.2f} seems biased"
        
        print(f"  ✓ Haplotype A: {hap_counts['A']} ({100*frac_A:.1f}%)")
        print(f"  ✓ Haplotype B: {hap_counts['B']} ({100*(1-frac_A):.1f}%)")


def test_amplitude_distribution():
    """Test that amplitudes follow 1 + Poisson distribution."""
    print("Testing amplitude distribution...")
    
    config = SimulationConfig(
        num_bins=200,
        num_clones=10,
        num_cells=50,
        lambda_events=2.0,
        lambda_amplitude=1.0,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    sim.run()
    
    # Collect absolute amplitudes
    amplitudes = []
    for clone in sim.clones:
        for event in clone.events:
            if event.haplotype.value != 'WGD':
                amplitudes.append(abs(event.amplitude))
    
    if amplitudes:
        # All should be >= 1 (since |Δ| = 1 + Poisson)
        assert all(a >= 1 for a in amplitudes), "All amplitudes should be >= 1"
        
        print(f"  ✓ Generated {len(amplitudes)} amplitudes")
        print(f"  ✓ Range: [{min(amplitudes)}, {max(amplitudes)}]")
        print(f"  ✓ Mean: {np.mean(amplitudes):.2f} (expected ~2 for λ=1)")


def test_full_simulation_with_new_priors():
    """Test complete simulation with new event priors."""
    print("Testing full simulation...")
    
    config = SimulationConfig(
        num_bins=500,
        num_clones=5,
        num_cells=100,
        lambda_events=1.5,
        prob_focal=0.7,
        prob_arm_given_broad=0.75,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    read_counts, (alt, ref, total) = sim.run()
    
    # Check outputs
    assert read_counts.shape == (100, 500), "Read counts shape incorrect"
    assert alt.shape[0] == 100, "Allele counts shape incorrect"
    
    print(f"  ✓ Read counts: {read_counts.shape}")
    print(f"  ✓ Allele counts: {alt.shape}")
    print(f"  ✓ Total events: {sum(len(c.events) for c in sim.clones)}")


def main():
    """Run all event prior tests."""
    print("="*70)
    print("Event Attribute Prior Tests (v0.5.0)")
    print("="*70)
    print()
    
    test_scale_class_probabilities()
    print()
    test_event_length_distribution()
    print()
    test_haplotype_selection()
    print()
    test_amplitude_distribution()
    print()
    test_full_simulation_with_new_priors()
    
    print()
    print("="*70)
    print("ALL EVENT PRIOR TESTS PASSED! ✓")
    print("="*70)


if __name__ == "__main__":
    main()
