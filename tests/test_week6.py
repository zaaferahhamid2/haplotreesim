"""
Week 6 Tests: CNA Events and Tree Structure
"""

import sys
sys.path.insert(0, '/Users/zaaferahhamid/Documents/haplotreesim/src')

import numpy as np
from haplotreesim import (
    SimulationConfig, 
    HaploTreeSimulator,
    EventGenerator,
    EventApplier,
    TreeBuilder,
    Clone,
    CNAEvent,
    Haplotype
)


def test_event_generator():
    """Test EventGenerator creates valid events."""
    print("Testing EventGenerator...")
    
    rng = np.random.default_rng(42)
    gen = EventGenerator(rng, num_bins=100, lambda_events=2.0)
    
    events = gen.generate_events_for_edge(allow_wgd=False)
    
    assert isinstance(events, list), "Events should be a list"
    assert all(isinstance(e, CNAEvent) for e in events), "All events should be CNAEvent"
    
    for event in events:
        assert 0 <= event.start_bin < 100, "Start bin out of range"
        assert 0 <= event.end_bin < 100, "End bin out of range"
        assert event.start_bin <= event.end_bin, "Start should be <= end"
        assert event.haplotype in [Haplotype.A, Haplotype.B], "Invalid haplotype"
        assert event.amplitude != 0, "Non-WGD events must have non-zero amplitude"
    
    print(f"  ✓ Generated {len(events)} valid events")


def test_event_applier():
    """Test EventApplier correctly modifies CN profiles."""
    print("\nTesting EventApplier...")
    
    # Create diploid clone
    cn_A = np.ones(100, dtype=int)
    cn_B = np.ones(100, dtype=int)
    
    clone = Clone(
        index=0,
        parent_index=None,
        cn_profile_A=cn_A,
        cn_profile_B=cn_B,
        events=[],
        is_root=True
    )
    
    # Create gain event on haplotype A
    event = CNAEvent(
        start_bin=10,
        end_bin=20,
        haplotype=Haplotype.A,
        amplitude=1  # Gain of 1
    )
    
    applier = EventApplier(max_copy_number=8)
    new_cn_A, new_cn_B = applier.apply_events(clone, [event])
    
    # Check that gain was applied
    assert new_cn_A[10] == 2, "Gain should increase CN to 2"
    assert new_cn_A[20] == 2, "Gain should affect end bin"
    assert new_cn_A[9] == 1, "Region before should be unchanged"
    assert new_cn_A[21] == 1, "Region after should be unchanged"
    assert np.all(new_cn_B == 1), "Haplotype B should be unchanged"
    
    print("  ✓ Events correctly modify CN profiles")


def test_wgd_event():
    """Test WGD event doubles both haplotypes."""
    print("\nTesting WGD event...")
    
    cn_A = np.ones(100, dtype=int)
    cn_B = np.ones(100, dtype=int)
    
    clone = Clone(
        index=0,
        parent_index=None,
        cn_profile_A=cn_A,
        cn_profile_B=cn_B,
        events=[],
        is_root=True
    )
    
    wgd_event = CNAEvent(
        start_bin=0,
        end_bin=99,
        haplotype=Haplotype.WGD,
        amplitude=0
    )
    
    applier = EventApplier(max_copy_number=8)
    new_cn_A, new_cn_B = applier.apply_events(clone, [wgd_event])
    
    assert np.all(new_cn_A == 2), "WGD should double haplotype A"
    assert np.all(new_cn_B == 2), "WGD should double haplotype B"
    
    print("  ✓ WGD correctly doubles both haplotypes")


def test_tree_builder():
    """Test TreeBuilder creates valid tree structures."""
    print("\nTesting TreeBuilder...")
    
    rng = np.random.default_rng(42)
    builder = TreeBuilder(rng, num_clones=5)
    
    # Test star tree
    star_tree = builder.build_tree(tree_type="star")
    assert star_tree[0] is None, "Root should have no parent"
    assert all(star_tree[i] == 0 for i in range(1, 5)), "All clones should be children of root"
    
    # Test linear tree
    linear_tree = builder.build_tree(tree_type="linear")
    assert linear_tree[0] is None, "Root should have no parent"
    assert linear_tree[1] == 0, "Clone 1 parent should be 0"
    assert linear_tree[2] == 1, "Clone 2 parent should be 1"
    
    print("  ✓ Tree structures are valid")


def test_week6_simulation():
    """Test complete Week 6 simulation with CNA events."""
    print("\nTesting Week 6 simulation...")
    
    config = SimulationConfig(
        num_bins=100,
        num_clones=5,
        num_cells=100,
        lambda_events=1.5,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    read_counts, alleles = sim.run()
    
    # Check clones have events
    total_events = sum(len(clone.events) for clone in sim.clones)
    assert total_events > 0, "Should have generated some events"
    
    # Check at least one clone is non-diploid
    non_diploid_clones = 0
    for clone in sim.clones:
        if not np.all(clone.total_cn() == 2):
            non_diploid_clones += 1
    
    assert non_diploid_clones > 0, "At least one clone should be non-diploid"
    
    print(f"  ✓ Generated {total_events} events across {len(sim.clones)} clones")
    print(f"  ✓ {non_diploid_clones} clones are non-diploid")


def main():
    """Run all Week 6 tests."""
    print("="*70)
    print("Week 6 Tests: CNA Events and Tree Structure")
    print("="*70)
    
    test_event_generator()
    test_event_applier()
    test_wgd_event()
    test_tree_builder()
    test_week6_simulation()
    
    print("\n" + "="*70)
    print("ALL WEEK 6 TESTS PASSED! ✓")
    print("="*70)


if __name__ == "__main__":
    main()
