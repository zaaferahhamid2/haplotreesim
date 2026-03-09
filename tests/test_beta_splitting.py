"""
Tests for Beta-Splitting Tree Implementation
"""

import sys
sys.path.insert(0, '/Users/zaaferahhamid/Documents/haplotreesim/src')

import numpy as np
from haplotreesim import (
    SimulationConfig,
    HaploTreeSimulator,
    BetaSplittingTreeBuilder
)


def test_beta_tree_structure():
    """Test that Beta-splitting creates correct tree structure."""
    print("Testing Beta-splitting tree structure...")
    
    rng = np.random.default_rng(42)
    builder = BetaSplittingTreeBuilder(
        rng=rng,
        num_clones=5,
        alpha_tree=0.5,
        beta_tree=0.3
    )
    
    nodes, parent_map = builder.build_tree()
    
    # Should have 2K - 1 nodes
    assert len(nodes) == 9, f"Expected 9 nodes, got {len(nodes)}"
    
    # Should have exactly K leaves
    leaves = [n for n in nodes if n.is_leaf]
    assert len(leaves) == 5, f"Expected 5 leaves, got {len(leaves)}"
    
    # Root should have parent_id = -1
    root = nodes[0]
    assert root.parent_id == -1, "Root should have parent_id = -1"
    assert root.is_leaf == False, "Root should not be a leaf"
    
    print(f"  ✓ Tree has {len(nodes)} nodes ({len(leaves)} leaves)")


def test_clone_proportions():
    """Test that clone proportions sum to 1 and are non-uniform."""
    print("\nTesting clone proportions...")
    
    rng = np.random.default_rng(42)
    builder = BetaSplittingTreeBuilder(
        rng=rng,
        num_clones=5,
        alpha_tree=0.5,
        beta_tree=0.3
    )
    
    nodes, _ = builder.build_tree()
    proportions = builder.get_clone_proportions(nodes)
    
    # Should sum to 1
    assert abs(proportions.sum() - 1.0) < 1e-10, \
        f"Proportions should sum to 1, got {proportions.sum()}"
    
    # Should be length K
    assert len(proportions) == 5, f"Expected 5 proportions, got {len(proportions)}"
    
    # With beta=0.3, should be imbalanced (not all equal)
    # Check that not all proportions are approximately equal
    std_dev = proportions.std()
    assert std_dev > 0.05, "Proportions should be imbalanced with beta=0.3"
    
    print(f"  ✓ Proportions sum to {proportions.sum():.6f}")
    print(f"  ✓ Proportions: {proportions}")
    print(f"  ✓ Std dev: {std_dev:.4f} (imbalanced)")


def test_branch_lengths():
    """Test that branch lengths are normalized."""
    print("\nTesting branch lengths...")
    
    rng = np.random.default_rng(42)
    builder = BetaSplittingTreeBuilder(
        rng=rng,
        num_clones=5,
        alpha_tree=0.5,
        beta_tree=0.3
    )
    
    nodes, _ = builder.build_tree()
    
    # Sum all branch lengths (excluding root's edge to itself)
    total_length = sum(node.edge_length for node in nodes)
    
    # Should sum to 1 (normalized)
    assert abs(total_length - 1.0) < 1e-10, \
        f"Branch lengths should sum to 1, got {total_length}"
    
    print(f"  ✓ Total branch length: {total_length:.6f}")


def test_events_tied_to_branch_length():
    """Test that events are tied to branch length (Equation 18)."""
    print("\nTesting events tied to branch length...")
    
    config = SimulationConfig(
        num_bins=200,
        num_clones=4,
        num_cells=50,
        lambda_events=2.0,  # High rate
        alpha_tree=0.5,
        beta_tree=0.3,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    sim.run()
    
    # Check that events vary by edge
    event_counts = [len(clone.events) for clone in sim.clones if not clone.is_root]
    
    # Should have variation in event counts
    assert len(set(event_counts)) > 1, "Event counts should vary across edges"
    
    # Total events should be reasonable
    total_events = sum(event_counts)
    assert total_events > 0, "Should have generated some events"
    
    print(f"  ✓ Event counts per edge: {event_counts}")
    print(f"  ✓ Total events: {total_events}")


def test_cell_sampling_from_proportions():
    """Test that cells are sampled according to Beta proportions."""
    print("\nTesting cell sampling from proportions...")
    
    config = SimulationConfig(
        num_bins=100,
        num_clones=5,
        num_cells=1000,  # Large sample for statistical test
        lambda_events=1.0,
        alpha_tree=0.5,
        beta_tree=0.3,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    sim.run()
    
    # Get true proportions
    true_props = sim._clone_proportions
    
    # Get observed proportions from cell assignments
    gt = sim.get_ground_truth()
    assignments = gt['clone_assignments']
    
    # Count assignments to leaf clones
    leaf_indices = sim._leaf_clone_indices
    observed_counts = np.zeros(len(leaf_indices))
    
    for i, leaf_idx in enumerate(leaf_indices):
        observed_counts[i] = (assignments == leaf_idx).sum()
    
    observed_props = observed_counts / config.num_cells
    
    # Check that observed is close to expected (within sampling variance)
    for i in range(len(true_props)):
        diff = abs(observed_props[i] - true_props[i])
        # With 1000 cells, should be within ~5% (generous bound)
        assert diff < 0.10, \
            f"Observed proportion {observed_props[i]:.3f} differs too much from expected {true_props[i]:.3f}"
    
    print(f"  ✓ Expected proportions: {true_props}")
    print(f"  ✓ Observed proportions: {observed_props}")
    print(f"  ✓ Max difference: {max(abs(observed_props - true_props)):.4f}")


def test_full_simulation_with_beta_tree():
    """Test complete simulation with Beta-splitting tree."""
    print("\nTesting full simulation...")
    
    config = SimulationConfig(
        num_bins=500,
        num_clones=5,
        num_cells=100,
        lambda_events=1.5,
        alpha_tree=0.5,
        beta_tree=0.3,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    read_counts, (alt, ref, total) = sim.run()
    
    # Check outputs
    assert read_counts.shape == (100, 500), "Read counts shape incorrect"
    assert alt.shape[0] == 100, "Allele counts shape incorrect"
    
    # Check that we have the tree structure
    assert len(sim.clones) == 9, "Should have 9 total clones (2*5-1)"
    assert len(sim._leaf_clone_indices) == 5, "Should have 5 leaf clones"
    
    print(f"  ✓ Read counts shape: {read_counts.shape}")
    print(f"  ✓ Allele counts shape: {alt.shape}")
    print(f"  ✓ Total clones: {len(sim.clones)} (leaves: {len(sim._leaf_clone_indices)})")


def main():
    """Run all Beta-splitting tests."""
    print("="*70)
    print("Beta-Splitting Tree Implementation Tests")
    print("="*70)
    
    test_beta_tree_structure()
    test_clone_proportions()
    test_branch_lengths()
    test_events_tied_to_branch_length()
    test_cell_sampling_from_proportions()
    test_full_simulation_with_beta_tree()
    
    print("\n" + "="*70)
    print("ALL BETA-SPLITTING TESTS PASSED! ✓")
    print("="*70)


if __name__ == "__main__":
    main()
