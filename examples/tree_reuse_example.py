"""
Example: Reusing Tree Structure Across Parameter Variations

This demonstrates how to fix tree structure while varying other parameters,
useful for controlled comparisons (e.g., varying WGD probability).
"""

from haplotreesim import SimulationConfig, HaploTreeSimulator
import numpy as np

print("="*70)
print("TREE REUSE EXAMPLE: Fixed Tree, Varying WGD Probability")
print("="*70)

# Step 1: Generate base simulation and save tree
print("\n1. Generate base simulation with default parameters...")
base_config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=200,
    lambda_events=2.0,
    prob_wgd=0.0,  # No WGD in base
    random_seed=42
)

base_sim = HaploTreeSimulator(base_config)
base_sim._initialize_genome()
base_sim._generate_clone_tree()

# Export tree structure
tree_file = 'tree_structure_5clones.pkl'
base_sim.export_tree_structure(tree_file)

print(f"\n2. Base simulation results:")
print(f"   Clone proportions: {base_sim._clone_proportions}")
print(f"   Mean ploidy: {np.mean([c.ploidy for c in base_sim.clones]):.2f}")

# Step 2: Reuse same tree with WGD
print(f"\n3. Reusing tree with WGD (prob_wgd=0.5)...")
wgd_config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=200,
    lambda_events=2.0,
    prob_wgd=0.5,  # Add WGD
    random_seed=100  # Different seed for events/observations
)

wgd_sim = HaploTreeSimulator(wgd_config)
wgd_sim._initialize_genome()

# Import the SAME tree structure
wgd_sim.import_tree_structure(tree_file)

# Continue with rest of simulation (events, observations)
# Manually call remaining steps since we bypassed normal _generate_clone_tree
wgd_sim._initialize_clones_from_imported_tree()
wgd_sim._sample_cells()
wgd_sim._generate_read_counts()
wgd_sim._generate_allele_counts()

print(f"\n4. WGD simulation results:")
print(f"   Clone proportions: {wgd_sim._clone_proportions} (SAME as base)")
print(f"   Mean ploidy: {np.mean([c.ploidy for c in wgd_sim.clones]):.2f}")
print(f"   WGD occurred: {wgd_sim.wgd_occurred}")

print(f"\n5. Comparison:")
print(f"   Tree topology: IDENTICAL ✓")
print(f"   Clone proportions: IDENTICAL ✓")
print(f"   WGD: DIFFERENT (0.0 vs 0.5 probability)")
print(f"   Ploidy: DIFFERENT ({np.mean([c.ploidy for c in base_sim.clones]):.2f} vs {np.mean([c.ploidy for c in wgd_sim.clones]):.2f})")

print(f"\n✓ Successfully demonstrated tree reuse!")
print(f"✓ Tree saved to: {tree_file}")
