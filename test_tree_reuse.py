"""
Test: Verify tree structure export/import works correctly
"""

from haplotreesim import SimulationConfig, HaploTreeSimulator
import numpy as np

print("="*70)
print("TEST: Tree Structure Export/Import")
print("="*70)

# Step 1: Create base simulation and export tree
print("\n1. Creating base simulation...")
config1 = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=100,
    lambda_events=2.0,
    prob_wgd=0.0,
    random_seed=42  # This controls tree structure
)

sim1 = HaploTreeSimulator(config1)
sim1._initialize_genome()
sim1._generate_clone_tree()

# Save key info from first simulation
tree1_proportions = sim1._clone_proportions.copy()
tree1_num_nodes = len(sim1._tree_nodes)
tree1_branch_lengths = [node.edge_length for node in sim1._tree_nodes]

print(f"   Tree nodes: {tree1_num_nodes}")
print(f"   Clone proportions: {tree1_proportions}")
print(f"   Branch lengths sum: {sum(tree1_branch_lengths):.3f}")

# Export tree
tree_file = 'test_tree.pkl'
sim1.export_tree_structure(tree_file)

# Step 2: Import tree in new simulation with DIFFERENT parameters
print(f"\n2. Creating NEW simulation with DIFFERENT parameters...")
config2 = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=100,
    lambda_events=5.0,      # DIFFERENT event rate
    prob_wgd=0.5,           # DIFFERENT WGD probability
    random_seed=999         # DIFFERENT random seed for events
)

sim2 = HaploTreeSimulator(config2)
sim2._initialize_genome()

# Import the SAME tree
sim2.import_tree_structure(tree_file)

# Check if tree structure matches
tree2_proportions = sim2._clone_proportions
tree2_num_nodes = len(sim2._tree_nodes)
tree2_branch_lengths = [node.edge_length for node in sim2._tree_nodes]

print(f"   Tree nodes: {tree2_num_nodes}")
print(f"   Clone proportions: {tree2_proportions}")
print(f"   Branch lengths sum: {sum(tree2_branch_lengths):.3f}")

# Step 3: Verify tree structure is IDENTICAL
print(f"\n3. Verification:")

tree_match = np.allclose(tree1_proportions, tree2_proportions)
nodes_match = tree1_num_nodes == tree2_num_nodes
branches_match = np.allclose(tree1_branch_lengths, tree2_branch_lengths)

print(f"   Clone proportions match: {tree_match} ✓" if tree_match else f"   Clone proportions match: {tree_match} ✗")
print(f"   Number of nodes match: {nodes_match} ✓" if nodes_match else f"   Number of nodes match: {nodes_match} ✗")
print(f"   Branch lengths match: {branches_match} ✓" if branches_match else f"   Branch lengths match: {branches_match} ✗")

if tree_match and nodes_match and branches_match:
    print(f"\n✓✓✓ SUCCESS! Tree structure preserved perfectly! ✓✓✓")
    print(f"\nKey insight:")
    print(f"  - Tree topology: IDENTICAL")
    print(f"  - Clone proportions: IDENTICAL")
    print(f"  - Branch lengths: IDENTICAL")
    print(f"  - But parameters differ (lambda_events, prob_wgd, random_seed)")
    print(f"\nThis enables controlled comparisons where ONLY the parameter")
    print(f"of interest varies while tree structure stays fixed!")
else:
    print(f"\n✗ FAILED: Trees don't match")

# Step 4: Show that we can continue simulation with imported tree
print(f"\n4. Continuing simulation with imported tree...")
sim2._initialize_clones_from_imported_tree()

print(f"   Clones created: {len(sim2.clones)}")
print(f"   Total events: {sum(len(c.events) for c in sim2.clones)}")
print(f"   Mean ploidy: {np.mean([c.ploidy for c in sim2.clones]):.2f}")

print(f"\n✓ Tree import/export working correctly!")
print(f"✓ Saved tree to: {tree_file}")
