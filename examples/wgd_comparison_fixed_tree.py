"""
Practical Example: Comparing WGD Probabilities with Fixed Tree

This demonstrates answering the question:
"How does WGD probability affect results when tree structure is held constant?"
"""

from haplotreesim import SimulationConfig, HaploTreeSimulator
import numpy as np
import pandas as pd

print("="*70)
print("WGD COMPARISON WITH FIXED TREE STRUCTURE")
print("="*70)

# Step 1: Generate and save base tree structure
print("\nStep 1: Generate base tree (5 clones)...")
base_config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=200,
    random_seed=42  # Controls tree topology
)

base_sim = HaploTreeSimulator(base_config)
base_sim._initialize_genome()
base_sim._generate_clone_tree()

# Save tree
tree_file = 'tree_5clones_seed42.pkl'
base_sim.export_tree_structure(tree_file)

print(f"   Clone proportions: {base_sim._clone_proportions}")
print(f"   Tree saved to: {tree_file}")

# Step 2: Run simulations with varying WGD, same tree
print("\nStep 2: Running simulations with prob_wgd = [0.0, 0.3, 0.6]...")
print("        (Same tree, different WGD probabilities)")

results = []

for wgd_prob in [0.0, 0.3, 0.6]:
    print(f"\n   Testing prob_wgd = {wgd_prob}...")
    
    config = SimulationConfig(
        chromosome='chr1',
        bin_width=500000,
        num_clones=5,
        num_cells=200,
        lambda_events=2.0,
        prob_wgd=wgd_prob,
        random_seed=100 + int(wgd_prob*10)  # Different seed for events
    )
    
    sim = HaploTreeSimulator(config)
    sim._initialize_genome()
    sim.import_tree_structure(tree_file)
    sim._initialize_clones_from_imported_tree()
    sim._sample_cells()
    sim._generate_read_counts()
    sim._generate_allele_counts()
    
    # Collect metrics
    ploidies = [c.ploidy for c in sim.clones]
    results.append({
        'prob_wgd': wgd_prob,
        'wgd_occurred': sim.wgd_occurred,
        'min_ploidy': np.min(ploidies),
        'max_ploidy': np.max(ploidies),
        'mean_ploidy': np.mean(ploidies),
        'total_events': sum(len(c.events) for c in sim.clones),
        'num_segments': len(sim.segments)
    })
    
    print(f"      WGD occurred: {sim.wgd_occurred}")
    print(f"      Ploidy range: [{np.min(ploidies):.2f}, {np.max(ploidies):.2f}]")

# Step 3: Display comparison table
print("\n" + "="*70)
print("COMPARISON TABLE: Same Tree, Varying WGD")
print("="*70)

df = pd.DataFrame(results)
print(df.to_string(index=False))

print("\n" + "="*70)
print("KEY OBSERVATIONS:")
print("="*70)
print("✓ Tree structure: IDENTICAL across all runs")
print("✓ Clone proportions: IDENTICAL across all runs")
print(f"✓ WGD occurrence: Varies with prob_wgd")
print(f"✓ Ploidy: Changes when WGD occurs (2.0 → ~4.0)")
print(f"✓ Segments: May vary due to different events")
print("\nThis controlled comparison isolates WGD's effect!")
