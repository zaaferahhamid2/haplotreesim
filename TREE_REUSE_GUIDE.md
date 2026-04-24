# Tree Structure Reuse Guide

## Overview

The simulator can save and reuse tree structures, enabling controlled parameter comparisons where tree topology remains fixed while other parameters vary.

## Use Cases

1. **WGD Probability Comparison**: Fix tree, vary `prob_wgd` from 0 to 0.6
2. **Event Rate Studies**: Fix tree, vary `lambda_events`  
3. **Coverage Analysis**: Fix tree and events, vary `mean_library_size`
4. **Doublet Sensitivity**: Fix tree, vary `prob_doublet`

## Example: Varying WGD While Fixing Tree
```python
from haplotreesim import SimulationConfig, HaploTreeSimulator

# Step 1: Generate base tree
base_config = SimulationConfig(
    chromosome='chr1',
    num_clones=5,
    random_seed=42  # Controls tree structure
)

sim = HaploTreeSimulator(base_config)
sim._initialize_genome()
sim._generate_clone_tree()

# Save tree
sim.export_tree_structure('my_tree.pkl')

# Step 2: Reuse tree with different parameters
for wgd_prob in [0.0, 0.3, 0.6]:
    config = SimulationConfig(
        chromosome='chr1',
        num_clones=5,
        prob_wgd=wgd_prob,
        random_seed=100 + int(wgd_prob*10)  # Different seed for events
    )
    
    sim = HaploTreeSimulator(config)
    sim._initialize_genome()
    sim.import_tree_structure('my_tree.pkl')
    sim._initialize_clones_from_imported_tree()
    sim._sample_cells()
    sim._generate_read_counts()
    sim._generate_allele_counts()
    
    print(f"WGD prob={wgd_prob}: ploidy={sim.clones[2].ploidy:.2f}")
```

## File Formats

- **`.pkl`**: Python pickle format (recommended, preserves types)
- **`.json`**: JSON format (human-readable, portable)

## What Gets Saved

- Node IDs and parent relationships
- Branch lengths  
- Clone proportions
- Leaf assignments

## What Gets Regenerated

- CNA events (controlled by `random_seed` and `lambda_events`)
- Observations (controlled by coverage parameters)
- WGD placement (controlled by `prob_wgd`)

## Important Notes

1. **Tree random seed**: Set `random_seed` when generating BASE tree
2. **Event random seed**: Change `random_seed` when reusing tree to get different events
3. **Same num_clones**: Imported tree must match config's `num_clones`

---

**See**: `examples/tree_reuse_example.py` for full working example
