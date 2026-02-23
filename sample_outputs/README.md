# Sample Simulation Outputs (Week 7 Update)

Generated with HaploTreeSim v0.3.0 including segment boundary detection.

## Configuration

- Genome bins: 500
- Number of clones: 5  
- Number of cells: 100
- CNA events: 13 total
- **Segments detected: 18** (from 19 breakpoints)
- Random seed: 42

## Files

### Visualizations

**`cn_profiles.png`** - Copy number profiles for all 5 clones
- Shows haplotype A (blue), haplotype B (orange), and total CN (black)
- Displays CNA events as changes in copy number
- Each clone subplot shows parent and event count

**`tree_structure.png`** - Clone tree diagram
- Green node = root clone
- Blue nodes = derived clones
- Numbers show event counts per clone

### Data Matrices (CSV)

**`read_counts.csv`** - 100 cells × 500 bins
- Sequencing read depth per cell per genomic bin
- Integer counts from Negative Binomial model

**`alternate_counts.csv`** - 100 cells × 18 segments
- Alternate allele counts per segment (Week 7: now per-segment!)
- From Beta-Binomial model with overdispersion

**`reference_counts.csv`** - 100 cells × 18 segments  
- Reference allele counts per segment

### Tree and Events (JSON)

**`tree_structure.json`**
- Clone tree topology (parent-child relationships)
- CNA events for each clone (position, haplotype, amplitude)
- Segment count: 18 segments

### Per-Clone CN Profiles (CSV)

**`clone_0_cn_profile.csv`** through **`clone_4_cn_profile.csv`**
- 500 bins × 3 columns (Haplotype_A, Haplotype_B, Total_CN)
- Integer copy numbers per genomic bin

## Week 7 Update

The key improvement in Week 7 is **segment boundary detection**:
- Previous: 1 segment covering entire genome
- Now: 18 segments based on CNA breakpoints
- Allele counts are now **per-segment** for finer resolution

## Loading Data
```python
import numpy as np
import json

# Load matrices
read_counts = np.loadtxt('read_counts.csv', delimiter=',')  # (100, 500)
alt_counts = np.loadtxt('alternate_counts.csv', delimiter=',')  # (100, 18)

# Load tree
with open('tree_structure.json') as f:
    tree = json.load(f)
    print(f"Segments: {tree['num_segments']}")  # 18
```

## Viewing Plots

The PNG files show:
1. CN profiles clearly show gains (CN > 2) and losses (CN < 2)
2. Different clones have different CN patterns (heterogeneity)
3. Haplotype-specific events visible (A vs B differences)
4. Tree structure shows all clones derived from root
