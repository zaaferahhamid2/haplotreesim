# HaploTreeSim

**A Controlled, End-to-End Benchmark with Haplotype-Resolved CNA and CNA-Phylogeny Ground Truth for Low-Pass scDNA-seq**

## Week 5 Deliverable: Repository + Configuration Skeleton

This repository contains the initial implementation of HaploTreeSim, a simulator for generating synthetic single-cell DNA sequencing (scDNA-seq) data with ground-truth haplotype-specific copy number alterations (CNAs) and clone phylogenies.

### Current Status (Week 5)

The minimal simulator is now functional and produces diploid output (no CNAs yet). The following components are implemented:

 yes - **Repository Structure**: Clean Python package with `src/`, `tests/`, and `configs/` directories  
 yes - **Data Models**: Complete class definitions for bins, segments, clones, events, and cells  
 yes - **Simulator Core**: Basic simulator that outputs diploid read-count and allele-count matrices  
 yes - **Configuration System**: Flexible configuration using dataclasses  
 yes - **Test Suite**: Verification tests for Week 5 deliverable

### Repository Structure

```
haplotreesim/
├── src/
│   └── haplotreesim/
│       ├── __init__.py           # Package initialization
│       ├── data_models.py        # Core data structures
│       └── simulator.py          # Main simulator class
├── tests/
│   └── test_week5.py             # Week 5 verification tests
├── configs/
│   └── example_configs.py        # Example configurations
├── docs/                         # Documentation (future)
├── setup.py                      # Package installation
├── README.md                     # This file
└── requirements.txt              # Python dependencies
```

### Data Models

The following data structures are implemented (see `src/haplotreesim/data_models.py`):

- **`Bin`**: Fixed-width genomic regions with GC/mappability annotations
- **`Segment`**: Contiguous intervals of bins with uniform copy number
- **`HaplotypeBlock`**: Groups of segments with shared phasing
- **`CNAEvent`**: Copy number alteration events (gains, losses, WGD)
- **`Clone`**: Tree nodes with haplotype-specific CN profiles
- **`Cell`**: Individual cells with clone assignments and observations
- **`SimulationConfig`**: Complete parameter configuration

### Installation

```bash
cd haplotreesim
pip install -e .
```

Or install with development dependencies:

```bash
pip install -e ".[dev]"
```

### Quick Start

```python
from haplotreesim import SimulationConfig, HaploTreeSimulator

# Create configuration
config = SimulationConfig(
    num_bins=1000,
    num_clones=5,
    num_cells=200,
    root_type="diploid",
    random_seed=42
)

# Initialize and run simulator
sim = HaploTreeSimulator(config)
read_counts, (alt_counts, ref_counts, total_counts) = sim.run()

# Get ground truth
ground_truth = sim.get_ground_truth()
```

### Running Tests

```bash
cd haplotreesim
python tests/test_week5.py
```

Expected output:
```
======================================================================
Week 5 Test: Minimal Diploid Simulator
======================================================================

Initializing genome...
Generating clone tree with 3 clones...
Sampling 50 cells...
Generating read-depth observations...
Generating allelic observations...
Simulation complete!

======================================================================
Verification
======================================================================
✓ Read counts shape: (50, 100)
✓ Allele counts shape: (50, 1)
✓ All 3 clones are diploid (1,1)
✓ Mean read count per bin: 25.34
✓ Total allelic reads: 75234

Week 5 Test PASSED! ✓
```

### Current Functionality (Week 5)

The simulator currently supports:

1. **Genome Initialization**
   - Fixed-width bins with configurable size
   - Single segment covering entire genome (breakpoints coming in Week 6)
   - Single haplotype block (phase switches coming later)

2. **Clone Tree Generation**
   - Root clone initialization (diploid or WGD)
   - Multiple clones (all copies of root for now)
   - Tree structure and CNA events coming in Week 6

3. **Cell Sampling**
   - Uniform sampling from clones
   - Log-normal library size and allelic coverage factors
   - Normal cells and doublets coming later

4. **Observation Models**
   - **Read counts**: Negative Binomial model with overdispersion
   - **Allele counts**: Beta-Binomial model with overdispersion
   - Realistic noise matching low-pass scDNA-seq

5. **Ground Truth Extraction**
   - Clone assignments (z_n)
   - Haplotype-specific CN profiles (c^(A), c^(B))
   - Clone tree structure
   - Event lists (empty for Week 5)

### Configuration Parameters

Key parameters in `SimulationConfig`:

```python
# Genome
num_bins: int = 10000              # Number of genomic bins
bin_length: int = 50000            # Bin size in bp (50kb default)

# Clones
num_clones: int = 5                # Number of clones
root_type: str = "diploid"         # "diploid" or "wgd"

# Cells
num_cells: int = 200               # Number of cells to simulate

# Read-depth model
mean_library_size: float = 100.0   # Mean coverage factor α_n
theta_x: float = 10.0              # Overdispersion parameter

# Allelic model
snp_density: float = 0.001         # Heterozygous SNPs per bp
mean_allelic_coverage: float = 50.0  # Mean allelic depth β_n
nu_a: float = 20.0                 # Beta-Binomial concentration
```

### Next Steps (Week 6)

The following features will be implemented in Week 6:

- [ ] Segment boundary detection from CNA breakpoints
- [ ] CNA event generation (gains, losses)
- [ ] Event application to clone CN profiles
- [ ] Tree structure generation
- [ ] Non-uniform clone frequencies

### Example Configurations

See `configs/example_configs.py` for three example configurations:

1. **`minimal_diploid`**: Simple 100-cell, 3-clone diploid simulation
2. **`low_pass_config`**: 200-cell simulation with low-pass coverage (0.01x)
3. **`wgd_root_config`**: WGD root simulation for tetraploid tumors

### References

This simulator implements the models described in the HaploTreeSim paper (sections 3.1-3.6).

### License

