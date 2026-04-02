# Week 8: Real Chromosomes and Global Segmentation (v0.6.0)

Implements real chromosome support using hg38 reference genome per professor feedback.

## Summary

Previously, the simulator used abstract genome coordinates. Now it uses **real human chromosomes from hg38** with actual chromosome lengths, centromere positions, and configurable bin widths.

## Key Features

### 1. Real Chromosome Data (hg38)

All 24 human chromosomes supported:
- **Autosomes**: chr1-22
- **Sex chromosomes**: chrX, chrY
- **Chromosome 1 example**: 248,956,422 bp (249 Mb)

Source: UCSC Genome Browser hg38/GRCh38 reference

### 2. Configurable Bin Width
```python
config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,  # 500 kb bins (configurable)
)
```

**Automatic computation:**
- `num_bins = ⌈chromosome_length / bin_width⌉`
- Chr1 at 500 kb → 498 bins
- Chr1 at 1 Mb → 249 bins

### 3. Real Chromosome Arms

**Centromere positions** from hg38:
- Chr1 p-arm: 0 - 122 Mb (244 bins)
- Chr1 centromere: 122 - 125 Mb
- Chr1 q-arm: 125 - 249 Mb (247 bins)

**Arm-level events** use real boundaries:
- Sample arm with probability ∝ arm length
- Cover entire p or q arm

### 4. Global Segmentation

**Section 3.6 implementation** (from Week 7):

Segments defined as maximal contiguous regions where **no extant clone** has a changepoint:
```
Breakpoint at (b, b+1) ⟺ ∃k: c_k,b^(A) ≠ c_k,b+1^(A) OR c_k,b^(B) ≠ c_k,b+1^(B)
```

**Example (chr1)**:
- 8 CNA events across 5 clones
- 15 breakpoints detected
- 14 segments created

## Usage

### Basic Simulation
```python
from haplotreesim import SimulationConfig, HaploTreeSimulator

config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=100,
    lambda_events=1.5
)

sim = HaploTreeSimulator(config)
read_counts, alleles = sim.run()

print(f"Simulated {config.chromosome}")
print(f"Bins: {config.num_bins}")
print(f"Segments: {len(sim.segments)}")
```

### Chromosome Information
```python
from haplotreesim import chromosome_data

# Get chromosome details
info = chromosome_data.describe_chromosome('chr1', bin_width=500000)

print(f"Length: {info['length_mb']:.1f} Mb")
print(f"Bins: {info['num_bins']}")
print(f"P arm: {info['p_arm_length_mb']:.1f} Mb")
print(f"Q arm: {info['q_arm_length_mb']:.1f} Mb")
```

### Visualize CN Profiles
```bash
python create_chr_plots.py
```

Creates `sample_outputs/chr1_cn_profiles.png` showing:
- All clones' CN profiles
- Genomic coordinates in Mb
- Centromere marked
- Haplotype-specific CN

## Implementation Details

### New Module: chromosome_data.py

**Functions**:
- `get_chromosome_length(chr)` - Length in bp
- `create_bins_for_chromosome(chr, bin_width)` - Compute bins
- `get_arm_boundaries(chr)` - Centromere position
- `describe_chromosome(chr, bin_width)` - Full info

**Data**: `HG38_CHROMOSOME_LENGTHS` dictionary with all 24 chromosomes

### Updated Components

**SimulationConfig**:
```python
chromosome: str = 'chr1'      # Which chromosome
bin_width: int = 500000       # Bin size in bp
num_bins: Optional[int]       # Auto-computed
```

**EventGenerator**:
- Uses real arm lengths for focal event max
- Samples arms with probability ∝ arm length
- Centromere-aware arm events

**Simulator**:
- Passes chromosome info to event generator
- Bins labeled with actual genomic coordinates

## Verification

**Test different chromosomes**:
```python
# Small chromosome
config = SimulationConfig(chromosome='chr21', bin_width=500000)
# → 94 bins (46.7 Mb)

# Large chromosome
config = SimulationConfig(chromosome='chr1', bin_width=500000)
# → 498 bins (249 Mb)

# Different bin width
config = SimulationConfig(chromosome='chr1', bin_width=1000000)
# → 249 bins (1 Mb each)
```

## Week 8 Deliverable

✅ **Plot showing CN profiles per clone across the genome**

File: `sample_outputs/chr1_cn_profiles.png`

Shows:
- 9 clones (5 leaves + 4 internal nodes)
- Haplotype A, B, and total CN
- Real genomic coordinates (0-249 Mb)
- Centromere position marked
- CNA events visible as CN changes

## Next Steps

- Week 9-10: Non-uniform clone frequencies, mirrored subclones
- Extend to multiple chromosomes
- Whole-genome simulations

## Files Changed

- `src/haplotreesim/chromosome_data.py` - **NEW**: hg38 data
- `src/haplotreesim/data_models.py` - Added chromosome params
- `src/haplotreesim/event_generator.py` - Real arm boundaries
- `src/haplotreesim/simulator.py` - Pass chromosome info
- `create_chr_plots.py` - **NEW**: Visualization script
- `sample_outputs/chr1_cn_profiles.png` - **NEW**: Week 8 deliverable
