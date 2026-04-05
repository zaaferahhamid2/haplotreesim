# Week 8: Real Chromosomes and Global Segmentation (v0.6.0)

Implements real chromosome support using the hg38 reference genome based on professor feedback.

## Summary

Previously, the simulator used abstract genome coordinates. It now supports real human chromosomes (hg38) with accurate chromosome lengths, centromere positions, and configurable bin widths.

## Key Features

- Supports chromosomes chr1–22, chrX, chrY  
- Uses real chromosome lengths and centromeres  
- Configurable bin width (default 500 kb)  
- Automatic bin computation  
- Global segmentation across clones  

## Usage

```python
from haplotreesim import SimulationConfig, HaploTreeSimulator

config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=100
)

sim = HaploTreeSimulator(config)
read_counts, alleles = sim.run()
