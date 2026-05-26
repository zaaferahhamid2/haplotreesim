# HaploTreeSim
**A Controlled, End-to-End Benchmark with Haplotype-Resolved CNA and CNA-Phylogeny Ground Truth for Low-Pass scDNA-seq**

[![Version](https://img.shields.io/badge/version-1.3.0-blue.svg)](https://github.com/zaaferahhamid2/haplotreesim/releases/tag/v1.3.0)
[![Paper Compliant](https://img.shields.io/badge/paper-compliant-brightgreen.svg)](#)

HaploTreeSim is a simulator for generating synthetic single-cell DNA sequencing (scDNA-seq) datasets with ground-truth haplotype-specific copy number alterations (CNAs) and clone phylogenies. It is designed to benchmark ASCN callers, CNA-tree methods, and joint phylogeny pipelines under controlled conditions.

---

## Table of Contents

- [Installation](#installation)
- [Running End-to-End](#running-end-to-end)
- [Current Pipeline Updates](#current-pipeline-updates)
- [SEACON Pipeline](#seacon-pipeline)
- [Parameter Configuration](#parameter-configuration)
- [Input Format](#input-format)
- [Expected Output](#expected-output)
- [For Evaluators](#for-evaluators)
- [Key Features](#key-features)
- [Roadmap](#roadmap)

---

## Installation

```bash
git clone https://github.com/zaaferahhamid2/haplotreesim.git
cd haplotreesim
pip install -e .
```

**Requirements:** Python 3.9+, numpy, scipy, pandas

---

## Running End-to-End

### Option 1: Benchmark CLI (recommended)

Run the full pipeline — simulate data, run a method, compute all metrics — with one command:

```bash
python3 -m haplotreesim.benchmark --config configs/default_benchmark.json
```

This will:
1. Generate a synthetic scDNA-seq dataset from the config
2. Run the specified method(s) on the simulated data
3. Compute all metrics against ground truth
4. Save results to `benchmark_results/`

### Option 2: Python API

```python
from haplotreesim import SimulationConfig, HaploTreeSimulator

# 1. Configure the simulation
config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,       # 500 kb bins
    num_clones=5,           # K clones
    num_cells=200,          # N cells
    lambda_events=1.5,      # mean CNA events per edge
    prob_wgd=0.3,           # WGD probability
    prob_doublet=0.05,      # doublet rate
    random_seed=42
)

# 2. Run the simulator
sim = HaploTreeSimulator(config)
read_counts, allele_counts = sim.run()
# read_counts:   np.ndarray shape (N, J) — read depth per cell per bin
# allele_counts: tuple of (alt_counts, ref_counts, total_counts) per segment

# 3. Access ground truth
ground_truth = sim.get_ground_truth()
# See "Expected Output" section for full ground truth contents

for clone in sim.clones:
    print(f"Clone {clone.index}: ploidy={clone.ploidy:.2f}")
```

### Option 3: Example script

```bash
python3 examples/run_benchmark_example.py
```

---

## Current Pipeline Updates

### Pipeline Additions

- File-based external-tool workflow.
- Dataset simulator writes reusable on-disk outputs.
- Tool adapters convert saved datasets into method-specific inputs.
- Method runners execute external callers such as SEACON.
- Evaluators compare method outputs against saved simulator truth.

### Simulator Changes

- Supports single-chromosome, multi-chromosome, and whole-genome simulation.
- Exports bins, segments, cells, read-count matrices, and RDR matrices.
- Exports allele-count matrices and segment-level allele-count summaries.
- Exports clone-level haplotype CN truth and cell-segment HSCN truth.
- Exports event tables, metadata, and reusable tree structure.

### Event Model Changes

- CNA events are generated in order along clone-tree edges.
- Gains and losses are validated against the current haplotype-specific CN state.
- Focal, arm-level, and chromosome-level events support whole-genome coordinates.
- Focal event length now uses `focal_length_min`.
- Focal event maximum length now uses `focal_length_max_fraction`.
- The old hardcoded focal length range has been removed.

### Metrics / Evaluation Changes

- Evaluation reads saved datasets directly.
- Truth is loaded from exported files instead of regenerated from a seed.
- SEACON calls can be parsed from segment-style or flat outputs.
- Predictions are mapped back to simulator bins and true segments.
- Reports HSCN error and LOH precision/recall/F1.
- Reports recurrent breakpoint metrics and breakpoint sensitivity curves.
- Reports clone assignment metrics.

---

## SEACON Pipeline

The current SEACON workflow is file-based. First simulate a general HaploTreeSim dataset, then convert it into SEACON input files, run SEACON, and evaluate the SEACON calls against the saved simulator truth.

### 1. Simulate a dataset

```bash
python scripts/simulate_dataset.py \
  --output-dir examples/whole_genome_simulation \
  --whole-genome \
  --num-clones 5 \
  --num-cells 100 \
  --lambda-events 5 \
  --lambda-amplitude 1 \
  --random-seed 42 \
  --overwrite
```

This writes the general dataset: bins, segments, cells, read counts, RDR, allele counts, clone CN truth, cell-segment HSCN truth, events, metadata, and tree structure.

### 2. Convert to SEACON input

```bash
python scripts/convert_to_seacon.py \
  --input-dir examples/whole_genome_simulation \
  --overwrite
```

This creates `examples/whole_genome_simulation/seacon_input/` with `RDR.tsv`, `readcounts.tsv`, `filtered_regions.bed`, `cells.txt`, and `precomputed_baf.tsv`.

### 3. Run SEACON

```bash
python scripts/run_seacon.py \
  --input-dir examples/whole_genome_simulation/seacon_input \
  --output-dir examples/whole_genome_simulation/seacon_output \
  --upper-filter 20 \
  --tolerance 0.05 \
  --max-wgd 1 \
  --max-cn 10 \
  --num-processors 1 \
  --overwrite
```

Use `--seacon-bin /path/to/seacon` if the `seacon` executable is not on your `PATH`.

### 4. Evaluate SEACON

```bash
python scripts/evaluate_seacon.py \
  --dataset-dir examples/whole_genome_simulation \
  --seacon-output-dir examples/whole_genome_simulation/seacon_output
```

The evaluator writes `evaluation_metrics.json` in the SEACON output directory and reports HSCN error, LOH precision/recall/F1, breakpoint metrics, and clone assignment metrics.

---

## Parameter Configuration

The config file (`configs/default_benchmark.json`) controls both the simulation and the benchmark run:

```json
{
  "sim_config": {
    "num_clones": 5,
    "num_cells": 100,
    "lambda_events": 1.5,
    "lambda_amplitude": 1.0,
    "prob_wgd": 0.0
  },
  "output_dir": "benchmark_results",
  "methods": ["dummy"],
  "breakpoint_tolerance": 2,
  "event_tolerance": 1
}
```

### Full SimulationConfig Parameters

| Parameter | Default | Description |
|---|---|---|
| `chromosome` | `"chr1"` | Chromosome to simulate (hg38) |
| `bin_width` | `500000` | Bin width in base pairs (500 kb) |
| `num_clones` | `5` | Number of extant clones (K) |
| `num_cells` | `200` | Number of cells (N) |
| `max_copy_number` | `8` | Maximum copy number (C_max) |
| `alpha_tree` | `0.5` | Beta-splitting α (tree balance) |
| `beta_tree` | `0.3` | Beta-splitting β (smaller = more imbalanced) |
| `lambda_events` | `1.5` | Mean CNA events per edge (λ_E) |
| `lambda_amplitude` | `1.0` | Event amplitude parameter (λ_Δ) |
| `prob_wgd` | `0.0` | Probability of WGD event (p_WGD) |
| `gain_prob` | `0.40` | Fraction of events that are gains (60% losses) |
| `prob_focal` | `0.7` | Probability of focal event (vs arm/chromosome) |
| `mean_library_size` | `50000` | Mean reads per cell (ā) |
| `snp_density` | `1e-3` | Heterozygous SNP density (ρ_SNP) |
| `prob_doublet` | `0.0` | Doublet rate (p_doublet) |
| `prob_phase_switch` | `0.01` | Haplotype phase switch rate (p_switch) |
| `random_seed` | `None` | Random seed for reproducibility |

### Benchmark-only Parameters

| Parameter | Default | Description |
|---|---|---|
| `output_dir` | `"benchmark_results"` | Directory to save results |
| `methods` | `["dummy"]` | Methods to run (currently: `"dummy"`) |
| `breakpoint_tolerance` | `2` | Tolerance window in bins for breakpoint matching |
| `event_tolerance` | `1` | Tolerance window for event matching |

---

## Input Format

HaploTreeSim generates all input data internally from the config — **no external input files are required**.

If you are integrating an external inference method and need to feed it simulated data, the simulator produces:

### Read counts matrix
- **Type:** `np.ndarray`, shape `(N, J)`
- **N:** number of cells, **J:** number of bins
- **Values:** integer read counts per cell per bin
- **Usage:** primary input for total CN callers

### Allele counts
- **Type:** tuple of three `np.ndarray`, each shape `(N, S)`
- **S:** number of segments
- **Contents:** `(alt_counts, ref_counts, total_counts)` — alternate allele reads, reference allele reads, and total SNP-overlapping reads per cell per segment
- **Usage:** input for allele-specific CN callers (ASCN/HSCN methods)

### Config file (JSON)
- Passed via `--config` flag
- Controls all simulation and benchmark parameters (see above)

---

## Expected Output

After running, the `benchmark_results/` directory contains:

### `benchmark_results.json`
Full results including config, per-method metrics, and runtime:
```json
{
  "config": { ... },
  "timestamp": "2026-05-19 10:00:00",
  "methods": {
    "dummy_oracle": {
      "metrics": {
        "hscn": {
          "hscn_error": 0.0,
          "loh_f1": 1.0,
          "msr": 1.0
        },
        "breakpoints": { "precision": 1.0, "recall": 1.0, "f1": 1.0 },
        "clone_assignment": { "ari": 1.0, "nmi": 1.0 },
        "tree": { ... }
      },
      "runtime_seconds": 0.01
    }
  }
}
```

### `ground_truth.json`
All ground truth data:

| Field | Type | Description |
|---|---|---|
| `segments` | list of `{start_bin, end_bin, index}` | True segment boundaries |
| `breakpoints` | list of int | True breakpoint bin indices |
| `clone_assignments` | list of int, length N | True clone label per cell |
| `cn_A` | (N, S) int array | True haplotype A copy number per cell per segment |
| `cn_B` | (N, S) int array | True haplotype B copy number per cell per segment |
| `events` | dict keyed by edge | CNA events per tree edge with start_bin, end_bin, haplotype, amplitude |
| `clone_proportions` | list of float | True mixing proportions (π_k) |

### Metrics computed

| Metric | Description |
|---|---|
| `hscn_error` | Swap-invariant haplotype CN MAE (Eq. 57) |
| `loh_f1` | Precision/recall/F1 for LOH detection |
| `msr` | Mirrored-subclone resolution rate (Eq. 60) |
| Breakpoint P/R/F1 | Within ±w bins tolerance window |
| `ari` | Adjusted Rand Index for clone assignment |
| `nmi` | Normalized Mutual Information for clone assignment |
| `nrf` | Normalized Robinson-Foulds tree distance |

---

## For Evaluators

This section describes what an inference method must accept and produce to plug into the benchmark.

### What the simulator provides to your method

```python
predictions = method.run(
    read_counts,    # np.ndarray (N, J) — integer read depths
    allele_counts,  # tuple: (alt_counts, ref_counts, total_counts), each (N, S)
    ground_truth    # dict — ground truth (only used by oracle; real methods should ignore)
)
```

### What your method must return

A dictionary with the following fields:

```python
{
    'cn_A':         np.ndarray,  # shape (N, S) int — inferred haplotype A copy number
    'cn_B':         np.ndarray,  # shape (N, S) int — inferred haplotype B copy number
    'clone_labels': np.ndarray,  # shape (N,) int — predicted clone assignment per cell
    'breakpoints':  np.ndarray,  # shape (B,) int — predicted breakpoint bin indices
    'tree_edges':   list,        # list of (parent, child) int tuples — optional
    'events':       list         # list of event dicts — optional
}
```

Fields `tree_edges` and `events` are optional. If empty, tree and event metrics are skipped.

### Adding a new method

Subclass `BaseMethod` in `benchmark.py`:

```python
class MyMethod(BaseMethod):
    def __init__(self):
        super().__init__("my_method")

    def run(self, read_counts, allele_counts, ground_truth):
        # ... your inference code here ...
        return {
            'cn_A': ...,
            'cn_B': ...,
            'clone_labels': ...,
            'breakpoints': ...,
            'tree_edges': [],
            'events': []
        }
```

Then register it in `BenchmarkRunner._initialize_methods()` and add `"my_method"` to the `methods` list in your config JSON.

---

## Key Features

- **Beta-splitting clone tree** with controlled imbalance
- **Focal/Arm/Chromosome CNA events** with difference-array overlap handling
- **Whole Genome Doubling (WGD)** with early-edge placement
- **Haplotype-specific** copy number per clone
- **Negative Binomial** read depth with library-size variation
- **Beta-Binomial** allelic counts with LOH floor, phase errors, sequencing noise
- **Doublet model** with CN-weighted allelic mixture
- **Full ground truth** export: segments, breakpoints, clone labels, CN profiles, events

---

## Roadmap

- [x] Weeks 5–10: Core simulator
- [x] Week 11: HSCN metrics (16 tests)
- [x] Week 12: Clone assignment metrics (11 tests)
- [x] Week 13: Tree metrics (9 tests)
- [x] Week 14: Benchmark pipeline with CLI (`v1.3.0`)
- [x] Week 15: SEACON integration
- [x] Week 16: CHISEL/Alleloscope integration (partial)
- [ ] Weeks 17–19: CNA-tree baselines + MEDICC2
- [ ] Weeks 20–24: Full experiment grid

---

## Contact

Repository: https://github.com/zaaferahhamid2/haplotreesim

**Version:** 1.3.0 | **Status:** ✅ Benchmark Pipeline Complete | **Updated:** May 2026
