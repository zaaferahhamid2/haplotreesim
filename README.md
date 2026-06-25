# HaploTreeSim
**A Controlled, End-to-End Benchmark with Haplotype-Resolved CNA and CNA-Phylogeny Ground Truth for Low-Pass scDNA-seq**

[![Version](https://img.shields.io/badge/version-2.1.0-blue.svg)](https://github.com/zaaferahhamid2/haplotreesim/releases/tag/v2.1.0)
[![Paper Compliant](https://img.shields.io/badge/paper-compliant-brightgreen.svg)](#)

HaploTreeSim is a simulator for generating synthetic single-cell DNA sequencing (scDNA-seq) datasets with ground-truth haplotype-specific copy number alterations (CNAs) and clone phylogenies. It is designed to benchmark ASCN callers, CNA-tree methods, and joint phylogeny pipelines under controlled conditions.

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Running the Full Pipeline](#running-the-full-pipeline)
- [SEACON Pipeline](#seacon-pipeline)
- [Alleloscope Pipeline](#alleloscope-pipeline)
- [SCICoNE Pipeline](#scicone-pipeline)
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

## Quick Start

```python
from haplotreesim import SimulationConfig, HaploTreeSimulator

config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=200,
    lambda_events=3.0,
    prob_wgd=0.0,
    random_seed=42
)

sim = HaploTreeSimulator(config)
read_counts, allele_counts = sim.run()
alt_counts, ref_counts, total_counts = allele_counts
ground_truth = sim.get_ground_truth()
```

---

## Running the Full Pipeline

All tools share the same simulated dataset. The workflow is:

### Step 1 — Simulate a dataset

```bash
python3 scripts/simulate_dataset.py \
  --output-dir examples/simulation_wgs \
  --whole-genome \
  --num-clones 5 \
  --num-cells 100 \
  --lambda-events 40 \
  --lambda-amplitude 2.0 \
  --prob-normal 0.1 \
  --random-seed 42 \
  --overwrite
```

This writes: bins.tsv, segments.tsv, cells.tsv, readcounts.tsv, rdr.tsv, allele_alt.tsv, allele_ref.tsv, snp_allele_alt.mtx, snp_allele_ref.mtx, snps.tsv, clone_cn_A.tsv, clone_cn_B.tsv, truth_cell_hscn_segments.tsv, events.tsv, tree_structure.json, metadata.json.

---

## SEACON Pipeline

```bash
python3 scripts/convert_to_seacon.py \
  --input-dir examples/simulation_wgs \
  --output-dir examples/simulation_wgs/seacon_input \
  --overwrite

conda activate SEACON
python3 scripts/run_seacon.py \
  --input-dir examples/simulation_wgs/seacon_input \
  --output-dir examples/simulation_wgs/seacon_output \
  --upper-filter 20 --tolerance 0.05 --max-wgd 1 --overwrite

conda deactivate
python3 scripts/evaluate_seacon.py \
  --dataset-dir examples/simulation_wgs \
  --seacon-output-dir examples/simulation_wgs/seacon_output \
  --metrics-out examples/simulation_wgs/seacon_metrics.json
```

**Results (whole genome, 100 cells, lambda_events=40):**

| Metric | Value |
|---|---|
| HSCN Error | 0.47 |
| LOH F1 | 0.88 |
| Clone ARI | 0.52 |
| Clone NMI | 0.63 |

---

## Alleloscope Pipeline

Requires R with Alleloscope installed: `remotes::install_github("seasoncloud/Alleloscope")`

```bash
Rscript scripts/convert_to_alleloscope.R \
  --input-dir examples/simulation_wgs \
  --output-dir examples/simulation_wgs/alleloscope_input

Rscript scripts/run_alleloscope.R \
  --input-dir examples/simulation_wgs/alleloscope_input \
  --output-dir examples/simulation_wgs/alleloscope_output

Rscript scripts/evaluate_alleloscope.R \
  --dataset-dir examples/simulation_wgs \
  --output-dir examples/simulation_wgs/alleloscope_output \
  --metrics-out examples/simulation_wgs/alleloscope_metrics.json
```

**Results (whole genome, 100 cells, lambda_events=40):**

| Metric | Value |
|---|---|
| HSCN Error (retained) | 0.16 |
| LOH F1 (retained) | 0.68 |
| Retention | 61/100 cells |

---

## SCICoNE Pipeline

Requires C++ compilation:
```bash
brew install cmake boost nlopt
cd ~/Documents/SCICoNE && mkdir -p build && cd build
cmake ../scicone/ && make -j4
```

```bash
python3 scripts/convert_to_scicone.py \
  --input-dir examples/simulation_wgs \
  --output-dir examples/simulation_wgs/scicone_input \
  --overwrite

python3 scripts/run_scicone.py \
  --input-dir examples/simulation_wgs/scicone_input \
  --output-dir examples/simulation_wgs/scicone_output \
  --scicone-build ~/Documents/SCICoNE/build \
  --overwrite

python3 scripts/evaluate_scicone.py \
  --dataset-dir examples/simulation_wgs \
  --scicone-output-dir examples/simulation_wgs/scicone_output \
  --metrics-out examples/simulation_wgs/scicone_metrics.json
```

**Results (whole genome, 100 cells, lambda_events=40):**

| Metric | Value |
|---|---|
| TCN MSE | 1.17 |
| Clone ARI | 0.74 |
| Clone NMI | 0.80 |
| Breakpoint F1 | 0.19 |

---

## CONET Pipeline

Requires C++ compilation and conet-py:

```bash
cd ~/Documents/CONET/src
sed -i "" "s|/usr/local|/opt/homebrew|g" Makefile
make
pip install -e ~/Documents/CONET/python/conet-py/
```

```bash
python3 scripts/convert_to_conet.py \
  --input-dir examples/simulation_wgs \
  --output-dir examples/simulation_wgs/conet_input \
  --overwrite

conda activate CONET
python3 scripts/run_conet.py \
  --input-dir examples/simulation_wgs/conet_input \
  --output-dir examples/simulation_wgs/conet_output \
  --counts-penalty 100 \
  --overwrite

python3 scripts/evaluate_conet.py \
  --dataset-dir examples/simulation_wgs \
  --conet-output-dir examples/simulation_wgs/conet_output \
  --metrics-out examples/simulation_wgs/conet_metrics.json
```

**Results (whole genome, 100 cells, lambda_events=40):**

| Metric | Value |
|---|---|
| TCN MSE | 2.42 |
| Breakpoint F1 | 0.45 (P=0.955, R=0.292) |
| nRF Distance | 0.40 |
| TreeCoverage | 0.83 |

**Note:** CONET infers CN values of 1-2 only on our dataset. This is a known limitation when true CN amplitudes exceed CONET's typical operating range.

---

## CNRein Pipeline

Requires CNRein installed via pip (in a conda env with Python 3.9):

```bash
conda create -n CONET python=3.9 -y
conda activate CONET
pip install CNRein scikit-learn
```

```bash
conda activate CONET
cd ~/Documents/haplotreesim

python3 scripts/convert_to_cnrein.py \
  --input-dir examples/simulation_wgs \
  --output-dir examples/simulation_wgs/cnrein_output \
  --overwrite

python3 scripts/run_cnrein.py \
  --input-dir examples/simulation_wgs/cnrein_output \
  --output-dir examples/simulation_wgs/cnrein_run \
  --overwrite

source venv/bin/activate
python3 scripts/evaluate_cnrein.py \
  --dataset-dir examples/simulation_wgs \
  --cnrein-output-dir examples/simulation_wgs/cnrein_run \
  --metrics-out examples/simulation_wgs/cnrein_metrics.json
```

**Results (whole genome, 100 cells, lambda_events=40):**

| Metric | Value |
|---|---|
| HSCN Error | 0.583 |
| LOH F1 | 0.000 |
| TCN MSE | 2.043 |
| Clone ARI | 0.465 |
| Clone NMI | 0.719 |
| Breakpoint F1 | 0.047 (P=0.024, R=0.889) |

**Note:** CNRein correctly detects CN values 1-11. High breakpoint recall (0.889) but low precision due to fine-grained segmentation. BAM processing step is bypassed by generating numpy inputs directly from simulator outputs.

---

## Parameter Configuration

| Parameter | Default | Description |
|---|---|---|
| `chromosome` | `"chr1"` | Chromosome to simulate (hg38) |
| `bin_width` | `500000` | Bin width in base pairs (500 kb) |
| `num_clones` | `5` | Number of extant clones (K) |
| `num_cells` | `200` | Number of cells (N) |
| `lambda_events` | `1.5` | Mean CNA events per edge |
| `lambda_amplitude` | `1.0` | Event amplitude parameter |
| `prob_wgd` | `0.0` | Probability of WGD event |
| `prob_normal` | `0.0` | Fraction of normal (diploid) cells |
| `prob_focal` | `0.7` | Probability of focal CNA events |
| `snp_density` | `1e-3` | Heterozygous SNP density per bp |
| `random_seed` | `None` | Random seed for reproducibility |

Use `--whole-genome` flag in `simulate_dataset.py` to simulate all autosomes (chr1-chr22).

---

## Input Format

- `readcounts.tsv` — (N, J) integer read counts per cell per bin
- `allele_alt.tsv`, `allele_ref.tsv` — (N, J) allele counts per cell per bin
- `snp_allele_alt.mtx`, `snp_allele_ref.mtx` — sparse SNP-level allele counts
- `snps.tsv` — SNP positions (chrom, start, end)
- `truth_cell_hscn_segments.tsv` — per-cell per-segment ground truth CN

---

## Expected Output

Each tool writes a metrics.json with:

**ASCN tools (SEACON, Alleloscope):** hscn_error, loh_f1, msr, clone_ari, clone_nmi

**CNA-tree tools (SCICoNE):** tcn_mse, clone_ari, clone_nmi, breakpoint_f1

---

## For Evaluators

To add a new tool create three scripts:
- `scripts/convert_to_TOOLNAME.py` — converts shared dataset to tool input
- `scripts/run_TOOLNAME.py` — runs the tool
- `scripts/evaluate_TOOLNAME.py` — evaluates against truth_cell_hscn_segments.tsv

---

## Key Features

- Beta-splitting clone tree with controlled imbalance
- Focal/Arm/Chromosome CNA events with difference-array overlap handling
- Whole Genome Doubling (WGD) with early-edge placement
- Haplotype-specific copy number per clone
- SNP-level allele counts for allele-specific tools
- Normal cell contamination via --prob-normal
- Negative Binomial read depth with library-size variation
- Beta-Binomial allelic counts with LOH floor, phase errors, sequencing noise
- Full ground truth export: segments, breakpoints, clone labels, CN profiles, events, tree

---

## Roadmap

- [x] Weeks 5-10: Core simulator
- [x] Week 11: HSCN metrics
- [x] Week 12: Clone assignment metrics
- [x] Week 13: Tree metrics
- [x] Week 14: Benchmark pipeline with CLI
- [x] Week 15: SEACON integration
- [x] Week 16: Alleloscope integration (whole-genome)
- [x] Week 17: SCICoNE integration (whole-genome)
- [x] Week 18: CONET integration (complete, CN range limitation noted)
- [x] Week 19: CNRein integration (complete, CNNaive step)
- [ ] Weeks 20-24: Full experiment grid

---

## Contact

Repository: https://github.com/zaaferahhamid2/haplotreesim

**Version:** 2.1.0 | **Status:** SCICoNE Complete | **Updated:** June 2026
