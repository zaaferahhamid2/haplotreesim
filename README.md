# HaploTreeSim
**A Controlled, End-to-End Benchmark with Haplotype-Resolved CNA and CNA-Phylogeny Ground Truth for Low-Pass scDNA-seq**

[![Version](https://img.shields.io/badge/version-3.0.0-blue.svg)](https://github.com/zaaferahhamid2/haplotreesim/releases)
[![Status](https://img.shields.io/badge/status-benchmark%20complete-brightgreen.svg)](#)

HaploTreeSim is a simulator for generating synthetic single-cell DNA sequencing (scDNA-seq) datasets with ground-truth haplotype-specific copy number alterations (CNAs) and clone phylogenies. It benchmarks five external methods — SEACON, Alleloscope, CNRein (CNNaive), SCICoNE, and CONET — under controlled, matched-ground-truth conditions.

**Full benchmark status:** 645 completed method-runs across the core 16-condition parameter grid plus four robustness experiments (sequencing coverage, phase-switch error, whole-genome doubling, doublet rate). See `results/master_summary.csv` for all parameters and metrics, `PROJECT_FILES.md` for a full file-by-file guide, and `plots/` for the manuscript figures and the script that generates them.

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Full Benchmark Results](#full-benchmark-results)
- [Project Structure](#project-structure)
- [Running the Full Pipeline](#running-the-full-pipeline)
- [SEACON Pipeline](#seacon-pipeline)
- [Alleloscope Pipeline](#alleloscope-pipeline)
- [SCICoNE Pipeline](#scicone-pipeline)
- [CONET Pipeline](#conet-pipeline)
- [CNRein Pipeline](#cnrein-pipeline)
- [Parameter Configuration](#parameter-configuration)
- [Input Format](#input-format)
- [For Evaluators](#for-evaluators)
- [Key Features](#key-features)
- [Known Limitations](#known-limitations)
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

## Full Benchmark Results

The complete benchmark grid is finished: 400 runs across the core 16-condition parameter sweep, plus four dedicated robustness experiments (coverage, phase-switch, whole-genome doubling, doublet rate), for 645 total method-runs across all five tools.

**Key findings:**
- Whole-genome doubling produced the largest, most consistent accuracy degradation of any factor tested, across every tool.
- Coverage shows a clear dose-response for SEACON and CNRein, plateauing around 0.05x, while Alleloscope fails outright below roughly 0.01x coverage (zero surviving SNPs after internal filtering) — a genuine, reproducible limitation rather than a bug.
- Phase-switch error selectively affects SEACON, with Alleloscope and CNRein comparatively robust.
- Doublet rate in the 0.02–0.05 range has minimal impact across all five methods.

Full aggregated results, including every parameter and metric for every run, are in `results/master_summary.csv`. Manuscript figures (coverage response, phase-switch sensitivity, and CNA phylogeny accuracy) are in `plots/`, along with the script used to generate them.

---

## Project Structure

See `PROJECT_FILES.md` for a description of every script and its path.

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

This writes: bins.tsv, segments.tsv, cells.tsv, readcounts.tsv, rdr.tsv, allele_alt.tsv, allele_ref.tsv, snp_allele_alt.mtx, snp_allele_ref.mtx, snps.tsv, clone_cn_A.tsv, clone_cn_B.tsv, truth_cell_hscn_segments.tsv, events.tsv, tree_structure.json, config.json.

### Running the full experiment grid

```bash
bash scripts/experiments/generate_experiments.sh   # generate all datasets
bash scripts/experiments/convert_all.sh             # convert to every tool's input format
python3 scripts/experiments/generate_slurm_jobs.py --tool SEACON --test   # test one job
python3 scripts/experiments/generate_slurm_jobs.py --tool SEACON --submit # submit full grid
```

Repeat the last two commands for each tool: `SEACON`, `Alleloscope`, `SCICoNE`, `CONET`, `CNRein`.

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

**Example single-run result** (whole genome, 100 cells, lambda_events=40; see `results/master_summary.csv` for the full 80+ run aggregate):

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

**Example single-run result:**

| Metric | Value |
|---|---|
| HSCN Error (retained) | 0.16 |
| LOH F1 (retained) | 0.68 |
| Retention | 61/100 cells |

**Note:** Alleloscope requires a minimum allelic coverage of roughly 0.01–0.02x to function; below that, zero SNPs survive its internal filtering and the run fails. This is a genuine finding, not a bug — see `plots/figure2_coverage.png`.

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

**Example single-run result:**

| Metric | Value |
|---|---|
| TCN MSE | 1.17 |
| Clone ARI | 0.74 |
| Clone NMI | 0.80 |
| Breakpoint F1 | 0.19 |

**Note:** SCICoNE's evaluation script does not currently compute tree-topology metrics (normalized RF distance, tree coverage) — see [Known Limitations](#known-limitations).

**Note on running multiple SCICoNE jobs concurrently:** each job must use a unique `postfix` value to avoid a temp-file collision in the underlying `pyscicone` library; `run_scicone.py` sets this automatically using the dataset name.

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

**Example single-run result:**

| Metric | Value |
|---|---|
| TCN MSE | 2.42 |
| Breakpoint F1 | 0.45 (P=0.955, R=0.292) |
| nRF Distance | 0.40 |
| TreeCoverage | 0.83 |

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

**Example single-run result:**

| Metric | Value |
|---|---|
| HSCN Error | 0.583 |
| LOH F1 | 0.000 |
| TCN MSE | 2.043 |
| Clone ARI | 0.465 |
| Clone NMI | 0.719 |
| Breakpoint F1 | 0.047 (P=0.024, R=0.889) |

**Note:** across the full benchmark grid, CNRein showed a consistent, near-total weakness in breakpoint detection (mean F1 = 0.01 across all 80 core-grid datasets) despite reasonable HSCN error overall.

---

## Parameter Configuration

| Parameter | Default | Description |
|---|---|---|
| `chromosome` | `"chr1"` | Chromosome to simulate (hg38) |
| `bin_width` | `500000` | Bin width in base pairs (500 kb) |
| `num_clones` | `4` | Number of extant clones (K) |
| `num_cells` | `200` | Number of cells (N) |
| `lambda_events` | `40.0` | Mean CNA events per edge |
| `lambda_amplitude` | `2.0` | Event amplitude parameter |
| `prob_wgd` | `0.0` | Probability of WGD event |
| `prob_normal` | `0.1` | Fraction of normal (diploid) cells |
| `prob_doublet` | `0.0` | Fraction of doublet cells |
| `prob_focal` | `0.7` | Probability of focal CNA events |
| `mean_library_size` | `100.0` | Mean read-depth coverage factor |
| `mean_allelic_coverage` | auto-calibrated | Mean per-SNP allelic coverage |
| `library_size_cv` | `0.3` | Coefficient of variation, library size |
| `allelic_coverage_cv` | `0.3` | Coefficient of variation, allelic coverage |
| `prob_phase_switch` | `0.01` | Phase-switch probability per chromosome-arm block |
| `snp_density` | `1e-3` | Heterozygous SNP density per bp |
| `random_seed` | `None` | Random seed for reproducibility |

Use `--whole-genome` flag in `simulate_dataset.py` to simulate all autosomes (chr1–chr22).

---

## Input Format

- `readcounts.tsv` — (N, J) integer read counts per cell per bin
- `allele_alt.tsv`, `allele_ref.tsv` — (N, J) allele counts per cell per bin
- `snp_allele_alt.mtx`, `snp_allele_ref.mtx` — sparse SNP-level allele counts
- `snps.tsv` — SNP positions (chrom, start, end)
- `truth_cell_hscn_segments.tsv` — per-cell per-segment ground truth CN

---

## For Evaluators

To add a new tool create three scripts:
- `scripts/convert_to_TOOLNAME.py` — converts shared dataset to tool input
- `scripts/run_TOOLNAME.py` — runs the tool
- `scripts/evaluate_TOOLNAME.py` — evaluates against truth_cell_hscn_segments.tsv

Then add the tool's resource configuration (conda env, CPUs, memory, time limit) to `TOOL_CONFIG` in `scripts/experiments/generate_slurm_jobs.py`.

---

## Key Features

- Beta-splitting clone tree with controlled imbalance
- Focal/Arm/Chromosome CNA events with difference-array overlap handling
- Whole Genome Doubling (WGD) with early-edge placement
- Haplotype-specific copy number per clone
- SNP-level allele counts for allele-specific tools
- Normal cell and doublet contamination
- Negative Binomial read depth with library-size variation
- Beta-Binomial allelic counts with LOH floor, phase errors, sequencing noise
- Full ground truth export: segments, breakpoints, clone labels, CN profiles, events, tree
- Automated SLURM benchmark pipeline with per-job runtime tracking

---

## Known Limitations

- **SCICoNE tree metrics:** `scripts/evaluate_scicone.py` computes clone assignment, total-CN error, and breakpoint metrics, but does not compute tree-topology metrics (normalized RF distance, tree coverage, cell-node match accuracy) that `evaluate_conet.py` does, despite both being CNA-phylogeny methods. See `PROJECT_FILES.md` for detail.
- **Alleloscope low-coverage failure:** Alleloscope cannot process datasets below roughly 0.01–0.02x allelic coverage, since zero SNPs survive its internal filtering step at that level. This is documented as a genuine finding, not a defect.

---

## Roadmap

- [x] Weeks 1–4: Literature review and background
- [x] Weeks 5–10: Core simulator
- [x] Week 11: HSCN metrics
- [x] Week 12: Clone assignment metrics
- [x] Week 13: Tree metrics
- [x] Week 14: Benchmark pipeline with CLI
- [x] SEACON, Alleloscope, SCICoNE, CONET, CNRein integration (all 5 tools)
- [x] Core 16-condition, 5-replicate, 5-tool benchmark grid (400 runs)
- [x] Coverage robustness experiment (100 runs, 4 non-default coverage levels)
- [x] Phase-switch robustness experiment (50 runs)
- [x] Whole-genome-doubling robustness experiment (50 runs)
- [x] Doublet-rate robustness experiment (50 runs)
- [x] Full results aggregation (`results/master_summary.csv`, 645 rows)
- [x] Manuscript figures (coverage response, phase-switch sensitivity, CNA phylogeny)
- [x] Full project file documentation (`PROJECT_FILES.md`)

---

## Contact

Repository: https://github.com/zaaferahhamid2/haplotreesim

**Version:** 3.0.0 | **Status:** Full benchmark complete (645 runs) | **Updated:** August 2026