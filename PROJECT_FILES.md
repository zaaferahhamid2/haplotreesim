# HaploTreeSim Project File Guide

Repository root: `/blue/lzhang.uwf/zh34.uwf/haplotreesim`

## Core simulator library (`src/haplotreesim/`)

| File | Description |
|---|---|
| `src/haplotreesim/__init__.py` | Package entry point; exports `SimulationConfig`, `HaploTreeSimulator`, and core data classes. |
| `src/haplotreesim/data_models.py` | Core data structures: `Bin`, `Segment`, `HaplotypeBlock`, `CNAEvent`, `Clone`, `Cell`, and `SimulationConfig` (all simulator parameters). |
| `src/haplotreesim/simulator.py` | Main simulator class (`HaploTreeSimulator`) that orchestrates genome init, tree building, event generation/application, and read/allele observation generation. |
| `src/haplotreesim/tree_builder.py` | Builds basic clone tree topologies. |
| `src/haplotreesim/beta_tree_builder.py` | Beta-splitting tree generator (mass-conserving) producing non-uniform clone proportions and branch lengths tied to event rate. |
| `src/haplotreesim/event_generator.py` | Samples CNA event attributes: scale class (focal/arm/chromosome), length, amplitude, haplotype. |
| `src/haplotreesim/event_applier.py` | Applies CNA events to clone copy-number profiles (difference-array based, handles overlapping events and WGD). |
| `src/haplotreesim/segment_detector.py` | Detects truth segment boundaries from the union of CNA breakpoints across all clones. |
| `src/haplotreesim/chromosome_data.py` | Real hg38 chromosome length reference data. |
| `src/haplotreesim/benchmark.py` | Benchmark pipeline scaffolding (Week 14): standardized runner interface for evaluating tools against simulated datasets. |

## Dataset generation and simulation-adjacent scripts (top level + `scripts/`)

| File | Description |
|---|---|
| `scripts/simulate_dataset.py` | **Primary dataset generator.** Creates one tool-neutral HaploTreeSim dataset (reads, allele counts, tree structure, ground truth) given CLI flags (num-clones, num-cells, lambda-events, beta-tree, mean-library-size, phase-switch-prob, prob-wgd, prob-doublet, etc.). Output feeds all tool-specific converters. |
| `scripts/dataset_io.py` | Shared file I/O helpers for reading/writing the tool-neutral dataset format. |
| `scripts/write_baf.py` | Standalone helper that builds `BAF.tsv` for SEACON by mapping `precomputed_baf.tsv` values onto `filtered_regions.bed` bin indices. |
| `setup.py` | Package installer (`pip install -e .`) for the `haplotreesim` src package. |
| `configs/example_configs.py` | Example `SimulationConfig` presets (e.g. minimal diploid config). |
| `test_tree_reuse.py` | Manual verification script confirming tree structure export/import (`--tree-structure` flag) works correctly — used to validate the "same tree per rep" requirement. |
| `verify_overlapping_events.py` | Manual trace/verification of how overlapping CNA events combine on tree edges. |
| `demonstrate_overlaps.py` | Small standalone demo of overlapping-event handling with a controlled example tree. |

## Plotting scripts (top level + `scripts/`)

| File | Description |
|---|---|
| `create_plots.py` | Generates CN profile plots for sample outputs. |
| `create_chr_plots.py` | Week 8: CN profile plots across real chromosomes. |
| `scripts/week8_cn_profiles_with_wgd.py` | Week 8: compares diploid vs WGD copy-number profiles. |
| `create_week9_plots.py` | Week 9: read-depth observation model sanity plots (mean/variance by CN, depth vs library size, library size distribution). |
| `scripts/week9_depth_with_wgd.py` | Week 9: read-depth vs CN plots with/without WGD (ploidy normalization check). |
| `create_week10_plots.py` | Week 10: BAF/mBAF distribution plots by CN state across coverage regimes. |
| `scripts/week10_baf_with_doublets.py` | Week 10: BAF/mBAF plots including doublet examples. |

## Per-tool pipeline scripts (`scripts/`) — pattern: convert → run → evaluate

Each of the 5 benchmarked tools (SEACON, Alleloscope, SCICoNE, CONET, CNRein) has 3 scripts following the same pattern: convert the shared dataset into that tool's native input format, run the tool, then evaluate its output against ground truth.

| File | Description |
|---|---|
| `scripts/convert_to_seacon.py` | Converts a HaploTreeSim dataset into SEACON's expected input files (read-depth ratio + BAF matrices, normal-cell labels). |
| `scripts/run_seacon.py` | Runs SEACON on a converted input directory (does not simulate or convert). |
| `scripts/evaluate_seacon.py` | Scores SEACON output against ground truth: HSCN error, LOH F1, breakpoint F1, clone ARI/NMI, MSR. |
| `scripts/convert_to_alleloscope.R` | Converts a dataset into Alleloscope's binned alt/ref count matrix format. |
| `scripts/run_alleloscope.R` | Runs Alleloscope; contains a runtime monkey-patch for an internal package bug in `Matrix_filter` (chromosomes losing all SNPs during filtering were being silently dropped, causing a crash). |
| `scripts/evaluate_alleloscope.R` | Scores Alleloscope output (all-cells and retained-cells metrics, since Alleloscope filters out some cells); includes a safe lookup fix for chromosomes missing from the output. |
| `scripts/convert_to_scicone.py` | Converts a dataset into SCICoNE's cells-by-bins read-count CSV format. |
| `scripts/run_scicone.py` | Runs SCICoNE via the `pyscicone` Python API; each job uses a unique `postfix` (dataset name) to avoid concurrent-job temp-file collisions. |
| `scripts/evaluate_scicone.py` | Scores SCICoNE output: TCN MSE, clone ARI/NMI, breakpoint F1. |
| `scripts/convert_to_conet.py` | Converts a dataset into CONET's corrected-counts matrix + breakpoint-candidate format. |
| `scripts/run_conet.py` | Runs CONET. |
| `scripts/evaluate_conet.py` | Scores CONET output: TCN MSE, breakpoint F1, clone ARI/NMI, nRF tree distance, tree coverage, cell-node match accuracy. |
| `scripts/convert_to_cnrein.py` | Converts a dataset into CNRein's RDR/BAF numpy array format. |
| `scripts/run_cnrein.py` | Runs CNRein's CNNaive step directly on numpy arrays (bypasses BAM processing). |
| `scripts/evaluate_cnrein.py` | Scores CNRein output: HSCN error, TCN MSE, LOH F1, breakpoint F1, clone ARI/NMI. |

## Experiment orchestration (`scripts/experiments/`)

| File | Description |
|---|---|
| `scripts/experiments/generate_experiments.sh` | Generates all experiment datasets (one-factor-at-a-time sweeps: clone number, tree balance, event rate, cell count, normal fraction, event amplitude, focal/broad ratio, plus coverage/phase-switch/WGD sweeps added later). Reuses the clone-4 tree per seed for every non-tree-structure experiment. |
| `scripts/experiments/convert_all.sh` | Runs all 5 tools' converters across every generated experiment dataset; skips datasets already converted. |
| `scripts/experiments/generate_slurm_jobs.py` | Generates and submits Slurm job scripts per tool/dataset. Usage: `--tool X --test --dataset Y` for a single test job, `--tool X --submit` for the full grid. Contains per-tool CPU/memory/time configs (`TOOL_CONFIG`). |
| `scripts/show_structure.sh` | Prints the repository directory structure. |

## Tests (`tests/`)

| File | Description |
|---|---|
| `tests/test_week5.py` | Verifies the minimal diploid simulator (bins, segments, diploid tree, basic sampling). |
| `tests/test_week6.py` | Verifies CNA event generation/application and tree structure (Week 6 deliverable). |
| `tests/test_week7.py` | Verifies segment boundary detection from CNA breakpoints. |
| `tests/test_beta_splitting.py` | Verifies the Beta-splitting tree model (topology, clone proportions, branch lengths). |
| `tests/test_event_priors.py` | Verifies event attribute priors (scale-class probabilities, focal length distribution, haplotype selection, amplitude distribution). |
| `tests/test_hscn_metrics.py` | Verifies haplotype-specific copy-number error metric calculations (Week 11). |
| `tests/test_clone_metrics.py` | Verifies clone-assignment metrics: ARI, NMI, Hungarian matching (Week 12). |
| `tests/test_tree_metrics.py` | Verifies tree-comparison metrics: Robinson-Foulds distance, ancestor-descendant accuracy, event matching (Week 13). |
| `tests/test_breakpoint_metrics.py` | Verifies breakpoint precision/recall/F1 calculations. |
| `tests/test_metrics_validation.py` | Sanity-check suite requested by professor: confirms metrics score perfect predictions as perfect, and correctly degrade on corrupted data. |
| `scripts/comprehensive_test.py` | End-to-end test verifying all Phase A + B simulator features together (full paper-compliance check). |

## Data and output directories (not tracked in git — see `.gitignore`)

| Directory | Contents |
|---|---|
| `experiments/` | Raw simulated datasets, one folder per condition (e.g. `clone_4_rep42`, `coverage_25_rep43`). Fully reproducible from code + seed + tree structure; never pushed to GitHub. |
| `SEACON/`, `Alleloscope/`, `SCICoNE/`, `CONET/`, `CNRein/` | Per-tool converted input files, one folder per dataset. |
| `results/` | Per-tool, per-dataset output: `metrics.json` (accuracy metrics), `metrics.json.runtime` (elapsed time log), plus raw tool output files. Also contains the aggregated summary CSVs (see below). |
| `results/runtime_summary.csv` | Aggregated elapsed-time-per-job across all tools/datasets (tracked in git — small file). |
| `results/metrics_summary_full.csv` | Aggregated accuracy metrics across all tools/datasets (tracked in git — small file). |
| `logs/` | Per-job Slurm stdout/stderr logs, one file per tool/dataset. |
| `slurm_scripts/` | Generated Slurm submission scripts, one per tool/dataset (written by `generate_slurm_jobs.py`). |

## Notes

- All tool environments are loaded via `module load conda && conda activate <env>` (e.g. `haplotreesim`, `SEACON`, `CONET`) — see `TOOL_CONFIG` in `generate_slurm_jobs.py` for which environment each tool needs.
- R scripts need `module load R/4.1` before running interactively.
- The core benchmark grid (80 datasets × 5 tools = 400 runs) follows the parameter design approved by the professor via email (num_clones=4/2/6/8, tree balance beta 0.1/0.3/0.5, lambda_events 10/40/80, cells 50/200/500, prob_normal 0/0.1/0.3, amplitude 1/2/4, focal 0.3/0.7/1.0), which differs from the illustrative grid in the manuscript's Supplementary Table S1.
- Coverage, phase-switch, WGD, and doublet-rate experiments were added later per direct professor request; parameter values used may need reconciliation with Table S1 (see email thread).
