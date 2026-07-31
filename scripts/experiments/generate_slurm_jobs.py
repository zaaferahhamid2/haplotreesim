#!/usr/bin/env python3
"""
Generate and optionally submit Slurm job scripts for the benchmark grid.
Usage:
  python3 scripts/experiments/generate_slurm_jobs.py --tool SEACON --test   # test job only
  python3 scripts/experiments/generate_slurm_jobs.py --tool SEACON --submit # submit all jobs
"""

import argparse
import os
import subprocess
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
ACCOUNT = "lzhang.uwf"
QOS = "lzhang.uwf"
BASE_DIR = Path("/blue/lzhang.uwf/zh34.uwf/haplotreesim")
EXP_DIR = BASE_DIR / "experiments"
CONV_DIR = BASE_DIR
OUT_DIR = BASE_DIR / "results"
LOG_DIR = BASE_DIR / "logs"
TOOLS_DIR = Path("/blue/lzhang.uwf/zh34.uwf/tools")
R_LIBS = "/blue/lzhang.uwf/zh34.uwf/R_libs"

TOOL_CONFIG = {
    "SEACON": {
        "conda_env": "SEACON",
        "eval_env": "haplotreesim",
        "cpus": 4,
        "mem": "32gb",
        "time": "06:00:00",
        "run_script": "scripts/run_seacon.py",
        "eval_script": "scripts/evaluate_seacon.py",
        "run_args": "--upper-filter 20 --tolerance 0.05 --max-wgd 1",
        "extra_modules": "",
    },
    "Alleloscope": {
        "conda_env": "haplotreesim",
        "eval_env": "haplotreesim",
        "cpus": 4,
        "mem": "96gb",
        "time": "03:00:00",
        "run_script": "scripts/run_alleloscope.R",
        "eval_script": "scripts/evaluate_alleloscope.R",
        "run_args": "",
        "extra_modules": "R/4.1",
    },
    "SCICoNE": {
        "conda_env": "haplotreesim",
        "eval_env": "haplotreesim",
        "cpus": 4,
        "mem": "16gb",
        "time": "02:00:00",
        "run_script": "scripts/run_scicone.py",
        "eval_script": "scripts/evaluate_scicone.py",
        "run_args": f"--scicone-build {TOOLS_DIR}/SCICoNE/build",
        "extra_modules": "",
    },
    "CONET": {
        "conda_env": "CONET",
        "eval_env": "haplotreesim",
        "cpus": 4,
        "mem": "16gb",
        "time": "00:30:00",
        "run_script": "scripts/run_conet.py",
        "eval_script": "scripts/evaluate_conet.py",
        "run_args": f"--conet-bin {TOOLS_DIR}/CONET/src/CONET --counts-penalty 100000",
        "extra_modules": "",
    },
    "CNRein": {
        "conda_env": "CONET",
        "eval_env": "haplotreesim",
        "cpus": 4,
        "mem": "16gb",
        "time": "03:00:00",
        "run_script": "scripts/run_cnrein.py",
        "eval_script": "scripts/evaluate_cnrein.py",
        "run_args": "",
        "extra_modules": "",
    },
}

def get_datasets():
    datasets = []
    for d in sorted(EXP_DIR.iterdir()):
        if d.is_dir() and any(d.name.startswith(p) for p in
                              ["clone_", "beta_", "events_", "cells_",
                               "normal_", "amplitude_", "focal_", "coverage_", "phaseswitch_", "wgd_", "doublet_"]):
            datasets.append(d.name)
    return datasets


def make_slurm_script(tool, dataset, cfg):
    conv_input = CONV_DIR / tool / dataset
    run_output = OUT_DIR / tool / dataset
    log_file = LOG_DIR / tool / f"{dataset}.log"
    metrics_out = run_output / "metrics.json"

    if tool == "Alleloscope":
        run_cmd = f"""
module load {cfg['extra_modules']}
export R_LIBS_USER={R_LIBS}

echo "=== CONVERT ==="
Rscript {BASE_DIR}/scripts/run_alleloscope.R \\
    --input-dir {conv_input} \\
    --output-dir {run_output} \\
    2>&1

echo "=== EVALUATE ==="
Rscript {BASE_DIR}/scripts/evaluate_alleloscope.R \\
    --dataset-dir {EXP_DIR}/{dataset} \\
    --output-dir {run_output} \\
    --metrics-out {metrics_out} \\
    2>&1
"""
    elif tool == "SEACON":
        run_cmd = f"""
conda activate {cfg['conda_env']}
python3 {BASE_DIR}/{cfg['run_script']} \\
    --input-dir {conv_input} \\
    --output-dir {conv_input} \\
    {cfg['run_args']} 2>&1

conda activate {cfg['eval_env']}
mkdir -p {run_output}
python3 {BASE_DIR}/{cfg['eval_script']} \\
    --dataset-dir {EXP_DIR}/{dataset} \\
    --seacon-output-dir {conv_input} \\
    --metrics-out {metrics_out} \\
    2>&1
"""
    elif tool == "SCICoNE":
        run_cmd = f"""
conda activate {cfg['conda_env']}
python3 {BASE_DIR}/{cfg['run_script']} \\
    --input-dir {conv_input} \\
    --output-dir {run_output} \\
    {cfg['run_args']} \\
    --overwrite 2>&1

python3 {BASE_DIR}/{cfg['eval_script']} \\
    --dataset-dir {EXP_DIR}/{dataset} \\
    --scicone-output-dir {run_output} \\
    --metrics-out {metrics_out} \\
    2>&1
"""
    elif tool == "CONET":
        run_cmd = f"""
export LD_LIBRARY_PATH=/blue/lzhang.uwf/zh34.uwf/.conda/envs/CONET/lib:$LD_LIBRARY_PATH
conda activate {cfg['conda_env']}
python3 {BASE_DIR}/{cfg['run_script']} \\
    --input-dir {conv_input} \\
    --output-dir {run_output} \\
    {cfg['run_args']} \\
    --overwrite 2>&1

conda activate {cfg['eval_env']}
python3 {BASE_DIR}/{cfg['eval_script']} \\
    --dataset-dir {EXP_DIR}/{dataset} \\
    --conet-output-dir {run_output} \\
    --metrics-out {metrics_out} \\
    2>&1
"""
    elif tool == "CNRein":
        run_cmd = f"""
conda activate {cfg['conda_env']}
python3 {BASE_DIR}/{cfg['run_script']} \\
    --input-dir {conv_input} \\
    --output-dir {run_output} \\
    --overwrite 2>&1

conda activate {cfg['eval_env']}
python3 {BASE_DIR}/{cfg['eval_script']} \\
    --dataset-dir {EXP_DIR}/{dataset} \\
    --cnrein-output-dir {run_output} \\
    --metrics-out {metrics_out} \\
    2>&1
"""

    script = f"""#!/bin/bash
#SBATCH --job-name={tool[:3]}_{dataset[:20]}
#SBATCH --account={ACCOUNT}
#SBATCH --qos={QOS}
#SBATCH --cpus-per-task={cfg['cpus']}
#SBATCH --mem={cfg['mem']}
#SBATCH --time={cfg['time']}
#SBATCH --output={log_file}
#SBATCH --error={log_file}

echo "Job started: $(date)"
echo "Dataset: {dataset}"
echo "Tool: {tool}"
START_TIME=$(date +%s)

cd {BASE_DIR}
module load conda

mkdir -p {run_output}
mkdir -p {log_file.parent}

{run_cmd}

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo "Job finished: $(date)"
echo "Elapsed time: ${{ELAPSED}} seconds"
echo "RUNTIME_SECONDS=${{ELAPSED}}" >> {metrics_out}.runtime
"""
    return script


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", required=True, choices=list(TOOL_CONFIG.keys()))
    parser.add_argument("--test", action="store_true", help="Generate and submit test job only")
    parser.add_argument("--dataset", default=None, help="Specific dataset to use for test job")
    parser.add_argument("--submit", action="store_true", help="Submit all jobs")
    parser.add_argument("--dry-run", action="store_true", help="Print scripts without submitting")
    args = parser.parse_args()

    tool = args.tool
    cfg = TOOL_CONFIG[tool]
    datasets = get_datasets()

    if not datasets:
        print("No datasets found. Run generate_experiments.sh first.")
        return

    slurm_dir = BASE_DIR / "slurm_scripts" / tool
    slurm_dir.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    (LOG_DIR / tool).mkdir(parents=True, exist_ok=True)

    if args.test:
        # Use specified dataset or first available
        test_ds = args.dataset if args.dataset else datasets[0]
        script = make_slurm_script(tool, test_ds, cfg)
        script_path = slurm_dir / f"test_{test_ds}.sh"
        script_path.write_text(script)
        print(f"Test script written: {script_path}")
        if not args.dry_run:
            result = subprocess.run(["sbatch", str(script_path)], capture_output=True, text=True)
            print(result.stdout)
            print(result.stderr)
        return

    if args.submit:
        submitted = 0
        for ds in datasets:
            script = make_slurm_script(tool, ds, cfg)
            script_path = slurm_dir / f"{ds}.sh"
            script_path.write_text(script)
            if not args.dry_run:
                result = subprocess.run(["sbatch", str(script_path)], capture_output=True, text=True)
                print(f"{ds}: {result.stdout.strip()}")
                submitted += 1
            else:
                print(f"[dry-run] Would submit: {script_path}")
        print(f"\nSubmitted {submitted} jobs for {tool}")


if __name__ == "__main__":
    main()
