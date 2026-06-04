"""
Run SEACON on a prepared SEACON input directory.

This script does not simulate data. Build a general dataset first with
simulate_dataset.py, convert it with convert_to_seacon.py, then run this script.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from dataset_io import ensure_clean_dir, require_files, write_json


SEACON_INPUT_FILES = [
    "cells.txt",
    "filtered_regions.bed",
    "readcounts.tsv",
    "RDR.tsv",
    "precomputed_baf.tsv",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing SEACON input files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory where SEACON should run and write outputs.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace output-dir if it exists.")
    parser.add_argument("--seacon-bin", default="seacon", help="SEACON executable name or path.")
    parser.add_argument("--upper-filter", type=int, default=20)
    parser.add_argument("--tolerance", type=float, default=0.05)
    parser.add_argument("--max-wgd", type=int, default=1)
    parser.add_argument("--max-cn", type=int, default=10)
    parser.add_argument("--num-processors", type=int, default=1)
    return parser.parse_args()


def stage_inputs(input_dir: Path, output_dir: Path) -> None:
    require_files(input_dir, SEACON_INPUT_FILES)
    for filename in SEACON_INPUT_FILES:
        shutil.copy2(input_dir / filename, output_dir / filename)


def run_command(command: list[str]) -> None:
    print("Running:", " ".join(map(str, command)))
    subprocess.run(command, check=True)


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir

    ensure_clean_dir(output_dir, overwrite=args.overwrite)
    stage_inputs(input_dir, output_dir)

    # Write BAF.tsv directly from precomputed_baf.tsv (flat cell x bin format)
    # This bypasses prep_baf which requires BAM files
    import numpy as np
    rdr = pd.read_csv(output_dir / "RDR.tsv", sep="\t", index_col=0)
    cell_names = rdr.index.tolist()
    n_bins = rdr.shape[1]
    baf_in = pd.read_csv(output_dir / "precomputed_baf.tsv", sep="\t")
    regions = pd.read_csv(output_dir / "filtered_regions.bed", sep="\t",
                          header=None, names=["chrom", "start", "end"])
    baf_mat = pd.DataFrame(0.5, index=cell_names, columns=list(range(n_bins)))
    for _, row in baf_in.iterrows():
        match = regions[(regions["chrom"] == row["chrom"]) &
                        ((regions["start"] == row["start"] + 1) |
                         (regions["start"] == row["start"]))]
        if len(match) == 0:
            continue
        bin_idx = match.index[0]
        if row["cell"] in baf_mat.index:
            baf_mat.loc[row["cell"], bin_idx] = row["BAF"]
    baf_mat.to_csv(output_dir / "BAF.tsv", sep="\t")
    print(f"Wrote BAF.tsv: {baf_mat.shape}")

    call_command = [
        args.seacon_bin,
        "call",
        "-o",
        str(output_dir),
        "--upper-filter",
        str(args.upper_filter),
        "--tolerance",
        str(args.tolerance),
        "--max-wgd",
        str(args.max_wgd),
        "--max-CN",
        str(args.max_cn),
    ]
    run_command(call_command)

    calls_file = output_dir / "calls.tsv"
    calls_flat_file = output_dir / "calls_flat.tsv"
    require_files(output_dir, ["BAF.tsv", "calls.tsv", "calls_flat.tsv"])

    calls = pd.read_csv(calls_file, sep="\t")
    manifest = {
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "seacon_bin": args.seacon_bin,
        "calls_rows": int(len(calls)),
        "calls_file": calls_file.name,
        "calls_flat_file": calls_flat_file.name,
    }
    write_json(output_dir / "run_manifest.json", manifest)

    print(f"SEACON complete. calls.tsv has {len(calls)} rows.")


if __name__ == "__main__":
    main()
