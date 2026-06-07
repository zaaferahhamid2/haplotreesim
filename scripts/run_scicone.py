"""
Run SCICoNE on converted HaploTreeSim input.
Steps: breakpoint detection -> segment counts -> tree inference
"""

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dataset_io import write_json


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--scicone-build", type=Path,
                        default=Path.home() / "Documents/SCICoNE/build",
                        help="Path to SCICoNE build directory")
    parser.add_argument("--n-iters", type=int, default=4000)
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--threshold", type=float, default=3.0)
    parser.add_argument("--ploidy", type=int, default=2)
    parser.add_argument("--copy-number-limit", type=int, default=4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def run(cmd, cwd=None):
    print("Running:", " ".join(str(c) for c in cmd))
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout[-500:])
    if result.returncode != 0:
        print("STDERR:", result.stderr[-500:])
        raise RuntimeError(f"Command failed: {cmd[0]}")
    return result


def main():
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    build = args.scicone_build
    scripts = build.parent / "scripts"

    import pandas as pd
    rc = pd.read_csv(input_dir / "readcounts.csv", header=None)
    n_cells, n_bins = rc.shape
    print(f"Input: {n_cells} cells x {n_bins} bins")

    rc_file = output_dir / "readcounts.csv"
    import shutil
    shutil.copy2(input_dir / "readcounts.csv", rc_file)

    # Step 1: Breakpoint detection
    print("\n--- Step 1: Breakpoint detection ---")
    run([
        build / "breakpoint_detection",
        f"--d_matrix_file={rc_file}",
        f"--n_bins={n_bins}",
        f"--n_cells={n_cells}",
        f"--window_size={args.window_size}",
        f"--threshold={args.threshold}",
        "--postfix=scicone"
    ], cwd=output_dir)

    # Step 2: Segment counts
    print("\n--- Step 2: Segment counts ---")
    segmented_file = output_dir / "readcounts_segmented_counts.txt"
    region_sizes_file = output_dir / "scicone_segmented_region_sizes.txt"

    # find the segmented regions file
    import glob
    seg_regions = list(output_dir.glob("*scicone_segmented_region_sizes.txt"))
    if not seg_regions:
        seg_regions = list(output_dir.glob("*segmented_region_sizes*"))

    run([
        sys.executable,
        scripts / "segment_counts.py",
        str(rc_file),
        str(output_dir / "scicone_segmented_region_sizes.txt")
    ], cwd=output_dir)

    # find segmented counts file
    seg_counts = list(output_dir.glob("*segmented_counts.txt"))
    if not seg_counts:
        raise FileNotFoundError("Segmented counts file not found")
    seg_counts_file = seg_counts[0]
    seg_sizes_file = output_dir / "scicone_segmented_region_sizes.txt"

    seg_rc = pd.read_csv(seg_counts_file, header=None, sep="\t")
    n_regions = seg_rc.shape[1]
    print(f"  Segmented into {n_regions} regions")

    # Step 3: Tree inference
    print("\n--- Step 3: Tree inference ---")
    run([
        build / "inference",
        f"--d_matrix_file={seg_counts_file}",
        f"--region_sizes_file={seg_sizes_file}",
        f"--n_regions={n_regions}",
        f"--n_cells={n_cells}",
        f"--ploidy={args.ploidy}",
        "--verbosity=1",
        f"--copy_number_limit={args.copy_number_limit}",
        f"--n_iters={args.n_iters}",
        f"--seed={args.seed}",
        "--postfix=scicone"
    ], cwd=output_dir)

    # Check outputs
    cnv_file = output_dir / "scicone_inferred_cnvs.csv"
    tree_file = output_dir / "scicone_tree_inferred.txt"

    if cnv_file.exists():
        cnvs = pd.read_csv(cnv_file, header=None)
        print(f"\n✓ CNV profiles: {cnvs.shape[0]} cells x {cnvs.shape[1]} regions")
    if tree_file.exists():
        print(f"✓ Tree saved to {tree_file}")

    write_json(output_dir / "manifest.json", {
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "n_cells": n_cells,
        "n_bins": n_bins,
        "n_regions": n_regions,
        "files": {
            "cnvs": str(cnv_file),
            "tree": str(tree_file),
            "segmented_counts": str(seg_counts_file),
            "segmented_region_sizes": str(seg_sizes_file)
        }
    })
    print(f"\nDone. Output in {output_dir}")


if __name__ == "__main__":
    main()
