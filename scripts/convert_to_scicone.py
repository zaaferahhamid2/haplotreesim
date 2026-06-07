"""
Convert HaploTreeSim dataset to SCICoNE input format.
SCICoNE takes a cells x bins read count matrix (CSV, no header).
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from dataset_io import load_dataset, ensure_clean_dir, write_json


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir or input_dir / "scicone_input"
    ensure_clean_dir(output_dir, overwrite=args.overwrite)

    print(f"Loading dataset from {input_dir}")
    dataset = load_dataset(input_dir)
    rc = dataset["readcounts"]
    bins = dataset["bins"]
    cells = dataset["cells"]

    n_cells, n_bins = rc.shape
    print(f"  {n_cells} cells x {n_bins} bins")

    # SCICoNE input: CSV no header, rows=cells, cols=bins, integer counts
    rc_int = rc.astype(int)
    rc_int.to_csv(output_dir / "readcounts.csv", header=False, index=False)
    print(f"  Wrote readcounts.csv")

    # Region sizes file: one integer per line = number of bins per region
    # For now use 1 bin per region (SCICoNE will do its own segmentation)
    region_sizes = [1] * n_bins
    with open(output_dir / "region_sizes.txt", "w") as f:
        for s in region_sizes:
            f.write(f"{s}\n")
    print(f"  Wrote region_sizes.txt ({n_bins} regions)")

    # Save cell and bin metadata for evaluator
    cells.to_csv(output_dir / "cells.tsv", sep="\t", index=False)
    bins.to_csv(output_dir / "bins.tsv", sep="\t", index=False)

    manifest = {
        "source_dataset": str(input_dir),
        "format": "scicone-input",
        "n_cells": n_cells,
        "n_bins": n_bins,
        "files": {
            "readcounts": "readcounts.csv",
            "region_sizes": "region_sizes.txt",
            "cells": "cells.tsv",
            "bins": "bins.tsv"
        }
    }
    write_json(output_dir / "manifest.json", manifest)
    print(f"Wrote SCICoNE input to {output_dir}")


if __name__ == "__main__":
    main()
