"""
Convert HaploTreeSim dataset to CONET input format.
CONET needs a corrected counts matrix: bins x cells with RDR values,
plus chr, start, end, width, candidate_brkp columns.
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
    output_dir = args.output_dir or input_dir / "conet_input"
    ensure_clean_dir(output_dir, overwrite=args.overwrite)

    print(f"Loading dataset from {input_dir}")
    dataset = load_dataset(input_dir)
    bins = dataset["bins"]
    cells = dataset["cells"]
    rdr = dataset["rdr"]

    n_bins, n_cells = len(bins), len(cells)
    cell_names = cells["cell"].tolist()
    print(f"  {n_cells} cells x {n_bins} bins")

    # Build corrected counts dataframe: bins x cells
    # CONET expects: chr, start, end, width, candidate_brkp, cell1, cell2, ...
    # chr must be integer (X=23, Y=24)
    chrom_to_int = {f"chr{i}": i for i in range(1, 23)}
    chrom_to_int["chrX"] = 23
    chrom_to_int["chrY"] = 24

    cc_rows = []
    prev_chrom = None
    for _, bin_row in bins.iterrows():
        chrom_str = bin_row["chrom"]
        chrom_int = chrom_to_int.get(chrom_str, int(chrom_str.replace("chr", "")))
        start = int(bin_row["start"])
        end = int(bin_row["end"])
        width = end - start

        # candidate_brkp: 1 if first bin of chromosome
        is_brkp = 1 if chrom_str != prev_chrom else 0
        prev_chrom = chrom_str

        bin_idx = int(bin_row["bin"])
        # RDR values for this bin across all cells
        rdr_vals = rdr.iloc[:, bin_idx].values if bin_idx < rdr.shape[1] else np.ones(n_cells)

        row = {
            "chr": chrom_int,
            "start": start,
            "end": end,
            "width": width,
            "candidate_brkp": is_brkp
        }
        for ci, cell in enumerate(cell_names):
            row[cell] = float(rdr_vals[ci])
        cc_rows.append(row)

    cc_df = pd.DataFrame(cc_rows)
    out_file = output_dir / "corrected_counts.csv"
    cc_df.to_csv(out_file, index=False)
    print(f"  Wrote corrected_counts.csv: {cc_df.shape}")

    # Save metadata
    cells.to_csv(output_dir / "cells.tsv", sep="\t", index=False)
    bins.to_csv(output_dir / "bins.tsv", sep="\t", index=False)
    write_json(output_dir / "manifest.json", {
        "source_dataset": str(input_dir),
        "format": "conet-input",
        "n_cells": n_cells,
        "n_bins": n_bins,
        "files": {"corrected_counts": "corrected_counts.csv"}
    })
    print(f"Wrote CONET input to {output_dir}")


if __name__ == "__main__":
    main()
