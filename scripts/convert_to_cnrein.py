"""
Convert HaploTreeSim dataset to CNRein input format.
CNRein needs RDR and BAF arrays in numpy format.
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
    output_dir = args.output_dir or input_dir / "cnrein_input"
    ensure_clean_dir(output_dir, overwrite=args.overwrite)

    print(f"Loading dataset from {input_dir}")
    dataset = load_dataset(input_dir)
    bins = dataset["bins"]
    cells = dataset["cells"]
    rdr = dataset["rdr"]
    alt = dataset["allele_alt"]
    ref = dataset["allele_ref"]

    n_cells, n_bins = rdr.shape
    print(f"  {n_cells} cells x {n_bins} bins")

    chrom_to_int = {f'chr{i}': i-1 for i in range(1, 23)}
    chr_arr = np.array([chrom_to_int.get(c, 0) for c in bins['chrom'].values])

    # RDR: (n_cells, n_bins)
    RDR = rdr.values.astype(float)

    # HAP: (n_cells, n_bins, 2)
    alt_arr = alt.values.astype(float)
    ref_arr = ref.values.astype(float)
    total = alt_arr + ref_arr
    total[total == 0] = 1
    baf = alt_arr / total
    HAP = np.stack([baf, 1 - baf], axis=2)

    # Save initial/ files
    initial_dir = output_dir / "initial"
    initial_dir.mkdir(parents=True, exist_ok=True)
    np.savez(initial_dir / "RDR_1M.npz", RDR)
    np.savez(initial_dir / "HAP_1M.npz", HAP)
    np.savez(initial_dir / "chr_1M.npz", chr_arr)
    np.savez(initial_dir / "cellNames.npz", np.array(cells["cell"].tolist()))
    np.savez(initial_dir / "allChr_100k.npz", chr_arr)

    (output_dir / "binScale").mkdir(exist_ok=True)
    (output_dir / "finalPrediction").mkdir(exist_ok=True)

    # Save metadata
    cells.to_csv(output_dir / "cells.tsv", sep="\t", index=False)
    bins.to_csv(output_dir / "bins.tsv", sep="\t", index=False)
    write_json(output_dir / "manifest.json", {
        "source_dataset": str(input_dir),
        "n_cells": n_cells, "n_bins": n_bins
    })
    print(f"Wrote CNRein input to {output_dir}")


if __name__ == "__main__":
    main()
