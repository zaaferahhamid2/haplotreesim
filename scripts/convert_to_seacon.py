"""
Convert a general HaploTreeSim dataset into SEACON input files.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from dataset_io import ensure_clean_dir, load_dataset, write_json


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="General HaploTreeSim dataset directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="SEACON input directory. Defaults to INPUT_DIR/seacon_input.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace output-dir if it exists.")
    return parser.parse_args()


def matrix_column(matrix: pd.DataFrame, bin_idx: int):
    """Return the column label used by a dataset matrix for a bin index."""
    return str(bin_idx) if str(bin_idx) in matrix.columns else bin_idx


def cell_clone_map(cells: pd.DataFrame) -> dict[str, str]:
    """Return cell-to-clone labels for clone-median BAF imputation."""
    if "clone_assignment" in cells.columns:
        clone_values = cells["clone_assignment"]
    elif "clone" in cells.columns:
        clone_values = cells["clone"]
    else:
        raise ValueError("cells.tsv must contain clone_assignment or clone for BAF imputation")

    return {
        str(cell): str(clone)
        for cell, clone in zip(cells["cell"], clone_values)
    }


def write_precomputed_baf(
    dataset: dict,
    output_dir: Path,
    bins: pd.DataFrame | None = None,
    missing_baf: int | None = None,
) -> tuple[int, int]:
    cells = dataset["cells"]
    if bins is None:
        bins = dataset["bins"]
        missing_baf = 0
    elif missing_baf is None:
        missing_baf = 0

    alt_counts = dataset["allele_alt"]
    ref_counts = dataset["allele_ref"]

    cell_names = cells["cell"].tolist()
    clones_by_cell = cell_clone_map(cells)
    baf_rows = []

    for _, bin_row in bins.iterrows():
        bin_idx = int(bin_row["bin"])
        bin_col = matrix_column(alt_counts, bin_idx)
        observed_by_clone: dict[str, list[float]] = {}

        for cell in cell_names:
            alt = alt_counts.loc[cell, bin_col]
            ref = ref_counts.loc[cell, bin_col]
            total = alt + ref
            if not np.isfinite(total) or total <= 0:
                continue
            raw_baf = alt / total
            baf = float(min(raw_baf, 1 - raw_baf))
            observed_by_clone.setdefault(clones_by_cell[str(cell)], []).append(baf)

        for cell in cell_names:
            alt = alt_counts.loc[cell, bin_col]
            ref = ref_counts.loc[cell, bin_col]
            total = alt + ref
            if not np.isfinite(total) or total <= 0:
                missing_baf += 1
                clone = clones_by_cell[str(cell)]
                clone_observed = observed_by_clone.get(clone, [])
                if not clone_observed:
                    baf = 0.5  # fallback: no depth in clone, assume balanced
                else:
                    baf = round(float(np.median(clone_observed)), 5)
            else:
                raw_baf = alt / total
                baf = round(float(min(raw_baf, 1 - raw_baf)), 5)

            baf_rows.append({
                "chrom": bin_row["chrom"],
                # SEACON adds 1 to precomputed BAF starts before matching filtered_regions.
                "start": int(bin_row["start"]) - 1,
                "end": int(bin_row["end"]),
                "cell": cell,
                "BAF": baf,
            })

    baf_df = pd.DataFrame(baf_rows, columns=["chrom", "start", "end", "cell", "BAF"])
    baf_df.to_csv(output_dir / "precomputed_baf.tsv", sep="\t", index=False)
    return len(baf_rows), missing_baf


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir or input_dir / "seacon_input"
    ensure_clean_dir(output_dir, overwrite=args.overwrite)

    dataset = load_dataset(input_dir)
    cells = dataset["cells"]
    bins = dataset["bins"]

    (output_dir / "cells.txt").write_text("\n".join(cells["cell"].tolist()) + "\n")

    regions = bins[["chrom", "start", "end"]].copy()
    regions.to_csv(output_dir / "filtered_regions.bed", sep="\t", header=False, index=False)

    dataset["readcounts"].to_csv(output_dir / "readcounts.tsv", sep="\t")
    dataset["rdr"].to_csv(output_dir / "RDR.tsv", sep="\t")

    observed_baf, missing_baf = write_precomputed_baf(dataset, output_dir, bins=bins)
    expected_baf_rows = len(cells) * len(bins)
    if observed_baf != expected_baf_rows:
        raise ValueError(
            f"precomputed_baf.tsv has {observed_baf} rows; expected {expected_baf_rows} "
            "for dense SEACON input."
        )

    # Write diploid.txt with normal cell names (clone_assignment == 0)
    normal_cells = cells[cells["clone_assignment"] == 0]["cell"].tolist()
    if normal_cells:
        (output_dir / "diploid.txt").write_text("\n".join(normal_cells) + "\n")
        print(f"  Wrote diploid.txt with {len(normal_cells)} normal cells")

    manifest = {
        "source_dataset": str(input_dir),
        "format": "seacon-input",
        "num_regions": int(len(bins)),
        "baf_imputation": "same_clone_median",
        "files": {
            "cells": "cells.txt",
            "regions": "filtered_regions.bed",
            "readcounts": "readcounts.tsv",
            "rdr": "RDR.tsv",
            "precomputed_baf": "precomputed_baf.tsv",
        },
        "observed_baf_rows": observed_baf,
        "imputed_bin_cell_baf": missing_baf,
    }
    write_json(output_dir / "manifest.json", manifest)

    if np.isnan(dataset["rdr"].to_numpy(dtype=float)).any():
        raise ValueError("RDR.tsv contains NaN values")

    print(f"Wrote SEACON input to {output_dir}")
    print(
        f"  BAF rows: {observed_baf}; "
        f"imputed zero-depth cell-bin BAF values: {missing_baf}"
    )


if __name__ == "__main__":
    main()
