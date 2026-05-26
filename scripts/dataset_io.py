"""
Shared helpers for file-based simulation datasets.

The dataset directory is intentionally tool-neutral. Tool-specific scripts
should read these files and write their own input/output folders.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd


REQUIRED_DATASET_FILES = [
    "metadata.json",
    "config.json",
    "bins.tsv",
    "segments.tsv",
    "cells.tsv",
    "readcounts.tsv",
    "rdr.tsv",
    "allele_alt.tsv",
    "allele_ref.tsv",
    "allele_total.tsv",
    "clone_cn_A.tsv",
    "clone_cn_B.tsv",
    "truth_cell_hscn_segments.tsv",
]

OPTIONAL_DATASET_FILES = {
    "segment_allele_alt": "segment_allele_alt.tsv",
    "segment_allele_ref": "segment_allele_ref.tsv",
    "segment_allele_total": "segment_allele_total.tsv",
}


def ensure_clean_dir(path: Path, overwrite: bool = False) -> None:
    """Create an output directory, optionally replacing existing contents."""
    import shutil

    if path.exists():
        if not overwrite:
            raise FileExistsError(f"{path} already exists; pass --overwrite to replace it")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


def require_files(directory: Path, filenames: Iterable[str]) -> None:
    """Raise FileNotFoundError if any required files are missing."""
    missing = [name for name in filenames if not (directory / name).exists()]
    if missing:
        joined = ", ".join(missing)
        raise FileNotFoundError(f"Missing required files in {directory}: {joined}")


def write_json(path: Path, payload: Dict) -> None:
    """Write JSON with stable indentation."""
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def read_json(path: Path) -> Dict:
    """Read JSON from disk."""
    with open(path) as handle:
        return json.load(handle)


def read_matrix(path: Path) -> pd.DataFrame:
    """Read a tab-separated matrix with row IDs in the first column."""
    return pd.read_csv(path, sep="\t", index_col=0)


def write_matrix(path: Path, data: np.ndarray, index: List[str], columns: List[int]) -> None:
    """Write a numeric matrix with stable row and column labels."""
    df = pd.DataFrame(data, index=index, columns=columns)
    df.index.name = None
    df.to_csv(path, sep="\t")


def load_dataset(dataset_dir: Path) -> Dict:
    """Load a general HaploTreeSim dataset directory."""
    require_files(dataset_dir, REQUIRED_DATASET_FILES)

    cells = pd.read_csv(dataset_dir / "cells.tsv", sep="\t")
    bins = pd.read_csv(dataset_dir / "bins.tsv", sep="\t")
    segments = pd.read_csv(dataset_dir / "segments.tsv", sep="\t")

    dataset = {
        "metadata": read_json(dataset_dir / "metadata.json"),
        "config": read_json(dataset_dir / "config.json"),
        "cells": cells,
        "bins": bins,
        "segments": segments,
        "readcounts": read_matrix(dataset_dir / "readcounts.tsv"),
        "rdr": read_matrix(dataset_dir / "rdr.tsv"),
        "allele_alt": read_matrix(dataset_dir / "allele_alt.tsv"),
        "allele_ref": read_matrix(dataset_dir / "allele_ref.tsv"),
        "allele_total": read_matrix(dataset_dir / "allele_total.tsv"),
        "clone_cn_A": read_matrix(dataset_dir / "clone_cn_A.tsv"),
        "clone_cn_B": read_matrix(dataset_dir / "clone_cn_B.tsv"),
        "truth_cell_hscn_segments": pd.read_csv(
            dataset_dir / "truth_cell_hscn_segments.tsv",
            sep="\t"
        ),
    }

    for key, filename in OPTIONAL_DATASET_FILES.items():
        path = dataset_dir / filename
        if path.exists():
            dataset[key] = read_matrix(path)

    return dataset
