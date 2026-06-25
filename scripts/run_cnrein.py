"""
Run CNRein (CNNaive step) on converted HaploTreeSim input.
Bypasses BAM processing by generating required numpy files from simulator data.
"""

import argparse
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-ploidy", type=int, default=10)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "binScale").mkdir(exist_ok=True)
    (out_dir / "finalPrediction").mkdir(exist_ok=True)
    (out_dir / "initial").mkdir(exist_ok=True)

    # Copy initial files from input
    import shutil
    for f in (args.input_dir / "initial").glob("*.npz"):
        shutil.copy2(f, out_dir / "initial" / f.name)

    bins_df = pd.read_csv(args.input_dir / "bins.tsv", sep="\t")
    cells_df = pd.read_csv(args.input_dir / "cells.tsv", sep="\t")
    n_bins = len(bins_df)
    n_cells = len(cells_df)
    cell_names = cells_df["cell"].tolist()

    chrom_to_int = {f'chr{i}': i-1 for i in range(1, 23)}
    chr_arr = np.array([chrom_to_int.get(c, 0) for c in bins_df["chrom"].values])

    # Run CNNaive via patched pipeline
    sys.path.insert(0, str(Path(__file__).parent))
    import importlib.util, site
    sp = site.getsitepackages()[0]
    scaler_path = Path(sp) / "CNRein" / "scaler.py"

    # Patch off-by-one bug in scaler.py
    content = scaler_path.read_text()
    old = "    for a in range(RDR_all.shape[0]):\n        errorList = divideError[a]"
    new = "    for a in range(divideError.shape[0]):\n        errorList = divideError[a]"
    if old in content:
        scaler_path.write_text(content.replace(old, new))
        print("  Patched scaler.py off-by-one bug")

    from CNRein.scaler import runNaiveCopy, findRegions, findDividers
    print("Running CNNaive Step 1/3: Finding regions...")
    chr_file = str(out_dir / "initial/chr_1M.npz")
    RDR_file = str(out_dir / "initial/RDR_1M.npz")
    HAP_file = str(out_dir / "initial/HAP_1M.npz")
    region_file = str(out_dir / "binScale/regions.npz")
    findRegions(RDR_file, HAP_file, chr_file, region_file)

    print("Running CNNaive Step 2/3: Finding scaling factors...")
    divider_file = str(out_dir / "binScale/dividers.npz")
    error_file = str(out_dir / "binScale/dividerError.npz")
    dividerList_file = str(out_dir / "binScale/dividerAll.npz")
    findDividers(RDR_file, HAP_file, chr_file, divider_file, error_file,
                 dividerList_file, region_file, maxPloidy=args.max_ploidy)

    print("Running CNNaive Step 3/3: Estimating CN profiles...")
    # Generate binScale files from regions
    regions = np.load(region_file)["arr_0"]
    n_regions_found = regions.shape[0]
    RDR = np.load(RDR_file)["arr_0"]
    HAP = np.load(HAP_file)["arr_0"]

    filtered_RDR = np.zeros((n_cells, n_regions_found))
    filtered_HAP = np.zeros((n_cells, n_regions_found, 2))
    filtered_RDR_noise = np.zeros((n_cells, n_regions_found))
    BAF_noise = np.zeros((n_cells, n_regions_found))
    chr_per_region = np.zeros(n_regions_found, dtype=int)

    for ri in range(n_regions_found):
        start, end = int(regions[ri, 1]), int(regions[ri, 2])
        end = min(end, n_bins)
        if start >= n_bins:
            continue
        filtered_RDR[:, ri] = np.mean(RDR[:, start:end], axis=1)
        filtered_RDR_noise[:, ri] = np.std(RDR[:, start:end], axis=1)
        filtered_HAP[:, ri, 0] = np.mean(HAP[:, start:end, 0], axis=1)
        filtered_HAP[:, ri, 1] = np.mean(HAP[:, start:end, 1], axis=1)
        BAF_noise[:, ri] = np.std(HAP[:, start:end, 0], axis=1)
        chr_per_region[ri] = chr_arr[start] if start < n_bins else 0

    # Fix NaNs
    filtered_RDR = np.nan_to_num(filtered_RDR, nan=1.0)
    filtered_HAP = np.nan_to_num(filtered_HAP, nan=0.5)
    filtered_RDR_noise = np.nan_to_num(filtered_RDR_noise, nan=0.0)
    BAF_noise = np.nan_to_num(BAF_noise, nan=0.0)

    bs = out_dir / "binScale"
    np.savez(bs / "filtered_RDR_avg.npz", filtered_RDR)
    np.savez(bs / "filtered_HAP_avg.npz", filtered_HAP)
    np.savez(bs / "filtered_RDR_noise.npz", filtered_RDR_noise)
    np.savez(bs / "BAF_noise.npz", BAF_noise)
    np.savez(bs / "chr_avg.npz", chr_per_region)

    from CNRein.scaler import findInitialCNA
    initialCNA_file = str(bs / "initialCNA.npz")
    initialUniqueCNA_file = str(bs / "initialUniqueCNA.npz")
    initialUniqueIndex_file = str(bs / "initialIndex.npz")
    findInitialCNA(
        str(bs / "filtered_RDR_avg.npz"),
        str(bs / "filtered_RDR_noise.npz"),
        str(bs / "filtered_HAP_avg.npz"),
        str(bs / "BAF_noise.npz"),
        str(bs / "chr_avg.npz"),
        str(bs / "dividers.npz"),
        str(bs / "dividerError.npz"),
        str(bs / "dividerAll.npz"),
        initialCNA_file,
        initialUniqueCNA_file,
        initialUniqueIndex_file
    )

    pred = np.load(initialCNA_file)["arr_0"]
    print(f"  initialCNA shape: {pred.shape}")

    # Write CSV output
    n_regions_pred = pred.shape[1]
    bins_per_region = n_bins // n_regions_pred
    rows = []
    for ci in range(n_cells):
        for ri in range(n_regions_pred):
            bin_start = ri * bins_per_region
            bin_end = min((ri + 1) * bins_per_region - 1, n_bins - 1)
            chrom = bins_df.iloc[bin_start]["chrom"].replace("chr", "")
            start = int(bins_df.iloc[bin_start]["start"])
            end = int(bins_df.iloc[bin_end]["end"])
            h1 = int(pred[ci, ri, 0])
            h2 = int(pred[ci, ri, 1])
            rows.append([cell_names[ci], chrom, start, end, h1, h2])

    out_file = out_dir / "finalPrediction/CNNaivePrediction.csv"
    pd.DataFrame(rows, columns=["Cell barcode", "Chromosome", "Start", "End",
                                 "Haplotype 1", "Haplotype 2"]).to_csv(out_file, index=False)
    print(f"  Wrote {out_file}: {len(rows)} rows")
    print(f"  TCN unique: {sorted(set((pred[:,:,0]+pred[:,:,1]).flatten().tolist()))}")
    print(f"\nDone. Output in {out_dir}")


if __name__ == "__main__":
    main()
