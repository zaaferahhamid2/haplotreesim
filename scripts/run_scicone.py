"""
Run SCICoNE on converted HaploTreeSim input using pyscicone Python API.
This gives access to cell_node_labels and proper tree objects.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--scicone-build", type=Path,
                        default=Path.home() / "Documents/SCICoNE/build")
    parser.add_argument("--n-reps", type=int, default=4)
    parser.add_argument("--n-iters", type=int, default=4000)
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--threshold", type=float, default=3.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    import scicone

    sci = scicone.SCICoNE(
        str(args.scicone_build) + "/",
        str(args.output_dir) + "/",
        verbose=True,
        postfix=args.input_dir.name
    )

    # Load read counts matrix (cells x bins)
    print(f"Loading data from {args.input_dir}")
    rc = pd.read_csv(args.input_dir / "readcounts.csv", header=None)
    d_mat = rc.values.astype(float)
    n_cells, n_bins = d_mat.shape
    print(f"  {n_cells} cells x {n_bins} bins")

    # Step 1: Breakpoint detection
    print("\n--- Breakpoint detection ---")
    bps = sci.detect_breakpoints(
        d_mat,
        window_size=args.window_size,
        threshold=args.threshold,
        bp_limit=300
    )
    n_regions = len(bps["segmented_region_sizes"])
    print(f"  Detected {n_regions} regions")

    # Save breakpoints
    np.savetxt(args.output_dir / "segmented_regions.txt",
               bps["segmented_regions"], fmt="%d")
    np.savetxt(args.output_dir / "segmented_region_sizes.txt",
               bps["segmented_region_sizes"], fmt="%d")

    # Step 2: Tree inference
    print("\n--- Tree inference ---")
    inferred_tree = sci.learn_tree(
        d_mat,
        bps["segmented_region_sizes"],
        n_reps=args.n_reps,
        seed=args.seed,
        cluster=True,
        full=True,
        max_scoring=True,
        robustness_thr=0.5
    )

    # Save outputs
    cell_node_labels = inferred_tree.cell_node_labels
    pd.Series(cell_node_labels).to_csv(
        args.output_dir / "cell_node_labels.csv", index=False, header=False)
    print(f"  Wrote cell_node_labels.csv")
    print(f"  Unique nodes: {sorted(set(cell_node_labels))}")

    cnvs = inferred_tree.outputs["inferred_cnvs"]
    np.savetxt(args.output_dir / "inferred_cnvs.csv", cnvs, delimiter=",", fmt="%d")
    print(f"  Wrote inferred_cnvs.csv: {cnvs.shape}")

    # Save tree structure from tree_str
    tree_str = getattr(inferred_tree, "tree_str", None)
    if tree_str:
        with open(args.output_dir / "tree_str.txt", "w") as f:
            f.write(tree_str)
        print(f"  Wrote tree_str.txt")
    
    # Parse edges from outputs if available
    edges = []
    tree_edges_file = list(Path(str(args.output_dir)).glob("*tree_inferred.txt"))
    if not tree_edges_file:
        tree_edges_file = list(Path(str(args.scicone_build)).parent.glob("*tree_inferred.txt"))
    if tree_edges_file:
        with open(tree_edges_file[0]) as f:
            for line in f:
                if line.startswith("node ") and "p_id:" in line:
                    node_id = int(line.split(":")[0].replace("node ","").strip())
                    pid_str = line.split("p_id:")[1].split(",")[0]
                    if pid_str != "NULL":
                        edges.append((int(pid_str), node_id))
        with open(args.output_dir / "tree_edges.txt", "w") as f:
            for p, c in edges:
                f.write(f"{p}\t{c}\n")
        print(f"  Wrote tree_edges.txt: {len(edges)} edges")

    # Save manifest
    with open(args.output_dir / "manifest.json", "w") as f:
        json.dump({
            "input_dir": str(args.input_dir),
            "output_dir": str(args.output_dir),
            "n_cells": n_cells,
            "n_bins": n_bins,
            "n_regions": n_regions,
            "unique_nodes": sorted(set(cell_node_labels))
        }, f, indent=2)

    print(f"\nDone. Output in {args.output_dir}")


if __name__ == "__main__":
    main()
