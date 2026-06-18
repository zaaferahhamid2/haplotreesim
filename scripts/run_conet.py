"""
Run CONET on converted HaploTreeSim input.
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--conet-bin", type=Path,
                        default=Path.home() / "Documents/CONET/src/CONET")
    parser.add_argument("--param-iters", type=int, default=3000)
    parser.add_argument("--pt-iters", type=int, default=20000)
    parser.add_argument("--counts-penalty", type=float, default=200000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--mixture-size", type=int, default=2)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Need to run from CONET env
    conet_py = Path.home() / "Documents/CONET/python/conet-py"
    sys.path.insert(0, str(conet_py))

    from conet import CorrectedCounts, DataConverter, CONET, CONETParameters, InferenceResult

    print(f"Loading corrected counts from {args.input_dir}")
    cc_file = args.input_dir / "corrected_counts.csv"
    df = pd.read_csv(cc_file, sep=',', header=0, low_memory=False)
    print(f"  Shape: {df.shape}")

    cc = CorrectedCounts(df)
    cc.add_chromosome_ends(neutral_cn=2, end_bin_length=500000)
    print("  CorrectedCounts created")

    # Write CONET input files to output dir
    data_dir = str(args.output_dir) + "/"
    DataConverter(event_length_normalizer=3095677412).create_CoNET_input_files(
        data_dir, corrected_counts=cc)
    print("  Input files written")

    # Run CONET
    print(f"\nRunning CONET (param_iters={args.param_iters}, pt_iters={args.pt_iters})...")
    output_dir = str(args.output_dir / "conet_output") + "/"
    os.makedirs(output_dir, exist_ok=True)

    conet = CONET(str(args.conet_bin), output_path=output_dir)
    params = CONETParameters(
        tree_structure_prior_k1=0.5,
        data_dir=data_dir,
        counts_penalty_s1=args.counts_penalty,
        counts_penalty_s2=args.counts_penalty,
        param_inf_iters=args.param_iters,
        seed=args.seed,
        mixture_size=args.mixture_size,
        pt_inf_iters=args.pt_iters,
        neutral_cn=2,
        output_dir=output_dir
    )
    conet.infer_tree(params)
    print("CONET inference done!")

    # Load results
    result = InferenceResult(output_dir, cc)
    inferred_cn = result.get_inferred_copy_numbers(neutral_cn=2)
    cn_df = pd.DataFrame(inferred_cn)
    cn_df.to_csv(args.output_dir / "inferred_cn.csv", index=False)
    print(f"  Wrote inferred_cn.csv: {cn_df.shape}")

    # Save tree edges
    tree = result.inferred_tree
    edges = list(tree.edges())
    with open(args.output_dir / "inferred_tree_edges.txt", "w") as f:
        for p, c in edges:
            f.write(f"{p}\t{c}\n")
    print(f"  Wrote inferred_tree_edges.txt: {len(edges)} edges")
    print(f"\nDone. Output in {args.output_dir}")


if __name__ == "__main__":
    main()
