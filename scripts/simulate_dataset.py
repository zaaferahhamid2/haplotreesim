"""
Simulate a tool-neutral HaploTreeSim dataset.

The output directory is the stable input for downstream tool adapters such as
SEACON, CHISEL, or future CNA-tree baselines.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmwrite

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from haplotreesim import HaploTreeSimulator, SimulationConfig
from haplotreesim.chromosome_data import AUTOSOMES, WHOLE_GENOME_CHROMOSOMES
from dataset_io import ensure_clean_dir, write_json, write_matrix


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("examples/small_simulation"),
        help="Directory to write the general dataset.",
    )
    parser.add_argument(
        "--tree-structure",
        type=Path,
        default=None,
        help="Optional exported tree structure to reuse while resampling events and observations.",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace output-dir if it exists.")
    parser.add_argument("--chromosome", default="chr1", help="Single chromosome for backward-compatible runs.")
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        default=None,
        help="Chromosome list to simulate, for example: --chromosomes chr1 chr2 chr3.",
    )
    parser.add_argument(
        "--whole-genome",
        action="store_true",
        help="Simulate all autosomes chr1-chr22. Overrides --chromosome/--chromosomes.",
    )
    parser.add_argument(
        "--include-sex-chromosomes",
        action="store_true",
        help="With --whole-genome, include chrX and chrY as well.",
    )
    parser.add_argument("--bin-width", type=int, default=500000)
    parser.add_argument("--num-bins", type=int, default=None, help="Explicit bin count for single-chromosome synthetic runs.")
    parser.add_argument("--num-clones", type=int, default=4)
    parser.add_argument("--num-cells", type=int, default=100)
    parser.add_argument("--lambda-events", type=float, default=5)
    parser.add_argument("--lambda-amplitude", type=float, default=1.0)
    parser.add_argument("--prob-wgd", type=float, default=0.0)
    parser.add_argument("--prob-normal", type=float, default=0.0, help="Fraction of normal (diploid) cells.")
    parser.add_argument("--prob-focal", type=float, default=0.7, help="Probability of focal CNA events.")
    parser.add_argument(
        "--prob-arm-given-broad",
        type=float,
        default=0.75,
        help="Probability of arm-level events conditional on drawing a broad event.",
    )
    parser.add_argument("--snp-density", type=float, default=0.001, help="Usable heterozygous SNP density per bp.")
    parser.add_argument("--mean-allelic-coverage", type=float, default=None, help="Mean per-SNP allelic coverage. Defaults to simulator auto-calibration.")
    parser.add_argument("--random-seed", type=int, default=42)
    return parser.parse_args()


def event_rows(sim: HaploTreeSimulator) -> list[dict]:
    rows = []
    for clone in sim.clones:
        for event_order, event in enumerate(clone.events):
            rows.append({
                "clone": f"clone{clone.index}",
                "clone_index": clone.index,
                "parent_clone_index": clone.parent_index,
                "event_order": event_order,
                "chrom": (
                    "whole_genome"
                    if event.haplotype.value == "WGD"
                    else sim.bins[event.start_bin].chromosome
                ),
                "start_bin": event.start_bin,
                "end_bin": event.end_bin,
                "haplotype": event.haplotype.value,
                "amplitude": event.amplitude,
                "event_time": event.event_time,
                "scale_class": event.scale_class,
                "event_id": event.event_id,
            })
    return rows


def aggregate_bin_counts_to_segments(counts: np.ndarray, segments) -> np.ndarray:
    segment_counts = np.zeros((counts.shape[0], len(segments)), dtype=counts.dtype)
    for segment in segments:
        segment_counts[:, segment.index] = counts[:, segment.start_bin:segment.end_bin + 1].sum(axis=1)
    return segment_counts


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    ensure_clean_dir(output_dir, overwrite=args.overwrite)

    if args.whole_genome:
        chromosomes = WHOLE_GENOME_CHROMOSOMES if args.include_sex_chromosomes else AUTOSOMES
    elif args.chromosomes is not None:
        chromosomes = args.chromosomes
    else:
        chromosomes = [args.chromosome]

    config = SimulationConfig(
        chromosome=chromosomes[0],
        chromosomes=chromosomes,
        bin_width=args.bin_width,
        num_bins=args.num_bins,
        num_clones=args.num_clones,
        num_cells=args.num_cells,
        lambda_events=args.lambda_events,
        lambda_amplitude=args.lambda_amplitude,
        prob_wgd=args.prob_wgd,
        prob_normal=args.prob_normal,
        prob_focal=args.prob_focal,
        prob_arm_given_broad=args.prob_arm_given_broad,
        snp_density=args.snp_density,
        mean_allelic_coverage=args.mean_allelic_coverage,
        random_seed=args.random_seed,
    )

    print("Simulating dataset...")
    sim = HaploTreeSimulator(config)
    read_counts, allele_counts = sim.run(
        tree_structure_path=str(args.tree_structure) if args.tree_structure else None
    )
    alt_counts, ref_counts, total_counts = allele_counts
    truth = sim.get_ground_truth()

    n_cells, n_bins = read_counts.shape
    cell_names = [f"cell{i}" for i in range(n_cells)]
    bin_columns = list(range(n_bins))
    segment_columns = list(range(len(sim.segments)))
    clone_names = [f"clone{clone.index}" for clone in sim.clones]
    segment_alt_counts = aggregate_bin_counts_to_segments(alt_counts, sim.segments)
    segment_ref_counts = aggregate_bin_counts_to_segments(ref_counts, sim.segments)
    segment_total_counts = aggregate_bin_counts_to_segments(total_counts, sim.segments)

    write_json(output_dir / "config.json", asdict(config))

    bins = pd.DataFrame([{
        "bin": bin_obj.index,
        "chrom": bin_obj.chromosome,
        "start": bin_obj.start + 1,
        "end": bin_obj.end,
        "length": bin_obj.length,
    } for bin_obj in sim.bins])
    bins.to_csv(output_dir / "bins.tsv", sep="\t", index=False)

    segments = pd.DataFrame([{
        "segment": segment.index,
        "chrom": sim.bins[segment.start_bin].chromosome,
        "start_bin": segment.start_bin,
        "end_bin": segment.end_bin,
        "length_bp": segment.length,
        "haplotype_block": segment.haplotype_block,
    } for segment in sim.segments])
    segments.to_csv(output_dir / "segments.tsv", sep="\t", index=False)

    cell_rows = []
    for cell in sim.cells:
        cell_rows.append({
            "cell": f"cell{cell.index}",
            "cell_index": cell.index,
            "clone_assignment": cell.clone_assignment,
            "clone": f"clone{cell.clone_assignment}",
            "library_size": cell.library_size,
            "allelic_coverage": cell.allelic_coverage,
            "is_doublet": cell.is_doublet,
            "doublet_clones": "" if cell.doublet_clones is None else ",".join(map(str, cell.doublet_clones)),
        })
    pd.DataFrame(cell_rows).to_csv(output_dir / "cells.tsv", sep="\t", index=False)

    write_matrix(output_dir / "readcounts.tsv", read_counts, cell_names, bin_columns)
    cell_avgs = pd.DataFrame(read_counts, index=cell_names, columns=bin_columns).mean(axis=1).replace(0, 1)
    rdr = pd.DataFrame(read_counts, index=cell_names, columns=bin_columns).div(cell_avgs, axis=0).round(5)
    rdr.to_csv(output_dir / "rdr.tsv", sep="\t")

    write_matrix(output_dir / "allele_alt.tsv", alt_counts, cell_names, bin_columns)
    write_matrix(output_dir / "allele_ref.tsv", ref_counts, cell_names, bin_columns)
    write_matrix(output_dir / "allele_total.tsv", total_counts, cell_names, bin_columns)
    write_matrix(output_dir / "segment_allele_alt.tsv", segment_alt_counts, cell_names, segment_columns)
    write_matrix(output_dir / "segment_allele_ref.tsv", segment_ref_counts, cell_names, segment_columns)
    write_matrix(output_dir / "segment_allele_total.tsv", segment_total_counts, cell_names, segment_columns)
    write_matrix(output_dir / "clone_cn_A.tsv", truth["cn_profiles_A"], clone_names, bin_columns)
    write_matrix(output_dir / "clone_cn_B.tsv", truth["cn_profiles_B"], clone_names, bin_columns)

    snps = pd.DataFrame({
        "snp": np.arange(len(sim.snp_bins), dtype=int),
        "chrom": sim.snp_chromosomes,
        "position": sim.snp_positions,
        "bin": sim.snp_bins,
        "bin_local_snp": sim.snp_bin_local_indices,
    })
    snps.to_csv(output_dir / "snps.tsv", sep="\t", index=False)
    pd.DataFrame({
        "bin": np.arange(len(sim.bin_snp_counts), dtype=int),
        "snp_count": sim.bin_snp_counts,
    }).to_csv(output_dir / "bin_snp_counts.tsv", sep="\t", index=False)
    mmwrite(output_dir / "snp_allele_alt.mtx", sim.snp_alt_counts)
    mmwrite(output_dir / "snp_allele_ref.mtx", sim.snp_ref_counts)

    clone_rows = []
    for clone in sim.clones:
        clone_rows.append({
            "clone": f"clone{clone.index}",
            "clone_index": clone.index,
            "parent_clone_index": clone.parent_index,
            "is_root": clone.is_root,
            "ploidy": clone.ploidy,
        })
    pd.DataFrame(clone_rows).to_csv(output_dir / "clones.tsv", sep="\t", index=False)

    truth_rows = []
    clone_labels = np.array(truth["clone_assignments"])
    for cell_idx, clone_idx in enumerate(clone_labels):
        clone = sim.clones[clone_idx]
        for segment in sim.segments:
            bins_in_segment = slice(segment.start_bin, segment.end_bin + 1)
            cn_A = int(round(float(np.mean(clone.cn_profile_A[bins_in_segment]))))
            cn_B = int(round(float(np.mean(clone.cn_profile_B[bins_in_segment]))))
            truth_rows.append({
                "cell": f"cell{cell_idx}",
                "cell_index": cell_idx,
                "clone_assignment": clone_idx,
                "segment": segment.index,
                "start_bin": segment.start_bin,
                "end_bin": segment.end_bin,
                "cn_A": cn_A,
                "cn_B": cn_B,
                "total_cn": cn_A + cn_B,
                "loh": min(cn_A, cn_B) == 0 and cn_A + cn_B >= 1,
            })
    pd.DataFrame(truth_rows).to_csv(
        output_dir / "truth_cell_hscn_segments.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(event_rows(sim)).to_csv(output_dir / "events.tsv", sep="\t", index=False)
    sim.export_tree_structure(str(output_dir / "tree_structure.json"))

    metadata = {
        "schema_version": 1,
        "format": "haplotreesim-general-dataset",
        "description": "Tool-neutral simulated low-pass scDNA-seq dataset.",
        "num_cells": n_cells,
        "num_bins": n_bins,
        "num_segments": len(sim.segments),
        "num_clones_total": len(sim.clones),
        "num_clone_leaves": sim.config.num_clones,
        "chromosomes": sim.config.chromosomes,
        "chromosome_bin_offsets": {
            chrom: [int(start), int(end)]
            for chrom, (start, end) in sim.chromosome_bin_offsets.items()
        },
        "interchromosome_breakpoints": [int(bp) for bp in sim.chromosome_boundary_bins],
        "allele_count_unit": "snp-derived-bin",
        "snp_count_unit": "snp",
        "coordinate_system": "1-based inclusive starts/ends in bins.tsv",
        "source_tree_structure": str(args.tree_structure) if args.tree_structure else None,
        "files": {
            "bins": "bins.tsv",
            "segments": "segments.tsv",
            "cells": "cells.tsv",
            "readcounts": "readcounts.tsv",
            "rdr": "rdr.tsv",
            "allele_alt": "allele_alt.tsv",
            "allele_ref": "allele_ref.tsv",
            "allele_total": "allele_total.tsv",
            "snps": "snps.tsv",
            "bin_snp_counts": "bin_snp_counts.tsv",
            "snp_allele_alt": "snp_allele_alt.mtx",
            "snp_allele_ref": "snp_allele_ref.mtx",
            "segment_allele_alt": "segment_allele_alt.tsv",
            "segment_allele_ref": "segment_allele_ref.tsv",
            "segment_allele_total": "segment_allele_total.tsv",
            "clone_cn_A": "clone_cn_A.tsv",
            "clone_cn_B": "clone_cn_B.tsv",
            "cell_hscn_truth": "truth_cell_hscn_segments.tsv",
            "events": "events.tsv",
            "tree": "tree_structure.json",
        },
    }
    write_json(output_dir / "metadata.json", metadata)

    print(f"Wrote general dataset to {output_dir}")


if __name__ == "__main__":
    main()
