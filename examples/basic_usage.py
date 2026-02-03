"""
Example usage of HaploTreeSim (Week 5 version).

This script demonstrates how to:
1. Configure and run a simulation
2. Access the generated data
3. Visualize basic statistics
"""

import sys
sys.path.insert(0, '/home/claude/haplotreesim/src')

import numpy as np
from haplotreesim import SimulationConfig, HaploTreeSimulator


def main():
    """Run a simple example simulation."""
    
    print("="*70)
    print("HaploTreeSim Example Usage (Week 5)")
    print("="*70)
    
    # Configure simulation
    print("\n1. Creating configuration...")
    config = SimulationConfig(
        # Genome
        num_bins=500,
        bin_length=100000,  # 100kb bins
        
        # Clones
        num_clones=4,
        root_type="diploid",
        
        # Cells
        num_cells=100,
        
        # Coverage (simulating ~0.05x per cell)
        mean_library_size=25.0,
        mean_allelic_coverage=15.0,
        
        # Random seed for reproducibility
        random_seed=2024
    )
    
    print(f"   • {config.num_bins} bins of {config.bin_length/1000:.0f}kb")
    print(f"   • {config.num_clones} clones")
    print(f"   • {config.num_cells} cells")
    print(f"   • ~{config.mean_library_size*2/100:.2f}x mean coverage")
    
    # Run simulation
    print("\n2. Running simulation...")
    sim = HaploTreeSimulator(config)
    read_counts, (alt_counts, ref_counts, total_counts) = sim.run()
    
    # Analyze results
    print("\n3. Analyzing results...")
    
    # Read count statistics
    print("\n   Read Counts:")
    print(f"   • Shape: {read_counts.shape} (cells × bins)")
    print(f"   • Mean per bin: {read_counts.mean():.2f}")
    print(f"   • Std per bin: {read_counts.std():.2f}")
    print(f"   • Total reads: {read_counts.sum():,}")
    
    # Allele count statistics
    print("\n   Allele Counts:")
    print(f"   • Shape: {alt_counts.shape} (cells × segments)")
    print(f"   • Total alternate: {alt_counts.sum():,}")
    print(f"   • Total reference: {ref_counts.sum():,}")
    mean_vaf = alt_counts.sum() / (total_counts.sum() + 1e-10)
    print(f"   • Mean VAF (diploid): {mean_vaf:.3f} (expect ~0.5)")
    
    # Ground truth
    print("\n4. Extracting ground truth...")
    gt = sim.get_ground_truth()
    
    print(f"   • Clone assignments: {gt['clone_assignments'].shape}")
    print(f"   • CN profiles (hap A): {gt['cn_profiles_A'].shape}")
    print(f"   • CN profiles (hap B): {gt['cn_profiles_B'].shape}")
    
    # Clone distribution
    clone_counts = np.bincount(gt['clone_assignments'])
    print("\n   Clone distribution:")
    for k, count in enumerate(clone_counts):
        pct = 100 * count / config.num_cells
        print(f"      Clone {k}: {count:3d} cells ({pct:5.1f}%)")
    
    # Per-cell statistics
    print("\n5. Per-cell statistics:")
    reads_per_cell = read_counts.sum(axis=1)
    print(f"   • Mean reads/cell: {reads_per_cell.mean():.0f}")
    print(f"   • Min reads/cell: {reads_per_cell.min():.0f}")
    print(f"   • Max reads/cell: {reads_per_cell.max():.0f}")
    
    allelic_per_cell = total_counts.sum(axis=1)
    print(f"   • Mean allelic reads/cell: {allelic_per_cell.mean():.0f}")
    print(f"   • Min allelic reads/cell: {allelic_per_cell.min():.0f}")
    print(f"   • Max allelic reads/cell: {allelic_per_cell.max():.0f}")
    
    # Coverage statistics
    print("\n6. Coverage statistics:")
    genome_size = config.num_bins * config.bin_length
    avg_coverage = reads_per_cell.mean() * config.bin_length / genome_size
    print(f"   • Genome size: {genome_size/1e6:.1f} Mb")
    print(f"   • Average coverage: {avg_coverage:.4f}x")
    
    print("\n" + "="*70)
    print("Example complete!")
    print("="*70)
    print("\nNext steps:")
    print("  • Inspect sim.bins, sim.segments, sim.clones, sim.cells")
    print("  • Use read_counts and allele_counts for downstream analysis")
    print("  • Compare inferred CNAs to gt['cn_profiles_A/B']")
    print("\nWeek 6 will add CNA events and tree structure!")
    
    return sim, read_counts, (alt_counts, ref_counts, total_counts)


if __name__ == "__main__":
    sim, reads, alleles = main()
