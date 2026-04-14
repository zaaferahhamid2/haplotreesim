"""
Week 10: BAF/mBAF Plots with Doublet Examples
Shows intermediate BAF values from doublets
"""
import numpy as np
import matplotlib.pyplot as plt
from haplotreesim import SimulationConfig, HaploTreeSimulator

def plot_baf_with_doublets():
    """Show BAF patterns including doublets."""
    
    print("Generating dataset with doublets...")
    config = SimulationConfig(
        chromosome='chr1',
        bin_width=500000,
        num_clones=5,
        num_cells=300,
        lambda_events=3.0,
        prob_doublet=0.15,  # 15% doublets for visibility
        mean_library_size=50000,
        snp_density=1e-3,  # 1 SNP per kb
        random_seed=400
    )
    
    sim = HaploTreeSimulator(config)
    _, (a, r, t) = sim.run()
    
    # Separate singlets and doublets
    singlet_cells = [c for c in sim.cells if not c.is_doublet]
    doublet_cells = [c for c in sim.cells if c.is_doublet]
    
    print(f"Singlets: {len(singlet_cells)}")
    print(f"Doublets: {len(doublet_cells)}")
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Plot 1: Singlets BAF
    ax = axes[0]
    singlet_baf = []
    singlet_tcn = []
    
    for cell in singlet_cells[:100]:  # Sample first 100 singlets
        clone = sim.clones[cell.clone_assignment]
        
        for s, segment in enumerate(sim.segments):
            if t[cell.index, s] > 5:  # Enough coverage
                baf = a[cell.index, s] / t[cell.index, s]
                tcn = clone.get_segment_cn(segment)[2]  # Total CN
                
                singlet_baf.append(baf)
                singlet_tcn.append(tcn)
    
    ax.scatter(singlet_tcn, singlet_baf, alpha=0.3, s=5, c='blue')
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Total Copy Number', fontsize=11)
    ax.set_ylabel('BAF (Alternate Allele Fraction)', fontsize=11)
    ax.set_title('A. Singlet Cells (Normal)', fontsize=12, fontweight='bold')
    ax.set_ylim([0, 1])
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Doublets BAF
    ax = axes[1]
    doublet_baf = []
    doublet_tcn_sum = []
    
    for cell in doublet_cells:
        k1, k2 = cell.doublet_clones
        clone1 = sim.clones[k1]
        clone2 = sim.clones[k2]
        
        for s, segment in enumerate(sim.segments):
            if t[cell.index, s] > 5:
                baf = a[cell.index, s] / t[cell.index, s]
                
                # Total CN is sum from both clones
                tcn1 = clone1.get_segment_cn(segment)[2]
                tcn2 = clone2.get_segment_cn(segment)[2]
                tcn_sum = tcn1 + tcn2
                
                doublet_baf.append(baf)
                doublet_tcn_sum.append(tcn_sum)
    
    ax.scatter(doublet_tcn_sum, doublet_baf, alpha=0.3, s=5, c='red')
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Total Copy Number (Sum)', fontsize=11)
    ax.set_ylabel('BAF (Alternate Allele Fraction)', fontsize=11)
    ax.set_title('B. Doublet Cells', fontsize=12, fontweight='bold')
    ax.set_ylim([0, 1])
    ax.grid(True, alpha=0.3)
    ax.text(0.05, 0.95, 'Note: Intermediate BAF\nfrom CN-weighted mixture', 
            transform=ax.transAxes, fontsize=10, 
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Plot 3: Comparison histogram
    ax = axes[2]
    
    # Convert to mirrored BAF (mBAF)
    singlet_mbaf = [min(b, 1-b) for b in singlet_baf]
    doublet_mbaf = [min(b, 1-b) for b in doublet_baf]
    
    ax.hist(singlet_mbaf, bins=50, alpha=0.5, label='Singlets', color='blue', density=True)
    ax.hist(doublet_mbaf, bins=50, alpha=0.5, label='Doublets', color='red', density=True)
    ax.set_xlabel('Mirrored BAF', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('C. mBAF Distribution', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sample_outputs/week10_baf_with_doublets.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: sample_outputs/week10_baf_with_doublets.png")
    plt.close()

if __name__ == '__main__':
    plot_baf_with_doublets()
