"""
Week 9: Read-Depth Sanity Plots with WGD
Shows importance of ploidy normalization
"""
import numpy as np
import matplotlib.pyplot as plt
from haplotreesim import SimulationConfig, HaploTreeSimulator

def plot_depth_vs_cn_with_wgd():
    """Show read depth vs CN with and without WGD."""
    
    print("Generating WGD example with events...")
    config = SimulationConfig(
        chromosome='chr1',
        bin_width=500000,
        num_clones=5,
        num_cells=200,
        lambda_events=3.0,
        prob_wgd=1.0,
        mean_library_size=50000,
        random_seed=300
    )
    
    sim = HaploTreeSimulator(config)
    read_counts, _ = sim.run()
    
    # Find diploid and tetraploid clones
    diploid_clones = [i for i, c in enumerate(sim.clones) if 1.8 < c.ploidy < 2.2]
    tetraploid_clones = [i for i, c in enumerate(sim.clones) if c.ploidy > 3.5]
    
    print(f"Diploid clones: {diploid_clones}")
    print(f"Tetraploid clones: {tetraploid_clones}")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Collect data for diploid clones
    if diploid_clones:
        diploid_cn = []
        diploid_depth = []
        
        for cell in sim.cells:
            if cell.clone_assignment in diploid_clones:
                clone = sim.clones[cell.clone_assignment]
                tcn = clone.cn_profile_A + clone.cn_profile_B
                
                # Sample 500 bins
                sample_indices = np.random.choice(len(tcn), size=min(500, len(tcn)), replace=False)
                
                for idx in sample_indices:
                    diploid_cn.append(tcn[idx])
                    diploid_depth.append(read_counts[cell.index, idx])
        
        # Plot diploid
        ax = axes[0]
        ax.scatter(diploid_cn, diploid_depth, alpha=0.3, s=10, c='blue')
        
        # Trend lines
        for cn in range(0, 6):
            mask = np.array(diploid_cn) == cn
            if np.sum(mask) > 10:
                mean_depth = np.mean(np.array(diploid_depth)[mask])
                ax.axhline(y=mean_depth, xmin=cn/6, xmax=(cn+1)/6, 
                          color='red', linewidth=2, alpha=0.7)
        
        ax.set_xlabel('Copy Number', fontsize=11)
        ax.set_ylabel('Read Depth', fontsize=11)
        ax.set_title('A. Diploid Clones (ploidy ~2.0)', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    # Collect data for tetraploid clones
    if tetraploid_clones:
        tetra_cn = []
        tetra_depth = []
        
        for cell in sim.cells:
            if cell.clone_assignment in tetraploid_clones:
                clone = sim.clones[cell.clone_assignment]
                tcn = clone.cn_profile_A + clone.cn_profile_B
                
                sample_indices = np.random.choice(len(tcn), size=min(500, len(tcn)), replace=False)
                
                for idx in sample_indices:
                    tetra_cn.append(tcn[idx])
                    tetra_depth.append(read_counts[cell.index, idx])
        
        # Plot tetraploid
        ax = axes[1]
        ax.scatter(tetra_cn, tetra_depth, alpha=0.3, s=10, c='green')
        
        # Trend lines
        for cn in range(0, 10):
            mask = np.array(tetra_cn) == cn
            if np.sum(mask) > 10:
                mean_depth = np.mean(np.array(tetra_depth)[mask])
                ax.axhline(y=mean_depth, xmin=cn/10, xmax=(cn+1)/10, 
                          color='red', linewidth=2, alpha=0.7)
        
        ax.set_xlabel('Copy Number', fontsize=11)
        ax.set_ylabel('Read Depth', fontsize=11)
        ax.set_title('B. Tetraploid Clones (ploidy ~4.0, WGD)', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.text(0.05, 0.95, 'Note: Baseline shifted up\ndue to WGD', 
                transform=ax.transAxes, fontsize=10, 
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig('sample_outputs/week9_depth_vs_cn_wgd.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: sample_outputs/week9_depth_vs_cn_wgd.png")
    plt.close()

if __name__ == '__main__':
    plot_depth_vs_cn_with_wgd()
