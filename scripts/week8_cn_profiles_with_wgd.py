"""
Week 8: CN Profile Visualization with WGD Examples
"""
import numpy as np
import matplotlib.pyplot as plt
from haplotreesim import SimulationConfig, HaploTreeSimulator

def plot_cn_profiles_comparison():
    """Compare diploid vs WGD CN profiles."""
    
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Simulation 1: No WGD (diploid baseline)
    print("Generating diploid example (no WGD)...")
    config1 = SimulationConfig(
        chromosome='chr1',
        bin_width=500000,
        num_clones=5,
        num_cells=100,
        lambda_events=2.5,
        prob_wgd=0.0,  # No WGD
        random_seed=100
    )
    
    sim1 = HaploTreeSimulator(config1)
    sim1.run()
    
    # Plot diploid profiles
    ax = axes[0]
    for i, clone in enumerate(sim1.clones[:6]):
        if i == 0:  # Root
            continue
        tcn = clone.cn_profile_A + clone.cn_profile_B
        ax.plot(tcn, label=f'Clone {i} (ploidy={clone.ploidy:.2f})', alpha=0.7, linewidth=1.5)
    
    ax.axhline(y=2, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='Diploid baseline')
    ax.set_xlabel('Bin Index', fontsize=11)
    ax.set_ylabel('Total Copy Number', fontsize=11)
    ax.set_title('A. Diploid Baseline (No WGD)', fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim([0, 6])
    ax.grid(True, alpha=0.3)
    
    # Simulation 2: With WGD
    print("Generating WGD example...")
    config2 = SimulationConfig(
        chromosome='chr1',
        bin_width=500000,
        num_clones=5,
        num_cells=100,
        lambda_events=2.5,
        prob_wgd=1.0,  # Force WGD
        random_seed=200
    )
    
    sim2 = HaploTreeSimulator(config2)
    sim2.run()
    
    # Plot WGD profiles
    ax = axes[1]
    for i, clone in enumerate(sim2.clones[:6]):
        if i == 0:  # Root
            continue
        
        tcn = clone.cn_profile_A + clone.cn_profile_B
        is_wgd = clone.ploidy > 3.5
        style = '-' if is_wgd else '--'
        marker = ' [WGD]' if is_wgd else ''
        
        ax.plot(tcn, style, label=f'Clone {i} (ploidy={clone.ploidy:.2f}){marker}', 
                alpha=0.7, linewidth=1.5)
    
    ax.axhline(y=2, color='gray', linestyle='--', linewidth=1, alpha=0.3, label='Diploid baseline')
    ax.axhline(y=4, color='red', linestyle='--', linewidth=1, alpha=0.5, label='Tetraploid baseline (WGD)')
    ax.set_xlabel('Bin Index', fontsize=11)
    ax.set_ylabel('Total Copy Number', fontsize=11)
    ax.set_title(f'B. With WGD (WGD on node {sim2.wgd_node})', fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim([0, 8])
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('sample_outputs/week8_cn_profiles_with_wgd.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: sample_outputs/week8_cn_profiles_with_wgd.png")
    plt.close()

if __name__ == '__main__':
    plot_cn_profiles_comparison()
