"""
Create CN profile plots across real chromosomes for Week 8.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from haplotreesim import SimulationConfig, HaploTreeSimulator, chromosome_data

print("Creating chr1 CN profile plots...")

# Run simulation
config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=100,
    lambda_events=1.5,
    random_seed=42
)

sim = HaploTreeSimulator(config)
sim.run()

# Get chromosome info
chr_info = chromosome_data.describe_chromosome('chr1', 500000)

print(f"\nSimulated {chr_info['chromosome']}: {chr_info['length_mb']:.1f} Mb")
print(f"Bins: {chr_info['num_bins']}")
print(f"Segments: {len(sim.segments)}")

# Create figure
fig, axes = plt.subplots(len(sim.clones), 1, figsize=(16, 3*len(sim.clones)))
if len(sim.clones) == 1:
    axes = [axes]

# Convert bins to genomic positions (in Mb)
bin_positions = np.arange(config.num_bins) * (config.bin_width / 1e6)

for i, clone in enumerate(sim.clones):
    ax = axes[i]
    
    # Plot CN profiles
    ax.plot(bin_positions, clone.cn_profile_A, label='Haplotype A', 
            linewidth=1.5, alpha=0.8, color='#1f77b4')
    ax.plot(bin_positions, clone.cn_profile_B, label='Haplotype B', 
            linewidth=1.5, alpha=0.8, color='#ff7f0e')
    ax.plot(bin_positions, clone.total_cn(), label='Total CN', 
            linewidth=2, color='black', alpha=0.7)
    
    # Mark centromere
    centro_p = chr_info['centromere_p_end'] / 1e6
    centro_q = chr_info['centromere_q_start'] / 1e6
    ax.axvspan(centro_p, centro_q, alpha=0.1, color='gray', label='Centromere')
    
    # Reference line at diploid
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.3, linewidth=1)
    
    # Labels
    ax.set_ylabel('Copy Number', fontsize=11)
    ax.set_ylim(-0.5, 5)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # Title
    is_leaf = i in sim._leaf_clone_indices
    leaf_str = " [LEAF]" if is_leaf else ""
    parent_str = f"parent: {clone.parent_index}" if clone.parent_index is not None else "root"
    ax.set_title(f"Clone {i} ({parent_str}, {len(clone.events)} events){leaf_str}", 
                 fontsize=12, fontweight='bold')
    
    if i == len(sim.clones) - 1:
        ax.set_xlabel(f'{chr_info["chromosome"]} position (Mb)', fontsize=11)

plt.suptitle(f'{chr_info["chromosome"]} Copy Number Profiles', 
             fontsize=14, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('sample_outputs/chr1_cn_profiles.png', dpi=150, bbox_inches='tight')
print("\n✓ Saved: sample_outputs/chr1_cn_profiles.png")

print("\n✓ All chr1 plots created!")
