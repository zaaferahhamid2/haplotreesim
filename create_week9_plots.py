"""
Week 9 Sanity Plots: Read-depth observation model validation.

Creates plots showing:
1. Mean vs variance by copy number
2. Read depth vs coverage (library size)
3. Distribution of library sizes
4. Bin bias distribution (if enabled)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from haplotreesim import SimulationConfig, HaploTreeSimulator

print("Creating Week 9 Sanity Plots...")
print("="*60)

# Run simulation WITHOUT bin bias
print("\nSimulation 1: No bin bias")
config1 = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=200,  # More cells for better statistics
    lambda_events=1.5,
    use_gc_bias=False,
    mean_library_size=100.0,
    library_size_cv=0.3,
    theta_x=10.0,
    random_seed=42
)

sim1 = HaploTreeSimulator(config1)
read_counts1, _ = sim1.run()

# Run simulation WITH bin bias
print("\nSimulation 2: With bin bias")
config2 = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=200,
    lambda_events=1.5,
    use_gc_bias=True,
    mean_library_size=100.0,
    library_size_cv=0.3,
    theta_x=10.0,
    random_seed=43
)

sim2 = HaploTreeSimulator(config2)
read_counts2, _ = sim2.run()

print("\n" + "="*60)
print("Creating Plots...")
print("="*60)

# Create figure with subplots
fig = plt.figure(figsize=(16, 12))

# Plot 1: Mean vs Variance by Copy Number
print("\n1. Mean vs Variance by CN...")
ax1 = plt.subplot(2, 3, 1)

# Collect mean and variance by CN for sim1
cn_to_stats = {}
for cell in sim1.cells:
    clone = sim1.clones[cell.clone_assignment]
    for b in range(config1.num_bins):
        tcn = clone.total_cn()[b]
        count = read_counts1[cell.index, b]
        
        if tcn not in cn_to_stats:
            cn_to_stats[tcn] = []
        cn_to_stats[tcn].append(count)

cns = sorted(cn_to_stats.keys())
means = [np.mean(cn_to_stats[cn]) for cn in cns]
variances = [np.var(cn_to_stats[cn]) for cn in cns]

ax1.scatter(means, variances, s=100, alpha=0.7, c=cns, cmap='viridis')
ax1.plot([0, max(means)], [0, max(means)], 'r--', alpha=0.5, label='Variance = Mean (Poisson)')

# Add theoretical Negative Binomial variance
theory_var = [m + m**2/config1.theta_x for m in means]
ax1.plot(means, theory_var, 'g--', alpha=0.7, label=f'NegBin (θ={config1.theta_x})')

ax1.set_xlabel('Mean Read Count', fontsize=11)
ax1.set_ylabel('Variance', fontsize=11)
ax1.set_title('Mean-Variance Relationship', fontsize=12, fontweight='bold')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Read Depth vs Library Size
print("2. Depth vs Library Size...")
ax2 = plt.subplot(2, 3, 2)

total_reads_per_cell = read_counts1.sum(axis=1)
library_sizes = [cell.library_size for cell in sim1.cells]

ax2.scatter(library_sizes, total_reads_per_cell, alpha=0.6, s=50)
ax2.set_xlabel('Library Size Factor (α)', fontsize=11)
ax2.set_ylabel('Total Reads per Cell', fontsize=11)
ax2.set_title('Read Depth vs Coverage', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)

# Add correlation
corr = np.corrcoef(library_sizes, total_reads_per_cell)[0, 1]
ax2.text(0.05, 0.95, f'Correlation: {corr:.3f}', 
         transform=ax2.transAxes, va='top', fontsize=10,
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Plot 3: Library Size Distribution
print("3. Library Size Distribution...")
ax3 = plt.subplot(2, 3, 3)

ax3.hist(library_sizes, bins=30, alpha=0.7, edgecolor='black')
ax3.axvline(config1.mean_library_size, color='r', linestyle='--', 
            linewidth=2, label=f'Mean={config1.mean_library_size}')
ax3.axvline(np.mean(library_sizes), color='g', linestyle='--', 
            linewidth=2, label=f'Observed={np.mean(library_sizes):.1f}')

ax3.set_xlabel('Library Size Factor (α)', fontsize=11)
ax3.set_ylabel('Number of Cells', fontsize=11)
ax3.set_title(f'Library Size Distribution (CV={config1.library_size_cv})', 
              fontsize=12, fontweight='bold')
ax3.legend()
ax3.grid(True, alpha=0.3, axis='y')

# Plot 4: Read Count Distribution by CN
print("4. Read Count Distribution by CN...")
ax4 = plt.subplot(2, 3, 4)

# Show distributions for different CNs
for cn in [0, 1, 2, 3, 4]:
    if cn in cn_to_stats and len(cn_to_stats[cn]) > 10:
        counts = cn_to_stats[cn]
        ax4.hist(counts, bins=30, alpha=0.5, label=f'CN={cn}', density=True)

ax4.set_xlabel('Read Count', fontsize=11)
ax4.set_ylabel('Density', fontsize=11)
ax4.set_title('Read Count Distributions by CN', fontsize=12, fontweight='bold')
ax4.legend()
ax4.grid(True, alpha=0.3, axis='y')

# Plot 5: Bin Bias Effect (if enabled)
print("5. Bin Bias Comparison...")
ax5 = plt.subplot(2, 3, 5)

if hasattr(sim2, 'bin_biases'):
    ax5.hist(sim2.bin_biases, bins=50, alpha=0.7, edgecolor='black', color='orange')
    ax5.axvline(1.0, color='r', linestyle='--', linewidth=2, label='Mean=1.0')
    ax5.axvline(np.mean(sim2.bin_biases), color='g', linestyle='--', 
                linewidth=2, label=f'Observed={np.mean(sim2.bin_biases):.3f}')
    ax5.set_xlabel('Bin Bias Factor (κ)', fontsize=11)
    ax5.set_ylabel('Number of Bins', fontsize=11)
    ax5.set_title('Bin Bias Distribution (GC/Mappability)', fontsize=12, fontweight='bold')
    ax5.legend()
    ax5.grid(True, alpha=0.3, axis='y')
else:
    ax5.text(0.5, 0.5, 'Bin bias not enabled', ha='center', va='center',
             transform=ax5.transAxes, fontsize=12)
    ax5.set_title('Bin Bias: Disabled', fontsize=12, fontweight='bold')

# Plot 6: Effect of Bin Bias on Read Counts
print("6. Bin Bias Effect on Reads...")
ax6 = plt.subplot(2, 3, 6)

# Compare mean reads per bin with and without bias
mean_reads_nobias = read_counts1.mean(axis=0)
mean_reads_withbias = read_counts2.mean(axis=0)

if hasattr(sim2, 'bin_biases'):
    # Normalize by expected (remove CN effect)
    # Get average CN per bin across all cells
    avg_cn_per_bin = np.zeros(config1.num_bins)
    for cell in sim1.cells:
        clone = sim1.clones[cell.clone_assignment]
        avg_cn_per_bin += clone.total_cn()
    avg_cn_per_bin /= len(sim1.cells)
    
    # Normalize reads by CN
    norm_nobias = mean_reads_nobias / (avg_cn_per_bin + 0.1)
    norm_withbias = mean_reads_withbias / (avg_cn_per_bin + 0.1)
    
    ax6.scatter(sim2.bin_biases, norm_withbias / norm_nobias.mean(), 
                alpha=0.5, s=20)
    ax6.plot([0.5, 1.5], [0.5, 1.5], 'r--', alpha=0.7, label='y=x')
    ax6.set_xlabel('Bin Bias Factor (κ)', fontsize=11)
    ax6.set_ylabel('Relative Read Depth', fontsize=11)
    ax6.set_title('Bin Bias Effect on Read Counts', fontsize=12, fontweight='bold')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
else:
    ax6.text(0.5, 0.5, 'No bias comparison', ha='center', va='center',
             transform=ax6.transAxes, fontsize=12)

plt.suptitle('Week 9: Read-Depth Model Sanity Checks', 
             fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('sample_outputs/week9_sanity_plots.png', dpi=150, bbox_inches='tight')
print("\n✓ Saved: sample_outputs/week9_sanity_plots.png")

print("\n" + "="*60)
print("Summary Statistics")
print("="*60)
print(f"\nLibrary Sizes:")
print(f"  Target mean: {config1.mean_library_size}")
print(f"  Observed mean: {np.mean(library_sizes):.2f}")
print(f"  Target CV: {config1.library_size_cv}")
print(f"  Observed CV: {np.std(library_sizes)/np.mean(library_sizes):.3f}")

print(f"\nNegative Binomial Overdispersion:")
print(f"  θ parameter: {config1.theta_x}")
print(f"  Expected variance/mean ratio: ~{1 + 100/config1.theta_x:.2f}")

if hasattr(sim2, 'bin_biases'):
    print(f"\nBin Biases:")
    print(f"  Mean: {np.mean(sim2.bin_biases):.3f} (target: 1.0)")
    print(f"  Std: {np.std(sim2.bin_biases):.3f}")
    print(f"  Range: [{sim2.bin_biases.min():.3f}, {sim2.bin_biases.max():.3f}]")

print("\n✓ Week 9 sanity plots complete!")
