"""
Week 10 Deliverable: BAF/mBAF distribution plots by CN state across coverage regimes.

BAF = B-allele frequency = alternate / (alternate + reference)
mBAF = mirrored BAF = |BAF - 0.5| + 0.5 (mirrors around 0.5)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from haplotreesim import SimulationConfig, HaploTreeSimulator

print("Creating Week 10: BAF/mBAF Distribution Plots")
print("="*60)

# Test across multiple coverage regimes
coverage_regimes = [
    {'name': 'Low (0.01x)', 'mean_allelic': 5.0, 'color': '#e74c3c'},
    {'name': 'Medium (0.05x)', 'mean_allelic': 25.0, 'color': '#3498db'},
    {'name': 'High (0.1x)', 'mean_allelic': 50.0, 'color': '#2ecc71'},
]

# Store results for each coverage
results = []

for regime in coverage_regimes:
    print(f"\nSimulating {regime['name']} coverage...")
    
    config = SimulationConfig(
        chromosome='chr1',
        bin_width=500000,
        num_clones=5,
        num_cells=200,  # More cells for better statistics
        lambda_events=2.0,  # More events for variety
        mean_allelic_coverage=regime['mean_allelic'],
        prob_phase_switch=0.05,
        random_seed=42
    )
    
    sim = HaploTreeSimulator(config)
    read_counts, (alt, ref, total) = sim.run()
    
    # Compute BAF and mBAF
    # BAF = a / t (where t > 0)
    mask = total > 0
    baf = np.zeros_like(alt, dtype=float)
    baf[mask] = alt[mask] / total[mask]
    
    # mBAF = |BAF - 0.5| + 0.5
    mbaf = np.abs(baf - 0.5) + 0.5
    
    # Collect by CN state
    cn_to_baf = {}
    cn_to_mbaf = {}
    
    for n in range(len(sim.cells)):
        cell = sim.cells[n]
        clone = sim.clones[cell.clone_assignment]
        
        for s, segment in enumerate(sim.segments):
            if total[n, s] == 0:
                continue
            
            cn_A, cn_B, tcn = clone.get_segment_cn(segment)
            
            # Store by total CN
            if tcn not in cn_to_baf:
                cn_to_baf[tcn] = []
                cn_to_mbaf[tcn] = []
            
            cn_to_baf[tcn].append(baf[n, s])
            cn_to_mbaf[tcn].append(mbaf[n, s])
    
    results.append({
        'regime': regime,
        'cn_to_baf': cn_to_baf,
        'cn_to_mbaf': cn_to_mbaf,
        'config': config,
        'sim': sim
    })

print("\n" + "="*60)
print("Creating Plots...")
print("="*60)

# Create figure with subplots
fig = plt.figure(figsize=(18, 12))

# Plot 1: BAF distributions by CN for each coverage regime
for idx, result in enumerate(results):
    ax = plt.subplot(2, 3, idx + 1)
    
    cn_to_baf = result['cn_to_baf']
    regime = result['regime']
    
    # Plot violin plots for each CN state
    cns = sorted(cn_to_baf.keys())
    data = [cn_to_baf[cn] for cn in cns]
    
    positions = np.arange(len(cns))
    parts = ax.violinplot(data, positions=positions, widths=0.7, 
                          showmeans=True, showmedians=True)
    
    # Color by CN
    for pc in parts['bodies']:
        pc.set_facecolor(regime['color'])
        pc.set_alpha(0.7)
    
    ax.set_xticks(positions)
    ax.set_xticklabels([f'CN={cn}' for cn in cns])
    ax.set_ylabel('BAF (B-Allele Frequency)', fontsize=11)
    ax.set_xlabel('Copy Number State', fontsize=11)
    ax.set_title(f'{regime["name"]} Coverage\nBAF by CN State', 
                 fontsize=12, fontweight='bold')
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.grid(True, alpha=0.3, axis='y')

# Plot 2: mBAF distributions by CN for each coverage regime
for idx, result in enumerate(results):
    ax = plt.subplot(2, 3, idx + 4)
    
    cn_to_mbaf = result['cn_to_mbaf']
    regime = result['regime']
    
    cns = sorted(cn_to_mbaf.keys())
    data = [cn_to_mbaf[cn] for cn in cns]
    
    positions = np.arange(len(cns))
    parts = ax.violinplot(data, positions=positions, widths=0.7,
                          showmeans=True, showmedians=True)
    
    for pc in parts['bodies']:
        pc.set_facecolor(regime['color'])
        pc.set_alpha(0.7)
    
    ax.set_xticks(positions)
    ax.set_xticklabels([f'CN={cn}' for cn in cns])
    ax.set_ylabel('mBAF (Mirrored BAF)', fontsize=11)
    ax.set_xlabel('Copy Number State', fontsize=11)
    ax.set_title(f'{regime["name"]} Coverage\nmBAF by CN State', 
                 fontsize=12, fontweight='bold')
    ax.set_ylim(0.45, 1.05)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    ax.grid(True, alpha=0.3, axis='y')

plt.suptitle('Week 10: BAF/mBAF Distributions by Copy Number Across Coverage Regimes', 
             fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('sample_outputs/week10_baf_distributions.png', dpi=150, bbox_inches='tight')
print("\n✓ Saved: sample_outputs/week10_baf_distributions.png")

# Create second figure: CN-specific BAF scatter plots
fig2 = plt.figure(figsize=(18, 10))

cn_states_to_plot = [1, 2, 3, 4]  # Focus on these CN states

for plot_idx, cn_state in enumerate(cn_states_to_plot):
    ax = plt.subplot(2, 4, plot_idx + 1)
    
    # Scatter BAF vs total depth for this CN state
    for result in results:
        regime = result['regime']
        sim = result['sim']
        alt, ref, total = sim.cells[0].allele_counts  # Get structure
        
        # Collect all data points for this CN
        bafs = []
        depths = []
        
        for n in range(len(sim.cells)):
            cell = sim.cells[n]
            clone = sim.clones[cell.clone_assignment]
            alt_n, ref_n, total_n = cell.allele_counts
            
            for s, segment in enumerate(sim.segments):
                if total_n[s] == 0:
                    continue
                
                cn_A, cn_B, tcn = clone.get_segment_cn(segment)
                
                if tcn == cn_state:
                    bafs.append(alt_n[s] / total_n[s])
                    depths.append(total_n[s])
        
        if bafs:
            ax.scatter(depths, bafs, alpha=0.4, s=10, 
                      color=regime['color'], label=regime['name'])
    
    ax.set_xlabel('Total Allelic Depth', fontsize=10)
    ax.set_ylabel('BAF', fontsize=10)
    ax.set_title(f'CN={cn_state}', fontsize=11, fontweight='bold')
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
    ax.grid(True, alpha=0.3)
    if plot_idx == 0:
        ax.legend(fontsize=8)

# mBAF scatter plots
for plot_idx, cn_state in enumerate(cn_states_to_plot):
    ax = plt.subplot(2, 4, plot_idx + 5)
    
    for result in results:
        regime = result['regime']
        sim = result['sim']
        
        mbafs = []
        depths = []
        
        for n in range(len(sim.cells)):
            cell = sim.cells[n]
            clone = sim.clones[cell.clone_assignment]
            alt_n, ref_n, total_n = cell.allele_counts
            
            for s, segment in enumerate(sim.segments):
                if total_n[s] == 0:
                    continue
                
                cn_A, cn_B, tcn = clone.get_segment_cn(segment)
                
                if tcn == cn_state:
                    baf_val = alt_n[s] / total_n[s]
                    mbaf_val = abs(baf_val - 0.5) + 0.5
                    mbafs.append(mbaf_val)
                    depths.append(total_n[s])
        
        if mbafs:
            ax.scatter(depths, mbafs, alpha=0.4, s=10,
                      color=regime['color'], label=regime['name'])
    
    ax.set_xlabel('Total Allelic Depth', fontsize=10)
    ax.set_ylabel('mBAF', fontsize=10)
    ax.set_title(f'CN={cn_state}', fontsize=11, fontweight='bold')
    ax.set_ylim(0.45, 1.05)
    ax.axhline(0.5, color='gray', linestyle='--', alpha=0.5)
    ax.grid(True, alpha=0.3)

plt.suptitle('Week 10: BAF/mBAF vs Depth by Copy Number State', 
             fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.savefig('sample_outputs/week10_baf_vs_depth.png', dpi=150, bbox_inches='tight')
print("✓ Saved: sample_outputs/week10_baf_vs_depth.png")

print("\n" + "="*60)
print("Summary Statistics")
print("="*60)

for result in results:
    regime = result['regime']
    cn_to_baf = result['cn_to_baf']
    
    print(f"\n{regime['name']}:")
    print(f"  Mean allelic coverage: {regime['mean_allelic']}")
    
    for cn in sorted(cn_to_baf.keys())[:5]:  # Show first 5 CN states
        bafs = cn_to_baf[cn]
        if bafs:
            print(f"  CN={cn}: {len(bafs)} observations, "
                  f"mean BAF={np.mean(bafs):.3f}, std={np.std(bafs):.3f}")

print("\n✓ Week 10 plots complete!")
