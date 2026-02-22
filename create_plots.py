"""
Create CN profile plots for sample outputs
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import json

print("Creating CN profile plots...")

# Load tree structure to get event info
with open('sample_outputs/tree_structure.json', 'r') as f:
    tree = json.load(f)

# Create figure with subplots
fig, axes = plt.subplots(5, 1, figsize=(14, 12))

for i in range(5):
    # Load CN profile
    cn = np.loadtxt(f'sample_outputs/clone_{i}_cn_profile.csv', 
                    delimiter=',', skiprows=1)
    
    ax = axes[i]
    
    # Plot
    ax.plot(cn[:, 0], label='Haplotype A', linewidth=1.5, alpha=0.8, color='#1f77b4')
    ax.plot(cn[:, 1], label='Haplotype B', linewidth=1.5, alpha=0.8, color='#ff7f0e')
    ax.plot(cn[:, 2], label='Total CN', linewidth=2, color='black', alpha=0.7)
    
    # Reference line at diploid
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.3, linewidth=1)
    
    # Labels and formatting
    ax.set_ylabel('Copy Number', fontsize=11)
    ax.set_ylim(-0.5, 5)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # Title with event info
    clone_info = tree['clones'][i]
    parent_str = f"parent: {clone_info['parent_id']}" if clone_info['parent_id'] is not None else "root"
    ax.set_title(f"Clone {i} ({parent_str}, {clone_info['num_events']} events)", 
                 fontsize=12, fontweight='bold')
    
    if i == 4:
        ax.set_xlabel('Genomic Bin', fontsize=11)

plt.tight_layout()
plt.savefig('sample_outputs/cn_profiles.png', dpi=150, bbox_inches='tight')
print("✓ Saved: sample_outputs/cn_profiles.png")

# Create tree diagram
print("\nCreating tree structure diagram...")

fig, ax = plt.subplots(figsize=(10, 8))

# Simple tree layout
positions = {
    0: (0.5, 0.9),   # Root at top
    1: (0.2, 0.6),   # Children spread out
    2: (0.4, 0.6),
    3: (0.6, 0.6),
    4: (0.8, 0.6)
}

# Draw edges
for clone_info in tree['clones']:
    if clone_info['parent_id'] is not None:
        parent_id = clone_info['parent_id']
        child_id = clone_info['clone_id']
        x1, y1 = positions[parent_id]
        x2, y2 = positions[child_id]
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=2, alpha=0.5)

# Draw nodes
for clone_info in tree['clones']:
    clone_id = clone_info['clone_id']
    x, y = positions[clone_id]
    
    # Color
    color = 'lightgreen' if clone_info['is_root'] else 'lightblue'
    
    # Circle
    circle = plt.Circle((x, y), 0.05, color=color, ec='black', linewidth=2, zorder=10)
    ax.add_patch(circle)
    
    # Label
    label = f"C{clone_id}\n{clone_info['num_events']} evt"
    ax.text(x, y, label, ha='center', va='center', fontsize=10, fontweight='bold', zorder=11)

ax.set_xlim(0, 1)
ax.set_ylim(0.3, 1.0)
ax.axis('off')
ax.set_title('Clone Tree Structure\n(green=root, blue=derived clones)', 
             fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('sample_outputs/tree_structure.png', dpi=150, bbox_inches='tight')
print("✓ Saved: sample_outputs/tree_structure.png")

print("\n✓ All plots created successfully!")
