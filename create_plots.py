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

num_clones = tree['num_clones']

# Create figure with subplots (one per clone)
fig, axes = plt.subplots(num_clones, 1, figsize=(14, 3*num_clones))
if num_clones == 1:
    axes = [axes]

for i in range(num_clones):
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
    
    # Title with clone info
    clone_info = tree['clones'][i]
    parent_str = f"parent: {clone_info['parent_id']}" if clone_info['parent_id'] is not None else "root"
    leaf_str = " [LEAF]" if clone_info['is_leaf'] else ""
    ax.set_title(f"Clone {i} ({parent_str}, {clone_info['num_events']} events){leaf_str}", 
                 fontsize=11, fontweight='bold')
    
    if i == num_clones - 1:
        ax.set_xlabel('Genomic Bin', fontsize=11)

plt.tight_layout()
plt.savefig('sample_outputs/cn_profiles.png', dpi=150, bbox_inches='tight')
print("✓ Saved: sample_outputs/cn_profiles.png")

# Create tree diagram (more sophisticated layout for binary tree)
print("\nCreating tree structure diagram...")

fig, ax = plt.subplots(figsize=(12, 10))

# Build tree structure from parent relationships
def build_layout(tree_data):
    """Create a hierarchical layout for the tree."""
    # Get depth for each node
    depths = {}
    for clone in tree_data['clones']:
        cid = clone['clone_id']
        depth = 0
        current = clone
        while current['parent_id'] is not None:
            depth += 1
            parent_id = current['parent_id']
            current = tree_data['clones'][parent_id]
        depths[cid] = depth
    
    max_depth = max(depths.values())
    
    # Group nodes by depth
    nodes_by_depth = {}
    for cid, depth in depths.items():
        if depth not in nodes_by_depth:
            nodes_by_depth[depth] = []
        nodes_by_depth[depth].append(cid)
    
    # Assign positions
    positions = {}
    y_spacing = 0.8 / (max_depth + 1) if max_depth > 0 else 0.5
    
    for depth in range(max_depth + 1):
        nodes = sorted(nodes_by_depth[depth])
        y = 0.9 - depth * y_spacing
        
        # Distribute nodes evenly at this depth
        n = len(nodes)
        if n == 1:
            x_positions = [0.5]
        else:
            x_positions = np.linspace(0.1, 0.9, n)
        
        for i, cid in enumerate(nodes):
            positions[cid] = (x_positions[i], y)
    
    return positions

positions = build_layout(tree)

# Draw edges
for clone in tree['clones']:
    if clone['parent_id'] is not None:
        parent_id = clone['parent_id']
        child_id = clone['clone_id']
        x1, y1 = positions[parent_id]
        x2, y2 = positions[child_id]
        ax.plot([x1, x2], [y1, y2], 'k-', linewidth=2, alpha=0.5)

# Draw nodes
for clone in tree['clones']:
    clone_id = clone['clone_id']
    x, y = positions[clone_id]
    
    # Color based on root or leaf
    if clone['is_root']:
        color = 'lightgreen'
    elif clone['is_leaf']:
        color = 'lightblue'
    else:
        color = 'lightgray'
    
    # Draw circle
    circle = plt.Circle((x, y), 0.04, color=color, ec='black', linewidth=2, zorder=10)
    ax.add_patch(circle)
    
    # Add label
    label = f"C{clone_id}\n{clone['num_events']} evt"
    ax.text(x, y, label, ha='center', va='center', fontsize=9, fontweight='bold', zorder=11)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1.0)
ax.axis('off')

title = f"Beta-Splitting Clone Tree\n"
title += f"(green=root, blue=leaves, gray=internal)\n"
title += f"{tree['num_leaves']} leaves, {tree['num_clones']} total nodes"
ax.set_title(title, fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('sample_outputs/tree_structure.png', dpi=150, bbox_inches='tight')
print("✓ Saved: sample_outputs/tree_structure.png")

print("\n✓ All plots created successfully!")
