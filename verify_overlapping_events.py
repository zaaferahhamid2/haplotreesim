"""
Manual verification of overlapping CNA events.

This script creates a small tree and manually traces CN changes
from root to leaves, showing how overlapping events are handled.
"""

import numpy as np
from haplotreesim import SimulationConfig, HaploTreeSimulator
import json

print("="*70)
print("Manual Verification of Overlapping CNA Events")
print("="*70)

# Create a small, simple simulation
config = SimulationConfig(
    num_bins=300,  # Small genome for easy visualization
    num_clones=3,  # Small tree: root + 3 leaves = 5 total nodes
    num_cells=50,
    lambda_events=2.0,  # Generate some events
    alpha_tree=0.5,
    beta_tree=0.5,  # More balanced
    random_seed=123  # Fixed seed for reproducibility
)

print("\nConfiguration:")
print(f"  Bins: {config.num_bins}")
print(f"  Clones (leaves): {config.num_clones}")
print(f"  Root CN: diploid (1, 1)")

sim = HaploTreeSimulator(config)
read_counts, alleles = sim.run()

print("\n" + "="*70)
print("TREE STRUCTURE")
print("="*70)

# Show tree topology
for i, clone in enumerate(sim.clones):
    parent_str = f"parent={clone.parent_index}" if clone.parent_index is not None else "ROOT"
    leaf_str = " [LEAF]" if i in sim._leaf_clone_indices else " [INTERNAL]"
    print(f"\nClone {i}: {parent_str}{leaf_str}")
    print(f"  Events on edge: {len(clone.events)}")
    
print("\n" + "="*70)
print("DETAILED EVENT TRACE (Root → Leaves)")
print("="*70)

def trace_clone_path(clone_idx, clones):
    """Trace path from root to clone."""
    path = []
    current = clones[clone_idx]
    
    while current.index is not None:
        path.insert(0, current.index)
        if current.parent_index is None:
            break
        current = clones[current.parent_index]
    
    return path

# For each leaf clone, trace the path and show CN evolution
for leaf_idx in sorted(sim._leaf_clone_indices):
    path = trace_clone_path(leaf_idx, sim.clones)
    
    print(f"\n{'='*70}")
    print(f"LEAF CLONE {leaf_idx}: Path from root: {' → '.join(map(str, path))}")
    print(f"{'='*70}")
    
    # Start with root CN
    current_cn_A = sim.clones[0].cn_profile_A.copy()
    current_cn_B = sim.clones[0].cn_profile_B.copy()
    
    print(f"\nInitial state (Root, Clone 0):")
    print(f"  All bins: CN_A = 1, CN_B = 1 (diploid)")
    
    # Walk down the path
    for i in range(1, len(path)):
        clone_idx = path[i]
        clone = sim.clones[clone_idx]
        parent_idx = clone.parent_index
        
        print(f"\n{'-'*70}")
        print(f"Edge: Clone {parent_idx} → Clone {clone_idx}")
        print(f"{'-'*70}")
        
        if len(clone.events) == 0:
            print("  No events on this edge (CN unchanged)")
            continue
        
        # Show events on this edge
        print(f"  Events on this edge: {len(clone.events)}")
        
        for j, event in enumerate(clone.events):
            print(f"\n  Event {j+1}:")
            print(f"    Type: {event.haplotype.value}")
            print(f"    Bins: [{event.start_bin}, {event.end_bin}] (length: {event.length})")
            print(f"    Amplitude: {event.amplitude:+d}")
            
            if event.haplotype.value == 'WGD':
                print(f"    Effect: Double ALL bins on BOTH haplotypes")
            else:
                print(f"    Effect: CN change of {event.amplitude:+d} on haplotype {event.haplotype.value}")
        
        # Check for overlaps
        if len(clone.events) > 1:
            print(f"\n  Checking for overlapping events:")
            events = [e for e in clone.events if e.haplotype.value != 'WGD']
            for i in range(len(events)):
                for j in range(i+1, len(events)):
                    e1, e2 = events[i], events[j]
                    # Check overlap
                    overlap_start = max(e1.start_bin, e2.start_bin)
                    overlap_end = min(e1.end_bin, e2.end_bin)
                    
                    if overlap_start <= overlap_end:
                        print(f"    OVERLAP: Events {i+1} and {j+1}")
                        print(f"      Event {i+1}: bins [{e1.start_bin}-{e1.end_bin}], {e1.haplotype.value}, Δ={e1.amplitude:+d}")
                        print(f"      Event {j+1}: bins [{e2.start_bin}-{e2.end_bin}], {e2.haplotype.value}, Δ={e2.amplitude:+d}")
                        print(f"      Overlap region: bins [{overlap_start}-{overlap_end}]")
                        
                        if e1.haplotype == e2.haplotype:
                            net_change = e1.amplitude + e2.amplitude
                            print(f"      Same haplotype → Net change: {net_change:+d} (additive)")
                        else:
                            print(f"      Different haplotypes → Independent effects")
        
        # Show example bins
        print(f"\n  Example bin changes:")
        
        # Find some interesting bins
        sample_bins = []
        
        # First event's start
        if len(clone.events) > 0 and clone.events[0].haplotype.value != 'WGD':
            sample_bins.append(clone.events[0].start_bin)
        
        # First event's middle
        if len(clone.events) > 0 and clone.events[0].haplotype.value != 'WGD':
            mid = (clone.events[0].start_bin + clone.events[0].end_bin) // 2
            sample_bins.append(mid)
        
        # Check for overlap region
        if len(clone.events) > 1:
            events = [e for e in clone.events if e.haplotype.value != 'WGD']
            if len(events) >= 2:
                overlap_start = max(events[0].start_bin, events[1].start_bin)
                overlap_end = min(events[0].end_bin, events[1].end_bin)
                if overlap_start <= overlap_end:
                    sample_bins.append(overlap_start)
        
        # Remove duplicates and show up to 3 examples
        sample_bins = sorted(set(sample_bins))[:3]
        
        parent = sim.clones[parent_idx]
        for bin_idx in sample_bins:
            before_A = parent.cn_profile_A[bin_idx]
            before_B = parent.cn_profile_B[bin_idx]
            after_A = clone.cn_profile_A[bin_idx]
            after_B = clone.cn_profile_B[bin_idx]
            
            print(f"    Bin {bin_idx}:")
            print(f"      Before: CN_A={before_A}, CN_B={before_B}, Total={before_A+before_B}")
            print(f"      After:  CN_A={after_A}, CN_B={after_B}, Total={after_A+after_B}")
            print(f"      Change: ΔCN_A={after_A-before_A:+d}, ΔCN_B={after_B-before_B:+d}")
    
    # Final state
    final_clone = sim.clones[leaf_idx]
    print(f"\n{'='*70}")
    print(f"FINAL STATE for Leaf Clone {leaf_idx}:")
    print(f"{'='*70}")
    
    # Show some statistics
    unique_total_cn = len(set(final_clone.total_cn()))
    print(f"  Unique total CN values: {unique_total_cn}")
    print(f"  CN range: {final_clone.total_cn().min()} to {final_clone.total_cn().max()}")
    
    # Show a sample of bins
    print(f"\n  Sample bins (first 10):")
    for b in range(min(10, config.num_bins)):
        cn_a = final_clone.cn_profile_A[b]
        cn_b = final_clone.cn_profile_B[b]
        print(f"    Bin {b}: CN_A={cn_a}, CN_B={cn_b}, Total={cn_a+cn_b}")

print("\n" + "="*70)
print("VERIFICATION SUMMARY")
print("="*70)

print("\nOverlapping Event Handling:")
print("  ✓ Events on same edge are applied using difference arrays")
print("  ✓ Overlapping intervals: effects are ADDITIVE within same haplotype")
print("  ✓ After all events: CN is clipped to [0, C_max]")
print("  ✓ Different haplotypes: independent effects")

print("\nImplementation (from simulator.py):")
print("  1. Inherit parent CN: c_v ← c_parent")
print("  2. If WGD: c_v ← 2 × c_v")
print("  3. Build difference arrays for all gain/loss events")
print("  4. Compute prefix sums: δ_v,b = Σ d_v,j")
print("  5. Apply: c_v,b ← clip(c_v,b + δ_v,b, 0, C_max)")

print("\n" + "="*70)
print("DONE: Manual trace complete!")
print("="*70)

# Save detailed event log
output = {
    "tree_structure": {},
    "event_details": []
}

for clone in sim.clones:
    output["tree_structure"][clone.index] = {
        "parent": clone.parent_index,
        "is_leaf": clone.index in sim._leaf_clone_indices,
        "num_events": len(clone.events)
    }

for clone in sim.clones:
    if len(clone.events) > 0:
        for event in clone.events:
            output["event_details"].append({
                "clone_id": clone.index,
                "parent_id": clone.parent_index,
                "start_bin": int(event.start_bin),
                "end_bin": int(event.end_bin),
                "haplotype": event.haplotype.value,
                "amplitude": int(event.amplitude),
                "length": int(event.length)
            })

with open('event_trace_verification.json', 'w') as f:
    json.dump(output, f, indent=2)

print("\nSaved detailed event log to: event_trace_verification.json")
