"""
Demonstrate overlapping event handling with a controlled example.
"""

import numpy as np
from haplotreesim.data_models import Clone, CNAEvent, Haplotype
from haplotreesim.event_applier import EventApplier

print("="*70)
print("Demonstrating Overlapping Event Handling")
print("="*70)

# Create a simple parent clone (diploid)
num_bins = 300
parent_cn_A = np.ones(num_bins, dtype=int)
parent_cn_B = np.ones(num_bins, dtype=int)

parent = Clone(
    index=0,
    parent_index=None,
    cn_profile_A=parent_cn_A,
    cn_profile_B=parent_cn_B,
    events=[],
    is_root=True
)

print("\nParent Clone (Root):")
print(f"  All bins: CN_A = 1, CN_B = 1 (diploid)")
print(f"  Total CN = 2 everywhere")

# Create overlapping events on same haplotype
print("\n" + "="*70)
print("EXAMPLE 1: Overlapping Events on SAME Haplotype")
print("="*70)

event1 = CNAEvent(
    start_bin=50,
    end_bin=150,
    haplotype=Haplotype.A,
    amplitude=2
)

event2 = CNAEvent(
    start_bin=100,
    end_bin=200,
    haplotype=Haplotype.A,
    amplitude=3
)

print("\nEvent 1 on Haplotype A:")
print(f"  Bins [50, 150], Amplitude: +2")

print("\nEvent 2 on Haplotype A:")
print(f"  Bins [100, 200], Amplitude: +3")

print("\nOverlap Analysis:")
print(f"  Overlap region: bins [100, 150]")
print(f"  In overlap: BOTH events apply (additive)")

# Apply events
applier = EventApplier(max_copy_number=10)
cn_A, cn_B = applier.apply_events(parent, [event1, event2])

print("\nResulting Copy Numbers:")

# Region 1: Before event1 (bins 0-49)
print(f"\n  Bins [0, 49] (before any event):")
print(f"    CN_A = {cn_A[0]}, CN_B = {cn_B[0]}, Total = {cn_A[0] + cn_B[0]}")
print(f"    Expected: CN_A = 1 (no change)")

# Region 2: Event1 only (bins 50-99)
print(f"\n  Bins [50, 99] (event1 only):")
print(f"    CN_A = {cn_A[50]}, CN_B = {cn_B[50]}, Total = {cn_A[50] + cn_B[50]}")
print(f"    Expected: CN_A = 1 + 2 = 3 ✓" if cn_A[50] == 3 else f"    ERROR: Expected 3, got {cn_A[50]}")

# Region 3: OVERLAP (bins 100-150)
print(f"\n  Bins [100, 150] (OVERLAP: both events):")
print(f"    CN_A = {cn_A[100]}, CN_B = {cn_B[100]}, Total = {cn_A[100] + cn_B[100]}")
print(f"    Expected: CN_A = 1 + 2 + 3 = 6 (ADDITIVE) ✓" if cn_A[100] == 6 else f"    ERROR: Expected 6, got {cn_A[100]}")
print(f"    This demonstrates: overlapping events ADD together!")

# Region 4: Event2 only (bins 151-200)
print(f"\n  Bins [151, 200] (event2 only):")
print(f"    CN_A = {cn_A[151]}, CN_B = {cn_B[151]}, Total = {cn_A[151] + cn_B[151]}")
print(f"    Expected: CN_A = 1 + 3 = 4 ✓" if cn_A[151] == 4 else f"    ERROR: Expected 4, got {cn_A[151]}")

# Region 5: After all events (bins 201-299)
print(f"\n  Bins [201, 299] (after all events):")
print(f"    CN_A = {cn_A[201]}, CN_B = {cn_B[201]}, Total = {cn_A[201] + cn_B[201]}")
print(f"    Expected: CN_A = 1 (no change)")

# Verify correctness
assert cn_A[0] == 1, "Bins before events should be unchanged"
assert cn_A[50] == 3, "Event1 only region should be 1+2=3"
assert cn_A[100] == 6, "Overlap region should be 1+2+3=6 (ADDITIVE)"
assert cn_A[151] == 4, "Event2 only region should be 1+3=4"
assert cn_A[201] == 1, "Bins after events should be unchanged"

print("\n✓ All assertions passed!")

# Example 2: Overlapping events on DIFFERENT haplotypes
print("\n" + "="*70)
print("EXAMPLE 2: Overlapping Events on DIFFERENT Haplotypes")
print("="*70)

event3 = CNAEvent(
    start_bin=50,
    end_bin=150,
    haplotype=Haplotype.A,
    amplitude=2
)

event4 = CNAEvent(
    start_bin=100,
    end_bin=200,
    haplotype=Haplotype.B,
    amplitude=-1
)

print("\nEvent 3 on Haplotype A:")
print(f"  Bins [50, 150], Amplitude: +2")

print("\nEvent 4 on Haplotype B:")
print(f"  Bins [100, 200], Amplitude: -1")

print("\nOverlap Analysis:")
print(f"  Overlap region: bins [100, 150]")
print(f"  In overlap: Events affect DIFFERENT haplotypes (independent)")

cn_A2, cn_B2 = applier.apply_events(parent, [event3, event4])

print("\nResulting Copy Numbers:")

print(f"\n  Bins [50, 99] (event3 only, haplotype A):")
print(f"    CN_A = {cn_A2[50]}, CN_B = {cn_B2[50]}, Total = {cn_A2[50] + cn_B2[50]}")
print(f"    Expected: CN_A = 3, CN_B = 1")

print(f"\n  Bins [100, 150] (OVERLAP: different haplotypes):")
print(f"    CN_A = {cn_A2[100]}, CN_B = {cn_B2[100]}, Total = {cn_A2[100] + cn_B2[100]}")
print(f"    Expected: CN_A = 3 (from event3), CN_B = 0 (from event4)")
print(f"    Total = 3 (not 2+2=4, because events are independent)")

print(f"\n  Bins [151, 200] (event4 only, haplotype B):")
print(f"    CN_A = {cn_A2[151]}, CN_B = {cn_B2[151]}, Total = {cn_A2[151] + cn_B2[151]}")
print(f"    Expected: CN_A = 1, CN_B = 0")

assert cn_A2[100] == 3 and cn_B2[100] == 0, "Different haplotypes should be independent"

print("\n✓ All assertions passed!")

# Example 3: Multiple overlaps on same haplotype
print("\n" + "="*70)
print("EXAMPLE 3: Three Overlapping Events on Same Haplotype")
print("="*70)

event5 = CNAEvent(start_bin=50, end_bin=150, haplotype=Haplotype.A, amplitude=1)
event6 = CNAEvent(start_bin=100, end_bin=200, haplotype=Haplotype.A, amplitude=2)
event7 = CNAEvent(start_bin=120, end_bin=180, haplotype=Haplotype.A, amplitude=1)

print("\nEvent 5: bins [50, 150], Δ = +1")
print("Event 6: bins [100, 200], Δ = +2")
print("Event 7: bins [120, 180], Δ = +1")

print("\nOverlap Analysis:")
print("  Bins [50, 99]:   Event 5 only → Δ = +1")
print("  Bins [100, 119]: Events 5+6 → Δ = +1+2 = +3")
print("  Bins [120, 150]: Events 5+6+7 → Δ = +1+2+1 = +4")
print("  Bins [151, 180]: Events 6+7 → Δ = +2+1 = +3")
print("  Bins [181, 200]: Event 6 only → Δ = +2")

cn_A3, cn_B3 = applier.apply_events(parent, [event5, event6, event7])

print("\nVerification:")
print(f"  Bins [50, 99]:   CN_A = {cn_A3[50]} (expected 2)")
print(f"  Bins [100, 119]: CN_A = {cn_A3[100]} (expected 4)")
print(f"  Bins [120, 150]: CN_A = {cn_A3[120]} (expected 5)")
print(f"  Bins [151, 180]: CN_A = {cn_A3[151]} (expected 4)")
print(f"  Bins [181, 200]: CN_A = {cn_A3[181]} (expected 3)")

assert cn_A3[50] == 2, f"Expected 2, got {cn_A3[50]}"
assert cn_A3[100] == 4, f"Expected 4, got {cn_A3[100]}"
assert cn_A3[120] == 5, f"Expected 5, got {cn_A3[120]}"
assert cn_A3[151] == 4, f"Expected 4, got {cn_A3[151]}"
assert cn_A3[181] == 3, f"Expected 3, got {cn_A3[181]}"

print("\n✓ All assertions passed!")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

print("\nOverlapping Event Handling Verified:")
print("  ✓ Same haplotype, overlapping regions: Effects are ADDITIVE")
print("  ✓ Different haplotypes: Effects are INDEPENDENT")
print("  ✓ Multiple overlaps: All effects sum correctly")
print("  ✓ Implementation uses difference arrays for O(B + M) efficiency")

print("\nAlgorithm (from event_applier.py):")
print("  1. Create difference arrays d_A and d_B")
print("  2. For each event (start, end, haplotype, Δ):")
print("       d[haplotype][start] += Δ")
print("       d[haplotype][end+1] -= Δ")
print("  3. Compute prefix sums: δ[b] = Σ d[j] for j=1 to b")
print("  4. Apply: CN[b] = clip(parent_CN[b] + δ[b], 0, C_max)")

print("\n" + "="*70)
print("✓ VERIFICATION COMPLETE")
print("="*70)
