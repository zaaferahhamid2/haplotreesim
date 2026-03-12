# Overlapping CNA Event Verification

This document provides manual verification that overlapping copy number events are handled correctly in HaploTreeSim.

## Implementation Overview

**Algorithm** (from `event_applier.py`):
1. Initialize child CN by copying parent: `c_child ← c_parent`
2. If WGD event: `c_child ← 2 × c_child`
3. Create difference arrays `d_A` and `d_B` (length B+1, initialized to 0)
4. For each gain/loss event `(start, end, haplotype, Δ)`:
   - `d[haplotype][start] += Δ`
   - `d[haplotype][end+1] -= Δ`
5. Compute prefix sums: `δ[b] = Σ(d[j])` for j=1 to b
6. Apply with clipping: `c[b] = clip(c[b] + δ[b], 0, C_max)`

**Key Properties**:
- ✅ **Additive on same haplotype**: Overlapping events sum their effects
- ✅ **Independent on different haplotypes**: Events affect each haplotype separately
- ✅ **Efficient**: O(B + M) time for B bins and M events per edge
- ✅ **Exact**: Matches the intended semantics from Section 3.5 of the paper

---

## Example 1: Overlapping Events on Same Haplotype

**Setup**:
- Parent: CN_A = 1, CN_B = 1 (diploid) across all 300 bins
- Event 1: Bins [50, 150], Haplotype A, Δ = +2
- Event 2: Bins [100, 200], Haplotype A, Δ = +3

**Overlap Region**: Bins [100, 150] — both events apply

**Expected Results**:
```
Bins [0, 49]:    CN_A = 1           (no events)
Bins [50, 99]:   CN_A = 1 + 2 = 3   (event 1 only)
Bins [100, 150]: CN_A = 1 + 2 + 3 = 6  (BOTH events, ADDITIVE)
Bins [151, 200]: CN_A = 1 + 3 = 4   (event 2 only)
Bins [201, 299]: CN_A = 1           (no events)
```

**Actual Results**:
```
Bins [0, 49]:    CN_A = 1 ✓
Bins [50, 99]:   CN_A = 3 ✓
Bins [100, 150]: CN_A = 6 ✓  (ADDITIVE VERIFIED)
Bins [151, 200]: CN_A = 4 ✓
Bins [201, 299]: CN_A = 1 ✓
```

**Verification**: ✅ Overlapping events on same haplotype sum correctly

---

## Example 2: Overlapping Events on Different Haplotypes

**Setup**:
- Parent: CN_A = 1, CN_B = 1 (diploid)
- Event 3: Bins [50, 150], Haplotype A, Δ = +2
- Event 4: Bins [100, 200], Haplotype B, Δ = -1

**Overlap Region**: Bins [100, 150] — events on different haplotypes

**Expected Results**:
```
Bins [50, 99]:   CN_A = 3, CN_B = 1    (event 3 only)
Bins [100, 150]: CN_A = 3, CN_B = 0    (independent effects)
Bins [151, 200]: CN_A = 1, CN_B = 0    (event 4 only)
```

**Actual Results**:
```
Bins [50, 99]:   CN_A = 3, CN_B = 1 ✓
Bins [100, 150]: CN_A = 3, CN_B = 0 ✓  (INDEPENDENT VERIFIED)
Bins [151, 200]: CN_A = 1, CN_B = 0 ✓
```

**Verification**: ✅ Events on different haplotypes are independent

---

## Example 3: Three Overlapping Events on Same Haplotype

**Setup**:
- Parent: CN_A = 1, CN_B = 1 (diploid)
- Event 5: Bins [50, 150], Haplotype A, Δ = +1
- Event 6: Bins [100, 200], Haplotype A, Δ = +2
- Event 7: Bins [120, 180], Haplotype A, Δ = +1

**Multiple Overlap Regions**:
```
Bins [50, 99]:   Event 5 only          → Δ = +1
Bins [100, 119]: Events 5 + 6          → Δ = +1 + 2 = +3
Bins [120, 150]: Events 5 + 6 + 7      → Δ = +1 + 2 + 1 = +4
Bins [151, 180]: Events 6 + 7          → Δ = +2 + 1 = +3
Bins [181, 200]: Event 6 only          → Δ = +2
```

**Expected & Actual Results**:
```
Bins [50, 99]:   CN_A = 1 + 1 = 2       ✓
Bins [100, 119]: CN_A = 1 + 3 = 4       ✓
Bins [120, 150]: CN_A = 1 + 4 = 5       ✓  (THREE events sum)
Bins [151, 180]: CN_A = 1 + 3 = 4       ✓
Bins [181, 200]: CN_A = 1 + 2 = 3       ✓
```

**Verification**: ✅ Multiple overlapping events sum correctly

---

## Real Simulation Example

Running a small simulation with 3 leaf clones (5 total nodes):

### Tree Structure
```
Clone 0 (Root) → Clone 1 [LEAF]
               → Clone 2 → Clone 3 [LEAF]
                         → Clone 4 [LEAF]
```

### Event Details by Edge

**Edge: Clone 0 → Clone 1**
- Event 1: Bins [0, 299], Haplotype B, Δ = +3
- Result: All bins have CN_A = 1, CN_B = 4, Total = 5

**Edge: Clone 0 → Clone 2**
- Event 1: Bins [59, 101], Haplotype A, Δ = +5
- Result: Bins [59-101] have CN_A = 6, rest have CN_A = 1

**Edge: Clone 2 → Clone 3**
- Event 1: Bins [120, 165], Haplotype B, Δ = -2
- Event 2: Bins [18, 63], Haplotype B, Δ = +2
- No overlap (different regions, same haplotype)
- Result: Combined effects create diverse CN profile

**Edge: Clone 2 → Clone 4**
- Event 1: Bins [8, 54], Haplotype A, Δ = -3
- Inherits parent's [59-101] gain, adds [8-54] loss
- Result: CN_A varies from 0 to 6 across genome

### Manual Trace: Root → Clone 3

**Step 1: Root (Clone 0)**
- All bins: CN_A = 1, CN_B = 1

**Step 2: Clone 0 → Clone 2**
- Inherit: CN_A = 1, CN_B = 1
- Apply Event: Bins [59, 101] get Δ = +5 on haplotype A
- Result: Bins [59-101]: CN_A = 6, CN_B = 1

**Step 3: Clone 2 → Clone 3**
- Inherit from Clone 2: CN_A = 6 (bins 59-101), else 1; CN_B = 1
- Apply Event 1: Bins [120, 165] get Δ = -2 on haplotype B → CN_B = 0 (clipped)
- Apply Event 2: Bins [18, 63] get Δ = +2 on haplotype B → CN_B = 3
- Final state verified: ✓

---

## Difference Array Mechanics

**Example**: Two events on Haplotype A:
- Event 1: bins [50, 150], Δ = +2
- Event 2: bins [100, 200], Δ = +3

**Difference Array Construction**:
```
d[50]  += 2    (event 1 starts)
d[151] -= 2    (event 1 ends)
d[100] += 3    (event 2 starts)
d[201] -= 3    (event 2 ends)
```

**Prefix Sum (δ)**:
```
δ[0-49]   = 0
δ[50-99]  = 2      (event 1 effect)
δ[100-150] = 2+3=5  (both events, ADDITIVE)
δ[151-200] = 3      (event 2 effect)
δ[201+]    = 0
```

**Final CN** = `clip(parent_CN + δ, 0, C_max)`

---

## Verification Scripts

Two verification scripts are provided:

1. **`demonstrate_overlaps.py`** - Controlled examples with assertions
2. **`verify_overlapping_events.py`** - Real simulation trace

Both produce output showing:
- Event details (start, end, haplotype, amplitude)
- Overlap detection and analysis
- Before/after CN values
- Manual trace from root to each leaf

---

## Conclusion

✅ **Overlapping events are handled correctly**:
- Same haplotype: Effects are **additive**
- Different haplotypes: Effects are **independent**
- Implementation matches Section 3.5 of the paper
- Uses efficient difference arrays (O(B + M) time)
- All verification tests pass

The implementation correctly handles:
- Single events
- Non-overlapping events
- Overlapping events (2+)
- Multiple overlaps (3+)
- Events across tree paths
- WGD + gain/loss combinations
- Clipping to [0, C_max]
