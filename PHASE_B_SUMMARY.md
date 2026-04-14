# Phase B Summary: WGD + Doublets Implementation

## Completed Features

### 1. WGD (Whole Genome Doubling)
**Implementation**: Section 3.4, Equations 18-19

- ✅ Probability `p_WGD` controls WGD occurrence
- ✅ WGD placed on ONE early edge (depth ≤ `d_WGD = 2`)
- ✅ Edge selected proportional to branch length
- ✅ WGD doubles ALL bins on BOTH haplotypes
- ✅ Applied BEFORE focal/arm/chr events
- ✅ Descendants inherit doubled state

**Verification**: Clones with WGD have ploidy ~4.0 (vs 2.0 baseline)

### 2. Doublet Allelic Model
**Implementation**: Equation 15

- ✅ Doublets sampled with probability `p_doublet`
- ✅ Two clones independently drawn per doublet
- ✅ CN-weighted allele fraction mixture
- ✅ Correct total depth (sum of both clones)
- ✅ Produces intermediate BAF for complementary states

**Verification**: Doublets show BAF ≈ 0.5 for complementary CN states

## Updated Graphs

### Week 8: CN Profiles with WGD
**File**: `sample_outputs/week8_cn_profiles_with_wgd.png`
- Panel A: Diploid baseline (no WGD)
- Panel B: With WGD showing ploidy ~4.0

### Week 9: Read Depth with WGD
**File**: `sample_outputs/week9_depth_vs_cn_wgd.png`
- Panel A: Diploid clones (depth ∝ CN/2)
- Panel B: Tetraploid clones (depth ∝ CN/4)
- Demonstrates need for ploidy normalization

### Week 10: BAF with Doublets
**File**: `sample_outputs/week10_baf_with_doublets.png`
- Panel A: Singlet BAF patterns
- Panel B: Doublet BAF patterns (intermediate values)
- Panel C: mBAF distribution comparison

## Parameters Added
```python
# WGD parameters
prob_wgd: float = 0.0      # Probability of WGD
d_wgd: int = 2             # Max depth for WGD placement

# Already existed
prob_doublet: float = 0.0  # Probability of doublet
```

## Test Results
```
WGD Test (prob_wgd=1.0):
  Clone 0: ploidy = 2.00 (root, no WGD)
  Clone 2: ploidy = 4.03 [WGD node]
  Clone 4: ploidy = 4.03 [descendant]
  Clone 5: ploidy = 4.44 [descendant]
  ✓ WGD working correctly

Doublet Test (prob_doublet=0.15):
  Singlets: 249/300
  Doublets: 51/300 (17.0%)
  ✓ Doublets sampled correctly
  ✓ CN-weighted mixture applied
```

## Version Tags

- **v0.9.0**: Phase A (parameter updates, LOH floor, haplotype blocks, ploidy)
- **v1.0.0**: Phase B (WGD + doublets) - **PAPER COMPLIANT**

## Next Steps

1. ✅ Phase A: Parameter updates
2. ✅ Phase B: WGD + Doublets  
3. ⏭️ Week 11: Core metrics (HSCN error, LOH F1, breakpoints, etc.)
4. ⏭️ Week 12-14: Full evaluation framework
5. ⏭️ Week 15+: Baseline integration

---

**Status**: Simulator is now fully compliant with paper specification ✓
