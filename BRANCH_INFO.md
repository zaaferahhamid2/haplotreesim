# Phase B Branch: WGD + Doublets (v1.0.0)

This branch contains Phase B implementation work for full paper compliance.

## Changes in this branch:
- **WGD (Whole Genome Doubling)** - Section 3.4, Equations 18-19
  - Sample WGD with probability `p_WGD`
  - Place on ONE early edge (depth ≤ `d_WGD = 2`)
  - Select edge proportional to branch length
  - Double ALL bins on BOTH haplotypes BEFORE other events
  - Descendants inherit doubled state (ploidy 2.0 → 4.0)

- **Doublet Allelic Model** - Equation 15
  - Sample doublets with probability `p_doublet`
  - Store two clone indices per doublet cell
  - Compute CN-weighted allele fraction mixture
  - Produces intermediate BAF for complementary states

## Updated Visualizations:
- Week 8: CN profiles with WGD comparison
- Week 9: Read depth with ploidy normalization
- Week 10: BAF patterns with doublets

## Files modified:
- `src/haplotreesim/simulator.py` (WGD placement, doublet sampling, allelic model)
- `src/haplotreesim/data_models.py` (ploidy field)
- `scripts/week8_cn_profiles_with_wgd.py`
- `scripts/week9_depth_with_wgd.py`
- `scripts/week10_baf_with_doublets.py`
- `scripts/comprehensive_test.py`

## Tag: v1.0.0
## Status: ✅ Completed and merged to main
## Date: April 2026
## Achievement: 🎉 Full paper compliance achieved!
