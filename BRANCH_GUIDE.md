# Branch Organization Guide

This repository uses branches to organize different phases of implementation.

## Active Branch
- **`main`** - Current working version with all completed features

## Phase Branches (Historical Documentation)

### Phase A: Parameter Updates (v0.9.0)
- **Branch**: `phase-a-parameter-updates`
- **Purpose**: Paper compliance - parameter defaults and model fixes
- **Key Features**:
  - gain_prob = 0.40
  - Pure Poisson event counts
  - Arm-level haplotype blocks
  - LOH floor + sequencing error
  - Ploidy output
- **Status**: ✅ Merged to main

### Phase B: WGD + Doublets (v1.0.0)
- **Branch**: `phase-b-wgd-doublets`
- **Purpose**: WGD model and doublet allelic model
- **Key Features**:
  - WGD placement (Equations 18-19)
  - Genome doubling (ploidy 2.0 → 4.0)
  - Doublet CN-weighted mixture
  - Updated visualizations
- **Status**: ✅ Merged to main

## How to Use

### To view Phase A work:
```bash
git checkout phase-a-parameter-updates
cat BRANCH_INFO.md
```

### To view Phase B work:
```bash
git checkout phase-b-wgd-doublets
cat BRANCH_INFO.md
```

### To return to current work:
```bash
git checkout main
```

---

**Current Status**: v1.0.0 - Paper Compliant ✅  
**Next**: Week 11 - Core Metrics Implementation
