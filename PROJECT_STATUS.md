# HaploTreeSim Project Status

**Current Version**: v1.0.0 (Paper-Compliant Simulator)  
**Last Updated**: April 14, 2026

---

## ✅ Completed Phases

### Phase 0: Setup (Weeks 1-4)
- [x] Biological fundamentals learned
- [x] scDNA-seq workflow understood
- [x] Paper reading completed

### Phase 1: Matrix-Mode Simulator (Weeks 5-10)
- [x] Week 5: Repository + configuration
- [x] Week 6: Beta-splitting clone tree
- [x] Week 7: CNA event generator + WGD placement
- [x] Week 8: Global segmentation + breakpoints
- [x] Week 9: Read-depth observation model
- [x] Week 10: Allelic observation model

### Phase A: Parameter Updates (v0.9.0)
- [x] Updated `gain_prob` to 0.40 (60% losses)
- [x] Removed +1 from event count (pure Poisson)
- [x] Fixed focal `L_max` to 10 Mb
- [x] Added `epsilon_seq`, `epsilon_floor` parameters
- [x] Fixed haplotype blocks to arm-level
- [x] Implemented LOH floor + sequencing error
- [x] Added ploidy computation and output

### Phase B: WGD + Doublets (v1.0.0)
- [x] WGD placement model (Equations 18-19)
- [x] WGD propagation to descendants
- [x] Doublet sampling
- [x] Doublet allelic model (Equation 15)
- [x] Updated visualization graphs

---

## 📊 Current Capabilities

### Simulator Features
1. **Clone Tree**: Beta-splitting with K leaves, controlled imbalance
2. **CNA Events**: Focal/arm/chromosome with configurable rates
3. **WGD**: Early placement, genome doubling, descendant propagation
4. **Haplotype Blocks**: Chromosome-arm level with phase switching
5. **Read Depth**: Negative Binomial with library size variation
6. **Allelic Counts**: Beta-Binomial with LOH floor, sequencing error
7. **Doublets**: CN-weighted allelic mixture
8. **Ploidy**: Mean genome-wide output for normalization

### Visualization
- Week 8: CN profiles (diploid vs WGD)
- Week 9: Read depth vs CN (ploidy normalization)
- Week 10: BAF patterns (singlets vs doublets)

---

## 🎯 Next Steps (Week 11+)

### Phase 2: Metrics + Evaluation (Weeks 11-14)

#### Week 11: Core Metrics ⏭️ NEXT
- [ ] HSCN swap-invariant error (Equation 57)
- [ ] Total CN MSE (Equation 58)
- [ ] LOH detection (precision/recall/F1)
- [ ] Mirrored-subclone resolution (Equation 60)

#### Week 12: Cluster/Clone Mapping
- [ ] Hungarian matching for K̂ ≠ K
- [ ] Merge/split reconciliation
- [ ] ARI and NMI metrics

#### Week 13: Tree Metrics
- [ ] Normalized Robinson-Foulds distance
- [ ] Event matching (interval overlap + haplotype)
- [ ] Ancestor-descendant accuracy

#### Week 14: Benchmark Pipeline
- [ ] Input/output standardization
- [ ] Run benchmark.py script
- [ ] Automated evaluation

---

## 📈 Implementation Progress
```
Weeks 1-4:   ████████████████████ 100% (Setup)
Weeks 5-10:  ████████████████████ 100% (Simulator Core)
Phase A+B:   ████████████████████ 100% (Paper Compliance)
Week 11:     ░░░░░░░░░░░░░░░░░░░░   0% (Metrics) ← YOU ARE HERE
Weeks 12-28: ░░░░░░░░░░░░░░░░░░░░   0% (Remaining)
```

---

## 🔧 Key Files

### Core Implementation
- `src/haplotreesim/simulator.py` - Main simulator
- `src/haplotreesim/data_models.py` - Configuration + data classes
- `src/haplotreesim/beta_tree_builder.py` - Beta-splitting tree
- `src/haplotreesim/event_generator.py` - CNA event sampling
- `src/haplotreesim/event_applier.py` - Event application
- `src/haplotreesim/segment_detector.py` - Breakpoint detection

### Documentation
- `README.md` - Project overview
- `PHASE_B_SUMMARY.md` - Phase B implementation details
- `PROJECT_STATUS.md` - This file

### Visualizations
- `scripts/week8_cn_profiles_with_wgd.py`
- `scripts/week9_depth_with_wgd.py`
- `scripts/week10_baf_with_doublets.py`
- `scripts/comprehensive_test.py`

---

## 📦 Version History

- **v0.1.0**: Initial simulator (basic functionality)
- **v0.9.0**: Phase A updates (parameter compliance)
- **v1.0.0**: Phase B updates (WGD + doublets) ✓ **CURRENT**

---

**Status**: ✅ Simulator complete and paper-compliant  
**Ready for**: Week 11 (Core Metrics Implementation)
