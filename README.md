# HaploTreeSim
**A Controlled, End-to-End Benchmark with Haplotype-Resolved CNA and CNA-Phylogeny Ground Truth for Low-Pass scDNA-seq**

[![Version](https://img.shields.io/badge/version-1.0.0-blue.svg)](https://github.com/zaaferahhamid2/haplotreesim/releases/tag/v1.0.0)
[![Paper Compliant](https://img.shields.io/badge/paper-compliant-brightgreen.svg)](#)

HaploTreeSim is a simulator for generating synthetic single-cell DNA sequencing (scDNA-seq) datasets with ground-truth haplotype-specific copy number alterations (CNAs) and clone phylogenies.

---

## 🎯 Current Status: v1.0.0 (Paper-Compliant)

**✅ Complete Implementation** of all features specified in the HaploTreeSim paper.

### Recent Updates
- **v1.0.0** (April 2026): Phase B - WGD placement model + doublet allelic model
- **v0.9.0** (April 2026): Phase A - Parameter updates for paper compliance

---

## 🚀 Key Features

### Clone Tree Generation
- **Beta-splitting tree** with controlled leaf proportions
- K extant clones with imbalanced topology
- Clone mixing proportions sum to 1

### CNA Evolution
- **Focal/Arm/Chromosome events** with configurable rates
- **Whole Genome Doubling (WGD)**: Early placement (Equations 18-19)
- **Haplotype-specific** copy number evolution
- **Overlapping events** via difference arrays

### Observation Models
- **Read-depth**: Negative Binomial with library size variation
- **Allelic counts**: Beta-Binomial with LOH floor and sequencing error
- **Doublet model**: CN-weighted allelic mixture (Equation 15)

### Ground Truth
- Haplotype-specific CN profiles per clone
- Clone tree topology and branch lengths
- CNA event history per edge
- **Ploidy per clone** for normalization
- Segment boundaries from breakpoints

---

## 📦 Installation
```bash
git clone https://github.com/zaaferahhamid2/haplotreesim.git
cd haplotreesim
pip install -e .
```

---

## 🎮 Quick Start
```python
from haplotreesim import SimulationConfig, HaploTreeSimulator

config = SimulationConfig(
    chromosome='chr1',
    bin_width=500000,
    num_clones=5,
    num_cells=200,
    lambda_events=2.0,
    prob_wgd=0.3,
    prob_doublet=0.05,
    random_seed=42
)

sim = HaploTreeSimulator(config)
read_counts, (alt_counts, ref_counts, total_counts) = sim.run()

# Access ground truth
for clone in sim.clones:
    print(f"Clone {clone.index}: ploidy={clone.ploidy:.2f}")
```

---

## 📊 Example Outputs

See `sample_outputs/` for visualization examples:
- Week 8: CN profiles with WGD
- Week 9: Read depth vs CN (ploidy normalization)
- Week 10: BAF patterns with doublets

---

## ⚙️ Key Parameters
```python
# Tree
num_clones: int = 5
alpha_tree: float = 0.5
beta_tree: float = 0.3

# Events
lambda_events: float = 1.5
gain_prob: float = 0.40        # 60% losses

# WGD & Contaminants
prob_wgd: float = 0.0
prob_doublet: float = 0.0
prob_phase_switch: float = 0.01

# Coverage
mean_library_size: float = 50000
snp_density: float = 1e-3
```

---

## 📖 Documentation

- [PHASE_B_SUMMARY.md](PHASE_B_SUMMARY.md) - WGD + doublet details
- [PROJECT_STATUS.md](PROJECT_STATUS.md) - Development roadmap
- [BRANCH_GUIDE.md](BRANCH_GUIDE.md) - Branch organization

---

## 🎯 Roadmap

- [x] Weeks 5-10: Core simulator
- [x] Phase A: Parameter updates
- [x] Phase B: WGD + doublets
- [ ] Week 11: Core metrics
- [ ] Weeks 12-28: Evaluation framework & baselines

---

## 📧 Contact

Repository: https://github.com/zaaferahhamid2/haplotreesim

---

**Version**: 1.0.0 | **Status**: ✅ Paper-Compliant | **Updated**: April 2026
