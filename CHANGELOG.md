# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned for Week 6
- CNA event generation (gains, losses, WGD)
- Segment boundary detection from breakpoints
- Non-trivial clone tree structure
- Non-uniform clone frequencies

## [0.1.0] - 2024-XX-XX

### Added (Week 5 Deliverable)
- Initial release of HaploTreeSim simulator
- Core data models:
  - `Bin`: Fixed-width genomic regions
  - `Segment`: Contiguous bin intervals
  - `HaplotypeBlock`: Groups of segments with shared phasing
  - `CNAEvent`: Copy number alteration events
  - `Clone`: Tree nodes with haplotype-specific CN profiles
  - `Cell`: Individual simulated cells
  - `SimulationConfig`: Complete parameter configuration
- `HaploTreeSimulator` class:
  - Genome initialization (bins, segments, haplotype blocks)
  - Diploid clone tree generation (no CNAs yet)
  - Cell sampling with realistic library sizes
  - Read-depth observation model (Negative Binomial)
  - Allelic observation model (Beta-Binomial)
  - Ground truth extraction
- Configuration system:
  - Example configurations for diploid and WGD roots
  - Low-pass coverage scenarios (0.01x)
- Test suite:
  - Integration tests for diploid simulation
  - WGD root validation
  - Output shape and correctness verification
- Documentation:
  - Comprehensive README with installation and usage
  - Example usage scripts
  - GitHub setup guide
- Repository structure:
  - Clean package layout (`src/`, `tests/`, `configs/`)
  - Installation via `setup.py`
  - Dependencies in `requirements.txt`

### Technical Details
- Implements observation models from HaploTreeSim paper (Sections 3.1-3.6)
- Generates realistic low-pass scDNA-seq data with proper noise
- Produces both read-count and allele-count matrices
- Supports diploid (1,1) and WGD (2,2) root initialization
- All clones currently diploid (CNAs coming in Week 6)

### Dependencies
- Python >= 3.8
- numpy >= 1.20.0
- scipy >= 1.7.0

## Version Roadmap

- **v0.1.0** (Week 5): Minimal diploid simulator ✅
- **v0.2.0** (Week 10): CNA events and tree structure
- **v0.3.0** (Week 15): Read mode with FASTQ generation
- **v0.4.0** (Week 20): Benchmark tracks and metrics
- **v1.0.0** (Week 28): Complete benchmark suite

---

[Unreleased]: https://github.com/YOUR_USERNAME/haplotreesim/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/YOUR_USERNAME/haplotreesim/releases/tag/v0.1.0
