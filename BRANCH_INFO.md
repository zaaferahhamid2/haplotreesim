# Phase A Branch: Parameter Updates (v0.9.0)

This branch contains Phase A implementation work from the paper compliance updates.

## Changes in this branch:
- Updated `gain_prob` to 0.40 (60% losses predominate in solid tumors)
- Removed +1 from event count (pure Poisson distribution)
- Fixed focal `L_max` to 10 Mb (was variable)
- Added `epsilon_seq` = 0.001 (sequencing error contamination)
- Added `epsilon_floor` = 1e-6 (LOH floor to prevent degenerate Beta)
- Fixed haplotype blocks to arm-level (was segment-level)
- Implemented LOH floor + sequencing error in allelic model
- Added ploidy computation and output for normalization

## Files modified:
- `src/haplotreesim/data_models.py`
- `src/haplotreesim/simulator.py`
- `src/haplotreesim/event_generator.py`

## Tag: v0.9.0
## Status: ✅ Completed and merged to main
## Date: April 2026
