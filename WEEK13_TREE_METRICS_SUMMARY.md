# Week 13: Tree Metrics

**Status**: ✅ COMPLETE  
**Date**: May 12, 2026

## Implemented Metrics

### 1. Robinson-Foulds (RF) Distance
Measures topological dissimilarity between two trees using bipartitions.
- Normalized RF ∈ [0, 1]
- 0 = identical topology, 1 = maximally different

### 2. Ancestor-Descendant Accuracy
Checks if pairwise ancestor-descendant relationships are preserved.
- Compares lowest common ancestor (LCA) relationships
- Returns fraction of correct relationships

### 3. Event Matching
Matches predicted CNA events to true events using:
- Interval overlap (with tolerance)
- Haplotype label agreement
- Returns precision/recall/F1

## Test Results: 9/9 passing ✓

## Next: Week 14 - Benchmark Pipeline
