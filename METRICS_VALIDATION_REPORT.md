# Metrics Validation Report

**Date**: May 12, 2026  
**Purpose**: Verify metrics work correctly before testing on real tools

## Test Setup

Generated synthetic ground truth data and tested two scenarios:
1. **Perfect Predictions**: Predictions = Ground Truth
2. **Corrupted Predictions**: Manually introduced errors

## Results

### 1. HSCN Metrics ✓
- **Perfect**: HSCN error = 0.0, TCN MSE = 0.0, LOH F1 = 1.0, MSR = 1.0
- **Corrupted** (10% CN changed): HSCN error = 0.058, LOH F1 = 0.98
- ✅ Correctly detects errors

### 2. Breakpoint Metrics ✓
- **Perfect**: Precision = 1.0, Recall = 1.0, F1 = 1.0
- **Corrupted** (2 missed, 1 FP): Precision = 0.75, Recall = 0.60, F1 = 0.67
- ✅ Correctly detects FP and FN

### 3. Clone Assignment Metrics ✓
- **Perfect**: ARI = 1.0, NMI = 1.0, Hungarian Acc = 1.0
- **Corrupted** (20 cells misassigned): ARI = 0.64, NMI = 0.69, Acc = 0.81
- ✅ Correctly detects misassignments

### 4. Tree Metrics ✓
- **Perfect**: RF = 0.0, Ancestor-Desc Acc = 1.0, Event F1 = 1.0
- **Corrupted** (different tree): Event Recall = 0.50, F1 = 0.67
- ✅ Correctly detects topology and event differences

## Conclusion

✅ All metrics validated and working correctly:
- Perfect data → perfect scores
- Corrupted data → appropriately degraded scores

**Status**: Ready to test on real baseline methods

## Next Steps

Week 14: Benchmark pipeline to run real methods (CHISEL, SEACON, etc.)
