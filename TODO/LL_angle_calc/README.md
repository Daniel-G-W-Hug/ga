# Lessons Learned: Angle Calculation Simplification

This directory contains documentation and test files from the refactoring of `active_bivt2dp.cpp` to simplify angle calculations from a two-angle system to a single `atan2`-based approach.

## Main Documentation

**[angle_calc.md](angle_calc.md)** - Complete documentation covering:

- Mathematical concepts (geometric lines vs. bivectors)
- Original two-angle implementation
- New single-angle implementation with `atan2`
- Edge case handling (±180° boundary)
- Comprehensive verification results
- Complete angle mapping tables

## Test Files

All test programs can be compiled from the project root:

```bash
cd build
clang++ -std=c++23 -I.. -I/opt/homebrew/include ../TODO/LL_angle_calc/<filename>.cpp -L/opt/homebrew/lib -lfmt -o <filename>
./<filename>
```

### Core Verification Tests

- **verify_complete_equivalence.cpp** - Most important: Compares actual p_from/p_to outputs between original and new implementations (73/73 tests passed ✓)
- **verify_mapping.cpp** - Verifies case classification mapping (25/25 tests passed ✓)
- **verify_edge_cases.cpp** - Tests critical edge cases at 0°, 90°, 180°, 270°

### Analysis Tools

- **angle_mapping.cpp** - Maps phi_x/phi_y values to line_angle for key angles
- **complete_angle_table.cpp** - Full table showing all angles every 15° from 0° to 360°
- **systematic_angles.cpp** - Systematic walkthrough with geometric interpretation
- **unique_line_angle.cpp** - Analysis of unique angles for geometric lines

### Edge Case Investigation

- **test_atan2_boundary.cpp** - Investigates ±180° boundary behavior
- **test_edge_conditions.cpp** - Maps original phi_x/phi_y conditions to atan2 ranges
- **test_bivt_orientation.cpp** - Tests bivector orientation with rotating control points

## Summary

**Result:** Successfully simplified from 2 GA `angle()` calls to 1 `std::atan2()` call with:

- ✅ 100% identical behavior (verified by comprehensive testing)
- ✅ Better performance (fewer function calls)
- ✅ Improved code clarity (single angle value)
- ✅ Platform-independent edge case handling

**Modified file:** `ga_view/src/active_bivt2dp.cpp` (lines 204-324)

**Date:** November 2025
