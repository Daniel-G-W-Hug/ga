# Angle Calculations for PGA Line Orientation

## Date: 2025-11-16

## Overview

This document explains the simplification of angle calculations in `active_bivt2dp.cpp`
from using two GA angle functions to a single `atan2` call, while maintaining identical
behavior for arrow rendering.

---

## Key Mathematical Concepts

### Geometric Lines vs. Bivectors

**Geometric Line (Undirected):**

- A line as a set of points has NO unique orientation
- The same geometric line appears twice in a full 360° rotation
- Example: horizontal line appears at q=0° and q=180°
- Can be uniquely identified by a single angle in range [0°, 180°]

**Bivector p∧q (Oriented):**

- Encodes BOTH the line AND an orientation (which direction the arrow points)
- The bivector `p∧q = -(q∧p)` flips when you swap the points
- Requires full orientation information for rendering arrows
- Needs angle range [-180°, 180°] or equivalent two-angle system

### Why We Need Orientation

For rendering in `active_bivt2dp.cpp`:

- Must draw an arrow from `p_from` to `p_to`
- The arrow head must consistently point toward the second control point `q`
- As the user drags control points, the arrow must maintain consistent orientation
- Therefore, we need the FULL bivector orientation, not just the geometric line

---

## Original Implementation (Two Angles)

### Calculation Method (orig.)

```cpp
auto const x_axis = bivec2dp{0, 1, 0};
auto const y_axis = bivec2dp{1, 0, 0};  // actually -y_axis_2dp
auto const phi_x = angle(x_axis, bvt);  // angle to x-axis
auto const phi_y = angle(y_axis, bvt);  // angle to y-axis
```

### Angle Ranges

- `phi_x ∈ [0°, 180°]` - angle between bivector and x-axis
- `phi_y ∈ [0°, 180°]` - angle between bivector and y-axis
- Both angles together uniquely determine the bivector orientation

### Example Angle Values

| q position | phi_x | phi_y | Geometric interpretation |
|------------|-------|-------|--------------------------|
| (3, 0)     | 0°    | 90°   | Horizontal line, arrow right |
| (0, 3)     | 90°   | 180°  | Vertical line, arrow up |
| (-3, 0)    | 180°  | 90°   | Horizontal line, arrow left |
| (0, -3)    | 90°   | 0°    | Vertical line, arrow down |

### Quadrant Detection (8 Cases)

The original code uses combinations of `phi_x` and `phi_y` to determine 8 cases:

1. **pos x-axis**: `phi_x ≈ 0°` AND `phi_y ≈ 90°`
2. **first quadrant**: `phi_x ∈ (0°, 90°)` AND `phi_y ∈ (90°, 180°)`
3. **pos y-axis**: `phi_x ≈ 90°` AND `phi_y ≈ 180°`
4. **second quadrant**: `phi_x ∈ (90°, 180°)` AND `phi_y ∈ (90°, 180°)`
5. **neg x-axis**: `phi_x ≈ 180°` AND `phi_y ≈ 90°`
6. **third quadrant**: `phi_x ∈ (90°, 180°)` AND `phi_y ∈ (0°, 90°)`
7. **neg y-axis**: `phi_x ≈ 90°` AND `phi_y ≈ 0°`
8. **fourth quadrant**: `phi_x ∈ (0°, 90°)` AND `phi_y ∈ (0°, 90°)`

### Intersection Calculation Conditions (orig.)

```cpp
// Calculate horizontal boundary intersections (left/right)
if (std::abs(phi_x - pi / 2.0) > eps) {  // NOT on ±y-axis
    p_from_24 = unitize(rwdg(l2, bvt));  // left boundary
    p_to_24 = unitize(rwdg(l4, bvt));    // right boundary
}

// Calculate vertical boundary intersections (top/bottom)
if (std::abs(phi_y - pi / 2.0) > eps) {  // NOT on ±x-axis
    p_from_13 = unitize(rwdg(l3, bvt));  // bottom boundary
    p_to_13 = unitize(rwdg(l1, bvt));    // top boundary
}
```

---

## New Implementation (Single atan2 Angle)

### Calculation Method (new)

```cpp
// For bivector (x,y,z), the line direction is perpendicular to normal (-y, x)
value_t line_angle = std::atan2(-bvt.y, bvt.x);
```

### Angle Range

- `line_angle ∈ [-180°, 180°]` - full orientation in atan2 range
- Single value uniquely determines the bivector orientation
- Simpler calculation: one `atan2` call instead of two GA `angle()` calls

### Complete Angle Mapping Table

| q angle | phi_x | phi_y | line_angle | Case | bvt components |
|---------|-------|-------|------------|------|----------------|
| 0°      | 0°    | 90°   | **-90°**   | pos x-axis | (0, +3, 0) |
| 15°     | 15°   | 105°  | **-105°**  | first quadrant | (-0.78, +2.90, 0) |
| 30°     | 30°   | 120°  | **-120°**  | first quadrant | (-1.50, +2.60, 0) |
| 45°     | 45°   | 135°  | **-135°**  | first quadrant | (-2.12, +2.12, 0) |
| 60°     | 60°   | 150°  | **-150°**  | first quadrant | (-2.60, +1.50, 0) |
| 75°     | 75°   | 165°  | **-165°**  | first quadrant | (-2.90, +0.78, 0) |
| 90°     | 90°   | 180°  | **±180°**  | pos y-axis | (-3, 0, 0) |
| 105°    | 105°  | 165°  | **+165°**  | second quadrant | (-2.90, -0.78, 0) |
| 120°    | 120°  | 150°  | **+150°**  | second quadrant | (-2.60, -1.50, 0) |
| 135°    | 135°  | 135°  | **+135°**  | second quadrant | (-2.12, -2.12, 0) |
| 150°    | 150°  | 120°  | **+120°**  | second quadrant | (-1.50, -2.60, 0) |
| 165°    | 165°  | 105°  | **+105°**  | second quadrant | (-0.78, -2.90, 0) |
| 180°    | 180°  | 90°   | **+90°**   | neg x-axis | (0, -3, 0) |
| 195°    | 165°  | 75°   | **+75°**   | third quadrant | (+0.78, -2.90, 0) |
| 210°    | 150°  | 60°   | **+60°**   | third quadrant | (+1.50, -2.60, 0) |
| 225°    | 135°  | 45°   | **+45°**   | third quadrant | (+2.12, -2.12, 0) |
| 240°    | 120°  | 30°   | **+30°**   | third quadrant | (+2.60, -1.50, 0) |
| 255°    | 105°  | 15°   | **+15°**   | third quadrant | (+2.90, -0.78, 0) |
| 270°    | 90°   | 0°    | **0°**     | neg y-axis | (+3, 0, 0) |
| 285°    | 75°   | 15°   | **-15°**   | fourth quadrant | (+2.90, +0.78, 0) |
| 300°    | 60°   | 30°   | **-30°**   | fourth quadrant | (+2.60, +1.50, 0) |
| 315°    | 45°   | 45°   | **-45°**   | fourth quadrant | (+2.12, +2.12, 0) |
| 330°    | 30°   | 60°   | **-60°**   | fourth quadrant | (+1.50, +2.60, 0) |
| 345°    | 15°   | 75°   | **-75°**   | fourth quadrant | (+0.78, +2.90, 0) |
| 360°    | 0°    | 90°   | **-90°**   | pos x-axis | (0, +3, 0) |

### Quadrant Detection (8 Cases with Single Angle)

Ordered from `line_angle = -180°` to `+180°`:

1. **pos y-axis**: `|abs(line_angle) - π| < eps` → `line_angle ≈ ±180°`
2. **first quadrant**: `line_angle ∈ (-180°, -90°)`
3. **pos x-axis**: `|line_angle + π/2| < eps` → `line_angle ≈ -90°`
4. **fourth quadrant**: `line_angle ∈ (-90°, 0°)`
5. **neg y-axis**: `|line_angle| < eps` → `line_angle ≈ 0°`
6. **third quadrant**: `line_angle ∈ (0°, +90°)`
7. **neg x-axis**: `|line_angle - π/2| < eps` → `line_angle ≈ +90°`
8. **second quadrant**: `line_angle ∈ (+90°, +180°)`

### Intersection Calculation Conditions (new)

```cpp
// Calculate horizontal boundary intersections (left/right)
// Skip when on ±y-axis
if (!(std::abs(std::abs(line_angle) - pi) < eps || std::abs(line_angle) < eps)) {
    p_from_24 = unitize(rwdg(l2, bvt));  // left boundary
    p_to_24 = unitize(rwdg(l4, bvt));    // right boundary
}

// Calculate vertical boundary intersections (top/bottom)
// Skip when on ±x-axis
if (!(std::abs(line_angle - pi / 2.0) < eps || std::abs(line_angle + pi / 2.0) < eps)) {
    p_from_13 = unitize(rwdg(l3, bvt));  // bottom boundary
    p_to_13 = unitize(rwdg(l1, bvt));    // top boundary
}
```

### Mapping from Original to New Conditions

| Original Condition | Maps to line_angle | Geometric Case |
|-------------------|-------------------|----------------|
| `phi_x ≈ 90°` | `line_angle ≈ ±180°` OR `line_angle ≈ 0°` | Skip p_24 (on ±y-axis) |
| `phi_y ≈ 90°` | `line_angle ≈ -90°` OR `line_angle ≈ +90°` | Skip p_13 (on ±x-axis) |

**Key insight for intersection skipping:**

- **Skip p_24** when `phi_x = 90°` → Skip when `line_angle ≈ ±180°` or `0°` (±y-axis)
- **Skip p_13** when `phi_y = 90°` → Skip when `line_angle ≈ ±90°` (±x-axis)

---

## Critical Edge Case: The ±180° Boundary

### The Problem

At `q = 90°` (positive y-axis), the bivector is `bvt = (-3, 0, 0)`, and:

```cpp
line_angle = atan2(-bvt.y, bvt.x) = atan2(0, -3)
```

The C++ standard allows `atan2(0, negative)` to return **either +180° or -180°**.

### The Solution

Use the condition:

```cpp
std::abs(std::abs(line_angle) - pi) < eps
```

This works for BOTH possible return values:

- If `atan2` returns `+π`: `|+π| - π = 0 < eps` ✓
- If `atan2` returns `-π`: `|-π| - π = 0 < eps` ✓

### Platform Behavior

On macOS ARM64 with clang++:

- `atan2(0, -3)` returns `+180.0°` (the positive value)
- Transition behavior:
  - `q = 89.95°` → `line_angle = -179.95°`
  - `q = 90.00°` → `line_angle = -180.00°` (wraps to negative π)
  - `q = 90.05°` → `line_angle = +179.95°`

The condition `|abs(line_angle) - π| < eps` is **platform-independent** because it treats
+π and -π as identical.

---

## Verification and Testing

### Comprehensive Test Results

**Test file:** `verify_complete_equivalence.cpp`

**Results:**

- ✅ 73/73 test cases passed (every 5° from 0° to 360°)
- ✅ All edge cases verified (0°, 90°, 180°, 270°)
- ✅ Perfect numerical match (coordinates identical to 1e-10 precision)
- ✅ Zero mismatches between original and new implementations

### Example Edge Case Verification

| Case | Original Output | New Output | Match |
|------|----------------|------------|-------|
| 0° (pos x-axis) | `(-5,0)→(+5,0)` | `(-5,0)→(+5,0)` | ✓ |
| 90° (pos y-axis) | `(0,-5)→(0,+5)` | `(0,-5)→(0,+5)` | ✓ |
| 180° (neg x-axis) | `(+5,0)→(-5,0)` | `(+5,0)→(-5,0)` | ✓ |
| 270° (neg y-axis) | `(0,+5)→(0,-5)` | `(0,+5)→(0,-5)` | ✓ |

---

## Advantages of the New Implementation

### Performance

- **One** `std::atan2()` call instead of **two** GA `angle()` calls
- Simpler, more efficient computation

### Code Clarity

- Single angle value is easier to reason about
- Direct geometric interpretation: angle of the bivector normal
- Fewer variables to track (`line_angle` instead of `phi_x` + `phi_y`)

### Correctness

- Mathematically equivalent to original (verified by testing)
- Platform-independent edge case handling
- No loss of precision or accuracy

---

## Important Notes

### What "Quadrants" Mean Here

The 8 cases are NOT classical mathematical quadrants! They are **boundary crossing
patterns**:

- **"first quadrant"**: line crosses from bottom-left to top-right
- **"second quadrant"**: line crosses from bottom-right to top-left
- **"third quadrant"**: line crosses from top-left to bottom-right
- **"fourth quadrant"**: line crosses from top-right to bottom-left

These names are legacy from the original implementation and represent which viewport
boundaries the line intersects, rotated 45° from classical quadrants.

### Unique Angle for Geometric Lines

If we only needed to identify the **undirected geometric line**:

- Could use `phi_x ∈ [0°, 180°]` as a unique identifier
- Each geometric line appears exactly once in this range

But for **rendering with arrows**:

- Need full bivector orientation (which direction the arrow points)
- Must distinguish `p∧q` from `q∧p`
- Requires either `(phi_x, phi_y)` or `line_angle ∈ [-180°, 180°]`

---

## Implementation Details

### File Modified

- `ga_view/src/active_bivt2dp.cpp` (lines 204-324)

### Key Changes

**Removed:**

```cpp
auto const x_axis = bivec2dp{0, 1, 0};
auto const y_axis = bivec2dp{1, 0, 0};
auto const phi_x = angle(x_axis, bvt);
auto const phi_y = angle(y_axis, bvt);
```

**Added:**

```cpp
// Determine the angle of the projective line using atan2
// For bivector (x,y,z), the line direction is (-y, x)
// This gives us the full orientation in range [-π, π]
value_t line_angle = std::atan2(-bvt.y, bvt.x);
```

**All 8 case conditions updated** to use `line_angle` instead of `(phi_x, phi_y)`
combinations.

**Intersection calculation conditions updated** with correct atan2 range mappings.

---

## Conclusion

The simplification from two angles to one atan2-based angle is:

✅ **Mathematically equivalent** - verified by comprehensive testing
✅ **More efficient** - fewer function calls
✅ **Clearer** - single angle easier to understand
✅ **Platform-independent** - robust edge case handling
✅ **Production-ready** - zero behavioral changes

This refactoring maintains 100% identical behavior while improving code clarity and
performance.

---

## References

### Test Files (in this directory)

All test files have been moved to `TODO/LL_angle_calc/` for archival purposes:

- `angle_mapping.cpp` - Maps phi_x/phi_y to line_angle
- `verify_mapping.cpp` - Verifies case classification mapping
- `verify_edge_cases.cpp` - Tests edge cases (0°, 90°, 180°, 270°)
- `test_atan2_boundary.cpp` - Investigates ±180° boundary behavior
- `verify_complete_equivalence.cpp` - **Complete output equivalence test**
- `complete_angle_table.cpp` - Full angle table with analysis
- `systematic_angles.cpp` - Systematic angle walkthrough
- `unique_line_angle.cpp` - Analysis of unique line angles
- `test_edge_conditions.cpp` - Maps original conditions to atan2 ranges
- `test_bivt_orientation.cpp` - Tests bivector orientation with rotation

To compile and run any test (from project root):

```bash
cd build
clang++ -std=c++23 -I.. -I/opt/homebrew/include ../TODO/LL_angle_calc/test_file.cpp -L/opt/homebrew/lib -lfmt -o test_file
./test_file
```

### Related Code

- Original implementation: `git show HEAD:ga_view/src/active_bivt2dp.cpp`
- Current implementation: `ga_view/src/active_bivt2dp.cpp` (lines 204-324)
