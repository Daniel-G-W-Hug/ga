# STA4D Extended Metric Calculation - Potential Issue

**Date**: 2026-01-06
**Status**: Needs Investigation
**Priority**: Low (current implementation works correctly, but may be overly complex)

## Current Situation

STA4D (spacetime algebra) requires **special-case handling** in the extended metric calculation that differs from all other algebras (EGA2D, EGA3D, PGA2DP, PGA3DP). This special case is hardcoded with a literal metric signature check.

**Code Location**: [ga_prdxpr_rule_generator.cpp:252-265](../ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp#L252-L265)

```cpp
if (config.metric_signature == std::vector<int>{+1, -1, -1, -1}) {
    // STA4D special case: check if g0 (index 0) is present
    bool has_g0 = std::find(indices.begin(), indices.end(), 0) != indices.end();

    if (indices.size() == 2) {
        // Bivectors: +1 if contains g0, -1 if only spatial indices
        metric_value = has_g0 ? +1 : -1;
    }
    else if (indices.size() == 3) {
        // Trivectors: -1 if contains g0, +1 if only spatial indices
        metric_value = has_g0 ? -1 : +1;
    }
}
```

## The Mathematical Problem

For most algebras (EGA, PGA), the extended metric follows a simple rule:

- **Product of constituent metrics**: `G(e_i ∧ e_j) = G(e_i) × G(e_j)`

### Example: EGA3D with metric `{+1, +1, +1}`

- All vectors: `e1² = e2² = e3² = +1`
- All bivectors: `e23² = e31² = e12² = +1` (product of two +1s)
- Pseudoscalar: `e123² = +1` (product of three +1s)

**This simple product rule works perfectly.**

### Example: STA4D G(1,3,0) with metric `{+1, -1, -1, -1}`

**Vectors:**

- `g0² = +1` (timelike)
- `g1² = g2² = g3² = -1` (spacelike)

**Bivectors** - the naive product rule gives **wrong answers**:

| Bivector | Naive Product          | Actual Value | Match?        |
|----------|------------------------|--------------|---------------|
| `g01²`   | `(+1) × (-1) = -1`     | `+1`         | ❌ **WRONG**  |
| `g02²`   | `(+1) × (-1) = -1`     | `+1`         | ❌ **WRONG**  |
| `g03²`   | `(+1) × (-1) = -1`     | `+1`         | ❌ **WRONG**  |
| `g23²`   | `(-1) × (-1) = +1`     | `-1`         | ❌ **WRONG**  |
| `g31²`   | `(-1) × (-1) = +1`     | `-1`         | ❌ **WRONG**  |
| `g12²`   | `(-1) × (-1) = +1`     | `-1`         | ❌ **WRONG**  |

The pattern that emerges:

- **Bivectors with g0 (timelike)**: extended metric = `+1`
- **Bivectors without g0 (spatial)**: extended metric = `-1`

**Trivectors** follow the opposite pattern:

- `g023², g031², g012²` (contain g0): extended metric = `-1`
- `g123²` (only spatial): extended metric = `+1`

**Pseudoscalar:**

- `g0123² = (+1) × (-1) × (-1) × (-1) = -1` ✓ **CORRECT**

## Current Solution Issues

1. **Hardcoded metric signature check**: Uses literal comparison `{+1, -1, -1, -1}`
2. **Index-based logic**: Assumes g0 is at index 0
3. **Special case handling**: Completely different code path for STA4D
4. **Not generalizable**: Won't work for other mixed-signature algebras

## Why This Matters

The extended metric is used in:

1. **Dot product calculations** - determines inner product values
2. **Regressive products** - affects dual operations
3. **Contraction operations** - left/right contractions depend on these values

Getting it wrong produces mathematically incorrect results, breaking physical interpretations in relativity applications.

## Questions to Investigate

1. **Is there a general mathematical rule** that covers STA4D without special cases?
   - Perhaps involving signed/unsigned components?
   - Some property of mixed signatures that we're missing?

2. **What's the actual mathematical definition** of extended metric for mixed-signature algebras?
   - Check Hestenes' work on spacetime algebra
   - Look at Dorst/Fontijne/Mann's GA references
   - Review Clifford algebra literature

3. **Can we detect this pattern automatically** instead of hardcoding?
   - Detect "mixed signature" algebras (mix of +1 and -1)
   - Identify which dimension is "special" (the one with opposite sign)
   - Apply appropriate transformation rules

4. **Would a full metric matrix solve this elegantly?**
   - Instead of just diagonal values, store full metric tensor
   - Extended metric becomes matrix operations
   - But: adds complexity for simple cases (EGA, PGA)

## Potential Solutions to Explore

### Option 1: Find the General Mathematical Rule

- Research the actual definition used in physics/GA literature
- Implement the mathematically correct formula
- Remove special cases entirely

### Option 2: Generalize the Detection

Instead of:

```cpp
if (config.metric_signature == std::vector<int>{+1, -1, -1, -1})
```

Do:

```cpp
if (has_mixed_signature(config.metric_signature)) {
    auto special_index = find_minority_signature_index(config.metric_signature);
    // Apply generalized transformation based on special_index presence
}
```

### Option 3: Configuration-Driven Approach

Add to `AlgebraConfig`:

```cpp
struct AlgebraConfig {
    // ... existing fields ...

    // Optional: custom extended metric calculation strategy
    enum class ExtendedMetricMode {
        Standard,        // Simple product rule (EGA, PGA)
        MixedSignature,  // STA-style mixed signature handling
        CustomMatrix     // Full metric matrix (future CGA)
    };
    ExtendedMetricMode extended_metric_mode = Standard;

    // For MixedSignature mode: which index triggers special handling
    std::optional<int> special_signature_index;
};
```

### Option 4: Accept the Complexity

- Current implementation is correct
- STA4D truly is a special case mathematically
- Document clearly and move on
- Wait for more mixed-signature algebras to see patterns

## Next Steps

1. **Literature review**: Find authoritative GA references on extended metrics for mixed signatures
2. **Mathematical analysis**: Derive the correct general formula (if it exists)
3. **Proof of concept**: Implement a cleaner solution and verify it produces identical results
4. **Decision**: Keep current approach vs. refactor based on findings

## References to Add

- [ ] Hestenes, D. "Space-Time Algebra"
- [ ] Doran & Lasenby "Geometric Algebra for Physicists"
- [ ] Dorst, Fontijne, Mann "Geometric Algebra for Computer Science"
- [ ] Relevant Clifford algebra papers on metric tensors

## Related Code Locations

- **Extended metric calculation**: [ga_prdxpr_rule_generator.cpp:187-280](../ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp#L187-L280)
- **Metric lookup**: [ga_prdxpr_metric_calc.cpp:139-140](../ga_prdxpr/src_prdxpr/ga_prdxpr_metric_calc.cpp#L139-L140)
- **STA4D basis definition**: [ga_prdxpr_sta4d.hpp:16-35](../ga_prdxpr/src_prdxpr/ga_prdxpr_sta4d.hpp#L16-L35)
- **Algebra config struct**: [ga_prdxpr_rule_generator.hpp:87](../ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.hpp#L87)

---

**Note**: The current implementation is **mathematically correct** and **fully tested**. This is an optimization/elegance concern, not a bug.
