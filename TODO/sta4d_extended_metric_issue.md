# STA4D Extended Metric Calculation - Potential Issue

**Date**: 2026-01-06
**Status**: Needs Investigation
**Priority**: Low (current implementation works correctly, but may be overly complex)

PLEASE REMEMBER: everything concerning STA currently is work in progress, not tested and may simply still be WRONG.

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

## Investigation Session 2026-01-06: Gram Determinant Analysis

### Discovery: Two Calculation Methods Produce Different Results

There are **two separate extended metric calculation functions**:

1. **`calculate_extended_metric()`** in `ga_prdxpr_rule_generator.cpp` (lines 188-281)
   - Diagonal-only calculation
   - Has hardcoded STA4D special case (lines 253-266)
   - **Produces correct results** for STA4D bivectors/trivectors

2. **`calculate_extended_metric_matrix_full()`** in `ga_prdxpr_metric_calc.cpp`
   - Full matrix calculation using Gram determinants
   - Mathematically rigorous approach: `G[i,j] = ⟨basis[i], basis[j]⟩ = det(Gram)`
   - **Produces WRONG results** for STA4D (signs flipped for bivectors)

### Actual Test Output for STA4D

**From diagonal calculation (correct):**

```text
g01 → +1 ✓
g02 → +1 ✓
g03 → +1 ✓
g23 → -1 ✓
g31 → -1 ✓
g12 → -1 ✓
g012 → -1 ✓
g013 → -1 ✓
g023 → -1 ✓
g123 → +1 ✓
g0123 → -1 ✓
```

**From full matrix calculation (WRONG - signs flipped for bivectors):**

```text
g01 → -1 ❌ (should be +1)
g02 → -1 ❌ (should be +1)
g03 → -1 ❌ (should be +1)
g23 → +1 ❌ (should be -1)
g31 → +1 ❌ (should be -1)
g12 → +1 ❌ (should be -1)
```

### The Gram Determinant Problem

For orthogonal bases, the Gram matrix is diagonal:

```text
Gram[i,j] = metric[indices[i]] if i==j, else 0
det(Gram) = ∏ metric[indices]
```

**EGA3D `e23`** (works correctly):

```text
Gram = [+1  0 ]  → det = +1 ✓
       [ 0 +1 ]
```

**STA4D `g01`** (wrong sign):

```text
Gram = [+1  0 ]  → det = -1 ❌ (should be +1)
       [ 0 -1 ]
```

**STA4D `g23`** (wrong sign):

```text
Gram = [-1  0 ]  → det = +1 ❌ (should be -1)
       [ 0 -1 ]
```

**Pattern**: For STA4D, the Gram determinant gives **exactly the opposite sign** from correct!

### Mathematical Observations

**Anticommutativity insight**: When squaring a blade, reordering is needed:

```text
g01² = (g0∧g1)² = (g0·g1)·(g0·g1)
     = g0·(g1·g0)·g1
     = g0·(-g0·g1)·g1    (anticommute: g1·g0 = -g0·g1)
     = -(g0·g0)·(g1·g1)
     = -(+1)·(-1) = +1 ✓
```

This suggests the formula: `G(blade) = (-1)^(k(k-1)/2) × ∏ metric[indices]`

But testing against EGA3D:

```text
e23: (-1)^(2×1/2) × (+1)×(+1) = -1 × 1 = -1 ❌
```

This would break EGA3D (which currently works), so this formula is wrong.

### Hypothesis: Minkowski Signature Requires Sign Flip

For **Minkowski-like signatures** (one timelike among spacelike, or vice versa):

- The Gram determinant produces the **negation** of the correct extended metric
- Need to detect this case and flip sign

But **why mathematically**? Possibilities:

1. Related to how signed volumes work in pseudo-Riemannian spaces?
2. Convention difference between Euclidean and Minkowski inner products?
3. Different interpretation of the exomorphism matrix for mixed signatures?

### Open Questions

1. **What is the mathematical principle** that explains why Gram determinants work for uniform signatures but fail (with sign flip) for mixed signatures?

2. **Reference needed**: Does the rigidgeometricalgebra.org definition of the exomorphism matrix explicitly address mixed signatures?

3. **Is this a convention issue** or a fundamental mathematical difference in how inner products work in pseudo-Riemannian spaces?

4. **Should we**:
   - Keep using `calculate_extended_metric()` (with STA4D special case) as the source of truth?
   - Fix `calculate_extended_metric_matrix_full()` to handle Minkowski signatures?
   - Replace Gram determinant approach entirely for mixed signatures?

### Immediate Action Taken

- Fixing inconsistency between diagonal and full matrix calculations
- Diagonal calculation (with STA4D special case) is considered correct
- Full matrix calculation needs to match diagonal results

### References for Further Investigation

- Hestenes "Space-Time Algebra" - Chapter on metric and inner products
- Doran & Lasenby "Geometric Algebra for Physicists" - Spacetime algebra metrics
- <https://rigidgeometricalgebra.org/wiki/index.php?title=Metrics> - Check for Minkowski case
- Clifford algebra literature on bilinear forms in pseudo-Riemannian spaces

---

## Investigation Session 2026-01-06 (Part 2): PDF Analysis - Exomorphism Matrix Theory

### Key Discoveries from "Geometric Algebra for Computer Graphics" (Lengyel)

**Source**: Chapter 2, Section 2.7 "Exomorphisms" and Section 2.8 "Metric Transformations"

#### 1. The Conforming Property (Equation 2.49)

> "For any basis bivector **eᵢⱼ**, Equation (2.48) requires that **Meᵢⱼ = (Meᵢ) ∧ (Meⱼ)**"

This is the **defining property** of the exomorphism matrix **M**. The extended metric matrix **G** must satisfy:

```text
G(a ∧ b) = (Ga) ∧ (Gb)  for all multivectors a, b
```

For bivector `eᵢⱼ = eᵢ ∧ eⱼ`:

```text
G(eᵢⱼ) = G(eᵢ) ∧ G(eⱼ) = metric[i]·eᵢ ∧ metric[j]·eⱼ
       = metric[i]·metric[j]·(eᵢ ∧ eⱼ)
```

**Therefore**: The extended metric value for `eᵢⱼ` should be `metric[i] × metric[j]`

**This is exactly what the Gram determinant computes!** Yet it fails for STA4D.

#### 2. The Antimetric Relationship (Equation 2.65)

```text
G · 𝔾 = det(g) · I
```

Where:

- **G** = metric exomorphism matrix
- **𝔾** = metric antiexomorphism matrix (related to complement)
- **g** = base metric tensor
- **I** = identity matrix

**Critical insight for STA4D**:

- det(g) = (+1)·(-1)·(-1)·(-1) = **-1**
- Therefore: **G · 𝔾 = -I**

For Euclidean algebras (EGA):

- det(g) = (+1)ⁿ = **+1**
- Therefore: **G · 𝔾 = +I**

**Hypothesis**: The negative determinant for Minkowski signatures might explain why the Gram approach produces opposite signs!

#### 3. Metric Transformation Matrix Structure (Figure 2.11, Equation 2.61)

The extended metric **G** is constructed as a **block diagonal matrix**:

- **C₀(m)** = 1 (scalar identity)
- **C₁(m)** = **m** (the base metric matrix)
- **C₂(m)** = wedge products of columns of **m** (bivector compound matrix)
- **Cₙ(m)** = determinant of **m** (pseudoscalar)

The **compound matrices Cₖ(m)** are constructed via wedge products, which must satisfy the conforming property.

### Mathematical Puzzle

**The Contradiction**:

1. **PDF states**: Extended metric follows conforming property → `G(eᵢⱼ) = metric[i] × metric[j]`
2. **Gram determinant computes**: Exactly this product for orthogonal bases
3. **Reality for STA4D**: This gives **wrong signs** (opposite of correct values)

**Possible Resolutions**:

1. **Convention difference**: Perhaps the PDF uses a different sign convention for Minkowski metrics?
2. **Gram determinant limitation**: The inner product definition `⟨blade, blade⟩ = det(Gram)` may not apply to pseudo-Riemannian spaces?
3. **Missing sign factor**: The conforming property might need an additional sign factor `(-1)^something` for mixed signatures?
4. **Antimetric relationship**: The equation `G·𝔾 = det(g)·I` with negative determinant might imply a sign flip in the extended metric calculation?

### Action Items

1. ✅ **Added extended metric tests** for EGA3D to verify correctness
   - Tests confirm EGA3D extended metric matches expected values
   - Conforming property verified: G(e23) = (+1)·(+1) = +1 ✓

2. **TODO**: Add similar tests for all algebras (EGA2D, PGA2DP, PGA3DP)

3. **TODO**: Research question - Look for literature on:
   - Extended metric (exomorphism matrix) for Minkowski/pseudo-Riemannian signatures
   - How Clifford algebra inner products work in mixed-signature spaces
   - Whether the Gram determinant formula needs modification for non-Euclidean metrics

4. **TODO**: Investigate the antimetric relationship `G·𝔾 = det(g)·I`
   - Can we compute 𝔾 for STA4D and verify this equation?
   - Does the negative determinant cause a sign flip somewhere in the calculation?

### Related Code Locations

- **Extended metric tests**: [ga_ega3d_test.hpp:2916-2984](../ga_test/src/ga_ega3d_test.hpp#L2916-L2984)
- **Stored metric matrices**: [ga_usr_consts.hpp:94-111](../ga/ga_usr_consts.hpp#L94-L111)

---

**Status Update**: The diagonal calculation is **mathematically correct** for STA4D. The Gram determinant approach produces results consistent with the **conforming property** from the PDF, yet these are **wrong for STA4D**. This suggests either:

1. A different sign convention for Minkowski metrics in GA literature, or
2. The Gram determinant formula for inner products doesn't directly apply to pseudo-Riemannian spaces

The negative metric determinant (det(g) = -1) and its relationship to the antimetric (`G·𝔾 = -I`) may be key to understanding this discrepancy.

---

## Investigation Session 2026-01-06 (Part 3): Antimetric Relationship Testing

### Objective

Test the antimetric relationship from equation 2.65 of the PDF:

```text
G · 𝔾 = det(g) · I
```

Where:

- **G** = metric exomorphism matrix (extended metric)
- **𝔾** = metric antiexomorphism matrix
- **det(g)** = determinant of the base metric tensor
- **I** = identity matrix

### Implementation Approach

Added helper functions to [ga_prdxpr_rule_generator_test.cpp:114-240](../ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator_test.cpp#L114-L240):

1. **`build_antimetric_matrix()`**: Attempts to construct 𝔾 from complement rules
   - Maps each basis element to its complement
   - Builds matrix where 𝔾[i,j] is the coefficient mapping basis[i] → complement(basis[i])

2. **`multiply_matrices()`**: Standard matrix multiplication for verification

3. **`investigate_antimetric_relationship()`**: Tests the equation `G·𝔾 = det(g)·I`
   - Computes extended metric matrix **G**
   - Constructs antiexomorphism matrix **𝔾** from complement rules
   - Multiplies **G·𝔾** and compares with **det(g)·I**

### Results: The Relationship Does NOT Hold

**Testing EGA3D (det(g) = +1, expecting G·𝔾 = +I):**

```text
=== Antimetric Relationship Investigation: ega3d ===
Metric signature: 1, 1, 1
det(g) = 1
Algebra dimension: 3 (odd) - using complement

Verifying: G · 𝔾 = det(g) · I
Expected diagonal value: 1
Mismatches found:
  [0,0]: got 0, expected 1    (scalar → e123)
  [0,7]: got 1, expected 0
  [1,1]: got 0, expected 1    (e1 → e23)
  [1,4]: got 1, expected 0
  [2,2]: got 0, expected 1    (e2 → e31)
  [2,5]: got 1, expected 0
  [3,3]: got 0, expected 1    (e3 → e12)
  [3,6]: got 1, expected 0
  [4,1]: got 1, expected 0    (e23 → e1)
  [4,4]: got 0, expected 1
  [5,2]: got 1, expected 0    (e31 → e2)
  [5,5]: got 0, expected 1
  [6,3]: got 1, expected 0    (e12 → e3)
  [6,6]: got 0, expected 1
  [7,0]: got 1, expected 0    (e123 → scalar)
  [7,7]: got 0, expected 1
✗ FAILED: G · 𝔾 ≠ det(g) · I for ega3d
```

**Pattern Analysis**: The result `G·𝔾` is a **permutation matrix** (off-diagonal 1s, diagonal 0s), not the identity matrix!

### Key Finding: Complement ≠ Antimetric

**The complement operation is NOT the antiexomorphism 𝔾 from the PDF.**

The complement rules map:

- EGA3D: `1 ↔ e123`, `e1 ↔ e23`, `e2 ↔ e31`, `e3 ↔ e12`

This creates a **permutation-like structure**, not an antiexomorphism that satisfies `G·𝔾 = det(g)·I`.

### Mathematical Interpretation

**What is the antiexomorphism 𝔾?**

From the PDF context (equations 2.63-2.65), the antiexomorphism should satisfy:

```text
𝔾(a ∧ b) = (𝔾b) ∧ (𝔾a)   (note reversed order)
```

The complement operation in GA is related to duality:

```text
complement(A) = A · I⁻¹   (or similar dual operation)
```

These are **different mathematical operations**:

- **Complement**: Hodge dual, maps k-vectors to (n-k)-vectors
- **Antiexomorphism 𝔾**: A specific linear transformation satisfying G·𝔾 = det(g)·I

### Implications

1. **Cannot verify antimetric relationship** using complement rules directly
2. **Need to understand** how to properly construct the antiexomorphism matrix 𝔾
3. **The PDF's 𝔾** likely involves metric-weighted duals or more complex transformations
4. **Complement is a permutation** in orthonormal bases, not a metric transformation

### Open Questions

1. **How to construct 𝔾?** What is the actual definition of the metric antiexomorphism in terms of basis elements?

2. **Is 𝔾 related to the inverse complement?** Perhaps `𝔾 = metric-weighted complement`?

3. **Does equation 2.65 help explain STA4D?** If we could compute the correct 𝔾, would the relationship `G·𝔾 = -I` (for STA4D) reveal why signs flip?

4. **Literature search needed**: Find explicit construction of the antiexomorphism matrix in GA references

### Conclusion

The investigation revealed that:

- The complement operation creates a permutation matrix structure
- This is fundamentally different from the antiexomorphism 𝔾 in equation 2.65
- Further mathematical research is needed to understand the proper construction of 𝔾
- The antimetric relationship remains an intriguing but unverified hypothesis for explaining STA4D behavior

**Related Code**: [ga_prdxpr_rule_generator_test.cpp:931-939](../ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator_test.cpp#L931-L939)

---

## Investigation Session 2026-01-06 (Part 4): Recursive Extraction Validation

### Objective

Validate that the recursive extraction approach (building extended metric from wedge products level-by-level) works correctly for all implemented algebras before attempting to apply it to STA4D.

### Methodology: Recursive Extraction from Algebra Operations

Instead of computing extended metric from Gram determinants (which fails for STA4D), we extract it directly from the algebra's own operations:

**Algorithm**:
1. **Level 0 (Scalar)**: Extended metric element is always 1 by definition
2. **Level 1 (Vectors)**: Extract from dot products: `g_ei = dot(ei, ei)`
3. **Level k (k≥2)**: Construct basis element from wedge products of lower-level elements, then extract: `g_basis = dot(basis_constructed, basis_constructed)`

**Key Insight**: This approach uses the *conforming property* of the exomorphism: G(a∧b) = (Ga)∧(Gb). We build higher-grade basis elements via wedge products, then extract their metric values.

### Implementation: Test Cases Added

#### EGA2D - G(2,0,0)
**File**: [ga_ega2d_test.hpp:2503-2528](../ga_test/src/ga_ega2d_test.hpp#L2503-L2528)

```cpp
TEST_CASE("G<2,0,0>: extended metric recursive extraction via wedge products")
{
    auto G = ega2d_metric_view();

    // Level 0: Scalar (always 1)
    CHECK(G[0, 0] == 1);

    // Level 1: Vectors (extract from dot products)
    value_t g_e1 = value_t(dot(e1_2d, e1_2d));    // = 1
    value_t g_e2 = value_t(dot(e2_2d, e2_2d));    // = 1
    CHECK(abs(g_e1 - value_t(G[1, 1])) < eps);
    CHECK(abs(g_e2 - value_t(G[2, 2])) < eps);

    // Level 2: Bivectors (extract from wedge products + dot)
    auto e12_constructed = wdg(e1_2d, e2_2d);
    value_t g_e12 = value_t(dot(e12_constructed, e12_constructed));  // = 1
    CHECK(abs(g_e12 - value_t(G[3, 3])) < eps);
}
```

**Result**: ✓ All checks pass

#### EGA3D - G(3,0,0)
**File**: [ga_ega3d_test.hpp:2964-3007](../ga_test/src/ga_ega3d_test.hpp#L2964-L3007)

```cpp
TEST_CASE("G<3,0,0>: extended metric recursive extraction via wedge products")
{
    auto G = ega3d_metric_view();

    // Level 0: Scalar (always 1)
    CHECK(G[0, 0] == 1);

    // Level 1: Vectors
    value_t g_e1 = value_t(dot(e1_3d, e1_3d));    // = 1
    value_t g_e2 = value_t(dot(e2_3d, e2_3d));    // = 1
    value_t g_e3 = value_t(dot(e3_3d, e3_3d));    // = 1
    CHECK(abs(g_e1 - value_t(G[1, 1])) < eps);
    CHECK(abs(g_e2 - value_t(G[2, 2])) < eps);
    CHECK(abs(g_e3 - value_t(G[3, 3])) < eps);

    // Level 2: Bivectors (constructed from wedge products)
    auto e23_constructed = wdg(e2_3d, e3_3d);
    auto e31_constructed = wdg(e3_3d, e1_3d);
    auto e12_constructed = wdg(e1_3d, e2_3d);
    value_t g_e23 = value_t(dot(e23_constructed, e23_constructed));  // = 1
    value_t g_e31 = value_t(dot(e31_constructed, e31_constructed));  // = 1
    value_t g_e12 = value_t(dot(e12_constructed, e12_constructed));  // = 1
    CHECK(abs(g_e23 - value_t(G[4, 4])) < eps);
    CHECK(abs(g_e31 - value_t(G[5, 5])) < eps);
    CHECK(abs(g_e12 - value_t(G[6, 6])) < eps);

    // Level 3: Trivector (constructed recursively)
    auto e123_constructed = wdg(wdg(e1_3d, e2_3d), e3_3d);
    value_t g_e123 = value_t(dot(e123_constructed, e123_constructed));  // = 1
    CHECK(abs(g_e123 - value_t(G[7, 7])) < eps);
}
```

**Result**: ✓ All checks pass

#### PGA2DP - G(2,0,1)
**File**: [ga_pga2dp_test.hpp:3518-3562](../ga_test/src/ga_pga2dp_test.hpp#L3518-L3562)

**Challenge**: Null dimension e3² = 0 propagates through wedge products

```cpp
TEST_CASE("G<2,0,1>: extended metric recursive extraction via wedge products")
{
    auto G = pga2dp_metric_view();

    // Level 0: Scalar
    CHECK(G[0, 0] == 1);

    // Level 1: Vectors (e3 is null)
    value_t g_e1 = value_t(dot(e1_2dp, e1_2dp));  // = 1
    value_t g_e2 = value_t(dot(e2_2dp, e2_2dp));  // = 1
    value_t g_e3 = value_t(dot(e3_2dp, e3_2dp));  // = 0 (null dimension)
    CHECK(abs(g_e1 - value_t(G[1, 1])) < eps);
    CHECK(abs(g_e2 - value_t(G[2, 2])) < eps);
    CHECK(abs(g_e3 - value_t(G[3, 3])) < eps);

    // Level 2: Bivectors
    auto e23_constructed = wdg(e2_2dp, e3_2dp);  // involves e3 → metric = 0
    auto e31_constructed = wdg(e3_2dp, e1_2dp);  // involves e3 → metric = 0
    auto e12_constructed = wdg(e1_2dp, e2_2dp);  // Euclidean → metric = 1
    value_t g_e23 = value_t(dot(e23_constructed, e23_constructed));  // = 0
    value_t g_e31 = value_t(dot(e31_constructed, e31_constructed));  // = 0
    value_t g_e12 = value_t(dot(e12_constructed, e12_constructed));  // = 1
    CHECK(abs(g_e23 - value_t(G[4, 4])) < eps);
    CHECK(abs(g_e31 - value_t(G[5, 5])) < eps);
    CHECK(abs(g_e12 - value_t(G[6, 6])) < eps);

    // Level 3: Trivector (pseudoscalar)
    auto e321_constructed = wdg(wdg(e3_2dp, e2_2dp), e1_2dp);  // involves e3 → metric = 0
    value_t g_e321 = value_t(dot(e321_constructed, e321_constructed));  // = 0
    CHECK(abs(g_e321 - value_t(G[7, 7])) < eps);
}
```

**Result**: ✓ All checks pass (null dimensions handled correctly)

#### PGA3DP - G(3,0,1)
**File**: [ga_pga3dp_test.hpp:4263-4328](../ga_test/src/ga_pga3dp_test.hpp#L4263-L4328)

**Challenge**: More complex propagation of null dimension e4² = 0

```cpp
TEST_CASE("G<3,0,1>: extended metric recursive extraction via wedge products")
{
    auto G = pga3dp_metric_view();

    // Level 0: Scalar
    CHECK(G[0, 0] == 1);

    // Level 1: Vectors (e4 is null)
    value_t g_e1 = value_t(dot(e1_3dp, e1_3dp));  // = 1
    value_t g_e2 = value_t(dot(e2_3dp, e2_3dp));  // = 1
    value_t g_e3 = value_t(dot(e3_3dp, e3_3dp));  // = 1
    value_t g_e4 = value_t(dot(e4_3dp, e4_3dp));  // = 0 (null)
    CHECK(abs(g_e1 - value_t(G[1, 1])) < eps);
    CHECK(abs(g_e2 - value_t(G[2, 2])) < eps);
    CHECK(abs(g_e3 - value_t(G[3, 3])) < eps);
    CHECK(abs(g_e4 - value_t(G[4, 4])) < eps);

    // Level 2: Bivectors (6 total)
    auto e41_constructed = wdg(e4_3dp, e1_3dp);  // involves e4 → 0
    auto e42_constructed = wdg(e4_3dp, e2_3dp);  // involves e4 → 0
    auto e43_constructed = wdg(e4_3dp, e3_3dp);  // involves e4 → 0
    auto e23_constructed = wdg(e2_3dp, e3_3dp);  // Euclidean → 1
    auto e31_constructed = wdg(e3_3dp, e1_3dp);  // Euclidean → 1
    auto e12_constructed = wdg(e1_3dp, e2_3dp);  // Euclidean → 1
    value_t g_e41 = value_t(dot(e41_constructed, e41_constructed));  // = 0
    value_t g_e42 = value_t(dot(e42_constructed, e42_constructed));  // = 0
    value_t g_e43 = value_t(dot(e43_constructed, e43_constructed));  // = 0
    value_t g_e23 = value_t(dot(e23_constructed, e23_constructed));  // = 1
    value_t g_e31 = value_t(dot(e31_constructed, e31_constructed));  // = 1
    value_t g_e12 = value_t(dot(e12_constructed, e12_constructed));  // = 1
    CHECK(abs(g_e41 - value_t(G[5, 5])) < eps);
    CHECK(abs(g_e42 - value_t(G[6, 6])) < eps);
    CHECK(abs(g_e43 - value_t(G[7, 7])) < eps);
    CHECK(abs(g_e23 - value_t(G[8, 8])) < eps);
    CHECK(abs(g_e31 - value_t(G[9, 9])) < eps);
    CHECK(abs(g_e12 - value_t(G[10, 10])) < eps);

    // Level 3: Trivectors (4 total)
    auto e423_constructed = wdg(wdg(e4_3dp, e2_3dp), e3_3dp);  // involves e4 → 0
    auto e431_constructed = wdg(wdg(e4_3dp, e3_3dp), e1_3dp);  // involves e4 → 0
    auto e412_constructed = wdg(wdg(e4_3dp, e1_3dp), e2_3dp);  // involves e4 → 0
    auto e321_constructed = wdg(wdg(e3_3dp, e2_3dp), e1_3dp);  // Euclidean → 1
    value_t g_e423 = value_t(dot(e423_constructed, e423_constructed));  // = 0
    value_t g_e431 = value_t(dot(e431_constructed, e431_constructed));  // = 0
    value_t g_e412 = value_t(dot(e412_constructed, e412_constructed));  // = 0
    value_t g_e321 = value_t(dot(e321_constructed, e321_constructed));  // = 1
    CHECK(abs(g_e423 - value_t(G[11, 11])) < eps);
    CHECK(abs(g_e431 - value_t(G[12, 12])) < eps);
    CHECK(abs(g_e412 - value_t(G[13, 13])) < eps);
    CHECK(abs(g_e321 - value_t(G[14, 14])) < eps);

    // Level 4: Pseudoscalar
    auto e1234_constructed = wdg(wdg(wdg(e1_3dp, e2_3dp), e3_3dp), e4_3dp);  // involves e4 → 0
    value_t g_e1234 = value_t(dot(e1234_constructed, e1234_constructed));  // = 0
    CHECK(abs(g_e1234 - value_t(G[15, 15])) < eps);
}
```

**Result**: ✓ All checks pass

### Test Results

**Test Execution**:
```bash
cd build/ga_test && ./ga_ega_test && ./ga_pga_test && cd ../..
```

**Output**:
```
[doctest] test cases:   76 |   76 passed | 0 failed | 0 skipped
[doctest] assertions: 2323 | 2323 passed | 0 failed |
...
[doctest] test cases:   84 |   84 passed | 0 failed | 0 skipped
[doctest] assertions: 2055 | 2055 passed | 0 failed |

Total: 4378 assertions passed ✓
```

### Key Findings

1. **Recursive extraction works perfectly** for all implemented algebras (EGA2D, EGA3D, PGA2DP, PGA3DP)

2. **Null dimensions are handled correctly**: The approach naturally propagates degenerate metrics through wedge products

3. **100% validation**: Every extracted extended metric value matches the pre-computed values in the algebra's metric matrix

4. **Implementation-agnostic**: The approach works regardless of how the algebra is implemented internally, since it only relies on the algebra's own wedge and dot product operations

### STA4D Status

**Cannot test STA4D yet**: The STA4D wedge and dot products are **not yet implemented** in the `ga/` library. They exist only in the `ga_prdxpr` code generator as text string rules, not as executable operations.

**File**: [ga_sta4d_test.hpp](../ga_test/src/ga_sta4d_test.hpp) contains only 4 basic test cases:
- Type checking test
- Dual composition properties test
- Complement composition properties test
- Comparison test

No product operations (wedge, dot, geometric) are available to construct basis elements for testing.

### Conclusion

**Validation Complete**: The recursive extraction approach is mathematically sound and has been proven to work for:
- ✓ Euclidean algebras (EGA2D, EGA3D)
- ✓ Projective algebras with null dimensions (PGA2DP, PGA3DP)
- ✓ All grade levels (0 through pseudoscalar)
- ✓ All 4378 test assertions pass

**Recommendation**: This approach should be adopted in `ga_prdxpr` for extended metric calculation, replacing the Gram determinant approach that fails for STA4D.

**Next Steps**:
1. Commit these validation tests as regression tests
2. Analyze why Gram determinants fail for STA4D
3. Implement recursive extraction in `ga_prdxpr_rule_generator.cpp` or `ga_prdxpr_rule_generator_test.cpp`

**Related Code**:
- [ga_ega2d_test.hpp:2503-2528](../ga_test/src/ga_ega2d_test.hpp#L2503-L2528)
- [ga_ega3d_test.hpp:2964-3007](../ga_test/src/ga_ega3d_test.hpp#L2964-L3007)
- [ga_pga2dp_test.hpp:3518-3562](../ga_test/src/ga_pga2dp_test.hpp#L3518-L3562)
- [ga_pga3dp_test.hpp:4263-4328](../ga_test/src/ga_pga3dp_test.hpp#L4263-L4328)
