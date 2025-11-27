# CGA Support: Matrix-Based Extended Metric Implementation

**Status**: Documentation only (implementation not started)
**Effort Estimate**: 3-4 hours of focused work
**Priority**: Medium (enables future CGA support)
**Risk**: Medium (core changes to metric system, extensive testing required)

## Executive Summary

Current implementation assumes **orthogonal basis** (diagonal metric). This works for:

- ✅ EGA (Euclidean GA) - orthonormal basis
- ✅ PGA (Projective GA) - orthogonal with degenerate directions
- ✅ STA (Space-Time Algebra) - orthogonal Minkowski signature

**Not supported**: CGA (Conformal GA) - has non-orthogonal null vectors with off-diagonal metric elements.

This document describes how to extend the system to support full metric matrices, enabling CGA and other non-orthogonal algebras.

## Current Architecture

### Key Data Structure

```cpp
// File: ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.hpp (line 85)
struct AlgebraConfig {
    std::vector<std::string> basis_vectors;     // {"e1", "e2", "e3"}
    std::vector<int> metric_signature;          // {+1, +1, +1} DIAGONAL ONLY!
    mvec_coeff multivector_basis;               // {"1", "e1", "e2", "e12", ...}
    std::string scalar_name = one_str();
    std::string basis_prefix = "e";
};
```

### Limitation

`metric_signature` only stores **diagonal elements** g_{ii}. For CGA, we need off-diagonal elements g_{ij} where i≠j.

**Example CGA3D metric** (what we need to support):

```text
        e1   e2   e3   e+   e-
e1   [  1    0    0    0    0  ]
e2   [  0    1    0    0    0  ]
e3   [  0    0    1    0    0  ]
e+   [  0    0    0    0   -1  ]  ← e+·e- = -1 (OFF-DIAGONAL!)
e-   [  0    0    0   -1    0  ]  ← e-·e+ = -1 (symmetric)
```

## Implementation Plan

### Phase 1: Extend AlgebraConfig (15 minutes)

**File**: `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.hpp`

**Option 1 - Backward Compatible** (recommended):

```cpp
struct AlgebraConfig {
    std::vector<std::string> basis_vectors;
    std::vector<int> metric_signature;  // Keep for existing algebras

    // NEW: Optional full metric matrix for non-orthogonal algebras
    std::optional<std::vector<std::vector<int>>> metric_matrix;  // n×n matrix

    mvec_coeff multivector_basis;
    std::string scalar_name = one_str();
    std::string basis_prefix = "e";

    // Helper: Check if metric matrix is provided
    bool has_full_metric() const { return metric_matrix.has_value(); }
};
```

**Option 2 - Breaking Change** (cleaner, but requires config updates):

```cpp
struct AlgebraConfig {
    std::vector<std::string> basis_vectors;
    std::vector<std::vector<int>> metric_matrix;  // Replace metric_signature
    // ... rest unchanged
};
```

**Recommendation**: Use Option 1 for backward compatibility.

### Phase 2: Rewrite Extended Metric Matrix Calculation (45 minutes)

**File**: `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp` (line 168)

**Current implementation** (diagonal only):

```cpp
std::vector<int> calculate_extended_metric_matrix(const AlgebraConfig& config)
{
    // ... initializes diagonal elements only
    // Off-diagonal elements remain 0 (line 265-287)
}
```

**New implementation** (computes all elements):

```cpp
std::vector<int> calculate_extended_metric_matrix(const AlgebraConfig& config)
{
    const auto& basis = config.multivector_basis;
    const size_t n = basis.size();

    // Flattened row-major storage
    std::vector<int> matrix_data(n * n, 0);
    std::mdspan G{matrix_data.data(), n, n};

    // Compute ALL elements G[i,j] = basis[i] · basis[j]
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {  // Upper triangular (symmetric)
            // Use existing dot product logic
            auto [result_element, sign] =
                compute_basis_dot_product(basis[i], basis[j], config);

            if (result_element == config.scalar_name) {
                G[i, j] = sign;
                G[j, i] = sign;  // Symmetric matrix
            } else {
                G[i, j] = 0;
                G[j, i] = 0;
            }
        }
    }

    return matrix_data;
}

// Helper function to compute basis element dot products
std::pair<std::string, int> compute_basis_dot_product(
    const std::string& a, const std::string& b, const AlgebraConfig& config)
{
    // Handle scalar cases
    if (a == config.scalar_name && b == config.scalar_name) {
        return {config.scalar_name, 1};
    }
    if (a == config.scalar_name || b == config.scalar_name) {
        return {zero_str(), 0};
    }

    // Use geometric product to compute: a·b = (ab + ba)/2
    // For basis elements: a·b is the grade-0 part of geometric product
    auto [product, sign] = multiply_basis_elements(a, b, config);

    // If product is scalar, this is the dot product
    if (product == config.scalar_name) {
        return {config.scalar_name, sign};
    }

    // If product is not scalar, dot product is zero (different grades)
    return {zero_str(), 0};
}
```

**Key insight**: The existing `multiply_basis_elements()` function already does most of the work!

### Phase 3: Update Dot Product Generation (30 minutes)

**File**: `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp` (line 536)

**Current limitation** (line 549):

```cpp
// Assumes diagonal metric - returns 0 for different basis elements
if (a != b) {
    return {zero_str(), 0};
}
```

**New implementation**:

```cpp
std::pair<std::string, int> multiply_basis_elements_dot(
    const std::string& a, const std::string& b, const AlgebraConfig& config)
{
    // Handle scalar cases
    if (a == config.scalar_name && b == config.scalar_name) {
        return {config.scalar_name, 1};
    }
    if (a == config.scalar_name || b == config.scalar_name) {
        return {zero_str(), 0};
    }

    // Get extended metric matrix (handles both diagonal and non-diagonal)
    auto matrix_data = calculate_extended_metric_matrix(config);
    size_t n = config.multivector_basis.size();
    std::mdspan G{matrix_data.data(), n, n};

    // Find indices of a and b in basis
    auto it_a = std::find(config.multivector_basis.begin(),
                          config.multivector_basis.end(), a);
    auto it_b = std::find(config.multivector_basis.begin(),
                          config.multivector_basis.end(), b);

    if (it_a == config.multivector_basis.end() ||
        it_b == config.multivector_basis.end()) {
        return {zero_str(), 0};
    }

    size_t idx_a = std::distance(config.multivector_basis.begin(), it_a);
    size_t idx_b = std::distance(config.multivector_basis.begin(), it_b);

    // Use matrix element (now handles off-diagonal!)
    int metric_value = G[idx_a, idx_b];

    return (metric_value == 0)
        ? std::make_pair(zero_str(), 0)
        : std::make_pair(config.scalar_name, metric_value);
}
```

### Phase 4: Update Geometric Product (30-45 minutes)

**File**: `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp` (line 89)

**Current code** (uses diagonal metric):

```cpp
sign *= config.metric_signature[metric_index];
```

**Updated code** (supports full matrix):

```cpp
// Apply metric when identical basis vectors appear: e_i * e_i
int metric_value;
if (config.has_full_metric()) {
    // Use full metric matrix diagonal element
    metric_value = (*config.metric_matrix)[metric_index][metric_index];
} else {
    // Fallback to diagonal-only metric (backward compatible)
    metric_value = config.metric_signature[metric_index];
}
sign *= metric_value;
```

**Important note**: Geometric product only needs **diagonal** of the basis vector metric. The off-diagonal elements only appear in the extended metric for higher-grade basis elements.

### Phase 5: Create CGA Configuration Example (20 minutes)

**New file**: `ga_prdxpr/src_prdxpr/ga_prdxpr_cga3d_config.cpp`

```cpp
#include "ga_prdxpr_cga3d.hpp"
#include "ga_prdxpr_rule_generator.hpp"

AlgebraConfig get_cga3d_algebra_config() {
    AlgebraConfig config;

    // Basis vectors: 3 Euclidean + 2 null vectors (e+, e-)
    config.basis_vectors = {"e1", "e2", "e3", "e+", "e-"};

    // Diagonal metric signature (for geometric product)
    config.metric_signature = {+1, +1, +1, 0, 0};

    // Full metric matrix (for extended metric)
    config.metric_matrix = {
        {1, 0, 0, 0, 0},   // e1·e1=1, others=0
        {0, 1, 0, 0, 0},   // e2·e2=1
        {0, 0, 1, 0, 0},   // e3·e3=1
        {0, 0, 0, 0, -1},  // e+·e+=0, e+·e-=-1 (OFF-DIAGONAL!)
        {0, 0, 0, -1, 0},  // e-·e-=0, e-·e+=-1 (symmetric)
    };

    // Complete basis (32 elements for 5D algebra)
    config.multivector_basis = mvCGA3d_basis;
    config.scalar_name = one_str();
    config.basis_prefix = "e";

    return config;
}

// Generate all rules
static auto cga3d_generated_rules = generate_algebra_rules(get_cga3d_algebra_config());

const prd_rules gpr_cga3d_rules = cga3d_generated_rules.geometric_product;
const prd_rules wdg_cga3d_rules = cga3d_generated_rules.wedge_product;
const prd_rules dot_cga3d_rules = cga3d_generated_rules.dot_product;
// ... etc
```

### Phase 6: Update Regressive Metric for CGA (Optional, Advanced)

**Current limitation**: Regressive metric uses simple inversion `Ḡ[i,i] = 1 - G[i,i]` which only works for diagonal matrices.

**For CGA**: Would need proper matrix transformation.

**File**: `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp` (line 418)

**Approach**: Implement via complement transformation:

```cpp
// Ḡ = C · G · C^T
// where C is the complement operation as a matrix
```

This is **advanced** and not strictly necessary for basic CGA support.

## Testing Strategy

### 1. Regression Testing (60 minutes)

Run all existing tests to ensure backward compatibility:

```bash
cd build
cmake --build . --target ga_prdxpr_rule_generator_test
./ga_prdxpr/ga_prdxpr_rule_generator_test
```

**Expected**: All existing algebras (EGA2D, EGA3D, PGA2DP, PGA3DP, STA3D) should produce **identical** results.

### 2. CGA Validation

Create test cases for CGA3D:

- Verify null vector properties: e+·e+ = 0, e-·e- = 0
- Verify non-diagonal: e+·e- = -1
- Test geometric product: e+ * e- = (e+ ∧ e-) + e+·e-
- Validate extended metric matrix has correct off-diagonal elements

### 3. Performance Testing

Ensure no significant performance regression:

- Time existing algebra generation
- Compare before/after

## Files to Modify

### Core System

1. **ga_prdxpr_rule_generator.hpp** (line 85)
   - Extend `AlgebraConfig` structure

2. **ga_prdxpr_rule_generator.cpp**
   - Line 168: `calculate_extended_metric_matrix()` - full matrix computation
   - Line 536: `multiply_basis_elements_dot()` - use matrix lookup
   - Line 89: `multiply_basis_elements()` - support full metric

### New Files (for CGA)

1. **ga_prdxpr_cga3d.hpp** (new)
   - CGA3D basis definitions and rule declarations

2. **ga_prdxpr_cga3d_config.cpp** (new)
   - CGA3D algebra configuration with full metric matrix

### Build System

1. **CMakeLists.txt**
   - Add CGA files to build

## Performance Analysis

**Current**: O(1) metric lookup via direct array access
**New**: O(1) metric lookup via 2D array access `G[i,j]`

**Impact**: Negligible (2D array access is still O(1))

**Caching**: Extended metric matrix is computed once during initialization, then reused.

## Migration Path

### Recommended Approach: Incremental

1. ✅ Add `metric_matrix` as **optional** field (backward compatible)
2. ✅ Update extended metric computation to calculate full matrix
3. ✅ Modify dot product to use matrix instead of diagonal check
4. ✅ Test with existing algebras (should be identical)
5. ✅ Add CGA as proof of concept
6. ⚠️ Deprecate diagonal-only approach later (if desired)

### Alternative: Breaking Change

1. Replace `metric_signature` with `metric_matrix` entirely
2. Update all existing configs (EGA, PGA, STA) to use matrix format
3. Remove all diagonal-only code paths

**Not recommended** - too risky for minimal benefit.

## Known Issues & Future Work

### Issue 1: Regressive Metric for CGA

Current regressive metric calculation assumes diagonal matrices. For CGA:

- Cannot use `Ḡ[i,i] = 1 - G[i,i]`
- Need matrix-based complement transformation
- **Workaround**: Only compute regressive metric for PGA (current behavior)

### Issue 2: Dual Rules for CGA

Dual calculations may need updates for non-diagonal metrics:

- Current: `dual(u) = complement(G × u)` assumes diagonal G
- May need full matrix multiplication

### Issue 3: Special Cases

Some algebras (like STA3D) have algebra-specific extended metric rules (line 247). Need to ensure these still work with matrix approach.

## References

### Key Code Locations

- **AlgebraConfig**: `ga_prdxpr_rule_generator.hpp:85`
- **Extended metric (diagonal only)**: `ga_prdxpr_rule_generator.cpp:168`
- **Dot product (diagonal assumption)**: `ga_prdxpr_rule_generator.cpp:536-575`
- **Geometric product (metric application)**: `ga_prdxpr_rule_generator.cpp:89`
- **Regressive metric**: `ga_prdxpr_rule_generator.cpp:418`

### Documentation

- **CGA limitations documented**:
  - `ga_prdxpr_rule_generator.hpp:62-71` (AlgebraConfig comments)
  - `ga_prdxpr_rule_generator.cpp:265-287` (extended metric comments)
  - `ga_prdxpr_rule_generator.cpp:395-400` (regressive metric comments)

### Mathematical Background

- Extended metric: G(u) where G[i,j] = e_i · e_j
- CGA null vectors: e+ = e₄ + e₅, e- = (e₄ - e₅)/2
- CGA metric: e+·e+ = 0, e-·e- = 0, e+·e- = -1

## Conclusion

**Effort**: 3-4 hours
**Benefit**: Enables CGA and all future non-orthogonal algebras
**Risk**: Medium (requires careful testing)
**Recommendation**: ✅ **Do it** - solid investment for future extensibility

The matrix-based approach is a **natural extension** that maintains backward compatibility while enabling powerful new capabilities.
