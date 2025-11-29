# Implementation Plan: Matrix-Based Extended Metric and Systematic Dual Generation

## Executive Summary

This implementation plan describes a systematic approach to:

1. **Compute the full extended metric matrix** (including off-diagonal elements) for all algebras
2. **Calculate the regressive metric matrix** using proper matrix transformations
3. **Generate dual rules systematically** using matrix-vector multiplication with the extended metric

**Key Goals:**

- ✅ 100% character-identical output to current implementation for all existing algebras (EGA2D, EGA3D, PGA2DP, PGA3DP, STA4D)
- ✅ Foundation for future CGA support (handles off-diagonal metric elements)
- ✅ Unified mathematical framework based on extended metric theory
- ✅ No performance regression (computation done once at initialization)
- ✅ Easy integration of calculated metrics into `ga_usr_consts.hpp`

---

## Current Implementation Issues

### Problem 1: Diagonal-Only Assumption

Current code at [ga_prdxpr_rule_generator.cpp:265-287](../ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp#L265-L287):

```cpp
// LIMITATION: Off-diagonal elements remain 0 (orthogonal basis assumption)
```

**Impact**: Cannot support CGA or any non-orthogonal algebra.

### Problem 2: Dual Calculation Uses Only Diagonal

Current code at [ga_prdxpr_rule_generator.cpp:733](../ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp#L733):

```cpp
int metric_value = G[i, i];  // Only uses diagonal!
```

**Impact**: Even if we compute full matrix, dual generation ignores off-diagonal elements.

### Problem 3: Regressive Metric Calculation

Current code at [ga_prdxpr_rule_generator.cpp:454](../ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp#L454):

```cpp
G_reg[i, i] = (G[i, i] == 0) ? 1 : 0;  // Simple inversion
```

**Impact**: Doesn't work for non-diagonal matrices.

---

## Mathematical Foundation

### Extended Metric Definition

For basis blades **B_i** and **B_j** of the same grade k:

```text
G[i,j] = ⟨B_i, B_j⟩ = det(Gram_Matrix(B_i, B_j))
```

Where the Gram matrix for blades **A** = **a₁** ∧ ... ∧ **aₖ** and **B** = **b₁** ∧ ... ∧ **bₖ** is:

```text
Gram[p,q] = aₚ · bᵩ   (using vector metric)
```

**Key Properties:**

1. **Block-diagonal by grade** (different grades → zero metric)
2. **Symmetric**: G[i,j] = G[j,i]
3. **Diagonal for orthogonal bases** (EGA, PGA, STA)
4. **Non-diagonal for null bases** (CGA with n₀, n_∞)

### Dual Operation

The dual operation is defined as:

```text
dual(u) = complement(G × u)
```

Where:

- **G** = extended metric matrix
- **×** = matrix-vector multiplication
- **complement** = wedge-based complement transformation

**Current issue**: Implementation only uses `G[i,i]` diagonal element!

### Regressive Metric

For projective algebras, the regressive metric **Ḡ** is defined via complement:

```text
Ḡ = C · G · C^T
```

Where **C** is the complement transformation matrix.

**Current issue**: Implementation uses `Ḡ[i,i] = 1 - G[i,i]` which only works for diagonal PGA.

---

## Implementation Roadmap

### Step 1: Add Full Metric Computation ✅ COMPLETED

**Files to modify:**

- `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp`
- `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.hpp`

**Tasks:**

1.1. **Add integer determinant computation** (30 min)

- Location: Add to `ga_prdxpr_rule_generator.cpp` before extended metric functions
- Function: `compute_integer_determinant_2x2()`, `compute_integer_determinant_3x3()`, `compute_integer_determinant_4x4()`, `compute_integer_determinant_5x5()`
- Reason: Support up to 5×5 for CGA3D bivector metrics
- Implementation: Direct formulas for 2×2, 3×3; cofactor expansion for 4×4, 5×5

1.2. **Add Gram matrix determinant computation** (30 min)

- Function: `compute_gram_determinant(indices_a, indices_b, config)`
- Build k×k Gram matrix where `Gram[p,q] = vector_metric(indices_a[p], indices_b[q])`
- Call appropriate determinant function based on k

1.3. **Add vector metric lookup** (20 min)

- Function: `get_vector_metric(index_i, index_j, config)`
- Handle 1-based vs 0-based indexing
- For diagonal: use `config.metric_signature[i]` when `i == j`
- For off-diagonal: use `config.metric_matrix[i][j]` if available (future CGA)
- Return 0 for orthogonal bases when `i ≠ j`

1.4. **Add blade inner product computation** (30 min)

- Function: `compute_blade_inner_product(blade_a, blade_b, config)`
- Handle special cases: scalars, different grades
- For grade 1: use `get_vector_metric()` directly
- For higher grades: use `compute_gram_determinant()`

1.5. **Update extended metric matrix calculation** (30 min)

- Keep existing `calculate_extended_metric_matrix()` function signature
- Replace diagonal-only logic with full matrix computation:

  ```cpp
  for (size_t i = 0; i < n; ++i) {
      for (size_t j = i; j < n; ++j) {  // Upper triangular
          int metric_value = compute_blade_inner_product(basis[i], basis[j], config);
          G[i, j] = metric_value;
          G[j, i] = metric_value;  // Symmetric
      }
  }
  ```

**Validation for Step 1:**

- Verify diagonal elements `G[i,i]` match current implementation
- Verify off-diagonal elements `G[i,j]` are zero for orthogonal bases (EGA, PGA, STA)
- Print matrix for each algebra to inspect structure
- Test: Create simple test that compares old diagonal with new matrix diagonal

**Deliverable:**

- Extended metric matrix with ALL elements computed correctly
- Diagonal elements identical to current implementation
- Off-diagonal elements zero for orthogonal bases (proving backward compatibility)

---

### Step 2: Update Dual Generation ✅ COMPLETED

**Files to modify:**

- `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp`

**Tasks:**

2.1. **Add helper functions** (20 min)

- `format_signed_element(coeff, elem)` - Format "coeff × elem" with proper signs
- `apply_complement(signed_elem, complement_rules)` - Apply complement with sign handling
- `parse_signed_basis_element(signed_elem, basis)` - Parse element and extract sign

2.2. **Update dual generation functions** (40 min)

- Modify `generate_left_dual_rules()`:
  - Change from `int metric_value = G[i, i]` to full row lookup
  - For each basis element i: compute `Σ G[i,j] × basis[j]` for all j
  - For diagonal matrices (current): sum has 1 term → same result as before
  - Apply complement to result

- Modify `generate_right_dual_rules()` similarly
- Modify `generate_dual_rules()` similarly

**Implementation pattern:**

```cpp
std::vector<std::pair<std::string, int>> weighted_sum;
for (size_t j = 0; j < n; ++j) {
    int g_ij = G[i, j];
    if (g_ij != 0) {
        weighted_sum.push_back({basis[j], g_ij});
    }
}

// For diagonal: weighted_sum.size() == 1 (at most)
// Apply complement to the single term (current behavior)
```

**Validation for Step 2:**

- Run all existing tests
- Character-by-character comparison of generated dual rules
- Verify: EGA2D, EGA3D, PGA2DP, PGA3DP, STA4D all produce identical output
- Test command: `diff old_dual_rules.txt new_dual_rules.txt` → no differences

**Deliverable:**

- Dual rules generated via full matrix-vector multiplication
- 100% character-identical output for all existing algebras

---

### Step 3: Implement Regressive Metric ✅ COMPLETED (as bonus step before Step 2)

**Files to modify:**

- `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp`

**Tasks:**

3.1. **Add matrix multiplication utilities** (30 min)

- `matrix_multiply(A, B, n)` - Standard matrix multiplication
- `matrix_transpose(M, n)` - Matrix transpose
- `matrix_triple_product(C, G, n)` - Compute C · G · C^T efficiently

3.2. **Build complement transformation matrix** (45 min)

- Function: `build_complement_matrix(config, complement_rules)`
- Create n×n matrix C where `C[i,j]` = coefficient of `basis[j]` in `complement(basis[i])`
- Parse complement rules to extract basis element and sign
- For orthogonal bases: C is permutation matrix with ±1 entries

3.3. **Update regressive metric calculation** (30 min)

- Keep existing function signature: `calculate_regressive_extended_metric_matrix(config)`
- Add new implementation path:
  - Generate complement rules
  - Build complement matrix C
  - Compute `Ḡ = C · G · C^T`
- For diagonal PGA: verify result matches `Ḡ[i,i] = 1 - G[i,i]`

3.4. **Validation** (15 min)

- Verify regressive metric is symmetric
- For PGA: verify diagonal matches simple inversion rule
- Test regressive product rules remain identical

**Validation for Step 3:**

- PGA2DP and PGA3DP regressive products unchanged
- Regressive metric matrix properties:
  - Symmetric: `Ḡ[i,j] = Ḡ[j,i]`
  - For PGA diagonal: `Ḡ[i,i] = 1 - G[i,i]`
  - Block-diagonal by grade

**Deliverable:**

- Regressive metric computed via proper matrix transformation
- 100% character-identical regressive product rules for PGA algebras

---

### Step 4: Integration with ga_usr_consts.hpp ✅ COMPLETED

**Goal**: Make calculated extended and regressive metrics easily exportable to `ga_usr_consts.hpp`

**Files to modify:**

- `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp`
- `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.hpp`

**Tasks:**

4.1. **Add metric export functionality** (30 min)

- Function: `export_extended_metric_to_code(algebra_name, basis, matrix)`
- Generates C++ code for constexpr arrays
- Output format:

  ```cpp
  // Extended metric for EGA3D
  inline constexpr std::array<int, 64> ega3d_extended_metric = {
      1, 0, 0, 0, ...  // Flattened row-major 8×8 matrix
  };

  // Helper: Access as G[i,j] = ega3d_extended_metric[i*8 + j]
  ```

4.2. **Add regressive metric export** (20 min)

- Function: `export_regressive_metric_to_code(algebra_name, basis, matrix)`
- Similar format for PGA regressive metrics
- Output format:

  ```cpp
  // Regressive extended metric for PGA3DP
  inline constexpr std::array<int, 256> pga3dp_regressive_metric = {
      0, 0, 0, 0, ...  // Flattened row-major 16×16 matrix
  };
  ```

4.3. **Add print/export command** (10 min)

- Add option to print metrics in format ready for copy-paste into `ga_usr_consts.hpp`
- Include basis element names as comments for clarity:

  ```cpp
  // Basis: 1, e1, e2, e3, e23, e31, e12, e123
  inline constexpr std::array<int, 64> ega3d_extended_metric = {
      // Row 0 (1):
      1, 0, 0, 0, 0, 0, 0, 0,
      // Row 1 (e1):
      0, 1, 0, 0, 0, 0, 0, 0,
      // ...
  };
  ```

**Validation for Step 4:**

- Copy exported metric into test file
- Verify constexpr compilation
- Verify values match computed matrix
- Test runtime access: `G[i*n + j]` gives correct value

**Deliverable:**

- Functions to export metrics as C++ constexpr arrays
- Format ready for integration into `ga_usr_consts.hpp`
- Documentation on how to use exported metrics

---

### Step 5: Backward Compatibility & Feature Flag (30 min)

**Files to modify:**

- `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.hpp`
- `ga_prdxpr/src_prdxpr/ga_prdxpr_rule_generator.cpp`

**Tasks:**

5.1. **Add feature flag** (10 min)

- Add compile-time option: `#define GA_PRDXPR_USE_FULL_METRIC_MATRIX 1`
- Allow toggling between old (diagonal-only) and new (full matrix) implementations
- Default: NEW implementation (full matrix)

5.2. **Conditional compilation** (10 min)

- Wrap old implementation in `#if !GA_PRDXPR_USE_FULL_METRIC_MATRIX`
- Keep old code for reference and emergency fallback
- Eventually remove after confidence in new implementation

5.3. **Add regression test** (10 min)

- Test that compares old vs new implementation output
- Automated verification of character-identical results
- Run as part of CI/test suite

**Deliverable:**

- Smooth migration path with fallback option
- Automated regression testing

---

### Step 6: Documentation (1 hour)

**Files to create/modify:**

- Update `CLAUDE.md` with new metric computation approach
- Update comments in `ga_prdxpr_rule_generator.cpp`

**Tasks:**

6.1. **Update CLAUDE.md** (30 min)

- Document new matrix-based approach
- Explain mathematical foundation
- Add section on extended metric computation
- Document export functionality for `ga_usr_consts.hpp`

6.2. **Code comments** (20 min)

- Add detailed comments to new functions
- Explain Gram matrix determinant calculation
- Document matrix-vector multiplication in dual generation
- Add references to mathematical sources

6.3. **Usage examples** (10 min)

- Show how to export metrics
- Example of integrating into `ga_usr_consts.hpp`
- Example of accessing matrix elements at runtime

**Deliverable:**

- Comprehensive documentation of new approach
- Clear usage examples
- Mathematical references

---

## Total Effort Estimate

| Step | Description | Time |
|------|-------------|------|
| 1 | Full metric computation | 2h |
| 2 | Update dual generation | 1h |
| 3 | Regressive metric | 2h |
| 4 | Integration with ga_usr_consts.hpp | 1h |
| 5 | Backward compatibility | 0.5h |
| 6 | Documentation | 1h |
| **Total** | | **7.5h** |

---

## Validation Strategy

### Requirement: 100% Character-Identical Output

**Critical**: All existing algebras must produce EXACTLY the same results.

**Mathematical Reason**: For orthogonal bases (EGA, PGA, STA):

- Extended metric is diagonal: `G[i,j] = 0` for `i ≠ j`
- Matrix-vector multiplication: `(G × e_i) = G[i,i] × e_i`
- This is IDENTICAL to current implementation!

**Test Plan:**

1. **Run existing tests**:

   ```bash
   cd build
   cmake --build . --target ga_prdxpr_rule_generator_test
   ./ga_prdxpr/ga_prdxpr_rule_generator_test
   ```

2. **Character-by-character comparison**:
   - Generate rules with old implementation → save to file
   - Generate rules with new implementation → save to file
   - Compare: `diff old_rules.txt new_rules.txt`
   - **Expected**: No differences!

3. **Verify for each algebra**:
   - EGA2D: All 4×4 = 16 basis elements
   - EGA3D: All 8×8 = 64 basis elements
   - PGA2DP: All 8×8 = 64 basis elements
   - PGA3DP: All 16×16 = 256 basis elements
   - STA4D: All 16×16 = 256 basis elements

### Mathematical Verification

**Property 1**: Extended metric is block-diagonal by grade

```cpp
for (i, j) where grade(basis[i]) ≠ grade(basis[j]):
    assert(G[i,j] == 0)
```

**Property 2**: Symmetry

```cpp
for all (i, j):
    assert(G[i,j] == G[j,i])
```

**Property 3**: Dual identity for EGA

```cpp
for all basis elements u:
    assert(dual(dual(u)) == ±u)  // Double dual = ± identity
```

**Property 4**: Regressive-Geometric duality for PGA

```cpp
for all u, v:
    assert(rgpr(u,v) == dual(gpr(dual(u), dual(v))))
```

---

## Benefits

### Immediate Benefits

1. ✅ **Mathematical rigor**: Proper implementation of extended metric theory
2. ✅ **Unified framework**: Single approach for all algebras
3. ✅ **Better documentation**: Clear connection to mathematical definitions
4. ✅ **Bug prevention**: Systematic computation reduces ad-hoc special cases
5. ✅ **Integration ready**: Easy export to `ga_usr_consts.hpp`

### Future Benefits

1. ✅ **CGA support**: Foundation ready for conformal geometric algebra
2. ✅ **Other non-orthogonal algebras**: Motor algebra, etc.
3. ✅ **Numerical approach option**: Can compute metric numerically if needed
4. ✅ **Symbolic conversion**: Matrix can be converted back to symbolic form
5. ✅ **Runtime access**: Metrics available as constexpr arrays

### No Downsides

- ✅ **Performance**: Same O(1) lookup after initial computation
- ✅ **Compatibility**: 100% backward compatible
- ✅ **Code size**: Minimal increase (better organized)

---

## Notes for Implementation

### Integer Arithmetic

- All metrics use integer values (no floating-point)
- Determinants computed with integer arithmetic
- Avoids numerical precision issues

### Basis Indexing

- Auto-detect 1-based vs 0-based indexing
- Handle both conventions transparently
- Same logic as current implementation

### Matrix Storage

- Row-major flattened storage: `matrix[i*n + j]`
- Use `std::mdspan` for convenient 2D access
- Compatible with constexpr arrays for export

### Export Format

- Generate valid C++ code ready for copy-paste
- Include basis names as comments
- Provide helper functions for matrix access
- Document indexing convention clearly

---

## References

### Code Locations

- **AlgebraConfig**: `ga_prdxpr_rule_generator.hpp:85`
- **Extended metric (current)**: `ga_prdxpr_rule_generator.cpp:168`
- **Dual generation (current)**: `ga_prdxpr_rule_generator.cpp:711`
- **Regressive metric (current)**: `ga_prdxpr_rule_generator.cpp:437`

### Mathematical References

- Extended metric definition: See `TODO/metric.md`
- CGA null basis structure: See `TODO/cga_mod.md`
- Gram matrix determinant: Standard linear algebra

### Related Documentation

- `CLAUDE.md`: Project architecture and GA fundamentals
- `TODO/metric.md`: Extended metric computation theory
- `TODO/cga_mod.md`: CGA requirements and matrix structure

---

## Success Criteria

✅ All existing tests pass
✅ Character-identical output for all algebras
✅ Extended metric matrix properly computed (diagonal + off-diagonal)
✅ Dual rules generated via matrix-vector multiplication
✅ Regressive metric computed via complement transformation
✅ Metrics exportable to `ga_usr_consts.hpp` format
✅ Documentation updated
✅ No performance regression

---

**Status**: Ready to implement
**Start with**: Step 1 - Full Metric Computation
**Next**: Proceed sequentially through steps 2-6
