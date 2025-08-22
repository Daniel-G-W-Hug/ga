# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Paths and Build System

**CRITICAL PATH INFORMATION:**

- **Project Root**: `/Users/hud3bh/prg/cpp/pj/ga/` (absolute path)
- **Build Directory**: `/Users/hud3bh/prg/cpp/pj/ga/build/` (absolute path)
- **Working Directory**: Always work from the build directory when running executables or when creating temporary files
- **Source Files**: Always in `/Users/hud3bh/prg/cpp/pj/ga/[module]/src/` or `/Users/hud3bh/prg/cpp/pj/ga/ga_prdxpr/src_trafo/`
- **File Organization**: Keep all temporary/debug files in the build directory (`build/`) rather than the source tree to maintain clean project organization.

**Build Directory Structure:**

```text
/Users/hud3bh/prg/cpp/pj/ga/build/          # Main build directory (working directory)
├── ga_lua/ga_lua                           # Lua interface executable  
├── ga_test/ga_ega_test                     # EGA test executable
├── ga_test/ga_pga_test                     # PGA test executable
├── ga_view/ga_view                         # Qt6 visualization executable
└── ga_prdxpr/ga_prdxpr                     # Main code generator
```

**File Path Rules for Code:**

- When reading source files from executables: `../[module]/src/filename` (relative to build dir)
- NEVER guess paths - use these established patterns
- All executables run from `/Users/hud3bh/prg/cpp/pj/ga/build/` directory

This project uses CMake (minimum version 3.28) with C++23 standard.

**IMPORTANT: Always use out-of-source builds exclusively in the `build/` directory.**

```bash
# Configure (from project root)
mkdir -p build
cd build
cmake ..

# Build (from build directory)  
cmake --build .

# Clean build (if needed)
cd .. && rm -rf build && mkdir build && cd build && cmake ..
```

**Critical Build Requirements:**

- ALL build artifacts must be generated in `build/` directory only
- NEVER run cmake directly in source directories (ga_lua/, ga_test/, etc.)
- NEVER run make/ninja commands outside of `build/` directory
- Build directory structure: `build/ga_lua/`, `build/ga_test/`, `build/ga_view/`, etc.
- Always run builds from the project root's `build/` directory

**Build Workflow:**

1. Always start from project root directory (`/Users/hud3bh/prg/cpp/pj/ga/`)
2. Ensure `build/` directory exists: `mkdir -p build`
3. Navigate to build directory: `cd build`
4. Configure: `cmake ..`
5. Build: `cmake --build .`
6. Run executables: `cd ga_lua && ./ga_lua && cd ..`, `cd ga_test && ./ga_ega_test &&
   cd..`, etc.
7. For refactoring with reference output data always match reference output 100% character-identical

## Running Tests

Execute test binaries from the build directory (ensure you're in `build/` directory):

```bash
cd ga_test && ./ga_ega_test && cd ..    # Tests for Euclidean GA
cd ga_test && ./ga_pga_test && cd ..    # Tests for Projective GA
```

## Running Applications

Execute applications from the build directory:

```bash
cd ga_lua && ./ga_lua && cd ..                    # Interactive Lua shell
cd ga_lua && ./ga_lua script.lua && cd ..         # Run Lua script
cd ga_view && ./ga_view && cd ..                  # Qt6 visualization tool
cd ga_prdxpr && ./ga_prdxpr && cd ..              # Code generator
```

## Project Architecture

This is a header-only geometric algebra library (`ga/`) with accompanying modules:

- **ga/**: Core header-only GA library supporting both Euclidean (EGA) and Projective (PGA) geometric algebra
- **ga_test/**: Test suite using doctest framework
- **ga_view/**: Qt6-based 2D visualization tool for GA objects and operations
- **ga_lua/**: Lua scripting interface for GA operations
- **ga_prdxpr/**: Code generator for GA product expressions (EGA2D, EGA3D, PGA2DP, PGA3DP)

### Library Usage Patterns

The library provides two main entry points:

- `#include "ga/ga_ega.hpp"` for Euclidean GA (2D, 3D, 4D)
- `#include "ga/ga_pga.hpp"` for Projective GA (2dp, 3dp)

Critical usage requirements:

1. All GA code must be wrapped in try-catch blocks to handle exceptions
2. Use namespace `hd::ga` and either `hd::ga::ega` or `hd::ga::pga`
3. **Always enclose left/right contractions (`<<` and `>>` operators) in parentheses due to operator precedence**
4. The base scalar type (float/double) is configured globally in `ga/ga_value_t.hpp`

### Type System Architecture

- Template types: `Scalar2d<T>`, `Vec2d<T>`, `BiVec2d<T>`, `MVec2d_E<T>`, etc.
- User convenience types: `scalar2d`, `vec2d`, `bivec2d`, `mvec2d`, etc. (based on `value_t`)
- Multivector types differentiate between even (`_E`), odd (`_U`), and full (`MVec`)
  multivectors
- Function call arguments for scalar and pscalar use `Scalar2d<T> a`, etc. directly
  instead of `Scalar2d<T> const& a` (cheap to copy), while for vector and bivector, etc. use
  `Vec2d<T> const&` for unmodifiable arguments
- ALWAYS use `Vec3dp<T> const&` instead of `const Vec3dp<T>&` etc. (i.e. use "east const" convention)

### Dependencies

Required: fmt library (header-only)
Optional: doctest (testing), Qt6 (visualization), Lua + sol2 (scripting), readline
(optional for scripting)

Install on macOS: `brew install fmt doctest qt6 lua`

### Development Notes

- Compiler definitions: `-D_HD_GA_EXTENDED_TEST_DIV_BY_ZERO` (extended testing), `-D_HD_GA_HAVE_READLINE`(readline support)
- The library supports switching between debug/release builds via CMAKE_BUILD_TYPE
- MSVC uses `/bigobj` flag due to template instantiation complexity
- Additional compiler definitions start with `_HD_GA_` in order to be consistent
  with already existing definitions
- Use the existing cmake infrastructure as much as possible to test changes

## GA Product Expression Generator (ga_prdxpr/)

The `ga_prdxpr/` directory contains a sophisticated **code generator** that produces optimized C++ expressions for geometric algebra operations. It supports four complete algebras: EGA2D, EGA3D, PGA2DP, and PGA3DP.

### Key Components

**Configuration System:**

- `ga_prdxpr_*.hpp`: User defined information on algebra, rules and coefficients
- `ga_prdxpr_*_config.cpp`: Algebra-specific configurations defining basis, coefficients,
  and product rules to be used for product generation
- **AlgebraConfig**: Complete algebra specification (basis elements, complement rules, etc.)
- **ProductDefinition**: Configures product types (geometric, wedge, contraction, sandwich, etc.)
- **ProductCase**: Individual expression cases with specific coefficient combinations

**Generation Engine:**

- `ga_prdxpr_generator.cpp`: Main generation logic with algebra-specific handlers
- **Dimensional dispatch**: 2D, 3D, 4D specialized implementations
- **Product-specific handlers**: Each ProductType gets specialized mathematical treatment

### Critical Mathematical Patterns

**Complement Transformations for regressive products:**

- for odd dimensinal algebras ega3d, pga2dp: `cmpl(operation(cmpl(A), cmpl(B)))`
- for even dimensinal algebras ega2d, pga3dp: `lcmpl(operation(rcmpl(A), rcmpl(B)))`

**Coefficient Usage Patterns:**

```cpp
A/B                   // Pure multivector × multivector operations
A_even/B_even         // Even part operations (rotors)
R_even/R_rev_even     // Rotor operations (EGA)
R_odd/R_rrev_odd      // Motor operations (PGA)
svps/svps             // Symmetric scalar/vector/pseudoscalar
svps1/svps2           // Asymmetric patterns (e.g., v1.x*v2.y)
```

**Sandwich Products (Two-Step Pattern):**

1. `rotor/motor * object → intermediate_result` (with parentheses via `brace_switch::use_braces`)
2. `intermediate_result * rev/rrev(rotor/motor) → final_result`

### Algebra-Specific Behaviors

**EGA2D/EGA3D (Euclidean):**

- Standard complement sequences
- Rotor-based sandwich products with `rev()` operations based on geometric product `gpr()`
- Coefficient pattern: R_even/R_rev_even

**PGA2DP/PGA3DP (Projective):**

- Regressive complement transformations
- Motor-based sandwich products with `rrev()` operations based on regressive geometric
  product `rgpr()`

**Success Criteria:**

- Character-perfect match with reference implementation
- All product types working (geometric, wedge, contraction, sandwich, etc.)
- Mathematical accuracy verified against known GA identities

### Common Pitfalls

1. **Coefficient Mismatches**: Each ProductType needs algebra-specific coefficient analysis
2. **Operator Precedence**: Always parenthesize left/right contractions (`<<`, `>>`) in GA expressions
3. **Output Format**: Case descriptions must match mathematical reality exactly

### Key Files for Modifications in ga_prdxpr subfolder

- `ga_prdxpr_generator.cpp`: Core generation logic and product-specific handlers
- `ga_prdxpr_*_config.cpp`: Algebra-specific configurations and coefficient definitions
- `ga_prdxpr_*.cpp`: Legacy reference implementations (don't modify, use for verification)

**Important**: The ga_prdxpr system is a **complete, production-ready geometric algebra code generator** that produces mathematically accurate, optimized C++ expressions for all supported algebras.

### Static Initialization Order Safety in ga_prdxpr

**Critical Lesson**: ga_prdxpr experienced static initialization order fiasco (segmentation faults) during complement rule generation implementation. Key safety principles learned:

**Root Causes of Static Initialization Issues:**

1. **Static string constants** in headers that depend on each other
2. **Circular references** between extern declarations and const initializations
3. **Cross-file static dependencies** during rule generation

**Required Safety Patterns:**

1. **Inline Functions for Constants**: Convert all static string constants to inline functions:

   ```cpp
   // UNSAFE: Static variables
   static const std::string one_str{"1"s};
   
   // SAFE: Inline functions  
   inline const std::string& one_str() { static const std::string s{"1"s}; return s; }
   ```

2. **Staged Initialization**: Generate all rules in single pass, then assign to const variables:

   ```cpp
   // Stage 1: Generate all rules together (no dependencies)
   static auto generated_rules = generate_algebra_rules(config);
   
   // Stage 2: Assign to const variables (safe, already generated)
   const prd_rules gpr_rules = generated_rules.geometric_product;
   const prd_rules complement_rules = generated_rules.complement;
   ```

3. **Explicit Values for Cross-Referenced Rules**: Never reference extern variables in const initializations:

   ```cpp
   // UNSAFE: References extern before initialization
   extern const prd_rules complement_rules;
   const prd_rules dual_rules = complement_rules;  // BAD!
   
   // SAFE: Explicit values
   const prd_rules dual_rules = {{"1", "e123"}, {"e1", "e23"}, ...};
   ```

**Future Dual Rule Implementation**: When adding dual rules, maintain these safety patterns to avoid recreating static initialization order issues.

## Congruence Testing in Geometric Algebra

**Key Mathematical Insight**: Elements are congruent if they represent the same geometric subspace up to scalar multiplication, **regardless of sign or magnitude**. The congruence test varies by grade:

**Congruence Rules by Grade:**

- **Grade 0 (Scalars)**: All non-zero scalars represent the same 0-dimensional subspace
- **Grade 1 (Vectors)**: Congruent if same subspace, e.g. parallel/antiparallel
- **Grade 2+ (Bivectors, etc.)**: Congruent if same subspace, e.g. parallel/antiparallel
- **Grade n (Pseudoscalars)**: All non-zero pseudoscalars represent the same n-dimensional subspace

**Implementation Pattern**: Use type-specific function overloads with `requires` clauses instead of generic templates, since wedge products between same-grade high-order elements are often undefined (e.g., `wdg(BiVec3d, BiVec3d)` → grade 4, invalid in 3D).

**Critical**: This logic applies to all algebras (EGA2D, EGA3D, PGA2DP, PGA3DP) but must account for algebra-specific grade structures and available operations.

- **EGA2D**: `is_congruent2d()` for `Scalar2d`, `Vec2d`, `PScalar2d`
- **EGA3D**: `is_congruent3d()` for `Scalar3d`, `Vec3d`, `BiVec3d`, `PScalar3d`
- **PGA2DP**: `is_congruent2dp()` for `Scalar2dp`, `Vec2dp`, `BiVec2dp`, `PScalar2dp`
- **PGA3DP**: `is_congruent3dp()` for `Scalar3dp`, `Vec3dp`, `BiVec3dp`, `TriVec3dp`, `PScalar3dp`

**Unified Algorithm**:

1. Component-wise zero detection (avoids issues with degenerate PGA metrics)
2. Find scale factor `k` from first non-zero component pair where `a = k*b`
3. Verify `a = k*b` holds for all components within tolerance
4. Special cases: scalars and pseudoscalars at max grade always congruent if non-zero

## Geometric Algebra Mathematical Foundations

### Core GA Concepts Learned

**Geometric Algebra Structure:**

- **G(p,q,r)** notation: p positive, q negative, r null (degenerate) basis vectors
- **Basis Elements**: Ordered combinations like `e1`, `e12`, `e123` with canonical ordering
- **Grade Structure**: scalar (grade 0), vectors (grade 1), bivectors (grade 2), trivectors (grade 3), etc.
- **Metric Signature**: Defines how basis vectors square (e.g., e1² = +1, e2² = -1, e3² = 0)

**Product Operations:**

1. **Geometric Product**: Full GA multiplication with metric effects
2. **Wedge Product**: Antisymmetric exterior product (zero for repeated indices)
3. **Dot Product**: Symmetric contraction using extended metric values

**Indexing Systems:**

- **1-based**: e1, e2, e3 (traditional mathematical notation)
- **0-based**: e0, g0, g1, g2 (often used in physics/spacetime)
- Algorithm must detect and handle both automatically

### Supported Algebra Types

**EGA2D - G(2,0,0)**: Euclidean 2D

- Basis: `{1, e1, e2, e12}`
- Metric: `{+1, +1}`
- Extended Metric: `{1, 1, 1, 1}`

**EGA3D - G(3,0,0)**: Euclidean 3D

- Basis: `{1, e1, e2, e3, e23, e31, e12, e123}`
- Metric: `{+1, +1, +1}`
- Extended Metric: `{1, 1, 1, 1, 1, 1, 1, 1}`

**PGA2DP - G(2,0,1)**: Projective 2D

- Basis: `{1, e1, e2, e3, e23, e31, e12, e321}`
- Metric: `{+1, +1, 0}`
- Extended Metric: `{1, 1, 1, 0, 0, 0, 1, 0}`

**PGA3DP - G(3,0,1)**: Projective 3D

- Basis: `{1, e1, e2, e3, e4, e41, e42, e43, e23, e31, e12, e423, e431, e412, e321, e1234}`
- Metric: `{+1, +1, +1, 0}`
- Extended Metric: `{1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0}`

**STA3D - G(1,3,0)**: Space-Time Algebra

- Basis: `{1, g0, g1, g2, g3, g01, g02, g03, g23, g31, g12, g023, g031, g012, g123, g0123}`
- Metric: `{+1, -1, -1, -1}` (timelike g0, spacelike g1,g2,g3)
- Extended Metric: Special mixed signature rules

### Extended Metric Calculation Rules

**Standard Algebras (EGA, PGA):**

- **Vectors**: Use base metric signature directly
- **Higher grades**: Product of constituent vector metrics
- **Pseudoscalar**: Determinant of metric tensor
- **Mixed signatures**: Apply conforming property G(a ∧ b) = G(a) ∧ G(b)

**Space-Time Algebra (STA3D) Special Rules:**

- **Bivectors with g0**: Extended metric = +1 (e.g., g01, g02, g03)
- **Bivectors without g0**: Extended metric = -1 (e.g., g12, g23, g31)
- **Trivectors with g0**: Extended metric = -1 (e.g., g012, g023, g031)
- **Trivectors without g0**: Extended metric = +1 (e.g., g123)
- **Pseudoscalar**: Extended metric = determinant = -1

### Automatic Rule Generation System

**Core Algorithm Components:**

- **AlgebraConfig**: Defines basis vectors, metric signature, canonical basis ordering
- **multiply_basis_elements()**: Implements geometric product with sign tracking
- **Extended Metric Calculator**: Handles dot product values for all basis elements
- **Validation System**: Character-by-character comparison with reference implementations

**Key Mathematical Insights:**

1. **Sign Calculation**: Count swaps needed for canonical ordering
2. **Metric Application**: Handle both 1-based and 0-based indexing automatically
3. **Canonical Forms**: Use user-provided basis ordering to avoid ambiguity
4. **Mixed Signatures**: Require algebra-specific extended metric rules

**Validation Success**: The automatic generation system produces character-identical results for all five tested algebras (EGA2D, EGA3D, PGA2DP, PGA3DP, STA3D) across geometric, wedge, and dot products.

## Congruence Implementation and Numerical Precision

### Congruence Testing for Geometric Algebra Elements

The library implements unified congruence testing using the `is_congruent*()` functions to determine when two GA elements represent the same subspace (i.e., are scalar multiples of each other).

**Implementation Pattern (Unified A = k*B Approach):**

```cpp
bool is_congruent(const T& a, const U& b, value_t tolerance = eps)
{
    // 1. Handle zero cases using component-wise checks
    // 2. Find scale factor k where a = k*b 
    // 3. Verify a = k*b for all components within tolerance
}
```

**Type-Specific Logic:**

- **Scalars/Pseudoscalars**: All non-zero elements of same grade are congruent (represent same-dimensional subspace)
- **Vectors/Bivectors/Trivectors**: Use component-wise proportionality testing with `A = k*B` relationship

### Critical Numerical Precision Lessons

**Problem**: Initial implementation used absolute tolerance (`eps ≈ 1.11e-15`) which
failed for floating-point calculations involving:

- Wedge products followed by division operations
- Complex geometric algebra expressions with accumulated numerical errors

**Solution**: Implemented **relative tolerance** scaling:

```cpp
// OLD: Absolute tolerance (too strict)
return (std::abs(a.x - k * b.x) < tolerance);

// NEW: Relative tolerance (numerically stable)
value_t rel_tol = tolerance * std::max({std::abs(a.x), std::abs(a.y), ..., value_t(1.0)});
return (std::abs(a.x - k * b.x) < rel_tol);
```

**Key Insight**: When working with GA operations that involve multiple floating-point calculations (wedge products, divisions, etc.), accumulated numerical errors can exceed extremely tight absolute tolerances. Relative tolerance provides better numerical stability while maintaining mathematical correctness.

**Applied Consistently Across All Algebras:**

- EGA2D: Vector functions
- EGA3D: Vector and bivector functions  
- PGA2DP: Vector and bivector functions
- PGA3DP: Vector, bivector, and trivector functions
