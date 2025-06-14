# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

This project uses CMake (minimum version 3.28) with C++23 standard. **IMPORTANT: Always use out-of-source builds exclusively in the `build/` directory.**

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
6. Run executables: `./ga_lua/ga_lua`, `./ga_test/ga_ega_test`, etc.

## Running Tests

Execute test binaries from the build directory (ensure you're in `build/` directory):

```bash
./ga_test/ga_ega_test    # Tests for Euclidean GA
./ga_test/ga_pga_test    # Tests for Projective GA
```

## Running Applications

Execute applications from the build directory:

```bash
./ga_lua/ga_lua                    # Interactive Lua shell
./ga_lua/ga_lua script.lua         # Run Lua script
./ga_view/ga_view                  # Qt6 visualization tool
./ga_prdxpr/ga_prdxpr              # Code generator
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
- Multivector types differentiate between even (`_E`), odd (`_U`), and full (`MVec`) algebras

### Dependencies

Required: fmt library (header-only)
Optional: doctest (testing), Qt6 (visualization), Lua + sol2 (scripting), readline
(optional for scripting)

Install on macOS: `brew install fmt doctest qt6 lua`

### Development Notes

- Compiler definitions: `-D_HD_GA_EXTENDED_TEST_DIV_BY_ZERO` (extended testing), `-D_HD_GA_PRINT_WITH_TYPE_INFO` (verbose printing)
- The library supports switching between debug/release builds via CMAKE_BUILD_TYPE
- MSVC uses `/bigobj` flag due to template instantiation complexity
- Additional compiler definitions should start with `_HD_GA_` in order to be consistent
  with already existing definitions
- Use the existing cmake infrastructure as much as possible to test changes

## GA Product Expression Generator (ga_prdxpr/)

The `ga_prdxpr/` directory contains a sophisticated **code generator** that produces optimized C++ expressions for geometric algebra operations. It supports four complete algebras: EGA2D, EGA3D, PGA2DP, and PGA3DP.

### Architecture Overview

**Hybrid Two-System Design:**
- **New System**: ProductDefinition-based configuration framework (recommended for new work)
- **Legacy System**: Individual function-based generators (maintained for compatibility)
- **Output**: Both systems generate identical 6928-line mathematical expressions

### Key Components

**Configuration System:**
- `ga_prdxpr_config_*.cpp`: Algebra-specific configurations defining basis, coefficients, and product rules
- **AlgebraConfig**: Complete algebra specification (basis elements, complement rules, etc.)
- **ProductDefinition**: Configures product types (geometric, wedge, contraction, sandwich, etc.)
- **ProductCase**: Individual expression cases with specific coefficient combinations

**Generation Engine:**
- `ga_prdxpr_generator.cpp`: Main generation logic with algebra-specific handlers
- **Dimensional dispatch**: 2D, 3D, 4D specialized implementations
- **Product-specific handlers**: Each ProductType gets specialized mathematical treatment

### Critical Mathematical Patterns

**Complement Transformations:**
- Standard: `cmpl(operation(cmpl(A), cmpl(B)))`
- PGA Regressive: `lcmpl(operation(rcmpl(A), rcmpl(B)))`
- **PGA3DP Exception**: Regressive sandwich uses asymmetric `rcmpl → gpr → lcmpl` sequence

**Coefficient Usage Patterns:**
```cpp
A/B                    // Pure multivector × multivector operations
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
- Rotor-based sandwich products with `rev()` operations
- Coefficient pattern: R_even/R_rev_even

**PGA2DP/PGA3DP (Projective):**
- Regressive complement transformations
- Motor-based sandwich products with `rrev()` operations
- **PGA3DP Critical**: Uses asymmetric complement sequence for regressive sandwich

### Usage & Verification

**Build and Run:**
```bash
cd build
./ga_prdxpr/ga_prdxpr          # Original reference implementation
./ga_prdxpr/ga_prdxpr_new      # New configuration-driven system
```

**Verification Commands:**
```bash
# Should generate exactly 6928 lines
wc -l output.txt

# Should show no differences (100% match required)
diff reference_output.txt current_output.txt
```

**Success Criteria:**
- Exactly 6928 lines of output
- Character-perfect match with reference implementation
- All product types working (geometric, wedge, contraction, sandwich, etc.)
- Mathematical accuracy verified against known GA identities

### Common Pitfalls

1. **Complement Confusion**: PGA algebras require specific complement sequences, especially PGA3DP
2. **Coefficient Mismatches**: Each ProductType needs algebra-specific coefficient analysis
3. **Operator Precedence**: Always parenthesize left/right contractions (`<<`, `>>`) in GA expressions
4. **Output Format**: Case descriptions must match mathematical reality exactly

### Key Files for Modifications
- `ga_prdxpr_generator.cpp`: Core generation logic and product-specific handlers
- `ga_prdxpr_config_*.cpp`: Algebra-specific configurations and coefficient definitions
- `ga_prdxpr_*.cpp`: Legacy reference implementations (use for verification)

**Important**: The ga_prdxpr system is a **complete, production-ready geometric algebra code generator** that produces mathematically accurate, optimized C++ expressions for all four supported algebras.
