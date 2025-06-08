# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build System

This project uses CMake (minimum version 3.28) with C++23 standard. Build from the `build/` directory:

```bash
# Configure
cd build
cmake ..

# Build
cmake --build .
```

## Running Tests

Execute test binaries from the build directory:

```bash
./ga_test/ga_ega_test    # Tests for Euclidean GA
./ga_test/ga_pga_test    # Tests for Projective GA
```

## Project Architecture

This is a header-only geometric algebra library (`ga/`) with accompanying modules:

- **ga/**: Core header-only GA library supporting both Euclidean (EGA) and Projective (PGA) geometric algebra
- **ga_test/**: Test suite using doctest framework
- **ga_view/**: Qt6-based 2D visualization tool for GA objects and operations
- **ga_lua/**: Lua scripting interface for GA operations
- **ga_prdxpr/**: Code generator for GA product expressions

### Library Usage Patterns

The library provides two main entry points:

- `#include "ga/ga_ega.hpp"` for Euclidean GA (2D, 3D, 4D)
- `#include "ga/ga_pga.hpp"` for Projective GA (2D+, 3D+)

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
Optional: doctest (testing), Qt6 (visualization), Lua + sol2 (scripting)

Install on macOS: `brew install fmt doctest qt6 lua`

### Development Notes

- Compiler definitions: `-D_HD_GA_EXTENDED_TEST_DIV_BY_ZERO` (extended testing), `-D_HD_GA_PRINT_WITH_TYPE_INFO` (verbose printing)
- The library supports switching between debug/release builds via CMAKE_BUILD_TYPE
- MSVC uses `/bigobj` flag due to template instantiation complexity
