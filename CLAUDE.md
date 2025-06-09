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
- **ga_prdxpr/**: Code generator for GA product expressions

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
