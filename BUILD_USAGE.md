# Cross-Platform Build Usage Guide

> **This is the public `ga` repository** — the header-only geometric-algebra library, the
> `ga_prdxpr` code generator, the `ga_lua` scripting shell, the `ga_py` Python wrapper, and
> the test suite. The Qt6 visualization viewer and the application bundles (wafer grinding,
> magnetic levitation) live in a private superset repo that embeds this one as a submodule,
> so **nothing here needs Qt6**.

## Quick Start

```bash
mkdir -p build && cd build
cmake ..
cmake --build .
```

This automatically:

- uses system dependencies when available (macOS: Homebrew, Windows: vcpkg),
- falls back to FetchContent for missing supported dependencies,
- prints a summary of what was found / downloaded.

## Build Options

```bash
# Force FetchContent for reproducible builds (CI/CD)
cmake .. -DGA_FORCE_FETCH_CONTENT=ON

# System-only mode (fail if a dependency is missing rather than fetching it)
cmake .. -DGA_USE_SYSTEM_DEPENDENCIES=ON -DGA_FORCE_FETCH_CONTENT=OFF

# Skip the Lua shell / its readline support
cmake .. -D_GA_USE_LUA=OFF
cmake .. -D_GA_USE_READLINE=OFF

# Build the Python wrapper (ga_py) as well — needs nanobind in the ga_py venv,
# see ga_py/README.md
cmake .. -D_GA_BUILD_PYTHON=ON

# Release build (recommended for the heavier test suites)
cmake .. -DCMAKE_BUILD_TYPE=Release
```

## Dependencies

Only **fmt** is needed to use the library (header-only; FetchContent fallback if it is not
installed). Everything else is optional and gates only its own target — a missing optional
dependency simply skips that target and the build stays green.

| Dependency        | Enables                                   | macOS (brew)     | Linux (apt)        | Fallback              |
|-------------------|-------------------------------------------|------------------|--------------------|-----------------------|
| fmt               | the library + every target                | `fmt`            | `libfmt-dev`       | FetchContent          |
| doctest           | the test suite (`ga_test`)                | `doctest`        | `libdoctest-dev`   | FetchContent          |
| Lua + sol2        | the Lua shell (`ga_lua`, `_GA_USE_LUA`)   | `lua`            | `lua5.1-dev`       | sol2 via FetchContent |
| readline          | nicer `ga_lua` interactive editing        | `readline`       | `libreadline-dev`  | optional              |
| Python + nanobind | the Python wrapper (`ga_py`, opt-in `-D_GA_BUILD_PYTHON=ON`) | `python` | `python3-dev` | see `ga_py/README.md` |

## Platform-Specific Setup

### macOS (Homebrew)

```bash
brew install fmt doctest lua readline
```

### Windows (vcpkg + system Lua)

```bash
vcpkg install fmt doctest sol2 readline-win32 --triplet x64-windows
# Lua is linked from a system install — see the Windows section of CLAUDE.md.
```

### Linux (Ubuntu/Debian)

```bash
sudo apt-get install libfmt-dev libdoctest-dev lua5.1-dev libreadline-dev
```

## Configuration Summary

After `cmake ..` you get a summary of what was resolved, e.g.:

```text
GA Project Configuration Summary:
=================================
Build type: Debug
C++ standard: 23
fmt: ✓ Available
doctest: ✓ Available
sol2: ✓ Available
Lua: ✓ Available
```

This shows at a glance what is enabled and what might need attention (an unavailable
optional dependency just means its target is skipped).
