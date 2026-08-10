# Cross-Platform Build Usage Guide

> **This is the public `ga` repository** — the header-only geometric-algebra library, the
> `ga_prdxpr` code generator, the `ga_lua` scripting shell, the `ga_py` Python wrapper, and
> the test suite. Application-specific bundles are not part of this repository.

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
# Default: prefer system installs, FetchContent only the missing supported deps
# (fmt / doctest / sol2). This is GA_USE_SYSTEM_DEPENDENCIES=ON (the default) with
# GA_FORCE_FETCH_CONTENT=OFF — nothing to pass.
cmake ..

# Force FetchContent for the supported deps even if they are installed system-wide
# (reproducible CI/CD builds)
cmake .. -DGA_FORCE_FETCH_CONTENT=ON

# Ignore system installs and FetchContent the supported deps instead (equivalent to
# the force flag for those deps)
cmake .. -DGA_USE_SYSTEM_DEPENDENCIES=OFF

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
