# Cross-Platform CMake Dependency Management

This directory contains modular CMake scripts for flexible dependency management across Windows, macOS, and Linux.

> This is the **public `ga` repository**: the header-only library, `ga_prdxpr`, `ga_lua`,
> `ga_py`, and the test suite. It needs **no Qt6** — the Qt6 viewer lives in a private
> superset repo that embeds this one as a submodule.

## Features

- **Three-Tier Resolution**: Automatically resolves dependencies in order: system → ../../include/ → FetchContent (readline prefers system)
- **Platform Awareness**: Different strategies for different platforms
- **User Control**: Options to force system-only or FetchContent-only approaches
- **Maximum Reuse**: Common dependency resolution logic shared across platforms
- **Local Development**: Supports local dependencies in ../../include/ for development workflows

## Usage Options

### Option 1: Default (Three-tier approach)

```cmake
# Uses system deps → ../../include/ deps → FetchContent for missing ones
cmake .. 
```

### Option 2: Force FetchContent

```cmake
# Forces FetchContent for all supported dependencies (fmt, doctest, sol2)
cmake .. -DGA_FORCE_FETCH_CONTENT=ON
```

### Option 3: FetchContent instead of system installs

```cmake
# Ignore system-installed fmt/doctest/sol2 and FetchContent them instead
# (equivalent to GA_FORCE_FETCH_CONTENT=ON for those deps)
cmake .. -DGA_USE_SYSTEM_DEPENDENCIES=OFF
```

### Option 4: Disable readline

```cmake
# Disables optional readline support
cmake .. -D_GA_USE_READLINE=OFF
```

## Dependencies by Category

### System-Installed (per target)

- **Lua 5.1+**: must be installed separately, for the `ga_lua` scripting shell
  (`_GA_USE_LUA`, ON by default; the target is skipped if Lua is absent)

### Flexible (System → ../../include/ → FetchContent)

- **fmt**: Formatting library (header-only preferred)
- **doctest**: Testing framework (header-only)
- **sol2**: Lua C++ binding (header-only)

### Optional (System → ../../include/ → skip)

- **readline**: Enhanced command-line editing (three-tier detection, no FetchContent fallback)

## Platform-Specific Installation

### macOS (Homebrew)

```bash
# For the Lua shell
brew install lua

# fmt/doctest/readline (if not using ../../include/ or FetchContent)
brew install fmt doctest readline
```

### Windows (vcpkg recommended)

```bash
# Lua is linked from a system install (see the Windows section of CLAUDE.md)
vcpkg install fmt doctest sol2 readline-win32 --triplet x64-windows
```

### Linux (Ubuntu/Debian)

```bash
# For the Lua shell
sudo apt-get install lua5.1-dev

# fmt/doctest/readline (if not using ../../include/ or FetchContent)
sudo apt-get install libfmt-dev libdoctest-dev libreadline-dev
```

## Architecture

- `dependencies.cmake`: Main coordinator, handles user options and three-tier resolution
- `find_dependencies.cmake`: Three-tier dependency detection (system → ../../include/ → FetchContent)
- `fetch_dependencies.cmake`: FetchContent setup with version pinning and configuration

## Three-Tier Dependency Resolution

1. **Tier 1: System Dependencies** - Uses `find_package()` to locate system-installed libraries
2. **Tier 2: Local Dependencies** - Checks `../../include/` for header-only libraries (development workflow)
3. **Tier 3: FetchContent** - Downloads and builds dependencies from source as fallback

## Benefits

1. **Developer Flexibility**: Works whether dependencies are system-installed, local, or downloaded
2. **Local Development**: Supports ../../include/ workflow for development with local dependencies
3. **CI/CD Friendly**: Can force FetchContent for reproducible builds
4. **Beginner Friendly**: Clear error messages with installation instructions
5. **Cross-Platform**: Handles platform differences automatically
6. **Maintainable**: Centralized dependency logic, easy to update versions
