# Cross-Platform CMake Dependency Management

This directory contains modular CMake scripts for flexible dependency management across Windows, macOS, and Linux.

## Features

- **Three-Tier Resolution**: Automatically resolves dependencies in order: system → ../../include/ → FetchContent
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

### Option 3: System Only
```cmake
# Only uses system dependencies, skips ../../include/ and FetchContent
cmake .. -DGA_USE_SYSTEM_DEPENDENCIES=ON -DGA_FORCE_FETCH_CONTENT=OFF
```

### Option 4: Disable readline
```cmake
# Disables optional readline support
cmake .. -DGA_USE_READLINE=OFF
```

## Dependencies by Category

### Always System-Installed (Required)
- **Qt6**: Must be installed separately (GUI framework)  
- **Lua 5.1+**: Must be installed separately (scripting runtime)

### Flexible (System → ../../include/ → FetchContent)
- **fmt**: Formatting library (header-only preferred)
- **doctest**: Testing framework (header-only)
- **sol2**: Lua C++ binding (header-only)

### Optional (System → ../../include/ → skip)
- **readline**: Enhanced command-line editing (three-tier detection, no FetchContent fallback)

## Platform-Specific Installation

### macOS (Homebrew)
```bash
# Required
brew install qt6 lua

# Optional (if not using ../../include/ or FetchContent)
brew install fmt doctest readline
```

### Windows (vcpkg recommended)
```bash
# Required - install Qt6 from qt.io and Lua from lua.org

# Optional (if not using ../../include/ or FetchContent)
vcpkg install fmt doctest
```

### Linux (Ubuntu/Debian)
```bash
# Required
sudo apt-get install qt6-base-dev qt6-tools-dev lua5.1-dev

# Optional (if not using ../../include/ or FetchContent)  
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