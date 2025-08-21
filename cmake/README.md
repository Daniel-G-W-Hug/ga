# Cross-Platform CMake Dependency Management

This directory contains modular CMake scripts for flexible dependency management across Windows, macOS, and Linux.

## Features

- **Flexible Dependencies**: Automatically uses system packages when available, falls back to FetchContent
- **Platform Awareness**: Different strategies for different platforms
- **User Control**: Options to force system-only or FetchContent-only approaches
- **Maximum Reuse**: Common dependency resolution logic shared across platforms

## Usage Options

### Option 1: Default (Hybrid approach)
```cmake
# Uses system dependencies when available, FetchContent for missing ones
cmake .. 
```

### Option 2: Force FetchContent
```cmake
# Forces FetchContent for all supported dependencies (fmt, doctest, sol2)
cmake .. -DGA_FORCE_FETCH_CONTENT=ON
```

### Option 3: System Only
```cmake
# Only uses system dependencies, fails if any are missing
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

### Flexible (System or FetchContent)
- **fmt**: Formatting library (header-only preferred)
- **doctest**: Testing framework (header-only)
- **sol2**: Lua C++ binding (header-only)

### Optional
- **readline**: Enhanced command-line editing (system only)

## Platform-Specific Installation

### macOS (Homebrew)
```bash
# Required
brew install qt6 lua

# Optional (if using system dependencies)
brew install fmt doctest readline
```

### Windows (vcpkg recommended)
```bash
# Required - install Qt6 from qt.io and Lua from lua.org

# Optional (if using system dependencies)
vcpkg install fmt doctest
```

### Linux (Ubuntu/Debian)
```bash
# Required
sudo apt-get install qt6-base-dev qt6-tools-dev lua5.1-dev

# Optional (if using system dependencies)  
sudo apt-get install libfmt-dev libdoctest-dev libreadline-dev
```

## Architecture

- `dependencies.cmake`: Main coordinator, handles user options
- `find_dependencies.cmake`: System dependency detection with helpful error messages
- `fetch_dependencies.cmake`: FetchContent setup with version pinning and configuration

## Benefits

1. **Developer Flexibility**: Works whether dependencies are system-installed or not
2. **CI/CD Friendly**: Can force FetchContent for reproducible builds
3. **Beginner Friendly**: Clear error messages with installation instructions
4. **Cross-Platform**: Handles platform differences automatically
5. **Maintainable**: Centralized dependency logic, easy to update versions