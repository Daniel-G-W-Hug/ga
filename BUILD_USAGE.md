# Cross-Platform Build Usage Guide

## Quick Start

### Default Configuration (Recommended)
```bash
mkdir -p build && cd build
cmake ..
cmake --build .
```

This automatically:
- Uses system dependencies when available (macOS: homebrew, Windows: vcpkg)
- Falls back to FetchContent for missing dependencies
- Shows clear summary of what was found/downloaded

### Build Options

```bash
# Force FetchContent for reproducible builds (CI/CD)
cmake .. -DGA_FORCE_FETCH_CONTENT=ON

# System-only mode (fail if dependencies missing)
cmake .. -DGA_USE_SYSTEM_DEPENDENCIES=ON -DGA_FORCE_FETCH_CONTENT=OFF

# Disable optional readline support
cmake .. -DGA_USE_READLINE=OFF
```

## Platform-Specific Setup

### macOS (Homebrew)
```bash
# Required dependencies
brew install qt6 lua

# Optional (for system-first approach)
brew install fmt doctest readline
```

### Windows (vcpkg recommended)
```bash
# Install Qt6 from qt.io installer
# Install Lua from lua.org

# Optional (if you want system dependencies)
vcpkg install fmt doctest
```

### Linux (Ubuntu/Debian)
```bash
# Required dependencies
sudo apt-get install qt6-base-dev qt6-tools-dev lua5.1-dev

# Optional dependencies
sudo apt-get install libfmt-dev libdoctest-dev libreadline-dev
```

## Build System Benefits

### ✅ Fixed Issues
- **Lua Linking**: Fixed scoping issues with CMake functions vs macros
- **sol2 Compatibility**: Updated to development branch for C++23 support
- **Cross-Platform**: Works consistently across macOS, Windows, Linux

### ✅ New Features
- **Flexible Dependencies**: System-first with FetchContent fallback
- **User Control**: Multiple build configuration options
- **Clear Feedback**: Detailed dependency resolution messages
- **Version Pinning**: Stable versions for FetchContent dependencies

### ✅ Developer Experience
- **Just Works**: Default configuration handles most scenarios
- **CI/CD Ready**: Force FetchContent for reproducible builds
- **Informative**: Clear error messages with installation instructions
- **Fast**: Reuses system packages when available

## Migration from Old System

The new system is **fully backward compatible**. Your existing build workflow continues to work:

```bash
cd build
cmake ..
cmake --build .
```

No changes needed to your existing scripts or documentation.

## Dependency Resolution Summary

| Dependency | macOS | Windows | Linux | Fallback |
|------------|-------|---------|-------|----------|
| Qt6 | brew | qt.io installer | apt-get | ❌ Required |
| Lua | brew | lua.org | apt-get | ❌ Required |
| fmt | brew | vcpkg | apt-get | ✅ FetchContent |
| doctest | brew | vcpkg | apt-get | ✅ FetchContent |
| sol2 | ❌ | ❌ | ❌ | ✅ FetchContent |
| readline | brew | ❌ | apt-get | ❌ Optional |

## Configuration Summary Display

After running `cmake ..`, you'll see a helpful summary:

```
GA Project Configuration Summary:
=================================
Build type: Debug
C++ standard: 23
fmt: ✓ Available
doctest: ✓ Available  
sol2: ✓ Available
Qt6: ✓ Available
Lua: ✓ Available
```

This immediately shows you what's working and what might need attention.