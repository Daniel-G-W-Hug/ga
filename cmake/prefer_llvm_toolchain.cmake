# Prefer an installed LLVM/clang over the platform default.
#
# MUST be included BEFORE project() -- that is where CMake resolves and caches the
# compiler; afterwards this is a no-op. Kept as a module rather than inline so an
# enclosing build, whose own project() call runs first, can apply the same rule before
# its project() line: `include(external/ga/cmake/prefer_llvm_toolchain.cmake)` (adjust
# the path to wherever this repo sits).
#
# Rationale: the vendor compiler shipped with the platform trails upstream LLVM by a
# release or more, and this project is built against the newer one -- the repo's own
# formatting baseline is clang-format 22, and the C++23 library features used here land
# earlier there. A bare `cmake ..` otherwise silently selects the older vendor compiler,
# which is invisible in the build log and only shows up as a missing feature or as
# results that do not compare against ones produced elsewhere.
#
# Deliberately conservative -- it only fills a gap, never overrides a choice:
#   * skipped entirely on a Windows host, where MSVC is the platform toolchain;
#   * skipped when the caller passed -DCMAKE_CXX_COMPILER=... or set CXX in the
#     environment, or when a cache already exists from an earlier configure;
#   * NO_DEFAULT_PATH, so it picks up a real LLVM installation and never re-selects
#     the vendor `clang++` sitting first on PATH;
#   * a no-op when no such installation is found (the platform default stands).
# Set LLVM_ROOT to point it at a specific installation.

if(NOT CMAKE_HOST_WIN32 AND NOT CMAKE_CXX_COMPILER AND NOT DEFINED ENV{CXX})
    find_program(_GA_LLVM_CXX
        NAMES clang++
        HINTS $ENV{LLVM_ROOT}/bin
              /opt/homebrew/opt/llvm/bin   # Homebrew, Apple silicon
              /usr/local/opt/llvm/bin      # Homebrew, Intel
              /opt/local/libexec/llvm/bin  # MacPorts
              /usr/lib/llvm/bin            # common Linux prefix
        NO_DEFAULT_PATH)
    if(_GA_LLVM_CXX)
        set(CMAKE_CXX_COMPILER "${_GA_LLVM_CXX}" CACHE FILEPATH
            "C++ compiler (LLVM preferred over the platform default)")
        message(STATUS "Preferring LLVM toolchain: ${_GA_LLVM_CXX}")
    endif()
endif()
