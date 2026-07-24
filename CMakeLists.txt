# Create tests for header-only files in namespace hd

##################################################################
# usage: execute "cd build && cmake .." from build directory
#
# to build:  "cmake --build ." from build directory
##################################################################
#]]
#[[ hint: backet comment will not be formated by cmake-format]]

cmake_minimum_required(VERSION 3.28)

#
# for trial of new features, libs or system values"
#
project(
  ga
  VERSION 0.1
  LANGUAGES CXX)

# Root of THIS (public `ga`) source tree. Use GA_ROOT -- not CMAKE_SOURCE_DIR -- for the
# repo-root include path and exporter paths in the sub-CMakeLists, so the build works both
# standalone AND when this repo is added as a subdirectory of a superproject (the private
# `ga_private` superset adds it as the `external/ga` submodule, where CMAKE_SOURCE_DIR is
# the superproject root, not this one). CMAKE_CURRENT_SOURCE_DIR here is always this repo's
# root. Propagates to all add_subdirectory descendants as a normal variable.
set(GA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

# set(CMAKE_VERBOSE_MAKEFILE ON)  // alternative: make VERBOSE=1
# add_definitions(-Dxxx)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS OFF)  # Disable GNU extensions (-std=c++23 instead of -std=gnu++23)

# Build type (Debug | Release | RelWithDebInfo | MinSizeRel).
#
# Single-config generators (Ninja, Unix Makefiles -- what VS Code's CMake Tools uses) read
# CMAKE_BUILD_TYPE. Defaulting it through the CACHE (not a plain set) lets the IDE / cmake-gui
# see and change it, and the STRINGS property gives the IDE a build-type dropdown, so Debug
# and Release are selectable directly from VS Code (status-bar variant or the CMake cache
# editor). An explicit -DCMAKE_BUILD_TYPE=... (e.g. from CMake Tools) overrides the default.
# Multi-config generators (Xcode, Visual Studio) ignore CMAKE_BUILD_TYPE and choose the
# configuration at build time (CMAKE_CONFIGURATION_TYPES), so the default is skipped there.
get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT _is_multi_config AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Build type: Debug | Release | RelWithDebInfo | MinSizeRel" FORCE)
  message(STATUS "Build type not specified: defaulting to Debug.")
endif()
# The STRINGS dropdown only applies to single-config generators, where the cache variable
# above exists. Multi-config generators (Visual Studio, Xcode) never create CMAKE_BUILD_TYPE
# as a cache entry, so guard the property to avoid "could not find CACHE variable".
if(NOT _is_multi_config)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
               Debug Release RelWithDebInfo MinSizeRel)
endif()

if(MSVC)
  # warning level 4 and all warnings as errors
  add_compile_options(/bigobj /Wall /W4 /WX /EHs /utf-8)
else()
  # lots of warnings and all warnings as errors
  #add_compile_options(-Wall -Wextra -pedantic -Werror)
  # add_compile_options(-Wall -Wextra -Wconversion -pedantic)
  # -Wshadow catches the same family MSVC's /W4 reports as C4456-4459 (local/parameter/
  # member/global shadowing) -- enabled so a shadowing bug surfaces on every platform
  # instead of only failing the MSVC /WX build.
  add_compile_options(-Wall -Wextra -Wshadow -pedantic)
endif()

### Extended testing for division by zero (optional, may cost speed)
add_definitions(-D_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)

#################################################################
# Cross-platform dependency management
#################################################################

# Add our custom cmake modules to the module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Optional lua wrapper that provides lua scripting and an interactive ga_lua shell
option(_GA_USE_LUA "Build the lua_wrapper (requires sol2 and optionally readline)" ON)

# Optional readline support for ga_lua interactive shell.
# Forced OFF (and hidden from the cache UI) when _GA_USE_LUA is OFF, since
# readline is only consumed by ga_lua/src/ga_lua.cpp.
include(CMakeDependentOption)
cmake_dependent_option(_GA_USE_READLINE
    "Enable readline support for interactive ga_lua shell" ON
    "_GA_USE_LUA" OFF)

# Optional Python wrapper (nanobind module + cross-check exporter).
# requires Python 3.10+ and nanobind installed in the
# Python environment that CMake's find_package(Python) discovers.
option(_GA_BUILD_PYTHON "Build the ga_py Python wrapper and its test-data exporter" ON)

# Include and setup dependencies
include(dependencies)
setup_ga_dependencies()

#################################################################
# Project modules
#################################################################

add_subdirectory(ga)                     # ga library (header only)
add_subdirectory(ga_test)                # tests for header only library

if(_GA_USE_LUA)
    add_subdirectory(ga_lua)             # lua scripting for ga
endif()

# ga_view is a private module (it moved to ga_private). The public repo has no
# ga_view/ of its own; the private superset supplies it at ${GA_PRIVATE_DIR}/ga_view
# and attaches it through the GA_PRIVATE_DIR overlay below. The local-directory
# branch is kept for a checkout that still carries ga_view/ (e.g. before the move).
#
# Adding the private viewer from HERE (this scope) rather than from the superset's
# own top-level CMake is deliberate: the Qt6/fmt/doctest/ga dependency targets are
# created by this project's dependency system in this scope and are not visible in
# the parent superset scope, so a sibling add_subdirectory there cannot link them.
#
# Skip on macOS with GCC due to a Qt6 header incompatibility (Qt6 uses _Pragma in
# function signatures which GCC rejects but Clang allows; may work on Linux).
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ga_view/CMakeLists.txt")
    set(_GA_VIEW_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ga_view")
elseif(GA_PRIVATE_DIR AND EXISTS "${GA_PRIVATE_DIR}/ga_view/CMakeLists.txt")
    set(_GA_VIEW_DIR "${GA_PRIVATE_DIR}/ga_view")
endif()
if(_GA_VIEW_DIR)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND APPLE)
        message(STATUS "ga_view: Skipped (Qt6 headers incompatible with GCC on macOS)")
        message(STATUS "         Use Clang to build ga_view, or build on Linux")
    else()
        # Binary dir is named so it works whether ga_view is in-tree or external.
        add_subdirectory("${_GA_VIEW_DIR}" ga_view)   # qt-based viewer
    endif()
endif()

add_subdirectory(ga_prdxpr)              # generate product expressions

if(_GA_BUILD_PYTHON)
    add_subdirectory(ga_py)              # nanobind Python wrapper
endif()

#################################################################
# Summary
#################################################################

# Print dependency summary
message(STATUS "")
message(STATUS "GA Project Configuration Summary:")
message(STATUS "=================================")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
if(TARGET fmt::fmt OR TARGET fmt::fmt-header-only)
    message(STATUS "fmt: ✓ Available")
else()
    message(STATUS "fmt: ✗ Missing")
endif()
if(TARGET doctest::doctest)
    message(STATUS "doctest: ✓ Available")
else()
    message(STATUS "doctest: ✗ Missing")
endif()
if(_GA_USE_LUA)
    if(TARGET sol2::sol2)
        message(STATUS "sol2: ✓ Available")
    else()
        message(STATUS "sol2: ✗ Missing")
    endif()
    if(Lua_FOUND OR LUA_FOUND OR GA_HAS_LUA)
        message(STATUS "Lua: ✓ Available")
    else()
        message(STATUS "Lua: ✗ Missing")
    endif()
else()
    message(STATUS "sol2: — Skipped (_GA_USE_LUA=OFF)")
    message(STATUS "Lua:  — Skipped (_GA_USE_LUA=OFF)")
endif()
if(TARGET Qt6::Core)
    message(STATUS "Qt6: ✓ Available")
else()
    message(STATUS "Qt6: ✗ Missing")
endif()
message(STATUS "")