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

# set(CMAKE_VERBOSE_MAKEFILE ON)  // alternative: make VERBOSE=1
# add_definitions(-Dxxx)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# set a default build type: Debug | RelWithDebInfo | Release | MinSizeRel
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
  message(STATUS "Build type not specified: use Debug mode by default.")
endif()

if(MSVC)
  # warning level 4 and all warnings as errors
  add_compile_options(/bigobj /Wall /W4 /WX /EHs)
else()
  # lots of warnings and all warnings as errors
  #add_compile_options(-Wall -Wextra -pedantic -Werror)
  # add_compile_options(-Wall -Wextra -Wconversion -pedantic)
  # add_compile_options(-Wall -Wextra -Wshadow -pedantic)
  add_compile_options(-Wall -Wextra -pedantic)
endif()

### Extended testing for division by zero (optional, may cost speed)
add_definitions(-D_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)

#################################################################
# Cross-platform dependency management
#################################################################

# Add our custom cmake modules to the module path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Include and setup dependencies
include(dependencies)
setup_ga_dependencies()

# Optional readline support for ga_lua interactive shell
option(GA_USE_READLINE "Enable readline support for ga_lua interactive shell" ON)

#################################################################
# Project modules
#################################################################
 
add_subdirectory(ga)                     # ga library (header only)
add_subdirectory(ga_test)                # tests for header only library
add_subdirectory(ga_lua)                 # lua scripting for ga
add_subdirectory(ga_view)                # qt-based viewer
add_subdirectory(ga_prdxpr)              # generate product expressions

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
if(TARGET sol2::sol2)
    message(STATUS "sol2: ✓ Available")
else()
    message(STATUS "sol2: ✗ Missing")
endif()
if(TARGET Qt6::Core)
    message(STATUS "Qt6: ✓ Available")
else()
    message(STATUS "Qt6: ✗ Missing")
endif()
if(TARGET ${LUA_LIBRARIES} OR Lua_FOUND)
    message(STATUS "Lua: ✓ Available")
else()
    message(STATUS "Lua: ✗ Missing")
endif()
if(GA_HAS_HD_UTILS)
    message(STATUS "hd utils: ✓ Available")
else()
    message(STATUS "hd utils: ✗ Missing")
endif()
message(STATUS "")