# Create tests for header-only files in namespace hd

##################################################################
# usage: execute "cmake .. && make" from build directory
#    or: execute "cmake .. -G"Ninja" && ninja from build directory
#
# alternative to build:  "cmake --build ." from build directory
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
  add_compile_options(-Wall -Wextra -pedantic)
endif()

### Default is using duality as defined by Macdonald: dual(A)=A*rev(I)
### in order to use duality as defined by Hestenes, Doran and Lasenby: dual(A)=I*A
### uncomment this definition or set it externally
#add_definitions(-D_HD_GA_HESTENES_DORAN_LASENBY_DUALITY)

### Extended testing for division by zero (optional, may cost speed)
add_definitions(-D_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)

### Print with type info for scalars, vectors, bivectors, etc.
add_definitions(-D_HD_GA_PRINT_WITH_TYPE_INFO)

add_subdirectory(ga)      # ga library (header only)
add_subdirectory(ga_test) # tests for header only library
add_subdirectory(ga_lua)  # lua scripting for ga
add_subdirectory(ga_view) # qt-based viewer