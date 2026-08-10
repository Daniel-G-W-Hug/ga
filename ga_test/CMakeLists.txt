# incrementally add test files needed

set(EXEC_NAME1 ga_ega_test)
set(EXEC_NAME2 ga_pga_test)
set(EXEC_NAME3 ga_appl2dp_test)
set(EXEC_NAME4 ga_appl3dp_test)
set(EXEC_NAME5 ga_sta_test)
set(EXEC_NAME6 ga_integrator_test)
set(EXEC_NAME7 ga_stencil_test)
# Application-specific bundles (formerly EXEC_NAME9/10) are not part of this
# repository; an enclosing build can declare its own through the GA_PRIVATE_DIR
# overlay hook at the end of this file.
# ga_export_python_cases now lives in python_utilities/ (localized EXPORT_NAME)
# set(EXEC_NAME5 ga_sta_test)  # Commented out - STA not yet implemented

set(EGA_FILES
    src/ga_ega_test.cpp)

set(PGA_FILES
    src/ga_pga_test.cpp)

set(APPL2DP_FILES
    src/ga_appl2dp_test.cpp)

set(APPL3DP_FILES
    src/ga_appl3dp_test.cpp)

set(STA_FILES
    src/ga_sta_test.cpp)

set(INTEGRATOR_FILES
    src/ga_integrator_test.cpp)

set(STENCIL_FILES
    src/ga_stencil_test.cpp)

add_executable(${EXEC_NAME1} ${EGA_FILES})        #dep: doctest, fmt, ga
add_executable(${EXEC_NAME2} ${PGA_FILES})        #dep: doctest, fmt, ga
add_executable(${EXEC_NAME3} ${APPL2DP_FILES})    #dep: doctest, fmt, ga
add_executable(${EXEC_NAME4} ${APPL3DP_FILES})    #dep: doctest, fmt, ga
add_executable(${EXEC_NAME5} ${STA_FILES})        #dep: doctest, fmt, ga
add_executable(${EXEC_NAME6} ${INTEGRATOR_FILES}) #dep: doctest, fmt, ga
add_executable(${EXEC_NAME7} ${STENCIL_FILES})    #dep: doctest, fmt, ga

target_include_directories(${EXEC_NAME1} PRIVATE ${GA_ROOT})
target_include_directories(${EXEC_NAME2} PRIVATE ${GA_ROOT})
target_include_directories(${EXEC_NAME3} PRIVATE ${GA_ROOT})
target_include_directories(${EXEC_NAME4} PRIVATE ${GA_ROOT})
target_include_directories(${EXEC_NAME5} PRIVATE ${GA_ROOT})
target_include_directories(${EXEC_NAME6} PRIVATE ${GA_ROOT})
target_include_directories(${EXEC_NAME7} PRIVATE ${GA_ROOT})

# Dependencies are handled by the main CMakeLists.txt dependency system
# Just link against the targets that should be available

target_link_libraries(${EXEC_NAME1} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME2} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME3} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME4} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME5} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME6} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME7} PRIVATE doctest::doctest ga)

# doctest anonymous-symbol counter + its c2y warning.
#
# DOCTEST_COUNTER=__COUNTER__: on clang >= 22 doctest disables __COUNTER__ and falls back
# to __LINE__ for the names it generates per TEST_SUITE/TEST_CASE. With several test
# headers #included into one TU that collides whenever two registrations land on the same
# line in different headers (e.g. ga_appl3dp_*: a TEST_SUITE at line 29 in both the
# 3d_cmp and physics headers -> "redefinition of getCurrentTestSuite"). The GA library is
# header-only and the tests use textual includes, so __COUNTER__ is monotonic within each
# TU -- restoring it gives every registration a unique name regardless of line layout.
#
# -Wno-c2y-extensions: silences the c2y-extension warning __COUNTER__ raises on clang.
foreach(target ${EXEC_NAME1} ${EXEC_NAME2} ${EXEC_NAME3} ${EXEC_NAME4} ${EXEC_NAME5}
               ${EXEC_NAME6} ${EXEC_NAME7})
    target_compile_definitions(${target} PRIVATE DOCTEST_COUNTER=__COUNTER__)
    target_compile_options(${target} PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-c2y-extensions>)
endforeach()
# target_link_libraries(${EXEC_NAME5} PRIVATE doctest::doctest ga)  # Commented out - STA not yet implemented

# Helper function to link fmt to a target
function(link_fmt_to_target target_name)
    if(TARGET fmt::fmt-header-only)
        target_link_libraries(${target_name} PRIVATE fmt::fmt-header-only)
    elseif(TARGET fmt::fmt)
        target_link_libraries(${target_name} PRIVATE fmt::fmt)
        # On MSVC, fmt needs to be compiled - define FMT_HEADER_ONLY=0 to use compiled version
        if(MSVC)
            target_compile_definitions(${target_name} PRIVATE FMT_HEADER_ONLY=0)
        endif()
    else()
        message(FATAL_ERROR "fmt library not found - dependency system may not be working correctly")
    endif()
endfunction()

# Link fmt to all executables
link_fmt_to_target(${EXEC_NAME1})
link_fmt_to_target(${EXEC_NAME2})
link_fmt_to_target(${EXEC_NAME3})
link_fmt_to_target(${EXEC_NAME4})
link_fmt_to_target(${EXEC_NAME5})
link_fmt_to_target(${EXEC_NAME6})
link_fmt_to_target(${EXEC_NAME7})

# ----------------------------------------------------------------------------
# Standalone benchmarks (utilities/).
#
# These force optimization (MSVC: /O2; others: -O3) regardless of the build
# config, which conflicts with the /RTC1 runtime checks the Debug config injects
# on MSVC (cl error D8016). They therefore live in their own subdirectory so the
# /RTC1 strip needed to make /O2 legal is scoped to the benchmarks alone and does
# not weaken the Debug runtime checks on the test executables above.
# ----------------------------------------------------------------------------
add_subdirectory(utilities)

# ----------------------------------------------------------------------------
# Application overlay (opt-in extension point).
#
# An enclosing build may set GA_PRIVATE_DIR and supply a
# ${GA_PRIVATE_DIR}/ga_test/private_targets.cmake, which is included here and can
# declare additional test/bench targets of its own. Standalone builds leave
# GA_PRIVATE_DIR unset -> nothing extra is added. Placed after
# link_fmt_to_target() (defined above) so the overlay can use it.
# ----------------------------------------------------------------------------
if(GA_PRIVATE_DIR AND EXISTS "${GA_PRIVATE_DIR}/ga_test/private_targets.cmake")
    include("${GA_PRIVATE_DIR}/ga_test/private_targets.cmake")
endif()

# ----------------------------------------------------------------------------
# Cross-check JSON exporter for the Python wrapper (python_utilities/).
#
# Standalone tool (ga + fmt, no doctest), only built with -D_GA_BUILD_PYTHON=ON.
# Lives in its own subdirectory, out of the doctest test targets. Depends on
# link_fmt_to_target (defined above), so this add_subdirectory stays after it.
# ----------------------------------------------------------------------------
add_subdirectory(python_utilities)