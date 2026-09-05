# Standalone GA benchmarks.
#
# These are NOT part of any test run --- they are built and run on demand only.
# Each is a plain executable (ga + fmt, no doctest). Build a single one with e.g.
#     cmake --build . --target ga_sta_bench_transform
# then run it directly, or use the convenience `run_<target>` custom target which
# builds AND runs it.
#
# Benchmarks are forced to -O3/NDEBUG (MSVC: /O2) regardless of CMAKE_BUILD_TYPE,
# because the relative cost of the measured variants depends strongly on the
# optimization level and a Debug build would be meaningless.
#
# MSVC gotcha: the Debug config injects /RTC1 (runtime checks), and /O2 + /RTC1
# is a hard cl error (D8016). /RTC1 cannot be removed per target, only per scope,
# so this CMakeLists lives in its own subdirectory and strips /RTC1 from the
# directory-local Debug flags below. The strip applies only to targets defined
# here --- the test executables in the parent directory keep their /RTC1.
if(MSVC)
    string(REGEX REPLACE "/RTC[1csu]+" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endif()

# Keep benchmark binaries next to the other ga_test artifacts (build/ga_test/)
# rather than build/ga_test/utilities/, so existing run/launch paths are unchanged.
set(_BENCH_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/..")

set(BENCH_TRANSFORM ga_sta_bench_transform)
add_executable(${BENCH_TRANSFORM} bench_sta4ds_transform.cpp)
target_include_directories(${BENCH_TRANSFORM} PRIVATE ${GA_ROOT})
target_link_libraries(${BENCH_TRANSFORM} PRIVATE ga)
link_fmt_to_target(${BENCH_TRANSFORM})
# exclude from the default "build all" so a normal build/test cycle never builds it
set_target_properties(${BENCH_TRANSFORM} PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    RUNTIME_OUTPUT_DIRECTORY "${_BENCH_OUTPUT_DIR}")
target_compile_definitions(${BENCH_TRANSFORM} PRIVATE NDEBUG)
if(MSVC)
    target_compile_options(${BENCH_TRANSFORM} PRIVATE /O2)
else()
    target_compile_options(${BENCH_TRANSFORM} PRIVATE -O3)
endif()

# Convenience: build AND run the benchmark in one step.
add_custom_target(run_${BENCH_TRANSFORM}
    COMMAND ${BENCH_TRANSFORM}
    DEPENDS ${BENCH_TRANSFORM}
    WORKING_DIRECTORY "${_BENCH_OUTPUT_DIR}"
    COMMENT "Running sta4ds transform benchmark"
    VERBATIM
)

# The division-guard bench PAIR: one source, two binaries, the second through
# ga_division_guard(... OFF). What the guard costs on this compiler, as a ratio of two
# ns/op figures -- see the comment block at the guard's definition in ga/CMakeLists.txt.
foreach(_arm "" "_off")
    set(_t ga_bench_division_guard${_arm})
    add_executable(${_t} bench_division_guard.cpp)
    target_include_directories(${_t} PRIVATE ${GA_ROOT})
    target_link_libraries(${_t} PRIVATE ga)
    link_fmt_to_target(${_t})
    set_target_properties(${_t} PROPERTIES
        EXCLUDE_FROM_ALL TRUE
        RUNTIME_OUTPUT_DIRECTORY "${_BENCH_OUTPUT_DIR}")
    target_compile_definitions(${_t} PRIVATE NDEBUG)
    if(MSVC)
        target_compile_options(${_t} PRIVATE /O2)
    else()
        target_compile_options(${_t} PRIVATE -O3)
    endif()
endforeach()
ga_division_guard(ga_bench_division_guard_off OFF)
