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
target_include_directories(${BENCH_TRANSFORM} PRIVATE ${CMAKE_SOURCE_DIR})
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

set(BENCH_DIPOLE ga_dipole_bench)
add_executable(${BENCH_DIPOLE} bench_dipole_ga_vs_classical.cpp)
target_include_directories(${BENCH_DIPOLE} PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(${BENCH_DIPOLE} PRIVATE ga)
link_fmt_to_target(${BENCH_DIPOLE})
set_target_properties(${BENCH_DIPOLE} PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    RUNTIME_OUTPUT_DIRECTORY "${_BENCH_OUTPUT_DIR}")
target_compile_definitions(${BENCH_DIPOLE} PRIVATE NDEBUG)
if(MSVC)
    target_compile_options(${BENCH_DIPOLE} PRIVATE /O2)
else()
    target_compile_options(${BENCH_DIPOLE} PRIVATE -O3)
endif()

add_custom_target(run_${BENCH_DIPOLE}
    COMMAND ${BENCH_DIPOLE}
    DEPENDS ${BENCH_DIPOLE}
    WORKING_DIRECTORY "${_BENCH_OUTPUT_DIR}"
    COMMENT "Running dipole GA-vs-classical benchmark"
    VERBATIM
)

# Maglev motion stability-envelope study (closed-loop sweeps; on-demand only).
set(MAGLEV_STABILITY ga_maglev_stability)
add_executable(${MAGLEV_STABILITY} bench_maglev_stability.cpp)
target_include_directories(${MAGLEV_STABILITY} PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(${MAGLEV_STABILITY} PRIVATE ga)
link_fmt_to_target(${MAGLEV_STABILITY})
set_target_properties(${MAGLEV_STABILITY} PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    RUNTIME_OUTPUT_DIRECTORY "${_BENCH_OUTPUT_DIR}")
target_compile_definitions(${MAGLEV_STABILITY} PRIVATE NDEBUG)
if(MSVC)
    target_compile_options(${MAGLEV_STABILITY} PRIVATE /O2)
else()
    target_compile_options(${MAGLEV_STABILITY} PRIVATE -O3)
endif()

add_custom_target(run_${MAGLEV_STABILITY}
    COMMAND ${MAGLEV_STABILITY}
    DEPENDS ${MAGLEV_STABILITY}
    WORKING_DIRECTORY "${_BENCH_OUTPUT_DIR}"
    COMMENT "Running maglev motion stability-envelope study"
    VERBATIM
)

# Headless replica of the ga_view maglev demo state machine: validate choreography /
# controller changes before shipping them to the scene (on-demand only).
set(MAGLEV_REPLICA maglev_scene_replica)
add_executable(${MAGLEV_REPLICA} maglev_scene_replica.cpp)
target_include_directories(${MAGLEV_REPLICA} PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(${MAGLEV_REPLICA} PRIVATE ga)
link_fmt_to_target(${MAGLEV_REPLICA})
set_target_properties(${MAGLEV_REPLICA} PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    RUNTIME_OUTPUT_DIRECTORY "${_BENCH_OUTPUT_DIR}")
target_compile_definitions(${MAGLEV_REPLICA} PRIVATE NDEBUG)
if(MSVC)
    target_compile_options(${MAGLEV_REPLICA} PRIVATE /O2)
else()
    target_compile_options(${MAGLEV_REPLICA} PRIVATE -O3)
endif()

add_custom_target(run_${MAGLEV_REPLICA}
    COMMAND ${MAGLEV_REPLICA}
    DEPENDS ${MAGLEV_REPLICA}
    WORKING_DIRECTORY "${_BENCH_OUTPUT_DIR}"
    COMMENT "Running maglev scene replica validation"
    VERBATIM
)
