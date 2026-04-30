# Micro-benchmarks for ga_py — quantify the Python dispatch tax.
#
# These targets are compiled with the same flags as the rest of the project
# (CMAKE_BUILD_TYPE flows through), so a Debug build of `gpr_overhead_bench`
# pairs with a Debug build of `_ga_py` for a fair head-to-head.

# C++ side of the geometric-product overhead benchmark. The matching Python
# script lives at gpr_overhead_bench.py — same calculation, three different
# call paths (ga_py bindings / pure Python / numpy vectorised).
add_executable(gpr_overhead_bench gpr_overhead_bench.cpp)

target_compile_features(gpr_overhead_bench PRIVATE cxx_std_23)

target_include_directories(gpr_overhead_bench PRIVATE
    "${GA_PROJECT_ROOT}"
)

target_compile_definitions(gpr_overhead_bench PRIVATE
    FMT_HEADER_ONLY=1
)

if(fmt_FOUND)
    target_link_libraries(gpr_overhead_bench PRIVATE fmt::fmt-header-only)
elseif(EXISTS "/opt/homebrew/include/fmt/format.h")
    target_include_directories(gpr_overhead_bench PRIVATE "/opt/homebrew/include")
endif()
