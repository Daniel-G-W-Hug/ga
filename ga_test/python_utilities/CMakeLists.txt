# Cross-check JSON exporter for the Python wrapper (python_utilities/).
#
# Builds a small executable that walks a representative set of GA operations,
# captures (op, inputs, expected output) tuples, and emits them as JSON to
# ga_py/tests/data/ga_test_cases.json. The Python test suite consumes that file
# to verify the wrapper produces bit-identical results to C++.
#
# Only built when the Python wrapper is enabled (-D_GA_BUILD_PYTHON=ON);
# otherwise the executable would be useless. The companion convenience target
# `regenerate_python_test_data` builds AND runs the exporter, refreshing the
# committed JSON snapshot.
#
# A standalone tool (ga + fmt, no doctest), kept in its own subdirectory rather
# than among the doctest test targets. Depends on link_fmt_to_target (defined in
# the parent), so add_subdirectory(python_utilities) must stay after it.
if(_GA_BUILD_PYTHON)
    set(EXPORT_NAME ga_export_python_cases)

    # Committed JSON snapshot consumed by the Python cross-check test suite.
    set(_GA_PY_CASES_JSON
        "${GA_ROOT}/ga_py/tests/data/ga_test_cases.json")

    add_executable(${EXPORT_NAME} export_python_cases.cpp)
    target_include_directories(${EXPORT_NAME} PRIVATE ${GA_ROOT})
    target_link_libraries(${EXPORT_NAME} PRIVATE ga)
    target_compile_options(${EXPORT_NAME} PRIVATE
        $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-c2y-extensions>)
    # Keep the binary next to the other ga_test artifacts (build/ga_test/).
    set_target_properties(${EXPORT_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/..")
    # Bake in the absolute default output path so the exporter works regardless
    # of the working directory it is launched from (e.g. as an IDE target). A
    # command-line argument still overrides.
    target_compile_definitions(${EXPORT_NAME} PRIVATE
        _GA_PY_CASES_JSON_DEFAULT="${_GA_PY_CASES_JSON}")
    link_fmt_to_target(${EXPORT_NAME})
    add_custom_target(regenerate_python_test_data
        COMMAND ${EXPORT_NAME} "${_GA_PY_CASES_JSON}"
        DEPENDS ${EXPORT_NAME}
        WORKING_DIRECTORY "${GA_ROOT}"
        COMMENT "Regenerating Python cross-check cases at ${_GA_PY_CASES_JSON}"
        VERBATIM
    )
endif()
