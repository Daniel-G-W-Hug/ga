# incrementally add test files needed

set(EXEC_NAME1 ga_ega_test)
set(EXEC_NAME2 ga_pga_test)
set(EXEC_NAME3 ga_appl2dp_test)
set(EXEC_NAME4 ga_appl3dp_test)
set(EXPORT_NAME ga_export_python_cases)
# set(EXEC_NAME5 ga_sta_test)  # Commented out - STA not yet implemented

set(EGA_FILES
    src/ga_ega_test.cpp)

set(PGA_FILES
    src/ga_pga_test.cpp)

set(APPL2DP_FILES
    src/ga_appl2dp_test.cpp)

set(APPL3DP_FILES
    src/ga_appl3dp_test.cpp)

# set(STA_FILES
#     src/ga_sta_test.cpp)  # Commented out - STA not yet implemented

add_executable(${EXEC_NAME1} ${EGA_FILES})      #dep: doctest, fmt, ga
add_executable(${EXEC_NAME2} ${PGA_FILES})      #dep: doctest, fmt, ga
add_executable(${EXEC_NAME3} ${APPL2DP_FILES})  #dep: doctest, fmt, ga
add_executable(${EXEC_NAME4} ${APPL3DP_FILES})  #dep: doctest, fmt, ga
# add_executable(${EXEC_NAME5} ${STA_FILES})    #dep: doctest, fmt, ga  # Commented out - STA not yet implemented

target_include_directories(${EXEC_NAME1} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME2} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME3} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME4} PRIVATE ${CMAKE_SOURCE_DIR})
# target_include_directories(${EXEC_NAME5} PRIVATE ${CMAKE_SOURCE_DIR})  # Commented out - STA not yet implemented

# Dependencies are handled by the main CMakeLists.txt dependency system
# Just link against the targets that should be available

target_link_libraries(${EXEC_NAME1} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME2} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME3} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME4} PRIVATE doctest::doctest ga)

# Suppress -Wc2y-extensions triggered by __COUNTER__ in doctest macros (third-party header)
foreach(target ${EXEC_NAME1} ${EXEC_NAME2} ${EXEC_NAME3} ${EXEC_NAME4})
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
# link_fmt_to_target(${EXEC_NAME5})  # Commented out - STA not yet implemented

# ----------------------------------------------------------------------------
# Cross-check JSON exporter for the Python wrapper.
#
# Builds a small executable that walks a representative set of GA operations,
# captures (op, inputs, expected output) tuples, and emits them as JSON to
# ga_py/tests/data/ga_test_cases.json. The Python test suite consumes that
# file to verify the wrapper produces bit-identical results to C++.
#
# Only built when the Python wrapper is enabled (-D_GA_BUILD_PYTHON=ON);
# otherwise the executable would be useless. The companion convenience target
# `regenerate_python_test_data` builds AND runs the exporter, refreshing the
# committed JSON snapshot.
#
# This target does NOT depend on doctest --- it's a standalone program that
# uses only ga + fmt.
# ----------------------------------------------------------------------------
if(_GA_BUILD_PYTHON)
    add_executable(${EXPORT_NAME} src/export_python_cases.cpp)
    target_include_directories(${EXPORT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
    target_link_libraries(${EXPORT_NAME} PRIVATE ga)
    target_compile_options(${EXPORT_NAME} PRIVATE
        $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-c2y-extensions>)
    link_fmt_to_target(${EXPORT_NAME})

    # Convenience: refresh the committed JSON snapshot in one step.
    set(_GA_PY_CASES_JSON
        "${CMAKE_SOURCE_DIR}/ga_py/tests/data/ga_test_cases.json")
    add_custom_target(regenerate_python_test_data
        COMMAND ${EXPORT_NAME} "${_GA_PY_CASES_JSON}"
        DEPENDS ${EXPORT_NAME}
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        COMMENT "Regenerating Python cross-check cases at ${_GA_PY_CASES_JSON}"
        VERBATIM
    )
endif()