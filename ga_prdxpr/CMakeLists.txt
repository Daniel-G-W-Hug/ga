# Copyright 2024-2025, Daniel Hug. All rights reserved.
set(EXEC_NAME ga_prdxpr)
set(EXEC2_NAME ga_prdxpr_viscmp)
set(EXEC3_NAME ga_prdxpr_parser_test)
set(EXEC4_NAME ga_prdxpr_rule_generator_test)

# Source files
set(SOURCES
    src_prdxpr/ga_prdxpr_main.cpp
    src_prdxpr/ga_prdxpr_generator.cpp
    src_prdxpr/ga_prdxpr_ega2d_config.cpp
    src_prdxpr/ga_prdxpr_ega3d_config.cpp
    src_prdxpr/ga_prdxpr_pga2dp_config.cpp
    src_prdxpr/ga_prdxpr_pga3dp_config.cpp
    # Reference mathematical computation files
    src_prdxpr/ga_prdxpr_common.cpp
    # Automatic rule generation system
    src_prdxpr/ga_prdxpr_rule_generator.cpp
    # Rule generation for each algebra
    src_prdxpr/ga_prdxpr_ega2d.cpp
    src_prdxpr/ga_prdxpr_ega3d.cpp
    src_prdxpr/ga_prdxpr_pga2dp.cpp
    src_prdxpr/ga_prdxpr_pga3dp.cpp
)

# Headers
set(HEADERS
    src_prdxpr/ga_prdxpr_config_types.hpp
    src_prdxpr/ga_prdxpr_generator.hpp
    # User configuration files
    src_prdxpr/ga_prdxpr_ega2d_config.hpp
    src_prdxpr/ga_prdxpr_ega3d_config.hpp
    src_prdxpr/ga_prdxpr_pga2dp_config.hpp
    src_prdxpr/ga_prdxpr_pga3dp_config.hpp
    # Reference mathematical headers
    src_prdxpr/ga_prdxpr_common.hpp
    src_prdxpr/ga_prdxpr_ega2d.hpp
    src_prdxpr/ga_prdxpr_ega3d.hpp
    src_prdxpr/ga_prdxpr_pga2dp.hpp
    src_prdxpr/ga_prdxpr_pga3dp.hpp
    src_prdxpr/ga_prdxpr_rule_generator.hpp
)

# Core transformation system
set(SOURCES_VISUAL_COMPARISON_TEST src_trafo/ga_prdxpr_trafo.cpp
            src_trafo/ga_prdxpr_trafo_expr_simplifier.cpp
            src_trafo/ga_prdxpr_trafo_sandwich_transformer.cpp
            src_trafo/ga_prdxpr_trafo_nary_expression.cpp
            src_trafo/ga_prdxpr_trafo_tests.cpp
            src_trafo/ga_prdxpr_trafo_viscmp.cpp)

set(SOURCES_PARSER_TEST src_trafo/ga_prdxpr_trafo_parser_test.cpp
            src_trafo/ga_prdxpr_trafo.cpp)

set(SOURCES_RULE_GENERATOR_TEST src_prdxpr/ga_prdxpr_rule_generator_test.cpp
            src_prdxpr/ga_prdxpr_common.cpp
            src_prdxpr/ga_prdxpr_rule_generator.cpp
            src_prdxpr/ga_prdxpr_ega2d_config.cpp
            src_prdxpr/ga_prdxpr_ega3d_config.cpp
            src_prdxpr/ga_prdxpr_pga2dp_config.cpp
            src_prdxpr/ga_prdxpr_pga3dp_config.cpp
            src_prdxpr/ga_prdxpr_ega2d.cpp
            src_prdxpr/ga_prdxpr_ega3d.cpp
            src_prdxpr/ga_prdxpr_pga2dp.cpp
            src_prdxpr/ga_prdxpr_pga3dp.cpp)

# Create executables
add_executable(${EXEC_NAME} ${SOURCES} ${HEADERS})  #dep: fmt
add_executable(${EXEC2_NAME} ${SOURCES_VISUAL_COMPARISON_TEST})  #dep: fmt
add_executable(${EXEC3_NAME} ${SOURCES_PARSER_TEST})  #dep: fmt, doctest
add_executable(${EXEC4_NAME} ${SOURCES_RULE_GENERATOR_TEST})  #dep: fmt

# Include directories
target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC2_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC3_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC4_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

# Dependencies are handled by the main CMakeLists.txt dependency system
# Just link against the targets that should be available

# Helper function to link fmt to a target
function(link_fmt_to_target target_name)
    if(TARGET fmt::fmt-header-only)
        target_link_libraries(${target_name} PRIVATE fmt::fmt-header-only)
    elseif(TARGET fmt::fmt)
        target_link_libraries(${target_name} PRIVATE fmt::fmt)
    else()
        message(FATAL_ERROR "fmt library not found - dependency system may not be working correctly")
    endif()
endfunction()

# Link fmt to all executables
link_fmt_to_target(${EXEC_NAME})
link_fmt_to_target(${EXEC2_NAME})
link_fmt_to_target(${EXEC4_NAME})

# For EXEC3_NAME, also link doctest
if(TARGET doctest::doctest)
    target_link_libraries(${EXEC3_NAME} PRIVATE doctest::doctest)
    link_fmt_to_target(${EXEC3_NAME})
else()
    message(FATAL_ERROR "doctest not found - dependency system may not be working correctly")
endif()

# Compiler-specific options
if(MSVC)
    target_compile_options(${EXEC_NAME} PRIVATE /bigobj)
endif()