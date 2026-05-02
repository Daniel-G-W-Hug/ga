# Copyright 2024-2025, Daniel Hug. All rights reserved.
set(EXEC_NAME ga_prdxpr)
set(EXEC2_NAME ga_prdxpr_viscmp)
set(EXEC3_NAME ga_prdxpr_parser_test)
set(EXEC4_NAME ga_prdxpr_rule_generator_test)

# Source files
set(SOURCES
    # Entry point
    src_prdxpr/ga_prdxpr_main.cpp
    # CLI option parsing
    src_prdxpr/cli/ga_prdxpr_options.cpp
    # Foundational utilities and shared types
    src_prdxpr/core/ga_prdxpr_common.cpp
    # Automatic rule generation engine + metric/dual math
    src_prdxpr/rules/ga_prdxpr_rule_generator.cpp
    src_prdxpr/rules/ga_prdxpr_metric_calc.cpp
    src_prdxpr/rules/ga_prdxpr_dual_calc.cpp
    src_prdxpr/rules/ga_prdxpr_metric_export.cpp
    # Per-algebra: *_config.cpp MUST come before *.cpp.
    # The *.cpp TUs initialize static rules by calling get_<algebra>_algebra_config(),
    # which reads namespace-scope const globals defined in headers (internal linkage,
    # one copy per TU). Listing the config TU first preserves the dynamic-initialization
    # order across translation units that the original layout depended on.
    src_prdxpr/algebras/ga_prdxpr_ega2d_config.cpp
    src_prdxpr/algebras/ga_prdxpr_ega3d_config.cpp
    src_prdxpr/algebras/ga_prdxpr_pga2dp_config.cpp
    src_prdxpr/algebras/ga_prdxpr_pga3dp_config.cpp
    src_prdxpr/algebras/ga_prdxpr_sta4d_config.cpp
    src_prdxpr/algebras/ga_prdxpr_ega2d.cpp
    src_prdxpr/algebras/ga_prdxpr_ega3d.cpp
    src_prdxpr/algebras/ga_prdxpr_pga2dp.cpp
    src_prdxpr/algebras/ga_prdxpr_pga3dp.cpp
    src_prdxpr/algebras/ga_prdxpr_sta4d.cpp
    # Sandwich product transformation pipeline
    src_prdxpr/sandwich/ga_prdxpr_sandwich_parser.cpp
    src_prdxpr/sandwich/ga_prdxpr_sandwich_simplifier.cpp
    src_prdxpr/sandwich/ga_prdxpr_sandwich_nary_expression.cpp
    src_prdxpr/sandwich/ga_prdxpr_sandwich_transformer.cpp
    # Generator orchestrator
    src_prdxpr/generator/ga_prdxpr_generator.cpp
    # Code generation system (--output=code)
    src_prdxpr/codegen/ga_codegen_types.cpp
    src_prdxpr/codegen/ga_codegen_emitter.cpp
)

# Headers
set(HEADERS
    # CLI
    src_prdxpr/cli/ga_prdxpr_options.hpp
    # Foundational
    src_prdxpr/core/ga_prdxpr_common.hpp
    src_prdxpr/core/ga_prdxpr_config_types.hpp
    # Rule generation engine + metric/dual math
    src_prdxpr/rules/ga_prdxpr_rule_generator.hpp
    src_prdxpr/rules/ga_prdxpr_metric_calc.hpp
    src_prdxpr/rules/ga_prdxpr_dual_calc.hpp
    src_prdxpr/rules/ga_prdxpr_metric_export.hpp
    # Per-algebra
    src_prdxpr/algebras/ga_prdxpr_ega2d.hpp
    src_prdxpr/algebras/ga_prdxpr_ega3d.hpp
    src_prdxpr/algebras/ga_prdxpr_pga2dp.hpp
    src_prdxpr/algebras/ga_prdxpr_pga3dp.hpp
    src_prdxpr/algebras/ga_prdxpr_sta4d.hpp
    src_prdxpr/algebras/ga_prdxpr_ega2d_config.hpp
    src_prdxpr/algebras/ga_prdxpr_ega3d_config.hpp
    src_prdxpr/algebras/ga_prdxpr_pga2dp_config.hpp
    src_prdxpr/algebras/ga_prdxpr_pga3dp_config.hpp
    src_prdxpr/algebras/ga_prdxpr_sta4d_config.hpp
    # Sandwich pipeline
    src_prdxpr/sandwich/ga_prdxpr_sandwich_parser.hpp
    src_prdxpr/sandwich/ga_prdxpr_sandwich_simplifier.hpp
    src_prdxpr/sandwich/ga_prdxpr_sandwich_nary_expression.hpp
    src_prdxpr/sandwich/ga_prdxpr_sandwich_transformer.hpp
    # Generator
    src_prdxpr/generator/ga_prdxpr_generator.hpp
    # Code generation system
    src_prdxpr/codegen/ga_codegen_types.hpp
    src_prdxpr/codegen/ga_codegen_emitter.hpp
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

set(SOURCES_RULE_GENERATOR_TEST src_prdxpr/tests/ga_prdxpr_rule_generator_test.cpp
            src_prdxpr/core/ga_prdxpr_common.cpp
            src_prdxpr/rules/ga_prdxpr_rule_generator.cpp
            src_prdxpr/rules/ga_prdxpr_metric_calc.cpp
            src_prdxpr/rules/ga_prdxpr_dual_calc.cpp
            src_prdxpr/algebras/ga_prdxpr_ega2d_config.cpp
            src_prdxpr/algebras/ga_prdxpr_ega3d_config.cpp
            src_prdxpr/algebras/ga_prdxpr_pga2dp_config.cpp
            src_prdxpr/algebras/ga_prdxpr_pga3dp_config.cpp
            src_prdxpr/algebras/ga_prdxpr_sta4d_config.cpp
            src_prdxpr/algebras/ga_prdxpr_ega2d.cpp
            src_prdxpr/algebras/ga_prdxpr_ega3d.cpp
            src_prdxpr/algebras/ga_prdxpr_pga2dp.cpp
            src_prdxpr/algebras/ga_prdxpr_pga3dp.cpp
            src_prdxpr/algebras/ga_prdxpr_sta4d.cpp)

# Create executables
add_executable(${EXEC_NAME} ${SOURCES} ${HEADERS})  #dep: fmt
add_executable(${EXEC2_NAME} ${SOURCES_VISUAL_COMPARISON_TEST})  #dep: fmt
add_executable(${EXEC3_NAME} ${SOURCES_PARSER_TEST})  #dep: fmt, doctest
add_executable(${EXEC4_NAME} ${SOURCES_RULE_GENERATOR_TEST})  #dep: fmt

# Include directories.
# Project root is needed for "ga/..." headers; src_prdxpr/ is needed for the
# subfolder-prefixed includes (e.g. "core/ga_prdxpr_common.hpp").
target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR}
                                                ${CMAKE_CURRENT_SOURCE_DIR}/src_prdxpr)
target_include_directories(${EXEC2_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC3_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC4_NAME} PRIVATE ${CMAKE_SOURCE_DIR}
                                                 ${CMAKE_CURRENT_SOURCE_DIR}/src_prdxpr)

# Dependencies are handled by the main CMakeLists.txt dependency system
# Just link against the targets that should be available

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

# mdspan configuration - enable square bracket operator
target_compile_definitions(${EXEC_NAME} PUBLIC MDSPAN_USE_BRACKET_OPERATOR=1 MDSPAN_CXX_STANDARD=23)
target_compile_definitions(${EXEC2_NAME} PUBLIC MDSPAN_USE_BRACKET_OPERATOR=1 MDSPAN_CXX_STANDARD=23)
target_compile_definitions(${EXEC3_NAME} PUBLIC MDSPAN_USE_BRACKET_OPERATOR=1 MDSPAN_CXX_STANDARD=23)
target_compile_definitions(${EXEC4_NAME} PUBLIC MDSPAN_USE_BRACKET_OPERATOR=1 MDSPAN_CXX_STANDARD=23)

# Suppress -Wc2y-extensions triggered by __COUNTER__ in doctest macros (third-party header)
target_compile_options(${EXEC3_NAME} PRIVATE $<$<CXX_COMPILER_ID:Clang,AppleClang>:-Wno-c2y-extensions>)

# Compiler-specific options
if(MSVC)
    target_compile_options(${EXEC_NAME} PRIVATE /bigobj)
endif()