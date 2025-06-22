# Copyright 2024-2025, Daniel Hug. All rights reserved.
set(EXEC_NAME ga_prdxpr)
set(EXEC2_NAME ga_prdxpr_viscmp)
set(EXEC3_NAME ga_prdxpr_parser_test)

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
)

# Core transformation system
set(SOURCES_VISUAL_COMPARISON_TEST src_trafo/ga_prdxpr_trafo.cpp
            src_trafo/ga_prdxpr_trafo_expr_simplifier.cpp
            src_trafo/ga_prdxpr_trafo_sandwich_transformer.cpp
            src_trafo/ga_prdxpr_trafo_nary_expression.cpp
            src_trafo/ga_prdxpr_trafo_viscmp.cpp)

set(SOURCES_PARSER_TEST src_trafo/ga_prdxpr_trafo_parser_test.cpp
            src_trafo/ga_prdxpr_trafo.cpp)

# Create executables
add_executable(${EXEC_NAME} ${SOURCES} ${HEADERS})  #dep: fmt
add_executable(${EXEC2_NAME} ${SOURCES_VISUAL_COMPARISON_TEST})  #dep: fmt
add_executable(${EXEC3_NAME} ${SOURCES_PARSER_TEST})  #dep: fmt, doctest

# Include directories
target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC2_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC3_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

find_package(fmt CONFIG REQUIRED)
find_package(doctest REQUIRED)

# Link with fmt library
target_link_libraries(${EXEC_NAME} PRIVATE fmt::fmt-header-only)
target_link_libraries(${EXEC2_NAME} PRIVATE fmt::fmt-header-only)
target_link_libraries(${EXEC3_NAME} PRIVATE fmt::fmt-header-only doctest::doctest)

# Compiler-specific options
if(MSVC)
    target_compile_options(${EXEC_NAME} PRIVATE /bigobj)
endif()