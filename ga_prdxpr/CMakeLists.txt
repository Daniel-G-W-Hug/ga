set(EXEC_NAME ga_prdxpr)
set(EXEC2_NAME ga_trafo)  # just temporary for testing

set(SOURCES_PRDXPR src/ga_prdxpr_main.cpp
            src/ga_prdxpr_common.cpp
            src/ga_prdxpr_ega2d.cpp
            src/ga_prdxpr_ega3d.cpp
            src/ga_prdxpr_pga2dp.cpp
            src/ga_prdxpr_pga3dp.cpp)

set(SOURCES_PRDXPR_TRANSFORMATION src/ga_prdxpr_transformation.cpp
            src/ga_prdxpr_transformation_test.cpp)

add_executable(${EXEC_NAME} ${SOURCES_PRDXPR})  #dep: fmt
add_executable(${EXEC2_NAME} ${SOURCES_PRDXPR_TRANSFORMATION})  #dep: fmt, doctest

target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC2_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

find_package(doctest REQUIRED)
find_package(fmt CONFIG REQUIRED)

target_link_libraries(${EXEC_NAME} PRIVATE fmt::fmt-header-only)
target_link_libraries(${EXEC2_NAME} PRIVATE doctest::doctest fmt::fmt-header-only)