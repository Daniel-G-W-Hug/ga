set(EXEC_NAME ga_prdxpr)
set(EXEC2_NAME trafo-s)  # just temporary for testing

set(SOURCES src/ga_prdxpr_main.cpp
            src/ga_prdxpr_common.cpp
            src/ga_prdxpr_transformer.cpp
            src/ga_prd_ega2d.cpp
            src/ga_prd_ega3d.cpp
            src/ga_prd_pga2dp.cpp
            src/ga_prd_pga3dp.cpp
            )

add_executable(${EXEC_NAME} ${SOURCES})  #dep: fmt
add_executable(${EXEC2_NAME} src/ga_prdxpr_transformer-simple.cpp)  #dep: fmt

target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC2_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

find_package(fmt CONFIG REQUIRED)

target_link_libraries(${EXEC_NAME} PRIVATE fmt::fmt-header-only)
target_link_libraries(${EXEC2_NAME} PRIVATE fmt::fmt-header-only)