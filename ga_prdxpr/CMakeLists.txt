set(EXEC_NAME ga_prdxpr)

set(SOURCES src/ga_prdxpr_main.cpp
            src/ga_prdxpr_common.cpp
            src/ga_prd_ega2d.cpp
            src/ga_prd_ega3d.cpp
            )

add_executable(${EXEC_NAME} ${SOURCES})  #dep: fmt

target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

find_package(fmt CONFIG REQUIRED)

target_link_libraries(${EXEC_NAME} PRIVATE fmt::fmt-header-only)