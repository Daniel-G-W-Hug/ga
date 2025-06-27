# incrementally add test files needed

set(EXEC_NAME1 ga_ega_test)
set(EXEC_NAME2 ga_pga_test)
set(EXEC_NAME3 ga_appl_test)

set(EGA_FILES
    src/ga_ega_test.cpp
    src/ga_ega2d_test.cpp 
    src/ga_ega3d_test.cpp
    src/ga_ega4d_test.cpp)

set(PGA_FILES
    src/ga_pga_test.cpp
    src/ga_pga2dp_test.cpp
    src/ga_pga3dp_test.cpp)

    set(APPL_FILES
    src/ga_appl_test.cpp
    src/ga_appl2dp_test.cpp)

add_executable(${EXEC_NAME1} ${EGA_FILES})   #dep: doctest, fmt, ga
add_executable(${EXEC_NAME2} ${PGA_FILES})   #dep: doctest, fmt, ga
add_executable(${EXEC_NAME3} ${APPL_FILES})  #dep: doctest, fmt, ga

target_include_directories(${EXEC_NAME1} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME2} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME3} PRIVATE ${CMAKE_SOURCE_DIR})

find_package(doctest REQUIRED)
find_package(fmt CONFIG REQUIRED)

target_link_libraries(${EXEC_NAME1} PRIVATE doctest::doctest fmt::fmt-header-only ga)
target_link_libraries(${EXEC_NAME2} PRIVATE doctest::doctest fmt::fmt-header-only ga)
target_link_libraries(${EXEC_NAME3} PRIVATE doctest::doctest fmt::fmt-header-only ga)