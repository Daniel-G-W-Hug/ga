# incrementally add test files needed

set(EXEC_NAME ga_ega_test)
set(EXEC_NAME2 ga_pga_test)

add_executable(${EXEC_NAME} src/ga_ega_test.cpp)  #dep: doctest, fmt, ga
add_executable(${EXEC_NAME2} src/ga_pga_test.cpp)  #dep: doctest, fmt, ga

target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME2} PRIVATE ${CMAKE_SOURCE_DIR})

find_package(doctest REQUIRED)
find_package(fmt CONFIG REQUIRED)

target_link_libraries(${EXEC_NAME} PRIVATE doctest::doctest fmt::fmt-header-only)
target_link_libraries(${EXEC_NAME2} PRIVATE doctest::doctest fmt::fmt-header-only)