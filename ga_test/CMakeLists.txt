# incrementally add test files needed

set(EXEC_NAME ga_test)

add_executable(${EXEC_NAME} ga_test.cpp)  #dep: doctest, fmt

target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

find_package(doctest REQUIRED)
find_package(fmt CONFIG REQUIRED)

target_link_libraries(${EXEC_NAME} PRIVATE doctest::doctest fmt::fmt-header-only)