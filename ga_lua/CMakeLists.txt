# incrementally add test files needed

set(EXEC_NAME ga_lua)

add_executable(${EXEC_NAME} ga_lua.cpp)   #dep: lua, sol2, ga

target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/../../include/sol2/include)

find_package(Lua 5.1 REQUIRED)
find_package(fmt CONFIG REQUIRED)

target_include_directories(${EXEC_NAME} PRIVATE ${LUA_INCLUDE_DIR})
target_link_libraries(${EXEC_NAME} PRIVATE ${LUA_LIBRARIES} fmt::fmt-header-only)
