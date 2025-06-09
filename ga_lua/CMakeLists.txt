# incrementally add test files needed

set(EXEC_NAME ga_lua)

add_executable(${EXEC_NAME} src/ga_lua.cpp)   #dep: lua, sol2, ga

target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/../../include/sol2/include)

find_package(Lua 5.1 REQUIRED)
find_package(fmt CONFIG REQUIRED)

# Optional readline support for better interactive experience
option(GA_USE_READLINE "Enable readline support for ga_lua interactive shell" ON)

if(GA_USE_READLINE)
    find_library(READLINE_LIBRARY readline)
    if(READLINE_LIBRARY)
        target_compile_definitions(${EXEC_NAME} PRIVATE _HD_GA_HAVE_READLINE)
        target_link_libraries(${EXEC_NAME} PRIVATE ${READLINE_LIBRARY})
        message(STATUS "ga_lua: Using readline for enhanced interactive experience")
    else()
        message(STATUS "ga_lua: readline not found, using basic input (install with: brew install readline / apt-get install libreadline-dev)")
    endif()
else()
    message(STATUS "ga_lua: readline support disabled via GA_USE_READLINE=OFF")
endif()

target_include_directories(${EXEC_NAME} PRIVATE ${LUA_INCLUDE_DIR})
target_link_libraries(${EXEC_NAME} PRIVATE ${LUA_LIBRARIES} fmt::fmt-header-only ga)
