# incrementally add test files needed

set(EXEC_NAME ga_lua)

add_executable(${EXEC_NAME} src/ga_lua.cpp)   #dep: lua, sol2, ga

target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})

# Dependencies are handled by the main CMakeLists.txt dependency system
# Just link against the targets that should be available

# Lua is required and should be found by the dependency system
if(NOT Lua_FOUND)
    message(FATAL_ERROR "Lua not found - please install Lua 5.1 or later on your system")
endif()

target_include_directories(${EXEC_NAME} PRIVATE ${LUA_INCLUDE_DIR})
target_link_libraries(${EXEC_NAME} PRIVATE ${LUA_LIBRARIES} ga)

# Link sol2 (should be available from dependency system)
if(TARGET sol2::sol2)
    target_link_libraries(${EXEC_NAME} PRIVATE sol2::sol2)
else()
    message(FATAL_ERROR "sol2 not found - dependency system may not be working correctly")
endif()

# Link fmt (prefer header-only if available, otherwise regular)
if(TARGET fmt::fmt-header-only)
    target_link_libraries(${EXEC_NAME} PRIVATE fmt::fmt-header-only)
elseif(TARGET fmt::fmt)
    target_link_libraries(${EXEC_NAME} PRIVATE fmt::fmt)
else()
    message(FATAL_ERROR "fmt library not found - dependency system may not be working correctly")
endif()

# Optional readline support for better interactive experience
if(GA_USE_READLINE AND TARGET readline::readline)
    target_compile_definitions(${EXEC_NAME} PRIVATE _HD_GA_HAVE_READLINE)
    target_link_libraries(${EXEC_NAME} PRIVATE readline::readline)
    message(STATUS "ga_lua: Using readline for enhanced interactive experience")
elseif(GA_USE_READLINE AND GA_HAS_READLINE)
    # Fallback to raw library if target not created
    target_compile_definitions(${EXEC_NAME} PRIVATE _HD_GA_HAVE_READLINE)
    target_link_libraries(${EXEC_NAME} PRIVATE ${READLINE_LIBRARY})
    message(STATUS "ga_lua: Using readline (fallback method)")
else()
    message(STATUS "ga_lua: readline not available or disabled")
endif()
