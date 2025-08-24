# Cross-platform dependency management for GA project
# Handles flexible dependency resolution between system packages and FetchContent

# Options for dependency management
option(GA_USE_SYSTEM_DEPENDENCIES "Use system-installed dependencies when available" ON)
option(GA_FORCE_FETCH_CONTENT "Force use of FetchContent for supported dependencies" OFF)

# Include platform-specific dependency handlers
include(${CMAKE_CURRENT_LIST_DIR}/find_dependencies.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/fetch_dependencies.cmake)

# Main dependency resolution function
function(setup_ga_dependencies)
    message(STATUS "Setting up GA dependencies...")
    
    # Check for local hd utility library (project-specific dependency)
    if(EXISTS "${CMAKE_SOURCE_DIR}/../../include/hd")
        message(STATUS "✓ Found local hd utility library at: ${CMAKE_SOURCE_DIR}/../../include/hd")
        set(GA_HAS_HD_UTILS TRUE)
    else()
        message(WARNING "hd utility library not found at expected location: ${CMAKE_SOURCE_DIR}/../../include/hd")
        message(STATUS "Some ga_view functionality may not work properly")
        set(GA_HAS_HD_UTILS FALSE)
    endif()
    
    # Find flexible dependencies using three-tier logic
    find_fmt_tiered()
    find_doctest_tiered() 
    find_sol2_tiered()
    
    # Find required system dependencies (Qt6, Lua) - must be system installed
    find_package(Qt6 COMPONENTS Core Gui Widgets QUIET)
    if(Qt6_FOUND)
        message(STATUS "✓ Found Qt6: ${Qt6_VERSION}")
    else()
        message(WARNING "Qt6 not found! Please install:")
        if(APPLE)
            message(STATUS "  macOS: brew install qt6")
        elseif(WIN32)
            message(STATUS "  Windows: Download from https://www.qt.io/download")
        else()
            message(STATUS "  Linux: sudo apt-get install qt6-base-dev qt6-tools-dev")
        endif()
    endif()
    
    find_package(Lua 5.1 QUIET)
    if(Lua_FOUND OR LUA_FOUND)
        message(STATUS "✓ Found Lua: ${LUA_VERSION_STRING}")
        # Make sure Lua variables are available in parent scope (for modules)
        set(Lua_FOUND ${Lua_FOUND} PARENT_SCOPE)
        set(LUA_FOUND ${LUA_FOUND} PARENT_SCOPE) 
        set(LUA_LIBRARIES ${LUA_LIBRARIES} PARENT_SCOPE)
        set(LUA_INCLUDE_DIR ${LUA_INCLUDE_DIR} PARENT_SCOPE)
        set(LUA_VERSION_STRING ${LUA_VERSION_STRING} PARENT_SCOPE)
        set(GA_HAS_LUA TRUE PARENT_SCOPE)
    else()
        message(WARNING "Lua not found! Please install:")
        if(APPLE)
            message(STATUS "  macOS: brew install lua")
        elseif(WIN32)
            message(STATUS "  Windows: Download from https://www.lua.org/download.html")
        else()
            message(STATUS "  Linux: sudo apt-get install lua5.1-dev")
        endif()
        set(GA_HAS_LUA FALSE PARENT_SCOPE)
    endif()
    
    # Setup optional dependencies using three-tier logic
    find_readline_tiered()
    
    # Handle flexible dependencies based on platform and options
    if(GA_FORCE_FETCH_CONTENT)
        message(STATUS "Forcing FetchContent for all supported dependencies")
        setup_all_fetchcontent_dependencies()
    elseif(GA_USE_SYSTEM_DEPENDENCIES)
        message(STATUS "Using system dependencies where available, FetchContent for missing ones")
        setup_missing_dependencies()
    else()
        message(STATUS "Using FetchContent for all supported dependencies")
        setup_all_fetchcontent_dependencies()
    endif()
    
endfunction()

