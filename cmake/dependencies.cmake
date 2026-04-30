# Cross-platform dependency management for GA project
# Handles flexible dependency resolution between system packages and FetchContent

# Options for dependency management
option(GA_USE_SYSTEM_DEPENDENCIES "Use system-installed dependencies when available" ON)
option(GA_FORCE_FETCH_CONTENT "Force use of FetchContent for supported dependencies" OFF)

# Include platform-specific dependency handlers
include(${CMAKE_CURRENT_LIST_DIR}/find_dependencies.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/fetch_dependencies.cmake)

# Detect vcpkg usage and include appropriate handler
if(WIN32 AND DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg")
    message(STATUS "Detected vcpkg toolchain on Windows")
    include(${CMAKE_CURRENT_LIST_DIR}/vcpkg_dependencies.cmake)
endif()

# Main dependency resolution function
function(setup_ga_dependencies)
    # Use vcpkg-based dependency management if vcpkg toolchain detected
    if(WIN32 AND DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg")
        setup_vcpkg_dependencies()
        return()
    endif()

    # Otherwise use traditional dependency management
    message(STATUS "Setting up GA dependencies...")

    # Find flexible dependencies using three-tier logic
    find_fmt_tiered()
    find_doctest_tiered()

    # sol2, Lua, readline are only needed when ga_lua is built. Skipping the
    # discovery (and the FetchContent fallbacks down below) when _GA_USE_LUA
    # is OFF means a Lua-less build doesn't drag in sol2, doesn't nag about
    # missing Lua, and doesn't hunt for readline.
    if(_GA_USE_LUA)
        find_sol2_tiered()

        # On macOS, prefer Lua 5.4 (sol2 is not compatible with Lua 5.5+)
        if(APPLE)
            set(LUA54_BREW_PREFIX "/opt/homebrew/opt/lua@5.4")
            if(EXISTS "${LUA54_BREW_PREFIX}")
                message(STATUS "ga_lua: Preferring Lua 5.4 from Homebrew (sol2 not compatible with Lua 5.5+)")
                list(INSERT CMAKE_PREFIX_PATH 0 "${LUA54_BREW_PREFIX}")
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
                message(STATUS "  macOS: brew install lua@5.4")
            elseif(WIN32)
                message(STATUS "  Windows: Download from https://www.lua.org/download.html")
            else()
                message(STATUS "  Linux: sudo apt-get install lua5.1-dev")
            endif()
            set(GA_HAS_LUA FALSE PARENT_SCOPE)
        endif()

        # Setup optional dependencies using three-tier logic
        find_readline_tiered()
    else()
        message(STATUS "Skipping sol2 / Lua / readline lookup (_GA_USE_LUA=OFF)")
        # Mark as "satisfied" so setup_missing_dependencies below doesn't try
        # to FetchContent-fetch sol2 — we don't need it for a Lua-less build.
        set(GA_HAS_SOL2 TRUE PARENT_SCOPE)
        set(GA_HAS_LUA FALSE PARENT_SCOPE)
        set(GA_HAS_READLINE FALSE PARENT_SCOPE)
    endif()

    # Find required system dependencies (Qt6) - must be system installed
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

