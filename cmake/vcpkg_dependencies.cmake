# vcpkg-based dependency management for Windows
# This file handles dependencies when using vcpkg toolchain

function(setup_vcpkg_dependencies)
    message(STATUS "Setting up dependencies using vcpkg...")

    # Required dependencies - these should always be available via vcpkg
    find_package(fmt CONFIG REQUIRED)
    if(fmt_FOUND)
        message(STATUS "✓ Found vcpkg fmt: ${fmt_VERSION}")
        set(GA_HAS_FMT TRUE PARENT_SCOPE)
    else()
        message(FATAL_ERROR "fmt not found via vcpkg - run: vcpkg install fmt")
    endif()

    find_package(doctest CONFIG REQUIRED)
    if(doctest_FOUND)
        message(STATUS "✓ Found vcpkg doctest")
        set(GA_HAS_DOCTEST TRUE PARENT_SCOPE)
    else()
        message(FATAL_ERROR "doctest not found via vcpkg - run: vcpkg install doctest")
    endif()

    # sol2 is header-only from vcpkg
    find_path(SOL2_INCLUDE_DIRS "sol/sol.hpp"
        PATHS ${CMAKE_BINARY_DIR}/vcpkg_installed/${VCPKG_TARGET_TRIPLET}/include
        NO_DEFAULT_PATH)
    if(SOL2_INCLUDE_DIRS)
        message(STATUS "✓ Found vcpkg sol2 headers at: ${SOL2_INCLUDE_DIRS}")
        set(GA_HAS_SOL2 TRUE PARENT_SCOPE)
        set(SOL2_INCLUDE_DIRS ${SOL2_INCLUDE_DIRS} PARENT_SCOPE)
    else()
        message(FATAL_ERROR "sol2 not found via vcpkg - run: vcpkg install sol2")
    endif()

    # Use ONLY system-installed Lua with static linking - no vcpkg conflicts
    message(STATUS "Using system-installed Lua with static linking (no vcpkg)")

    # Direct static linking to system Lua installation
    set(LUA_INCLUDE_DIR "C:/Users/$ENV{USERNAME}/AppData/Local/Programs/Lua/include")
    set(LUA_LIBRARY "C:/Users/$ENV{USERNAME}/AppData/Local/Programs/Lua/lib/lua54.lib")

    # Verify the files exist
    if(EXISTS "${LUA_INCLUDE_DIR}/lua.h" AND EXISTS "${LUA_LIBRARY}")
        message(STATUS "✓ Found static system lua: headers at ${LUA_INCLUDE_DIR}, lib at ${LUA_LIBRARY}")
        set(LUA_LIBRARIES ${LUA_LIBRARY})
        set(LUA_FOUND TRUE)
        set(LUA_VERSION_STRING "5.4.6")
        set(GA_HAS_LUA TRUE PARENT_SCOPE)
        # Make variables available in parent scope AND cache for subdirectories
        set(LUA_FOUND ${LUA_FOUND} PARENT_SCOPE)
        set(LUA_LIBRARIES ${LUA_LIBRARIES} PARENT_SCOPE)
        set(LUA_INCLUDE_DIR ${LUA_INCLUDE_DIR} PARENT_SCOPE)
        set(LUA_VERSION_STRING ${LUA_VERSION_STRING} PARENT_SCOPE)
        # Cache variables for subdirectories
        set(LUA_FOUND ${LUA_FOUND} CACHE INTERNAL "")
        set(LUA_LIBRARIES ${LUA_LIBRARIES} CACHE INTERNAL "")
        set(LUA_INCLUDE_DIR ${LUA_INCLUDE_DIR} CACHE INTERNAL "")
        set(LUA_VERSION_STRING ${LUA_VERSION_STRING} CACHE INTERNAL "")
        set(GA_HAS_LUA TRUE CACHE INTERNAL "")
    else()
        message(FATAL_ERROR "System Lua not found at expected location: ${LUA_INCLUDE_DIR} / ${LUA_LIBRARY}")
        message(STATUS "Please ensure Lua 5.4 is installed at C:/Users/$ENV{USERNAME}/AppData/Local/Programs/Lua/")
        set(GA_HAS_LUA FALSE PARENT_SCOPE)
    endif()

    # Check for local hd utility library (project-specific dependency)
    if(EXISTS "${CMAKE_SOURCE_DIR}/../../include/hd")
        message(STATUS "✓ Found local hd utility library at: ${CMAKE_SOURCE_DIR}/../../include/hd")
        set(GA_HAS_HD_UTILS TRUE PARENT_SCOPE)
    else()
        message(WARNING "hd utility library not found at expected location: ${CMAKE_SOURCE_DIR}/../../include/hd")
        message(STATUS "Some ga_view functionality may not work properly")
        set(GA_HAS_HD_UTILS FALSE PARENT_SCOPE)
    endif()

    # Optional readline support (Windows-compatible implementation)
    if(GA_USE_READLINE)
        find_package(readline-win32 CONFIG QUIET)
        if(readline-win32_FOUND)
            message(STATUS "✓ Found vcpkg readline-win32")
            # Create readline::readline target for compatibility with existing code
            if(NOT TARGET readline::readline)
                add_library(readline::readline ALIAS readline-win32::readline-win32)
            endif()
            set(GA_HAS_READLINE TRUE PARENT_SCOPE)
        else()
            message(STATUS "✗ readline-win32 not found")
            message(STATUS "  To enable readline support, run: vcpkg install readline-win32")
            set(GA_HAS_READLINE FALSE PARENT_SCOPE)
        endif()
    else()
        message(STATUS "readline support disabled (GA_USE_READLINE=OFF)")
        set(GA_HAS_READLINE FALSE PARENT_SCOPE)
    endif()

    # System dependencies that must be installed separately
    find_package(Qt6 COMPONENTS Core Gui Widgets QUIET)
    if(Qt6_FOUND)
        message(STATUS "✓ Found Qt6: ${Qt6_VERSION}")
        set(GA_HAS_QT6 TRUE PARENT_SCOPE)
    else()
        message(WARNING "Qt6 not found! Please install Qt6 separately")
        message(STATUS "  Windows: Download from https://www.qt.io/download")
        set(GA_HAS_QT6 FALSE PARENT_SCOPE)
    endif()


endfunction()