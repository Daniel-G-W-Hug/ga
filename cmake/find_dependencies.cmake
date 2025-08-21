# System dependency finding for GA project
# Handles finding system-installed packages

# Enhanced system dependency finding with better error messages
function(find_system_dependencies)
    message(STATUS "Searching for system dependencies...")
    
    # Platform-specific package manager hints
    if(APPLE)
        set(PKG_MANAGER_HINT "Install with: brew install <package>")
        set(PREFIX_HINTS /opt/homebrew /usr/local)
    elseif(WIN32)
        set(PKG_MANAGER_HINT "Install with vcpkg or download manually")
        set(PREFIX_HINTS "C:/vcpkg/installed/${VCPKG_TARGET_TRIPLET}")
    else()
        set(PKG_MANAGER_HINT "Install with your system package manager")
        set(PREFIX_HINTS /usr /usr/local)
    endif()
    
    # fmt - required for formatting
    find_package(fmt CONFIG QUIET)
    if(fmt_FOUND)
        message(STATUS "✓ Found system fmt: ${fmt_VERSION}")
        set(GA_HAS_SYSTEM_FMT TRUE PARENT_SCOPE)
    else()
        message(STATUS "✗ fmt not found. ${PKG_MANAGER_HINT}")
        set(GA_HAS_SYSTEM_FMT FALSE PARENT_SCOPE)
    endif()
    
    # doctest - required for testing
    find_package(doctest QUIET)
    if(doctest_FOUND)
        message(STATUS "✓ Found system doctest")
        set(GA_HAS_SYSTEM_DOCTEST TRUE PARENT_SCOPE)
    else()
        message(STATUS "✗ doctest not found. ${PKG_MANAGER_HINT}")
        set(GA_HAS_SYSTEM_DOCTEST FALSE PARENT_SCOPE)
    endif()
    
    # sol2 - header-only Lua binding (special handling)
    find_sol2_system()
    
    # Qt6 - must be system installed
    find_package(Qt6 COMPONENTS Core Gui Widgets QUIET)
    if(Qt6_FOUND)
        message(STATUS "✓ Found Qt6: ${Qt6_VERSION}")
        set(GA_HAS_QT6 TRUE PARENT_SCOPE)
    else()
        message(WARNING "Qt6 not found! Please install:")
        if(APPLE)
            message(STATUS "  macOS: brew install qt6")
        elseif(WIN32)
            message(STATUS "  Windows: Download from https://www.qt.io/download")
        else()
            message(STATUS "  Linux: sudo apt-get install qt6-base-dev qt6-tools-dev")
        endif()
        set(GA_HAS_QT6 FALSE PARENT_SCOPE)
    endif()
    
    # Lua - must be system installed
    find_package(Lua 5.1 QUIET)
    if(Lua_FOUND)
        message(STATUS "✓ Found Lua: ${LUA_VERSION_STRING}")
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
endfunction()

# Special handling for sol2 (header-only library)
function(find_sol2_system)
    # Look for sol2 in common locations
    find_path(SOL2_INCLUDE_DIR 
              NAMES sol/sol.hpp
              PATHS ${PREFIX_HINTS}
              PATH_SUFFIXES include sol2/include
              DOC "sol2 header-only library include directory")
    
    if(SOL2_INCLUDE_DIR)
        # Create imported target for system sol2
        if(NOT TARGET sol2::sol2)
            add_library(sol2::sol2 INTERFACE IMPORTED)
            target_include_directories(sol2::sol2 INTERFACE ${SOL2_INCLUDE_DIR})
        endif()
        message(STATUS "✓ Found system sol2 at: ${SOL2_INCLUDE_DIR}")
        set(GA_HAS_SYSTEM_SOL2 TRUE PARENT_SCOPE)
    else()
        message(STATUS "✗ sol2 not found. Will use FetchContent.")
        set(GA_HAS_SYSTEM_SOL2 FALSE PARENT_SCOPE)
    endif()
endfunction()

# Optional: Enhanced readline detection
function(find_readline_optional)
    if(GA_USE_READLINE)
        find_library(READLINE_LIBRARY readline)
        find_path(READLINE_INCLUDE_DIR readline/readline.h)
        
        if(READLINE_LIBRARY AND READLINE_INCLUDE_DIR)
            message(STATUS "✓ Found readline: ${READLINE_LIBRARY}")
            if(NOT TARGET readline::readline)
                add_library(readline::readline UNKNOWN IMPORTED)
                set_target_properties(readline::readline PROPERTIES
                    IMPORTED_LOCATION ${READLINE_LIBRARY}
                    INTERFACE_INCLUDE_DIRECTORIES ${READLINE_INCLUDE_DIR})
            endif()
            set(GA_HAS_READLINE TRUE PARENT_SCOPE)
        else()
            message(STATUS "✗ readline not found (optional)")
            if(APPLE)
                message(STATUS "  Install with: brew install readline")
            endif()
            set(GA_HAS_READLINE FALSE PARENT_SCOPE)
        endif()
    endif()
endfunction()