# System dependency finding for GA project
# Handles finding system-installed packages

# Enhanced system dependency finding with three-tier resolution
function(find_system_dependencies)
    message(STATUS "Searching for dependencies (system -> ../../include/ -> FetchContent)...")
    
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
    
    # fmt - required for formatting (tier 1: system, tier 2: ../../include/, tier 3: FetchContent)
    find_fmt_tiered()
    
    # doctest - required for testing (tier 1: system, tier 2: ../../include/, tier 3: FetchContent)  
    find_doctest_tiered()
    
    # sol2 - header-only Lua binding (tier 1: system, tier 2: ../../include/, tier 3: FetchContent)
    find_sol2_tiered()
    
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
    if(Lua_FOUND OR LUA_FOUND)
        message(STATUS "✓ Found Lua: ${LUA_VERSION_STRING}")
        set(GA_HAS_LUA TRUE PARENT_SCOPE)
        # Make sure the variables are available in parent scope
        set(Lua_FOUND ${Lua_FOUND} PARENT_SCOPE)
        set(LUA_FOUND ${LUA_FOUND} PARENT_SCOPE)
        set(LUA_LIBRARIES ${LUA_LIBRARIES} PARENT_SCOPE)
        set(LUA_INCLUDE_DIR ${LUA_INCLUDE_DIR} PARENT_SCOPE)
        set(LUA_VERSION_STRING ${LUA_VERSION_STRING} PARENT_SCOPE)
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

# Three-tier fmt finding: system -> ../../include/ -> FetchContent
function(find_fmt_tiered)
    # Tier 1: System packages
    find_package(fmt CONFIG QUIET)
    if(fmt_FOUND)
        message(STATUS "✓ Found system fmt: ${fmt_VERSION}")
        set(GA_HAS_SYSTEM_FMT TRUE PARENT_SCOPE)
        set(GA_HAS_LOCAL_FMT FALSE PARENT_SCOPE)
        set(GA_HAS_FMT TRUE PARENT_SCOPE)  # Overall flag
        return()
    endif()
    
    # Tier 2: Local ../../include/ directory
    set(LOCAL_FMT_DIR "${CMAKE_SOURCE_DIR}/../../include/fmt")
    if(EXISTS "${LOCAL_FMT_DIR}/include/fmt/core.h")
        message(STATUS "✓ Found local fmt at: ${LOCAL_FMT_DIR}")
        if(NOT TARGET fmt::fmt)
            add_library(fmt::fmt INTERFACE IMPORTED)
            target_include_directories(fmt::fmt INTERFACE "${LOCAL_FMT_DIR}/include")
        endif()
        set(GA_HAS_SYSTEM_FMT FALSE PARENT_SCOPE)
        set(GA_HAS_LOCAL_FMT TRUE PARENT_SCOPE)
        set(GA_HAS_FMT TRUE PARENT_SCOPE)  # Overall flag
        return()
    endif()
    
    # Tier 3: Will use FetchContent (handled elsewhere)
    message(STATUS "✗ fmt not found in system or ../../include/. Will use FetchContent.")
    set(GA_HAS_SYSTEM_FMT FALSE PARENT_SCOPE)
    set(GA_HAS_LOCAL_FMT FALSE PARENT_SCOPE)
    set(GA_HAS_FMT FALSE PARENT_SCOPE)  # Overall flag
endfunction()

# Three-tier doctest finding: system -> ../../include/ -> FetchContent  
function(find_doctest_tiered)
    # Tier 1: System packages
    find_package(doctest QUIET)
    if(doctest_FOUND)
        message(STATUS "✓ Found system doctest")
        set(GA_HAS_SYSTEM_DOCTEST TRUE PARENT_SCOPE)
        set(GA_HAS_LOCAL_DOCTEST FALSE PARENT_SCOPE)
        set(GA_HAS_DOCTEST TRUE PARENT_SCOPE)  # Overall flag
        return()
    endif()
    
    # Tier 2: Local ../../include/ directory  
    set(LOCAL_DOCTEST_DIR "${CMAKE_SOURCE_DIR}/../../include/doctest")
    if(EXISTS "${LOCAL_DOCTEST_DIR}/doctest/doctest.h")
        message(STATUS "✓ Found local doctest at: ${LOCAL_DOCTEST_DIR}")
        if(NOT TARGET doctest::doctest)
            add_library(doctest::doctest INTERFACE IMPORTED)
            target_include_directories(doctest::doctest INTERFACE "${LOCAL_DOCTEST_DIR}")
        endif()
        set(GA_HAS_SYSTEM_DOCTEST FALSE PARENT_SCOPE)
        set(GA_HAS_LOCAL_DOCTEST TRUE PARENT_SCOPE)
        set(GA_HAS_DOCTEST TRUE PARENT_SCOPE)  # Overall flag
        return()
    endif()
    
    # Tier 3: Will use FetchContent (handled elsewhere)
    message(STATUS "✗ doctest not found in system or ../../include/. Will use FetchContent.")
    set(GA_HAS_SYSTEM_DOCTEST FALSE PARENT_SCOPE)
    set(GA_HAS_LOCAL_DOCTEST FALSE PARENT_SCOPE)
    set(GA_HAS_DOCTEST FALSE PARENT_SCOPE)  # Overall flag
endfunction()

# Three-tier sol2 finding: system -> ../../include/ -> FetchContent
function(find_sol2_tiered)
    # Tier 1: System packages
    find_path(SOL2_SYSTEM_INCLUDE_DIR 
              NAMES sol/sol.hpp
              PATHS ${PREFIX_HINTS}
              PATH_SUFFIXES include sol2/include
              DOC "sol2 header-only library include directory")
    
    if(SOL2_SYSTEM_INCLUDE_DIR)
        message(STATUS "✓ Found system sol2 at: ${SOL2_SYSTEM_INCLUDE_DIR}")
        if(NOT TARGET sol2::sol2)
            add_library(sol2::sol2 INTERFACE IMPORTED)
            target_include_directories(sol2::sol2 INTERFACE ${SOL2_SYSTEM_INCLUDE_DIR})
        endif()
        set(GA_HAS_SYSTEM_SOL2 TRUE PARENT_SCOPE)
        set(GA_HAS_LOCAL_SOL2 FALSE PARENT_SCOPE)
        set(GA_HAS_SOL2 TRUE PARENT_SCOPE)  # Overall flag
        return()
    endif()
    
    # Tier 2: Local ../../include/ directory
    set(LOCAL_SOL2_DIR "${CMAKE_SOURCE_DIR}/../../include/sol2")
    if(EXISTS "${LOCAL_SOL2_DIR}/include/sol/sol.hpp")
        message(STATUS "✓ Found local sol2 at: ${LOCAL_SOL2_DIR}")
        if(NOT TARGET sol2::sol2)
            add_library(sol2::sol2 INTERFACE IMPORTED)
            target_include_directories(sol2::sol2 INTERFACE "${LOCAL_SOL2_DIR}/include")
        endif()
        set(GA_HAS_SYSTEM_SOL2 FALSE PARENT_SCOPE)
        set(GA_HAS_LOCAL_SOL2 TRUE PARENT_SCOPE)
        set(GA_HAS_SOL2 TRUE PARENT_SCOPE)  # Overall flag
        return()
    endif()
    
    # Tier 3: Will use FetchContent (handled elsewhere)
    message(STATUS "✗ sol2 not found in system or ../../include/. Will use FetchContent.")
    set(GA_HAS_SYSTEM_SOL2 FALSE PARENT_SCOPE)
    set(GA_HAS_LOCAL_SOL2 FALSE PARENT_SCOPE)
    set(GA_HAS_SOL2 FALSE PARENT_SCOPE)  # Overall flag
endfunction()

# Three-tier readline finding: local compiled -> system -> skip (optional dependency)
function(find_readline_tiered)
    if(NOT GA_USE_READLINE)
        message(STATUS "readline support disabled")
        set(GA_HAS_SYSTEM_READLINE FALSE PARENT_SCOPE)
        set(GA_HAS_LOCAL_READLINE FALSE PARENT_SCOPE)
        set(GA_HAS_READLINE FALSE PARENT_SCOPE)
        return()
    endif()
    
    # Tier 1: Local compiled libraries (prioritized for development)
    set(LOCAL_READLINE_DIR "${CMAKE_SOURCE_DIR}/../../include/readline")
    if(EXISTS "${LOCAL_READLINE_DIR}/readline.h")
        # Enhanced search for compiled libraries in common build output directories
        find_library(READLINE_LOCAL_LIBRARY 
                     NAMES libreadline.dylib libreadline.so libreadline.a readline
                     PATHS ${LOCAL_READLINE_DIR}
                     PATH_SUFFIXES . shlib .libs lib
                     NO_DEFAULT_PATH)
        
        if(READLINE_LOCAL_LIBRARY)
            message(STATUS "✓ Found local compiled readline: ${READLINE_LOCAL_LIBRARY}")
            message(STATUS "  Using local version instead of system (for development)")
            
            # Find required terminal capabilities library for local readline
            find_library(TERMCAP_LIBRARY NAMES ncurses curses termcap)
            if(TERMCAP_LIBRARY)
                message(STATUS "  Found terminal library: ${TERMCAP_LIBRARY}")
            else()
                message(STATUS "  Warning: No terminal library found - linking may fail")
            endif()
            
            if(NOT TARGET readline::readline)
                add_library(readline::readline UNKNOWN IMPORTED)
                if(TERMCAP_LIBRARY)
                    set_target_properties(readline::readline PROPERTIES
                        IMPORTED_LOCATION ${READLINE_LOCAL_LIBRARY}
                        INTERFACE_INCLUDE_DIRECTORIES ${LOCAL_READLINE_DIR}
                        INTERFACE_LINK_LIBRARIES ${TERMCAP_LIBRARY})
                else()
                    set_target_properties(readline::readline PROPERTIES
                        IMPORTED_LOCATION ${READLINE_LOCAL_LIBRARY}
                        INTERFACE_INCLUDE_DIRECTORIES ${LOCAL_READLINE_DIR})
                endif()
            endif()
            set(GA_HAS_SYSTEM_READLINE FALSE PARENT_SCOPE)
            set(GA_HAS_LOCAL_READLINE TRUE PARENT_SCOPE)
            set(GA_HAS_READLINE TRUE PARENT_SCOPE)
            return()
        else()
            # Headers found but no compiled library - continue to check system
            message(STATUS "✓ Found local readline headers at: ${LOCAL_READLINE_DIR}")
            message(STATUS "  No compiled library found - checking system version")
        endif()
    endif()
    
    # Tier 2: System packages (fallback)
    find_library(READLINE_SYSTEM_LIBRARY readline)
    find_path(READLINE_SYSTEM_INCLUDE_DIR readline/readline.h)
    
    if(READLINE_SYSTEM_LIBRARY AND READLINE_SYSTEM_INCLUDE_DIR)
        message(STATUS "✓ Found system readline: ${READLINE_SYSTEM_LIBRARY}")
        if(NOT TARGET readline::readline)
            add_library(readline::readline UNKNOWN IMPORTED)
            set_target_properties(readline::readline PROPERTIES
                IMPORTED_LOCATION ${READLINE_SYSTEM_LIBRARY}
                INTERFACE_INCLUDE_DIRECTORIES ${READLINE_SYSTEM_INCLUDE_DIR})
        endif()
        set(GA_HAS_SYSTEM_READLINE TRUE PARENT_SCOPE)
        set(GA_HAS_LOCAL_READLINE FALSE PARENT_SCOPE)
        set(GA_HAS_READLINE TRUE PARENT_SCOPE)
        return()
    endif()
    
    # Tier 3: Headers-only fallback (if local headers exist but no system library)
    if(EXISTS "${LOCAL_READLINE_DIR}/readline.h")
        message(STATUS "✓ Using local readline headers with no library fallback")
        message(STATUS "  Warning: May have linking issues without compiled library")
        if(NOT TARGET readline::readline)
            add_library(readline::readline INTERFACE IMPORTED)
            target_include_directories(readline::readline INTERFACE ${LOCAL_READLINE_DIR})
        endif()
        set(GA_HAS_SYSTEM_READLINE FALSE PARENT_SCOPE)
        set(GA_HAS_LOCAL_READLINE TRUE PARENT_SCOPE)
        set(GA_HAS_READLINE TRUE PARENT_SCOPE)
        return()
    endif()
    
    # Tier 4: Not found (optional dependency - don't use FetchContent)
    message(STATUS "✗ readline not found (optional)")
    if(APPLE)
        message(STATUS "  Install with: brew install readline")
    endif()
    set(GA_HAS_SYSTEM_READLINE FALSE PARENT_SCOPE)
    set(GA_HAS_LOCAL_READLINE FALSE PARENT_SCOPE)
    set(GA_HAS_READLINE FALSE PARENT_SCOPE)
endfunction()