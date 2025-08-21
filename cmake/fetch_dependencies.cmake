# FetchContent dependency management for GA project
# Handles downloading and building dependencies from source

include(FetchContent)

# Configure FetchContent behavior
set(FETCHCONTENT_QUIET OFF)  # Show download progress
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)  # Don't update on every build

# fmt library - header-only formatting
function(setup_fmt_fetchcontent)
    message(STATUS "Setting up fmt via FetchContent...")
    
    FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG        10.2.1  # Use a stable release tag
        GIT_SHALLOW    TRUE
    )
    
    # Configure fmt options
    set(FMT_DOC OFF CACHE BOOL "Generate the doc target.")
    set(FMT_INSTALL OFF CACHE BOOL "Generate the install target.")
    set(FMT_TEST OFF CACHE BOOL "Generate the test target.")
    
    FetchContent_MakeAvailable(fmt)
    
    message(STATUS "✓ fmt configured via FetchContent")
endfunction()

# doctest library - lightweight testing framework
function(setup_doctest_fetchcontent)
    message(STATUS "Setting up doctest via FetchContent...")
    
    FetchContent_Declare(
        doctest
        GIT_REPOSITORY https://github.com/doctest/doctest.git
        GIT_TAG        v2.4.11  # Use a stable release tag
        GIT_SHALLOW    TRUE
    )
    
    # Configure doctest options
    set(DOCTEST_WITH_TESTS OFF CACHE BOOL "Build tests/examples")
    set(DOCTEST_WITH_MAIN_IN_STATIC_LIB OFF CACHE BOOL "Build a static lib (cmake target) with a default main entry point")
    
    FetchContent_MakeAvailable(doctest)
    
    message(STATUS "✓ doctest configured via FetchContent")
endfunction()

# sol2 library - Lua binding (header-only)
function(setup_sol2_fetchcontent)
    message(STATUS "Setting up sol2 via FetchContent...")
    
    FetchContent_Declare(
        sol2
        GIT_REPOSITORY https://github.com/ThePhD/sol2.git
        GIT_TAG        v3.3.0  # Use a stable release tag
        GIT_SHALLOW    TRUE
    )
    
    # Configure sol2 options
    set(SOL2_BUILD_LUA OFF CACHE BOOL "Build Lua from source")
    set(SOL2_ENABLE_INSTALL OFF CACHE BOOL "Enable installation of Sol2")
    set(SOL2_SINGLE_HEADER ON CACHE BOOL "Generate single header")
    
    FetchContent_MakeAvailable(sol2)
    
    message(STATUS "✓ sol2 configured via FetchContent")
endfunction()

# Master function to setup all FetchContent dependencies
function(setup_all_fetchcontent_dependencies)
    message(STATUS "Setting up all dependencies via FetchContent...")
    
    setup_fmt_fetchcontent()
    setup_doctest_fetchcontent()
    setup_sol2_fetchcontent()
    
    message(STATUS "All FetchContent dependencies configured")
endfunction()

# Fallback setup for missing system dependencies
function(setup_missing_dependencies)
    # Only setup dependencies that weren't found on the system
    if(NOT GA_HAS_SYSTEM_FMT)
        setup_fmt_fetchcontent()
    endif()
    
    if(NOT GA_HAS_SYSTEM_DOCTEST)
        setup_doctest_fetchcontent()
    endif()
    
    if(NOT GA_HAS_SYSTEM_SOL2)
        setup_sol2_fetchcontent()
    endif()
endfunction()

# Platform-specific FetchContent optimizations
function(optimize_fetchcontent_for_platform)
    if(WIN32)
        # Windows-specific optimizations
        set(CMAKE_POLICY_DEFAULT_CMP0077 NEW CACHE STRING "")
        
        # Use faster generators if available
        if(CMAKE_GENERATOR MATCHES "Visual Studio")
            # Enable parallel builds in Visual Studio
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP" CACHE STRING "")
        endif()
        
    elseif(APPLE)
        # macOS-specific optimizations
        # Use system compiler optimizations
        set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE STRING "")
        
    else()
        # Linux-specific optimizations
        # Enable parallel builds
        include(ProcessorCount)
        ProcessorCount(N)
        if(NOT N EQUAL 0)
            set(CMAKE_BUILD_PARALLEL_LEVEL ${N} CACHE STRING "")
        endif()
    endif()
endfunction()