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
    
    # First, attempt to find system dependencies
    find_system_dependencies()
    
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
    
    # Setup optional dependencies
    find_readline_optional()
endfunction()

