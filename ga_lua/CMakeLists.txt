# incrementally add test files needed

set(EXEC_NAME ga_lua)

add_executable(${EXEC_NAME} src/ga_lua.cpp)   #dep: lua, sol2, ga

target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/../../include)

# Dependencies are handled by the main CMakeLists.txt dependency system
# Just link against the targets that should be available

# Lua is required and should be found by the dependency system
if(NOT GA_HAS_LUA)
    message(FATAL_ERROR "Lua not found - please install Lua 5.1 or later on your system")
endif()

# Debug: Check what Lua variables we received
message(STATUS "ga_lua: DEBUG - LUA_INCLUDE_DIR=${LUA_INCLUDE_DIR}")
message(STATUS "ga_lua: DEBUG - LUA_LIBRARIES=${LUA_LIBRARIES}")

# CRITICAL: Lua include directory must come FIRST so sol2 finds the right lua.h
target_include_directories(${EXEC_NAME} PRIVATE ${LUA_INCLUDE_DIR})
target_link_libraries(${EXEC_NAME} PRIVATE ${LUA_LIBRARIES} ga)

# sol2 is header-only - add AFTER Lua includes to ensure correct lua.h is found
if(SOL2_INCLUDE_DIRS)
    target_include_directories(${EXEC_NAME} PRIVATE ${SOL2_INCLUDE_DIRS})
    message(STATUS "ga_lua: Using sol2 headers from ${SOL2_INCLUDE_DIRS}")
    message(STATUS "ga_lua: Include order: Lua (${LUA_INCLUDE_DIR}) then sol2 (${SOL2_INCLUDE_DIRS})")
else()
    # Fallback: try to find sol2::sol2 target for backwards compatibility
    if(TARGET sol2::sol2)
        target_link_libraries(${EXEC_NAME} PRIVATE sol2::sol2)
        message(STATUS "ga_lua: Using sol2::sol2 target")
    else()
        message(FATAL_ERROR "sol2 not found - dependency system may not be working correctly")
    endif()
endif()

# Helper function to link fmt to a target
function(link_fmt_to_target target_name)
    if(TARGET fmt::fmt-header-only)
        target_link_libraries(${target_name} PRIVATE fmt::fmt-header-only)
    elseif(TARGET fmt::fmt)
        target_link_libraries(${target_name} PRIVATE fmt::fmt)
        # On MSVC, fmt needs to be compiled - define FMT_HEADER_ONLY=0 to use compiled version
        if(MSVC)
            target_compile_definitions(${target_name} PRIVATE FMT_HEADER_ONLY=0)
        endif()
    else()
        message(FATAL_ERROR "fmt library not found - dependency system may not be working correctly")
    endif()
endfunction()

# Link fmt to executable
link_fmt_to_target(${EXEC_NAME})

# Note: DLL copying not needed - using static Lua library (lua54.lib) for clean deployment

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

# Fix for macOS C++ runtime linking - sol2 requires proper hash symbol resolution
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(STATUS "ga_lua: Configuring robust C++ runtime linking for sol2 compatibility")
    
    # Comprehensive Homebrew LLVM detection
    set(NEEDS_EXPLICIT_LIBCXX FALSE)
    
    # Check 1: Compiler path explicitly contains Homebrew
    if(CMAKE_CXX_COMPILER MATCHES "/opt/homebrew")
        message(STATUS "ga_lua: Detected Homebrew LLVM compiler: ${CMAKE_CXX_COMPILER}")
        set(NEEDS_EXPLICIT_LIBCXX TRUE)
    endif()
    
    # Check 2: Homebrew paths in CMAKE_PREFIX_PATH
    if(CMAKE_PREFIX_PATH MATCHES "/opt/homebrew")
        message(STATUS "ga_lua: Detected Homebrew in CMAKE_PREFIX_PATH")
        set(NEEDS_EXPLICIT_LIBCXX TRUE)
    endif()
    
    # Check 3: Check if LUA_LIBRARIES comes from Homebrew
    if(LUA_LIBRARIES MATCHES "/opt/homebrew")
        message(STATUS "ga_lua: Detected Homebrew Lua libraries")
        set(NEEDS_EXPLICIT_LIBCXX TRUE)
    endif()
    
    # Check 4: Check various CMake variables for Homebrew paths
    get_cmake_property(_variableNames VARIABLES)
    foreach(_variableName ${_variableNames})
        if("${${_variableName}}" MATCHES "/opt/homebrew")
            message(STATUS "ga_lua: Detected Homebrew in ${_variableName}")
            set(NEEDS_EXPLICIT_LIBCXX TRUE)
            break()
        endif()
    endforeach()
    
    if(NEEDS_EXPLICIT_LIBCXX)
        message(STATUS "ga_lua: Applying Homebrew LLVM-compatible linking strategy")
        # Essential for Homebrew LLVM: explicit stdlib specification
        target_compile_options(${EXEC_NAME} PRIVATE -stdlib=libc++)
        target_link_options(${EXEC_NAME} PRIVATE -stdlib=libc++)
        
        # Critical: Add Homebrew LLVM C++ library path to ensure we link against the right libc++
        target_link_options(${EXEC_NAME} PRIVATE -L/opt/homebrew/opt/llvm/lib/c++)
        target_link_options(${EXEC_NAME} PRIVATE -Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++)
        
        # Link the complete C++ runtime to resolve __hash_memory symbols using explicit paths
        target_link_libraries(${EXEC_NAME} PRIVATE /opt/homebrew/opt/llvm/lib/c++/libc++.dylib)
        target_link_libraries(${EXEC_NAME} PRIVATE /opt/homebrew/opt/llvm/lib/c++/libc++abi.dylib)
    else()
        # Safety fallback: always ensure C++ runtime is linked for sol2's hash operations
        message(STATUS "ga_lua: Applying standard C++ runtime linking for sol2 compatibility")
        target_link_libraries(${EXEC_NAME} PRIVATE c++)
    endif()
endif()

