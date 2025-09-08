# incrementally add test files needed

set(EXEC_NAME1 ga_ega_test)
set(EXEC_NAME2 ga_pga_test)
set(EXEC_NAME3 ga_appl_test)

set(EGA_FILES
    src/ga_ega_test.cpp
    src/ga_ega2d_test.cpp 
    src/ga_ega3d_test.cpp
    src/ga_ega4d_test.cpp)

set(PGA_FILES
    src/ga_pga_test.cpp
    src/ga_pga2dp_test.cpp
    src/ga_pga3dp_test.cpp)

    set(APPL_FILES
    src/ga_appl_test.cpp
    src/ga_appl2dp_test.cpp
    src/ga_2d_2dp_compare_test.cpp)

add_executable(${EXEC_NAME1} ${EGA_FILES})   #dep: doctest, fmt, ga
add_executable(${EXEC_NAME2} ${PGA_FILES})   #dep: doctest, fmt, ga
add_executable(${EXEC_NAME3} ${APPL_FILES})  #dep: doctest, fmt, ga

target_include_directories(${EXEC_NAME1} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME2} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME3} PRIVATE ${CMAKE_SOURCE_DIR})

# Dependencies are handled by the main CMakeLists.txt dependency system
# Just link against the targets that should be available

target_link_libraries(${EXEC_NAME1} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME2} PRIVATE doctest::doctest ga)
target_link_libraries(${EXEC_NAME3} PRIVATE doctest::doctest ga)

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

# Link fmt to all executables
link_fmt_to_target(${EXEC_NAME1})
link_fmt_to_target(${EXEC_NAME2})
link_fmt_to_target(${EXEC_NAME3})