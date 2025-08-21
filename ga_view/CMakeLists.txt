#[[
# directory structure:
# .
# ├── build
# ├── CMakeLists.txt
#     └── *.[cpp|cxx]
#
# ##############################################################################
# to configure execute from build directory:
#                cmake ..
#             or cmake .. -G "Ninja"
#             or cmake .. -G "Visual Studio 17 2022" -T host=x64 -A x64
#
# to build execute from build directory:
#
#                cmake --build .
#             or ninja
#
# ##############################################################################
#]]

set(EXEC_NAME ga_view)

#
# configure qt
#
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

set(SOURCES src/main.cpp
            src/w_mainwindow.cpp
            src/coordsys.cpp
            src/coordsys_model.cpp
            src/w_common.cpp
            src/w_coordsys.cpp
            src/w_statusbar.cpp
            src/active_common.cpp
            src/active_pt2d.cpp
            src/active_vt2d.cpp
            src/active_bivt2d.cpp
            src/active_projection.cpp
            src/active_reflection.cpp
            src/active_bivt2dp.cpp
            src/item_pt2d.cpp
            src/item_pt2dp.cpp
            src/item_ln2d.cpp
            src/item_cln2dp.cpp
            src/item_vt2d.cpp
            src/item_bivt2dp.cpp
            src/active_reflectionp.cpp
            src/active_kinematics2dp.cpp)

add_executable(${EXEC_NAME} ${SOURCES} )

# target link libraries have to be added AFTER add_executable or add_lnibrary!
target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/../../include)

# Dependencies are handled by the main CMakeLists.txt dependency system
# Just link against the targets that should be available

# Qt6 is required and should be found by the dependency system
if(NOT TARGET Qt6::Core)
    message(FATAL_ERROR "Qt6 not found - please install Qt6 on your system")
endif()

target_link_libraries(${EXEC_NAME} PRIVATE Qt6::Core Qt6::Gui Qt6::Widgets ga)

# Link fmt (prefer header-only if available, otherwise regular)
if(TARGET fmt::fmt-header-only)
    target_link_libraries(${EXEC_NAME} PRIVATE fmt::fmt-header-only)
elseif(TARGET fmt::fmt)
    target_link_libraries(${EXEC_NAME} PRIVATE fmt::fmt)
else()
    message(FATAL_ERROR "fmt library not found - dependency system may not be working correctly")
endif()
