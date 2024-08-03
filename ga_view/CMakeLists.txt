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
            src/active_pt.cpp
            src/active_vec.cpp
            src/active_bivec.cpp
            src/active_projection.cpp
            src/active_reflection.cpp
            src/item_pt2d.cpp
            src/item_ln2d.cpp
            src/item_vt2d.cpp)

add_executable(${EXEC_NAME} ${SOURCES} )

# target link libraries have to be added AFTER add_executable or add_lnibrary!
target_include_directories(${EXEC_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
target_include_directories(${EXEC_NAME}
                           PRIVATE ${CMAKE_SOURCE_DIR}/../../include)

#
# Use the Widgets module from Qt 6.
#
# Find the QtWidgets library
find_package(
  Qt6
  COMPONENTS Core Gui Widgets 
  REQUIRED)

  # make fmt available
find_package(fmt CONFIG REQUIRED)

target_link_libraries(${EXEC_NAME} PRIVATE Qt6::Core
                                           Qt6::Gui
                                           Qt6::Widgets
                                           fmt::fmt-header-only)
