# Overview on modules

- ga: header-only library for geometric algebra in 2d and 3d
- ga_test: tests for ga header only library
- ga_lua: lua scripting based on ga
- ga_view: qt-based visualisation of 2d ga objects and operations

## Dependencies

- ga: fmt (must be included by every user - header only lib!)
- ga_test: fmt, doctest
- ga_lua: fmt, lua, sol2
- ga_view: fmt, qt6

## Installing dependencies on MacOS

    - brew install fmt
    - brew istalll doctest
    - brew install lua
    - brew install qt6
    - sol2 to be installed from: <https://github.com/ThePhD/sol2>

## Installing dependencies on Windows

    - fmt to be installed from: https://github.com/fmtlib/fmt
    - doctest to be installed from: https://github.com/doctest/doctest
    - winget install "Lua for Windows"
    - qt6 to be installed from official installer
    - sol2 to be installed from: <https://github.com/ThePhD/sol2>
