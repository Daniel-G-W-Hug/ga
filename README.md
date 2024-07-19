## Overview on modules

- ga: header-only library for geometric algebra in 2d and 3d
- ga_test: tests for ga header only library
- ga_lua: lua scripting based on ga
- ga_view: qt-based visualisation of 2d ga objects and operations

### Dependencies

- ga: fmt (must be included by every user - header only lib!)
- ga_test: fmt, doctest
- ga_lua: fmt, lua, sol2
- ga_view: fmt, qt6

#### installing dependencies on macos

    - brew install fmt
    - brew istalll doctest
    - brew install lua
    - brew install qt6

 sol2 to be installed from: <https://github.com/ThePhD/sol2>