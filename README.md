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

    - brew install fmt         // always required
    - brew istalll doctest     // only required for ga_test
    - brew install qt6         // only required for ga_view
    - brew install lua         // lua and sol2 only required for ga_lua
    - sol2 to be installed from: <https://github.com/ThePhD/sol2>

## Installing dependencies on Windows

    - fmt to be installed from: https://github.com/fmtlib/fmt
    - doctest to be installed from: https://github.com/doctest/doctest
    - winget install lua (you need a x64 binary otherwise linking will not work)
    - qt6 to be installed from official installer
    - sol2 to be installed from: <https://github.com/ThePhD/sol2>

## Usage of ga library

1) Make all available that is needed for GA, e.g. by

    ```
    using namespace hd::ga;
    using namespace hd::ga::ega;
    ```

2) Create a try catch block to catch any exceptions thrown by GA code

    ```
    try {

        // call GA stuff here

    } catch (const std::exception &e) {

        // whichever alternative suits your needs
        qDebug() << "Something went wrong: " << e.what();
        std::cout << "Something went wrong: " << e.what();
        fmt::print("Something went wrong: {}", e.what());

        // ... or similar

    }
    ```
