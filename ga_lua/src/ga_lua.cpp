// author: Daniel Hug, 2024 & 2025

#include "ga_lua.hpp"

#include <exception> // std::exception
#include <iostream>  // std::cerr

// int main(int argc, char* argv[])
int main()
{

    try {

        sol::state lua;
        lua.open_libraries();

        register_2d_types(lua);
        register_3d_types(lua);
        register_functions(lua);
        register_constants(lua);

        // run the lua script
#if _MSC_VER >= 1900
        // on the msvc compiler there will be debug and release subfolders
        // (thus we need one level more to get back)
        lua.safe_script_file("../../../ga_lua/input/ga_lua.lua");
#else
        lua.safe_script_file("../../ga_lua/input/ga_lua.lua");
#endif
    }
    catch (sol::error const& e) {

        // some lua error occured
        std::cerr << e.what() << "\n";
        return -3;
    }
    catch (std::exception const& e) {
        // some non-lua error occured
        std::cerr << e.what() << "\n";
        return -2;
    }
    catch (...) {
        std::cerr << "ERROR: other, unexpected stuff happened.\n";
        return -1;
    }

    return 0;
}