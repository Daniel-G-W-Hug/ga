// author: Daniel Hug, 2024

#include "ga_lua.hpp"

#include <iostream> // cout, cerr

int main(int argc, char* argv[])
{

    try {

        sol::state lua;
        lua.open_libraries();

        register_2d_types(lua);
        register_3d_types(lua);
        register_functions(lua);
        register_constants(lua);

        // run the lua script
        lua.safe_script_file("../../ga_lua/input/ga_lua.lua");
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