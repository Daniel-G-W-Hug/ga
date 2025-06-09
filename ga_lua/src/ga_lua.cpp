// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_lua.hpp"

#include "fmt/format.h" // fmt::print, fmt::println
#include <exception>    // std::exception
#include <filesystem>   // std::filesystem::exists
#include <iostream>     // std::cerr, std::cin
#include <string>       // std::string

#ifdef _HD_GA_HAVE_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

void print_usage(const char* program_name)
{
    fmt::println("Usage: {} [options] [script_file]\n", program_name);
    fmt::println("Options:");
    fmt::println("  -h, --help     Show this help message");
    fmt::println(
        "  -i, --interactive  Start interactive shell (default if no arguments)");
    fmt::println("  -d, --default  Run the default demo script\n");
    fmt::println("Arguments:");
    fmt::println("  script_file    Lua script file to execute\n");
    fmt::println("Examples:");
    fmt::println("  {}                     # Interactive shell", program_name);
    fmt::println("  {} -i                  # Interactive shell", program_name);
    fmt::println("  {} -d                  # Run default demo script", program_name);
    fmt::println("  {} my_script.lua       # Run custom script file", program_name);
}

void run_interactive_shell(sol::state& lua)
{
    fmt::println("GA Lua Interactive Shell");
    fmt::println("Type 'quit' or 'exit' to leave, 'help' for GA functions");
    fmt::println("All GA types and functions are available (e.g., vec2d.new(1,2))");
#ifdef _HD_GA_HAVE_READLINE
    fmt::println("Use arrow keys for history, Tab for completion\n");
#else
    fmt::println();
#endif

    std::string input;
    int line_number = 1;

    while (true) {
#ifdef _HD_GA_HAVE_READLINE
        char* line = readline(fmt::format("ga_lua[{}]> ", line_number).c_str());
        if (!line) {
            // EOF (Ctrl+D on Unix, Ctrl+Z on Windows)
            fmt::println("\nExiting...");
            break;
        }
        input = line;
        
        // Add non-empty lines to history
        if (!input.empty()) {
            add_history(line);
        }
        
        free(line);
#else
        fmt::print("ga_lua[{}]> ", line_number);
        
        if (!std::getline(std::cin, input)) {
            // EOF (Ctrl+D on Unix, Ctrl+Z on Windows)
            fmt::println("\nExiting...");
            break;
        }
#endif

        // Trim whitespace and check for special commands
        if (input.empty()) {
            continue;
        }

        if (input == "quit" || input == "exit") {
            break;
        }

        if (input == "help") {
            fmt::println("Available GA functions and types:");
            fmt::println("  2D: scalar2d, vec2d, pscalar2d, mvec2d_e, mvec2d");
            fmt::println(
                "  3D: scalar3d, vec3d, bivec3d, pscalar3d, mvec3d_e, mvec3d_u, mvec3d");
            fmt::println("  Functions: dot, wdg, nrm, normalize, inv, dual, etc.");
            fmt::println(
                "  Constants: e1_2d, e2_2d, I_2d, e1_3d, e2_3d, e3_3d, I_3d, etc.");
            fmt::println("  Example: v = vec2d.new(1, 2); print(v); print(nrm(v))\n");
            continue;
        }

        try {
            // Execute the input as Lua code
            auto result = lua.safe_script(input, sol::script_pass_on_error);
            if (!result.valid()) {
                sol::error err = result;
                std::cerr << "Error: " << err.what() << "\n";
            }
        }
        catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << "\n";
        }

        line_number++;
    }
}

std::string get_default_script_path()
{
#if _MSC_VER >= 1900
    // on the msvc compiler there will be debug and release subfolders
    return "../../../ga_lua/input/ga_lua.lua";
#else
    return "../../ga_lua/input/ga_lua.lua";
#endif
}

int main(int argc, char* argv[])
{
    bool interactive_mode = false;
    bool use_default_script = false;
    std::string script_file;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        }
        else if (arg == "-i" || arg == "--interactive") {
            interactive_mode = true;
        }
        else if (arg == "-d" || arg == "--default") {
            use_default_script = true;
        }
        else if (arg.size() > 0 && arg[0] == '-') {
            std::cerr << "Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
        else {
            // Assume it's a script file
            if (!script_file.empty()) {
                std::cerr << "Multiple script files specified\n";
                print_usage(argv[0]);
                return 1;
            }
            script_file = arg;
        }
    }

    // Determine execution mode
    if (use_default_script) {
        script_file = get_default_script_path();
        interactive_mode = false;
    }
    else if (script_file.empty() && !interactive_mode) {
        // Default to interactive mode if no script file specified and no explicit mode
        interactive_mode = true;
    }

    try {
        sol::state lua;
        lua.open_libraries();

        register_2d_types(lua);
        register_3d_types(lua);
        register_2dp_types(lua);
        register_3dp_types(lua);
        register_functions(lua);
        register_constants(lua);

        if (interactive_mode) {
            run_interactive_shell(lua);
        }
        else {
            // Check if script file exists
            if (!std::filesystem::exists(script_file)) {
                if (use_default_script) {
                    std::cerr << "Error: Default script '" << script_file
                              << "' not found\n";
                    std::cerr << "Make sure you're running from the correct directory\n";
                }
                else {
                    std::cerr << "Error: Script file '" << script_file << "' not found\n";
                }
                return 1;
            }

            if (use_default_script) {
                fmt::println("Running default demo script: {}\n", script_file);
            }
            else {
                fmt::println("Running script: {}\n", script_file);
            }
            lua.safe_script_file(script_file);
        }
    }
    catch (sol::error const& e) {
        std::cerr << "Lua error: " << e.what() << "\n";
        return -3;
    }
    catch (std::exception const& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return -2;
    }
    catch (...) {
        std::cerr << "ERROR: Unexpected error occurred\n";
        return -1;
    }

    return 0;
}