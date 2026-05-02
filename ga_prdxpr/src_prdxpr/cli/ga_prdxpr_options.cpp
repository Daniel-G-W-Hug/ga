// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "cli/ga_prdxpr_options.hpp"
#include <algorithm>
#include <fmt/core.h>
#include <sstream>

namespace configurable {

///////////////////////////////////////////////////////////////////////////////
// GeneratorOptions implementation
///////////////////////////////////////////////////////////////////////////////

bool GeneratorOptions::should_generate_algebra(std::string const& algebra_name) const
{
    return algebras.empty() || algebras.count(algebra_name) > 0;
}

bool GeneratorOptions::should_generate_product(std::string const& product_name) const
{
    return products.empty() || products.count(product_name) > 0;
}

bool GeneratorOptions::should_show_coeffs() const
{
    return output_types.empty() || output_types.count(OutputType::Coeffs) > 0;
}

bool GeneratorOptions::should_show_tables() const
{
    return output_types.empty() || output_types.count(OutputType::Tables) > 0;
}

bool GeneratorOptions::should_show_metrics() const
{
    return output_types.empty() || output_types.count(OutputType::Metrics) > 0;
}

bool GeneratorOptions::should_show_code() const
{
    // Opt-in: never on by default (would inflate the standard golden output).
    // Only emitted when the user explicitly asks for it.
    return output_types.count(OutputType::Code) > 0;
}

bool GeneratorOptions::is_filtering_active() const
{
    return !algebras.empty() || !products.empty() || !output_types.empty();
}

///////////////////////////////////////////////////////////////////////////////
// ArgumentParser implementation
///////////////////////////////////////////////////////////////////////////////

ArgumentParser::ArgumentParser(int argc, char const* argv[])
{
    // Skip program name (argv[0])
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--help" || arg == "-h") {
            show_help_ = true;
            return;
        }

        if (arg == "--list" || arg == "-l") {
            show_list_ = true;
            return;
        }

        parse_argument(arg);

        if (has_error_) {
            return;
        }
    }
}

void ArgumentParser::parse_argument(std::string const& arg)
{
    // Expect format: --key=value
    if (arg.substr(0, 2) != "--") {
        has_error_ = true;
        error_message_ = "Invalid argument format: " + arg + " (expected --key=value)";
        return;
    }

    auto equals_pos = arg.find('=');
    if (equals_pos == std::string::npos) {
        has_error_ = true;
        error_message_ = "Invalid argument format: " + arg + " (expected --key=value)";
        return;
    }

    std::string key = arg.substr(2, equals_pos - 2);
    std::string value = arg.substr(equals_pos + 1);

    if (key == "algebra" || key == "a") {
        parse_algebra_list(value);
    }
    else if (key == "products" || key == "product" || key == "p") {
        parse_product_list(value);
    }
    else if (key == "output" || key == "o") {
        parse_output_list(value);
    }
    else if (key == "show-sym-asym") {
        // Parse boolean value
        std::string lower_val = to_lower(trim(value));
        if (lower_val == "true" || lower_val == "1" || lower_val == "yes") {
            options_.show_symmetric_asymmetric = true;
        }
        else if (lower_val == "false" || lower_val == "0" || lower_val == "no") {
            options_.show_symmetric_asymmetric = false;
        }
        else {
            has_error_ = true;
            error_message_ =
                "Invalid value for --show-sym-asym: " + value + " (expected true/false)";
        }
    }
    else {
        has_error_ = true;
        error_message_ = "Unknown option: --" + key;
    }
}

void ArgumentParser::parse_algebra_list(std::string const& value)
{
    auto names = split(value, ',');
    for (auto& name : names) {
        name = trim(to_lower(name));

        // Handle special keywords
        if (name == "all") {
            // Insert all valid algebras
            options_.algebras.insert("ega2d");
            options_.algebras.insert("ega3d");
            options_.algebras.insert("pga2dp");
            options_.algebras.insert("pga3dp");
            options_.algebras.insert("sta4d");
            continue;
        }
        else if (name == "none") {
            // Clear all algebras
            options_.algebras.clear();
            continue;
        }

        if (!is_valid_algebra(name)) {
            has_error_ = true;
            error_message_ = "Invalid algebra name: " + name;
            return;
        }
        options_.algebras.insert(name);
    }
}

void ArgumentParser::parse_product_list(std::string const& value)
{
    auto names = split(value, ',');
    for (auto& name : names) {
        name = trim(name);

        // Handle special keywords (case-insensitive for keywords only)
        std::string lower_name = to_lower(name);
        if (lower_name == "all") {
            // Clear filter to show all products (empty set = all)
            options_.products.clear();
            continue;
        }
        else if (lower_name == "none") {
            // Set to special marker that matches nothing
            // Clear and insert a product name that will never match
            options_.products.clear();
            options_.products.insert("__NONE__");
            continue;
        }

        // Product names are case-sensitive and can contain special characters
        // We'll validate them when we try to match against actual product configs
        options_.products.insert(name);
    }
}

void ArgumentParser::parse_output_list(std::string const& value)
{
    auto types = split(value, ',');
    for (auto& type : types) {
        type = trim(to_lower(type));

        // Handle special keywords
        if (type == "all") {
            // Insert all output types
            options_.output_types.insert(OutputType::Coeffs);
            options_.output_types.insert(OutputType::Tables);
            options_.output_types.insert(OutputType::Metrics);
            options_.output_types.insert(OutputType::Code);
            continue;
        }
        else if (type == "none") {
            // Clear all output types
            options_.output_types.clear();
            continue;
        }

        if (!is_valid_output_type(type)) {
            has_error_ = true;
            error_message_ = "Invalid output type: " + type;
            return;
        }

        if (type == "coeffs") {
            options_.output_types.insert(OutputType::Coeffs);
        }
        else if (type == "tables") {
            options_.output_types.insert(OutputType::Tables);
        }
        else if (type == "metrics") {
            options_.output_types.insert(OutputType::Metrics);
        }
        else if (type == "code") {
            options_.output_types.insert(OutputType::Code);
        }
    }
}

bool ArgumentParser::is_valid_algebra(std::string const& name) const
{
    static std::set<std::string> const valid_algebras = {
        "ega2d", "ega3d", "pga2dp", "pga3dp", "sta4d", "all", "none"};
    return valid_algebras.count(name) > 0;
}

bool ArgumentParser::is_valid_product(std::string const& name) const
{
    // We'll validate against actual product names at runtime
    // For now, just check it's not empty
    return !name.empty();
}

bool ArgumentParser::is_valid_output_type(std::string const& name) const
{
    static std::set<std::string> const valid_types = {"coeffs", "tables", "metrics",
                                                      "code", "all", "none"};
    return valid_types.count(name) > 0;
}

std::vector<std::string> ArgumentParser::split(std::string const& str, char delimiter)
{
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string item;

    while (std::getline(ss, item, delimiter)) {
        if (!item.empty()) {
            result.push_back(item);
        }
    }

    return result;
}

std::string ArgumentParser::trim(std::string const& str)
{
    auto start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) {
        return "";
    }

    auto end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
}

std::string ArgumentParser::to_lower(std::string const& str)
{
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return result;
}

void ArgumentParser::print_help() const
{
    fmt::println("GA Product Expression Generator");
    fmt::println("================================");
    fmt::println("");
    fmt::println("Usage: ga_prdxpr [OPTIONS]");
    fmt::println("");
    fmt::println("Options:");
    fmt::println("  --help, -h              Show this help message");
    fmt::println("  --list, -l              List available algebras and product types");
    fmt::println("");
    fmt::println(
        "  --algebra=ALGEBRAS      Generate only specified algebras (comma-separated)");
    fmt::println(
        "                          Special: 'all' (all algebras), 'none' (no algebras)");
    fmt::println("                          Example: --algebra=ega2d,ega3d");
    fmt::println("                          Example: --algebra=all");
    fmt::println("");
    fmt::println(
        "  --products=PRODUCTS     Generate only specified products (comma-separated)");
    fmt::println(
        "                          Special: 'all' (all products), 'none' (no products)");
    fmt::println("                          Example: --products=gpr,wdg,dot");
    fmt::println("                          Example: --products=all");
    fmt::println("");
    fmt::println(
        "  --output=TYPES          Show only specified output types (comma-separated)");
    fmt::println("                          Types: coeffs, tables, metrics, code");
    fmt::println(
        "                          Special: 'all' (all types), 'none' (no types)");
    fmt::println("                          Note: 'code' is opt-in only");
    fmt::println("                          Example: --output=coeffs,tables");
    fmt::println("                          Example: --output=code");
    fmt::println("                          Example: --output=all");
    fmt::println("");
    fmt::println(
        "  --show-sym-asym=BOOL    Show symmetric/asymmetric parts of basis tables");
    fmt::println("                          Values: true/false (default: false)");
    fmt::println("                          Example: --show-sym-asym=true");
    fmt::println("");
    fmt::println("Common Workflows:");
    fmt::println(
        "  # View only product tables for one algebra (no coefficient expressions)");
    fmt::println("  ga_prdxpr --algebra=ega2d --output=tables");
    fmt::println("");
    fmt::println("  # View only coefficient expressions for one product (no tables)");
    fmt::println("  ga_prdxpr --algebra=ega3d --products=gpr --output=coeffs");
    fmt::println("");
    fmt::println("  # Compare one product across algebras");
    fmt::println("  ga_prdxpr --algebra=ega2d,ega3d --products=gpr --output=coeffs");
    fmt::println("");
    fmt::println("More Examples:");
    fmt::println(
        "  ga_prdxpr                              # Generate everything (default)");
    fmt::println("  ga_prdxpr --algebra=ega2d              # Only EGA2D algebra");
    fmt::println(
        "  ga_prdxpr --output=tables              # Only basis tables (all algebras)");
    fmt::println("  ga_prdxpr --output=tables --show-sym-asym=true");
    fmt::println("");
    fmt::println("With no arguments, generates complete output (searchable reference).");
}

void ArgumentParser::print_list() const
{
    fmt::println("Available Algebras:");
    fmt::println("  ega2d   - Euclidean Geometric Algebra 2D");
    fmt::println("  ega3d   - Euclidean Geometric Algebra 3D");
    fmt::println("  pga2dp  - Projective Geometric Algebra 2D");
    fmt::println("  pga3dp  - Projective Geometric Algebra 3D");
    fmt::println("  sta4d   - Space-Time Algebra 4D");
    fmt::println("");
    fmt::println("Common Product Types (all algebras):");
    fmt::println("  gpr                 - Geometric product");
    fmt::println("  wdg                 - Wedge (outer) product");
    fmt::println("  dot                 - Inner (dot) product");
    fmt::println("  cmt                 - Commutator product");
    fmt::println("  rgpr                - Regressive geometric product");
    fmt::println("  rwdg                - Regressive wedge product");
    fmt::println("  rcmt                - Regressive commutator product");
    fmt::println("  rdot                - Regressive inner product");
    fmt::println("");
    fmt::println("EGA Product Types (ega2d, ega3d):");
    fmt::println("  l_contract          - Left contraction");
    fmt::println("  r_contract          - Right contraction");
    fmt::println("  l_expand            - Left expansion");
    fmt::println("  r_expand            - Right expansion");
    fmt::println("  rgpr (alternative)  - Alternative regressive geometric product");
    fmt::println("  rtwdg1              - Regressive triple wedge (first variant)");
    fmt::println("  sandwich_gpr        - Sandwich product (geometric)");
    fmt::println("");
    fmt::println("STA4D Product Types (sta4d):");
    fmt::println("  l_contract          - Left contraction");
    fmt::println("  r_contract          - Right contraction");
    fmt::println("  l_expand            - Left expansion");
    fmt::println("  r_expand            - Right expansion");
    fmt::println("  sandwich_gpr        - Sandwich product (geometric)");
    fmt::println("");
    fmt::println("PGA Product Types (pga2dp, pga3dp):");
    fmt::println("  l_bulk_contract     - Left bulk contraction");
    fmt::println("  r_bulk_contract     - Right bulk contraction");
    fmt::println("  l_weight_contract   - Left weight contraction");
    fmt::println("  r_weight_contract   - Right weight contraction");
    fmt::println("  l_bulk_expand       - Left bulk expansion");
    fmt::println("  r_bulk_expand       - Right bulk expansion");
    fmt::println("  l_weight_expand     - Left weight expansion");
    fmt::println("  r_weight_expand     - Right weight expansion");
    fmt::println("  rgpr (alternative)  - Alternative regressive geometric product");
    fmt::println("  rtwdg1              - Regressive triple wedge (first variant)");
    fmt::println("  sandwich_rgpr       - Sandwich product (regressive)");
    fmt::println("");
    fmt::println("Output Types:");
    fmt::println("  coeffs              - Coefficient expressions");
    fmt::println("  tables              - Basis multiplication tables");
    fmt::println("  metrics             - Extended metric exports (C++ arrays)");
    fmt::println("  code                - C++ function implementations (opt-in only)");
}

} // namespace configurable
