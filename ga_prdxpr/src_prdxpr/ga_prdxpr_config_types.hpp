#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <map>
#include <string>
#include <vector>

// Mathematical computation functions are unchanged - these are just the common types
// that both the mathematical core and configuration system need to know about
#include "ga_prdxpr_common.hpp"

namespace configurable {

// Simple configuration for individual output cases
struct OutputCase {
    std::string case_name;         // "mv * mv -> mv"
    std::string left_coeff_name;   // "A" (maps to existing mv2d_coeff_A)
    std::string right_coeff_name;  // "B" (maps to existing mv2d_coeff_B)
    std::string left_filter_name;  // "mv" (maps to existing filter_2d::mv)
    std::string right_filter_name; // "mv" (maps to existing filter_2d::mv)

    // For sandwich products - specifies if this is a two-step operation
    bool is_two_step = false;
    std::string intermediate_name = ""; // For sandwich: intermediate result name

    // Constructor for simple cases (5 parameters)
    OutputCase(const std::string& case_n, const std::string& left_coeff,
               const std::string& right_coeff, const std::string& left_filter,
               const std::string& right_filter) :
        case_name(case_n), left_coeff_name(left_coeff), right_coeff_name(right_coeff),
        left_filter_name(left_filter), right_filter_name(right_filter)
    {
    }

    // Constructor for sandwich products (7 parameters)
    OutputCase(const std::string& case_n, const std::string& left_coeff,
               const std::string& right_coeff, const std::string& left_filter,
               const std::string& right_filter, bool two_step,
               const std::string& inter_name) :
        case_name(case_n), left_coeff_name(left_coeff), right_coeff_name(right_coeff),
        left_filter_name(left_filter), right_filter_name(right_filter),
        is_two_step(two_step), intermediate_name(inter_name)
    {
    }
};

// Configuration for a complete product type
struct ProductConfig {
    std::string product_name;      // "gpr", "wdg", "lcontract", etc.
    std::string description;       // "ega2d geometric product"
    std::string display_name;      // "geometric product" for output headers
    std::vector<OutputCase> cases; // What expressions to generate

    // Special handling flags
    bool is_sandwich_product = false;
    bool uses_brace_switch = false;
    bool show_basis_table = true; // Whether to show basis product table
};

// Minimal algebra data - just maps names to existing mathematical objects
struct AlgebraData {
    std::string name; // "ega2d", "ega3d", etc.

    // Maps from string names to existing coefficient objects
    std::map<std::string, mvec_coeff> coefficients;

    // Maps from string names to existing filter enums
    std::map<std::string, filter_2d> filters_2d;
    std::map<std::string, filter_3d> filters_3d;
    std::map<std::string, filter_4d> filters_4d;

    // Reference to existing basis (unchanged from reference implementation)
    mvec_coeff basis;

    // Dimension for filter selection
    int dimension;
};

} // namespace configurable