#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <set>
#include <string>
#include <vector>

namespace configurable {

// Output type flags
enum class OutputType { Products, Tables, Metrics };

// Generator runtime options
struct GeneratorOptions {
    // Algebra filters (empty = all algebras)
    std::set<std::string> algebras;

    // Product filters (empty = all products)
    std::set<std::string> products;

    // Output type filters (empty = all types)
    std::set<OutputType> output_types;

    // Table output options
    bool show_symmetric_asymmetric = false; // Show sym/asym parts of basis tables

    // Helper functions
    bool should_generate_algebra(std::string const& algebra_name) const;
    bool should_generate_product(std::string const& product_name) const;
    bool should_show_products() const;
    bool should_show_tables() const;
    bool should_show_metrics() const;
    bool is_filtering_active() const;
};

// Command-line argument parser
class ArgumentParser {
  public:

    ArgumentParser(int argc, char const* argv[]);

    GeneratorOptions const& get_options() const { return options_; }

    bool should_show_help() const { return show_help_; }
    bool should_show_list() const { return show_list_; }
    bool has_error() const { return has_error_; }
    std::string const& get_error_message() const { return error_message_; }

    // Display functions
    void print_help() const;
    void print_list() const;

  private:

    GeneratorOptions options_;
    bool show_help_ = false;
    bool show_list_ = false;
    bool has_error_ = false;
    std::string error_message_;

    // Parsing helpers
    void parse_argument(std::string const& arg);
    void parse_algebra_list(std::string const& value);
    void parse_product_list(std::string const& value);
    void parse_output_list(std::string const& value);

    // Validation
    bool is_valid_algebra(std::string const& name) const;
    bool is_valid_product(std::string const& name) const;
    bool is_valid_output_type(std::string const& name) const;

    // String utilities
    static std::vector<std::string> split(std::string const& str, char delimiter);
    static std::string trim(std::string const& str);
    static std::string to_lower(std::string const& str);
};

} // namespace configurable
