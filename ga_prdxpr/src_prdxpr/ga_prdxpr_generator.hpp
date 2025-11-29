#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"
#include "ga_prdxpr_config_types.hpp"
#include "ga_prdxpr_options.hpp"
#include "ga_prdxpr_sandwich_transformer.hpp"

namespace configurable {

class ConfigurableGenerator {
  public:

    // Main generation function - replaces all the individual generate_and_print_*
    // functions
    void generate_product_expressions(AlgebraData const& algebra,
                                      ProductConfig const& config,
                                      GeneratorOptions const& options);

  private:

    // Generate a single output case
    void generate_single_case(AlgebraData const& algebra, ProductConfig const& config,
                              OutputCase const& case_def, prd_table const& basis_tab);

    // Generate two-step sandwich product case
    void generate_sandwich_case(AlgebraData const& algebra, ProductConfig const& config,
                                OutputCase const& case_def, prd_table const& basis_tab);

    // Get the appropriate basis table for a product type
    // This uses the EXISTING mathematical functions from reference implementation
    prd_table get_basis_table_for_product(AlgebraData const& algebra,
                                          std::string const& product_name);

    // Helper to get filter based on algebra dimension
    filter_2d get_filter_2d(AlgebraData const& algebra, std::string const& filter_name);
    filter_3d get_filter_3d(AlgebraData const& algebra, std::string const& filter_name);
    filter_4d get_filter_4d(AlgebraData const& algebra, std::string const& filter_name);

    // Print formatted output matching reference implementation exactly
    void print_product_header(AlgebraData const& algebra, ProductConfig const& config);
    void print_basis_table(AlgebraData const& algebra, ProductConfig const& config,
                           prd_table const& basis_tab, bool show_sym_asym);
    void print_case_header(AlgebraData const& algebra, ProductConfig const& config,
                           std::string const& case_name);
    void print_case_result(mvec_coeff const& result, mvec_coeff const& basis);

    // Transform and print sandwich result in consistent multivector format
    void print_transformed_result(mvec_coeff const& result, mvec_coeff const& basis,
                                  AlgebraData const& algebra,
                                  ProductConfig const& config);

    // Helper function to apply coefficient alignment to transformation expressions
    void apply_coefficient_alignment(mvec_coeff& expressions,
                                     std::string const& algebra_name);
};

} // namespace configurable