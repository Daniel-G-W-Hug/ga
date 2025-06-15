#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"
#include "ga_prdxpr_config_types.hpp"

namespace configurable {

class ConfigurableGenerator {
  public:

    // Main generation function - replaces all the individual generate_and_print_*
    // functions
    void generate_product_expressions(const AlgebraData& algebra,
                                      const ProductConfig& config);

  private:

    // Generate a single output case
    void generate_single_case(const AlgebraData& algebra, const ProductConfig& config,
                              const OutputCase& case_def, const prd_table& basis_tab);

    // Generate two-step sandwich product case
    void generate_sandwich_case(const AlgebraData& algebra, const ProductConfig& config,
                                const OutputCase& case_def, const prd_table& basis_tab);

    // Get the appropriate basis table for a product type
    // This uses the EXISTING mathematical functions from reference implementation
    prd_table get_basis_table_for_product(const AlgebraData& algebra,
                                          const std::string& product_name);

    // Helper to get filter based on algebra dimension
    filter_2d get_filter_2d(const AlgebraData& algebra, const std::string& filter_name);
    filter_3d get_filter_3d(const AlgebraData& algebra, const std::string& filter_name);
    filter_4d get_filter_4d(const AlgebraData& algebra, const std::string& filter_name);

    // Print formatted output matching reference implementation exactly
    void print_product_header(const AlgebraData& algebra, const ProductConfig& config);
    void print_basis_table(const AlgebraData& algebra, const ProductConfig& config,
                           const prd_table& basis_tab);
    void print_case_header(const AlgebraData& algebra, const ProductConfig& config,
                           const std::string& case_name);
    void print_case_result(const mvec_coeff& result, const mvec_coeff& basis);
};

} // namespace configurable