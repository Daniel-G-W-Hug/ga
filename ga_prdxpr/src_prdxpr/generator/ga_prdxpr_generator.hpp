#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "cli/ga_prdxpr_options.hpp"
#include "core/ga_prdxpr_common.hpp"
#include "core/ga_prdxpr_config_types.hpp"
#include "sandwich/ga_prdxpr_sandwich_transformer.hpp"

namespace configurable {

// Print a summary of all OutputCase validation warnings emitted during the
// run, then clear the buffer. Writes to stderr to match the inline warnings.
// Stays silent (no header, no output) if zero warnings were collected.
void print_validation_summary();

class ConfigurableGenerator {
  public:

    // Main generation function - replaces all the individual generate_and_print_*
    // functions
    void generate_product_expressions(AlgebraData const& algebra,
                                      ProductConfig const& config,
                                      GeneratorOptions const& options);

    // Emit the unary complements / duals (l_cmpl, r_cmpl, l_dual, r_dual, and the
    // PGA bulk/weight variants) as ready-to-paste C++ for --output=code. These are
    // single-operand rule-table ops (not binary ProductConfigs), so they have their
    // own emission path: scatter each graded input type's components through the
    // complement/dual rule, then emit via codegen::emit_unary_function. Aggregate
    // inputs (mv/mv_e/mv_u) are skipped -- they stay hand-written delegations.
    void emit_unary_products_code(AlgebraData const& algebra,
                                  GeneratorOptions const& options);

  private:

    // Generate a single output case
    void generate_single_case(AlgebraData const& algebra, ProductConfig const& config,
                              OutputCase const& case_def, prd_table const& basis_tab);

    // Generate two-step sandwich product case
    void generate_sandwich_case(AlgebraData const& algebra, ProductConfig const& config,
                                OutputCase const& case_def, prd_table const& basis_tab);

    // Emit a single non-sandwich case as ready-to-paste C++ code.
    // Sandwich cases are intentionally skipped (handled separately).
    void emit_single_case_code(AlgebraData const& algebra, ProductConfig const& config,
                               OutputCase const& case_def, prd_table const& basis_tab);

    // Get the appropriate basis table for a product type
    // This uses the EXISTING mathematical functions from reference implementation
    prd_table get_basis_table_for_product(AlgebraData const& algebra,
                                          std::string const& product_name);

    // Non-orthogonal metrics: the geometric product of basis blades is
    // multi-term, so its basis table is an mt_table and case generation goes
    // through get_mv_from_mt_tab. These mirror the single-term path.
    static bool uses_mt_basis_table(AlgebraData const& algebra,
                                    std::string const& product_name);
    mt_table get_mt_basis_table_for_product(AlgebraData const& algebra,
                                            std::string const& product_name);
    void generate_single_case_mt(AlgebraData const& algebra, ProductConfig const& config,
                                 OutputCase const& case_def, mt_table const& mt_tab);
    void emit_single_case_code_mt(AlgebraData const& algebra, ProductConfig const& config,
                                  OutputCase const& case_def, mt_table const& mt_tab);
    void generate_sandwich_case_mt(AlgebraData const& algebra,
                                   ProductConfig const& config, mt_table const& mt_tab);

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