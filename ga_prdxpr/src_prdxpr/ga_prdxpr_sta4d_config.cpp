// Copyright 2024-2025, Daniel Hug. All rights reserved.

/*
 * GA PRODUCT RULE GENERATION SYSTEM GUIDE
 * =======================================
 *
 * This configuration file defines two types of product rule generation:
 *
 * 1. AUTOMATIC GENERATION (AlgebraConfig): Complement rules generated mathematically
 * 2. MANUAL CASES (ProductConfig): Expression cases configured manually
 *
 * AUTOMATIC RULE GENERATION:
 * --------------------------
 * The get_sta4d_algebra_config() function provides mathematical algebra parameters:
 * - basis_vectors: Basis vector names (e.g., {"e1", "e2", "e3", "e4"})
 * - metric_signature: Quadratic form values (e.g., {+1, +1, +1, 0})
 * - multivector_basis: Complete basis element names (e.g., {"1", "e1", "e2", "e3", "e4",
 * "e41", "e42", "e43", "e23", "e31", "e12", "e423", "e431", "e412", "e321", "e1234"})
 * - scalar_name: Name for scalar element (typically "1")
 * - basis_prefix: Prefix for basis elements (typically "e")
 *
 * From this configuration, the system automatically generates:
 * - Geometric product rules (geometric multiplication with metric)
 * - Wedge product rules (antisymmetric exterior product)
 * - Dot product rules (symmetric contraction with extended metric)
 * - Complement rules (computed from wedge product tables)
 *
 * COMPLEMENT RULE GENERATION:
 * The complement rules are now AUTOMATICALLY GENERATED using mathematical algorithms:
 * - Even algebras (EGA2D, PGA3DP): Generate left_complement and right_complement
 * - Odd algebras (EGA3D, PGA2DP): Generate single complement
 * - Algorithm: Search wedge product table for pseudoscalar relationships
 * - Result: Mathematically verified complement rules with zero transcription errors
 *
 * MANUAL CASE CONFIGURATION:
 * --------------------------
 * Each ProductConfig contains a .cases array with mathematical operation descriptions.
 * Cases define which coefficient combinations and type filters to use for code
 * generation.
 *
 * STANDARD FORMAT (5 parameters):
 * {"case_name", "left_coeff", "right_coeff", "left_filter", "right_filter"}
 *
 * SANDBOX FORMAT (7 parameters):
 * {"case_name", "left_coeff", "right_coeff", "left_filter", "right_filter", is_two_step,
 * "intermediate_name"}
 *
 * PARAMETER DESCRIPTIONS:
 * ----------------------
 * 1. case_name:         Mathematical operation description
 *    Examples: "mv * mv -> mv", "vec ^ s -> vec", "dot(bivec,bivec) -> s"
 *    Format: "left_type OPERATOR right_type -> result_type"
 *    Special: Use function notation for non-infix operators: "dot(A,B) -> result"
 *
 * 2. left_coeff_name:   Coefficient name for left operand (see COEFFICIENT SYSTEM below)
 * 3. right_coeff_name:  Coefficient name for right operand (see COEFFICIENT SYSTEM below)
 * 4. left_filter_name:  Type filter for left operand (see AVAILABLE FILTERS below)
 * 5. right_filter_name: Type filter for right operand (see AVAILABLE FILTERS below)
 *
 * SANDWICH PRODUCTS (7 parameters):
 * 6. is_two_step:       Always true for sandwich products
 * 7. intermediate_name: Name for intermediate result (e.g., "vec_tmp")
 *
 * COEFFICIENT SYSTEM LINKAGE:
 * ---------------------------
 * Coefficient names in config files map to mvec_coeff definitions in corresponding header
 * files. The mapping is established in create_[algebra]_algebra_data() via
 * AlgebraData.coefficients:
 *
 * Config Name -> Header Definition:
 * "A"         -> mv2d_coeff_A / mv2dp_coeff_A / mv3d_coeff_A / mv3dp_coeff_A
 * "B"         -> mv2d_coeff_B / mv2dp_coeff_B / mv3d_coeff_B / mv3dp_coeff_B
 * "A_even"    -> mv2d_coeff_A_even / mv2dp_coeff_A_even / etc.
 * "R_even"    -> mv2d_coeff_R_even / mv2dp_coeff_R_even / etc.
 * "svps"      -> mv2d_coeff_svps (EGA2D/3D naming)
 * "svBps"     -> mv2dp_coeff_svBps (PGA2DP/3DP naming)
 *
 * Each mvec_coeff contains component strings matching the algebra's basis:
 * - EGA2D: {"1", "e1", "e2", "e12"} -> {"A.c0", "A.c1", "A.c2", "A.c3"}
 * - PGA3DP: {"1", "e1", "e2", "e3", "e01", "e02", "e03", "e23", "e31", "e12", "e032",
 * "e013", "e021", "e123", "e0123"} -> {"A.c0"..."A.c15"}
 *
 * TO EXTEND THE SYSTEM:
 * 1. Define new mvec_coeff in the algebra header file (e.g., ga_prdxpr_pga3dp.hpp)
 * 2. Add mapping in create_[algebra]_algebra_data() coefficients map
 * 3. Use the coefficient name in .cases arrays
 *
 * Example extension for new coefficient "C":
 * Header: const mvec_coeff mv3dp_coeff_C = {"C.c0", "C.c1", ..., "C.c15"};
 * Config: pga3dp.coefficients = {..."C", mv3dp_coeff_C}...
 * Usage:  {"mv * mv -> mv", "C", "B", "mv", "mv"}
 *
 * AVAILABLE COEFFICIENTS (PGA3DP):
 * - General: A, B, M, M1, M2 (full multivectors)
 * - Even/Odd: A_even, B_even, A_odd, B_odd, M_even, M_odd (grade-filtered)
 * - Motors: R_even, R_odd, R_rev_even, R_rev_odd, R_rrev_even, R_rrev_odd
 * - Symmetric: svBtps, svBtps1, svBtps2 (scalar, vector, bivector, trivector,
 * pseudoscalar patterns)
 *
 * AVAILABLE FILTERS (PGA3DP):
 * - s (scalar), vec (vector), bivec (bivector), trivec (trivector), ps (pseudoscalar)
 * - mv_e (even), mv_u (odd), mv (multivector)
 *
 * NAMING CONVENTIONS:
 * - Use lowercase for all type names: s, vec, bivec, trivec, ps, mv, mv_e, mv_u
 * - Result types should match geometric algebra conventions
 * - Zero results: Use "0" as result type for operations that yield zero
 */

#include "ga_prdxpr_sta4d_config.hpp"

// Automatic rule generation configuration for sta4d
AlgebraConfig get_sta4d_algebra_config()
{
    // Extract basis prefix from vector basis and validate consistency
    std::string const prefix = extract_basis_prefix(mvsta4d_basis_kvec[1]);
    validate_basis_consistency(mvsta4d_basis, mvsta4d_basis_kvec, prefix, one_str());

    return {.basis_vectors = mvsta4d_basis_kvec[1],       // Use vector basis from header
            .metric_signature = mvsta4d_metric_signature, // Use metric from header
            .multivector_basis = mvsta4d_basis,           // Use mvsta4d_basis from header
            .scalar_name = one_str(),
            .basis_prefix = prefix}; // Use extracted and validated prefix
}

////////////////////////////////////////////////////////////////////////////////
// ALGEBRA CONFIGURATION - sta4d
//
// This file contains pure configuration for the sta4d geometric algebra.
// Product rule generation is handled in ga_prdxpr_sta4d.cpp.
//
// Configuration includes:
// - Algebra parameters (basis vectors, metric signature)
// - Coefficient definitions
// - Product case specifications
////////////////////////////////////////////////////////////////////////////////

namespace configurable {

AlgebraData create_sta4d_algebra_data()
{
    AlgebraData sta4d;
    sta4d.name = "sta4d";
    sta4d.basis = mvsta4d_basis; // Use existing basis
    sta4d.dimension = 4;         // g0: time, g1,g2,g3: space

    // Map coefficient names to existing coefficient objects
    sta4d.coefficients = {{"A", mvsta4d_coeff_A},
                          {"B", mvsta4d_coeff_B},
                          {"M", mvsta4d_coeff_M},
                          {"M1", mvsta4d_coeff_M1},
                          {"M2", mvsta4d_coeff_M2},
                          {"M_even", mvsta4d_coeff_M_even},
                          {"M_odd", mvsta4d_coeff_M_odd},
                          {"A_even", mvsta4d_coeff_A_even},
                          {"B_even", mvsta4d_coeff_B_even},
                          {"A_odd", mvsta4d_coeff_A_odd},
                          {"B_odd", mvsta4d_coeff_B_odd},
                          {"R_even", mvsta4d_coeff_R_even},
                          {"R_odd", mvsta4d_coeff_R_odd},
                          {"R_rev_even", mvsta4d_coeff_R_rev_even},
                          {"R_rev_odd", mvsta4d_coeff_R_rev_odd},
                          {"R_rrev_even", mvsta4d_coeff_R_rrev_even},
                          {"R_rrev_odd", mvsta4d_coeff_R_rrev_odd},
                          {"svBtps", mvsta4d_coeff_svBtps},
                          {"svBtps1", mvsta4d_coeff_svBtps1},
                          {"svBtps2", mvsta4d_coeff_svBtps2}};

    // Map filter names to existing filter enums
    sta4d.filters_4d = {{"s", filter_4d::s},         {"vec", filter_4d::vec},
                        {"bivec", filter_4d::bivec}, {"trivec", filter_4d::trivec},
                        {"ps", filter_4d::ps},       {"mv_e", filter_4d::mv_e},
                        {"mv_u", filter_4d::mv_u},   {"mv", filter_4d::mv}};

    return sta4d;
}

ProductConfig get_sta4d_gpr_config()
{
    return {.product_name = "gpr",
            .description = "sta4d geometric product",
            .display_name = "geometric product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_cmt_config()
{
    return {.product_name = "cmt",
            .description = "sta4d commutator product",
            .display_name = "commutator product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_wdg_config()
{
    return {.product_name = "wdg",
            .description = "sta4d wedge product",
            .display_name = "wedge product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_dot_config()
{
    return {.product_name = "dot",
            .description = "sta4d inner product",
            .display_name = "inner product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_left_contract_config()
{
    return {.product_name = "left_contract",
            .description = "sta4d left contraction",
            .display_name = "left contraction",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_right_contract_config()
{
    return {.product_name = "right_contract",
            .description = "sta4d right contraction",
            .display_name = "right contraction",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_left_expand_config()
{
    return {.product_name = "left_expand",
            .description = "sta4d left expansion",
            .display_name = "left expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_right_expand_config()
{
    return {.product_name = "right_expand",
            .description = "sta4d right expansion",
            .display_name = "right expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_rgpr_config()
{
    return {.product_name = "rgpr",
            .description = "sta4d regressive geometric product",
            .display_name = "regressive geometric product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_rcmt_config()
{
    return {.product_name = "rcmt",
            .description = "sta4d regressive commutator product",
            .display_name = "regressive commutator product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_rwdg_config()
{
    return {.product_name = "rwdg",
            .description = "sta4d regressive wedge product",
            .display_name = "regressive wedge product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_rdot_config()
{
    return {.product_name = "rdot",
            .description = "sta4d regressive inner product",
            .display_name = "regressive inner product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4d_sandwich_gpr_config()
{
    return {
        .product_name = "sandwich_gpr",
        .description = "sta4d sandwich product",
        .display_name = "sandwich product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter", is_two_step, "intermediate"}
        .cases =
            {// Single case that triggers sandwich product behavior - motor
             // operations
             {"dummy", "dummy", "dummy", "dummy", "dummy", true, "vec_tmp"}},
        .is_sandwich_product = true,
        .uses_brace_switch = true, // true needed for sandwich products
        .show_basis_table = true   // Reference shows basis table for sandwich product
    };
}

} // namespace configurable