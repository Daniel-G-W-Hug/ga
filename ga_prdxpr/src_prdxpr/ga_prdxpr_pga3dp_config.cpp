// Copyright 2024-2025, Daniel Hug. All rights reserved.

/*
 * CASE DESCRIPTION FORMAT GUIDE
 * =============================
 *
 * Each ProductConfig contains a .cases array with mathematical operation descriptions.
 * Cases define which coefficient combinations and type filters to use for code
 * generation.
 *
 * STANDARD FORMAT (5 parameters):
 * {"case_name", "left_coeff", "right_coeff", "left_filter", "right_filter"}
 *
 * SANDWICH FORMAT (7 parameters):
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

#include "ga_prdxpr_pga3dp_config.hpp"

// Automatic rule generation configuration for PGA3DP
AlgebraConfig get_pga3dp_algebra_config()
{
    return {
        .basis_vectors = {"e1", "e2", "e3", "e4"},
        .metric_signature = {+1, +1, +1, 0}, // G(3,0,1) - e1²=+1, e2²=+1, e3²=+1, e4²=0
        .multivector_basis = {"1", "e1", "e2", "e3", "e4", "e41", "e42", "e43", "e23",
                              "e31", "e12", "e423", "e431", "e412", "e321",
                              "e1234"}, // Copy from mv3dp_basis
        .scalar_name = one_str,
        .basis_prefix = "e"};
}

////////////////////////////////////////////////////////////////////////////////
// ALGEBRA CONFIGURATION - PGA3DP
//
// This file contains pure configuration for the PGA3DP geometric algebra.
// Product rule generation is handled in ga_prdxpr_pga3dp.cpp.
//
// Configuration includes:
// - Algebra parameters (basis vectors, metric signature)
// - Coefficient definitions
// - Product case specifications
////////////////////////////////////////////////////////////////////////////////

namespace configurable {

AlgebraData create_pga3dp_algebra_data()
{
    AlgebraData pga3dp;
    pga3dp.name = "pga3dp";
    pga3dp.basis = mv3dp_basis; // Use existing basis
    pga3dp.dimension = 4;       // 3d projective = 4d homogeneous

    // Map coefficient names to existing coefficient objects
    pga3dp.coefficients = {{"A", mv3dp_coeff_A},
                           {"B", mv3dp_coeff_B},
                           {"M", mv3dp_coeff_M},
                           {"M1", mv3dp_coeff_M1},
                           {"M2", mv3dp_coeff_M2},
                           {"M_even", mv3dp_coeff_M_even},
                           {"M_odd", mv3dp_coeff_M_odd},
                           {"A_even", mv3dp_coeff_A_even},
                           {"B_even", mv3dp_coeff_B_even},
                           {"A_odd", mv3dp_coeff_A_odd},
                           {"B_odd", mv3dp_coeff_B_odd},
                           {"R_even", mv3dp_coeff_R_even},
                           {"R_odd", mv3dp_coeff_R_odd},
                           {"R_rev_even", mv3dp_coeff_R_rev_even},
                           {"R_rev_odd", mv3dp_coeff_R_rev_odd},
                           {"R_rrev_even", mv3dp_coeff_R_rrev_even},
                           {"R_rrev_odd", mv3dp_coeff_R_rrev_odd},
                           {"svBtps", mv3dp_coeff_svBtps},
                           {"svBtps1", mv3dp_coeff_svBtps1},
                           {"svBtps2", mv3dp_coeff_svBtps2}};

    // Map filter names to existing filter enums
    pga3dp.filters_4d = {{"s", filter_4d::s},         {"vec", filter_4d::vec},
                         {"bivec", filter_4d::bivec}, {"trivec", filter_4d::trivec},
                         {"ps", filter_4d::ps},       {"mv_e", filter_4d::mv_e},
                         {"mv_u", filter_4d::mv_u},   {"mv", filter_4d::mv}};

    return pga3dp;
}

ProductConfig get_pga3dp_gpr_config()
{
    return {
        .product_name = "gpr",
        .description = "pga3dp geometric product",
        .display_name = "geometric product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases = {{"mv * mv -> mv", "A", "B", "mv", "mv"},
                  {"mv * mv_e -> mv", "A", "B_even", "mv", "mv_e"},
                  {"mv_e * mv -> mv", "A_even", "B", "mv_e", "mv"},
                  {"mv * mv_u -> mv", "A", "B_odd", "mv", "mv_u"},
                  {"mv_u * mv -> mv", "A_odd", "B", "mv_u", "mv"},
                  {"mv * ps -> mv", "A", "svBtps", "mv", "ps"},
                  {"ps * mv -> mv", "svBtps", "B", "ps", "mv"},
                  {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                  {"mv_u * mv_u -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
                  {"mv_e * mv_u -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
                  {"mv_u * mv_e -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
                  {"mv_e * ps -> mv_e", "A_even", "svBtps", "mv_e", "ps"},
                  {"ps * mv_e -> mv_e", "svBtps", "B_even", "ps", "mv_e"},
                  {"mv_u * ps -> mv_u", "A_odd", "svBtps", "mv_u", "ps"},
                  {"ps * mv_u -> mv_u", "svBtps", "B_odd", "ps", "mv_u"},
                  {"mv_e * trivec -> mv_u", "M_even", "svBtps", "mv_e", "trivec"},
                  {"trivec * mv_e -> mv_u", "svBtps", "M_even", "trivec", "mv_e"},
                  {"mv_e * bivec -> mv_e", "M_even", "svBtps", "mv_e", "bivec"},
                  {"bivec * mv_e -> mv_e", "svBtps", "M_even", "bivec", "mv_e"},
                  {"mv_u * bivec -> mv_u", "M_odd", "svBtps", "mv_u", "bivec"},
                  {"bivec * mv_u -> mv_u", "svBtps", "M_odd", "bivec", "mv_u"},
                  {"mv_e * vec -> mv_u", "A_even", "svBtps", "mv_e", "vec"},
                  {"vec * mv_e -> mv_u", "svBtps", "B_even", "vec", "mv_e"},
                  {"ps * s -> ps", "svBtps", "svBtps", "ps", "s"},
                  {"s * ps -> ps", "svBtps", "svBtps", "s", "ps"},
                  {"trivec * trivec -> mv_e", "svBtps1", "svBtps2", "trivec", "trivec"},
                  {"trivec * bivec -> mv_u", "svBtps", "svBtps", "trivec", "bivec"},
                  {"bivec * trivec -> mv_u", "svBtps", "svBtps", "bivec", "trivec"},
                  {"trivec * vec -> mv_e", "svBtps", "svBtps", "trivec", "vec"},
                  {"vec * trivec -> mv_e", "svBtps", "svBtps", "vec", "trivec"},
                  {"trivec * s -> trivec", "svBtps", "svBtps", "trivec", "s"},
                  {"s * trivec -> trivec", "svBtps", "svBtps", "s", "trivec"},
                  {"bivec * bivec -> mv_e", "svBtps1", "svBtps2", "bivec", "bivec"},
                  {"bivec * vec -> mv_u", "svBtps", "svBtps", "bivec", "vec"},
                  {"vec * bivec -> mv_u", "svBtps", "svBtps", "vec", "bivec"},
                  {"bivec * s -> bivec", "svBtps", "svBtps", "bivec", "s"},
                  {"s * bivec -> bivec", "svBtps", "svBtps", "s", "bivec"},
                  {"vec * vec -> mv_e", "svBtps1", "svBtps2", "vec", "vec"},
                  {"vec * s -> vec", "svBtps", "svBtps", "vec", "s"},
                  {"s * vec -> vec", "svBtps", "svBtps", "s", "vec"},
                  {"s * s -> s", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga3dp_cmt_config()
{
    return {
        .product_name = "cmt",
        .description = "pga3dp commutator product",
        .display_name = "commutator product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases = {{"cmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
                  {"cmt(trivec,trivec) -> bivec", "svBtps1", "svBtps2", "trivec",
                   "trivec"},
                  {"cmt(trivec,bivec) -> trivec", "svBtps", "svBtps", "trivec", "bivec"},
                  {"cmt(bivec,trivec) -> trivec", "svBtps", "svBtps", "bivec", "trivec"},
                  {"cmt(bivec,bivec) -> bivec", "svBtps1", "svBtps2", "bivec", "bivec"},
                  {"cmt(bivec,vec) -> vec", "svBtps", "svBtps", "bivec", "vec"},
                  {"cmt(vec,bivec) -> vec", "svBtps", "svBtps", "vec", "bivec"},
                  {"cmt(vec,vec) -> bivec", "svBtps1", "svBtps2", "vec", "vec"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga3dp_wdg_config()
{
    return {.product_name = "wdg",
            .description = "pga3dp wedge product",
            .display_name = "wedge product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv ^ mv -> mv", "A", "B", "mv", "mv"},
                      {"ps ^ s -> ps", "svBtps", "svBtps", "ps", "s"},
                      {"s ^ ps -> ps", "svBtps", "svBtps", "s", "ps"},
                      {"trivec ^ vec -> ps", "svBtps", "svBtps", "trivec", "vec"},
                      {"vec ^ trivec -> ps", "svBtps", "svBtps", "vec", "trivec"},
                      {"trivec ^ s -> trivec", "svBtps", "svBtps", "trivec", "s"},
                      {"s ^ trivec -> trivec", "svBtps", "svBtps", "s", "trivec"},
                      {"bivec ^ bivec -> ps", "svBtps1", "svBtps2", "bivec", "bivec"},
                      {"bivec ^ vec -> trivec", "svBtps", "svBtps", "bivec", "vec"},
                      {"vec ^ bivec -> trivec", "svBtps", "svBtps", "vec", "bivec"},
                      {"bivec ^ s -> bivec", "svBtps", "svBtps", "bivec", "s"},
                      {"s ^ bivec -> bivec", "svBtps", "svBtps", "s", "bivec"},
                      {"vec ^ vec -> bivec", "svBtps1", "svBtps2", "vec", "vec"},
                      {"vec ^ s -> vec", "svBtps", "svBtps", "vec", "s"},
                      {"s ^ vec -> vec", "svBtps", "svBtps", "s", "vec"},
                      {"s ^ s -> s", "svBtps1", "svBtps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_dot_config()
{
    return {
        .product_name = "dot",
        .description = "pga3dp inner product",
        .display_name = "inner product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases = {{"dot(mv,mv) -> s", "A", "B", "mv", "mv"},
                  {"dot(ps,ps) -> s", "svBtps1", "svBtps2", "ps", "ps"},
                  {"dot(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec", "trivec"},
                  {"dot(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec"},
                  {"dot(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec"},
                  {"dot(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga3dp_rwdg_config()
{
    return {
        .product_name = "rwdg",
        .description = "pga3dp regressive wedge product",
        .display_name = "regressive wedge product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases = {{"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv"},
                  {"rwdg(ps,trivec) -> trivec", "svBtps", "svBtps", "ps", "trivec"},
                  {"rwdg(trivec,ps) -> trivec", "svBtps", "svBtps", "trivec", "ps"},
                  {"rwdg(trivec,trivec) -> bivec", "svBtps1", "svBtps2", "trivec",
                   "trivec"},
                  {"rwdg(trivec,bivec) -> vec", "svBtps", "svBtps", "trivec", "bivec"},
                  {"rwdg(bivec,trivec) -> vec", "svBtps", "svBtps", "bivec", "trivec"},
                  {"rwdg(trivec,vec) -> s", "svBtps", "svBtps", "trivec", "vec"},
                  {"rwdg(vec,trivec) -> s", "svBtps", "svBtps", "vec", "trivec"},
                  {"rwdg(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga3dp_rdot_config()
{
    return {
        .product_name = "rdot",
        .description = "pga3dp regressive inner product",
        .display_name = "regressive inner product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases = {{"rdot(mv,mv) -> ps", "A", "B", "mv", "mv"},
                  {"rdot(ps,ps) -> ps", "svBtps1", "svBtps2", "ps", "ps"},
                  {"rdot(trivec,trivec) -> ps", "svBtps1", "svBtps2", "trivec", "trivec"},
                  {"rdot(bivec,bivec) -> ps", "svBtps1", "svBtps2", "bivec", "bivec"},
                  {"rdot(vec,vec) -> ps", "svBtps1", "svBtps2", "vec", "vec"},
                  {"rdot(s,s) -> ps", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga3dp_rgpr_config()
{
    return {
        .product_name = "rgpr",
        .description = "pga3dp regressive geometric product",
        .display_name = "regressive geometric product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases = {{"rgpr(mv,mv) -> mv", "A", "B", "mv", "mv"},
                  {"rgpr(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                  {"rgpr(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
                  {"rgpr(mv_e,trivec) -> mv_u", "M_even", "svBtps", "mv_e", "trivec"},
                  {"rgpr(mv_u,trivec) -> mv_e", "M_odd", "svBtps", "mv_u", "trivec"},
                  {"rgpr(mv_e,bivec) -> mv_e", "M_even", "svBtps", "mv_e", "bivec"},
                  {"rgpr(mv_u,bivec) -> mv_u", "M_odd", "svBtps", "mv_u", "bivec"},
                  {"rgpr(mv_e,vec) -> mv_u", "M_even", "svBtps", "mv_e", "vec"},
                  {"rgpr(trivec,trivec) -> mv_e", "svBtps1", "svBtps2", "trivec",
                   "trivec"},
                  {"rgpr(trivec,bivec) -> mv_u", "svBtps", "svBtps", "trivec", "bivec"},
                  {"rgpr(trivec,vec) -> mv_e", "svBtps", "svBtps", "trivec", "vec"},
                  {"rgpr(bivec,bivec) -> mv_e", "svBtps1", "svBtps2", "bivec", "bivec"},
                  {"rgpr(vec,vec) -> mv_e", "svBtps1", "svBtps2", "vec", "vec"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga3dp_rcmt_config()
{
    return {
        .product_name = "rcmt",
        .description = "pga3dp regressive commutator product",
        .display_name = "regressive commutator product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases = {{"rcmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
                  {"rcmt(trivec,trivec) -> bivec", "svBtps1", "svBtps2", "trivec",
                   "trivec"},
                  {"rcmt(trivec,bivec) -> trivec", "svBtps", "svBtps", "trivec", "bivec"},
                  {"rcmt(bivec,trivec) -> trivec", "svBtps", "svBtps", "bivec", "trivec"},
                  {"rcmt(bivec,bivec) -> bivec", "svBtps1", "svBtps2", "bivec", "bivec"},
                  {"rcmt(bivec,vec) -> vec", "svBtps", "svBtps", "bivec", "vec"},
                  {"rcmt(vec,bivec) -> vec", "svBtps", "svBtps", "vec", "bivec"},
                  {"rcmt(vec,vec) -> bivec", "svBtps1", "svBtps2", "vec", "vec"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga3dp_right_bulk_contract_config()
{
    return {
        .product_name = "right_bulk_contract",
        .description = "pga3dp right bulk contraction",
        .display_name = "right bulk contraction",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {{"right_bulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"right_bulk_contract(ps,ps) -> 0", "svBtps1", "svBtps2", "ps", "ps"},
             {"right_bulk_contract(ps,trivec) -> vec", "svBtps", "svBtps", "ps",
              "trivec"},
             {"right_bulk_contract(trivec,ps) -> 0", "svBtps", "svBtps", "trivec", "ps"},
             {"right_bulk_contract(ps,bivec) -> bivec", "svBtps", "svBtps", "ps",
              "bivec"},
             {"right_bulk_contract(bivec,ps) -> 0", "svBtps", "svBtps", "bivec", "ps"},
             {"right_bulk_contract(ps,vec) -> trivec", "svBtps", "svBtps", "ps", "vec"},
             {"right_bulk_contract(vec,ps) -> 0", "svBtps", "svBtps", "vec", "ps"},
             {"right_bulk_contract(ps,s) -> ps", "svBtps", "svBtps", "ps", "s"},
             {"right_bulk_contract(s,ps) -> 0", "svBtps", "svBtps", "s", "ps"},
             {"right_bulk_contract(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec",
              "trivec"},
             {"right_bulk_contract(trivec,bivec) -> vec", "svBtps", "svBtps", "trivec",
              "bivec"},
             {"right_bulk_contract(bivec,trivec) -> 0", "svBtps", "svBtps", "bivec",
              "trivec"},
             {"right_bulk_contract(trivec,vec) -> bivec", "svBtps", "svBtps", "trivec",
              "vec"},
             {"right_bulk_contract(vec,trivec) -> 0", "svBtps", "svBtps", "vec",
              "trivec"},
             {"right_bulk_contract(trivec,s) -> trivec", "svBtps", "svBtps", "trivec",
              "s"},
             {"right_bulk_contract(s,trivec) -> 0", "svBtps", "svBtps", "s", "trivec"},
             {"right_bulk_contract(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec",
              "bivec"},
             {"right_bulk_contract(bivec,vec) -> vec", "svBtps", "svBtps", "bivec",
              "vec"},
             {"right_bulk_contract(vec,bivec) -> 0", "svBtps", "svBtps", "vec", "bivec"},
             {"right_bulk_contract(bivec,s) -> bivec", "svBtps", "svBtps", "bivec", "s"},
             {"right_bulk_contract(s,bivec) -> 0", "svBtps", "svBtps", "s", "bivec"},
             {"right_bulk_contract(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec"},
             {"right_bulk_contract(vec,s) -> vec", "svBtps", "svBtps", "vec", "s"},
             {"right_bulk_contract(s,vec) -> 0", "svBtps", "svBtps", "s", "vec"},
             {"right_bulk_contract(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}


ProductConfig get_pga3dp_right_weight_contract_config()
{
    return {.product_name = "right_weight_contract",
            .description = "pga3dp right weight contraction",
            .display_name = "right weight contraction",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {}, // For brevity, the weight contractions and expansions show basis
                         // tables but no specific
                         // cases yet
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_left_bulk_contract_config()
{
    return {
        .product_name = "left_bulk_contract",
        .description = "pga3dp left bulk contraction",
        .display_name = "left bulk contraction",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {{"left_bulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"left_bulk_contract(ps,ps) -> 0", "svBtps1", "svBtps2", "ps", "ps"},
             {"left_bulk_contract(ps,trivec) -> 0", "svBtps", "svBtps", "ps", "trivec"},
             {"left_bulk_contract(trivec,ps) -> vec", "svBtps", "svBtps", "trivec", "ps"},
             {"left_bulk_contract(ps,bivec) -> 0", "svBtps", "svBtps", "ps", "bivec"},
             {"left_bulk_contract(bivec,ps) -> bivec", "svBtps", "svBtps", "bivec", "ps"},
             {"left_bulk_contract(ps,vec) -> 0", "svBtps", "svBtps", "ps", "vec"},
             {"left_bulk_contract(vec,ps) -> trivec", "svBtps", "svBtps", "vec", "ps"},
             {"left_bulk_contract(ps,s) -> 0", "svBtps", "svBtps", "ps", "s"},
             {"left_bulk_contract(s,ps) -> ps", "svBtps", "svBtps", "s", "ps"},
             {"left_bulk_contract(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec",
              "trivec"},
             {"left_bulk_contract(trivec,bivec) -> 0", "svBtps", "svBtps", "trivec",
              "bivec"},
             {"left_bulk_contract(bivec,trivec) -> vec", "svBtps", "svBtps", "bivec",
              "trivec"},
             {"left_bulk_contract(trivec,vec) -> 0", "svBtps", "svBtps", "trivec", "vec"},
             {"left_bulk_contract(vec,trivec) -> bivec", "svBtps", "svBtps", "vec",
              "trivec"},
             {"left_bulk_contract(trivec,s) -> 0", "svBtps", "svBtps", "trivec", "s"},
             {"left_bulk_contract(s,trivec) -> trivec", "svBtps", "svBtps", "s",
              "trivec"},
             {"left_bulk_contract(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec",
              "bivec"},
             {"left_bulk_contract(bivec,vec) -> 0", "svBtps", "svBtps", "bivec", "vec"},
             {"left_bulk_contract(vec,bivec) -> vec", "svBtps", "svBtps", "vec", "bivec"},
             {"left_bulk_contract(bivec,s) -> 0", "svBtps", "svBtps", "bivec", "s"},
             {"left_bulk_contract(s,bivec) -> bivec", "svBtps", "svBtps", "s", "bivec"},
             {"left_bulk_contract(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec"},
             {"left_bulk_contract(vec,s) -> 0", "svBtps", "svBtps", "vec", "s"},
             {"left_bulk_contract(s,vec) -> vec", "svBtps", "svBtps", "s", "vec"},
             {"left_bulk_contract(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga3dp_left_weight_contract_config()
{
    return {.product_name = "left_weight_contract",
            .description = "pga3dp left weight contraction",
            .display_name = "left weight contraction",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_right_bulk_expand_config()
{
    return {.product_name = "right_bulk_expand",
            .description = "pga3dp right bulk expansion",
            .display_name = "right bulk expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_right_weight_expand_config()
{
    return {.product_name = "right_weight_expand",
            .description = "pga3dp right weight expansion",
            .display_name = "right weight expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_left_bulk_expand_config()
{
    return {.product_name = "left_bulk_expand",
            .description = "pga3dp left bulk expansion",
            .display_name = "left bulk expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_left_weight_expand_config()
{
    return {.product_name = "left_weight_expand",
            .description = "pga3dp left weight expansion",
            .display_name = "left weight expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_sandwich_rgpr_config()
{
    return {
        .product_name = "sandwich_rgpr",
        .description = "pga3dp regressive sandwich product",
        .display_name = "regressive sandwich product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter", is_two_step, "intermediate"}
        .cases =
            {// Single case that triggers regressive sandwich product behavior - motor
             // operations
             {"dummy", "dummy", "dummy", "dummy", "dummy", true, "vec_tmp"}},
        .is_sandwich_product = true,
        .uses_brace_switch = true, // true needed for sandwich products
        .show_basis_table = true   // Reference shows basis table for sandwich product
    };
}

} // namespace configurable