// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

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
 * The get_sta4ds_algebra_config() function provides mathematical algebra parameters:
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
 * - Even algebras (EGA2D, PGA3DP): Generate l_cmpl and r_cmpl
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
 * Usage:  {"gpr(mv,mv) -> mv", "C", "B", "mv", "mv"}
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

#include "algebras/ga_prdxpr_sta4ds_config.hpp"

// Automatic rule generation configuration for sta4ds
AlgebraConfig get_sta4ds_algebra_config()
{
    // Extract basis prefix from vector basis and validate consistency
    std::string const prefix = extract_basis_prefix(mvsta4ds_basis_kvec[1]);
    validate_basis_consistency(mvsta4ds_basis, mvsta4ds_basis_kvec, prefix, one_str());

    return {.basis_vectors = mvsta4ds_basis_kvec[1],       // Use vector basis from header
            .metric_signature = mvsta4ds_metric_signature, // Use metric from header
            .multivector_basis = mvsta4ds_basis, // Use mvsta4ds_basis from header
            .scalar_name = one_str(),
            .basis_prefix = prefix}; // Use extracted and validated prefix
}

////////////////////////////////////////////////////////////////////////////////
// ALGEBRA CONFIGURATION - sta4ds
//
// This file contains pure configuration for the sta4ds geometric algebra.
// Product rule generation is handled in ga_prdxpr_sta4ds.cpp.
//
// Configuration includes:
// - Algebra parameters (basis vectors, metric signature)
// - Coefficient definitions
// - Product case specifications
////////////////////////////////////////////////////////////////////////////////

namespace configurable {

AlgebraData create_sta4ds_algebra_data()
{
    AlgebraData sta4ds;
    sta4ds.name = "sta4ds";
    sta4ds.basis = mvsta4ds_basis; // Use existing basis
    sta4ds.dimension = 4;          // g0: time, g1,g2,g3: space

    // Map coefficient names to existing coefficient objects
    sta4ds.coefficients = {{"A", mvsta4ds_coeff_A},
                           {"B", mvsta4ds_coeff_B},
                           {"M", mvsta4ds_coeff_M},
                           {"M1", mvsta4ds_coeff_M1},
                           {"M2", mvsta4ds_coeff_M2},
                           {"M_even", mvsta4ds_coeff_M_even},
                           {"M_odd", mvsta4ds_coeff_M_odd},
                           {"A_even", mvsta4ds_coeff_A_even},
                           {"B_even", mvsta4ds_coeff_B_even},
                           {"A_odd", mvsta4ds_coeff_A_odd},
                           {"B_odd", mvsta4ds_coeff_B_odd},
                           {"R_even", mvsta4ds_coeff_R_even},
                           {"R_odd", mvsta4ds_coeff_R_odd},
                           {"R_rev_even", mvsta4ds_coeff_R_rev_even},
                           {"R_rev_odd", mvsta4ds_coeff_R_rev_odd},
                           {"R_rrev_even", mvsta4ds_coeff_R_rrev_even},
                           {"R_rrev_odd", mvsta4ds_coeff_R_rrev_odd},
                           {"svBtps", mvsta4ds_coeff_svBtps},
                           {"svBtps1", mvsta4ds_coeff_svBtps1},
                           {"svBtps2", mvsta4ds_coeff_svBtps2}};

    // Map filter names to existing filter enums
    sta4ds.filters_4d = {{"s", filter_4d::s},         {"vec", filter_4d::vec},
                         {"bivec", filter_4d::bivec}, {"trivec", filter_4d::trivec},
                         {"ps", filter_4d::ps},       {"mv_e", filter_4d::mv_e},
                         {"mv_u", filter_4d::mv_u},   {"mv", filter_4d::mv}};

    return sta4ds;
}

ProductConfig get_sta4ds_gpr_config()
{
    return {.product_name = "gpr",
            .description = "sta4ds geometric product",
            .display_name = "geometric product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases =
                {// mv
                 {"gpr(mv,mv) -> mv", "A", "B", "mv", "mv"},
                 {"gpr(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
                 {"gpr(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
                 {"gpr(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
                 {"gpr(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
                 {"gpr(mv,ps) -> mv", "M", "svBtps", "mv", "ps"},
                 {"gpr(ps,mv) -> mv", "svBtps", "M", "ps", "mv"},
                 {"gpr(mv,trivec) -> mv", "M", "svBtps", "mv", "trivec"},
                 {"gpr(trivec,mv) -> mv", "svBtps", "M", "trivec", "mv"},
                 {"gpr(mv,bivec) -> mv", "M", "svBtps", "mv", "bivec"},
                 {"gpr(bivec,mv) -> mv", "svBtps", "M", "bivec", "mv"},
                 {"gpr(mv,vec) -> mv", "M", "svBtps", "mv", "vec"},
                 {"gpr(vec,mv) -> mv", "svBtps", "M", "vec", "mv"},
                 {"gpr(mv,s) -> mv", "M", "svBtps", "mv", "s"},
                 {"gpr(s,mv) -> mv", "svBtps", "M", "s", "mv"},
                 // mv_e
                 {"gpr(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                 {"gpr(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
                 {"gpr(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
                 {"gpr(mv_e,ps) -> mv_e", "A_even", "svBtps", "mv_e", "ps"},
                 {"gpr(ps,mv_e) -> mv_e", "svBtps", "B_even", "ps", "mv_e"},
                 {"gpr(mv_e,trivec) -> mv_u", "M_even", "svBtps", "mv_e", "trivec"},
                 {"gpr(trivec,mv_e) -> mv_u", "svBtps", "M_even", "trivec", "mv_e"},
                 {"gpr(mv_e,bivec) -> mv_e", "M_even", "svBtps", "mv_e", "bivec"},
                 {"gpr(bivec,mv_e) -> mv_e", "svBtps", "M_even", "bivec", "mv_e"},
                 {"gpr(mv_e,vec) -> mv_u", "A_even", "svBtps", "mv_e", "vec"},
                 {"gpr(vec,mv_e) -> mv_u", "svBtps", "B_even", "vec", "mv_e"},
                 {"gpr(mv_e,s) -> mv_e", "M_even", "svBtps", "mv_e", "s"},
                 {"gpr(s,mv_e) -> mv_e", "svBtps", "M_even", "s", "mv_e"},
                 // mv_u
                 {"gpr(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
                 {"gpr(mv_u,ps) -> mv_u", "A_odd", "svBtps", "mv_u", "ps"},
                 {"gpr(ps,mv_u) -> mv_u", "svBtps", "B_odd", "ps", "mv_u"},
                 {"gpr(mv_u,trivec) -> mv_e", "M_odd", "svBtps", "mv_u", "trivec"},
                 {"gpr(trivec,mv_u) -> mv_e", "svBtps", "M_odd", "trivec", "mv_u"},
                 {"gpr(mv_u,bivec) -> mv_u", "M_odd", "svBtps", "mv_u", "bivec"},
                 {"gpr(bivec,mv_u) -> mv_u", "svBtps", "M_odd", "bivec", "mv_u"},
                 {"gpr(mv_u,vec) -> mv_e", "M_odd", "svBtps", "mv_u", "vec"},
                 {"gpr(vec,mv_u) -> mv_e", "svBtps", "M_odd", "vec", "mv_u"},
                 {"gpr(mv_u,s) -> mv_u", "M_odd", "svBtps", "mv_u", "s"},
                 {"gpr(s,mv_u) -> mv_u", "svBtps", "M_odd", "s", "mv_u"},
                 // ps
                 {"gpr(ps,ps) -> s", "svBtps1", "svBtps2", "ps", "ps"},
                 {"gpr(ps,trivec) -> vec", "svBtps", "svBtps", "ps", "trivec"},
                 {"gpr(trivec,ps) -> vec", "svBtps", "svBtps", "trivec", "ps"},
                 {"gpr(ps,bivec) -> bivec", "svBtps", "svBtps", "ps", "bivec"},
                 {"gpr(bivec,ps) -> bivec", "svBtps", "svBtps", "bivec", "ps"},
                 {"gpr(ps,vec) -> trivec", "svBtps", "svBtps", "ps", "vec"},
                 {"gpr(vec,ps) -> trivec", "svBtps", "svBtps", "vec", "ps"},
                 {"gpr(ps,s) -> ps", "svBtps", "svBtps", "ps", "s"},
                 {"gpr(s,ps) -> ps", "svBtps", "svBtps", "s", "ps"},
                 // trivec
                 {"gpr(trivec,trivec) -> mv_e", "svBtps1", "svBtps2", "trivec", "trivec"},
                 {"gpr(trivec,bivec) -> mv_u", "svBtps", "svBtps", "trivec", "bivec"},
                 {"gpr(bivec,trivec) -> mv_u", "svBtps", "svBtps", "bivec", "trivec"},
                 {"gpr(trivec,vec) -> mv_e", "svBtps", "svBtps", "trivec", "vec"},
                 {"gpr(vec,trivec) -> mv_e", "svBtps", "svBtps", "vec", "trivec"},
                 {"gpr(trivec,s) -> trivec", "svBtps", "svBtps", "trivec", "s"},
                 {"gpr(s,trivec) -> trivec", "svBtps", "svBtps", "s", "trivec"},
                 // bivec
                 {"gpr(bivec,bivec) -> mv_e", "svBtps1", "svBtps2", "bivec", "bivec"},
                 {"gpr(bivec,vec) -> mv_u", "svBtps", "svBtps", "bivec", "vec"},
                 {"gpr(vec,bivec) -> mv_u", "svBtps", "svBtps", "vec", "bivec"},
                 {"gpr(bivec,s) -> bivec", "svBtps", "svBtps", "bivec", "s"},
                 {"gpr(s,bivec) -> bivec", "svBtps", "svBtps", "s", "bivec"},
                 // vec
                 {"gpr(vec,vec) -> mv_e", "svBtps1", "svBtps2", "vec", "vec"},
                 {"gpr(vec,s) -> vec", "svBtps", "svBtps", "vec", "s"},
                 {"gpr(s,vec) -> vec", "svBtps", "svBtps", "s", "vec"},
                 // s
                 {"gpr(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4ds_cmt_config()
{
    return {
        .product_name = "cmt",
        .description = "sta4ds commutator product",
        .display_name = "commutator product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {// mv
             {"cmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
             // ps
             {"cmt(ps,ps) -> 0", "svBtps1", "svBtps2", "ps", "ps"},
             {"cmt(ps,trivec) -> vec", "svBtps", "svBtps", "ps", "trivec"},
             {"cmt(trivec,ps) -> vec", "svBtps", "svBtps", "trivec", "ps"},
             {"cmt(ps,bivec) -> 0", "svBtps", "svBtps", "ps", "bivec"},
             {"cmt(bivec,ps) -> 0", "svBtps", "svBtps", "bivec", "ps"},
             {"cmt(ps,vec) -> trivec", "svBtps", "svBtps", "ps", "vec"},
             {"cmt(vec,ps) -> trivec", "svBtps", "svBtps", "vec", "ps"},
             {"cmt(ps,s) -> 0 ps", "svBtps", "svBtps", "ps", "s"},
             {"cmt(s,ps) -> 0 ps", "svBtps", "svBtps", "s", "ps"},
             // trivec
             {"cmt(trivec,trivec) -> bivec", "svBtps1", "svBtps2", "trivec", "trivec"},
             {"cmt(trivec,bivec) -> trivec", "svBtps", "svBtps", "trivec", "bivec"},
             {"cmt(bivec,trivec) -> trivec", "svBtps", "svBtps", "bivec", "trivec"},
             {"cmt(trivec,vec) -> ps", "svBtps", "svBtps", "trivec", "vec"},
             {"cmt(vec,trivec) -> ps", "svBtps", "svBtps", "vec", "trivec"},
             {"cmt(trivec,s) -> 0", "svBtps", "svBtps", "trivec", "s"},
             {"cmt(s,trivec) -> 0", "svBtps", "svBtps", "s", "trivec"},
             // bivec
             {"cmt(bivec,bivec) -> bivec", "svBtps1", "svBtps2", "bivec", "bivec"},
             {"cmt(bivec,vec) -> vec", "svBtps", "svBtps", "bivec", "vec"},
             {"cmt(vec,bivec) -> vec", "svBtps", "svBtps", "vec", "bivec"},
             {"cmt(bivec,s) -> 0", "svBtps", "svBtps", "bivec", "s"},
             {"cmt(s,bivec) -> 0", "svBtps", "svBtps", "s", "bivec"},
             // vec
             {"cmt(vec,vec) -> bivec", "svBtps1", "svBtps2", "vec", "vec"},
             {"cmt(vec,s) -> 0", "svBtps", "svBtps", "vec", "s"},
             {"cmt(s,vec) -> 0", "svBtps", "svBtps", "s", "vec"},
             // s
             {"cmt(s,s) -> 0", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_sta4ds_wdg_config()
{
    return {.product_name = "wdg",
            .description = "sta4ds wedge product",
            .display_name = "wedge product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases =
                {// mv
                 {"wdg(mv,mv) -> mv", "A", "B", "mv", "mv"},
                 {"wdg(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
                 {"wdg(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
                 {"wdg(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
                 {"wdg(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
                 {"wdg(mv,ps) -> ps", "M", "svBtps", "mv", "ps"},
                 {"wdg(ps,mv) -> ps", "svBtps", "M", "ps", "mv"},
                 {"wdg(mv,trivec) -> mv", "M", "svBtps", "mv", "trivec"},
                 {"wdg(trivec,mv) -> mv", "svBtps", "M", "trivec", "mv"},
                 {"wdg(mv,bivec) -> mv", "M", "svBtps", "mv", "bivec"},
                 {"wdg(bivec,mv) -> mv", "svBtps", "M", "bivec", "mv"},
                 {"wdg(mv,vec) -> mv", "M", "svBtps", "mv", "vec"},
                 {"wdg(vec,mv) -> mv", "svBtps", "M", "vec", "mv"},
                 {"wdg(mv,s) -> mv", "M", "svBtps", "mv", "s"},
                 {"wdg(s,mv) -> mv", "svBtps", "M", "s", "mv"},
                 // mv_e
                 {"wdg(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                 {"wdg(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
                 {"wdg(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
                 {"wdg(mv_e,ps) -> ps", "M_even", "svBtps", "mv_e", "ps"},
                 {"wdg(ps,mv_e) -> ps", "svBtps", "M_even", "ps", "mv_e"},
                 {"wdg(mv_e,trivec) -> trivec", "M_even", "svBtps", "mv_e", "trivec"},
                 {"wdg(trivec,mv_e) -> trivec", "svBtps", "M_even", "trivec", "mv_e"},
                 {"wdg(mv_e,bivec) -> mv_e", "M_even", "svBtps", "mv_e", "bivec"},
                 {"wdg(bivec,mv_e) -> mv_e", "svBtps", "M_even", "bivec", "mv_e"},
                 {"wdg(mv_e,vec) -> mv_u", "M_even", "svBtps", "mv_e", "vec"},
                 {"wdg(vec,mv_e) -> mv_u", "svBtps", "M_even", "vec", "mv_e"},
                 {"wdg(mv_e,s) -> mv_e", "M_even", "svBtps", "mv_e", "s"},
                 {"wdg(s,mv_e) -> mv_e", "svBtps", "M_even", "s", "mv_e"},
                 // mv_u
                 {"wdg(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
                 {"wdg(mv_u,ps) -> 0", "M_odd", "svBtps", "mv_u", "ps"},
                 {"wdg(ps,mv_u) -> 0", "svBtps", "M_odd", "ps", "mv_u"},
                 {"wdg(mv_u,trivec) -> ps", "M_odd", "svBtps", "mv_u", "trivec"},
                 {"wdg(trivec,mv_u) -> ps", "svBtps", "M_odd", "trivec", "mv_u"},
                 {"wdg(mv_u,bivec) -> trivec", "M_odd", "svBtps", "mv_u", "bivec"},
                 {"wdg(bivec,mv_u) -> trivec", "svBtps", "M_odd", "bivec", "mv_u"},
                 {"wdg(mv_u,vec) -> mv_e", "M_odd", "svBtps", "mv_u", "vec"},
                 {"wdg(vec,mv_u) -> mv_e", "svBtps", "M_odd", "vec", "mv_u"},
                 {"wdg(mv_u,s) -> mv_u", "M_odd", "svBtps", "mv_u", "s"},
                 {"wdg(s,mv_u) -> mv_u", "svBtps", "M_odd", "s", "mv_u"},
                 // ps
                 {"wdg(ps,ps) -> 0", "svBtps1", "svBtps2", "ps", "ps"},
                 {"wdg(ps,trivec) -> 0", "svBtps", "svBtps", "ps", "trivec"},
                 {"wdg(trivec,ps) -> 0", "svBtps", "svBtps", "trivec", "ps"},
                 {"wdg(ps,bivec) -> 0", "svBtps", "svBtps", "ps", "bivec"},
                 {"wdg(bivec,ps) -> 0", "svBtps", "svBtps", "bivec", "ps"},
                 {"wdg(ps,vec) -> 0", "svBtps", "svBtps", "ps", "vec"},
                 {"wdg(vec,ps) -> 0", "svBtps", "svBtps", "vec", "ps"},
                 {"wdg(ps,s) -> ps", "svBtps", "svBtps", "ps", "s"},
                 {"wdg(s,ps) -> ps", "svBtps", "svBtps", "s", "ps"},
                 // trivec
                 {"wdg(trivec,trivec) -> 0", "svBtps1", "svBtps2", "trivec", "trivec"},
                 {"wdg(trivec,bivec) -> 0", "svBtps", "svBtps", "trivec", "bivec"},
                 {"wdg(bivec,trivec) -> 0", "svBtps", "svBtps", "bivec", "trivec"},
                 {"wdg(trivec,vec) -> ps", "svBtps", "svBtps", "trivec", "vec"},
                 {"wdg(vec,trivec) -> ps", "svBtps", "svBtps", "vec", "trivec"},
                 {"wdg(trivec,s) -> trivec", "svBtps", "svBtps", "trivec", "s"},
                 {"wdg(s,trivec) -> trivec", "svBtps", "svBtps", "s", "trivec"},
                 // bivec
                 {"wdg(bivec,bivec) -> ps", "svBtps1", "svBtps2", "bivec", "bivec"},
                 {"wdg(bivec,vec) -> trivec", "svBtps", "svBtps", "bivec", "vec"},
                 {"wdg(vec,bivec) -> trivec", "svBtps", "svBtps", "vec", "bivec"},
                 {"wdg(bivec,s) -> bivec", "svBtps", "svBtps", "bivec", "s"},
                 {"wdg(s,bivec) -> bivec", "svBtps", "svBtps", "s", "bivec"},
                 // vec
                 {"wdg(vec,vec) -> bivec", "svBtps1", "svBtps2", "vec", "vec"},
                 {"wdg(vec,s) -> vec", "svBtps", "svBtps", "vec", "s"},
                 {"wdg(s,vec) -> vec", "svBtps", "svBtps", "s", "vec"},
                 // s
                 {"wdg(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4ds_dot_config()
{
    return {.product_name = "dot",
            .description = "sta4ds inner product",
            .display_name = "inner product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases =
                {// mv
                 {"dot(mv,mv) -> s", "A", "B", "mv", "mv"},
                 {"dot(mv_e,mv_e) -> s", "A_even", "B_even", "mv_e", "mv_e"},
                 {"dot(mv_u,mv_u) -> s", "A_odd", "B_odd", "mv_u", "mv_u"},
                 // ps
                 {"dot(ps,ps) -> s", "svBtps1", "svBtps2", "ps", "ps"},
                 // trivec
                 {"dot(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec", "trivec"},
                 // bivec
                 {"dot(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec"},
                 // vec
                 {"dot(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec"},
                 // s
                 {"dot(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4ds_l_contract_config()
{
    return {
        .product_name = "l_contract",
        .description = "sta4ds left contraction",
        .display_name = "left contraction",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {// mv
             {"l_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"l_contract(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
             {"l_contract(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
             {"l_contract(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
             {"l_contract(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
             {"l_contract(mv,ps) -> mv", "M", "svBtps", "mv", "ps"},
             {"l_contract(ps,mv) -> mv", "svBtps", "M", "ps", "mv"},
             {"l_contract(mv,trivec) -> mv", "M", "svBtps", "mv", "trivec"},
             {"l_contract(trivec,mv) -> mv", "svBtps", "M", "trivec", "mv"},
             {"l_contract(mv,bivec) -> mv", "M", "svBtps", "mv", "bivec"},
             {"l_contract(bivec,mv) -> mv", "svBtps", "M", "bivec", "mv"},
             {"l_contract(mv,vec) -> mv", "M", "svBtps", "mv", "vec"},
             {"l_contract(vec,mv) -> mv", "svBtps", "M", "vec", "mv"},
             {"l_contract(mv,s) -> mv", "M", "svBtps", "mv", "s"},
             {"l_contract(s,mv) -> mv", "svBtps", "M", "s", "mv"},
             // mv_e
             {"l_contract(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
             {"l_contract(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
             {"l_contract(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
             {"l_contract(mv_e,ps) -> mv_e", "A_even", "svBtps", "mv_e", "ps"},
             {"l_contract(ps,mv_e) -> mv_e", "svBtps", "B_even", "ps", "mv_e"},
             {"l_contract(mv_e,trivec) -> mv_u", "M_even", "svBtps", "mv_e", "trivec"},
             {"l_contract(trivec,mv_e) -> mv_u", "svBtps", "M_even", "trivec", "mv_e"},
             {"l_contract(mv_e,bivec) -> mv_e", "M_even", "svBtps", "mv_e", "bivec"},
             {"l_contract(bivec,mv_e) -> mv_e", "svBtps", "M_even", "bivec", "mv_e"},
             {"l_contract(mv_e,vec) -> mv_u", "A_even", "svBtps", "mv_e", "vec"},
             {"l_contract(vec,mv_e) -> mv_u", "svBtps", "B_even", "vec", "mv_e"},
             {"l_contract(mv_e,s) -> mv_e", "M_even", "svBtps", "mv_e", "s"},
             {"l_contract(s,mv_e) -> mv_e", "svBtps", "M_even", "s", "mv_e"},
             // mv_u
             {"l_contract(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"l_contract(mv_u,ps) -> mv_u", "A_odd", "svBtps", "mv_u", "ps"},
             {"l_contract(ps,mv_u) -> 0", "svBtps", "B_odd", "ps", "mv_u"},
             {"l_contract(mv_u,trivec) -> mv_e", "M_odd", "svBtps", "mv_u", "trivec"},
             {"l_contract(trivec,mv_u) -> mv_e", "svBtps", "M_odd", "trivec", "mv_u"},
             {"l_contract(mv_u,bivec) -> mv_u", "M_odd", "svBtps", "mv_u", "bivec"},
             {"l_contract(bivec,mv_u) -> mv_u", "svBtps", "M_odd", "bivec", "mv_u"},
             {"l_contract(mv_u,vec) -> mv_e", "M_odd", "svBtps", "mv_u", "vec"},
             {"l_contract(vec,mv_u) -> mv_e", "svBtps", "M_odd", "vec", "mv_u"},
             {"l_contract(mv_u,s) -> 0", "M_odd", "svBtps", "mv_u", "s"},
             {"l_contract(s,mv_u) -> mv_u", "svBtps", "M_odd", "s", "mv_u"},
             // ps
             {"l_contract(ps,ps) -> s", "svBtps1", "svBtps2", "ps", "ps"},
             {"l_contract(ps,trivec) -> 0", "svBtps", "svBtps", "ps", "trivec"},
             {"l_contract(trivec,ps) -> vec", "svBtps", "svBtps", "trivec", "ps"},
             {"l_contract(ps,bivec) -> 0", "svBtps", "svBtps", "ps", "bivec"},
             {"l_contract(bivec,ps) -> bivec", "svBtps", "svBtps", "bivec", "ps"},
             {"l_contract(ps,vec) -> 0", "svBtps", "svBtps", "ps", "vec"},
             {"l_contract(vec,ps) -> trivec", "svBtps", "svBtps", "vec", "ps"},
             {"l_contract(ps,s) -> 0", "svBtps", "svBtps", "ps", "s"},
             {"l_contract(s,ps) -> ps", "svBtps", "svBtps", "s", "ps"},
             // trivec
             {"l_contract(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec", "trivec"},
             {"l_contract(trivec,bivec) -> 0", "svBtps", "svBtps", "trivec", "bivec"},
             {"l_contract(bivec,trivec) -> vec", "svBtps", "svBtps", "bivec", "trivec"},
             {"l_contract(trivec,vec) -> 0", "svBtps", "svBtps", "trivec", "vec"},
             {"l_contract(vec,trivec) -> bivec", "svBtps", "svBtps", "vec", "trivec"},
             {"l_contract(trivec,s) -> 0", "svBtps", "svBtps", "trivec", "s"},
             {"l_contract(s,trivec) -> trivec", "svBtps", "svBtps", "s", "trivec"},
             // bivec
             {"l_contract(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec"},
             {"l_contract(bivec,vec) -> 0", "svBtps", "svBtps", "bivec", "vec"},
             {"l_contract(vec,bivec) -> vec", "svBtps", "svBtps", "vec", "bivec"},
             {"l_contract(bivec,s) -> 0", "svBtps", "svBtps", "bivec", "s"},
             {"l_contract(s,bivec) -> bivec", "svBtps", "svBtps", "s", "bivec"},
             // vec
             {"l_contract(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec"},
             {"l_contract(vec,s) -> 0", "svBtps", "svBtps", "vec", "s"},
             {"l_contract(s,vec) -> vec", "svBtps", "svBtps", "s", "vec"},
             // s
             {"l_contract(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_sta4ds_r_contract_config()
{
    return {
        .product_name = "r_contract",
        .description = "sta4ds right contraction",
        .display_name = "right contraction",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {// mv
             {"r_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"r_contract(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
             {"r_contract(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
             {"r_contract(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
             {"r_contract(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
             {"r_contract(mv,ps) -> mv", "M", "svBtps", "mv", "ps"},
             {"r_contract(ps,mv) -> mv", "svBtps", "M", "ps", "mv"},
             {"r_contract(mv,trivec) -> mv", "M", "svBtps", "mv", "trivec"},
             {"r_contract(trivec,mv) -> mv", "svBtps", "M", "trivec", "mv"},
             {"r_contract(mv,bivec) -> mv", "M", "svBtps", "mv", "bivec"},
             {"r_contract(bivec,mv) -> mv", "svBtps", "M", "bivec", "mv"},
             {"r_contract(mv,vec) -> mv", "M", "svBtps", "mv", "vec"},
             {"r_contract(vec,mv) -> mv", "svBtps", "M", "vec", "mv"},
             {"r_contract(mv,s) -> mv", "M", "svBtps", "mv", "s"},
             {"r_contract(s,mv) -> mv", "svBtps", "M", "s", "mv"},
             // mv_e
             {"r_contract(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
             {"r_contract(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
             {"r_contract(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
             {"r_contract(mv_e,ps) -> mv_e", "A_even", "svBtps", "mv_e", "ps"},
             {"r_contract(ps,mv_e) -> mv_e", "svBtps", "B_even", "ps", "mv_e"},
             {"r_contract(mv_e,trivec) -> mv_u", "M_even", "svBtps", "mv_e", "trivec"},
             {"r_contract(trivec,mv_e) -> mv_u", "svBtps", "M_even", "trivec", "mv_e"},
             {"r_contract(mv_e,bivec) -> mv_e", "M_even", "svBtps", "mv_e", "bivec"},
             {"r_contract(bivec,mv_e) -> mv_e", "svBtps", "M_even", "bivec", "mv_e"},
             {"r_contract(mv_e,vec) -> mv_u", "A_even", "svBtps", "mv_e", "vec"},
             {"r_contract(vec,mv_e) -> mv_u", "svBtps", "B_even", "vec", "mv_e"},
             {"r_contract(mv_e,s) -> mv_e", "M_even", "svBtps", "mv_e", "s"},
             {"r_contract(s,mv_e) -> mv_e", "svBtps", "M_even", "s", "mv_e"},
             // mv_u
             {"r_contract(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"r_contract(mv_u,ps) -> 0", "A_odd", "svBtps", "mv_u", "ps"},
             {"r_contract(ps,mv_u) -> mv_u", "svBtps", "B_odd", "ps", "mv_u"},
             {"r_contract(mv_u,trivec) -> mv_e", "M_odd", "svBtps", "mv_u", "trivec"},
             {"r_contract(trivec,mv_u) -> mv_e", "svBtps", "M_odd", "trivec", "mv_u"},
             {"r_contract(mv_u,bivec) -> mv_u", "M_odd", "svBtps", "mv_u", "bivec"},
             {"r_contract(bivec,mv_u) -> mv_u", "svBtps", "M_odd", "bivec", "mv_u"},
             {"r_contract(mv_u,vec) -> mv_e", "M_odd", "svBtps", "mv_u", "vec"},
             {"r_contract(vec,mv_u) -> mv_e", "svBtps", "M_odd", "vec", "mv_u"},
             {"r_contract(mv_u,s) -> mv_u", "M_odd", "svBtps", "mv_u", "s"},
             {"r_contract(s,mv_u) -> 0", "svBtps", "M_odd", "s", "mv_u"},
             // ps
             {"r_contract(ps,ps) -> s", "svBtps1", "svBtps2", "ps", "ps"},
             {"r_contract(ps,trivec) -> vec", "svBtps", "svBtps", "ps", "trivec"},
             {"r_contract(trivec,ps) -> 0", "svBtps", "svBtps", "trivec", "ps"},
             {"r_contract(ps,bivec) -> bivec", "svBtps", "svBtps", "ps", "bivec"},
             {"r_contract(bivec,ps) -> 0", "svBtps", "svBtps", "bivec", "ps"},
             {"r_contract(ps,vec) -> trivec", "svBtps", "svBtps", "ps", "vec"},
             {"r_contract(vec,ps) -> 0", "svBtps", "svBtps", "vec", "ps"},
             {"r_contract(ps,s) -> ps", "svBtps", "svBtps", "ps", "s"},
             {"r_contract(s,ps) -> 0", "svBtps", "svBtps", "s", "ps"},
             // trivec
             {"r_contract(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec", "trivec"},
             {"r_contract(trivec,bivec) -> vec", "svBtps", "svBtps", "trivec", "bivec"},
             {"r_contract(bivec,trivec) -> 0", "svBtps", "svBtps", "bivec", "trivec"},
             {"r_contract(trivec,vec) -> bivec", "svBtps", "svBtps", "trivec", "vec"},
             {"r_contract(vec,trivec) -> 0", "svBtps", "svBtps", "vec", "trivec"},
             {"r_contract(trivec,s) -> trivec", "svBtps", "svBtps", "trivec", "s"},
             {"r_contract(s,trivec) -> 0", "svBtps", "svBtps", "s", "trivec"},
             // bivec
             {"r_contract(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec"},
             {"r_contract(bivec,vec) -> vec", "svBtps", "svBtps", "bivec", "vec"},
             {"r_contract(vec,bivec) -> 0", "svBtps", "svBtps", "vec", "bivec"},
             {"r_contract(bivec,s) -> bivec", "svBtps", "svBtps", "bivec", "s"},
             {"r_contract(s,bivec) -> 0", "svBtps", "svBtps", "s", "bivec"},
             // vec
             {"r_contract(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec"},
             {"r_contract(vec,s) -> vec", "svBtps", "svBtps", "vec", "s"},
             {"r_contract(s,vec) -> 0", "svBtps", "svBtps", "s", "vec"},
             // s
             {"r_contract(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_sta4ds_l_expand_config()
{
    return {
        .product_name = "l_expand",
        .description = "sta4ds left expansion",
        .display_name = "left expansion",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        // Only entries producing vec/bivec/trivec results listed.
        .cases =
            {// ps
             {"l_expand(ps,trivec) -> trivec", "svBtps", "svBtps", "ps", "trivec"},
             {"l_expand(ps,bivec) -> bivec", "svBtps", "svBtps", "ps", "bivec"},
             {"l_expand(ps,vec) -> vec", "svBtps", "svBtps", "ps", "vec"},
             // trivec
             {"l_expand(trivec,bivec) -> trivec", "svBtps", "svBtps", "trivec", "bivec"},
             {"l_expand(trivec,vec) -> bivec", "svBtps", "svBtps", "trivec", "vec"},
             {"l_expand(trivec,s) -> vec", "svBtps", "svBtps", "trivec", "s"},
             // bivec
             {"l_expand(bivec,vec) -> trivec", "svBtps", "svBtps", "bivec", "vec"},
             {"l_expand(bivec,s) -> bivec", "svBtps", "svBtps", "bivec", "s"},
             // vec
             {"l_expand(vec,s) -> trivec", "svBtps", "svBtps", "vec", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_sta4ds_r_expand_config()
{
    return {
        .product_name = "r_expand",
        .description = "sta4ds right expansion",
        .display_name = "right expansion",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        // Only entries producing vec/bivec/trivec results listed.
        .cases =
            {// trivec
             {"r_expand(trivec,ps) -> trivec", "svBtps", "svBtps", "trivec", "ps"},
             // bivec
             {"r_expand(bivec,trivec) -> trivec", "svBtps", "svBtps", "bivec", "trivec"},
             {"r_expand(bivec,ps) -> bivec", "svBtps", "svBtps", "bivec", "ps"},
             // vec
             {"r_expand(vec,bivec) -> trivec", "svBtps", "svBtps", "vec", "bivec"},
             {"r_expand(vec,trivec) -> bivec", "svBtps", "svBtps", "vec", "trivec"},
             {"r_expand(vec,ps) -> vec", "svBtps", "svBtps", "vec", "ps"},
             // s
             {"r_expand(s,vec) -> trivec", "svBtps", "svBtps", "s", "vec"},
             {"r_expand(s,bivec) -> bivec", "svBtps", "svBtps", "s", "bivec"},
             {"r_expand(s,trivec) -> vec", "svBtps", "svBtps", "s", "trivec"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_sta4ds_rgpr_config()
{
    return {.product_name = "rgpr",
            .description = "sta4ds regressive geometric product",
            .display_name = "regressive geometric product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4ds_rcmt_config()
{
    return {.product_name = "rcmt",
            .description = "sta4ds regressive commutator product",
            .display_name = "regressive commutator product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4ds_rwdg_config()
{
    return {
        .product_name = "rwdg",
        .description = "sta4ds regressive wedge product",
        .display_name = "regressive wedge product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {// mv
             {"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"rwdg(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
             {"rwdg(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
             {"rwdg(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
             {"rwdg(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
             {"rwdg(mv,ps) -> mv", "M", "svBtps", "mv", "ps"},
             {"rwdg(ps,mv) -> mv", "svBtps", "M", "ps", "mv"},
             {"rwdg(mv,trivec) -> mv", "M", "svBtps", "mv", "trivec"},
             {"rwdg(trivec,mv) -> mv", "svBtps", "M", "trivec", "mv"},
             {"rwdg(mv,bivec) -> mv", "M", "svBtps", "mv", "bivec"},
             {"rwdg(bivec,mv) -> mv", "svBtps", "M", "bivec", "mv"},
             {"rwdg(mv,vec) -> mv", "M", "svBtps", "mv", "vec"},
             {"rwdg(vec,mv) -> mv", "svBtps", "M", "vec", "mv"},
             {"rwdg(mv,s) -> s", "M", "svBtps", "mv", "s"},
             {"rwdg(s,mv) -> s", "svBtps", "M", "s", "mv"},
             // mv_e
             {"rwdg(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
             {"rwdg(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
             {"rwdg(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
             {"rwdg(mv_e,ps) -> mv_e", "M_even", "svBtps", "mv_e", "ps"},
             {"rwdg(ps,mv_e) -> mv_e", "svBtps", "M_even", "ps", "mv_e"},
             {"rwdg(mv_e,trivec) -> mv_u", "M_even", "svBtps", "mv_e", "trivec"},
             {"rwdg(trivec,mv_e) -> mv_u", "svBtps", "M_even", "trivec", "mv_e"},
             {"rwdg(mv_e,bivec) -> mv_e", "M_even", "svBtps", "mv_e", "bivec"},
             {"rwdg(bivec,mv_e) -> mv_e", "svBtps", "M_even", "bivec", "mv_e"},
             {"rwdg(mv_e,vec) -> vec", "M_even", "svBtps", "mv_e", "vec"},
             {"rwdg(vec,mv_e) -> vec", "svBtps", "M_even", "vec", "mv_e"},
             {"rwdg(mv_e,s) -> s", "M_even", "svBtps", "mv_e", "s"},
             {"rwdg(s,mv_e) -> s", "svBtps", "M_even", "s", "mv_e"},
             // mv_u
             {"rwdg(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"rwdg(mv_u,ps) -> mv_u", "M_odd", "svBtps", "mv_u", "ps"},
             {"rwdg(ps,mv_u) -> mv_u", "svBtps", "M_odd", "ps", "mv_u"},
             {"rwdg(mv_u,trivec) -> mv_e", "M_odd", "svBtps", "mv_u", "trivec"},
             {"rwdg(trivec,mv_u) -> mv_e", "svBtps", "M_odd", "trivec", "mv_u"},
             {"rwdg(mv_u,bivec) -> vec", "M_odd", "svBtps", "mv_u", "bivec"},
             {"rwdg(bivec,mv_u) -> vec", "svBtps", "M_odd", "bivec", "mv_u"},
             {"rwdg(mv_u,vec) -> s", "M_odd", "svBtps", "mv_u", "vec"},
             {"rwdg(vec,mv_u) -> s", "svBtps", "M_odd", "vec", "mv_u"},
             {"rwdg(mv_u,s) -> 0", "M_odd", "svBtps", "mv_u", "s"},
             {"rwdg(s,mv_u) -> 0", "svBtps", "M_odd", "s", "mv_u"},
             // ps
             {"rwdg(ps,ps) -> ps", "svBtps1", "svBtps2", "ps", "ps"},
             {"rwdg(ps,trivec) -> trivec", "svBtps", "svBtps", "ps", "trivec"},
             {"rwdg(trivec,ps) -> trivec", "svBtps", "svBtps", "trivec", "ps"},
             {"rwdg(ps,bivec) -> bivec", "svBtps", "svBtps", "ps", "bivec"},
             {"rwdg(bivec,ps) -> bivec", "svBtps", "svBtps", "bivec", "ps"},
             {"rwdg(ps,vec) -> vec", "svBtps", "svBtps", "ps", "vec"},
             {"rwdg(vec,ps) -> vec", "svBtps", "svBtps", "vec", "ps"},
             {"rwdg(ps,s) -> s", "svBtps", "svBtps", "ps", "s"},
             {"rwdg(s,ps) -> s", "svBtps", "svBtps", "s", "ps"},
             // trivec
             {"rwdg(trivec,trivec) -> bivec", "svBtps1", "svBtps2", "trivec", "trivec"},
             {"rwdg(trivec,bivec) -> vec", "svBtps", "svBtps", "trivec", "bivec"},
             {"rwdg(bivec,trivec) -> vec", "svBtps", "svBtps", "bivec", "trivec"},
             {"rwdg(trivec,vec) -> s", "svBtps", "svBtps", "trivec", "vec"},
             {"rwdg(vec,trivec) -> s", "svBtps", "svBtps", "vec", "trivec"},
             {"rwdg(trivec,s) -> 0", "svBtps", "svBtps", "trivec", "s"},
             {"rwdg(s,trivec) -> 0", "svBtps", "svBtps", "s", "trivec"},
             // bivec
             {"rwdg(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec"},
             {"rwdg(bivec,vec) -> 0", "svBtps", "svBtps", "bivec", "vec"},
             {"rwdg(vec,bivec) -> 0", "svBtps", "svBtps", "vec", "bivec"},
             {"rwdg(bivec,s) -> 0", "svBtps", "svBtps", "bivec", "s"},
             {"rwdg(s,bivec) -> 0", "svBtps", "svBtps", "s", "bivec"},
             // vec
             {"rwdg(vec,vec) -> 0", "svBtps1", "svBtps2", "vec", "vec"},
             {"rwdg(vec,s) -> 0", "svBtps", "svBtps", "vec", "s"},
             {"rwdg(s,vec) -> 0", "svBtps", "svBtps", "s", "vec"},
             // s
             {"rwdg(s,s) -> 0", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_sta4ds_rdot_config()
{
    return {.product_name = "rdot",
            .description = "sta4ds regressive inner product",
            .display_name = "regressive inner product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_sta4ds_sandwich_gpr_config()
{
    return {
        .product_name = "sandwich_gpr",
        .description = "sta4ds sandwich product",
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