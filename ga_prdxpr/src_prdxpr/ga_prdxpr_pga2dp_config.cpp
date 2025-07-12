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
 * - PGA2DP: {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e321"} -> {"A.c0"..."A.c7"}
 *
 * TO EXTEND THE SYSTEM:
 * 1. Define new mvec_coeff in the algebra header file (e.g., ga_prdxpr_pga2dp.hpp)
 * 2. Add mapping in create_[algebra]_algebra_data() coefficients map
 * 3. Use the coefficient name in .cases arrays
 *
 * Example extension for new coefficient "C":
 * Header: const mvec_coeff mv2dp_coeff_C = {"C.c0", "C.c1", "C.c2", "C.c3", "C.c4",
 * "C.c5", "C.c6", "C.c7"}; Config: pga2dp.coefficients = {..."C", mv2dp_coeff_C}...
 * Usage:  {"mv * mv -> mv", "C", "B", "mv", "mv"}
 *
 * AVAILABLE COEFFICIENTS (PGA2DP):
 * - General: A, B, M, M1, M2 (full multivectors)
 * - Even/Odd: A_even, B_even, A_odd, B_odd, M_even, M_odd (grade-filtered)
 * - Motors: R_even, R_rrev_even (even elements for regressive operations)
 * - Symmetric: svBps, svBps1, svBps2 (scalar, vector, bivector, pseudoscalar patterns)
 *
 * AVAILABLE FILTERS (PGA2DP):
 * - s (scalar), vec (vector), bivec (bivector), ps (pseudoscalar)
 * - mv_e (even), mv_u (odd), mv (multivector)
 *
 * NAMING CONVENTIONS:
 * - Use lowercase for all type names: s, vec, bivec, ps, mv, mv_e, mv_u
 * - Result types should match geometric algebra conventions
 * - Zero results: Use "0" as result type for operations that yield zero
 */

#include "ga_prdxpr_pga2dp_config.hpp"
#include "ga_prdxpr_pga2dp.hpp"

namespace configurable {

AlgebraData create_pga2dp_algebra_data()
{
    AlgebraData pga2dp;
    pga2dp.name = "pga2dp";
    pga2dp.basis = mv2dp_basis; // Use existing basis
    pga2dp.dimension = 3;       // PGA2DP: 2D algebra embedded in 3D space

    // Map coefficient names to existing coefficient objects from header
    pga2dp.coefficients = {{"A", mv2dp_coeff_A},
                           {"B", mv2dp_coeff_B},
                           {"M", mv2dp_coeff_M},
                           {"M1", mv2dp_coeff_M1},
                           {"M2", mv2dp_coeff_M2},
                           {"A_even", mv2dp_coeff_A_even},
                           {"B_even", mv2dp_coeff_B_even},
                           {"A_odd", mv2dp_coeff_A_odd},
                           {"B_odd", mv2dp_coeff_B_odd},
                           {"M_even", mv2dp_coeff_M_even},
                           {"M_odd", mv2dp_coeff_M_odd},
                           {"R_even", mv2dp_coeff_R_even},
                           {"R_rrev_even", mv2dp_coeff_R_rrev_even},
                           {"svBps", mv2dp_coeff_svBps},
                           {"svBps1", mv2dp_coeff_svBps1},
                           {"svBps2", mv2dp_coeff_svBps2}};

    // Map filter names to existing filter enums for 3D (PGA2DP uses 3D filters)
    pga2dp.filters_3d = {{"s", filter_3d::s},         {"vec", filter_3d::vec},
                         {"bivec", filter_3d::bivec}, {"ps", filter_3d::ps},
                         {"mv_e", filter_3d::mv_e},   {"mv_u", filter_3d::mv_u},
                         {"mv", filter_3d::mv}};

    return pga2dp;
}

ProductConfig get_pga2dp_gpr_config()
{
    return {.product_name = "gpr",
            .description = "pga2dp geometric product",
            .display_name = "geometric product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv * mv -> mv", "A", "B", "mv", "mv"},
                      {"mv * mv_e -> mv", "A", "B_even", "mv", "mv_e"},
                      {"mv_e * mv -> mv", "A_even", "B", "mv_e", "mv"},
                      {"mv * mv_u -> mv", "A", "B_odd", "mv", "mv_u"},
                      {"mv_u * mv -> mv", "A_odd", "B", "mv_u", "mv"},
                      {"mv * ps -> mv", "A", "svBps", "mv", "ps"},
                      {"ps * mv -> mv", "svBps", "B", "ps", "mv"},
                      {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                      {"mv_u * mv_u -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
                      {"mv_e * mv_u -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
                      {"mv_u * mv_e -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
                      {"mv_e * ps -> mv_u", "A_even", "svBps", "mv_e", "ps"},
                      {"ps * mv_e -> mv_u", "svBps", "B_even", "ps", "mv_e"},
                      {"mv_u * ps -> mv_e", "A_odd", "svBps", "mv_u", "ps"},
                      {"ps * mv_u -> mv_e", "svBps", "B_odd", "ps", "mv_u"},
                      {"mv_e * bivec -> mv_e", "M_even", "svBps", "mv_e", "bivec"},
                      {"bivec * mv_e -> mv_e", "svBps", "M_even", "bivec", "mv_e"},
                      {"mv_u * bivec -> mv_u", "M_odd", "svBps", "mv_u", "bivec"},
                      {"bivec * mv_u -> mv_u", "svBps", "M_odd", "bivec", "mv_u"},
                      {"mv_e * vec -> mv_u", "A_even", "svBps", "mv_e", "vec"},
                      {"vec * mv_e -> mv_u", "svBps", "B_even", "vec", "mv_e"},
                      {"ps * ps -> 0", "svBps1", "svBps2", "ps", "ps"},
                      {"ps * bivec -> vec", "svBps", "svBps", "ps", "bivec"},
                      {"bivec * ps -> vec", "svBps", "svBps", "bivec", "ps"},
                      {"ps * vec -> bivec", "svBps", "svBps", "ps", "vec"},
                      {"vec * ps -> bivec", "svBps", "svBps", "vec", "ps"},
                      {"ps * s -> ps", "svBps", "svBps", "ps", "s"},
                      {"s * ps -> ps", "svBps", "svBps", "s", "ps"},
                      {"bivec * bivec -> mv_e", "svBps1", "svBps2", "bivec", "bivec"},
                      {"bivec * vec -> mv_u", "svBps", "svBps", "bivec", "vec"},
                      {"vec * bivec -> mv_u", "svBps", "svBps", "vec", "bivec"},
                      {"bivec * s -> bivec", "svBps", "svBps", "bivec", "s"},
                      {"s * bivec -> bivec", "svBps", "svBps", "s", "bivec"},
                      {"vec * vec -> mv_e", "svBps1", "svBps2", "vec", "vec"},
                      {"vec * s -> vec", "svBps", "svBps", "vec", "s"},
                      {"s * vec -> vec", "svBps", "svBps", "s", "vec"},
                      {"s * s -> s", "svBps1", "svBps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_cmt_config()
{
    return {.product_name = "cmt",
            .description = "pga2dp commutator product",
            .display_name = "commutator product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"cmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"cmt(bivec,bivec) -> bivec", "svBps1", "svBps2", "bivec", "bivec"},
                      {"cmt(bivec,vec) -> vec", "svBps", "svBps", "bivec", "vec"},
                      {"cmt(vec,bivec) -> vec", "svBps", "svBps", "vec", "bivec"},
                      {"cmt(vec,vec) -> bivec", "svBps1", "svBps2", "vec", "vec"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_wdg_config()
{
    return {.product_name = "wdg",
            .description = "pga2dp wedge product",
            .display_name = "wedge product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv ^ mv -> mv", "A", "B", "mv", "mv"},
                      {"ps ^ ps -> 0", "svBps1", "svBps2", "ps", "ps"},
                      {"ps ^ bivec -> 0", "svBps", "svBps", "ps", "bivec"},
                      {"bivec ^ ps -> 0", "svBps", "svBps", "bivec", "ps"},
                      {"ps ^ vec -> 0", "svBps", "svBps", "ps", "vec"},
                      {"vec ^ ps -> 0", "svBps", "svBps", "vec", "ps"},
                      {"ps ^ s -> ps", "svBps", "svBps", "ps", "s"},
                      {"s ^ ps -> ps", "svBps", "svBps", "s", "ps"},
                      {"bivec ^ bivec -> 0", "svBps1", "svBps2", "bivec", "bivec"},
                      {"bivec ^ vec -> ps", "svBps", "svBps", "bivec", "vec"},
                      {"vec ^ bivec -> ps", "svBps", "svBps", "vec", "bivec"},
                      {"bivec ^ s -> bivec", "svBps", "svBps", "bivec", "s"},
                      {"s ^ bivec -> bivec", "svBps", "svBps", "s", "bivec"},
                      {"vec ^ vec -> bivec", "svBps1", "svBps2", "vec", "vec"},
                      {"vec ^ s -> vec", "svBps", "svBps", "vec", "s"},
                      {"s ^ vec -> vec", "svBps", "svBps", "s", "vec"},
                      {"s ^ s -> s", "svBps1", "svBps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_dot_config()
{
    return {.product_name = "dot",
            .description = "pga2dp inner product",
            .display_name = "inner product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"dot(mv,mv) -> s", "A", "B", "mv", "mv"},
                      {"dot(ps,ps) -> 0", "svBps1", "svBps2", "ps", "ps"},
                      {"dot(bivec,bivec) -> s", "svBps1", "svBps2", "bivec", "bivec"},
                      {"dot(vec,vec) -> s", "svBps1", "svBps2", "vec", "vec"},
                      {"dot(s,s) -> s", "svBps1", "svBps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_rwdg_config()
{
    return {.product_name = "rwdg",
            .description = "pga2dp regressive wedge product",
            .display_name = "regressive wedge product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"rwdg(ps,bivec) -> bivec", "svBps", "svBps", "ps", "bivec"},
                      {"rwdg(bivec,ps) -> bivec", "svBps", "svBps", "bivec", "ps"},
                      {"rwdg(bivec,bivec) -> vec", "svBps1", "svBps2", "bivec", "bivec"},
                      {"rwdg(bivec,vec) -> s", "svBps", "svBps", "bivec", "vec"},
                      {"rwdg(vec,bivec) -> s", "svBps", "svBps", "vec", "bivec"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_rdot_config()
{
    return {.product_name = "rdot",
            .description = "pga2dp regressive inner product",
            .display_name = "regressive inner product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"rdot(ps,ps) -> ps", "svBps1", "svBps2", "ps", "ps"},
                      {"rdot(bivec,bivec) -> ps", "svBps1", "svBps2", "bivec", "bivec"},
                      {"rdot(vec,vec) -> ps", "svBps1", "svBps2", "vec", "vec"},
                      {"rdot(s,s) -> 0", "svBps1", "svBps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_rgpr_config()
{
    return {.product_name = "rgpr",
            .description = "pga2dp regressive geometric product",
            .display_name = "regressive geometric product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"rgpr(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"rgpr(mv_e,mv_e) -> mv_u", "A_even", "B_even", "mv_e", "mv_e"},
                      {"rgpr(mv_u,mv_u) -> mv_u", "A_odd", "B_odd", "mv_u", "mv_u"},
                      {"rgpr(mv_e,mv_u) -> mv_e", "A_even", "B_odd", "mv_e", "mv_u"},
                      {"rgpr(mv_u,mv_e) -> mv_e", "A_odd", "B_even", "mv_u", "mv_e"},
                      {"rgpr(mv_e,bivec) -> mv_u", "M_even", "svBps", "mv_e", "bivec"},
                      {"rgpr(mv_u,bivec) -> mv_e", "M_odd", "svBps", "mv_u", "bivec"},
                      {"rgpr(mv_e,vec) -> mv_e", "M_even", "svBps", "mv_e", "vec"},
                      {"rgpr(mv_u,vec) -> mv_u", "M_odd", "svBps", "mv_u", "vec"},
                      {"rgpr(bivec,bivec) -> mv_u", "svBps1", "svBps2", "bivec", "bivec"},
                      {"rgpr(bivec,vec) -> mv_e", "svBps", "svBps", "bivec", "vec"},
                      {"rgpr(vec,vec) -> mv_u", "svBps1", "svBps2", "vec", "vec"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_rcmt_config()
{
    return {.product_name = "rcmt",
            .description = "pga2dp regressive commutator product",
            .display_name = "regressive commutator product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"rcmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"rcmt(bivec,bivec) -> vec", "svBps1", "svBps2", "bivec", "bivec"},
                      {"rcmt(bivec,vec) -> bivec", "svBps", "svBps", "bivec", "vec"},
                      {"rcmt(vec,bivec) -> bivec", "svBps", "svBps", "vec", "bivec"},
                      {"rcmt(vec,vec) -> vec", "svBps1", "svBps2", "vec", "vec"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_right_bulk_contract_config()
{
    return {
        .product_name = "right_bulk_contract",
        .description = "pga2dp right bulk contraction",
        .display_name = "right bulk contraction",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {{"right_bulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"right_bulk_contract(ps,ps) -> 0", "svBps1", "svBps2", "ps", "ps"},
             {"right_bulk_contract(ps,bivec) -> vec", "svBps", "svBps", "ps", "bivec"},
             {"right_bulk_contract(bivec,ps) -> 0", "svBps", "svBps", "bivec", "ps"},
             {"right_bulk_contract(vec,ps) -> 0", "svBps", "svBps", "vec", "ps"},
             {"right_bulk_contract(ps,s) -> ps", "svBps", "svBps", "ps", "s"},
             {"right_bulk_contract(s,ps) -> 0", "svBps", "svBps", "s", "ps"},
             {"right_bulk_contract(bivec,bivec) -> s", "svBps1", "svBps2", "bivec",
              "bivec"},
             {"right_bulk_contract(bivec,vec) -> vec", "svBps", "svBps", "bivec", "vec"},
             {"right_bulk_contract(vec,bivec) -> 0", "svBps", "svBps", "vec", "bivec"},
             {"right_bulk_contract(bivec,s) -> bivec", "svBps", "svBps", "bivec", "s"},
             {"right_bulk_contract(s,bivec) -> 0", "svBps", "svBps", "s", "bivec"},
             {"right_bulk_contract(vec,vec) -> s", "svBps1", "svBps2", "vec", "vec"},
             {"right_bulk_contract(vec,s) -> vec", "svBps", "svBps", "vec", "s"},
             {"right_bulk_contract(s,vec) -> 0", "svBps", "svBps", "s", "vec"},
             {"right_bulk_contract(s,s) -> s", "svBps1", "svBps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga2dp_right_weight_contract_config()
{
    return {.product_name = "right_weight_contract",
            .description = "pga2dp right weight contraction",
            .display_name = "right weight contraction",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_left_bulk_contract_config()
{
    return {
        .product_name = "left_bulk_contract",
        .description = "pga2dp left bulk contraction",
        .display_name = "left bulk contraction",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {{"left_bulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"left_bulk_contract(ps,ps) -> 0", "svBps1", "svBps2", "ps", "ps"},
             {"left_bulk_contract(ps,bivec) -> 0", "svBps", "svBps", "ps", "bivec"},
             {"left_bulk_contract(bivec,ps) -> vec", "svBps", "svBps", "bivec", "ps"},
             {"left_bulk_contract(ps,vec) -> 0", "svBps", "svBps", "ps", "vec"},
             {"left_bulk_contract(vec,ps) -> bivec", "svBps", "svBps", "vec", "ps"},
             {"left_bulk_contract(ps,s) -> 0", "svBps", "svBps", "ps", "s"},
             {"left_bulk_contract(s,ps) -> ps", "svBps", "svBps", "s", "ps"},
             {"left_bulk_contract(bivec,bivec) -> s", "svBps1", "svBps2", "bivec",
              "bivec"},
             {"left_bulk_contract(bivec,vec) -> 0", "svBps", "svBps", "bivec", "vec"},
             {"left_bulk_contract(vec,bivec) -> vec", "svBps", "svBps", "vec", "bivec"},
             {"left_bulk_contract(bivec,s) -> 0", "svBps", "svBps", "bivec", "s"},
             {"left_bulk_contract(s,bivec) -> bivec", "svBps", "svBps", "s", "bivec"},
             {"left_bulk_contract(vec,vec) -> s", "svBps1", "svBps2", "vec", "vec"},
             {"left_bulk_contract(vec,s) -> 0", "svBps", "svBps", "vec", "s"},
             {"left_bulk_contract(s,vec) -> vec", "svBps", "svBps", "s", "vec"},
             {"left_bulk_contract(s,s) -> s", "svBps1", "svBps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga2dp_left_weight_contract_config()
{
    return {.product_name = "left_weight_contract",
            .description = "pga2dp left weight contraction",
            .display_name = "left weight contraction",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_right_bulk_expand_config()
{
    return {.product_name = "right_bulk_expand",
            .description = "pga2dp right bulk expansion",
            .display_name = "right bulk expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_right_weight_expand_config()
{
    return {.product_name = "right_weight_expand",
            .description = "pga2dp right weight expansion",
            .display_name = "right weight expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_left_bulk_expand_config()
{
    return {.product_name = "left_bulk_expand",
            .description = "pga2dp left bulk expansion",
            .display_name = "left bulk expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_left_weight_expand_config()
{
    return {.product_name = "left_weight_expand",
            .description = "pga2dp left weight expansion",
            .display_name = "left weight expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga2dp_sandwich_rgpr_config()
{
    return {.product_name = "sandwich_rgpr",
            .description = "pga2dp regressive sandwich product",
            .display_name = "regressive sandwich product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter", is_two_step, "intermediate"}
            .cases =
                {// Single case that triggers regressive sandwich product behavior - motor
                 // operations
                 {"dummy", "dummy", "dummy", "dummy", "dummy", true, "vec_tmp"}},
            .is_sandwich_product = true,
            .uses_brace_switch = true,
            .show_basis_table = true};
}

} // namespace configurable