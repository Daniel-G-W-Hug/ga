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
 * 1. Define new mvec_coeff in the algebra header file (e.g., ga_prdxpr_ega2d.hpp)
 * 2. Add mapping in create_[algebra]_algebra_data() coefficients map
 * 3. Use the coefficient name in .cases arrays
 *
 * Example extension for new coefficient "C":
 * Header: const mvec_coeff mv2d_coeff_C = {"C.c0", "C.c1", "C.c2", "C.c3"};
 * Config: ega2d.coefficients = {..."C", mv2d_coeff_C}...
 * Usage:  {"mv * mv -> mv", "C", "B", "mv", "mv"}
 *
 * AVAILABLE COEFFICIENTS (EGA2D):
 * - General: A, B, M, M1, M2 (full multivectors)
 * - Even: A_even, B_even, M_even, R_even, R_rev_even (scalar + pseudoscalar)
 * - Symmetric: svps, svps1, svps2 (scalar, vector, pseudoscalar patterns)
 *
 * AVAILABLE FILTERS (EGA2D):
 * - s (scalar), vec (vector), ps (pseudoscalar), mv_e (even), mv (multivector)
 *
 * NAMING CONVENTIONS:
 * - Use lowercase for all type names: s, vec, ps, mv, mv_e
 * - Result types should match geometric algebra conventions
 * - Zero results: Use "0" as result type for operations that yield zero
 */

#include "ga_prdxpr_ega2d_config.hpp"

namespace configurable {

AlgebraData create_ega2d_algebra_data()
{
    AlgebraData ega2d;
    ega2d.name = "ega2d";
    ega2d.basis = mv2d_basis; // Use existing basis
    ega2d.dimension = 2;

    // Map coefficient names to existing coefficient objects
    ega2d.coefficients = {
        {"A", mv2d_coeff_A},           {"B", mv2d_coeff_B},
        {"M", mv2d_coeff_M},           {"M1", mv2d_coeff_M1},
        {"M2", mv2d_coeff_M2},         {"A_even", mv2d_coeff_A_even},
        {"B_even", mv2d_coeff_B_even}, {"M_even", mv2d_coeff_M_even},
        {"R_even", mv2d_coeff_R_even}, {"R_rev_even", mv2d_coeff_R_rev_even},
        {"svps", mv2d_coeff_svps},     {"svps1", mv2d_coeff_svps1},
        {"svps2", mv2d_coeff_svps2}};

    // Map filter names to existing filter enums
    ega2d.filters_2d = {{"s", filter_2d::s},
                        {"vec", filter_2d::vec},
                        {"ps", filter_2d::ps},
                        {"mv_e", filter_2d::mv_e},
                        {"mv", filter_2d::mv}};

    return ega2d;
}

ProductConfig get_ega2d_gpr_config()
{
    return {.product_name = "gpr",
            .description = "ega2d geometric product",
            .display_name = "geometric product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv * mv -> mv", "A", "B", "mv", "mv"},
                      {"mv * mv_e -> mv", "A", "B_even", "mv", "mv_e"},
                      {"mv_e * mv -> mv", "A_even", "B", "mv_e", "mv"},
                      {"mv * ps -> mv", "A", "svps", "mv", "ps"},
                      {"ps * mv -> mv", "svps", "B", "ps", "mv"},
                      {"mv * vec -> mv", "A", "svps", "mv", "vec"},
                      {"vec * mv -> mv", "svps", "B", "vec", "mv"},
                      {"mv * s -> mv", "A", "svps", "mv", "s"},
                      {"s * mv -> mv", "svps", "B", "s", "mv"},
                      {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                      {"mv_e * ps -> mv_e", "A_even", "svps", "mv_e", "ps"},
                      {"ps * mv_e -> mv_e", "svps", "B_even", "ps", "mv_e"},
                      {"mv_e * vec -> vec", "A_even", "svps", "mv_e", "vec"},
                      {"vec * mv_e -> vec", "svps", "B_even", "vec", "mv_e"},
                      {"mv_e * s -> mv_e", "A_even", "svps", "mv_e", "s"},
                      {"s * mv_e -> mv_e", "svps", "B_even", "s", "mv_e"},
                      {"ps * ps -> s", "svps1", "svps2", "ps", "ps"},
                      {"ps * vec -> vec", "svps", "svps", "ps", "vec"},
                      {"vec * ps -> vec", "svps", "svps", "vec", "ps"},
                      {"ps * s -> ps", "svps", "svps", "ps", "s"},
                      {"s * ps -> ps", "svps", "svps", "s", "ps"},
                      {"vec * vec -> mv_e", "svps1", "svps2", "vec", "vec"},
                      {"vec * s -> vec", "svps", "svps", "vec", "s"},
                      {"s * vec -> vec", "svps", "svps", "s", "vec"},
                      {"s * s -> s", "svps1", "svps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_ega2d_cmt_config()
{
    return {.product_name = "cmt",
            .description = "ega2d commutator product",
            .display_name = "commutator product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"cmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"cmt(ps,vec) -> vec", "svps", "svps", "ps", "vec"},
                      {"cmt(vec,ps) -> vec", "svps", "svps", "vec", "ps"},
                      {"cmt(vec,vec) -> ps", "svps1", "svps2", "vec", "vec"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_ega2d_wdg_config()
{
    return {.product_name = "wdg",
            .description = "ega2d wedge product",
            .display_name = "wedge product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv ^ mv -> mv", "A", "B", "mv", "mv"},
                      {"mv ^ mv_e -> mv", "A", "B_even", "mv", "mv_e"},
                      {"mv_e ^ mv -> mv", "A_even", "B", "mv_e", "mv"},
                      {"mv ^ ps -> ps", "M", "svps", "mv", "ps"},
                      {"ps ^ mv -> ps", "svps", "M", "ps", "mv"},
                      {"mv ^ vec -> mv", "M", "svps", "mv", "vec"},
                      {"vec ^ mv -> mv", "svps", "M", "vec", "mv"},
                      {"mv ^ s -> mv", "M", "svps", "mv", "s"},
                      {"s ^ mv -> mv", "svps", "M", "s", "mv"},
                      {"mv_e ^ mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                      {"mv_e ^ ps -> ps", "M_even", "svps", "mv_e", "ps"},
                      {"ps ^ mv_e -> ps", "svps", "M_even", "ps", "mv_e"},
                      {"mv_e ^ vec -> vec", "M_even", "svps", "mv_e", "vec"},
                      {"vec ^ mv_e -> vec", "svps", "M_even", "vec", "mv_e"},
                      {"mv_e ^ s -> mv_e", "M_even", "svps", "mv_e", "s"},
                      {"s ^ mv_e -> mv_e", "svps", "M_even", "s", "mv_e"},
                      {"ps ^ ps -> 0", "svps1", "svps2", "ps", "ps"},
                      {"ps ^ vec -> 0", "svps", "svps", "ps", "vec"},
                      {"vec ^ ps -> 0", "svps", "svps", "vec", "ps"},
                      {"ps ^ s -> ps", "svps", "svps", "ps", "s"},
                      {"s ^ ps -> ps", "svps", "svps", "s", "ps"},
                      {"vec ^ vec -> ps", "svps1", "svps2", "vec", "vec"},
                      {"vec ^ s -> vec", "svps", "svps", "vec", "s"},
                      {"s ^ vec -> vec", "svps", "svps", "s", "vec"},
                      {"s ^ s -> s", "svps1", "svps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_ega2d_dot_config()
{
    return {.product_name = "dot",
            .description = "ega2d inner product",
            .display_name = "inner product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"dot(mv,mv) -> s", "A", "B", "mv", "mv"},
                      {"dot(ps,ps) -> s", "svps1", "svps2", "ps", "ps"},
                      {"dot(vec,vec) -> s", "svps1", "svps2", "vec", "vec"},
                      {"dot(s,s) -> s", "svps1", "svps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_ega2d_dot_alt_config()
{
    return {.product_name = "dot (alternative)",
            .description = "ega2d inner product (defined via regressive outer product)",
            .display_name = "inner product (alternative definition, only main diagonale "
                            "valid, i.e. for arguments of same grade)",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_ega2d_rwdg_config()
{
    return {.product_name = "rwdg",
            .description = "ega2d regressive wedge product",
            .display_name = "regressive wedge product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"rwdg(ps,ps) -> ps", "svps1", "svps2", "ps", "ps"},
                      {"rwdg(ps,vec) -> vec", "svps", "svps", "ps", "vec"},
                      {"rwdg(vec,ps) -> vec", "svps", "svps", "vec", "ps"},
                      {"rwdg(ps,s) -> s", "svps", "svps", "ps", "s"},
                      {"rwdg(s,ps) -> s", "svps", "svps", "s", "ps"},
                      {"rwdg(vec,vec) -> s", "svps1", "svps2", "vec", "vec"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_ega2d_lcontract_config()
{
    return {.product_name = "lcontract",
            .description = "ega2d left contraction",
            .display_name = "left contraction",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv << mv -> mv", "A", "B", "mv", "mv"},
                      {"mv << mv_e -> mv", "A", "B_even", "mv", "mv_e"},
                      {"mv_e << mv -> mv", "A_even", "B", "mv_e", "mv"},
                      {"mv << ps -> mv", "M", "svps", "mv", "ps"},
                      {"ps << mv -> s", "svps", "M", "ps", "mv"},
                      {"mv << vec -> mv", "M", "svps", "mv", "vec"},
                      {"vec << mv -> mv", "svps", "M", "vec", "mv"},
                      {"mv << s -> s", "M", "svps", "mv", "s"},
                      {"s << mv -> mv", "svps", "M", "s", "mv"},
                      {"mv_e << mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                      {"mv_e << ps -> mv_e", "M_even", "svps", "mv_e", "ps"},
                      {"ps << mv_e -> s", "svps", "M_even", "ps", "mv_e"},
                      {"mv_e << vec -> vec", "M_even", "svps", "mv_e", "vec"},
                      {"vec << mv_e -> vec", "svps", "M_even", "vec", "mv_e"},
                      {"mv_e << s -> s", "M_even", "svps", "mv_e", "s"},
                      {"s << mv_e -> mv_e", "svps", "M_even", "s", "mv_e"},
                      {"ps << ps -> s", "svps1", "svps2", "ps", "ps"},
                      {"ps << vec -> 0", "svps", "svps", "ps", "vec"},
                      {"vec << ps -> vec", "svps", "svps", "vec", "ps"},
                      {"ps << s -> 0", "svps", "svps", "ps", "s"},
                      {"s << ps -> ps", "svps", "svps", "s", "ps"},
                      {"vec << vec -> s", "svps1", "svps2", "vec", "vec"},
                      {"vec << s -> 0", "svps", "svps", "vec", "s"},
                      {"s << vec -> vec", "svps", "svps", "s", "vec"},
                      {"s << s -> s", "svps1", "svps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_ega2d_rcontract_config()
{
    return {.product_name = "rcontract",
            .description = "ega2d right contraction",
            .display_name = "right contraction",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv >> mv -> mv", "A", "B", "mv", "mv"},
                      {"mv >> mv_e -> mv", "A", "B_even", "mv", "mv_e"},
                      {"mv_e >> mv -> mv", "A_even", "B", "mv_e", "mv"},
                      {"mv >> ps -> s", "M", "svps", "mv", "ps"},
                      {"ps >> mv -> mv", "svps", "M", "ps", "mv"},
                      {"mv >> vec -> mv", "M", "svps", "mv", "vec"},
                      {"vec >> mv -> mv", "svps", "M", "vec", "mv"},
                      {"mv >> s -> mv", "M", "svps", "mv", "s"},
                      {"s >> mv -> s", "svps", "M", "s", "mv"},
                      {"mv_e >> mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                      {"mv_e >> ps -> s", "M_even", "svps", "mv_e", "ps"},
                      {"ps >> mv_e -> mv_e", "svps", "M_even", "ps", "mv_e"},
                      {"mv_e >> vec -> vec", "M_even", "svps", "mv_e", "vec"},
                      {"vec >> mv_e -> vec", "svps", "M_even", "vec", "mv_e"},
                      {"mv_e >> s -> mv_e", "M_even", "svps", "mv_e", "s"},
                      {"s >> mv_e -> s", "svps", "M_even", "s", "mv_e"},
                      {"ps >> ps -> s", "svps1", "svps2", "ps", "ps"},
                      {"ps >> vec -> vec", "svps", "svps", "ps", "vec"},
                      {"vec >> ps -> 0", "svps", "svps", "vec", "ps"},
                      {"ps >> s -> ps", "svps", "svps", "ps", "s"},
                      {"s >> ps -> 0", "svps", "svps", "s", "ps"},
                      {"vec >> vec -> s", "svps1", "svps2", "vec", "vec"},
                      {"vec >> s -> vec", "svps", "svps", "vec", "s"},
                      {"s >> vec -> 0", "svps", "svps", "s", "vec"},
                      {"s >> s -> s", "svps1", "svps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_ega2d_lexpand_config()
{
    return {
        .product_name = "lexpand",
        .description = "ega2d left expansion",
        .display_name = "left expansion",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {
                // Expansions show basis product tables but no specific cases for ega2d
            },
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_ega2d_rexpand_config()
{
    return {
        .product_name = "rexpand",
        .description = "ega2d right expansion",
        .display_name = "right expansion",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {
                // Expansions show basis product tables but no specific cases for ega2d
            },
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_ega2d_sandwich_gpr_config()
{
    return {
        .product_name = "sandwich_gpr",
        .description = "ega2d sandwich product",
        .display_name = "sandwich product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter", is_two_step, "intermediate"}
        .cases =
            {// Single case that triggers sandwich product behavior - reference only does
             // vec case
             {"dummy", "dummy", "dummy", "dummy", "dummy", true, "vec_tmp"}},
        .is_sandwich_product = true,
        .uses_brace_switch = true, // true is important for sandwich products
        .show_basis_table = true   // Reference shows basis table for sandwich product
    };
}

} // namespace configurable