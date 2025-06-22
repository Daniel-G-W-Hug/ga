// Copyright 2024-2025, Daniel Hug. All rights reserved.

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