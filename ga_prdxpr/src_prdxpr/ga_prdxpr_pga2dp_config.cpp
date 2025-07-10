// Copyright 2024-2025, Daniel Hug. All rights reserved.

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

ProductConfig get_pga2dp_right_bulk_contract_config()
{
    return {
        .product_name = "right_bulk_contract",
        .description = "pga2dp right bulk contraction",
        .display_name = "right bulk contraction",
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
            .cases =
                {// Single case that triggers regressive sandwich product behavior - motor
                 // operations
                 {"dummy", "dummy", "dummy", "dummy", "dummy", true, "vec_tmp"}},
            .is_sandwich_product = true,
            .uses_brace_switch = true,
            .show_basis_table = true};
}

} // namespace configurable