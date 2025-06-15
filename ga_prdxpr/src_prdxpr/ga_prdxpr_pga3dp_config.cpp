// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_pga3dp_config.hpp"

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

ProductConfig get_pga3dp_rbulk_contract_config()
{
    return {
        .product_name = "rbulk_contract",
        .description = "pga3dp right bulk contraction",
        .display_name = "right bulk contraction",
        .cases =
            {{"rbulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"rbulk_contract(ps,ps) -> 0", "svBtps1", "svBtps2", "ps", "ps"},
             {"rbulk_contract(ps,trivec) -> vec", "svBtps", "svBtps", "ps", "trivec"},
             {"rbulk_contract(trivec,ps) -> 0", "svBtps", "svBtps", "trivec", "ps"},
             {"rbulk_contract(ps,bivec) -> bivec", "svBtps", "svBtps", "ps", "bivec"},
             {"rbulk_contract(bivec,ps) -> 0", "svBtps", "svBtps", "bivec", "ps"},
             {"rbulk_contract(ps,vec) -> trivec", "svBtps", "svBtps", "ps", "vec"},
             {"rbulk_contract(vec,ps) -> 0", "svBtps", "svBtps", "vec", "ps"},
             {"rbulk_contract(ps,s) -> ps", "svBtps", "svBtps", "ps", "s"},
             {"rbulk_contract(s,ps) -> 0", "svBtps", "svBtps", "s", "ps"},
             {"rbulk_contract(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec",
              "trivec"},
             {"rbulk_contract(trivec,bivec) -> vec", "svBtps", "svBtps", "trivec",
              "bivec"},
             {"rbulk_contract(bivec,trivec) -> 0", "svBtps", "svBtps", "bivec", "trivec"},
             {"rbulk_contract(trivec,vec) -> bivec", "svBtps", "svBtps", "trivec", "vec"},
             {"rbulk_contract(vec,trivec) -> 0", "svBtps", "svBtps", "vec", "trivec"},
             {"rbulk_contract(trivec,s) -> trivec", "svBtps", "svBtps", "trivec", "s"},
             {"rbulk_contract(s,trivec) -> 0", "svBtps", "svBtps", "s", "trivec"},
             {"rbulk_contract(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec"},
             {"rbulk_contract(bivec,vec) -> vec", "svBtps", "svBtps", "bivec", "vec"},
             {"rbulk_contract(vec,bivec) -> 0", "svBtps", "svBtps", "vec", "bivec"},
             {"rbulk_contract(bivec,s) -> bivec", "svBtps", "svBtps", "bivec", "s"},
             {"rbulk_contract(s,bivec) -> 0", "svBtps", "svBtps", "s", "bivec"},
             {"rbulk_contract(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec"},
             {"rbulk_contract(vec,s) -> vec", "svBtps", "svBtps", "vec", "s"},
             {"rbulk_contract(s,vec) -> 0", "svBtps", "svBtps", "s", "vec"},
             {"rbulk_contract(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}


ProductConfig get_pga3dp_rweight_contract_config()
{
    return {.product_name = "rweight_contract",
            .description = "pga3dp right weight contraction",
            .display_name = "right weight contraction",
            .cases = {}, // For brevity, the weight contractions and expansions show basis
                         // tables but no specific
                         // cases yet
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_lbulk_contract_config()
{
    return {
        .product_name = "lbulk_contract",
        .description = "pga3dp left bulk contraction",
        .display_name = "left bulk contraction",
        .cases =
            {{"lbulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"lbulk_contract(ps,ps) -> 0", "svBtps1", "svBtps2", "ps", "ps"},
             {"lbulk_contract(ps,trivec) -> 0", "svBtps", "svBtps", "ps", "trivec"},
             {"lbulk_contract(trivec,ps) -> vec", "svBtps", "svBtps", "trivec", "ps"},
             {"lbulk_contract(ps,bivec) -> 0", "svBtps", "svBtps", "ps", "bivec"},
             {"lbulk_contract(bivec,ps) -> bivec", "svBtps", "svBtps", "bivec", "ps"},
             {"lbulk_contract(ps,vec) -> 0", "svBtps", "svBtps", "ps", "vec"},
             {"lbulk_contract(vec,ps) -> trivec", "svBtps", "svBtps", "vec", "ps"},
             {"lbulk_contract(ps,s) -> 0", "svBtps", "svBtps", "ps", "s"},
             {"lbulk_contract(s,ps) -> ps", "svBtps", "svBtps", "s", "ps"},
             {"lbulk_contract(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec",
              "trivec"},
             {"lbulk_contract(trivec,bivec) -> 0", "svBtps", "svBtps", "trivec", "bivec"},
             {"lbulk_contract(bivec,trivec) -> vec", "svBtps", "svBtps", "bivec",
              "trivec"},
             {"lbulk_contract(trivec,vec) -> 0", "svBtps", "svBtps", "trivec", "vec"},
             {"lbulk_contract(vec,trivec) -> bivec", "svBtps", "svBtps", "vec", "trivec"},
             {"lbulk_contract(trivec,s) -> 0", "svBtps", "svBtps", "trivec", "s"},
             {"lbulk_contract(s,trivec) -> trivec", "svBtps", "svBtps", "s", "trivec"},
             {"lbulk_contract(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec"},
             {"lbulk_contract(bivec,vec) -> 0", "svBtps", "svBtps", "bivec", "vec"},
             {"lbulk_contract(vec,bivec) -> vec", "svBtps", "svBtps", "vec", "bivec"},
             {"lbulk_contract(bivec,s) -> 0", "svBtps", "svBtps", "bivec", "s"},
             {"lbulk_contract(s,bivec) -> bivec", "svBtps", "svBtps", "s", "bivec"},
             {"lbulk_contract(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec"},
             {"lbulk_contract(vec,s) -> 0", "svBtps", "svBtps", "vec", "s"},
             {"lbulk_contract(s,vec) -> vec", "svBtps", "svBtps", "s", "vec"},
             {"lbulk_contract(s,s) -> s", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_pga3dp_lweight_contract_config()
{
    return {.product_name = "lweight_contract",
            .description = "pga3dp left weight contraction",
            .display_name = "left weight contraction",
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_rbulk_expand_config()
{
    return {.product_name = "rbulk_expand",
            .description = "pga3dp right bulk expansion",
            .display_name = "right bulk expansion",
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_rweight_expand_config()
{
    return {.product_name = "rweight_expand",
            .description = "pga3dp right weight expansion",
            .display_name = "right weight expansion",
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_lbulk_expand_config()
{
    return {.product_name = "lbulk_expand",
            .description = "pga3dp left bulk expansion",
            .display_name = "left bulk expansion",
            .cases = {},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}

ProductConfig get_pga3dp_lweight_expand_config()
{
    return {.product_name = "lweight_expand",
            .description = "pga3dp left weight expansion",
            .display_name = "left weight expansion",
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