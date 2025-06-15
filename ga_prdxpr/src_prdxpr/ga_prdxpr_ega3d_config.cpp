// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_ega3d_config.hpp"

namespace configurable {

AlgebraData create_ega3d_algebra_data()
{
    AlgebraData ega3d;
    ega3d.name = "ega3d";
    ega3d.basis = mv3d_basis; // Use existing basis
    ega3d.dimension = 3;

    // Map coefficient names to existing coefficient objects from header
    ega3d.coefficients = {{"A", mv3d_coeff_A},
                          {"B", mv3d_coeff_B},
                          {"M", mv3d_coeff_M},
                          {"M1", mv3d_coeff_M1},
                          {"M2", mv3d_coeff_M2},
                          {"A_even", mv3d_coeff_A_even},
                          {"B_even", mv3d_coeff_B_even},
                          {"A_odd", mv3d_coeff_A_odd},
                          {"B_odd", mv3d_coeff_B_odd},
                          {"M_even", mv3d_coeff_M_even},
                          {"M_odd", mv3d_coeff_M_odd},
                          {"R_even", mv3d_coeff_R_even},
                          {"R_rev_even", mv3d_coeff_R_rev_even},
                          {"svBps", mv3d_coeff_svBps},
                          {"svBps1", mv3d_coeff_svBps1},
                          {"svBps2", mv3d_coeff_svBps2}};

    // Map filter names to existing filter enums
    ega3d.filters_3d = {{"s", filter_3d::s},         {"vec", filter_3d::vec},
                        {"bivec", filter_3d::bivec}, {"ps", filter_3d::ps},
                        {"mv_e", filter_3d::mv_e},   {"mv_u", filter_3d::mv_u},
                        {"mv", filter_3d::mv}};

    return ega3d;
}

ProductConfig get_ega3d_gpr_config()
{
    return {.product_name = "gpr",
            .description = "ega3d geometric product",
            .display_name = "geometric product",
            .cases = {{"mv * mv -> mv", "A", "B", "mv", "mv"},
                      {"mv * mv_e -> mv", "A", "B_even", "mv", "mv_e"},
                      {"mv_e * mv -> mv", "A_even", "B", "mv_e", "mv"},
                      {"mv * mv_u -> mv", "A", "B_odd", "mv", "mv_u"},
                      {"mv_u * mv -> mv", "A_odd", "B", "mv_u", "mv"},
                      {"mv * ps -> mv", "A", "svBps", "mv", "ps"},
                      {"ps * mv -> mv", "svBps", "B", "ps", "mv"},
                      {"mv * bivec -> mv", "M", "svBps", "mv", "bivec"},
                      {"bivec * mv -> mv", "svBps", "M", "bivec", "mv"},
                      {"mv * vec -> mv", "M", "svBps", "mv", "vec"},
                      {"vec * mv -> mv", "svBps", "M", "vec", "mv"},
                      {"mv * s -> mv", "M", "svBps", "mv", "s"},
                      {"s * mv -> mv", "svBps", "M", "s", "mv"},
                      {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                      {"mv_u * mv_u -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
                      {"mv_e * mv_u -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
                      {"mv_u * mv_e -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
                      {"mv_e * ps -> mv_u", "M_even", "svBps", "mv_e", "ps"},
                      {"ps * mv_e -> mv_u", "svBps", "M_even", "ps", "mv_e"},
                      {"mv_u * ps -> mv_e", "M_odd", "svBps", "mv_u", "ps"},
                      {"ps * mv_u -> mv_e", "svBps", "M_odd", "ps", "mv_u"},
                      {"mv_e * bivec -> mv_e", "M_even", "svBps", "mv_e", "bivec"},
                      {"bivec * mv_e -> mv_e", "svBps", "M_even", "bivec", "mv_e"},
                      {"mv_u * bivec -> mv_u", "M_odd", "svBps", "mv_u", "bivec"},
                      {"bivec * mv_u -> mv_u", "svBps", "M_odd", "bivec", "mv_u"},
                      {"mv_e * vec -> mv_u", "M_even", "svBps", "mv_e", "vec"},
                      {"vec * mv_e -> mv_u", "svBps", "M_even", "vec", "mv_e"},
                      {"mv_u * vec -> mv_e", "M_odd", "svBps", "mv_u", "vec"},
                      {"vec * mv_u -> mv_e", "svBps", "M_odd", "vec", "mv_u"},
                      {"mv_e * s -> mv_e", "M_even", "svBps", "mv_e", "s"},
                      {"s * mv_e -> mv_e", "svBps", "M_even", "s", "mv_e"},
                      {"mv_u * s -> mv_u", "M_odd", "svBps", "mv_u", "s"},
                      {"s * mv_u -> mv_u", "svBps", "M_odd", "s", "mv_u"},
                      {"ps * ps -> s", "svBps1", "svBps2", "ps", "ps"},
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
                      {"s * s -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_cmt_config()
{
    return {.product_name = "cmt",
            .description = "ega3d commutator product",
            .display_name = "commutator product",
            .cases = {{"cmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"cmt(bivec,bivec) -> bivec", "svBps1", "svBps2", "bivec", "bivec"},
                      {"cmt(bivec,vec) -> vec", "svBps", "svBps", "bivec", "vec"},
                      {"cmt(vec,bivec) -> vec", "svBps", "svBps", "vec", "bivec"},
                      {"cmt(vec,vec) -> bivec", "svBps1", "svBps2", "vec", "vec"}}};
}

ProductConfig get_ega3d_wdg_config()
{
    return {.product_name = "wdg",
            .description = "ega3d wedge product",
            .display_name = "wedge product",
            .cases = {{"mv ^ mv -> mv", "A", "B", "mv", "mv"},
                      {"mv ^ bivec -> mv", "M", "svBps", "mv", "bivec"},
                      {"bivec ^ mv -> mv", "svBps", "M", "bivec", "mv"},
                      {"mv ^ vec -> mv", "M", "svBps", "mv", "vec"},
                      {"vec ^ mv -> mv", "svBps", "M", "vec", "mv"},
                      {"mv ^ s -> mv", "M", "svBps", "mv", "s"},
                      {"s ^ mv -> mv", "svBps", "M", "s", "mv"},
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
                      {"s ^ s -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_dot_config()
{
    return {.product_name = "dot",
            .description = "ega3d inner product",
            .display_name = "inner product",
            .cases = {{"dot(mv,mv) -> s", "A", "B", "mv", "mv"},
                      {"dot(ps,ps) -> s", "svBps1", "svBps2", "ps", "ps"},
                      {"dot(bivec,bivec) -> s", "svBps1", "svBps2", "bivec", "bivec"},
                      {"dot(vec,vec) -> s", "svBps1", "svBps2", "vec", "vec"},
                      {"dot(s,s) -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_rwdg_config()
{
    return {.product_name = "rwdg",
            .description = "ega3d regressive wedge product",
            .display_name = "regressive wedge product",
            .cases = {{"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"rwdg(ps,ps) -> ps", "svBps1", "svBps2", "ps", "ps"},
                      {"rwdg(ps,bivec) -> bivec", "svBps", "svBps", "ps", "bivec"},
                      {"rwdg(bivec,ps) -> bivec", "svBps", "svBps", "bivec", "ps"},
                      {"rwdg(ps,vec) -> vec", "svBps", "svBps", "ps", "vec"},
                      {"rwdg(vec,ps) -> vec", "svBps", "svBps", "vec", "ps"},
                      {"rwdg(ps,s) -> s", "svBps", "svBps", "ps", "s"},
                      {"rwdg(s,ps) -> s", "svBps", "svBps", "s", "ps"},
                      {"rwdg(bivec,bivec) -> vec", "svBps1", "svBps2", "bivec", "bivec"},
                      {"rwdg(bivec,vec) -> s", "svBps", "svBps", "bivec", "vec"},
                      {"rwdg(vec,bivec) -> s", "svBps", "svBps", "vec", "bivec"}}};
}

ProductConfig get_ega3d_lcontract_config()
{
    return {.product_name = "lcontract",
            .description = "ega3d left contraction",
            .display_name = "left contraction",
            .cases = {{"mv << mv -> mv", "A", "B", "mv", "mv"},
                      {"bivec << mv -> mv", "svBps", "M", "bivec", "mv"},
                      {"vec << mv -> mv", "svBps", "M", "vec", "mv"},
                      {"s << mv -> mv", "svBps", "M", "s", "mv"},
                      {"s << mv_e -> mv_e", "svBps", "M_even", "s", "mv_e"},
                      {"s << mv_u -> mv_u", "svBps", "M_odd", "s", "mv_u"},
                      {"ps << ps -> s", "svBps1", "svBps2", "ps", "ps"},
                      {"ps << bivec -> 0", "svBps", "svBps", "ps", "bivec"},
                      {"bivec << ps -> vec", "svBps", "svBps", "bivec", "ps"},
                      {"ps << vec -> 0", "svBps", "svBps", "ps", "vec"},
                      {"vec << ps -> bivec", "svBps", "svBps", "vec", "ps"},
                      {"ps << s -> 0", "svBps", "svBps", "ps", "s"},
                      {"s << ps -> ps", "svBps", "svBps", "s", "ps"},
                      {"bivec << bivec -> s", "svBps1", "svBps2", "bivec", "bivec"},
                      {"bivec << vec -> 0", "svBps", "svBps", "bivec", "vec"},
                      {"vec << bivec -> vec", "svBps", "svBps", "vec", "bivec"},
                      {"bivec << s -> 0", "svBps", "svBps", "bivec", "s"},
                      {"s << bivec -> bivec", "svBps", "svBps", "s", "bivec"},
                      {"vec << vec -> s", "svBps1", "svBps2", "vec", "vec"},
                      {"vec << s -> 0", "svBps", "svBps", "vec", "s"},
                      {"s << vec -> vec", "svBps", "svBps", "s", "vec"},
                      {"s << s -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_rcontract_config()
{
    return {.product_name = "rcontract",
            .description = "ega3d right contraction",
            .display_name = "right contraction",
            .cases = {{"mv >> mv -> mv", "A", "B", "mv", "mv"},
                      {"mv >> bivec -> mv", "M", "svBps", "mv", "bivec"},
                      {"mv >> vec -> mv", "M", "svBps", "mv", "vec"},
                      {"mv >> s -> mv", "M", "svBps", "mv", "s"},
                      {"mv_e >> s -> mv_e", "M_even", "svBps", "mv_e", "s"},
                      {"mv_u >> s -> mv_u", "M_odd", "svBps", "mv_u", "s"},
                      {"ps >> ps -> s", "svBps1", "svBps2", "ps", "ps"},
                      {"ps >> bivec -> vec", "svBps", "svBps", "ps", "bivec"},
                      {"bivec >> ps -> 0", "svBps", "svBps", "bivec", "ps"},
                      {"ps >> vec -> bivec", "svBps", "svBps", "ps", "vec"},
                      {"vec >> ps -> 0", "svBps", "svBps", "vec", "ps"},
                      {"ps >> s -> ps", "svBps", "svBps", "ps", "s"},
                      {"s >> ps -> 0", "svBps", "svBps", "s", "ps"},
                      {"bivec >> bivec -> s", "svBps1", "svBps2", "bivec", "bivec"},
                      {"bivec >> vec -> vec", "svBps", "svBps", "bivec", "vec"},
                      {"vec >> bivec -> 0", "svBps", "svBps", "vec", "bivec"},
                      {"bivec >> s -> bivec", "svBps", "svBps", "bivec", "s"},
                      {"s >> bivec -> 0", "svBps", "svBps", "s", "bivec"},
                      {"vec >> vec -> s", "svBps1", "svBps2", "vec", "vec"},
                      {"vec >> s -> vec", "svBps", "svBps", "vec", "s"},
                      {"s >> vec -> 0", "svBps", "svBps", "s", "vec"},
                      {"s >> s -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_lexpand_config()
{
    return {.product_name = "lexpand",
            .description = "ega3d left expansion",
            .display_name = "left expansion",
            .cases = {{"lexpand(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"lexpand(bivec,vec) -> bivec", "svBps", "svBps", "bivec", "vec"}}};
}

ProductConfig get_ega3d_rexpand_config()
{
    return {.product_name = "rexpand",
            .description = "ega3d right expansion",
            .display_name = "right expansion",
            .cases = {{"rexpand(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"rexpand(vec,bivec) -> bivec", "svBps", "svBps", "vec", "bivec"}}};
}

ProductConfig get_ega3d_sandwich_gpr_config()
{
    return {.product_name = "sandwich_gpr",
            .description = "ega3d rotor geometric product",
            .display_name = "sandwich product",
            .cases = {{"mv_e * vec -> mv_u", "R_even", "svBps", "mv_e", "vec"},
                      {"mv_e * bivec -> mv_e", "R_even", "svBps", "mv_e", "bivec"}},
            .is_sandwich_product = true};
}

} // namespace configurable