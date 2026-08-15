// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "algebras/ga_prdxpr_cga2dc_config.hpp"
#include "algebras/ga_prdxpr_cga2dc.hpp"

////////////////////////////////////////////////////////////////////////////////
// cga2dc algebra configuration
//
// Basis vectors: e1, e2 Euclidean; e3 projective (w); e4 round (u).
// The null pair satisfies e3*e3 = e4*e4 = 0 and e3.e4 = -1, i.e. the metric
// has an off-diagonal coupling — carried by metric_matrix (metric_signature
// holds its diagonal, as required by the validation).
//
// Diagonalization (exact, dyadic coefficients): internal orthogonal basis
// (e1, e2, e-, e+) with signature (+1, +1, -1, +1) and
//
//     e3 = 0.5*(e- - e+)          e- =  e3 + 0.5*e4
//     e4 = e- + e+                e+ = -e3 + 0.5*e4
//
// The normalized point embedding q = x*e1 + y*e2 + e3 + 0.5*(x²+y²)*e4 is
// null: q ⟑ q = x² + y² + 2*(0.5*(x²+y²))*(e3.e4) = 0.
////////////////////////////////////////////////////////////////////////////////

AlgebraConfig get_cga2dc_algebra_config()
{
    // Extract basis prefix from vector basis and validate consistency
    std::string const prefix = extract_basis_prefix(mv2dc_basis_kvec[1]);
    validate_basis_consistency(mv2dc_basis, mv2dc_basis_kvec, prefix, one_str());

    AlgebraConfig config;
    config.basis_vectors = {"e1", "e2", "e3", "e4"};
    config.metric_signature = {1, 1, 0, 0}; // diagonal of metric_matrix
    config.multivector_basis = mv2dc_basis;
    config.scalar_name = one_str();
    config.basis_prefix = prefix;

    config.metric_matrix = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, -1}, {0, 0, -1, 0}};
    config.basis_change = {{1.0, 0.0, 0.0, 0.0},
                           {0.0, 1.0, 0.0, 0.0},
                           {0.0, 0.0, 0.5, -0.5},
                           {0.0, 0.0, 1.0, 1.0}};
    config.basis_change_inv = {{1.0, 0.0, 0.0, 0.0},
                               {0.0, 1.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.5},
                               {0.0, 0.0, -1.0, 0.5}};
    config.diag_signature = {1, 1, -1, 1};

    return config;
}

////////////////////////////////////////////////////////////////////////////////
// Generator-facing algebra data and product configurations
////////////////////////////////////////////////////////////////////////////////

namespace configurable {

AlgebraData create_cga2dc_algebra_data()
{
    AlgebraData cga2dc;
    cga2dc.name = "cga2dc";
    cga2dc.basis = mv2dc_basis; // Use existing basis
    cga2dc.dimension = 4;       // e1, e2: Euclidean; e3: projective; e4: round

    // Map coefficient names to existing coefficient objects
    cga2dc.coefficients = {{"A", mv2dc_coeff_A},
                           {"B", mv2dc_coeff_B},
                           {"M", mv2dc_coeff_M},
                           {"M1", mv2dc_coeff_M1},
                           {"M2", mv2dc_coeff_M2},
                           {"M_even", mv2dc_coeff_M_even},
                           {"M_odd", mv2dc_coeff_M_odd},
                           {"A_even", mv2dc_coeff_A_even},
                           {"B_even", mv2dc_coeff_B_even},
                           {"A_odd", mv2dc_coeff_A_odd},
                           {"B_odd", mv2dc_coeff_B_odd},
                           {"R_even", mv2dc_coeff_R_even},
                           {"R_odd", mv2dc_coeff_R_odd},
                           {"R_rev_even", mv2dc_coeff_R_rev_even},
                           {"R_rev_odd", mv2dc_coeff_R_rev_odd},
                           {"R_rrev_even", mv2dc_coeff_R_rrev_even},
                           {"R_rrev_odd", mv2dc_coeff_R_rrev_odd},
                           {"svBtps", mv2dc_coeff_svBtps},
                           {"svBtps1", mv2dc_coeff_svBtps1},
                           {"svBtps2", mv2dc_coeff_svBtps2}};

    // Map filter names to existing filter enums
    cga2dc.filters_4d = {{"s", filter_4d::s},         {"vec", filter_4d::vec},
                         {"bivec", filter_4d::bivec}, {"trivec", filter_4d::trivec},
                         {"ps", filter_4d::ps},       {"mv_e", filter_4d::mv_e},
                         {"mv_u", filter_4d::mv_u},   {"mv", filter_4d::mv}};

    return cga2dc;
}

ProductConfig get_cga2dc_gpr_config()
{
    return {.product_name = "gpr",
            .description = "cga2dc geometric product",
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

ProductConfig get_cga2dc_wdg_config()
{
    return {.product_name = "wdg",
            .description = "cga2dc wedge product",
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

ProductConfig get_cga2dc_dot_config()
{
    return {.product_name = "dot",
            .description = "cga2dc inner product",
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


ProductConfig get_cga2dc_cmt_config()
{
    return {
        .product_name = "cmt",
        .description = "cga2dc commutator product",
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

ProductConfig get_cga2dc_l_contract_config()
{
    return {
        .product_name = "l_contract",
        .description = "cga2dc left contraction",
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

ProductConfig get_cga2dc_r_contract_config()
{
    return {
        .product_name = "r_contract",
        .description = "cga2dc right contraction",
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

ProductConfig get_cga2dc_l_expand_config()
{
    return {
        .product_name = "l_expand",
        .description = "cga2dc left expansion",
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

ProductConfig get_cga2dc_r_expand_config()
{
    return {
        .product_name = "r_expand",
        .description = "cga2dc right expansion",
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

ProductConfig get_cga2dc_rgpr_config()
{
    return {
        .product_name = "rgpr",
        .description = "cga2dc regressive geometric product",
        .display_name = "regressive geometric product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {// mv
             {"rgpr(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"rgpr(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
             {"rgpr(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
             {"rgpr(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
             {"rgpr(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
             {"rgpr(mv,ps) -> mv", "M", "svBtps", "mv", "ps"},
             {"rgpr(ps,mv) -> mv", "svBtps", "M", "ps", "mv"},
             {"rgpr(mv,trivec) -> mv", "M", "svBtps", "mv", "trivec"},
             {"rgpr(trivec,mv) -> mv", "svBtps", "M", "trivec", "mv"},
             {"rgpr(mv,bivec) -> mv", "M", "svBtps", "mv", "bivec"},
             {"rgpr(bivec,mv) -> mv", "svBtps", "M", "bivec", "mv"},
             {"rgpr(mv,vec) -> mv", "M", "svBtps", "mv", "vec"},
             {"rgpr(vec,mv) -> mv", "svBtps", "M", "vec", "mv"},
             {"rgpr(mv,s) -> mv", "M", "svBtps", "mv", "s"},
             {"rgpr(s,mv) -> mv", "svBtps", "M", "s", "mv"},
             // mv_e
             {"rgpr(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
             {"rgpr(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
             {"rgpr(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
             {"rgpr(mv_e,ps) -> mv_e", "A_even", "svBtps", "mv_e", "ps"},
             {"rgpr(ps,mv_e) -> mv_e", "svBtps", "B_even", "ps", "mv_e"},
             {"rgpr(mv_e,trivec) -> mv_u", "M_even", "svBtps", "mv_e", "trivec"},
             {"rgpr(trivec,mv_e) -> mv_u", "svBtps", "M_even", "trivec", "mv_e"},
             {"rgpr(mv_e,bivec) -> mv_e", "M_even", "svBtps", "mv_e", "bivec"},
             {"rgpr(bivec,mv_e) -> mv_e", "svBtps", "M_even", "bivec", "mv_e"},
             {"rgpr(mv_e,vec) -> mv_u", "M_even", "svBtps", "mv_e", "vec"},
             {"rgpr(vec,mv_e) -> mv_u", "svBtps", "M_even", "vec", "mv_e"},
             {"rgpr(mv_e,s) -> mv_e", "M_even", "svBtps", "mv_e", "s"},
             {"rgpr(s,mv_e) -> mv_e", "svBtps", "M_even", "s", "mv_e"},
             // mv_u
             {"rgpr(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"rgpr(mv_u,ps) -> mv_u", "A_odd", "svBtps", "mv_u", "ps"},
             {"rgpr(ps,mv_u) -> mv_u", "svBtps", "B_odd", "ps", "mv_u"},
             {"rgpr(mv_u,trivec) -> mv_e", "M_odd", "svBtps", "mv_u", "trivec"},
             {"rgpr(trivec,mv_u) -> mv_e", "svBtps", "M_odd", "trivec", "mv_u"},
             {"rgpr(mv_u,bivec) -> mv_u", "M_odd", "svBtps", "mv_u", "bivec"},
             {"rgpr(bivec,mv_u) -> mv_u", "svBtps", "M_odd", "bivec", "mv_u"},
             {"rgpr(mv_u,vec) -> mv_e", "M_odd", "svBtps", "mv_u", "vec"},
             {"rgpr(vec,mv_u) -> mv_e", "svBtps", "M_odd", "vec", "mv_u"},
             {"rgpr(mv_u,s) -> mv_u", "M_odd", "svBtps", "mv_u", "s"},
             {"rgpr(s,mv_u) -> mv_u", "svBtps", "M_odd", "s", "mv_u"},
             // ps
             {"rgpr(ps,ps) -> ps", "svBtps1", "svBtps2", "ps", "ps"},
             {"rgpr(ps,trivec) -> trivec", "svBtps", "svBtps", "ps", "trivec"},
             {"rgpr(trivec,ps) -> trivec", "svBtps", "svBtps", "trivec", "ps"},
             {"rgpr(ps,bivec) -> bivec", "svBtps", "svBtps", "ps", "bivec"},
             {"rgpr(bivec,ps) -> bivec", "svBtps", "svBtps", "bivec", "ps"},
             {"rgpr(ps,vec) -> vec", "svBtps", "svBtps", "ps", "vec"},
             {"rgpr(vec,ps) -> vec", "svBtps", "svBtps", "vec", "ps"},
             {"rgpr(ps,s) -> s", "svBtps", "svBtps", "ps", "s"},
             {"rgpr(s,ps) -> s", "svBtps", "svBtps", "s", "ps"},
             // trivec
             {"rgpr(trivec,trivec) -> mv_e", "svBtps1", "svBtps2", "trivec", "trivec"},
             {"rgpr(trivec,bivec) -> mv_u", "svBtps", "svBtps", "trivec", "bivec"},
             {"rgpr(bivec,trivec) -> mv_u", "svBtps", "svBtps", "bivec", "trivec"},
             {"rgpr(trivec,vec) -> mv_e", "svBtps", "svBtps", "trivec", "vec"},
             {"rgpr(vec,trivec) -> mv_e", "svBtps", "svBtps", "vec", "trivec"},
             {"rgpr(trivec,s) -> vec", "svBtps", "svBtps", "trivec", "s"},
             {"rgpr(s,trivec) -> vec", "svBtps", "svBtps", "s", "trivec"},
             // bivec
             {"rgpr(bivec,bivec) -> mv_e", "svBtps1", "svBtps2", "bivec", "bivec"},
             {"rgpr(bivec,vec) -> mv_u", "svBtps", "svBtps", "bivec", "vec"},
             {"rgpr(vec,bivec) -> mv_u", "svBtps", "svBtps", "vec", "bivec"},
             {"rgpr(bivec,s) -> bivec", "svBtps", "svBtps", "bivec", "s"},
             {"rgpr(s,bivec) -> bivec", "svBtps", "svBtps", "s", "bivec"},
             // vec
             {"rgpr(vec,vec) -> mv_e", "svBtps1", "svBtps2", "vec", "vec"},
             {"rgpr(vec,s) -> trivec", "svBtps", "svBtps", "vec", "s"},
             {"rgpr(s,vec) -> trivec", "svBtps", "svBtps", "s", "vec"},
             // s
             {"rgpr(s,s) -> ps", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga2dc_rcmt_config()
{
    return {
        .product_name = "rcmt",
        .description = "cga2dc regressive commutator product",
        .display_name = "regressive commutator product",
        // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
        // "left_filter", "right_filter"}
        .cases =
            {// mv
             {"rcmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
             // ps
             {"rcmt(ps,ps) -> 0", "svBtps1", "svBtps2", "ps", "ps"},
             {"rcmt(ps,trivec) -> 0", "svBtps", "svBtps", "ps", "trivec"},
             {"rcmt(trivec,ps) -> 0", "svBtps", "svBtps", "trivec", "ps"},
             {"rcmt(ps,bivec) -> 0", "svBtps", "svBtps", "ps", "bivec"},
             {"rcmt(bivec,ps) -> 0", "svBtps", "svBtps", "bivec", "ps"},
             {"rcmt(ps,vec) -> 0", "svBtps", "svBtps", "ps", "vec"},
             {"rcmt(vec,ps) -> 0", "svBtps", "svBtps", "vec", "ps"},
             {"rcmt(ps,s) -> 0", "svBtps", "svBtps", "ps", "s"},
             {"rcmt(s,ps) -> 0", "svBtps", "svBtps", "s", "ps"},
             // trivec
             {"rcmt(trivec,trivec) -> bivec", "svBtps1", "svBtps2", "trivec", "trivec"},
             {"rcmt(trivec,bivec) -> trivec", "svBtps", "svBtps", "trivec", "bivec"},
             {"rcmt(bivec,trivec) -> trivec", "svBtps", "svBtps", "bivec", "trivec"},
             {"rcmt(trivec,vec) -> s", "svBtps", "svBtps", "trivec", "vec"},
             {"rcmt(vec,trivec) -> s", "svBtps", "svBtps", "vec", "trivec"},
             {"rcmt(trivec,s) -> vec", "svBtps", "svBtps", "trivec", "s"},
             {"rcmt(s,trivec) -> vec", "svBtps", "svBtps", "s", "trivec"},
             // bivec
             {"rcmt(bivec,bivec) -> bivec", "svBtps1", "svBtps2", "bivec", "bivec"},
             {"rcmt(bivec,vec) -> vec", "svBtps", "svBtps", "bivec", "vec"},
             {"rcmt(vec,bivec) -> vec", "svBtps", "svBtps", "vec", "bivec"},
             {"rcmt(bivec,s) -> 0", "svBtps", "svBtps", "bivec", "s"},
             {"rcmt(s,bivec) -> 0", "svBtps", "svBtps", "s", "bivec"},
             // vec
             {"rcmt(vec,vec) -> bivec", "svBtps1", "svBtps2", "vec", "vec"},
             {"rcmt(vec,s) -> trivec", "svBtps", "svBtps", "vec", "s"},
             {"rcmt(s,vec) -> trivec", "svBtps", "svBtps", "s", "vec"},
             // s
             {"rcmt(s,s) -> 0", "svBtps1", "svBtps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga2dc_rwdg_config()
{
    return {
        .product_name = "rwdg",
        .description = "cga2dc regressive wedge product",
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

ProductConfig get_cga2dc_rdot_config()
{
    return {.product_name = "rdot",
            .description = "cga2dc regressive inner product",
            .display_name = "regressive inner product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases =
                {// mv
                 {"rdot(mv,mv) -> ps", "A", "B", "mv", "mv"},
                 {"rdot(mv_e,mv_e) -> ps", "A_even", "B_even", "mv_e", "mv_e"},
                 {"rdot(mv_u,mv_u) -> ps", "A_odd", "B_odd", "mv_u", "mv_u"},
                 // ps
                 {"rdot(ps,ps) -> ps", "svBtps1", "svBtps2", "ps", "ps"},
                 // trivec
                 {"rdot(trivec,trivec) -> ps", "svBtps1", "svBtps2", "trivec", "trivec"},
                 // bivec
                 {"rdot(bivec,bivec) -> ps", "svBtps1", "svBtps2", "bivec", "bivec"},
                 // vec
                 {"rdot(vec,vec) -> ps", "svBtps1", "svBtps2", "vec", "vec"},
                 // s
                 {"rdot(s,s) -> ps", "svBtps1", "svBtps2", "s", "s"}},
            .is_sandwich_product = false,
            .uses_brace_switch = false,
            .show_basis_table = true};
}


ProductConfig get_cga2dc_sandwich_rgpr_config()
{
    return {
        .product_name = "sandwich_rgpr",
        .description = "cga2dc regressive sandwich product",
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
