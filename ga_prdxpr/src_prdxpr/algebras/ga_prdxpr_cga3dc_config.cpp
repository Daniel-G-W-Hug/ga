// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "algebras/ga_prdxpr_cga3dc_config.hpp"
#include "algebras/ga_prdxpr_cga3dc.hpp"

////////////////////////////////////////////////////////////////////////////////
// cga3dc algebra configuration
//
// Basis vectors: e1, e2, e3 Euclidean; e4 projective (w); e5 round (u).
// The null pair satisfies e4*e4 = e5*e5 = 0 and e4.e5 = -1, i.e. the metric
// has an off-diagonal coupling — carried by metric_matrix (metric_signature
// holds its diagonal, as required by the validation).
//
// Diagonalization (exact, dyadic coefficients): internal orthogonal basis
// (e1, e2, e3, e-, e+) with signature (+1, +1, +1, -1, +1) and
//
//     e4 = 0.5*(e- - e+)          e- =  e4 + 0.5*e5
//     e5 = e- + e+                e+ = -e4 + 0.5*e5
//
// The normalized point embedding q = x*e1 + y*e2 + z*e3 + e4
// + 0.5*(x²+y²+z²)*e5 is null: q ⟑ q = 0.
////////////////////////////////////////////////////////////////////////////////

AlgebraConfig get_cga3dc_algebra_config()
{
    // Extract basis prefix from vector basis and validate consistency
    std::string const prefix = extract_basis_prefix(mv3dc_basis_kvec[1]);
    validate_basis_consistency(mv3dc_basis, mv3dc_basis_kvec, prefix, one_str());

    AlgebraConfig config;
    config.basis_vectors = {"e1", "e2", "e3", "e4", "e5"};
    config.metric_signature = {1, 1, 1, 0, 0}; // diagonal of metric_matrix
    config.multivector_basis = mv3dc_basis;
    config.scalar_name = one_str();
    config.basis_prefix = prefix;

    config.metric_matrix = {{1, 0, 0, 0, 0},
                            {0, 1, 0, 0, 0},
                            {0, 0, 1, 0, 0},
                            {0, 0, 0, 0, -1},
                            {0, 0, 0, -1, 0}};
    config.basis_change = {{1.0, 0.0, 0.0, 0.0, 0.0},
                           {0.0, 1.0, 0.0, 0.0, 0.0},
                           {0.0, 0.0, 1.0, 0.0, 0.0},
                           {0.0, 0.0, 0.0, 0.5, -0.5},
                           {0.0, 0.0, 0.0, 1.0, 1.0}};
    config.basis_change_inv = {{1.0, 0.0, 0.0, 0.0, 0.0},
                               {0.0, 1.0, 0.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.0, 0.0},
                               {0.0, 0.0, 0.0, 1.0, 0.5},
                               {0.0, 0.0, 0.0, -1.0, 0.5}};
    config.diag_signature = {1, 1, 1, -1, 1};

    return config;
}

////////////////////////////////////////////////////////////////////////////////
// Generator-facing algebra data and product configurations
////////////////////////////////////////////////////////////////////////////////

namespace configurable {

AlgebraData create_cga3dc_algebra_data()
{
    AlgebraData cga3dc;
    cga3dc.name = "cga3dc";
    cga3dc.basis = mv3dc_basis; // Use existing basis
    cga3dc.dimension = 5;       // e1, e2, e3: Euclidean; e4: projective; e5: round

    // Map coefficient names to existing coefficient objects
    cga3dc.coefficients = {{"A", mv3dc_coeff_A},
                           {"B", mv3dc_coeff_B},
                           {"M", mv3dc_coeff_M},
                           {"M1", mv3dc_coeff_M1},
                           {"M2", mv3dc_coeff_M2},
                           {"M_even", mv3dc_coeff_M_even},
                           {"M_odd", mv3dc_coeff_M_odd},
                           {"A_even", mv3dc_coeff_A_even},
                           {"B_even", mv3dc_coeff_B_even},
                           {"A_odd", mv3dc_coeff_A_odd},
                           {"B_odd", mv3dc_coeff_B_odd},
                           {"R_even", mv3dc_coeff_R_even},
                           {"R_odd", mv3dc_coeff_R_odd},
                           {"R_rev_even", mv3dc_coeff_R_rev_even},
                           {"R_rev_odd", mv3dc_coeff_R_rev_odd},
                           {"R_rrev_even", mv3dc_coeff_R_rrev_even},
                           {"R_rrev_odd", mv3dc_coeff_R_rrev_odd},
                           {"svBtQps", mv3dc_coeff_svBtQps},
                           {"svBtQps1", mv3dc_coeff_svBtQps1},
                           {"svBtQps2", mv3dc_coeff_svBtQps2}};

    // Map filter names to existing filter enums
    cga3dc.filters_5d = {{"s", filter_5d::s},
                         {"vec", filter_5d::vec},
                         {"bivec", filter_5d::bivec},
                         {"trivec", filter_5d::trivec},
                         {"quadvec", filter_5d::quadvec},
                         {"ps", filter_5d::ps},
                         {"mv_e", filter_5d::mv_e},
                         {"mv_u", filter_5d::mv_u},
                         {"mv", filter_5d::mv}};

    return cga3dc;
}

ProductConfig get_cga3dc_gpr_config()
{
    return {
        .product_name = "gpr",
        .description = "cga3dc geometric product",
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
             {"gpr(mv,ps) -> mv", "M", "svBtQps", "mv", "ps"},
             {"gpr(ps,mv) -> mv", "svBtQps", "M", "ps", "mv"},
             {"gpr(mv,quadvec) -> mv", "M", "svBtQps", "mv", "quadvec"},
             {"gpr(quadvec,mv) -> mv", "svBtQps", "M", "quadvec", "mv"},
             {"gpr(mv,trivec) -> mv", "M", "svBtQps", "mv", "trivec"},
             {"gpr(trivec,mv) -> mv", "svBtQps", "M", "trivec", "mv"},
             {"gpr(mv,bivec) -> mv", "M", "svBtQps", "mv", "bivec"},
             {"gpr(bivec,mv) -> mv", "svBtQps", "M", "bivec", "mv"},
             {"gpr(mv,vec) -> mv", "M", "svBtQps", "mv", "vec"},
             {"gpr(vec,mv) -> mv", "svBtQps", "M", "vec", "mv"},
             {"gpr(mv,s) -> mv", "M", "svBtQps", "mv", "s"},
             {"gpr(s,mv) -> mv", "svBtQps", "M", "s", "mv"},
             // mv_e
             {"gpr(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
             {"gpr(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
             {"gpr(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
             {"gpr(mv_e,ps) -> mv_u", "M_even", "svBtQps", "mv_e", "ps"},
             {"gpr(ps,mv_e) -> mv_u", "svBtQps", "M_even", "ps", "mv_e"},
             {"gpr(mv_e,quadvec) -> mv_e", "M_even", "svBtQps", "mv_e", "quadvec"},
             {"gpr(quadvec,mv_e) -> mv_e", "svBtQps", "M_even", "quadvec", "mv_e"},
             {"gpr(mv_e,trivec) -> mv_u", "M_even", "svBtQps", "mv_e", "trivec"},
             {"gpr(trivec,mv_e) -> mv_u", "svBtQps", "M_even", "trivec", "mv_e"},
             {"gpr(mv_e,bivec) -> mv_e", "M_even", "svBtQps", "mv_e", "bivec"},
             {"gpr(bivec,mv_e) -> mv_e", "svBtQps", "M_even", "bivec", "mv_e"},
             {"gpr(mv_e,vec) -> mv_u", "M_even", "svBtQps", "mv_e", "vec"},
             {"gpr(vec,mv_e) -> mv_u", "svBtQps", "M_even", "vec", "mv_e"},
             {"gpr(mv_e,s) -> mv_e", "M_even", "svBtQps", "mv_e", "s"},
             {"gpr(s,mv_e) -> mv_e", "svBtQps", "M_even", "s", "mv_e"},
             // mv_u
             {"gpr(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"gpr(mv_u,ps) -> mv_e", "M_odd", "svBtQps", "mv_u", "ps"},
             {"gpr(ps,mv_u) -> mv_e", "svBtQps", "M_odd", "ps", "mv_u"},
             {"gpr(mv_u,quadvec) -> mv_u", "M_odd", "svBtQps", "mv_u", "quadvec"},
             {"gpr(quadvec,mv_u) -> mv_u", "svBtQps", "M_odd", "quadvec", "mv_u"},
             {"gpr(mv_u,trivec) -> mv_e", "M_odd", "svBtQps", "mv_u", "trivec"},
             {"gpr(trivec,mv_u) -> mv_e", "svBtQps", "M_odd", "trivec", "mv_u"},
             {"gpr(mv_u,bivec) -> mv_u", "M_odd", "svBtQps", "mv_u", "bivec"},
             {"gpr(bivec,mv_u) -> mv_u", "svBtQps", "M_odd", "bivec", "mv_u"},
             {"gpr(mv_u,vec) -> mv_e", "M_odd", "svBtQps", "mv_u", "vec"},
             {"gpr(vec,mv_u) -> mv_e", "svBtQps", "M_odd", "vec", "mv_u"},
             {"gpr(mv_u,s) -> mv_u", "M_odd", "svBtQps", "mv_u", "s"},
             {"gpr(s,mv_u) -> mv_u", "svBtQps", "M_odd", "s", "mv_u"},
             // ps
             // quadvec
             // trivec
             // bivec
             // vec
             // s
             // ps
             {"gpr(ps,ps) -> s", "svBtQps1", "svBtQps2", "ps", "ps"},
             {"gpr(ps,quadvec) -> vec", "svBtQps", "svBtQps", "ps", "quadvec"},
             {"gpr(quadvec,ps) -> vec", "svBtQps", "svBtQps", "quadvec", "ps"},
             {"gpr(ps,trivec) -> bivec", "svBtQps", "svBtQps", "ps", "trivec"},
             {"gpr(trivec,ps) -> bivec", "svBtQps", "svBtQps", "trivec", "ps"},
             {"gpr(ps,bivec) -> trivec", "svBtQps", "svBtQps", "ps", "bivec"},
             {"gpr(bivec,ps) -> trivec", "svBtQps", "svBtQps", "bivec", "ps"},
             {"gpr(ps,vec) -> quadvec", "svBtQps", "svBtQps", "ps", "vec"},
             {"gpr(vec,ps) -> quadvec", "svBtQps", "svBtQps", "vec", "ps"},
             {"gpr(ps,s) -> ps", "svBtQps", "svBtQps", "ps", "s"},
             {"gpr(s,ps) -> ps", "svBtQps", "svBtQps", "s", "ps"},
             // quadvec
             {"gpr(quadvec,quadvec) -> mv_e", "svBtQps1", "svBtQps2", "quadvec",
              "quadvec"},
             {"gpr(quadvec,trivec) -> mv_u", "svBtQps", "svBtQps", "quadvec", "trivec"},
             {"gpr(trivec,quadvec) -> mv_u", "svBtQps", "svBtQps", "trivec", "quadvec"},
             {"gpr(quadvec,bivec) -> mv_e", "svBtQps", "svBtQps", "quadvec", "bivec"},
             {"gpr(bivec,quadvec) -> mv_e", "svBtQps", "svBtQps", "bivec", "quadvec"},
             {"gpr(quadvec,vec) -> mv_u", "svBtQps", "svBtQps", "quadvec", "vec"},
             {"gpr(vec,quadvec) -> mv_u", "svBtQps", "svBtQps", "vec", "quadvec"},
             {"gpr(quadvec,s) -> quadvec", "svBtQps", "svBtQps", "quadvec", "s"},
             {"gpr(s,quadvec) -> quadvec", "svBtQps", "svBtQps", "s", "quadvec"},
             // trivec
             {"gpr(trivec,trivec) -> mv_e", "svBtQps1", "svBtQps2", "trivec", "trivec"},
             {"gpr(trivec,bivec) -> mv_u", "svBtQps", "svBtQps", "trivec", "bivec"},
             {"gpr(bivec,trivec) -> mv_u", "svBtQps", "svBtQps", "bivec", "trivec"},
             {"gpr(trivec,vec) -> mv_e", "svBtQps", "svBtQps", "trivec", "vec"},
             {"gpr(vec,trivec) -> mv_e", "svBtQps", "svBtQps", "vec", "trivec"},
             {"gpr(trivec,s) -> trivec", "svBtQps", "svBtQps", "trivec", "s"},
             {"gpr(s,trivec) -> trivec", "svBtQps", "svBtQps", "s", "trivec"},
             // bivec
             {"gpr(bivec,bivec) -> mv_e", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             {"gpr(bivec,vec) -> mv_u", "svBtQps", "svBtQps", "bivec", "vec"},
             {"gpr(vec,bivec) -> mv_u", "svBtQps", "svBtQps", "vec", "bivec"},
             {"gpr(bivec,s) -> bivec", "svBtQps", "svBtQps", "bivec", "s"},
             {"gpr(s,bivec) -> bivec", "svBtQps", "svBtQps", "s", "bivec"},
             // vec
             {"gpr(vec,vec) -> mv_e", "svBtQps1", "svBtQps2", "vec", "vec"},
             {"gpr(vec,s) -> vec", "svBtQps", "svBtQps", "vec", "s"},
             {"gpr(s,vec) -> vec", "svBtQps", "svBtQps", "s", "vec"},
             // s
             {"gpr(s,s) -> s", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_wdg_config()
{
    return {
        .product_name = "wdg",
        .description = "cga3dc wedge product",
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
             {"wdg(mv,ps) -> mv", "M", "svBtQps", "mv", "ps"},
             {"wdg(ps,mv) -> mv", "svBtQps", "M", "ps", "mv"},
             {"wdg(mv,quadvec) -> mv", "M", "svBtQps", "mv", "quadvec"},
             {"wdg(quadvec,mv) -> mv", "svBtQps", "M", "quadvec", "mv"},
             {"wdg(mv,trivec) -> mv", "M", "svBtQps", "mv", "trivec"},
             {"wdg(trivec,mv) -> mv", "svBtQps", "M", "trivec", "mv"},
             {"wdg(mv,bivec) -> mv", "M", "svBtQps", "mv", "bivec"},
             {"wdg(bivec,mv) -> mv", "svBtQps", "M", "bivec", "mv"},
             {"wdg(mv,vec) -> mv", "M", "svBtQps", "mv", "vec"},
             {"wdg(vec,mv) -> mv", "svBtQps", "M", "vec", "mv"},
             {"wdg(mv,s) -> mv", "M", "svBtQps", "mv", "s"},
             {"wdg(s,mv) -> mv", "svBtQps", "M", "s", "mv"},
             // mv_e
             {"wdg(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
             {"wdg(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
             {"wdg(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
             {"wdg(mv_e,ps) -> mv_u", "M_even", "svBtQps", "mv_e", "ps"},
             {"wdg(ps,mv_e) -> mv_u", "svBtQps", "M_even", "ps", "mv_e"},
             {"wdg(mv_e,quadvec) -> mv_e", "M_even", "svBtQps", "mv_e", "quadvec"},
             {"wdg(quadvec,mv_e) -> mv_e", "svBtQps", "M_even", "quadvec", "mv_e"},
             {"wdg(mv_e,trivec) -> mv_u", "M_even", "svBtQps", "mv_e", "trivec"},
             {"wdg(trivec,mv_e) -> mv_u", "svBtQps", "M_even", "trivec", "mv_e"},
             {"wdg(mv_e,bivec) -> mv_e", "M_even", "svBtQps", "mv_e", "bivec"},
             {"wdg(bivec,mv_e) -> mv_e", "svBtQps", "M_even", "bivec", "mv_e"},
             {"wdg(mv_e,vec) -> mv_u", "M_even", "svBtQps", "mv_e", "vec"},
             {"wdg(vec,mv_e) -> mv_u", "svBtQps", "M_even", "vec", "mv_e"},
             {"wdg(mv_e,s) -> mv_e", "M_even", "svBtQps", "mv_e", "s"},
             {"wdg(s,mv_e) -> mv_e", "svBtQps", "M_even", "s", "mv_e"},
             // mv_u
             {"wdg(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"wdg(mv_u,ps) -> 0", "M_odd", "svBtQps", "mv_u", "ps"},
             {"wdg(ps,mv_u) -> 0", "svBtQps", "M_odd", "ps", "mv_u"},
             {"wdg(mv_u,quadvec) -> mv_u", "M_odd", "svBtQps", "mv_u", "quadvec"},
             {"wdg(quadvec,mv_u) -> mv_u", "svBtQps", "M_odd", "quadvec", "mv_u"},
             {"wdg(mv_u,trivec) -> mv_e", "M_odd", "svBtQps", "mv_u", "trivec"},
             {"wdg(trivec,mv_u) -> mv_e", "svBtQps", "M_odd", "trivec", "mv_u"},
             {"wdg(mv_u,bivec) -> mv_u", "M_odd", "svBtQps", "mv_u", "bivec"},
             {"wdg(bivec,mv_u) -> mv_u", "svBtQps", "M_odd", "bivec", "mv_u"},
             {"wdg(mv_u,vec) -> mv_e", "M_odd", "svBtQps", "mv_u", "vec"},
             {"wdg(vec,mv_u) -> mv_e", "svBtQps", "M_odd", "vec", "mv_u"},
             {"wdg(mv_u,s) -> mv_u", "M_odd", "svBtQps", "mv_u", "s"},
             {"wdg(s,mv_u) -> mv_u", "svBtQps", "M_odd", "s", "mv_u"},
             // ps
             // quadvec
             // trivec
             // bivec
             // vec
             // s
             // ps
             {"wdg(ps,ps) -> 0", "svBtQps1", "svBtQps2", "ps", "ps"},
             {"wdg(ps,quadvec) -> 0", "svBtQps", "svBtQps", "ps", "quadvec"},
             {"wdg(quadvec,ps) -> 0", "svBtQps", "svBtQps", "quadvec", "ps"},
             {"wdg(ps,trivec) -> 0", "svBtQps", "svBtQps", "ps", "trivec"},
             {"wdg(trivec,ps) -> 0", "svBtQps", "svBtQps", "trivec", "ps"},
             {"wdg(ps,bivec) -> 0", "svBtQps", "svBtQps", "ps", "bivec"},
             {"wdg(bivec,ps) -> 0", "svBtQps", "svBtQps", "bivec", "ps"},
             {"wdg(ps,vec) -> 0", "svBtQps", "svBtQps", "ps", "vec"},
             {"wdg(vec,ps) -> 0", "svBtQps", "svBtQps", "vec", "ps"},
             {"wdg(ps,s) -> ps", "svBtQps", "svBtQps", "ps", "s"},
             {"wdg(s,ps) -> ps", "svBtQps", "svBtQps", "s", "ps"},
             // quadvec
             {"wdg(quadvec,quadvec) -> 0", "svBtQps1", "svBtQps2", "quadvec", "quadvec"},
             {"wdg(quadvec,trivec) -> 0", "svBtQps", "svBtQps", "quadvec", "trivec"},
             {"wdg(trivec,quadvec) -> 0", "svBtQps", "svBtQps", "trivec", "quadvec"},
             {"wdg(quadvec,bivec) -> 0", "svBtQps", "svBtQps", "quadvec", "bivec"},
             {"wdg(bivec,quadvec) -> 0", "svBtQps", "svBtQps", "bivec", "quadvec"},
             {"wdg(quadvec,vec) -> ps", "svBtQps", "svBtQps", "quadvec", "vec"},
             {"wdg(vec,quadvec) -> ps", "svBtQps", "svBtQps", "vec", "quadvec"},
             {"wdg(quadvec,s) -> quadvec", "svBtQps", "svBtQps", "quadvec", "s"},
             {"wdg(s,quadvec) -> quadvec", "svBtQps", "svBtQps", "s", "quadvec"},
             // trivec
             {"wdg(trivec,trivec) -> 0", "svBtQps1", "svBtQps2", "trivec", "trivec"},
             {"wdg(trivec,bivec) -> ps", "svBtQps", "svBtQps", "trivec", "bivec"},
             {"wdg(bivec,trivec) -> ps", "svBtQps", "svBtQps", "bivec", "trivec"},
             {"wdg(trivec,vec) -> quadvec", "svBtQps", "svBtQps", "trivec", "vec"},
             {"wdg(vec,trivec) -> quadvec", "svBtQps", "svBtQps", "vec", "trivec"},
             {"wdg(trivec,s) -> trivec", "svBtQps", "svBtQps", "trivec", "s"},
             {"wdg(s,trivec) -> trivec", "svBtQps", "svBtQps", "s", "trivec"},
             // bivec
             {"wdg(bivec,bivec) -> quadvec", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             {"wdg(bivec,vec) -> trivec", "svBtQps", "svBtQps", "bivec", "vec"},
             {"wdg(vec,bivec) -> trivec", "svBtQps", "svBtQps", "vec", "bivec"},
             {"wdg(bivec,s) -> bivec", "svBtQps", "svBtQps", "bivec", "s"},
             {"wdg(s,bivec) -> bivec", "svBtQps", "svBtQps", "s", "bivec"},
             // vec
             {"wdg(vec,vec) -> bivec", "svBtQps1", "svBtQps2", "vec", "vec"},
             {"wdg(vec,s) -> vec", "svBtQps", "svBtQps", "vec", "s"},
             {"wdg(s,vec) -> vec", "svBtQps", "svBtQps", "s", "vec"},
             // s
             {"wdg(s,s) -> s", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_dot_config()
{
    return {
        .product_name = "dot",
        .description = "cga3dc inner product",
        .display_name = "inner product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// mv
             {"dot(mv,mv) -> s", "A", "B", "mv", "mv"},
             {"dot(mv_e,mv_e) -> s", "A_even", "B_even", "mv_e", "mv_e"},
             {"dot(mv_u,mv_u) -> s", "A_odd", "B_odd", "mv_u", "mv_u"},
             // ps
             {"dot(ps,ps) -> s", "svBtQps1", "svBtQps2", "ps", "ps"},
             // quadvec
             {"dot(quadvec,quadvec) -> s", "svBtQps1", "svBtQps2", "quadvec", "quadvec"},
             // trivec
             {"dot(trivec,trivec) -> s", "svBtQps1", "svBtQps2", "trivec", "trivec"},
             // bivec
             {"dot(bivec,bivec) -> s", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             // vec
             {"dot(vec,vec) -> s", "svBtQps1", "svBtQps2", "vec", "vec"},
             // s
             {"dot(s,s) -> s", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_cmt_config()
{
    return {
        .product_name = "cmt",
        .description = "cga3dc commutator product",
        .display_name = "commutator product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// mv
             {"cmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
             // ps
             {"cmt(ps,ps) -> 0", "svBtQps1", "svBtQps2", "ps", "ps"},
             {"cmt(ps,quadvec) -> 0", "svBtQps", "svBtQps", "ps", "quadvec"},
             {"cmt(quadvec,ps) -> 0", "svBtQps", "svBtQps", "quadvec", "ps"},
             {"cmt(ps,trivec) -> 0", "svBtQps", "svBtQps", "ps", "trivec"},
             {"cmt(trivec,ps) -> 0", "svBtQps", "svBtQps", "trivec", "ps"},
             {"cmt(ps,bivec) -> 0", "svBtQps", "svBtQps", "ps", "bivec"},
             {"cmt(bivec,ps) -> 0", "svBtQps", "svBtQps", "bivec", "ps"},
             {"cmt(ps,vec) -> 0", "svBtQps", "svBtQps", "ps", "vec"},
             {"cmt(vec,ps) -> 0", "svBtQps", "svBtQps", "vec", "ps"},
             {"cmt(ps,s) -> 0", "svBtQps", "svBtQps", "ps", "s"},
             {"cmt(s,ps) -> 0", "svBtQps", "svBtQps", "s", "ps"},
             // quadvec
             {"cmt(quadvec,quadvec) -> bivec", "svBtQps1", "svBtQps2", "quadvec",
              "quadvec"},
             {"cmt(quadvec,trivec) -> vec", "svBtQps", "svBtQps", "quadvec", "trivec"},
             {"cmt(trivec,quadvec) -> vec", "svBtQps", "svBtQps", "trivec", "quadvec"},
             {"cmt(quadvec,bivec) -> quadvec", "svBtQps", "svBtQps", "quadvec", "bivec"},
             {"cmt(bivec,quadvec) -> quadvec", "svBtQps", "svBtQps", "bivec", "quadvec"},
             {"cmt(quadvec,vec) -> trivec", "svBtQps", "svBtQps", "quadvec", "vec"},
             {"cmt(vec,quadvec) -> trivec", "svBtQps", "svBtQps", "vec", "quadvec"},
             {"cmt(quadvec,s) -> 0", "svBtQps", "svBtQps", "quadvec", "s"},
             {"cmt(s,quadvec) -> 0", "svBtQps", "svBtQps", "s", "quadvec"},
             // trivec
             {"cmt(trivec,trivec) -> bivec", "svBtQps1", "svBtQps2", "trivec", "trivec"},
             {"cmt(trivec,bivec) -> trivec", "svBtQps", "svBtQps", "trivec", "bivec"},
             {"cmt(bivec,trivec) -> trivec", "svBtQps", "svBtQps", "bivec", "trivec"},
             {"cmt(trivec,vec) -> quadvec", "svBtQps", "svBtQps", "trivec", "vec"},
             {"cmt(vec,trivec) -> quadvec", "svBtQps", "svBtQps", "vec", "trivec"},
             {"cmt(trivec,s) -> 0", "svBtQps", "svBtQps", "trivec", "s"},
             {"cmt(s,trivec) -> 0", "svBtQps", "svBtQps", "s", "trivec"},
             // bivec
             {"cmt(bivec,bivec) -> bivec", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             {"cmt(bivec,vec) -> vec", "svBtQps", "svBtQps", "bivec", "vec"},
             {"cmt(vec,bivec) -> vec", "svBtQps", "svBtQps", "vec", "bivec"},
             {"cmt(bivec,s) -> 0", "svBtQps", "svBtQps", "bivec", "s"},
             {"cmt(s,bivec) -> 0", "svBtQps", "svBtQps", "s", "bivec"},
             // vec
             {"cmt(vec,vec) -> bivec", "svBtQps1", "svBtQps2", "vec", "vec"},
             {"cmt(vec,s) -> 0", "svBtQps", "svBtQps", "vec", "s"},
             {"cmt(s,vec) -> 0", "svBtQps", "svBtQps", "s", "vec"},
             // s
             {"cmt(s,s) -> 0", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_l_contract_config()
{
    return {
        .product_name = "l_contract",
        .description = "cga3dc left contraction",
        .display_name = "left contraction",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// mv
             {"l_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"l_contract(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
             {"l_contract(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
             {"l_contract(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
             {"l_contract(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
             {"l_contract(mv,ps) -> mv", "M", "svBtQps", "mv", "ps"},
             {"l_contract(ps,mv) -> mv", "svBtQps", "M", "ps", "mv"},
             {"l_contract(mv,quadvec) -> mv", "M", "svBtQps", "mv", "quadvec"},
             {"l_contract(quadvec,mv) -> mv", "svBtQps", "M", "quadvec", "mv"},
             {"l_contract(mv,trivec) -> mv", "M", "svBtQps", "mv", "trivec"},
             {"l_contract(trivec,mv) -> mv", "svBtQps", "M", "trivec", "mv"},
             {"l_contract(mv,bivec) -> mv", "M", "svBtQps", "mv", "bivec"},
             {"l_contract(bivec,mv) -> mv", "svBtQps", "M", "bivec", "mv"},
             {"l_contract(mv,vec) -> mv", "M", "svBtQps", "mv", "vec"},
             {"l_contract(vec,mv) -> mv", "svBtQps", "M", "vec", "mv"},
             {"l_contract(mv,s) -> mv", "M", "svBtQps", "mv", "s"},
             {"l_contract(s,mv) -> mv", "svBtQps", "M", "s", "mv"},
             // mv_e
             {"l_contract(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
             {"l_contract(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
             {"l_contract(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
             {"l_contract(mv_e,ps) -> mv_u", "M_even", "svBtQps", "mv_e", "ps"},
             {"l_contract(ps,mv_e) -> 0", "svBtQps", "M_even", "ps", "mv_e"},
             {"l_contract(mv_e,quadvec) -> mv_e", "M_even", "svBtQps", "mv_e", "quadvec"},
             {"l_contract(quadvec,mv_e) -> mv_e", "svBtQps", "M_even", "quadvec", "mv_e"},
             {"l_contract(mv_e,trivec) -> mv_u", "M_even", "svBtQps", "mv_e", "trivec"},
             {"l_contract(trivec,mv_e) -> mv_u", "svBtQps", "M_even", "trivec", "mv_e"},
             {"l_contract(mv_e,bivec) -> mv_e", "M_even", "svBtQps", "mv_e", "bivec"},
             {"l_contract(bivec,mv_e) -> mv_e", "svBtQps", "M_even", "bivec", "mv_e"},
             {"l_contract(mv_e,vec) -> mv_u", "M_even", "svBtQps", "mv_e", "vec"},
             {"l_contract(vec,mv_e) -> mv_u", "svBtQps", "M_even", "vec", "mv_e"},
             {"l_contract(mv_e,s) -> mv_e", "M_even", "svBtQps", "mv_e", "s"},
             {"l_contract(s,mv_e) -> mv_e", "svBtQps", "M_even", "s", "mv_e"},
             // mv_u
             {"l_contract(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"l_contract(mv_u,ps) -> mv_e", "M_odd", "svBtQps", "mv_u", "ps"},
             {"l_contract(ps,mv_u) -> mv_e", "svBtQps", "M_odd", "ps", "mv_u"},
             {"l_contract(mv_u,quadvec) -> mv_u", "M_odd", "svBtQps", "mv_u", "quadvec"},
             {"l_contract(quadvec,mv_u) -> mv_u", "svBtQps", "M_odd", "quadvec", "mv_u"},
             {"l_contract(mv_u,trivec) -> mv_e", "M_odd", "svBtQps", "mv_u", "trivec"},
             {"l_contract(trivec,mv_u) -> mv_e", "svBtQps", "M_odd", "trivec", "mv_u"},
             {"l_contract(mv_u,bivec) -> mv_u", "M_odd", "svBtQps", "mv_u", "bivec"},
             {"l_contract(bivec,mv_u) -> mv_u", "svBtQps", "M_odd", "bivec", "mv_u"},
             {"l_contract(mv_u,vec) -> mv_e", "M_odd", "svBtQps", "mv_u", "vec"},
             {"l_contract(vec,mv_u) -> mv_e", "svBtQps", "M_odd", "vec", "mv_u"},
             {"l_contract(mv_u,s) -> 0", "M_odd", "svBtQps", "mv_u", "s"},
             {"l_contract(s,mv_u) -> mv_u", "svBtQps", "M_odd", "s", "mv_u"},
             // ps
             {"l_contract(ps,ps) -> s", "svBtQps1", "svBtQps2", "ps", "ps"},
             {"l_contract(ps,quadvec) -> 0", "svBtQps", "svBtQps", "ps", "quadvec"},
             {"l_contract(quadvec,ps) -> vec", "svBtQps", "svBtQps", "quadvec", "ps"},
             {"l_contract(ps,trivec) -> 0", "svBtQps", "svBtQps", "ps", "trivec"},
             {"l_contract(trivec,ps) -> bivec", "svBtQps", "svBtQps", "trivec", "ps"},
             {"l_contract(ps,bivec) -> 0", "svBtQps", "svBtQps", "ps", "bivec"},
             {"l_contract(bivec,ps) -> trivec", "svBtQps", "svBtQps", "bivec", "ps"},
             {"l_contract(ps,vec) -> 0", "svBtQps", "svBtQps", "ps", "vec"},
             {"l_contract(vec,ps) -> quadvec", "svBtQps", "svBtQps", "vec", "ps"},
             {"l_contract(ps,s) -> 0", "svBtQps", "svBtQps", "ps", "s"},
             {"l_contract(s,ps) -> ps", "svBtQps", "svBtQps", "s", "ps"},
             // quadvec
             {"l_contract(quadvec,quadvec) -> s", "svBtQps1", "svBtQps2", "quadvec",
              "quadvec"},
             {"l_contract(quadvec,trivec) -> 0", "svBtQps", "svBtQps", "quadvec",
              "trivec"},
             {"l_contract(trivec,quadvec) -> vec", "svBtQps", "svBtQps", "trivec",
              "quadvec"},
             {"l_contract(quadvec,bivec) -> 0", "svBtQps", "svBtQps", "quadvec", "bivec"},
             {"l_contract(bivec,quadvec) -> bivec", "svBtQps", "svBtQps", "bivec",
              "quadvec"},
             {"l_contract(quadvec,vec) -> 0", "svBtQps", "svBtQps", "quadvec", "vec"},
             {"l_contract(vec,quadvec) -> trivec", "svBtQps", "svBtQps", "vec",
              "quadvec"},
             {"l_contract(quadvec,s) -> 0", "svBtQps", "svBtQps", "quadvec", "s"},
             {"l_contract(s,quadvec) -> quadvec", "svBtQps", "svBtQps", "s", "quadvec"},
             // trivec
             {"l_contract(trivec,trivec) -> s", "svBtQps1", "svBtQps2", "trivec",
              "trivec"},
             {"l_contract(trivec,bivec) -> 0", "svBtQps", "svBtQps", "trivec", "bivec"},
             {"l_contract(bivec,trivec) -> vec", "svBtQps", "svBtQps", "bivec", "trivec"},
             {"l_contract(trivec,vec) -> 0", "svBtQps", "svBtQps", "trivec", "vec"},
             {"l_contract(vec,trivec) -> bivec", "svBtQps", "svBtQps", "vec", "trivec"},
             {"l_contract(trivec,s) -> 0", "svBtQps", "svBtQps", "trivec", "s"},
             {"l_contract(s,trivec) -> trivec", "svBtQps", "svBtQps", "s", "trivec"},
             // bivec
             {"l_contract(bivec,bivec) -> s", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             {"l_contract(bivec,vec) -> 0", "svBtQps", "svBtQps", "bivec", "vec"},
             {"l_contract(vec,bivec) -> vec", "svBtQps", "svBtQps", "vec", "bivec"},
             {"l_contract(bivec,s) -> 0", "svBtQps", "svBtQps", "bivec", "s"},
             {"l_contract(s,bivec) -> bivec", "svBtQps", "svBtQps", "s", "bivec"},
             // vec
             {"l_contract(vec,vec) -> s", "svBtQps1", "svBtQps2", "vec", "vec"},
             {"l_contract(vec,s) -> 0", "svBtQps", "svBtQps", "vec", "s"},
             {"l_contract(s,vec) -> vec", "svBtQps", "svBtQps", "s", "vec"},
             // s
             {"l_contract(s,s) -> s", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_r_contract_config()
{
    return {
        .product_name = "r_contract",
        .description = "cga3dc right contraction",
        .display_name = "right contraction",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// mv
             {"r_contract(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"r_contract(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
             {"r_contract(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
             {"r_contract(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
             {"r_contract(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
             {"r_contract(mv,ps) -> mv", "M", "svBtQps", "mv", "ps"},
             {"r_contract(ps,mv) -> mv", "svBtQps", "M", "ps", "mv"},
             {"r_contract(mv,quadvec) -> mv", "M", "svBtQps", "mv", "quadvec"},
             {"r_contract(quadvec,mv) -> mv", "svBtQps", "M", "quadvec", "mv"},
             {"r_contract(mv,trivec) -> mv", "M", "svBtQps", "mv", "trivec"},
             {"r_contract(trivec,mv) -> mv", "svBtQps", "M", "trivec", "mv"},
             {"r_contract(mv,bivec) -> mv", "M", "svBtQps", "mv", "bivec"},
             {"r_contract(bivec,mv) -> mv", "svBtQps", "M", "bivec", "mv"},
             {"r_contract(mv,vec) -> mv", "M", "svBtQps", "mv", "vec"},
             {"r_contract(vec,mv) -> mv", "svBtQps", "M", "vec", "mv"},
             {"r_contract(mv,s) -> mv", "M", "svBtQps", "mv", "s"},
             {"r_contract(s,mv) -> mv", "svBtQps", "M", "s", "mv"},
             // mv_e
             {"r_contract(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
             {"r_contract(mv_e,mv_u) -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
             {"r_contract(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
             {"r_contract(mv_e,ps) -> 0", "M_even", "svBtQps", "mv_e", "ps"},
             {"r_contract(ps,mv_e) -> mv_u", "svBtQps", "M_even", "ps", "mv_e"},
             {"r_contract(mv_e,quadvec) -> mv_e", "M_even", "svBtQps", "mv_e", "quadvec"},
             {"r_contract(quadvec,mv_e) -> mv_e", "svBtQps", "M_even", "quadvec", "mv_e"},
             {"r_contract(mv_e,trivec) -> mv_u", "M_even", "svBtQps", "mv_e", "trivec"},
             {"r_contract(trivec,mv_e) -> mv_u", "svBtQps", "M_even", "trivec", "mv_e"},
             {"r_contract(mv_e,bivec) -> mv_e", "M_even", "svBtQps", "mv_e", "bivec"},
             {"r_contract(bivec,mv_e) -> mv_e", "svBtQps", "M_even", "bivec", "mv_e"},
             {"r_contract(mv_e,vec) -> mv_u", "M_even", "svBtQps", "mv_e", "vec"},
             {"r_contract(vec,mv_e) -> mv_u", "svBtQps", "M_even", "vec", "mv_e"},
             {"r_contract(mv_e,s) -> mv_e", "M_even", "svBtQps", "mv_e", "s"},
             {"r_contract(s,mv_e) -> mv_e", "svBtQps", "M_even", "s", "mv_e"},
             // mv_u
             {"r_contract(mv_u,mv_u) -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"r_contract(mv_u,ps) -> mv_e", "M_odd", "svBtQps", "mv_u", "ps"},
             {"r_contract(ps,mv_u) -> mv_e", "svBtQps", "M_odd", "ps", "mv_u"},
             {"r_contract(mv_u,quadvec) -> mv_u", "M_odd", "svBtQps", "mv_u", "quadvec"},
             {"r_contract(quadvec,mv_u) -> mv_u", "svBtQps", "M_odd", "quadvec", "mv_u"},
             {"r_contract(mv_u,trivec) -> mv_e", "M_odd", "svBtQps", "mv_u", "trivec"},
             {"r_contract(trivec,mv_u) -> mv_e", "svBtQps", "M_odd", "trivec", "mv_u"},
             {"r_contract(mv_u,bivec) -> mv_u", "M_odd", "svBtQps", "mv_u", "bivec"},
             {"r_contract(bivec,mv_u) -> mv_u", "svBtQps", "M_odd", "bivec", "mv_u"},
             {"r_contract(mv_u,vec) -> mv_e", "M_odd", "svBtQps", "mv_u", "vec"},
             {"r_contract(vec,mv_u) -> mv_e", "svBtQps", "M_odd", "vec", "mv_u"},
             {"r_contract(mv_u,s) -> mv_u", "M_odd", "svBtQps", "mv_u", "s"},
             {"r_contract(s,mv_u) -> 0", "svBtQps", "M_odd", "s", "mv_u"},
             // ps
             {"r_contract(ps,ps) -> s", "svBtQps1", "svBtQps2", "ps", "ps"},
             {"r_contract(ps,quadvec) -> vec", "svBtQps", "svBtQps", "ps", "quadvec"},
             {"r_contract(quadvec,ps) -> 0", "svBtQps", "svBtQps", "quadvec", "ps"},
             {"r_contract(ps,trivec) -> bivec", "svBtQps", "svBtQps", "ps", "trivec"},
             {"r_contract(trivec,ps) -> 0", "svBtQps", "svBtQps", "trivec", "ps"},
             {"r_contract(ps,bivec) -> trivec", "svBtQps", "svBtQps", "ps", "bivec"},
             {"r_contract(bivec,ps) -> 0", "svBtQps", "svBtQps", "bivec", "ps"},
             {"r_contract(ps,vec) -> quadvec", "svBtQps", "svBtQps", "ps", "vec"},
             {"r_contract(vec,ps) -> 0", "svBtQps", "svBtQps", "vec", "ps"},
             {"r_contract(ps,s) -> ps", "svBtQps", "svBtQps", "ps", "s"},
             {"r_contract(s,ps) -> 0", "svBtQps", "svBtQps", "s", "ps"},
             // quadvec
             {"r_contract(quadvec,quadvec) -> s", "svBtQps1", "svBtQps2", "quadvec",
              "quadvec"},
             {"r_contract(quadvec,trivec) -> vec", "svBtQps", "svBtQps", "quadvec",
              "trivec"},
             {"r_contract(trivec,quadvec) -> 0", "svBtQps", "svBtQps", "trivec",
              "quadvec"},
             {"r_contract(quadvec,bivec) -> bivec", "svBtQps", "svBtQps", "quadvec",
              "bivec"},
             {"r_contract(bivec,quadvec) -> 0", "svBtQps", "svBtQps", "bivec", "quadvec"},
             {"r_contract(quadvec,vec) -> trivec", "svBtQps", "svBtQps", "quadvec",
              "vec"},
             {"r_contract(vec,quadvec) -> 0", "svBtQps", "svBtQps", "vec", "quadvec"},
             {"r_contract(quadvec,s) -> quadvec", "svBtQps", "svBtQps", "quadvec", "s"},
             {"r_contract(s,quadvec) -> 0", "svBtQps", "svBtQps", "s", "quadvec"},
             // trivec
             {"r_contract(trivec,trivec) -> s", "svBtQps1", "svBtQps2", "trivec",
              "trivec"},
             {"r_contract(trivec,bivec) -> vec", "svBtQps", "svBtQps", "trivec", "bivec"},
             {"r_contract(bivec,trivec) -> 0", "svBtQps", "svBtQps", "bivec", "trivec"},
             {"r_contract(trivec,vec) -> bivec", "svBtQps", "svBtQps", "trivec", "vec"},
             {"r_contract(vec,trivec) -> 0", "svBtQps", "svBtQps", "vec", "trivec"},
             {"r_contract(trivec,s) -> trivec", "svBtQps", "svBtQps", "trivec", "s"},
             {"r_contract(s,trivec) -> 0", "svBtQps", "svBtQps", "s", "trivec"},
             // bivec
             {"r_contract(bivec,bivec) -> s", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             {"r_contract(bivec,vec) -> vec", "svBtQps", "svBtQps", "bivec", "vec"},
             {"r_contract(vec,bivec) -> 0", "svBtQps", "svBtQps", "vec", "bivec"},
             {"r_contract(bivec,s) -> bivec", "svBtQps", "svBtQps", "bivec", "s"},
             {"r_contract(s,bivec) -> 0", "svBtQps", "svBtQps", "s", "bivec"},
             // vec
             {"r_contract(vec,vec) -> s", "svBtQps1", "svBtQps2", "vec", "vec"},
             {"r_contract(vec,s) -> vec", "svBtQps", "svBtQps", "vec", "s"},
             {"r_contract(s,vec) -> 0", "svBtQps", "svBtQps", "s", "vec"},
             // s
             {"r_contract(s,s) -> s", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_l_expand_config()
{
    return {
        .product_name = "l_expand",
        .description = "cga3dc left expansion",
        .display_name = "left expansion",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// ps
             {"l_expand(ps,quadvec) -> quadvec", "svBtQps", "svBtQps", "ps", "quadvec"},
             {"l_expand(ps,trivec) -> trivec", "svBtQps", "svBtQps", "ps", "trivec"},
             {"l_expand(ps,bivec) -> bivec", "svBtQps", "svBtQps", "ps", "bivec"},
             {"l_expand(ps,vec) -> vec", "svBtQps", "svBtQps", "ps", "vec"},
             // quadvec
             {"l_expand(quadvec,trivec) -> quadvec", "svBtQps", "svBtQps", "quadvec",
              "trivec"},
             {"l_expand(quadvec,bivec) -> trivec", "svBtQps", "svBtQps", "quadvec",
              "bivec"},
             {"l_expand(quadvec,vec) -> bivec", "svBtQps", "svBtQps", "quadvec", "vec"},
             {"l_expand(quadvec,s) -> vec", "svBtQps", "svBtQps", "quadvec", "s"},
             // trivec
             {"l_expand(trivec,bivec) -> quadvec", "svBtQps", "svBtQps", "trivec",
              "bivec"},
             {"l_expand(trivec,vec) -> trivec", "svBtQps", "svBtQps", "trivec", "vec"},
             {"l_expand(trivec,s) -> bivec", "svBtQps", "svBtQps", "trivec", "s"},
             // bivec
             {"l_expand(bivec,vec) -> quadvec", "svBtQps", "svBtQps", "bivec", "vec"},
             {"l_expand(bivec,s) -> trivec", "svBtQps", "svBtQps", "bivec", "s"},
             // vec
             {"l_expand(vec,s) -> quadvec", "svBtQps", "svBtQps", "vec", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_r_expand_config()
{
    return {
        .product_name = "r_expand",
        .description = "cga3dc right expansion",
        .display_name = "right expansion",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// quadvec
             {"r_expand(quadvec,ps) -> quadvec", "svBtQps", "svBtQps", "quadvec", "ps"},
             // trivec
             {"r_expand(trivec,ps) -> trivec", "svBtQps", "svBtQps", "trivec", "ps"},
             {"r_expand(trivec,quadvec) -> quadvec", "svBtQps", "svBtQps", "trivec",
              "quadvec"},
             // bivec
             {"r_expand(bivec,ps) -> bivec", "svBtQps", "svBtQps", "bivec", "ps"},
             {"r_expand(bivec,quadvec) -> trivec", "svBtQps", "svBtQps", "bivec",
              "quadvec"},
             {"r_expand(bivec,trivec) -> quadvec", "svBtQps", "svBtQps", "bivec",
              "trivec"},
             // vec
             {"r_expand(vec,ps) -> vec", "svBtQps", "svBtQps", "vec", "ps"},
             {"r_expand(vec,quadvec) -> bivec", "svBtQps", "svBtQps", "vec", "quadvec"},
             {"r_expand(vec,trivec) -> trivec", "svBtQps", "svBtQps", "vec", "trivec"},
             {"r_expand(vec,bivec) -> quadvec", "svBtQps", "svBtQps", "vec", "bivec"},
             // s
             {"r_expand(s,quadvec) -> vec", "svBtQps", "svBtQps", "s", "quadvec"},
             {"r_expand(s,trivec) -> bivec", "svBtQps", "svBtQps", "s", "trivec"},
             {"r_expand(s,bivec) -> trivec", "svBtQps", "svBtQps", "s", "bivec"},
             {"r_expand(s,vec) -> quadvec", "svBtQps", "svBtQps", "s", "vec"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_rgpr_config()
{
    return {
        .product_name = "rgpr",
        .description = "cga3dc regressive geometric product",
        .display_name = "regressive geometric product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// mv
             {"rgpr(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"rgpr(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
             {"rgpr(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
             {"rgpr(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
             {"rgpr(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
             {"rgpr(mv,ps) -> mv", "M", "svBtQps", "mv", "ps"},
             {"rgpr(ps,mv) -> mv", "svBtQps", "M", "ps", "mv"},
             {"rgpr(mv,quadvec) -> mv", "M", "svBtQps", "mv", "quadvec"},
             {"rgpr(quadvec,mv) -> mv", "svBtQps", "M", "quadvec", "mv"},
             {"rgpr(mv,trivec) -> mv", "M", "svBtQps", "mv", "trivec"},
             {"rgpr(trivec,mv) -> mv", "svBtQps", "M", "trivec", "mv"},
             {"rgpr(mv,bivec) -> mv", "M", "svBtQps", "mv", "bivec"},
             {"rgpr(bivec,mv) -> mv", "svBtQps", "M", "bivec", "mv"},
             {"rgpr(mv,vec) -> mv", "M", "svBtQps", "mv", "vec"},
             {"rgpr(vec,mv) -> mv", "svBtQps", "M", "vec", "mv"},
             {"rgpr(mv,s) -> mv", "M", "svBtQps", "mv", "s"},
             {"rgpr(s,mv) -> mv", "svBtQps", "M", "s", "mv"},
             // mv_e
             {"rgpr(mv_e,mv_e) -> mv_u", "A_even", "B_even", "mv_e", "mv_e"},
             {"rgpr(mv_e,mv_u) -> mv_e", "A_even", "B_odd", "mv_e", "mv_u"},
             {"rgpr(mv_u,mv_e) -> mv_e", "A_odd", "B_even", "mv_u", "mv_e"},
             {"rgpr(mv_e,ps) -> mv_e", "M_even", "svBtQps", "mv_e", "ps"},
             {"rgpr(ps,mv_e) -> mv_e", "svBtQps", "M_even", "ps", "mv_e"},
             {"rgpr(mv_e,quadvec) -> mv_u", "M_even", "svBtQps", "mv_e", "quadvec"},
             {"rgpr(quadvec,mv_e) -> mv_u", "svBtQps", "M_even", "quadvec", "mv_e"},
             {"rgpr(mv_e,trivec) -> mv_e", "M_even", "svBtQps", "mv_e", "trivec"},
             {"rgpr(trivec,mv_e) -> mv_e", "svBtQps", "M_even", "trivec", "mv_e"},
             {"rgpr(mv_e,bivec) -> mv_u", "M_even", "svBtQps", "mv_e", "bivec"},
             {"rgpr(bivec,mv_e) -> mv_u", "svBtQps", "M_even", "bivec", "mv_e"},
             {"rgpr(mv_e,vec) -> mv_e", "M_even", "svBtQps", "mv_e", "vec"},
             {"rgpr(vec,mv_e) -> mv_e", "svBtQps", "M_even", "vec", "mv_e"},
             {"rgpr(mv_e,s) -> mv_u", "M_even", "svBtQps", "mv_e", "s"},
             {"rgpr(s,mv_e) -> mv_u", "svBtQps", "M_even", "s", "mv_e"},
             // mv_u
             {"rgpr(mv_u,mv_u) -> mv_u", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"rgpr(mv_u,ps) -> mv_u", "M_odd", "svBtQps", "mv_u", "ps"},
             {"rgpr(ps,mv_u) -> mv_u", "svBtQps", "M_odd", "ps", "mv_u"},
             {"rgpr(mv_u,quadvec) -> mv_e", "M_odd", "svBtQps", "mv_u", "quadvec"},
             {"rgpr(quadvec,mv_u) -> mv_e", "svBtQps", "M_odd", "quadvec", "mv_u"},
             {"rgpr(mv_u,trivec) -> mv_u", "M_odd", "svBtQps", "mv_u", "trivec"},
             {"rgpr(trivec,mv_u) -> mv_u", "svBtQps", "M_odd", "trivec", "mv_u"},
             {"rgpr(mv_u,bivec) -> mv_e", "M_odd", "svBtQps", "mv_u", "bivec"},
             {"rgpr(bivec,mv_u) -> mv_e", "svBtQps", "M_odd", "bivec", "mv_u"},
             {"rgpr(mv_u,vec) -> mv_u", "M_odd", "svBtQps", "mv_u", "vec"},
             {"rgpr(vec,mv_u) -> mv_u", "svBtQps", "M_odd", "vec", "mv_u"},
             {"rgpr(mv_u,s) -> mv_e", "M_odd", "svBtQps", "mv_u", "s"},
             {"rgpr(s,mv_u) -> mv_e", "svBtQps", "M_odd", "s", "mv_u"},
             // ps
             {"rgpr(ps,ps) -> ps", "svBtQps1", "svBtQps2", "ps", "ps"},
             {"rgpr(ps,quadvec) -> quadvec", "svBtQps", "svBtQps", "ps", "quadvec"},
             {"rgpr(quadvec,ps) -> quadvec", "svBtQps", "svBtQps", "quadvec", "ps"},
             {"rgpr(ps,trivec) -> trivec", "svBtQps", "svBtQps", "ps", "trivec"},
             {"rgpr(trivec,ps) -> trivec", "svBtQps", "svBtQps", "trivec", "ps"},
             {"rgpr(ps,bivec) -> bivec", "svBtQps", "svBtQps", "ps", "bivec"},
             {"rgpr(bivec,ps) -> bivec", "svBtQps", "svBtQps", "bivec", "ps"},
             {"rgpr(ps,vec) -> vec", "svBtQps", "svBtQps", "ps", "vec"},
             {"rgpr(vec,ps) -> vec", "svBtQps", "svBtQps", "vec", "ps"},
             {"rgpr(ps,s) -> s", "svBtQps", "svBtQps", "ps", "s"},
             {"rgpr(s,ps) -> s", "svBtQps", "svBtQps", "s", "ps"},
             // quadvec
             {"rgpr(quadvec,quadvec) -> mv_u", "svBtQps1", "svBtQps2", "quadvec",
              "quadvec"},
             {"rgpr(quadvec,trivec) -> mv_e", "svBtQps", "svBtQps", "quadvec", "trivec"},
             {"rgpr(trivec,quadvec) -> mv_e", "svBtQps", "svBtQps", "trivec", "quadvec"},
             {"rgpr(quadvec,bivec) -> mv_u", "svBtQps", "svBtQps", "quadvec", "bivec"},
             {"rgpr(bivec,quadvec) -> mv_u", "svBtQps", "svBtQps", "bivec", "quadvec"},
             {"rgpr(quadvec,vec) -> mv_e", "svBtQps", "svBtQps", "quadvec", "vec"},
             {"rgpr(vec,quadvec) -> mv_e", "svBtQps", "svBtQps", "vec", "quadvec"},
             {"rgpr(quadvec,s) -> vec", "svBtQps", "svBtQps", "quadvec", "s"},
             {"rgpr(s,quadvec) -> vec", "svBtQps", "svBtQps", "s", "quadvec"},
             // trivec
             {"rgpr(trivec,trivec) -> mv_u", "svBtQps1", "svBtQps2", "trivec", "trivec"},
             {"rgpr(trivec,bivec) -> mv_e", "svBtQps", "svBtQps", "trivec", "bivec"},
             {"rgpr(bivec,trivec) -> mv_e", "svBtQps", "svBtQps", "bivec", "trivec"},
             {"rgpr(trivec,vec) -> mv_u", "svBtQps", "svBtQps", "trivec", "vec"},
             {"rgpr(vec,trivec) -> mv_u", "svBtQps", "svBtQps", "vec", "trivec"},
             {"rgpr(trivec,s) -> bivec", "svBtQps", "svBtQps", "trivec", "s"},
             {"rgpr(s,trivec) -> bivec", "svBtQps", "svBtQps", "s", "trivec"},
             // bivec
             {"rgpr(bivec,bivec) -> mv_u", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             {"rgpr(bivec,vec) -> mv_e", "svBtQps", "svBtQps", "bivec", "vec"},
             {"rgpr(vec,bivec) -> mv_e", "svBtQps", "svBtQps", "vec", "bivec"},
             {"rgpr(bivec,s) -> trivec", "svBtQps", "svBtQps", "bivec", "s"},
             {"rgpr(s,bivec) -> trivec", "svBtQps", "svBtQps", "s", "bivec"},
             // vec
             {"rgpr(vec,vec) -> mv_u", "svBtQps1", "svBtQps2", "vec", "vec"},
             {"rgpr(vec,s) -> quadvec", "svBtQps", "svBtQps", "vec", "s"},
             {"rgpr(s,vec) -> quadvec", "svBtQps", "svBtQps", "s", "vec"},
             // s
             {"rgpr(s,s) -> ps", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_rcmt_config()
{
    return {
        .product_name = "rcmt",
        .description = "cga3dc regressive commutator product",
        .display_name = "regressive commutator product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// mv
             {"rcmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
             // ps
             {"rcmt(ps,ps) -> 0", "svBtQps1", "svBtQps2", "ps", "ps"},
             {"rcmt(ps,quadvec) -> 0", "svBtQps", "svBtQps", "ps", "quadvec"},
             {"rcmt(quadvec,ps) -> 0", "svBtQps", "svBtQps", "quadvec", "ps"},
             {"rcmt(ps,trivec) -> 0", "svBtQps", "svBtQps", "ps", "trivec"},
             {"rcmt(trivec,ps) -> 0", "svBtQps", "svBtQps", "trivec", "ps"},
             {"rcmt(ps,bivec) -> 0", "svBtQps", "svBtQps", "ps", "bivec"},
             {"rcmt(bivec,ps) -> 0", "svBtQps", "svBtQps", "bivec", "ps"},
             {"rcmt(ps,vec) -> 0", "svBtQps", "svBtQps", "ps", "vec"},
             {"rcmt(vec,ps) -> 0", "svBtQps", "svBtQps", "vec", "ps"},
             {"rcmt(ps,s) -> 0", "svBtQps", "svBtQps", "ps", "s"},
             {"rcmt(s,ps) -> 0", "svBtQps", "svBtQps", "s", "ps"},
             // quadvec
             {"rcmt(quadvec,quadvec) -> trivec", "svBtQps1", "svBtQps2", "quadvec",
              "quadvec"},
             {"rcmt(quadvec,trivec) -> quadvec", "svBtQps", "svBtQps", "quadvec",
              "trivec"},
             {"rcmt(trivec,quadvec) -> quadvec", "svBtQps", "svBtQps", "trivec",
              "quadvec"},
             {"rcmt(quadvec,bivec) -> vec", "svBtQps", "svBtQps", "quadvec", "bivec"},
             {"rcmt(bivec,quadvec) -> vec", "svBtQps", "svBtQps", "bivec", "quadvec"},
             {"rcmt(quadvec,vec) -> bivec", "svBtQps", "svBtQps", "quadvec", "vec"},
             {"rcmt(vec,quadvec) -> bivec", "svBtQps", "svBtQps", "vec", "quadvec"},
             {"rcmt(quadvec,s) -> 0", "svBtQps", "svBtQps", "quadvec", "s"},
             {"rcmt(s,quadvec) -> 0", "svBtQps", "svBtQps", "s", "quadvec"},
             // trivec
             {"rcmt(trivec,trivec) -> trivec", "svBtQps1", "svBtQps2", "trivec",
              "trivec"},
             {"rcmt(trivec,bivec) -> bivec", "svBtQps", "svBtQps", "trivec", "bivec"},
             {"rcmt(bivec,trivec) -> bivec", "svBtQps", "svBtQps", "bivec", "trivec"},
             {"rcmt(trivec,vec) -> vec", "svBtQps", "svBtQps", "trivec", "vec"},
             {"rcmt(vec,trivec) -> vec", "svBtQps", "svBtQps", "vec", "trivec"},
             {"rcmt(trivec,s) -> 0", "svBtQps", "svBtQps", "trivec", "s"},
             {"rcmt(s,trivec) -> 0", "svBtQps", "svBtQps", "s", "trivec"},
             // bivec
             {"rcmt(bivec,bivec) -> trivec", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             {"rcmt(bivec,vec) -> quadvec", "svBtQps", "svBtQps", "bivec", "vec"},
             {"rcmt(vec,bivec) -> quadvec", "svBtQps", "svBtQps", "vec", "bivec"},
             {"rcmt(bivec,s) -> 0", "svBtQps", "svBtQps", "bivec", "s"},
             {"rcmt(s,bivec) -> 0", "svBtQps", "svBtQps", "s", "bivec"},
             // vec
             {"rcmt(vec,vec) -> trivec", "svBtQps1", "svBtQps2", "vec", "vec"},
             {"rcmt(vec,s) -> 0", "svBtQps", "svBtQps", "vec", "s"},
             {"rcmt(s,vec) -> 0", "svBtQps", "svBtQps", "s", "vec"},
             // s
             {"rcmt(s,s) -> 0", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_rwdg_config()
{
    return {
        .product_name = "rwdg",
        .description = "cga3dc regressive wedge product",
        .display_name = "regressive wedge product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// mv
             {"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv"},
             {"rwdg(mv,mv_e) -> mv", "A", "B_even", "mv", "mv_e"},
             {"rwdg(mv_e,mv) -> mv", "A_even", "B", "mv_e", "mv"},
             {"rwdg(mv,mv_u) -> mv", "A", "B_odd", "mv", "mv_u"},
             {"rwdg(mv_u,mv) -> mv", "A_odd", "B", "mv_u", "mv"},
             {"rwdg(mv,ps) -> mv", "M", "svBtQps", "mv", "ps"},
             {"rwdg(ps,mv) -> mv", "svBtQps", "M", "ps", "mv"},
             {"rwdg(mv,quadvec) -> mv", "M", "svBtQps", "mv", "quadvec"},
             {"rwdg(quadvec,mv) -> mv", "svBtQps", "M", "quadvec", "mv"},
             {"rwdg(mv,trivec) -> mv", "M", "svBtQps", "mv", "trivec"},
             {"rwdg(trivec,mv) -> mv", "svBtQps", "M", "trivec", "mv"},
             {"rwdg(mv,bivec) -> mv", "M", "svBtQps", "mv", "bivec"},
             {"rwdg(bivec,mv) -> mv", "svBtQps", "M", "bivec", "mv"},
             {"rwdg(mv,vec) -> mv", "M", "svBtQps", "mv", "vec"},
             {"rwdg(vec,mv) -> mv", "svBtQps", "M", "vec", "mv"},
             {"rwdg(mv,s) -> mv", "M", "svBtQps", "mv", "s"},
             {"rwdg(s,mv) -> mv", "svBtQps", "M", "s", "mv"},
             // mv_e
             {"rwdg(mv_e,mv_e) -> mv_u", "A_even", "B_even", "mv_e", "mv_e"},
             {"rwdg(mv_e,mv_u) -> mv_e", "A_even", "B_odd", "mv_e", "mv_u"},
             {"rwdg(mv_u,mv_e) -> mv_e", "A_odd", "B_even", "mv_u", "mv_e"},
             {"rwdg(mv_e,ps) -> mv_e", "M_even", "svBtQps", "mv_e", "ps"},
             {"rwdg(ps,mv_e) -> mv_e", "svBtQps", "M_even", "ps", "mv_e"},
             {"rwdg(mv_e,quadvec) -> mv_u", "M_even", "svBtQps", "mv_e", "quadvec"},
             {"rwdg(quadvec,mv_e) -> mv_u", "svBtQps", "M_even", "quadvec", "mv_e"},
             {"rwdg(mv_e,trivec) -> mv_e", "M_even", "svBtQps", "mv_e", "trivec"},
             {"rwdg(trivec,mv_e) -> mv_e", "svBtQps", "M_even", "trivec", "mv_e"},
             {"rwdg(mv_e,bivec) -> mv_u", "M_even", "svBtQps", "mv_e", "bivec"},
             {"rwdg(bivec,mv_e) -> mv_u", "svBtQps", "M_even", "bivec", "mv_e"},
             {"rwdg(mv_e,vec) -> mv_e", "M_even", "svBtQps", "mv_e", "vec"},
             {"rwdg(vec,mv_e) -> mv_e", "svBtQps", "M_even", "vec", "mv_e"},
             {"rwdg(mv_e,s) -> 0", "M_even", "svBtQps", "mv_e", "s"},
             {"rwdg(s,mv_e) -> 0", "svBtQps", "M_even", "s", "mv_e"},
             // mv_u
             {"rwdg(mv_u,mv_u) -> mv_u", "A_odd", "B_odd", "mv_u", "mv_u"},
             {"rwdg(mv_u,ps) -> mv_u", "M_odd", "svBtQps", "mv_u", "ps"},
             {"rwdg(ps,mv_u) -> mv_u", "svBtQps", "M_odd", "ps", "mv_u"},
             {"rwdg(mv_u,quadvec) -> mv_e", "M_odd", "svBtQps", "mv_u", "quadvec"},
             {"rwdg(quadvec,mv_u) -> mv_e", "svBtQps", "M_odd", "quadvec", "mv_u"},
             {"rwdg(mv_u,trivec) -> mv_u", "M_odd", "svBtQps", "mv_u", "trivec"},
             {"rwdg(trivec,mv_u) -> mv_u", "svBtQps", "M_odd", "trivec", "mv_u"},
             {"rwdg(mv_u,bivec) -> mv_e", "M_odd", "svBtQps", "mv_u", "bivec"},
             {"rwdg(bivec,mv_u) -> mv_e", "svBtQps", "M_odd", "bivec", "mv_u"},
             {"rwdg(mv_u,vec) -> mv_u", "M_odd", "svBtQps", "mv_u", "vec"},
             {"rwdg(vec,mv_u) -> mv_u", "svBtQps", "M_odd", "vec", "mv_u"},
             {"rwdg(mv_u,s) -> mv_e", "M_odd", "svBtQps", "mv_u", "s"},
             {"rwdg(s,mv_u) -> mv_e", "svBtQps", "M_odd", "s", "mv_u"},
             // ps
             {"rwdg(ps,ps) -> ps", "svBtQps1", "svBtQps2", "ps", "ps"},
             {"rwdg(ps,quadvec) -> quadvec", "svBtQps", "svBtQps", "ps", "quadvec"},
             {"rwdg(quadvec,ps) -> quadvec", "svBtQps", "svBtQps", "quadvec", "ps"},
             {"rwdg(ps,trivec) -> trivec", "svBtQps", "svBtQps", "ps", "trivec"},
             {"rwdg(trivec,ps) -> trivec", "svBtQps", "svBtQps", "trivec", "ps"},
             {"rwdg(ps,bivec) -> bivec", "svBtQps", "svBtQps", "ps", "bivec"},
             {"rwdg(bivec,ps) -> bivec", "svBtQps", "svBtQps", "bivec", "ps"},
             {"rwdg(ps,vec) -> vec", "svBtQps", "svBtQps", "ps", "vec"},
             {"rwdg(vec,ps) -> vec", "svBtQps", "svBtQps", "vec", "ps"},
             {"rwdg(ps,s) -> s", "svBtQps", "svBtQps", "ps", "s"},
             {"rwdg(s,ps) -> s", "svBtQps", "svBtQps", "s", "ps"},
             // quadvec
             {"rwdg(quadvec,quadvec) -> trivec", "svBtQps1", "svBtQps2", "quadvec",
              "quadvec"},
             {"rwdg(quadvec,trivec) -> bivec", "svBtQps", "svBtQps", "quadvec", "trivec"},
             {"rwdg(trivec,quadvec) -> bivec", "svBtQps", "svBtQps", "trivec", "quadvec"},
             {"rwdg(quadvec,bivec) -> vec", "svBtQps", "svBtQps", "quadvec", "bivec"},
             {"rwdg(bivec,quadvec) -> vec", "svBtQps", "svBtQps", "bivec", "quadvec"},
             {"rwdg(quadvec,vec) -> s", "svBtQps", "svBtQps", "quadvec", "vec"},
             {"rwdg(vec,quadvec) -> s", "svBtQps", "svBtQps", "vec", "quadvec"},
             {"rwdg(quadvec,s) -> 0", "svBtQps", "svBtQps", "quadvec", "s"},
             {"rwdg(s,quadvec) -> 0", "svBtQps", "svBtQps", "s", "quadvec"},
             // trivec
             {"rwdg(trivec,trivec) -> vec", "svBtQps1", "svBtQps2", "trivec", "trivec"},
             {"rwdg(trivec,bivec) -> s", "svBtQps", "svBtQps", "trivec", "bivec"},
             {"rwdg(bivec,trivec) -> s", "svBtQps", "svBtQps", "bivec", "trivec"},
             {"rwdg(trivec,vec) -> 0", "svBtQps", "svBtQps", "trivec", "vec"},
             {"rwdg(vec,trivec) -> 0", "svBtQps", "svBtQps", "vec", "trivec"},
             {"rwdg(trivec,s) -> 0", "svBtQps", "svBtQps", "trivec", "s"},
             {"rwdg(s,trivec) -> 0", "svBtQps", "svBtQps", "s", "trivec"},
             // bivec
             {"rwdg(bivec,bivec) -> 0", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             {"rwdg(bivec,vec) -> 0", "svBtQps", "svBtQps", "bivec", "vec"},
             {"rwdg(vec,bivec) -> 0", "svBtQps", "svBtQps", "vec", "bivec"},
             {"rwdg(bivec,s) -> 0", "svBtQps", "svBtQps", "bivec", "s"},
             {"rwdg(s,bivec) -> 0", "svBtQps", "svBtQps", "s", "bivec"},
             // vec
             {"rwdg(vec,vec) -> 0", "svBtQps1", "svBtQps2", "vec", "vec"},
             {"rwdg(vec,s) -> 0", "svBtQps", "svBtQps", "vec", "s"},
             {"rwdg(s,vec) -> 0", "svBtQps", "svBtQps", "s", "vec"},
             // s
             {"rwdg(s,s) -> 0", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

ProductConfig get_cga3dc_rdot_config()
{
    return {
        .product_name = "rdot",
        .description = "cga3dc regressive inner product",
        .display_name = "regressive inner product",
        // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
        // "right_filter"}
        .cases =
            {// mv
             {"rdot(mv,mv) -> ps", "A", "B", "mv", "mv"},
             {"rdot(mv_e,mv_e) -> ps", "A_even", "B_even", "mv_e", "mv_e"},
             {"rdot(mv_u,mv_u) -> ps", "A_odd", "B_odd", "mv_u", "mv_u"},
             // ps
             {"rdot(ps,ps) -> ps", "svBtQps1", "svBtQps2", "ps", "ps"},
             // quadvec
             {"rdot(quadvec,quadvec) -> ps", "svBtQps1", "svBtQps2", "quadvec",
              "quadvec"},
             // trivec
             {"rdot(trivec,trivec) -> ps", "svBtQps1", "svBtQps2", "trivec", "trivec"},
             // bivec
             {"rdot(bivec,bivec) -> ps", "svBtQps1", "svBtQps2", "bivec", "bivec"},
             // vec
             {"rdot(vec,vec) -> ps", "svBtQps1", "svBtQps2", "vec", "vec"},
             // s
             {"rdot(s,s) -> ps", "svBtQps1", "svBtQps2", "s", "s"}},
        .is_sandwich_product = false,
        .uses_brace_switch = false,
        .show_basis_table = true};
}

} // namespace configurable
