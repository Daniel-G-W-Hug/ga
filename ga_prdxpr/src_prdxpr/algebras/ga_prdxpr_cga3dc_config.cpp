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

} // namespace configurable
