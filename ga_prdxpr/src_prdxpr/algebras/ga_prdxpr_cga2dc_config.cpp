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

} // namespace configurable
