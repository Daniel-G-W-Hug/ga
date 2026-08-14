#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "core/ga_prdxpr_common.hpp"

////////////////////////////////////////////////////////////////////////////////
// cga2dc: conformal geometric algebra of 2D Euclidean space
//
// 4D algebra over the basis vectors e1, e2 (Euclidean), e3 (projective, w),
// e4 (round, u). e3 and e4 form a NULL PAIR: e3*e3 = e4*e4 = 0, e3.e4 = -1 —
// a non-orthogonal metric carried by AlgebraConfig::metric_matrix.
//
// Consequence: the geometric product of basis blades is MULTI-TERM
// (e.g. e3 * e4 = -1 + e34), so gpr rules live in prd_rules_mt; wedge, dot,
// complements and duals stay single-term (see the rule generator docs).
////////////////////////////////////////////////////////////////////////////////

// multivector basis: spatial components first, then the projective dimension
// (w = e3), then the round dimension (u = e4, always written trailing).
// Bivectors: the pga2dp bivectors {e31, e32, e12} first (the LIBRARY's pga2dp
// embeds verbatim), then grade-1 ∧ e4. Trivectors: pga2dp bivectors ∧ e4 in
// pga2dp order (a pga2dp line's components land verbatim in fields x,y,z),
// then e321 (the pga2dp pseudoscalar) last.
const mvec_coeff mv2dc_basis = {"1",    "e1",   "e2",   "e3",   "e4",  "e31",
                                "e32",  "e12",  "e14",  "e24",  "e34", "e314",
                                "e324", "e124", "e321", "e1234"};

const std::vector<mvec_coeff> mv2dc_basis_kvec = {
    {"1"},
    {"e1", "e2", "e3", "e4"},
    {"e31", "e32", "e12", "e14", "e24", "e34"},
    {"e314", "e324", "e124", "e321"},
    {"e1234"}};

////////////////////////////////////////////////////////////////////////////////
// product rules (generated at static-init time in ga_prdxpr_cga2dc.cpp)
////////////////////////////////////////////////////////////////////////////////

// geometric product: MULTI-TERM (null-pair metric) — a * b
extern const prd_rules_mt gpr_cga2dc_rules_mt;

// single-term products
extern const prd_rules wdg_cga2dc_rules; // Wedge product: a ^ b
extern const prd_rules dot_cga2dc_rules; // Dot product: a · b

// complements (even-dimensional algebra: left and right)
extern const prd_rules l_cmpl_cga2dc_rules;
extern const prd_rules r_cmpl_cga2dc_rules;

// duals (non-degenerate metric: regular left/right duals, no bulk/weight split)
extern const prd_rules l_dual_cga2dc_rules;
extern const prd_rules r_dual_cga2dc_rules;
