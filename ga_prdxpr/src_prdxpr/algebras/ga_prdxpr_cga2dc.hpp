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

////////////////////////////////////////////////////////////////////////////////
// coefficients needed to create a multivector = [coeff]^T [mv2dc_basis]
//
// component layout (16): [0] scalar, [1..4] vector (x,y,w,u roles in fields
// x,y,z,w), [5..10] bivector (vx,vy,vz,mx,my,mz), [11..14] trivector
// (x,y,z,w fields), [15] pseudoscalar. Even part: 0,5..10,15; odd: 1..4,11..14.
////////////////////////////////////////////////////////////////////////////////

const mvec_coeff mv2dc_coeff_A = {"A.c0",  "A.c1",  "A.c2",  "A.c3", "A.c4",  "A.c5",
                                  "A.c6",  "A.c7",  "A.c8",  "A.c9", "A.c10", "A.c11",
                                  "A.c12", "A.c13", "A.c14", "A.c15"};
const mvec_coeff mv2dc_coeff_B = {"B.c0",  "B.c1",  "B.c2",  "B.c3", "B.c4",  "B.c5",
                                  "B.c6",  "B.c7",  "B.c8",  "B.c9", "B.c10", "B.c11",
                                  "B.c12", "B.c13", "B.c14", "B.c15"};
const mvec_coeff mv2dc_coeff_M = {"M.c0",  "M.c1",  "M.c2",  "M.c3", "M.c4",  "M.c5",
                                  "M.c6",  "M.c7",  "M.c8",  "M.c9", "M.c10", "M.c11",
                                  "M.c12", "M.c13", "M.c14", "M.c15"};
const mvec_coeff mv2dc_coeff_M1 = {
    "M1.c0", "M1.c1", "M1.c2",  "M1.c3",  "M1.c4",  "M1.c5",  "M1.c6",  "M1.c7",
    "M1.c8", "M1.c9", "M1.c10", "M1.c11", "M1.c12", "M1.c13", "M1.c14", "M1.c15"};
const mvec_coeff mv2dc_coeff_M2 = {
    "M2.c0", "M2.c1", "M2.c2",  "M2.c3",  "M2.c4",  "M2.c5",  "M2.c6",  "M2.c7",
    "M2.c8", "M2.c9", "M2.c10", "M2.c11", "M2.c12", "M2.c13", "M2.c14", "M2.c15"};
const mvec_coeff mv2dc_coeff_M_even = {"M.c0", "0",    "0",    "0",    "0",    "M.c1",
                                       "M.c2", "M.c3", "M.c4", "M.c5", "M.c6", "0",
                                       "0",    "0",    "0",    "M.c7"};
const mvec_coeff mv2dc_coeff_M_odd = {"0",    "M.c0", "M.c1", "M.c2", "M.c3", "0",
                                      "0",    "0",    "0",    "0",    "0",    "M.c4",
                                      "M.c5", "M.c6", "M.c7", "0"};
const mvec_coeff mv2dc_coeff_A_even = {"A.c0", "0",    "0",    "0",    "0",    "A.c1",
                                       "A.c2", "A.c3", "A.c4", "A.c5", "A.c6", "0",
                                       "0",    "0",    "0",    "A.c7"};
const mvec_coeff mv2dc_coeff_B_even = {"B.c0", "0",    "0",    "0",    "0",    "B.c1",
                                       "B.c2", "B.c3", "B.c4", "B.c5", "B.c6", "0",
                                       "0",    "0",    "0",    "B.c7"};
const mvec_coeff mv2dc_coeff_A_odd = {"0",    "A.c0", "A.c1", "A.c2", "A.c3", "0",
                                      "0",    "0",    "0",    "0",    "0",    "A.c4",
                                      "A.c5", "A.c6", "A.c7", "0"};
const mvec_coeff mv2dc_coeff_B_odd = {"0",    "B.c0", "B.c1", "B.c2", "B.c3", "0",
                                      "0",    "0",    "0",    "0",    "0",    "B.c4",
                                      "B.c5", "B.c6", "B.c7", "0"};
const mvec_coeff mv2dc_coeff_R_even = {"R.c0", "0",    "0",    "0",    "0",    "R.c1",
                                       "R.c2", "R.c3", "R.c4", "R.c5", "R.c6", "0",
                                       "0",    "0",    "0",    "R.c7"};
const mvec_coeff mv2dc_coeff_R_rev_even = {
    "R.c0",  "0",     "0",     "0", "0", "-R.c1", "-R.c2", "-R.c3",
    "-R.c4", "-R.c5", "-R.c6", "0", "0", "0",     "0",     "R.c7"};
const mvec_coeff mv2dc_coeff_R_rrev_even = {
    "R.c0",  "0",     "0",     "0", "0", "-R.c1", "-R.c2", "-R.c3",
    "-R.c4", "-R.c5", "-R.c6", "0", "0", "0",     "0",     "R.c7"};
const mvec_coeff mv2dc_coeff_R_odd = {"0",    "R.c0", "R.c1", "R.c2", "R.c3", "0",
                                      "0",    "0",    "0",    "0",    "0",    "R.c4",
                                      "R.c5", "R.c6", "R.c7", "0"};
const mvec_coeff mv2dc_coeff_R_rev_odd = {
    "0", "R.c0", "R.c1", "R.c2",  "R.c3",  "0",     "0",     "0",
    "0", "0",    "0",    "-R.c4", "-R.c5", "-R.c6", "-R.c7", "0"};
const mvec_coeff mv2dc_coeff_R_rrev_odd = {
    "0", "-R.c0", "-R.c1", "-R.c2", "-R.c3", "0",    "0",    "0",
    "0", "0",     "0",     "R.c4",  "R.c5",  "R.c6", "R.c7", "0"};
const mvec_coeff mv2dc_coeff_svBtps = {"s",    "v.x",  "v.y",  "v.z",  "v.w",  "B.vx",
                                       "B.vy", "B.vz", "B.mx", "B.my", "B.mz", "t.x",
                                       "t.y",  "t.z",  "t.w",  "ps"};
const mvec_coeff mv2dc_coeff_svBtps1 = {
    "s1",    "v1.x",  "v1.y",  "v1.z", "v1.w", "B1.vx", "B1.vy", "B1.vz",
    "B1.mx", "B1.my", "B1.mz", "t1.x", "t1.y", "t1.z",  "t1.w",  "ps1"};
const mvec_coeff mv2dc_coeff_svBtps2 = {
    "s2",    "v2.x",  "v2.y",  "v2.z", "v2.w", "B2.vx", "B2.vy", "B2.vz",
    "B2.mx", "B2.my", "B2.mz", "t2.x", "t2.y", "t2.z",  "t2.w",  "ps2"};
