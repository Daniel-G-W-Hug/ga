#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "core/ga_prdxpr_common.hpp"

////////////////////////////////////////////////////////////////////////////////
// cga3dc: conformal geometric algebra of 3D Euclidean space
//
// 5D algebra over the basis vectors e1, e2, e3 (Euclidean), e4 (projective,
// w), e5 (round, u). e4 and e5 form a NULL PAIR: e4*e4 = e5*e5 = 0,
// e4.e5 = -1 — a non-orthogonal metric carried by AlgebraConfig::metric_matrix.
//
// Consequence: the geometric product of basis blades is MULTI-TERM
// (e.g. e4 * e5 = -1 + e45), so gpr rules live in prd_rules_mt; wedge, dot,
// complements and duals stay single-term (see the rule generator docs).
//
// NOTE (5D): reversion and regressive reversion coincide in five dimensions
// (the sign exponents agree for every grade), so the R_rev_* and R_rrev_*
// coefficient sets below are identical by construction.
////////////////////////////////////////////////////////////////////////////////

// multivector basis: spatial components first, then the projective dimension
// (w = e4), then the round dimension (u = e5, always written trailing).
// Grades 0-2 = [pga3dp grade-k blades, pga3dp order] + [pga3dp grade-(k-1)
// blades ^ e5, pga3dp order]. Grade 3 = the conformal LINE's blades first
// (v = e415,e425,e435; m = e235,e315,e125), then p = e423,e431,e412 and the
// pure-Euclidean blade e321 LAST. Grade 4 = label order x,y,z,w,u by omitted
// basis vector (= pga3dp trivectors ^ e5, then e1234).
const mvec_coeff mv3dc_basis = {
    "1",    "e1",   "e2",    "e3",    "e4",    "e5",    "e41",   "e42",
    "e43",  "e23",  "e31",   "e12",   "e15",   "e25",   "e35",   "e45",
    "e415", "e425", "e435",  "e235",  "e315",  "e125",  "e423",  "e431",
    "e412", "e321", "e4235", "e4315", "e4125", "e3215", "e1234", "e12345"};

const std::vector<mvec_coeff> mv3dc_basis_kvec = {
    {"1"},
    {"e1", "e2", "e3", "e4", "e5"},
    {"e41", "e42", "e43", "e23", "e31", "e12", "e15", "e25", "e35", "e45"},
    {"e415", "e425", "e435", "e235", "e315", "e125", "e423", "e431", "e412", "e321"},
    {"e4235", "e4315", "e4125", "e3215", "e1234"},
    {"e12345"}};

////////////////////////////////////////////////////////////////////////////////
// product rules (generated at static-init time in ga_prdxpr_cga3dc.cpp)
////////////////////////////////////////////////////////////////////////////////

// geometric product: MULTI-TERM (null-pair metric) — a * b
extern const prd_rules_mt gpr_cga3dc_rules_mt;

// single-term products
extern const prd_rules wdg_cga3dc_rules; // Wedge product: a ^ b
extern const prd_rules dot_cga3dc_rules; // Dot product: a · b

// complements (odd-dimensional algebra: left and right complements coincide,
// but the l_/r_ pair is kept so the generator's even/odd-dimension machinery
// stays uniform across algebras)
extern const prd_rules l_cmpl_cga3dc_rules;
extern const prd_rules r_cmpl_cga3dc_rules;

// duals (non-degenerate metric: regular left/right duals, no bulk/weight split)
extern const prd_rules l_dual_cga3dc_rules;
extern const prd_rules r_dual_cga3dc_rules;

// antiduals (dual w.r.t. the anti-exomorphism = the negated duals;
// the identity is verified at rule-generation time)
extern const prd_rules l_antidual_cga3dc_rules;
extern const prd_rules r_antidual_cga3dc_rules;

////////////////////////////////////////////////////////////////////////////////
// coefficients needed to create a multivector = [coeff]^T [mv3dc_basis]
//
// component layout (32): [0] scalar, [1..5] vector (x,y,z,w,u), [6..15]
// bivector (vx,vy,vz,mx,my,mz,px,py,pz,pw), [16..25] trivector (same field
// names), [26..30] quadvector (x,y,z,w,u), [31] pseudoscalar.
// Even part: 0, 6..15, 26..30; odd part: 1..5, 16..25, 31.
////////////////////////////////////////////////////////////////////////////////

const mvec_coeff mv3dc_coeff_A = {
    "A.c0",  "A.c1",  "A.c2",  "A.c3",  "A.c4",  "A.c5",  "A.c6",  "A.c7",
    "A.c8",  "A.c9",  "A.c10", "A.c11", "A.c12", "A.c13", "A.c14", "A.c15",
    "A.c16", "A.c17", "A.c18", "A.c19", "A.c20", "A.c21", "A.c22", "A.c23",
    "A.c24", "A.c25", "A.c26", "A.c27", "A.c28", "A.c29", "A.c30", "A.c31"};
const mvec_coeff mv3dc_coeff_B = {
    "B.c0",  "B.c1",  "B.c2",  "B.c3",  "B.c4",  "B.c5",  "B.c6",  "B.c7",
    "B.c8",  "B.c9",  "B.c10", "B.c11", "B.c12", "B.c13", "B.c14", "B.c15",
    "B.c16", "B.c17", "B.c18", "B.c19", "B.c20", "B.c21", "B.c22", "B.c23",
    "B.c24", "B.c25", "B.c26", "B.c27", "B.c28", "B.c29", "B.c30", "B.c31"};
const mvec_coeff mv3dc_coeff_M = {
    "M.c0",  "M.c1",  "M.c2",  "M.c3",  "M.c4",  "M.c5",  "M.c6",  "M.c7",
    "M.c8",  "M.c9",  "M.c10", "M.c11", "M.c12", "M.c13", "M.c14", "M.c15",
    "M.c16", "M.c17", "M.c18", "M.c19", "M.c20", "M.c21", "M.c22", "M.c23",
    "M.c24", "M.c25", "M.c26", "M.c27", "M.c28", "M.c29", "M.c30", "M.c31"};
const mvec_coeff mv3dc_coeff_M1 = {
    "M1.c0",  "M1.c1",  "M1.c2",  "M1.c3",  "M1.c4",  "M1.c5",  "M1.c6",  "M1.c7",
    "M1.c8",  "M1.c9",  "M1.c10", "M1.c11", "M1.c12", "M1.c13", "M1.c14", "M1.c15",
    "M1.c16", "M1.c17", "M1.c18", "M1.c19", "M1.c20", "M1.c21", "M1.c22", "M1.c23",
    "M1.c24", "M1.c25", "M1.c26", "M1.c27", "M1.c28", "M1.c29", "M1.c30", "M1.c31"};
const mvec_coeff mv3dc_coeff_M2 = {
    "M2.c0",  "M2.c1",  "M2.c2",  "M2.c3",  "M2.c4",  "M2.c5",  "M2.c6",  "M2.c7",
    "M2.c8",  "M2.c9",  "M2.c10", "M2.c11", "M2.c12", "M2.c13", "M2.c14", "M2.c15",
    "M2.c16", "M2.c17", "M2.c18", "M2.c19", "M2.c20", "M2.c21", "M2.c22", "M2.c23",
    "M2.c24", "M2.c25", "M2.c26", "M2.c27", "M2.c28", "M2.c29", "M2.c30", "M2.c31"};
const mvec_coeff mv3dc_coeff_M_even = {
    "M.c0", "0",    "0",     "0",     "0",     "0",     "M.c1",  "M.c2",
    "M.c3", "M.c4", "M.c5",  "M.c6",  "M.c7",  "M.c8",  "M.c9",  "M.c10",
    "0",    "0",    "0",     "0",     "0",     "0",     "0",     "0",
    "0",    "0",    "M.c11", "M.c12", "M.c13", "M.c14", "M.c15", "0"};
const mvec_coeff mv3dc_coeff_M_odd = {
    "0",     "M.c0",  "M.c1", "M.c2", "M.c3", "M.c4",  "0",     "0",
    "0",     "0",     "0",    "0",    "0",    "0",     "0",     "0",
    "M.c5",  "M.c6",  "M.c7", "M.c8", "M.c9", "M.c10", "M.c11", "M.c12",
    "M.c13", "M.c14", "0",    "0",    "0",    "0",     "0",     "M.c15"};
const mvec_coeff mv3dc_coeff_A_even = {
    "A.c0", "0",    "0",     "0",     "0",     "0",     "A.c1",  "A.c2",
    "A.c3", "A.c4", "A.c5",  "A.c6",  "A.c7",  "A.c8",  "A.c9",  "A.c10",
    "0",    "0",    "0",     "0",     "0",     "0",     "0",     "0",
    "0",    "0",    "A.c11", "A.c12", "A.c13", "A.c14", "A.c15", "0"};
const mvec_coeff mv3dc_coeff_B_even = {
    "B.c0", "0",    "0",     "0",     "0",     "0",     "B.c1",  "B.c2",
    "B.c3", "B.c4", "B.c5",  "B.c6",  "B.c7",  "B.c8",  "B.c9",  "B.c10",
    "0",    "0",    "0",     "0",     "0",     "0",     "0",     "0",
    "0",    "0",    "B.c11", "B.c12", "B.c13", "B.c14", "B.c15", "0"};
const mvec_coeff mv3dc_coeff_A_odd = {
    "0",     "A.c0",  "A.c1", "A.c2", "A.c3", "A.c4",  "0",     "0",
    "0",     "0",     "0",    "0",    "0",    "0",     "0",     "0",
    "A.c5",  "A.c6",  "A.c7", "A.c8", "A.c9", "A.c10", "A.c11", "A.c12",
    "A.c13", "A.c14", "0",    "0",    "0",    "0",     "0",     "A.c15"};
const mvec_coeff mv3dc_coeff_B_odd = {
    "0",     "B.c0",  "B.c1", "B.c2", "B.c3", "B.c4",  "0",     "0",
    "0",     "0",     "0",    "0",    "0",    "0",     "0",     "0",
    "B.c5",  "B.c6",  "B.c7", "B.c8", "B.c9", "B.c10", "B.c11", "B.c12",
    "B.c13", "B.c14", "0",    "0",    "0",    "0",     "0",     "B.c15"};
const mvec_coeff mv3dc_coeff_R_even = {
    "R.c0", "0",    "0",     "0",     "0",     "0",     "R.c1",  "R.c2",
    "R.c3", "R.c4", "R.c5",  "R.c6",  "R.c7",  "R.c8",  "R.c9",  "R.c10",
    "0",    "0",    "0",     "0",     "0",     "0",     "0",     "0",
    "0",    "0",    "R.c11", "R.c12", "R.c13", "R.c14", "R.c15", "0"};
const mvec_coeff mv3dc_coeff_R_rev_even = {
    "R.c0",  "0",     "0",     "0",     "0",     "0",     "-R.c1", "-R.c2",
    "-R.c3", "-R.c4", "-R.c5", "-R.c6", "-R.c7", "-R.c8", "-R.c9", "-R.c10",
    "0",     "0",     "0",     "0",     "0",     "0",     "0",     "0",
    "0",     "0",     "R.c11", "R.c12", "R.c13", "R.c14", "R.c15", "0"};
const mvec_coeff mv3dc_coeff_R_rrev_even = {
    "R.c0",  "0",     "0",     "0",     "0",     "0",     "-R.c1", "-R.c2",
    "-R.c3", "-R.c4", "-R.c5", "-R.c6", "-R.c7", "-R.c8", "-R.c9", "-R.c10",
    "0",     "0",     "0",     "0",     "0",     "0",     "0",     "0",
    "0",     "0",     "R.c11", "R.c12", "R.c13", "R.c14", "R.c15", "0"};
const mvec_coeff mv3dc_coeff_R_odd = {
    "0",     "R.c0",  "R.c1", "R.c2", "R.c3", "R.c4",  "0",     "0",
    "0",     "0",     "0",    "0",    "0",    "0",     "0",     "0",
    "R.c5",  "R.c6",  "R.c7", "R.c8", "R.c9", "R.c10", "R.c11", "R.c12",
    "R.c13", "R.c14", "0",    "0",    "0",    "0",     "0",     "R.c15"};
const mvec_coeff mv3dc_coeff_R_rev_odd = {
    "0",      "R.c0",   "R.c1",  "R.c2",  "R.c3",  "R.c4",   "0",      "0",
    "0",      "0",      "0",     "0",     "0",     "0",      "0",      "0",
    "-R.c5",  "-R.c6",  "-R.c7", "-R.c8", "-R.c9", "-R.c10", "-R.c11", "-R.c12",
    "-R.c13", "-R.c14", "0",     "0",     "0",     "0",      "0",      "R.c15"};
const mvec_coeff mv3dc_coeff_R_rrev_odd = {
    "0",      "R.c0",   "R.c1",  "R.c2",  "R.c3",  "R.c4",   "0",      "0",
    "0",      "0",      "0",     "0",     "0",     "0",      "0",      "0",
    "-R.c5",  "-R.c6",  "-R.c7", "-R.c8", "-R.c9", "-R.c10", "-R.c11", "-R.c12",
    "-R.c13", "-R.c14", "0",     "0",     "0",     "0",      "0",      "R.c15"};
const mvec_coeff mv3dc_coeff_svBtQps = {
    "s",    "v.x",  "v.y",  "v.z",  "v.w",  "v.u",  "B.vx", "B.vy",
    "B.vz", "B.mx", "B.my", "B.mz", "B.px", "B.py", "B.pz", "B.pw",
    "t.vx", "t.vy", "t.vz", "t.mx", "t.my", "t.mz", "t.px", "t.py",
    "t.pz", "t.pw", "Q.x",  "Q.y",  "Q.z",  "Q.w",  "Q.u",  "ps"};
const mvec_coeff mv3dc_coeff_svBtQps1 = {
    "s1",    "v1.x",  "v1.y",  "v1.z",  "v1.w",  "v1.u",  "B1.vx", "B1.vy",
    "B1.vz", "B1.mx", "B1.my", "B1.mz", "B1.px", "B1.py", "B1.pz", "B1.pw",
    "t1.vx", "t1.vy", "t1.vz", "t1.mx", "t1.my", "t1.mz", "t1.px", "t1.py",
    "t1.pz", "t1.pw", "Q1.x",  "Q1.y",  "Q1.z",  "Q1.w",  "Q1.u",  "ps1"};
const mvec_coeff mv3dc_coeff_svBtQps2 = {
    "s2",    "v2.x",  "v2.y",  "v2.z",  "v2.w",  "v2.u",  "B2.vx", "B2.vy",
    "B2.vz", "B2.mx", "B2.my", "B2.mz", "B2.px", "B2.py", "B2.pz", "B2.pw",
    "t2.vx", "t2.vy", "t2.vz", "t2.mx", "t2.my", "t2.mz", "t2.px", "t2.py",
    "t2.pz", "t2.pw", "Q2.x",  "Q2.y",  "Q2.z",  "Q2.w",  "Q2.u",  "ps2"};
