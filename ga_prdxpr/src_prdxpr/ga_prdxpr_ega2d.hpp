#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"

// ega2d: euclidean geometric algebra with
//        e1*e1 = +1, e2*e2 = +1
// G(2,0,0)

// multivector basis components
const mvec_coeff mv2d_basis = {"1", "e1", "e2", "e12"};

// rules to create produc tables (must contain space_str around operator_str on each side)
const prd_rules gpr_ega2d_rules = {
    {"1 * 1", "1"},     {"1 * e1", "e1"},    {"1 * e2", "e2"},   {"1 * e12", "e12"},
    {"e1 * 1", "e1"},   {"e1 * e1", "1"},    {"e1 * e2", "e12"}, {"e1 * e12", "e2"},
    {"e2 * 1", "e2"},   {"e2 * e1", "-e12"}, {"e2 * e2", "1"},   {"e2 * e12", "-e1"},
    {"e12 * 1", "e12"}, {"e12 * e1", "-e2"}, {"e12 * e2", "e1"}, {"e12 * e12", "-1"}};

const prd_rules wdg_ega2d_rules = {
    {"1 ^ 1", "1"},     {"1 ^ e1", "e1"},    {"1 ^ e2", "e2"},   {"1 ^ e12", "e12"},
    {"e1 ^ 1", "e1"},   {"e1 ^ e1", "0"},    {"e1 ^ e2", "e12"}, {"e1 ^ e12", "0"},
    {"e2 ^ 1", "e2"},   {"e2 ^ e1", "-e12"}, {"e2 ^ e2", "0"},   {"e2 ^ e12", "0"},
    {"e12 ^ 1", "e12"}, {"e12 ^ e1", "0"},   {"e12 ^ e2", "0"},  {"e12 ^ e12", "0"}};

const prd_rules dot_ega2d_rules = {
    {"1 * 1", "1"},   {"1 * e1", "0"},   {"1 * e2", "0"},   {"1 * e12", "0"},
    {"e1 * 1", "0"},  {"e1 * e1", "1"},  {"e1 * e2", "0"},  {"e1 * e12", "0"},
    {"e2 * 1", "0"},  {"e2 * e1", "0"},  {"e2 * e2", "1"},  {"e2 * e12", "0"},
    {"e12 * 1", "0"}, {"e12 * e1", "0"}, {"e12 * e2", "0"}, {"e12 * e12", "1"}};


// complement in 2d: lcmpl(rcmpl(u)) = u
// lcmpl:  lcmpl(u) ^ u = e12
const prd_rules lcmpl_ega2d_rules = {
    {"1", "e12"}, {"e1", "-e2"}, {"e2", "e1"}, {"e12", "1"}};

// rcmpl:  u ^ rcmpl(u) = e12
const prd_rules rcmpl_ega2d_rules = {
    {"1", "e12"}, {"e1", "e2"}, {"e2", "-e1"}, {"e12", "1"}};

// complement and dual are identical in ega2d (the metric is the identity matrix)
const prd_rules left_dual_ega2d_rules = lcmpl_ega2d_rules;
const prd_rules right_dual_ega2d_rules = rcmpl_ega2d_rules;

// coefficients needed to create a multivector = [coeff]^T [mv2d_basis]
const mvec_coeff mv2d_coeff_A = {"A.c0", "A.c1", "A.c2", "A.c3"};
const mvec_coeff mv2d_coeff_B = {"B.c0", "B.c1", "B.c2", "B.c3"};
const mvec_coeff mv2d_coeff_M = {"M.c0", "M.c1", "M.c2", "M.c3"};
const mvec_coeff mv2d_coeff_M1 = {"M1.c0", "M1.c1", "M1.c2", "M1.c3"};
const mvec_coeff mv2d_coeff_M2 = {"M2.c0", "M2.c1", "M2.c2", "M2.c3"};
const mvec_coeff mv2d_coeff_A_even = {"A.c0", "0", "0", "A.c1"};
const mvec_coeff mv2d_coeff_B_even = {"B.c0", "0", "0", "B.c1"};
const mvec_coeff mv2d_coeff_M_even = {"M.c0", "0", "0", "M.c1"};
const mvec_coeff mv2d_coeff_R_even = {"R.c0", "0", "0", "R.c1"};
const mvec_coeff mv2d_coeff_R_rev_even = {"R.c0", "0", "0", "-R.c1"};
const mvec_coeff mv2d_coeff_svps = {"s", "v.x", "v.y", "ps"};
const mvec_coeff mv2d_coeff_svps1 = {"s1", "v1.x", "v1.y", "ps1"};
const mvec_coeff mv2d_coeff_svps2 = {"s2", "v2.x", "v2.y", "ps2"};
