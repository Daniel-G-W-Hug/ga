#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"
#include "ga_prdxpr_rule_generator.hpp"

// ega2d: euclidean geometric algebra with
//        e1*e1 = +1, e2*e2 = +1
// G(2,0,0)

// multivector basis components
const mvec_coeff mv2d_basis = {"1", "e1", "e2", "e12"};

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - EGA2D
//
// Product tables are generated automatically using the rule generation system:
// - Configuration: ga_prdxpr_ega2d_config.cpp (defines algebra parameters)
// - Generation: ga_prdxpr_ega2d.cpp (creates product rules at startup)
// - Interface: This file (declares available rules)
//
// This ensures mathematically correct, consistent rules without manual errors.
////////////////////////////////////////////////////////////////////////////////
extern const prd_rules gpr_ega2d_rules; // Geometric product: a * b
extern const prd_rules wdg_ega2d_rules; // Wedge product: a ^ b
extern const prd_rules dot_ega2d_rules; // Dot product: a · b


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
