#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"
#include "ga_prdxpr_rule_generator.hpp"

// ega3d: euclidean geometric algebra with
//        e1*e1 = +1, e2*e2 = +1, e3*e3 = +1
// G(3,0,0)

// multivector basis components
const mvec_coeff mv3d_basis = {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e123"};

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - EGA3D
//
// Product tables are generated automatically using the rule generation system:
// - Configuration: ga_prdxpr_ega3d_config.cpp (defines algebra parameters)
// - Generation: ga_prdxpr_ega3d.cpp (creates product rules at startup)
// - Interface: This file (declares available rules)
//
// This ensures mathematically correct, consistent rules without manual errors.
////////////////////////////////////////////////////////////////////////////////
extern const prd_rules gpr_ega3d_rules; // Geometric product: a * b
extern const prd_rules wdg_ega3d_rules; // Wedge product: a ^ b
extern const prd_rules dot_ega3d_rules; // Dot product: a · b

// complement in 3d: cmpl(cmpl(u)) = u
// lcmpl(u) ^ u = u ^ rcmpl(u) = cmpl(u) ^ u = u ^ cmpl(u) = e123
extern const prd_rules cmpl_ega3d_rules;

// complement and dual are identical in ega3d (the metric is the identity matrix)
const prd_rules dual_ega3d_rules = {{"1", "e123"}, {"e1", "e23"}, {"e2", "e31"},
                                    {"e3", "e12"}, {"e23", "e1"}, {"e31", "e2"},
                                    {"e12", "e3"}, {"e123", "1"}};

// coefficients needed to create a multivector = [coeff]^T [mv3d_basis]
const mvec_coeff mv3d_coeff_A = {"A.c0", "A.c1", "A.c2", "A.c3",
                                 "A.c4", "A.c5", "A.c6", "A.c7"};
const mvec_coeff mv3d_coeff_B = {"B.c0", "B.c1", "B.c2", "B.c3",
                                 "B.c4", "B.c5", "B.c6", "B.c7"};
const mvec_coeff mv3d_coeff_M = {"M.c0", "M.c1", "M.c2", "M.c3",
                                 "M.c4", "M.c5", "M.c6", "M.c7"};
const mvec_coeff mv3d_coeff_M1 = {"M1.c0", "M1.c1", "M1.c2", "M1.c3",
                                  "M1.c4", "M1.c5", "M1.c6", "M1.c7"};
const mvec_coeff mv3d_coeff_M2 = {"M2.c0", "M2.c1", "M2.c2", "M2.c3",
                                  "M2.c4", "M2.c5", "M2.c6", "M2.c7"};
const mvec_coeff mv3d_coeff_A_even = {"A.c0", "0", "0", "0", "A.c1", "A.c2", "A.c3", "0"};
const mvec_coeff mv3d_coeff_B_even = {"B.c0", "0", "0", "0", "B.c1", "B.c2", "B.c3", "0"};
const mvec_coeff mv3d_coeff_A_odd = {"0", "A.c0", "A.c1", "A.c2", "0", "0", "0", "A.c3"};
const mvec_coeff mv3d_coeff_B_odd = {"0", "B.c0", "B.c1", "B.c2", "0", "0", "0", "B.c3"};
const mvec_coeff mv3d_coeff_M_even = {"M.c0", "0", "0", "0", "M.c1", "M.c2", "M.c3", "0"};
const mvec_coeff mv3d_coeff_M_odd = {"0", "M.c0", "M.c1", "M.c2", "0", "0", "0", "M.c3"};
const mvec_coeff mv3d_coeff_R_even = {"R.c0", "0", "0", "0", "R.c1", "R.c2", "R.c3", "0"};
const mvec_coeff mv3d_coeff_R_rev_even = {"R.c0",  "0",     "0",     "0",
                                          "-R.c1", "-R.c2", "-R.c3", "0"};
const mvec_coeff mv3d_coeff_svBps = {"s", "v.x", "v.y", "v.z", "B.x", "B.y", "B.z", "ps"};
const mvec_coeff mv3d_coeff_svBps1 = {"s1",   "v1.x", "v1.y", "v1.z",
                                      "B1.x", "B1.y", "B1.z", "ps1"};
const mvec_coeff mv3d_coeff_svBps2 = {"s2",   "v2.x", "v2.y", "v2.z",
                                      "B2.x", "B2.y", "B2.z", "ps2"};
