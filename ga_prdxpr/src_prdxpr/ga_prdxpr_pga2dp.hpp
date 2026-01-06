#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"
#include "ga_prdxpr_rule_generator.hpp"

// pga2dp: projective geometric algebra with
//         e1*e1 = +1, e2*e2 = +1, e3*e3 = 0
// G(2,0,1)

// metric signature: e1²=+1, e2²=+1, e3²=0
const std::vector<int> mv2dp_metric_signature = {+1, +1, 0};

// multivector basis components
const mvec_coeff mv2dp_basis = {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e321"};
// const mvec_coeff mv2dp_basis = {"1", "e1", "e2", "e3", "e31", "e32", "e12", "e321"};
// const mvec_coeff mv2dp_basis = {"1", "e1", "e2", "e3", "e31", "e32", "e12", "e123"};

// multivector basis components as vector for each order k
// order k=0: scalar
// order k=1: vectors
// order k=2: bivectors
// order k=3: trivector (=pseudoscalar for pga2dp)
const std::vector<mvec_coeff> mv2dp_basis_kvec = {
    {"1"}, {"e1", "e2", "e3"}, {"e23", "e31", "e12"}, {"e321"}};
// const std::vector<mvec_coeff> mv2dp_basis_kvec = {
//     {"1"}, {"e1", "e2", "e3"}, {"e31", "e32", "e12"}, {"e321"}};
// const std::vector<mvec_coeff> mv2dp_basis_kvec = {
//     {"1"}, {"e1", "e2", "e3"}, {"e31", "e32", "e12"}, {"e123"}};

////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - PGA2DP
//
// Product tables are generated automatically using the rule generation system:
// - Configuration: ga_prdxpr_pga2dp_config.cpp (defines algebra parameters)
// - Generation: ga_prdxpr_pga2dp.cpp (creates product rules at startup)
// - Interface: This file (declares available rules)
//
// This ensures mathematically correct, consistent rules without manual errors.
////////////////////////////////////////////////////////////////////////////////
extern const prd_rules gpr_pga2dp_rules; // Geometric product: a * b
extern const prd_rules wdg_pga2dp_rules; // Wedge product: a ^ b
extern const prd_rules dot_pga2dp_rules; // Dot product: a · b

// complement in 2dp: cmpl(cmpl(u)) = u
// lcmpl(u) ^ u = u ^ rcmpl(u) = cmpl(u) ^ u = u ^ cmpl(u) = e321
extern const prd_rules cmpl_pga2dp_rules;

// dual rules: generated automatically from complement rules + extended metric
// PGA2DP (odd-dimensional) has bulk_dual and weight_dual
// bulk_dual(u) = complement(G × u) where G is the extended metric matrix
// weight_dual(u) = complement(G × complement(u)) uses the regressive metric
// Note: bulk_dual is equivalent to the Hodge dual for PGA2DP
extern const prd_rules bulk_dual_pga2dp_rules;
extern const prd_rules weight_dual_pga2dp_rules;

// coefficients needed to create a multivector = [coeff]^T [mv2dp_basis]
const mvec_coeff mv2dp_coeff_A = {"A.c0", "A.c1", "A.c2", "A.c3",
                                  "A.c4", "A.c5", "A.c6", "A.c7"};
const mvec_coeff mv2dp_coeff_B = {"B.c0", "B.c1", "B.c2", "B.c3",
                                  "B.c4", "B.c5", "B.c6", "B.c7"};
const mvec_coeff mv2dp_coeff_M = {"M.c0", "M.c1", "M.c2", "M.c3",
                                  "M.c4", "M.c5", "M.c6", "M.c7"};
const mvec_coeff mv2dp_coeff_M1 = {"M1.c0", "M1.c1", "M1.c2", "M1.c3",
                                   "M1.c4", "M1.c5", "M1.c6", "M1.c7"};
const mvec_coeff mv2dp_coeff_M2 = {"M2.c0", "M2.c1", "M2.c2", "M2.c3",
                                   "M2.c4", "M2.c5", "M2.c6", "M2.c7"};
const mvec_coeff mv2dp_coeff_M_even = {"M.c0", "0",    "0",    "0",
                                       "M.c1", "M.c2", "M.c3", "0"};
const mvec_coeff mv2dp_coeff_M_odd = {"0", "M.c0", "M.c1", "M.c2", "0", "0", "0", "M.c3"};
const mvec_coeff mv2dp_coeff_A_even = {"A.c0", "0",    "0",    "0",
                                       "A.c1", "A.c2", "A.c3", "0"};
const mvec_coeff mv2dp_coeff_B_even = {"B.c0", "0",    "0",    "0",
                                       "B.c1", "B.c2", "B.c3", "0"};
const mvec_coeff mv2dp_coeff_A_odd = {"0", "A.c0", "A.c1", "A.c2", "0", "0", "0", "A.c3"};
const mvec_coeff mv2dp_coeff_B_odd = {"0", "B.c0", "B.c1", "B.c2", "0", "0", "0", "B.c3"};
const mvec_coeff mv2dp_coeff_R_even = {"R.c0", "0",    "0",    "0",
                                       "R.c1", "R.c2", "R.c3", "0"};
const mvec_coeff mv2dp_coeff_R_odd = {"0", "R.c0", "R.c1", "R.c2", "0", "0", "0", "R.c3"};
const mvec_coeff mv2dp_coeff_R_rev_even = {"R.c0",  "0",     "0",     "0",
                                           "-R.c1", "-R.c2", "-R.c3", "0"};
const mvec_coeff mv2dp_coeff_R_rev_odd = {"0", "R.c0", "R.c1", "R.c2",
                                          "0", "0",    "0",    "-R.c3"};
const mvec_coeff mv2dp_coeff_R_rrev_even = {"-R.c0", "0",    "0",    "0",
                                            "R.c1",  "R.c2", "R.c3", "0"};
const mvec_coeff mv2dp_coeff_R_rrev_odd = {"0", "-R.c0", "-R.c1", "-R.c2",
                                           "0", "0",     "0",     "R.c3"};
const mvec_coeff mv2dp_coeff_svBps = {"s",   "v.x", "v.y", "v.z",
                                      "B.x", "B.y", "B.z", "ps"};
const mvec_coeff mv2dp_coeff_svBps1 = {"s1",   "v1.x", "v1.y", "v1.z",
                                       "B1.x", "B1.y", "B1.z", "ps1"};
const mvec_coeff mv2dp_coeff_svBps2 = {"s2",   "v2.x", "v2.y", "v2.z",
                                       "B2.x", "B2.y", "B2.z", "ps2"};
