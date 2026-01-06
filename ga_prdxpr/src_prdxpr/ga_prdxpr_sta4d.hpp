#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"
#include "ga_prdxpr_rule_generator.hpp"

// sta4d: space-time geometric algebra with (g is short for gamma)
//        g0*g0 = +1, g1*g1 = -1, g2*g2 = -1, g3*g3 = -1
//
//        g0         ... time-like dimension
//        g1, g2, g3 ... space-like dimensions
//
// G(1,3,0)

// metric signature: g0²=+1, g1²=-1, g2²=-1, g3²=-1
const std::vector<int> mvsta4d_metric_signature = {+1, -1, -1, -1};

// multivector basis components
const mvec_coeff mvsta4d_basis = {"1",    "g0",   "g1",   "g2",   "g3",  "g01",
                                  "g02",  "g03",  "g23",  "g31",  "g12", "g023",
                                  "g031", "g012", "g123", "g0123"};

// multivector basis components as vector for each order k
// order k=0: scalar
// order k=1: vectors
// order k=2: bivectors
// order k=3: trivectors
// order k=4: quadvector (=pseudoscalar for sta4d)
const std::vector<mvec_coeff> mvsta4d_basis_kvec = {
    {"1"},
    {"g0", "g1", "g2", "g3"},
    {"g01", "g02", "g03", "g23", "g31", "g12"},
    {"g023", "g031", "g012", "g123"},
    {"g0123"}};
////////////////////////////////////////////////////////////////////////////////
// AUTOMATIC RULE GENERATION - sta4d
//
// Product tables are generated automatically using the rule generation system:
// - Configuration: ga_prdxpr_sta4d_config.cpp (defines algebra parameters)
// - Generation: ga_prdxpr_sta4d.cpp (creates product rules at startup)
// - Interface: This file (declares available rules)
//
// This ensures mathematically correct, consistent rules without manual errors.
////////////////////////////////////////////////////////////////////////////////
extern const prd_rules gpr_sta4d_rules; // Geometric product: a * b
extern const prd_rules wdg_sta4d_rules; // Wedge product: a ^ b
extern const prd_rules dot_sta4d_rules; // Dot product: a · b

// complement in sta4d: lcmpl(rcmpl(u)) = u
// lcmpl:  lcmpl(u) ^ u        = g0123
// rcmpl:  u        ^ rcmpl(u) = g0123
extern const prd_rules lcmpl_sta4d_rules;
extern const prd_rules rcmpl_sta4d_rules;

// dual rules: generated automatically from complement rules + extended metric
// dual(u) = complement(G × u) where G is the extended metric matrix
// left_dual and right_dual are derived from left_complement and right_complement
// respectively Note: right_dual is equivalent to the Hodge dual for STA4D
extern const prd_rules left_dual_sta4d_rules;
extern const prd_rules right_dual_sta4d_rules;

// coefficients needed to create a multivector = [coeff]^T [mvsta4d_basis]
const mvec_coeff mvsta4d_coeff_A = {"A.c0",  "A.c1",  "A.c2",  "A.c3", "A.c4",  "A.c5",
                                    "A.c6",  "A.c7",  "A.c8",  "A.c9", "A.c10", "A.c11",
                                    "A.c12", "A.c13", "A.c14", "A.c15"};
const mvec_coeff mvsta4d_coeff_B = {"B.c0",  "B.c1",  "B.c2",  "B.c3", "B.c4",  "B.c5",
                                    "B.c6",  "B.c7",  "B.c8",  "B.c9", "B.c10", "B.c11",
                                    "B.c12", "B.c13", "B.c14", "B.c15"};
const mvec_coeff mvsta4d_coeff_M = {"M.c0",  "M.c1",  "M.c2",  "M.c3", "M.c4",  "M.c5",
                                    "M.c6",  "M.c7",  "M.c8",  "M.c9", "M.c10", "M.c11",
                                    "M.c12", "M.c13", "M.c14", "M.c15"};
const mvec_coeff mvsta4d_coeff_M1 = {
    "M1.c0", "M1.c1", "M1.c2",  "M1.c3",  "M1.c4",  "M1.c5",  "M1.c6",  "M1.c7",
    "M1.c8", "M1.c9", "M1.c10", "M1.c11", "M1.c12", "M1.c13", "M1.c14", "M1.c15"};
const mvec_coeff mvsta4d_coeff_M2 = {
    "M2.c0", "M2.c1", "M2.c2",  "M2.c3",  "M2.c4",  "M2.c5",  "M2.c6",  "M2.c7",
    "M2.c8", "M2.c9", "M2.c10", "M2.c11", "M2.c12", "M2.c13", "M2.c14", "M2.c15"};
const mvec_coeff mvsta4d_coeff_M_even = {"M.c0", "0",    "0",    "0",    "0",    "M.c1",
                                         "M.c2", "M.c3", "M.c4", "M.c5", "M.c6", "0",
                                         "0",    "0",    "0",    "M.c7"};
const mvec_coeff mvsta4d_coeff_M_odd = {"0",    "M.c0", "M.c1", "M.c2", "M.c3", "0",
                                        "0",    "0",    "0",    "0",    "0",    "M.c4",
                                        "M.c5", "M.c6", "M.c7", "0"};
const mvec_coeff mvsta4d_coeff_A_even = {"A.c0", "0",    "0",    "0",    "0",    "A.c1",
                                         "A.c2", "A.c3", "A.c4", "A.c5", "A.c6", "0",
                                         "0",    "0",    "0",    "A.c7"};
const mvec_coeff mvsta4d_coeff_B_even = {"B.c0", "0",    "0",    "0",    "0",    "B.c1",
                                         "B.c2", "B.c3", "B.c4", "B.c5", "B.c6", "0",
                                         "0",    "0",    "0",    "B.c7"};
const mvec_coeff mvsta4d_coeff_R_even = {"R.c0", "0",    "0",    "0",    "0",    "R.c1",
                                         "R.c2", "R.c3", "R.c4", "R.c5", "R.c6", "0",
                                         "0",    "0",    "0",    "R.c7"};
const mvec_coeff mvsta4d_coeff_R_rev_even = {
    "R.c0",  "0",     "0",     "0", "0", "-R.c1", "-R.c2", "-R.c3",
    "-R.c4", "-R.c5", "-R.c6", "0", "0", "0",     "0",     "R.c7"};
const mvec_coeff mvsta4d_coeff_R_rrev_even = {
    "R.c0",  "0",     "0",     "0", "0", "-R.c1", "-R.c2", "-R.c3",
    "-R.c4", "-R.c5", "-R.c6", "0", "0", "0",     "0",     "R.c7"};
const mvec_coeff mvsta4d_coeff_A_odd = {"0",    "A.c0", "A.c1", "A.c2", "A.c3", "0",
                                        "0",    "0",    "0",    "0",    "0",    "A.c4",
                                        "A.c5", "A.c6", "A.c7", "0"};
const mvec_coeff mvsta4d_coeff_B_odd = {"0",    "B.c0", "B.c1", "B.c2", "B.c3", "0",
                                        "0",    "0",    "0",    "0",    "0",    "B.c4",
                                        "B.c5", "B.c6", "B.c7", "0"};
const mvec_coeff mvsta4d_coeff_R_odd = {"0",    "R.c0", "R.c1", "R.c2", "R.c3", "0",
                                        "0",    "0",    "0",    "0",    "0",    "R.c4",
                                        "R.c5", "R.c6", "R.c7", "0"};
const mvec_coeff mvsta4d_coeff_R_rev_odd = {
    "0", "R.c0", "R.c1", "R.c2",  "R.c3",  "0",     "0",     "0",
    "0", "0",    "0",    "-R.c4", "-R.c5", "-R.c6", "-R.c7", "0"};
const mvec_coeff mvsta4d_coeff_R_rrev_odd = {
    "0", "-R.c0", "-R.c1", "-R.c2", "-R.c3", "0",    "0",    "0",
    "0", "0",     "0",     "R.c4",  "R.c5",  "R.c6", "R.c7", "0"};
const mvec_coeff mvsta4d_coeff_svBtps = {"s",    "v.x",  "v.y",  "v.z",  "v.w",  "B.vx",
                                         "B.vy", "B.vz", "B.mx", "B.my", "B.mz", "t.x",
                                         "t.y",  "t.z",  "t.w",  "ps"};
const mvec_coeff mvsta4d_coeff_svBtps1 = {
    "s1",    "v1.x",  "v1.y",  "v1.z", "v1.w", "B1.vx", "B1.vy", "B1.vz",
    "B1.mx", "B1.my", "B1.mz", "t1.x", "t1.y", "t1.z",  "t1.w",  "ps1"};
const mvec_coeff mvsta4d_coeff_svBtps2 = {
    "s2",    "v2.x",  "v2.y",  "v2.z", "v2.w", "B2.vx", "B2.vy", "B2.vz",
    "B2.mx", "B2.my", "B2.mz", "t2.x", "t2.y", "t2.z",  "t2.w",  "ps2"};
