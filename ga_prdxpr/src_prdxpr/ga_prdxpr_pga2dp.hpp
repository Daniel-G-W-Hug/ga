#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"

// pga2dp: projective geometric algebra with
//         e1*e1 = +1, e2*e2 = +1, e3*e3 = 0
// G(2,0,1)

// multivector basis components
const mvec_coeff mv2dp_basis = {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e321"};

// rules to create produc tables (must contain space_str around operator_str on each side)
const prd_rules gpr_pga2dp_rules = {
    {"1 * 1", "1"},        {"1 * e1", "e1"},      {"1 * e2", "e2"},
    {"1 * e3", "e3"},      {"1 * e23", "e23"},    {"1 * e31", "e31"},
    {"1 * e12", "e12"},    {"1 * e321", "e321"},

    {"e1 * 1", "e1"},      {"e1 * e1", "1"},      {"e1 * e2", "e12"},
    {"e1 * e3", "-e31"},   {"e1 * e23", "-e321"}, {"e1 * e31", "-e3"},
    {"e1 * e12", "e2"},    {"e1 * e321", "-e23"},

    {"e2 * 1", "e2"},      {"e2 * e1", "-e12"},   {"e2 * e2", "1"},
    {"e2 * e3", "e23"},    {"e2 * e23", "e3"},    {"e2 * e31", "-e321"},
    {"e2 * e12", "-e1"},   {"e2 * e321", "-e31"},

    {"e3 * 1", "e3"},      {"e3 * e1", "e31"},    {"e3 * e2", "-e23"},
    {"e3 * e3", "0"},      {"e3 * e23", "0"},     {"e3 * e31", "0"},
    {"e3 * e12", "-e321"}, {"e3 * e321", "0"},

    {"e23 * 1", "e23"},    {"e23 * e1", "-e321"}, {"e23 * e2", "-e3"},
    {"e23 * e3", "0"},     {"e23 * e23", "0"},    {"e23 * e31", "0"},
    {"e23 * e12", "e31"},  {"e23 * e321", "0"},

    {"e31 * 1", "e31"},    {"e31 * e1", "e3"},    {"e31 * e2", "-e321"},
    {"e31 * e3", "0"},     {"e31 * e23", "0"},    {"e31 * e31", "0"},
    {"e31 * e12", "-e23"}, {"e31 * e321", "0"},

    {"e12 * 1", "e12"},    {"e12 * e1", "-e2"},   {"e12 * e2", "e1"},
    {"e12 * e3", "-e321"}, {"e12 * e23", "-e31"}, {"e12 * e31", "e23"},
    {"e12 * e12", "-1"},   {"e12 * e321", "e3"},

    {"e321 * 1", "e321"},  {"e321 * e1", "-e23"}, {"e321 * e2", "-e31"},
    {"e321 * e3", "0"},    {"e321 * e23", "0"},   {"e321 * e31", "0"},
    {"e321 * e12", "e3"},  {"e321 * e321", "0"}};

const prd_rules wdg_pga2dp_rules = {
    {"1 ^ 1", "1"},        {"1 ^ e1", "e1"},      {"1 ^ e2", "e2"},
    {"1 ^ e3", "e3"},      {"1 ^ e23", "e23"},    {"1 ^ e31", "e31"},
    {"1 ^ e12", "e12"},    {"1 ^ e321", "e321"},

    {"e1 ^ 1", "e1"},      {"e1 ^ e1", "0"},      {"e1 ^ e2", "e12"},
    {"e1 ^ e3", "-e31"},   {"e1 ^ e23", "-e321"}, {"e1 ^ e31", "0"},
    {"e1 ^ e12", "0"},     {"e1 ^ e321", "0"},

    {"e2 ^ 1", "e2"},      {"e2 ^ e1", "-e12"},   {"e2 ^ e2", "0"},
    {"e2 ^ e3", "e23"},    {"e2 ^ e23", "0"},     {"e2 ^ e31", "-e321"},
    {"e2 ^ e12", "0"},     {"e2 ^ e321", "0"},

    {"e3 ^ 1", "e3"},      {"e3 ^ e1", "e31"},    {"e3 ^ e2", "-e23"},
    {"e3 ^ e3", "0"},      {"e3 ^ e23", "0"},     {"e3 ^ e31", "0"},
    {"e3 ^ e12", "-e321"}, {"e3 ^ e321", "0"},

    {"e23 ^ 1", "e23"},    {"e23 ^ e1", "-e321"}, {"e23 ^ e2", "0"},
    {"e23 ^ e3", "0"},     {"e23 ^ e23", "0"},    {"e23 ^ e31", "0"},
    {"e23 ^ e12", "0"},    {"e23 ^ e321", "0"},

    {"e31 ^ 1", "e31"},    {"e31 ^ e1", "0"},     {"e31 ^ e2", "-e321"},
    {"e31 ^ e3", "0"},     {"e31 ^ e23", "0"},    {"e31 ^ e31", "0"},
    {"e31 ^ e12", "0"},    {"e31 ^ e321", "0"},

    {"e12 ^ 1", "e12"},    {"e12 ^ e1", "0"},     {"e12 ^ e2", "0"},
    {"e12 ^ e3", "-e321"}, {"e12 ^ e23", "0"},    {"e12 ^ e31", "0"},
    {"e12 ^ e12", "0"},    {"e12 ^ e321", "0"},

    {"e321 ^ 1", "e321"},  {"e321 ^ e1", "0"},    {"e321 ^ e2", "0"},
    {"e321 ^ e3", "0"},    {"e321 ^ e23", "0"},   {"e321 ^ e31", "0"},
    {"e321 ^ e12", "0"},   {"e321 ^ e321", "0"}};

const prd_rules dot_pga2dp_rules = {
    {"1 * 1", "1"},      {"1 * e1", "0"},     {"1 * e2", "0"},     {"1 * e3", "0"},
    {"1 * e23", "0"},    {"1 * e31", "0"},    {"1 * e12", "0"},    {"1 * e321", "0"},

    {"e1 * 1", "0"},     {"e1 * e1", "1"},    {"e1 * e2", "0"},    {"e1 * e3", "0"},
    {"e1 * e23", "0"},   {"e1 * e31", "0"},   {"e1 * e12", "0"},   {"e1 * e321", "0"},

    {"e2 * 1", "0"},     {"e2 * e1", "0"},    {"e2 * e2", "1"},    {"e2 * e3", "0"},
    {"e2 * e23", "0"},   {"e2 * e31", "0"},   {"e2 * e12", "0"},   {"e2 * e321", "0"},

    {"e3 * 1", "0"},     {"e3 * e1", "0"},    {"e3 * e2", "0"},    {"e3 * e3", "0"},
    {"e3 * e23", "0"},   {"e3 * e31", "0"},   {"e3 * e12", "0"},   {"e3 * e321", "0"},

    {"e23 * 1", "0"},    {"e23 * e1", "0"},   {"e23 * e2", "0"},   {"e23 * e3", "0"},
    {"e23 * e23", "0"},  {"e23 * e31", "0"},  {"e23 * e12", "0"},  {"e23 * e321", "0"},

    {"e31 * 1", "0"},    {"e31 * e1", "0"},   {"e31 * e2", "0"},   {"e31 * e3", "0"},
    {"e31 * e23", "0"},  {"e31 * e31", "0"},  {"e31 * e12", "0"},  {"e31 * e321", "0"},

    {"e12 * 1", "0"},    {"e12 * e1", "0"},   {"e12 * e2", "0"},   {"e12 * e3", "0"},
    {"e12 * e23", "0"},  {"e12 * e31", "0"},  {"e12 * e12", "1"},  {"e12 * e321", "0"},

    {"e321 * 1", "0"},   {"e321 * e1", "0"},  {"e321 * e2", "0"},  {"e321 * e3", "0"},
    {"e321 * e23", "0"}, {"e321 * e31", "0"}, {"e321 * e12", "0"}, {"e321 * e321", "0"}};


// complement in 2dp: cmpl(cmpl(u)) = u
// lcmpl(u) ^ u = u ^ rcmpl(u) = cmpl(u) ^ u = u ^ cmpl(u) = e321
const prd_rules cmpl_pga2dp_rules = {{"1", "e321"},  {"e1", "-e23"}, {"e2", "-e31"},
                                     {"e3", "-e12"}, {"e23", "-e1"}, {"e31", "-e2"},
                                     {"e12", "-e3"}, {"e321", "1"}};

// rules needed to derive contractions and expansions using bulk_dual and weight_dual
// same for left hand and right hand side
// (derived from the complement -> bulk dual is equivalent to the Hodge dual)
const prd_rules bulk_dual_pga2dp_rules = {{"1", "e321"},  {"e1", "-e23"}, {"e2", "-e31"},
                                          {"e3", "0"},    {"e23", "0"},   {"e31", "0"},
                                          {"e12", "-e3"}, {"e321", "0"}};
const prd_rules weight_dual_pga2dp_rules = {
    {"1", "0"},     {"e1", "0"},    {"e2", "0"},  {"e3", "-e12"},
    {"e23", "-e1"}, {"e31", "-e2"}, {"e12", "0"}, {"e321", "1"}};

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
