#pragma once

// author: Daniel Hug, 2024

#include "ga_prdxpr_common.hpp"

// multivector basis components
const mvec_coeff mv3d_basis = {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e123"};

// multivector basis products
const prd_table gpr_ega3d_basis_tab = {
    {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e123"},
    {"e1", "1", "e12", "-e31", "e123", "-e3", "e2", "e23"},
    {"e2", "-e12", "1", "e23", "e3", "e123", "-e1", "e31"},
    {"e3", "e31", "-e23", "1", "-e2", "e1", "e123", "e12"},
    {"e23", "e123", "-e3", "e2", "-1", "-e12", "e31", "-e1"},
    {"e31", "e3", "e123", "-e1", "e12", "-1", "-e23", "-e2"},
    {"e12", "-e2", "e1", "e123", "-e31", "e23", "-1", "-e3"},
    {"e123", "e23", "e31", "e12", "-e1", "-e2", "-e3", "-1"}};

const prd_table wdg_ega3d_basis_tab = {
    {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e123"},
    {"e1", "0", "e12", "-e31", "e123", "0", "0", "0"},
    {"e2", "-e12", "0", "e23", "0", "e123", "0", "0"},
    {"e3", "e31", "-e23", "0", "0", "0", "e123", "0"},
    {"e23", "e123", "0", "0", "0", "0", "0", "0"},
    {"e31", "0", "e123", "0", "0", "0", "0", "0"},
    {"e12", "0", "0", "e123", "0", "0", "0", "0"},
    {"e123", "0", "0", "0", "0", "0", "0", "0"}};

const prd_table lcontr_ega3d_basis_tab = {
    {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e123"},
    {"0", "1", "0", "0", "0", "-e3", "e2", "e23"},
    {"0", "0", "1", "0", "e3", "0", "-e1", "e31"},
    {"0", "0", "0", "1", "-e2", "e1", "0", "e12"},
    {"0", "0", "0", "0", "-1", "0", "0", "-e1"},
    {"0", "0", "0", "0", "0", "-1", "0", "-e2"},
    {"0", "0", "0", "0", "0", "0", "-1", "-e3"},
    {"0", "0", "0", "0", "0", "0", "0", "-1"}};

const prd_table rcontr_ega3d_basis_tab = {
    {"1", "0", "0", "0", "0", "0", "0", "0"},
    {"e1", "1", "0", "0", "0", "0", "0", "0"},
    {"e2", "0", "1", "0", "0", "0", "0", "0"},
    {"e3", "0", "0", "1", "0", "0", "0", "0"},
    {"e23", "0", "-e3", "e2", "-1", "0", "0", "0"},
    {"e31", "e3", "0", "-e1", "0", "-1", "0", "0"},
    {"e12", "-e2", "e1", "0", "0", "0", "-1", "0"},
    {"e123", "e23", "e31", "e12", "-e1", "-e2", "-e3", "-1"}};

const prd_table dot_ega3d_basis_tab = {
    {"1", "0", "0", "0", "0", "0", "0", "0"},  {"0", "1", "0", "0", "0", "0", "0", "0"},
    {"0", "0", "1", "0", "0", "0", "0", "0"},  {"0", "0", "0", "1", "0", "0", "0", "0"},
    {"0", "0", "0", "0", "-1", "0", "0", "0"}, {"0", "0", "0", "0", "0", "-1", "0", "0"},
    {"0", "0", "0", "0", "0", "0", "-1", "0"}, {"0", "0", "0", "0", "0", "0", "0", "-1"}};

// coefficients need to create a multivector = [coeff]^T [basis]

// const mvec_coeff mv_ega3d_lcoeff = {"sl",  "vxl", "vyl", "vzl",
//                                     "Bxl", "Byl", "Bzl", "psl"};
// const mvec_coeff mv_ega3d_rcoeff = {"sr",  "vxr", "vyr", "vzr",
//                                     "Bxr", "Byr", "Bzr", "psr"};

const mvec_coeff mv3d_coeff_A = {"A.c0", "A.c1", "A.c2", "A.c3",
                                 "A.c4", "A.c5", "A.c6", "A.c7"};
const mvec_coeff mv3d_coeff_B = {"B.c0", "B.c1", "B.c2", "B.c3",
                                 "B.c4", "B.c5", "B.c6", "B.c7"};
const mvec_coeff mv3d_coeff_M1 = {"M1.c0", "M1.c1", "M1.c2", "M1.c3",
                                  "M1.c4", "M1.c5", "M1.c6", "M1.c7"};
const mvec_coeff mv3d_coeff_M2 = {"M2.c0", "M2.c1", "M2.c2", "M2.c3",
                                  "M2.c4", "M2.c5", "M2.c6", "M2.c7"};
const mvec_coeff mv3d_coeff_A_even = {"A.c0", "0", "0", "0", "A.c1", "A.c2", "A.c3", "0"};
const mvec_coeff mv3d_coeff_B_even = {"B.c0", "0", "0", "0", "B.c1", "B.c2", "B.c3", "0"};
const mvec_coeff mv3d_coeff_A_uneven = {"0", "A.c0", "A.c1", "A.c2",
                                        "0", "0",    "0",    "A.c3"};
const mvec_coeff mv3d_coeff_B_uneven = {"0", "B.c0", "B.c1", "B.c2",
                                        "0", "0",    "0",    "B.c3"};
const mvec_coeff mv3d_coeff_svBps = {"s", "v.x", "v.y", "v.z", "B.x", "B.y", "B.z", "ps"};
const mvec_coeff mv3d_coeff_svBps1 = {"s1",   "v1.x", "v1.y", "v1.z",
                                      "B1.x", "B1.y", "B1.z", "ps1"};
const mvec_coeff mv3d_coeff_svBps2 = {"s2",   "v2.x", "v2.y", "v2.z",
                                      "B2.x", "B2.y", "B2.z", "ps2"};


void generate_and_print_ega3d_gpr();    // geometric product
void generate_and_print_ega3d_wdg();    // wedge product
void generate_and_print_ega3d_lcontr(); // left contraction
void generate_and_print_ega3d_rcontr(); // right contraction
void generate_and_print_ega3d_dot();    // scalar product