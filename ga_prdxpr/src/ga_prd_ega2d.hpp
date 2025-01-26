#pragma once

// author: Daniel Hug, 2024

#include "ga_prdxpr_common.hpp"

// multivector basis components
const mvec_coeff mv2d_basis = {"1", "e1", "e2", "e12"};

// multivector basis products
const prd_table gpr_ega2d_basis_tab = {{"1", "e1", "e2", "e12"},
                                       {"e1", "1", "e12", "e2"},
                                       {"e2", "-e12", "1", "-e1"},
                                       {"e12", "-e2", "e1", "-1"}};

const prd_table wdg_ega2d_basis_tab = {{"1", "e1", "e2", "e12"},
                                       {"e1", "0", "e12", "0"},
                                       {"e2", "-e12", "0", "0"},
                                       {"e12", "0", "0", "0"}};

const prd_table lcontr_ega2d_basis_tab = {{"1", "e1", "e2", "e12"},
                                          {"0", "1", "0", "e2"},
                                          {"0", "0", "1", "-e1"},
                                          {"0", "0", "0", "-1"}};

const prd_table rcontr_ega2d_basis_tab = {{"1", "0", "0", "0"},
                                          {"e1", "1", "0", "0"},
                                          {"e2", "0", "1", "0"},
                                          {"e12", "-e2", "e1", "-1"}};

const prd_table dot_ega2d_basis_tab = {{"1", "0", "0", "0"},
                                       {"0", "1", "0", "0"},
                                       {"0", "0", "1", "0"},
                                       {"0", "0", "0", "-1"}};

// coefficients need to create a multivector = [coeff]^T [basis]

// const mvec_coeff mv2d_lcoeff = {"sl", "vxl", "vyl", "psl"};
// const mvec_coeff mv2d_rcoeff = {"sr", "vxr", "vyr", "psr"};

const mvec_coeff mv2d_coeff_A = {"A.c0", "A.c1", "A.c2", "A.c3"};
const mvec_coeff mv2d_coeff_B = {"B.c0", "B.c1", "B.c2", "B.c3"};
const mvec_coeff mv2d_coeff_M1 = {"M1.c0", "M1.c1", "M1.c2", "M1.c3"};
const mvec_coeff mv2d_coeff_M2 = {"M2.c0", "M2.c1", "M2.c2", "M2.c3"};
const mvec_coeff mv2d_coeff_A_even = {"A.c0", "0", "0", "A.c1"};
const mvec_coeff mv2d_coeff_B_even = {"B.c0", "0", "0", "B.c1"};
const mvec_coeff mv2d_coeff_svps = {"s", "v.x", "v.y", "ps"};
const mvec_coeff mv2d_coeff_svps1 = {"s1", "v1.x", "v1.y", "ps1"};
const mvec_coeff mv2d_coeff_svps2 = {"s2", "v2.x", "v2.y", "ps2"};


void generate_and_print_ega2d_gpr();    // geometric product
void generate_and_print_ega2d_wdg();    // wedge product
void generate_and_print_ega2d_lcontr(); // left contraction
void generate_and_print_ega2d_rcontr(); // right contraction
void generate_and_print_ega2d_dot();    // scalar product