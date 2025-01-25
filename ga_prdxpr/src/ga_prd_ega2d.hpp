#pragma once

// author: Daniel Hug, 2024

#include "ga_prdxpr_common.hpp"

const mvec_coeff mv_ega2d_basis = {"1", "e1", "e2", "e12"};
const prd_table gpr_ega2d_basis_tab = {{"1", "e1", "e2", "e12"},
                                       {"e1", "1", "e12", "e2"},
                                       {"e2", "-e12", "1", "-e1"},
                                       {"e12", "-e2", "e1", "-1"}};

// const mvec_coeff mv_ega2d_basis = {"1", "e_1", "e_2", "e_{12}"};
// const prd_table gpr_ega2d_basis_tab = {{"1", "e_1", "e_2", "e_{12}"},
//                                        {"e_1", "1", "e_{12}", "e_2"},
//                                        {"e_2", "-e_{12}", "1", "-e_1"},
//                                        {"e_{12}", "-e_2", "e_1", "-1"}};

const prd_table wdg_ega2d_basis_tab = {{"1", "e1", "e2", "e12"},
                                       {"e1", "0", "e12", "0"},
                                       {"e2", "-e12", "0", "0"},
                                       {"e12", "0", "0", "0"}};

const prd_table dot_ega2d_basis_tab = {{"1", "0", "0", "0"},
                                       {"0", "1", "0", "0"},
                                       {"0", "0", "1", "0"},
                                       {"0", "0", "0", "-1"}};

const prd_table lcontr_ega2d_basis_tab = {{"1", "e1", "e2", "e12"},
                                          {"0", "1", "0", "e2"},
                                          {"0", "0", "1", "-e1"},
                                          {"0", "0", "0", "-1"}};

const prd_table rcontr_ega2d_basis_tab = {{"1", "0", "0", "0"},
                                          {"e1", "1", "0", "0"},
                                          {"e2", "0", "1", "0"},
                                          {"e12", "-e2", "e1", "-1"}};

// coefficients need to create a multivector = [coeff]^T [basis]

// const mvec_coeff mv_ega2d_lcoeff = {"sl", "vxl", "vyl", "psl"};
// const mvec_coeff mv_ega2d_rcoeff = {"sr", "vxr", "vyr", "psr"};

// const mvec_coeff mv_ega2d_lcoeff = {"^l", "v_x^l", "v_x^l", "ps^l"};
// const mvec_coeff mv_ega2d_rcoeff = {"^r", "v_x^r", "v_x^r", "ps^r"};

const mvec_coeff mv_ega2d_lcoeff = {"A.c0", "A.c1", "A.c2", "A.c3"};
const mvec_coeff mv_ega2d_rcoeff = {"B.c0", "B.c1", "B.c2", "B.c3"};
