// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_ega3d.hpp"
#include "ga_prdxpr_common.hpp"


void generate_and_print_ega3d_gpr()
{

    // ega3d geometric product
    std::string prd_name = "ega3d geometric product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str), gpr_ega3d_rules);
    auto basis = mv3d_basis;

    // fmt::println("equal = {}", basis_tab == gpr_ega3d_basis_tab);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv * mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * mv_u -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_odd);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_u * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * ps -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * bivec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * vec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * s -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_odd);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_odd);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_u * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * ps -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * ps -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e * s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * s -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * ps -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * ps -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * ps -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s * s -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega3d_cmt()
{

    // the commutator product is the asymmetric part of the geometric product
    std::string prd_name = "ega3d commutator product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str), gpr_ega3d_rules);
    auto basis = mv3d_basis;

    basis_tab = get_prd_tab_asym(basis_tab); // use the asymmetric part only

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "cmt(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(ps,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(ps,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(ps,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(vec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(vec,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(ps,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(ps,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(ps,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "cmt(bivec,bivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "cmt(bivec,vec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "cmt(vec,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "cmt(vec,vec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega3d_wdg()
{

    // ega3d wedge product
    std::string prd_name = "ega3d wedge product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, wdg_str), wdg_ega3d_rules);
    auto basis = mv3d_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv ^ mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ bivec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ vec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ s -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv_e -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv_u -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ bivec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv_e -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ bivec -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv_u -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ vec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ s -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ vec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ bivec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s ^ s -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_ega3d_lcontract()
{

    std::string prd_name = "ega3d left contraction";

    // (A << B) = rwdg( ldual(A),B ) -> in ega3d = rwdg(cmpl(A),B)
    //
    // with rwdg(A,B) = cmpl( wdg( cmpl(A),cmpl(B) ) ) in ega3d

    // create the complement from the input multivector (=forward transformation)
    auto basis = mv3d_basis;

    // lhs (dualize as inner function before complement)
    auto lhs =
        apply_rules_to_mv(apply_rules_to_mv(basis, cmpl_ega3d_rules), cmpl_ega3d_rules);
    // fmt::println("lhs:");
    // print_mvec(mv3d_coeff_svBps, lhs);
    // fmt::println("");

    // rhs
    auto rhs = apply_rules_to_mv(basis, cmpl_ega3d_rules);
    // fmt::println("rhs:");
    // print_mvec(mv3d_coeff_svBps, rhs);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega3d_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "mv << mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e << mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u << mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << ps -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps << mv -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << s -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u << mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << ps -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps << mv_e -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << ps -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps << mv_u -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << bivec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec << mv_e -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << bivec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec << mv_u -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec << mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << vec -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec << mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e << s -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << s -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << ps -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec << ps -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << ps -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << s -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec << bivec -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec << vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec << s -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec << vec -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec << s -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s << s -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega3d_rcontract()
{

    std::string prd_name = "ega3d right contraction";

    // (A >> B) = rwdg(A, rdual(B)) -> in ega3d = rwdg(A,cmpl(B))
    //
    // with rwdg(A,B) = cmpl( wdg(cmpl(A),cmpl(B)) ) in ega3d

    // create the complement from the input multivector (=forward transformation)
    auto basis = mv3d_basis;

    // lhs
    auto lhs = apply_rules_to_mv(basis, cmpl_ega3d_rules);
    // fmt::println("lhs:");
    // print_mvec(mv3d_coeff_svBps, lhs);
    // fmt::println("");

    // lhs (dualize as inner function before complement)
    auto rhs =
        apply_rules_to_mv(apply_rules_to_mv(basis, cmpl_ega3d_rules), cmpl_ega3d_rules);
    // fmt::println("rhs:");
    // print_mvec(mv3d_coeff_svBps, rhs);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega3d_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "mv >> mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv >> mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv >> mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv >> ps -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> bivec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> vec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> s -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s >> mv -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u >> mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u >> mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps >> mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u >> ps -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps >> mv_u -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> bivec -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec >> mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u >> bivec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec >> mv_u -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec >> mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u >> vec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec >> mv_u -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s >> mv_e -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u >> s -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M_odd, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s >> mv_u -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> ps -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec >> bivec -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec >> vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec >> s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec >> vec -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec >> s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s >> s -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_ega3d_lexpand()
{

    std::string prd_name = "ega3d left expansion";
    auto basis = mv3d_basis;

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_ega3d_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3d_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis, wdg_str), wdg_ega3d_rules);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "lexpand(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(s,mv) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_e,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_u,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_e,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_u,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_u,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_e,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(bivec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_u,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(bivec,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_e,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_u,vec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,mv_u) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_e,s) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(s,mv_e) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_u,s) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(s,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(bivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(bivec,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lexpand(bivec,vec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(bivec,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,s) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(s,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega3d_rexpand()
{

    std::string prd_name = "ega3d right expansion";
    auto basis = mv3d_basis;

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_ega3d_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3d_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis, basis_cmpl_func, wdg_str), wdg_ega3d_rules);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "rexpand(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,mv) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(s,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_u,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_u,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_u,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,mv_u) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,bivec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(bivec,mv_e) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_u,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(bivec,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_u,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(s,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_u,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(s,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(bivec,ps) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(s,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(bivec,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(bivec,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rexpand(vec,bivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(s,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(s,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(s,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_ega3d_dot()
{

    std::string prd_name = "ega3d inner product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str), dot_ega3d_rules);
    auto basis = mv3d_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv,mv) -> s");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "dot(mv_e,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "dot(mv_u,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(ps,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(bivec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega3d_rwdg()
{

    // Checking the roundtrip cmpl(cmpl(mvec)) == mvec
    //
    // fmt::println("mv_basis:");
    // print_mvec(mv3d_coeff_svBps, mv3d_basis);
    // fmt::println("");
    // auto double_cmpl_func = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
    // fmt::println("double_cmpl 1st:");
    // print_mvec(mv3d_coeff_svBps, double_cmpl_func);
    // fmt::println("");
    // double_cmpl_func = apply_rules_to_mv(double_cmpl_func, cmpl_ega3d_rules);
    // fmt::println("double_cmpl 2nd:");
    // print_mvec(mv3d_coeff_svBps, double_cmpl_func);
    // fmt::println("");

    std::string prd_name = "ega3d regressive wedge product";
    auto basis = mv3d_basis;

    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv3d_coeff_svBps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_ega3d_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3d_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str),
        wdg_ega3d_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "rwdg(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_M, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_A_odd, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_B_odd);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(ps,ps) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(ps,bivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(bivec,ps) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(ps,vec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(vec,ps) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(ps,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(s,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(bivec,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(bivec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(vec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega3d_rotor_gpr()
{

    std::string prd_name = "ega3d sandwich product";
    fmt::println("ega3d sandwich product: rotor * object * rev(rotor):");
    fmt::println("");

    auto basis = mv3d_basis;
    // fmt::println("mv_basis for sandwich product:");
    // print_mvec(mv3d_coeff_svps, basis);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str), gpr_ega3d_rules);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // first product between multivectors in basis_tab (R * v)
    fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u_tmp");
    auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_R_even, mv3d_coeff_svBps);
    // fmt::println("{} - product table with coeffs:", prd_name);
    // print_prd_tab(prd_tab);
    // fmt::println("");

    auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec,
                                        brace_switch::use_braces);
    fmt::println("mv_u_tmp:");
    print_mvec(mv_u_tmp, basis);
    fmt::println("");

    // second product between multivectors for the product v * rev(R)
    fmt::println("{}:", prd_name + space_str + "mv_u_tmp * rev(mv_e) -> mv_u_res");
    auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv3d_coeff_R_rev_even);
    // fmt::println("prd_tab_v:");
    // print_prd_tab(prd_tab_v);
    // fmt::println("");

    auto mv_u_res_v =
        get_mv_from_prd_tab(prd_tab_v, basis, filter_3d::mv_u, filter_3d::mv_e);
    print_mvec(mv_u_res_v, basis);
    fmt::println("");

    ////

    // first product between multivectors in basis_tab (R * B)
    fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e_tmp");
    auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec,
                                        brace_switch::use_braces);
    fmt::println("mv_e_tmp:");
    print_mvec(mv_e_tmp, basis);
    fmt::println("");

    // second product between multivectors for the product B * rev(R)
    fmt::println("{}:", prd_name + space_str + "mv_e_tmp * rev(mv_e) -> mv_e_res");
    auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv3d_coeff_R_rev_even);
    // fmt::println("prd_tab_B:");
    // print_prd_tab(prd_tab_B);
    // fmt::println("");

    auto mv_e_res_B =
        get_mv_from_prd_tab(prd_tab_B, basis, filter_3d::mv_e, filter_3d::mv_e);
    print_mvec(mv_e_res_B, basis);
    fmt::println("");

    // TODO: (requires string parsing and string manipulation for calculation product
    // expressions by expanding brackets, simplifying and regrouping expressions)
    //
    // manipulate the expression automatically in order to regroup the terms
    // and show by regrouping that
    //   - the pseudoscalar part vanishes (and doesn't need to be calculated)
    //   - that the terms can be rearranged automatically such that a matrix
    //     transformation transformed_vector = matrix * vector can be derived from the
    //     rotor. This is useful when the same transformation should be done for many
    //     objects. This will allow for the transformation to be computed with the minimum
    //     number of required operations (in contrast to somewhat more effort when doing
    //     this with the orgininal rotor expression directly).
    //     Despite this need for optimization to get minimum computational effort, the
    //     rotor formulation is still very meaningful, because of ease of derivation, ease
    //     of rotation interpolation and helpful geometric interpretation when compared to
    //     mere matrix based modeling.
    //
    fmt::println("-------------------------------------------------------------------\n");

    return;
}