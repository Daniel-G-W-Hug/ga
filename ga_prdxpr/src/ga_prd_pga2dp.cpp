// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prd_pga2dp.hpp"
#include "ga_prdxpr_common.hpp"

void generate_and_print_pga2dp_gpr()
{

    // pga2dp geometric product
    std::string prd_name = "pga2dp geometric product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str), gpr_pga2dp_rules);
    auto basis = mv2dp_basis;

    // fmt::println("equal = {}", basis_tab == gpr_pga2dp_basis_tab);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv * mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * mv_u -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_u * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * ps -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv * bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv * vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv * s -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_u * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * ps -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * ps -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u * vec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec * mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e * s -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s * mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u * s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s * mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * ps -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * ps -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s * s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2dp_cmt()
{

    // the commutator product is the asymmetric part of the geometric product
    std::string prd_name = "pga2dp commutator product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str), gpr_pga2dp_rules);
    auto basis = mv2dp_basis;

    basis_tab = get_prd_tab_asym(basis_tab); // use the asymmetric part only

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "cmt(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(ps,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(ps,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(ps,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(vec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(vec,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(mv_u,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(ps,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(ps,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(ps,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "cmt(bivec,bivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "cmt(bivec,vec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "cmt(vec,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "cmt(vec,vec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "cmt(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "cmt(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga2dp_wdg()
{

    // pga2dp wedge product
    std::string prd_name = "pga2dp wedge product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, wdg_str), wdg_pga2dp_rules);
    auto basis = mv2dp_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv ^ mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ s -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv_e -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv_u -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ bivec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv_e -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ bivec -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv_u -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ vec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ s -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ vec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ bivec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s ^ s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_dot()
{

    std::string prd_name = "pga2dp inner product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str), dot_pga2dp_rules);
    auto basis = mv2dp_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv,mv) -> s");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "dot(mv_e,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "dot(mv_u,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(ps,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(bivec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga2dp_rdot()
{

    std::string prd_name = "pga2dp regressive inner product";
    auto basis = mv2dp_basis;

    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
        dot_pga2dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(ps,ps) -> ps");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(bivec,bivec) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(vec,vec) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(s,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_rwdg()
{

    // Checking the roundtrip cmpl(cmpl(mvec)) == mvec
    //
    // fmt::println("mv_basis:");
    // print_mvec(mv2dp_coeff_svBps, mv2dp_basis);
    // fmt::println("");
    // auto double_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
    // fmt::println("double_cmpl 1st:");
    // print_mvec(mv2dp_coeff_svBps, double_cmpl_func);
    // fmt::println("");
    // double_cmpl_func = apply_rules_to_mv(double_cmpl_func, cmpl_pga2dp_rules);
    // fmt::println("double_cmpl 2nd:");
    // print_mvec(mv2dp_coeff_svBps, double_cmpl_func);
    // fmt::println("");

    std::string prd_name = "pga2dp regressive wedge product";
    auto basis = mv2dp_basis;

    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str),
        wdg_pga2dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "rwdg(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(ps,bivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(bivec,ps) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(bivec,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(bivec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(vec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_rgpr()
{

    std::string prd_name = "pga2dp regressive geometric product";
    auto basis = mv2dp_basis;

    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
        gpr_pga2dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "rgpr(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,mv_e) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,mv_u) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,mv_u) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,mv_e) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,bivec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,bivec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,vec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,vec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,s) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,ps) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(bivec,bivec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(bivec,vec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(vec,vec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,s) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_rbulk_contract()
{

    std::string prd_name = "pga2dp right bulk contraction";

    // (A >> B) = rwdg(A, rdual(B)) -> in pga2dp = rwdg(A, bulk_dual(B))
    //
    // with rwdg(A,B) = cmpl( wdg(cmpl(A),cmpl(B)) ) in pga2dp

    auto basis = mv2dp_basis;
    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    // check on sequence of applicationn of dualization and complement:

    // // complement first and then dualize -> wrong!
    // auto rhs_cmpl = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("rhs_cmpl (cmpl applied):");
    // print_mvec(mv2dp_coeff_svBps, rhs_cmpl);
    // fmt::println("");

    // auto rhs_cmpl_dual = apply_rules_to_mv(rhs_cmpl, bulk_dual_pga2dp_rules);
    // fmt::println("rhs_mv (cmpl + bulk_dualized):");
    // print_mvec(mv2dp_coeff_svBps, rhs_cmpl_dual);
    // fmt::println("");

    // // dualized fist and then complement -> right!
    // auto rhs_dual = apply_rules_to_mv(basis, bulk_dual_pga2dp_rules);
    // fmt::println("rhs_dual (bulk_dualized):");
    // print_mvec(mv2dp_coeff_svBps, rhs_dual);
    // fmt::println("");

    // auto rhs_dual_cmpl = apply_rules_to_mv(rhs_dual, cmpl_pga2dp_rules);
    // fmt::println("rhs_dual_cmpl (bulk_dualized + cmpl):");
    // print_mvec(mv2dp_coeff_svBps, rhs_dual_cmpl);
    // fmt::println("");

    // => dualization must be the inner operation:
    //    cmpl( wdg( cmpl(a) , cmpl(b) ) ) = rwdg(a, b)
    //
    //    thus: rwdg(a, bulk_dual(b) ) = cmpl( wdg( cmpl(a) , cmpl( bulk_dual(b) ) ) )


    // create the complement from the input multivector (=forward transformation)
    // lhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // rhs (dualize as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, bulk_dual_pga2dp_rules), cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_dual_cmpl_func, wdg_str),
        wdg_pga2dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,bivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,s) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,s) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,vec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,bivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,s) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,bivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,s) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,vec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_rweight_contract()
{

    std::string prd_name = "pga2dp right weight contraction";
    auto basis = mv2dp_basis;
    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    // check on sequence of application of dualization and complement:

    // // complement first and then dualize -> wrong!
    // auto rhs_cmpl = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("rhs_cmpl (cmpl applied):");
    // print_mvec(mv2dp_coeff_svBps, rhs_cmpl);
    // fmt::println("");

    // auto rhs_cmpl_dual = apply_rules_to_mv(rhs_cmpl, weight_dual_pga2dp_rules);
    // fmt::println("rhs_mv (cmpl + weight_dualized):");
    // print_mvec(mv2dp_coeff_svBps, rhs_cmpl_dual);
    // fmt::println("");

    // // dualized fist and then complement -> right!
    // auto rhs_dual = apply_rules_to_mv(basis, weight_dual_pga2dp_rules);
    // fmt::println("rhs_dual (weight_dualized):");
    // print_mvec(mv2dp_coeff_svBps, rhs_dual);
    // fmt::println("");

    // auto rhs_dual_cmpl = apply_rules_to_mv(rhs_dual, cmpl_pga2dp_rules);
    // fmt::println("rhs_dual_cmpl (weight_dualized + cmpl):");
    // print_mvec(mv2dp_coeff_svBps, rhs_dual_cmpl);
    // fmt::println("");

    // => dualization must be the inner operation:
    //    cmpl( wdg( cmpl(a) , cmpl(b) ) ) = rwdg(a, b)
    //
    //    thus: rwdg(a, weight_dual(b) ) = cmpl( wdg( cmpl(a) , cmpl( weight_dual(b) )
    //    ) )


    // create the complement from the input multivector (=forward transformation)
    // lhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // rhs (dualize as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, weight_dual_pga2dp_rules), cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_dual_cmpl_func, wdg_str),
        wdg_pga2dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(vec,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(s,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_e,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_u,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_e,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_u,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(ps,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_u,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(ps,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_e,bivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(bivec,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_u,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(bivec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_e,vec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(vec,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_u,vec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(vec,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv_u,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(s,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(ps,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(ps,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(bivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(ps,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(bivec,bivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(bivec,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(vec,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(vec,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_contract(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_contract(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_lbulk_contract()
{

    std::string prd_name = "pga2dp left bulk contraction";

    // (A << B) = rwdg( ldual(A),B ) -> in pga2dp = rwdg(bulk_dual(A),B)
    //
    // with rwdg(A,B) = cmpl( wdg( cmpl(A),cmpl(B) ) ) in pga2dp

    auto basis = mv2dp_basis;
    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    // check on sequence of applicationn of dualization and complement:

    // // complement first and then dualize -> wrong!
    // auto rhs_cmpl = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("rhs_cmpl (cmpl applied):");
    // print_mvec(mv2dp_coeff_svBps, rhs_cmpl);
    // fmt::println("");

    // auto rhs_cmpl_dual = apply_rules_to_mv(rhs_cmpl, bulk_dual_pga2dp_rules);
    // fmt::println("rhs_mv (cmpl + bulk_dualized):");
    // print_mvec(mv2dp_coeff_svBps, rhs_cmpl_dual);
    // fmt::println("");

    // // dualized fist and then complement -> right!
    // auto rhs_dual = apply_rules_to_mv(basis, bulk_dual_pga2dp_rules);
    // fmt::println("rhs_dual (bulk_dualized):");
    // print_mvec(mv2dp_coeff_svBps, rhs_dual);
    // fmt::println("");

    // auto rhs_dual_cmpl = apply_rules_to_mv(rhs_dual, cmpl_pga2dp_rules);
    // fmt::println("rhs_dual_cmpl (bulk_dualized + cmpl):");
    // print_mvec(mv2dp_coeff_svBps, rhs_dual_cmpl);
    // fmt::println("");

    // => dualization must be the inner operation:
    //    cmpl( wdg( cmpl(a) , cmpl(b) ) ) = rwdg(a, b)
    //
    //    thus: rwdg(a, bulk_dual(b) ) = cmpl( wdg( cmpl(a) , cmpl( bulk_dual(b) ) ) )


    // create the complement from the input multivector (=forward transformation)

    // lhs (dualize as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, bulk_dual_pga2dp_rules), cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // rhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_dual_cmpl_func, basis_cmpl_func, wdg_str),
        wdg_pga2dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,bivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,ps) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps, vec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,ps) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,ps) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,vec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,bivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,vec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_lweight_contract()
{

    std::string prd_name = "pga2dp left weight contraction";
    auto basis = mv2dp_basis;
    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    // check on sequence of application of dualization and complement:

    // // complement first and then dualize -> wrong!
    // auto rhs_cmpl = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("rhs_cmpl (cmpl applied):");
    // print_mvec(mv2dp_coeff_svBps, rhs_cmpl);
    // fmt::println("");

    // auto rhs_cmpl_dual = apply_rules_to_mv(rhs_cmpl, weight_dual_pga2dp_rules);
    // fmt::println("rhs_mv (cmpl + weight_dualized):");
    // print_mvec(mv2dp_coeff_svBps, rhs_cmpl_dual);
    // fmt::println("");

    // // dualized fist and then complement -> right!
    // auto rhs_dual = apply_rules_to_mv(basis, weight_dual_pga2dp_rules);
    // fmt::println("rhs_dual (weight_dualized):");
    // print_mvec(mv2dp_coeff_svBps, rhs_dual);
    // fmt::println("");

    // auto rhs_dual_cmpl = apply_rules_to_mv(rhs_dual, cmpl_pga2dp_rules);
    // fmt::println("rhs_dual_cmpl (weight_dualized + cmpl):");
    // print_mvec(mv2dp_coeff_svBps, rhs_dual_cmpl);
    // fmt::println("");

    // => dualization must be the inner operation:
    //    cmpl( wdg( cmpl(a) , cmpl(b) ) ) = rwdg(a, b)
    //
    //    thus: rwdg(a, weight_dual(b) ) = cmpl( wdg( cmpl(a) , cmpl( weight_dual(b) )
    //    ) )


    // create the complement from the input multivector (=forward transformation)

    // lhs (dualize as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, weight_dual_pga2dp_rules), cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // rhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_dual_cmpl_func, basis_cmpl_func, wdg_str),
        wdg_pga2dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(ps,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(s,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_e,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_u,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_e,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_u,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_e,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(ps,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_u,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(ps,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_e,bivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(bivec,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_u,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(bivec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_e,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(vec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_u,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(vec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(mv_u,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(s,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(ps,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(ps,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(bivec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(ps,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(vec,ps) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(bivec,bivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(bivec,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(vec,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(vec,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_contract(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_contract(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_rbulk_expand()
{

    std::string prd_name = "pga2dp right bulk expansion";
    auto basis = mv2dp_basis;

    auto lhs = basis;
    // fmt::println("mv_basis for rbulk_expand:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    auto rhs = apply_rules_to_mv(basis, bulk_dual_pga2dp_rules);
    // fmt::println("rhs for rbulk_expand:");
    // print_mvec(mv2dp_coeff_svBps, rhs);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(ps,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(bivec,mv) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(s,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv_e,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(ps,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv_e,bivec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(bivec,mv_e) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(mv_e,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(s,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(ps,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(ps,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(bivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(bivec,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(bivec,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(vec,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(s,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(vec,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(s,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_expand(s,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_rweight_expand()
{

    std::string prd_name = "pga2dp right weight expansion";
    auto basis = mv2dp_basis;

    auto lhs = basis;
    // fmt::println("mv_basis for rweight_expand:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    auto rhs = apply_rules_to_mv(basis, weight_dual_pga2dp_rules);
    // fmt::println("rhs for rweight_expand:");
    // print_mvec(mv2dp_coeff_svBps, rhs);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(ps,mv) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(s,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_e,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_u,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_e,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_u,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(ps,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_u,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(ps,mv_u) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_e,bivec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(bivec,mv_e) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_u,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(bivec,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_e,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(vec,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_u,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(vec,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(s,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(mv_u,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(s,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(ps,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(bivec,ps) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(ps,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(vec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(s,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(bivec,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(bivec,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(vec,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(s,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(vec,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rweight_expand(s,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rweight_expand(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga2dp_lbulk_expand()
{

    std::string prd_name = "pga2dp left bulk expansion";
    auto basis = mv2dp_basis;

    auto lhs = apply_rules_to_mv(basis, bulk_dual_pga2dp_rules);
    // fmt::println("rhs for lbulk_expand:");
    // print_mvec(mv2dp_coeff_svBps, rhs);
    // fmt::println("");

    auto rhs = basis;
    // fmt::println("mv_basis for lbulk_expand:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(ps,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(s,mv) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv_e,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(ps,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv_e,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(bivec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(mv_e,s) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(s,mv_e) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(ps,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(ps,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(bivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(bivec,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(bivec,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(vec,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(bivec,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(vec,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(vec,s) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_expand(s,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_lweight_expand()
{

    std::string prd_name = "pga2dp left weight expansion";
    auto basis = mv2dp_basis;

    auto lhs = apply_rules_to_mv(basis, weight_dual_pga2dp_rules);
    // fmt::println("rhs for lweight_expand:");
    // print_mvec(mv2dp_coeff_svBps, rhs);
    // fmt::println("");

    auto rhs = basis;
    // fmt::println("mv_basis for lweight_expand:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(s,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_e,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_u,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_e,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_u,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_u,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_e,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(bivec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_u,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(bivec,mv_u) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_e,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(vec,mv_e) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_u,vec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(vec,mv_u) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_e,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(mv_u,s) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(s,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(ps,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(bivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(ps,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(ps,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(bivec,bivec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(bivec,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(vec,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(bivec,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(vec,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(vec,s) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lweight_expand(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lweight_expand(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_motor_rgpr()
{

    std::string prd_name = "pga2dp regressive sandwich product";
    fmt::println("pga2dp regressive sandwich product - motor * object * rrev(motor):");
    fmt::println("");

    auto basis = mv2dp_basis;
    // fmt::println("mv_basis for regressive sandwich product:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
        gpr_pga2dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");


    // first product between multivectors in basis_tab (R * v)
    fmt::println("{}:", prd_name + space_str + "rgpr(mv_u, vec) -> mv_u_tmp");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_R_uneven, mv2dp_coeff_svBps);
    // fmt::println("{} - product table with coeffs:", prd_name);
    // print_prd_tab(prd_tab);
    // fmt::println("");

    auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec,
                                        brace_switch::use_braces);
    fmt::println("mv_u_tmp:");
    print_mvec(mv_u_tmp, basis);
    fmt::println("");

    // second product between multivectors for the product v * rrev(R)
    fmt::println("{}:", prd_name + space_str + "mv_u_tmp * rrev(mv_u) -> mv_u_res");
    auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv2dp_coeff_R_rrev_uneven);
    // fmt::println("prd_tab_v:");
    // print_prd_tab(prd_tab_v);
    // fmt::println("");

    auto mv_u_res_v =
        get_mv_from_prd_tab(prd_tab_v, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    print_mvec(mv_u_res_v, basis);
    fmt::println("");

    ////

    // first product between multivectors in basis_tab (R * B)
    fmt::println("{}:", prd_name + space_str + "mv_u * bivec -> mv_e_tmp");
    auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u,
                                        filter_2dp::bivec, brace_switch::use_braces);
    fmt::println("mv_e_tmp:");
    print_mvec(mv_e_tmp, basis);
    fmt::println("");

    // second product between multivectors for the product B * rrev(R)
    fmt::println("{}:", prd_name + space_str + "mv_e_tmp * rrev(mv_u) -> mv_e_res");
    auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv2dp_coeff_R_rrev_uneven);
    // fmt::println("prd_tab_B:");
    // print_prd_tab(prd_tab_B);
    // fmt::println("");

    auto mv_e_res_B =
        get_mv_from_prd_tab(prd_tab_B, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    print_mvec(mv_e_res_B, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}