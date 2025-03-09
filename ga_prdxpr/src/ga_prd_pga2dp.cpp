// author: Daniel Hug, 2025

#include "ga_prd_pga2dp.hpp"
#include "ga_prdxpr_common.hpp"

using fmt::println;

void generate_and_print_pga2dp_gpr()
{

    // pga2dp geometric product
    std::string prd_name = "pga2dp geometric product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str), gpr_pga2dp_rules);
    auto basis = mv2dp_basis;

    // println("equal = {}", basis_tab == gpr_pga2dp_basis_tab);

    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    println("{}:", prd_name + space_str + "mv * mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv * mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "mv_e * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv * mv_u -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "mv_u * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv * ps -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "ps * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "mv * bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv * vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv * s -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "mv_e * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv_u * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv_e * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "mv_u * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv_e * ps -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "ps * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv_u * ps -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "ps * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "bivec * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv_u * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "bivec * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "mv_e * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "vec * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "mv_u * vec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec * mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_e * s -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s * mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u * s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s * mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "ps * ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps * bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "bivec * ps -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps * vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "vec * ps -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps * s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s * ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "vec * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec * s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s * bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "vec * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "vec * s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "s * s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga2dp_wdg()
{

    // pga2dp wedge product
    std::string prd_name = "pga2dp wedge product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, wdg_str), wdg_pga2dp_rules);
    auto basis = mv2dp_basis;

    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    println("{}:", prd_name + space_str + "mv ^ mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "mv ^ mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_e ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv ^ mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_u ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "ps ^ mv -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv ^ bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv ^ vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv ^ s -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_u ^ mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "ps ^ mv_e -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u ^ ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "ps ^ mv_u -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e ^ bivec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec ^ mv_e -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u ^ bivec -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec ^ mv_u -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e ^ vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec ^ mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u ^ vec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e ^ s -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u ^ s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "ps ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps ^ bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "bivec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps ^ vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "vec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps ^ s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec ^ bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec ^ vec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "vec ^ bivec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec ^ s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s ^ bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "vec ^ vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "vec ^ s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s ^ vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "s ^ s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga2dp_lcontr()
{

    // pga2dp left contraction
    std::string prd_name = "pga2dp left contraction";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, lcontr_str),
        lcontr_pga2dp_rules);
    auto basis = mv2dp_basis;

    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    // println("{}:", prd_name + space_str + "mv << mv -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv << mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_e << mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv << mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_u << mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv << ps -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "ps << mv -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv << bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec << mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv << vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec << mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv << s -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s << mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e << mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u << mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e << mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_u << mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e << ps -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "ps << mv_e -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u << ps -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "ps << mv_u -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e << bivec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec << mv_e -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u << bivec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec << mv_u -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e << vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec << mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u << vec -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec << mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_e << s -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s << mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u << s -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s << mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "ps << ps -> 0");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps << bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "bivec << ps -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps << vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "vec << ps -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps << s -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s << ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec << bivec -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec << vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "vec << bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec << s -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s << bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "vec << vec -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "vec << s -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s << vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "s << s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga2dp_rcontr()
{

    // pga2dp right contraction
    std::string prd_name = "pga2dp right contraction";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, rcontr_str),
        rcontr_pga2dp_rules);
    auto basis = mv2dp_basis;

    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    // println("{}:", prd_name + space_str + "mv >> mv -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv >> mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_e >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv >> mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_u >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv >> ps -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "ps >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv >> bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv >> vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv >> s -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s >> mv -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e >> mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u >> mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e >> mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "mv_u >> mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e >> ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "ps >> mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u >> ps -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "ps >> mv_u -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e >> bivec -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec >> mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u >> bivec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bivec >> mv_u -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e >> vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec >> mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u >> vec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "vec >> mv_u -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_e >> s -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s >> mv_e -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "mv_u >> s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "s >> mv_u -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "ps >> ps -> 0");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps >> bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "bivec >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps >> vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "vec >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "ps >> s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec >> bivec -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec >> vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "vec >> bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "bivec >> s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s >> bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "vec >> vec -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "vec >> s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "s >> vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "s >> s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga2dp_dot()
{

    // pga2dp scalar product
    std::string prd_name = "pga2dp scalar product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str), dot_pga2dp_rules);
    auto basis = mv2dp_basis;

    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    println("{}:", prd_name + space_str + "dot(mv,mv) -> s");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "dot(mv_e,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "dot(mv_u,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "dot(ps,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "dot(bivec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "dot(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "dot(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");

    println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga2dp_rdot()
{

    std::string prd_name = "pga2dp regressive dot product (rdot)";
    auto basis = mv2dp_basis;

    // println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
        dot_pga2dp_rules);
    // println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");


    println("{}:", prd_name + space_str + "rdot(ps,ps) -> ps");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "rdot(bivec,bivec) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "rdot(vec,vec) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "rdot(s,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    print_mvec(prd_mv, basis);
    println("");


    println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_rwdg()
{

    // Checking the roundtrip cmpl(cmpl(mvec)) == mvec
    //
    // println("mv_basis:");
    // print_mvec(mv2dp_coeff_svBps, mv2dp_basis);
    // println("");
    // auto double_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
    // println("double_cmpl 1st:");
    // print_mvec(mv2dp_coeff_svBps, double_cmpl_func);
    // println("");
    // double_cmpl_func = apply_rules_to_mv(double_cmpl_func, cmpl_pga2dp_rules);
    // println("double_cmpl 2nd:");
    // print_mvec(mv2dp_coeff_svBps, double_cmpl_func);
    // println("");

    std::string prd_name = "pga2dp regressive wedge product";
    auto basis = mv2dp_basis;

    // println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str),
        wdg_pga2dp_rules);
    // println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    // create specific products
    // println("{}:", prd_name + space_str + "rwdg(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(s,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_e,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_u,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_e,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(mv_u,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_u,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_e,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(bivec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_u,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(bivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_e,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(vec,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_u,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(vec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(mv_u,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(s,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "rwdg(ps,bivec) -> bivec");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "rwdg(bivec,ps) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "rwdg(ps,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(vec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(ps,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(s,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "rwdg(bivec,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    println("{}:", prd_name + space_str + "rwdg(bivec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    println("{}:", prd_name + space_str + "rwdg(vec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "rwdg(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(vec,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rwdg(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rwdg(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_rgpr()
{

    std::string prd_name = "pga2dp regressive geometric product";
    auto basis = mv2dp_basis;

    // println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBps, basis);
    // println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // println("cmpl:");
    // print_mvec(mv2dp_coeff_svBps, basis_cmpl_func);
    // println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
        gpr_pga2dp_rules);
    // println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    // create specific products
    println("{}:", prd_name + space_str + "rgpr(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "rgpr(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(s,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "rgpr(mv_e,mv_e) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "rgpr(mv_u,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv_e,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    println("{}:", prd_name + space_str + "rgpr(mv_u,mv_e) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "rgpr(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv_u,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "rgpr(mv_e,bivec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");

    // println("{}:", prd_name + space_str + "rgpr(bivec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv_u,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(bivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "rgpr(mv_e,vec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");

    // println("{}:", prd_name + space_str + "rgpr(vec,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv_u,vec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(vec,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(mv_e,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(s,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(mv_u,s) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(s,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(ps,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(bivec,ps) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(ps,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(vec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(ps,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(s,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "rgpr(bivec,bivec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "rgpr(bivec,vec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(vec,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(bivec,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(s,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");


    println("{}:", prd_name + space_str + "rgpr(vec,vec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    print_mvec(prd_mv, basis);
    println("");


    // println("{}:", prd_name + space_str + "rgpr(vec,s) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "rgpr(s,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps, mv2dp_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "rgpr(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBps1, mv2dp_coeff_svBps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_bulk_contraction()
{

    // provides the (right) bulk_contraction

    std::string prd_name = "pga2dp bulk_contraction";
    auto basis = mv2dp_basis;
    // println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBtps, basis);
    // println("");

    // check on sequence of applicationn of dualization and complement:

    // // complement first and then dualize -> wrong!
    // auto rhs_cmpl = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // println("rhs_cmpl (cmpl applied):");
    // print_mvec(mv2dp_coeff_svBtps, rhs_cmpl);
    // println("");

    // auto rhs_cmpl_dual = apply_rules_to_mv(rhs_cmpl, bulk_dual_pga2dp_rules);
    // println("rhs_mv (cmpl + bulk_dualized):");
    // print_mvec(mv2dp_coeff_svBtps, rhs_cmpl_dual);
    // println("");

    // // dualized fist and then complement -> right!
    // auto rhs_dual = apply_rules_to_mv(basis, bulk_dual_pga2dp_rules);
    // println("rhs_dual (bulk_dualized):");
    // print_mvec(mv2dp_coeff_svBtps, rhs_dual);
    // println("");

    // auto rhs_dual_cmpl = apply_rules_to_mv(rhs_dual, cmpl_pga2dp_rules);
    // println("rhs_dual_cmpl (bulk_dualized + cmpl):");
    // print_mvec(mv2dp_coeff_svBtps, rhs_dual_cmpl);
    // println("");

    // => dualization must be the inner operation:
    //    cmpl( wdg( cmpl(a) , cmpl(b) ) ) = rwdg(a, b)
    //
    //    thus: rwdg(a, bulk_dual(b) ) = cmpl( wdg( cmpl(a) , cmpl( bulk_dual(b) ) ) )


    // create the complement from the input multivector (=forward transformation)
    // lhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // println("cmpl:");
    // print_mvec(mv2dp_coeff_svBtps, basis_cmpl_func);
    // println("");

    // rhs (dualize as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, bulk_dual_pga2dp_rules), cmpl_pga2dp_rules);
    // println("cmpl:");
    // print_mvec(mv2dp_coeff_svBtps, basis_cmpl_func);
    // println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_dual_cmpl_func, wdg_str),
        wdg_pga2dp_rules);
    // println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    // create specific products
    // println("{}:", prd_name + space_str + "bulk_contraction(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(mv,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(s,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(mv_e,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(mv_e,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(bivec,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(mv_e,s) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(s,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(ps,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps1, mv2dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(ps,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(bivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(ps,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(bivec,bivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps1, mv2dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(bivec,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(vec,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(bivec,s) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(vec,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps1, mv2dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(vec,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "bulk_contraction(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "bulk_contraction(s,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps1, mv2dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_weight_contraction()
{

    // provides the (right) weight_contraction

    std::string prd_name = "pga2dp weight_contraction";
    auto basis = mv2dp_basis;
    // println("mv_basis for rwdg:");
    // print_mvec(mv2dp_coeff_svBtps, basis);
    // println("");

    // check on sequence of applicationn of dualization and complement:

    // // complement first and then dualize -> wrong!
    // auto rhs_cmpl = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // println("rhs_cmpl (cmpl applied):");
    // print_mvec(mv2dp_coeff_svBtps, rhs_cmpl);
    // println("");

    // auto rhs_cmpl_dual = apply_rules_to_mv(rhs_cmpl, weight_dual_pga2dp_rules);
    // println("rhs_mv (cmpl + weight_dualized):");
    // print_mvec(mv2dp_coeff_svBtps, rhs_cmpl_dual);
    // println("");

    // // dualized fist and then complement -> right!
    // auto rhs_dual = apply_rules_to_mv(basis, weight_dual_pga2dp_rules);
    // println("rhs_dual (weight_dualized):");
    // print_mvec(mv2dp_coeff_svBtps, rhs_dual);
    // println("");

    // auto rhs_dual_cmpl = apply_rules_to_mv(rhs_dual, cmpl_pga2dp_rules);
    // println("rhs_dual_cmpl (weight_dualized + cmpl):");
    // print_mvec(mv2dp_coeff_svBtps, rhs_dual_cmpl);
    // println("");

    // => dualization must be the inner operation:
    //    cmpl( wdg( cmpl(a) , cmpl(b) ) ) = rwdg(a, b)
    //
    //    thus: rwdg(a, weight_dual(b) ) = cmpl( wdg( cmpl(a) , cmpl( weight_dual(b) )
    //    ) )


    // create the complement from the input multivector (=forward transformation)
    // lhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, cmpl_pga2dp_rules);
    // println("cmpl:");
    // print_mvec(mv2dp_coeff_svBtps, basis_cmpl_func);
    // println("");

    // rhs (dualize as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, weight_dual_pga2dp_rules), cmpl_pga2dp_rules);
    // println("cmpl:");
    // print_mvec(mv2dp_coeff_svBtps, basis_cmpl_func);
    // println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_dual_cmpl_func, wdg_str),
        wdg_pga2dp_rules);
    // println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    // create specific products
    // println("{}:", prd_name + space_str + "weight_contraction(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv,trivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::trivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(trivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::trivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(vec,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(s,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_e,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_u,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_e,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(mv_u,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_e,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_u,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_e,trivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_even, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::trivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(trivec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::trivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_u,trivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_M_uneven, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::trivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(trivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::trivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_e,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(bivec,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_u,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(bivec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_e,vec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(vec,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_u,vec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(vec,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_even, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_e);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(mv_u,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_A_uneven, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_u, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(s,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::mv_u);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(ps,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps1, mv2dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(ps,trivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::trivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(trivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::trivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(ps,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(bivec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(ps,vec) -> trivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::ps, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::ps);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(trivec,trivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps1, mv2dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::trivec,
    // filter_2dp::trivec); print_mvec(prd_mv, basis); println("");


    // println("{}:", prd_name + space_str + "weight_contraction(trivec,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::trivec,
    // filter_2dp::bivec); print_mvec(prd_mv, basis); println("");

    // println("{}:", prd_name + space_str + "weight_contraction(bivec,trivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec,
    // filter_2dp::trivec); print_mvec(prd_mv, basis); println("");


    // println("{}:", prd_name + space_str + "weight_contraction(trivec,vec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::trivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(vec,trivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::trivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(trivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::trivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(s,trivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::trivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(bivec,bivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps1, mv2dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(bivec,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(vec,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::bivec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::bivec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(vec,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps1, mv2dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::vec, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    // println("{}:", prd_name + space_str + "weight_contraction(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps, mv2dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::vec);
    // print_mvec(prd_mv, basis);
    // println("");


    // println("{}:", prd_name + space_str + "weight_contraction(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_svBtps1, mv2dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::s, filter_2dp::s);
    // print_mvec(prd_mv, basis);
    // println("");

    println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_bulk_expansion()
{

    // provides the (right) bulk_expansion

    std::string prd_name = "pga2dp bulk_expansion";
    auto basis = mv2dp_basis;

    auto lhs = basis;
    // println("mv_basis for bulk_expansion:");
    // print_mvec(mv2dp_coeff_svBtps, basis);
    // println("");

    auto rhs = apply_rules_to_mv(basis, bulk_dual_pga2dp_rules);
    // println("rhs for bulk_expansion:");
    // print_mvec(mv2dp_coeff_svBtps, rhs);
    // println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga2dp_rules);
    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga2dp_weight_expansion()
{

    // provides the (right) weight_expansion

    std::string prd_name = "pga2dp weight_expansion";
    auto basis = mv2dp_basis;

    auto lhs = basis;
    // println("mv_basis for weight_expansion:");
    // print_mvec(mv2dp_coeff_svBtps, basis);
    // println("");

    auto rhs = apply_rules_to_mv(basis, weight_dual_pga2dp_rules);
    // println("rhs for weight_expansion:");
    // print_mvec(mv2dp_coeff_svBtps, rhs);
    // println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga2dp_rules);
    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga2dp_rotor_gpr()
{

    std::string prd_name = "pga2dp sandwich product:";
    println("pga2dp sandwich product: rotor * object * rev(rotor):");
    println("");

    auto basis = mv2dp_basis;
    // println("mv_basis for sandwich product:");
    // print_mvec(mv2dp_coeff_svps, basis);
    // println("");

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str), gpr_pga2dp_rules);

    println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    println("");

    // first product between multivectors in basis_tab (R * v)
    println("{}:", prd_name + space_str + "mv_e * vec -> mv_u_tmp");
    auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_R_even, mv2dp_coeff_svBps);
    // println("{} - product table with coeffs:", prd_name);
    // print_prd_tab(prd_tab);
    // println("");

    auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e, filter_2dp::vec,
                                        brace_switch::use_braces);
    println("mv_u_tmp:");
    print_mvec(mv_u_tmp, basis);
    println("");

    // second product between multivectors for the product v * rev(R)
    println("{}:", prd_name + space_str + "mv_u_tmp * rev(mv_e) -> mv_u_res");
    auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv2dp_coeff_R_rev_even);
    // println("prd_tab_v:");
    // print_prd_tab(prd_tab_v);
    // println("");

    auto mv_u_res_v =
        get_mv_from_prd_tab(prd_tab_v, basis, filter_2dp::mv_u, filter_2dp::mv_e);
    print_mvec(mv_u_res_v, basis);
    println("");

    ////

    // first product between multivectors in basis_tab (R * B)
    println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e_tmp");
    auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, basis, filter_2dp::mv_e,
                                        filter_2dp::bivec, brace_switch::use_braces);
    println("mv_e_tmp:");
    print_mvec(mv_e_tmp, basis);
    println("");

    // second product between multivectors for the product B * rev(R)
    println("{}:", prd_name + space_str + "mv_e_tmp * rev(mv_e) -> mv_e_res");
    auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv2dp_coeff_R_rev_even);
    // println("prd_tab_B:");
    // print_prd_tab(prd_tab_B);
    // println("");

    auto mv_e_res_B =
        get_mv_from_prd_tab(prd_tab_B, basis, filter_2dp::mv_e, filter_2dp::mv_e);
    print_mvec(mv_e_res_B, basis);
    println("");

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
    println("-------------------------------------------------------------------\n");

    return;
}