// author: Daniel Hug, 2024

#include "ga_prd_ega3d.hpp"
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
    auto prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * mv_u -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_u * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * ps -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * bivec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * vec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * s -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_u * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * ps -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * ps -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e * s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * s -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * ps -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * ps -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * ps -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * s -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * ps -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * s -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * s -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s * s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

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
    auto prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv -> ps");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ bivec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ vec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ s -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv_e -> ps");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv_u -> 0");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ bivec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv_e -> bivec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ bivec -> ps");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv_u -> ps");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ vec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ s -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ vec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ s -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ vec -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ bivec -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ s -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ s -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s ^ s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega3d_lcontr()
{

    // ega3d left contraction
    std::string prd_name = "ega3d left contraction";

    auto basis_tab =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, lcontr_str),
                           lcontr_ega3d_rules);
    auto basis = mv3d_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv << mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e << mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u << mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << ps -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps << mv -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv << s -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u << mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << ps -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps << mv_e -> 0");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << ps -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps << mv_u -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << bivec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec << mv_e -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << bivec -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec << mv_u -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e << vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec << mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << vec -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec << mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e << s -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u << s -> 0");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << ps -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec << ps -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << vec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << ps -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << s -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << ps -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec << bivec -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec << vec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec << s -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec << vec -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec << s -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s << s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega3d_rcontr()
{

    // ega3d right contraction
    std::string prd_name = "ega3d right contraction";

    auto basis_tab =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, rcontr_str),
                           rcontr_ega3d_rules);
    auto basis = mv3d_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv >> mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv >> mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv >> mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv >> ps -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> bivec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> vec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec >> mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> s -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_M, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s >> mv -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u >> mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u >> mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps >> mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u >> ps -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps >> mv_u -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> bivec -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec >> mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u >> bivec -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec >> mv_u -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e >> vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec >> mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u >> vec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec >> mv_u -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s >> mv_e -> s");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u >> s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s >> mv_u -> 0");
    // prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> ps -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> bivec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> s -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec >> bivec -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec >> vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec >> s -> bivec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> bivec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec >> vec -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec >> s -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> vec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s >> s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega3d_dot()
{

    // ega3d scalar product
    std::string prd_name = "ega3d scalar product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str), dot_ega3d_rules);
    auto basis = mv3d_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv,mv) -> s");
    auto prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,mv_e) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv_e,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,mv_u) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv_u,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(ps,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(bivec,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(vec,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,s) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(s,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv_e,mv_e) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv_u,mv_u) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv_u,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(ps,mv_u) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv_e,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(bivec,mv_e) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::bivec, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv_u,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_uneven, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_u, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(vec,mv_u) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv_e,s) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_A_even, mv3d_coeff_svBps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::mv_e, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(s,mv_e) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps, mv3d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(ps,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::ps, filter_3d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::vec, filter_3d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv3d_coeff_svBps1, mv3d_coeff_svBps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3d::s, filter_3d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}