// author: Daniel Hug, 2025

#include "ga_prd_pga3dp.hpp"
#include "ga_prdxpr_common.hpp"

void generate_and_print_pga3dp_gpr()
{

    std::string prd_name = "pga3dp geometric product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, mul_str), gpr_pga3dp_rules);
    auto basis = mv3dp_basis;

    // fmt::println("equal = {}", basis_tab == gpr_pga3dp_basis_tab);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv * mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv * mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv * ps -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv * trivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "trivec * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv * bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv * vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv * s -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s * mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * mv_u -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_u * mv_e -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e * ps -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps * mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u * ps -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps * mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * trivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "trivec * mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u * trivec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "trivec * mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec * mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_u * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * mv_u -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_u);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec * mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u * vec -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec * mv_u -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e * s -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s * mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u * s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s * mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "ps * ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "ps * trivec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "trivec * ps -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "ps * bivec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec * ps -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "ps * vec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec * ps -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "trivec * trivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "trivec * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "bivec * trivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "trivec * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * trivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "trivec * s -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * trivec -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * bivec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * vec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * bivec -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec * s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s * s -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga3dp_wdg()
{

    std::string prd_name = "pga3dp wedge product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, wdg_str), wdg_pga3dp_rules);
    auto basis = mv3dp_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv ^ mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ mv_e -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ mv_u -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ trivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "trivec ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ bivec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ vec -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv ^ s -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "mv_u ^ mv_e -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ ps -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv_e -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "ps ^ mv_u -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ trivec -> trivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "trivec ^ mv_e -> trivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ trivec -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "trivec ^ mv_u -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ bivec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv_e -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ bivec -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ mv_u -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ vec -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv_e -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ vec -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ mv_u -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_e ^ s -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv_e -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "mv_u ^ s -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "s ^ mv_u -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "ps ^ ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "ps ^ trivec -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "trivec ^ ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "ps ^ bivec -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "ps ^ vec -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "vec ^ ps -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "trivec ^ trivec -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec,
    // filter_3dp::trivec); print_mvec(prd_mv, basis); fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "trivec ^ bivec -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec,
    // filter_3dp::bivec); print_mvec(prd_mv, basis); fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "bivec ^ trivec -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec,
    // filter_3dp::trivec); print_mvec(prd_mv, basis); fmt::println("");


    fmt::println("{}:", prd_name + space_str + "trivec ^ vec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ trivec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "trivec ^ s -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ trivec -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ bivec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ vec -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ bivec -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "bivec ^ s -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ bivec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ vec -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s ^ s -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_dot()
{

    std::string prd_name = "pga3dp inner product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, mul_str), dot_pga3dp_rules);
    auto basis = mv3dp_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv,mv) -> s");
    auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "dot(mv_e,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "dot(mv_u,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(ps,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(trivec,trivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(bivec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga3dp_rdot()
{

    std::string prd_name = "pga3dp regressive inner product";
    auto basis = mv3dp_basis;

    // fmt::println("mv_basis for rdot:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
        dot_pga3dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(mv,mv) -> ps");
    auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(ps,ps) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(trivec,trivec) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(bivec,bivec) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(vec,vec) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rdot(s,s) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_rwdg()
{

    // Checking the roundtrip lcmpl(rcmpl(mvec)) == mvec
    //
    // fmt::println("mv_basis:");
    // print_mvec(mv3dp_coeff_svBtps, mv3dp_basis);
    // fmt::println("");
    // auto double_cmpl_func = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
    // fmt::println("double_cmpl 1st:");
    // print_mvec(mv3dp_coeff_svBtps, double_cmpl_func);
    // fmt::println("");
    // double_cmpl_func = apply_rules_to_mv(double_cmpl_func, lcmpl_pga3dp_rules);
    // fmt::println("double_cmpl 2nd:");
    // print_mvec(mv3dp_coeff_svBtps, double_cmpl_func);
    // fmt::println("");

    std::string prd_name = "pga3dp regressive wedge product";
    auto basis = mv3dp_basis;

    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str),
        wdg_pga3dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "rwdg(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,trivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(trivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,trivec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(trivec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,trivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(trivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_u,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(ps,trivec) -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(trivec,ps) -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,ps) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(trivec,trivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(trivec,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(bivec,trivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(trivec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(vec,trivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(trivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,trivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(bivec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(bivec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,bivec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_rgpr()
{

    std::string prd_name = "pga3dp regressive geometric product";
    auto basis = mv3dp_basis;

    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
        gpr_pga3dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "rgpr(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,trivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(trivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,mv_e) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,mv_e) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,trivec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(trivec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,trivec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(trivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,bivec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,bivec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,vec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,vec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_e,s) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(mv_u,s) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,trivec) -> trivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(trivec,ps) -> trivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,ps) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(ps,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(trivec,trivec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(trivec,bivec) -> mv_u");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,trivec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec,
    // filter_3dp::trivec); print_mvec(prd_mv, basis); fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(trivec,vec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,trivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(trivec,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,trivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(bivec,bivec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,vec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,bivec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(bivec,s) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,bivec) -> bivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rgpr(vec,vec) -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(vec,s) -> trivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rgpr(s,vec) -> trivec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rgpr(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_rbulk_contract()
{

    std::string prd_name = "pga3dp right bulk contraction";

    // (A >> B) = rwdg(A, rdual(B)) -> in pga3dp = rwdg(A, rbulk_dual(B))
    //
    // with rwdg(A,B) = lcmpl( wdg(rcmpl(A),rcmpl(B)) ) in pga3dp

    auto basis = mv3dp_basis;
    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    // check on sequence of applicationn of dualization and complement:

    // // complement first and then dualize -> wrong!
    // auto rhs_cmpl = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("rhs_cmpl (rcmpl applied):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_cmpl);
    // fmt::println("");

    // auto rhs_cmpl_dual = apply_rules_to_mv(rhs_cmpl, right_bulk_dual_pga3dp_rules);
    // fmt::println("rhs_mv (rcmpl + bulk_dualized):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_cmpl_dual);
    // fmt::println("");

    // // dualized fist and then complement -> right!
    // auto rhs_dual = apply_rules_to_mv(basis, right_bulk_dual_pga3dp_rules);
    // fmt::println("rhs_dual (bulk_dualized):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_dual);
    // fmt::println("");

    // auto rhs_dual_cmpl = apply_rules_to_mv(rhs_dual, rcmpl_pga3dp_rules);
    // fmt::println("rhs_dual_cmpl (bulk_dualized + rcmpl):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_dual_cmpl);
    // fmt::println("");

    // => dualization must be the inner operation:
    //    lcmpl( wdg( rcmpl(a) , rcmpl(b) ) ) = rwdg(a, b)
    //
    //    thus: rwdg(a, bulk_dual(b) ) = lcmpl( wdg( rcmpl(a) , rcmpl( bulk_dual(b) ) ) )


    // create the complement from the input multivector (=forward transformation)
    // lhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");

    // rhs (dualize as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, right_bulk_dual_pga3dp_rules), rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_dual_cmpl_func, wdg_str),
        wdg_pga3dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,trivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(trivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_u,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_u,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_u,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,trivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(trivec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_u,trivec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(trivec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_u,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,vec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_u,vec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_e,s) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(mv_u,s) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,trivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(trivec,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,bivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,vec) -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(ps,s) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(trivec,trivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(trivec,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,trivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(trivec,vec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,trivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(trivec,s) -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,trivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,vec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,bivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(bivec,s) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,bivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(vec,s) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,vec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rbulk_contract(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_rweight_contract()
{

    std::string prd_name = "pga3dp right weight contraction";
    auto basis = mv3dp_basis;
    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    // check on sequence of applicationn of dualization and complement:

    // // complement first and then dualize -> wrong!
    // auto rhs_cmpl = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("rhs_cmpl (rcmpl applied):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_cmpl);
    // fmt::println("");

    // auto rhs_cmpl_dual = apply_rules_to_mv(rhs_cmpl, right_weight_dual_pga3dp_rules);
    // fmt::println("rhs_mv (rcmpl + weight_dualized):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_cmpl_dual);
    // fmt::println("");

    // // dualized fist and then complement -> right!
    // auto rhs_dual = apply_rules_to_mv(basis, right_weight_dual_pga3dp_rules);
    // fmt::println("rhs_dual (weight_dualized):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_dual);
    // fmt::println("");

    // auto rhs_dual_cmpl = apply_rules_to_mv(rhs_dual, rcmpl_pga3dp_rules);
    // fmt::println("rhs_dual_cmpl (weight_dualized + rcmpl):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_dual_cmpl);
    // fmt::println("");

    // => dualization must be the inner operation:
    //    lcmpl( wdg( rcmpl(a) , rcmpl(b) ) ) = rwdg(a, b)
    //
    //    thus: rwdg(a, weight_dual(b) ) = lcmpl( wdg( rcmpl(a) , rcmpl( weight_dual(b) )
    //    ) )


    // create the complement from the input multivector (=forward transformation)
    // lhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");

    // rhs (dualize as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, right_weight_dual_pga3dp_rules), rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_dual_cmpl_func, wdg_str),
        wdg_pga3dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga3dp_lweight_contract()
{

    std::string prd_name = "pga3dp left weight contraction";
    auto basis = mv3dp_basis;
    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    // check on sequence of applicationn of dualization and complement:

    // // complement first and then dualize -> wrong!
    // auto rhs_cmpl = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("rhs_cmpl (rcmpl applied):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_cmpl);
    // fmt::println("");

    // auto rhs_cmpl_dual = apply_rules_to_mv(rhs_cmpl, right_weight_dual_pga3dp_rules);
    // fmt::println("rhs_mv (rcmpl + weight_dualized):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_cmpl_dual);
    // fmt::println("");

    // // dualized fist and then complement -> right!
    // auto rhs_dual = apply_rules_to_mv(basis, right_weight_dual_pga3dp_rules);
    // fmt::println("rhs_dual (weight_dualized):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_dual);
    // fmt::println("");

    // auto rhs_dual_cmpl = apply_rules_to_mv(rhs_dual, rcmpl_pga3dp_rules);
    // fmt::println("rhs_dual_cmpl (weight_dualized + rcmpl):");
    // print_mvec(mv3dp_coeff_svBtps, rhs_dual_cmpl);
    // fmt::println("");

    // => dualization must be the inner operation:
    //    lcmpl( wdg( rcmpl(a) , rcmpl(b) ) ) = rwdg(a, b)
    //
    //    thus: rwdg(a, weight_dual(b) ) = lcmpl( wdg( rcmpl(a) , rcmpl( weight_dual(b) )
    //    ) )


    // create the complement from the input multivector (=forward transformation)
    // lhs (dualize as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, left_weight_dual_pga3dp_rules), rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");


    // rhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");


    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_dual_cmpl_func, basis_cmpl_func, wdg_str),
        wdg_pga3dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "rweight_contract(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_lbulk_contract()
{

    std::string prd_name = "pga3dp left bulk contraction";

    // (A << B)  = rwdg( ldual(A),B ) -> in pga3dp= rwdg(lbulk_dual(A),B)
    //
    // with rwdg(A,B) = lcmpl( wdg( rcmpl(A),rcmpl(B) ) ) in pga3dp

    // create the complement from the input multivector (=forward transformation)
    auto basis = mv3dp_basis;
    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    // lhs (dualize with left_bulk_dual as inner function before complement)
    auto basis_dual_cmpl_func = apply_rules_to_mv(
        apply_rules_to_mv(basis, left_bulk_dual_pga3dp_rules), rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");

    // rhs
    auto basis_cmpl_func = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");


    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_dual_cmpl_func, basis_cmpl_func, wdg_str),
        wdg_pga3dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,mv_u) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_u,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,mv) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,trivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(trivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,bivec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_u,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_u,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,mv_e) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_u,ps) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,mv_u) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,trivec) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(trivec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_u,trivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_M_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::trivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(trivec,mv_u) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_M_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,bivec) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_u,bivec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::bivec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,mv_u) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,mv_e) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_u,vec) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,mv_u) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_e,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_even, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(mv_u,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_A_uneven, mv3dp_coeff_svBtps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_u, filter_3dp::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,mv_u) -> mv_u");
    // prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_B_uneven);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::mv_u);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,trivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(trivec,ps) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,bivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,ps) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,vec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,ps) -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(ps,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::ps, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,ps) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(trivec,trivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(trivec,bivec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,trivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(trivec,vec) ->  0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,trivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(trivec,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::trivec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,trivec) -> trivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::trivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,bivec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,vec) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,bivec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(bivec,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::bivec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,bivec) -> bivec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::bivec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(vec,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::vec, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,vec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps, mv3dp_coeff_svBtps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "lbulk_contract(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_svBtps1, mv3dp_coeff_svBtps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::s, filter_3dp::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_rbulk_expand()
{

    std::string prd_name = "pga3dp right bulk expansion";
    auto basis = mv3dp_basis;

    auto lhs = basis;
    // fmt::println("mv_basis for rbulk_expand:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    auto rhs = apply_rules_to_mv(basis, right_bulk_dual_pga3dp_rules);
    // fmt::println("rhs for rbulk_expand:");
    // print_mvec(mv3dp_coeff_svBtps, rhs);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_rweight_expand()
{

    std::string prd_name = "pga3dp right weight expansion";
    auto basis = mv3dp_basis;

    auto lhs = basis;
    // fmt::println("mv_basis for rweight_expand:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    auto rhs = apply_rules_to_mv(basis, right_weight_dual_pga3dp_rules);
    // fmt::println("rhs for rweight_expand:");
    // print_mvec(mv3dp_coeff_svBtps, rhs);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_pga3dp_lbulk_expand()
{

    std::string prd_name = "pga3dp left bulk expansion";
    auto basis = mv3dp_basis;

    auto lhs = apply_rules_to_mv(basis, left_bulk_dual_pga3dp_rules);
    // fmt::println("rhs for lbulk_expand:");
    // print_mvec(mv3dp_coeff_svBtps, rhs);
    // fmt::println("");

    auto rhs = basis;
    // fmt::println("mv_basis for lbulk_expand:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_lweight_expand()
{

    std::string prd_name = "pga3dp left weight expansion";
    auto basis = mv3dp_basis;

    auto lhs = apply_rules_to_mv(basis, left_weight_dual_pga3dp_rules);
    // fmt::println("rhs for lweight_expand:");
    // print_mvec(mv3dp_coeff_svBtps, rhs);
    // fmt::println("");

    auto rhs = basis;
    // fmt::println("mv_basis for lweight_expand:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                        wdg_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_pga3dp_motor_rgpr()
{

    std::string prd_name = "pga3dp regressive sandwich product";
    fmt::println("pga3dp regressive sandwich product - motor * object * rrev(motor):");
    fmt::println("");

    auto basis = mv3dp_basis;

    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv3dp_coeff_svBtps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, rcmpl_pga3dp_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv3dp_coeff_svBtps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
        gpr_pga3dp_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");


    // first product between multivectors in basis_tab (R * v)
    fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u_tmp");
    auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_R_even, mv3dp_coeff_svBtps);
    // fmt::println("{} - product table with coeffs:", prd_name);
    // print_prd_tab(prd_tab);
    // fmt::println("");

    auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e, filter_3dp::vec,
                                        brace_switch::use_braces);
    fmt::println("mv_u_tmp:");
    print_mvec(mv_u_tmp, basis);
    fmt::println("");

    // second product between multivectors for the product v * rrev(R)
    fmt::println("{}:", prd_name + space_str + "mv_u_tmp * rrev(mv_e) -> mv_u_res");
    auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv3dp_coeff_R_rrev_even);
    // fmt::println("prd_tab_v:");
    // print_prd_tab(prd_tab_v);
    // fmt::println("");

    auto mv_u_res_v =
        get_mv_from_prd_tab(prd_tab_v, basis, filter_3dp::mv_u, filter_3dp::mv_e);
    print_mvec(mv_u_res_v, basis);
    fmt::println("");

    //// bivector

    // first product between multivectors in basis_tab (R * B)
    fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e_tmp");
    auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e,
                                        filter_3dp::bivec, brace_switch::use_braces);
    fmt::println("mv_e_tmp:");
    print_mvec(mv_e_tmp, basis);
    fmt::println("");

    // second product between multivectors for the product B * rrev(R)
    fmt::println("{}:", prd_name + space_str + "mv_e_tmp * rrev(mv_e) -> mv_e_res");
    auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv3dp_coeff_R_rrev_even);
    // fmt::println("prd_tab_B:");
    // print_prd_tab(prd_tab_B);
    // fmt::println("");

    auto mv_e_res_B =
        get_mv_from_prd_tab(prd_tab_B, basis, filter_3dp::mv_e, filter_3dp::mv_e);
    print_mvec(mv_e_res_B, basis);
    fmt::println("");

    //// trivector

    // first product between multivectors in basis_tab (R * B)
    fmt::println("{}:", prd_name + space_str + "mv_e * trivec -> mv_u_tmp_t");
    auto mv_u_tmp_t = get_mv_from_prd_tab(prd_tab, basis, filter_3dp::mv_e,
                                          filter_3dp::trivec, brace_switch::use_braces);
    fmt::println("mv_u_tmp_t:");
    print_mvec(mv_u_tmp_t, basis);
    fmt::println("");

    // second product between multivectors for the product B * rrev(R)
    fmt::println("{}:", prd_name + space_str + "mv_u_tmp_t * rrev(mv_e) -> mv_u_res_t");
    auto prd_tab_t = get_prd_tab(basis_tab, mv_u_tmp_t, mv3dp_coeff_R_rrev_even);
    // fmt::println("prd_tab_B:");
    // print_prd_tab(prd_tab_B);
    // fmt::println("");

    auto mv_u_res_t =
        get_mv_from_prd_tab(prd_tab_t, basis, filter_3dp::mv_u, filter_3dp::mv_e);
    print_mvec(mv_u_res_t, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}