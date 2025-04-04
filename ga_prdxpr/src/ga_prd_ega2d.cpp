// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prd_ega2d.hpp"
#include "ga_prdxpr_common.hpp"


void generate_and_print_ega2d_gpr()
{

    // ega2d geometric product
    std::string prd_name = "ega2d geometric product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str), gpr_ega2d_rules);
    auto basis = mv2d_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv * mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * ps -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * vec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * s -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * ps -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv_e -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * ps -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * ps -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s * s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2d_wdg()
{

    // ega2d wedge product
    std::string prd_name = "ega2d wedge product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, wdg_str), wdg_ega2d_rules);
    auto basis = mv2d_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv ^ mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps ^ mv -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ vec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ s -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps ^ mv_e -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e ^ vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ mv_e -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e ^ s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ vec -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s ^ s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2d_lcontract()
{

    std::string prd_name = "ega2d left contraction";

    // (A << B) = rwdg( ldual(A),B ) -> in ega2d = rwdg(lcmpl(A),B)
    //
    // with rwdg(A,B) = lcmpl( wdg( rcmpl(A),rcmpl(B) ) ) in ega2d

    // create the complement from the input multivector (=forward transformation)
    auto basis = mv2d_basis;

    // lhs (dualize with lcmpl as inner function before doing the rcmpl for rwdg)
    auto lhs =
        apply_rules_to_mv(apply_rules_to_mv(basis, lcmpl_ega2d_rules), rcmpl_ega2d_rules);
    // fmt::println("lhs:");
    // print_mvec(mv2d_coeff_svps, lhs);
    // fmt::println("");

    // rhs
    auto rhs = apply_rules_to_mv(basis, rcmpl_ega2d_rules);
    // fmt::println("rhs:");
    // print_mvec(mv2d_coeff_svps, rhs);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega2d_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_ega2d_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv << mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv << mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv << ps -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps << mv -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv << vec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv << s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e << mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e << ps -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps << mv_e -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e << vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << mv_e -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e << s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << ps -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << ps -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << s -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << ps -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec << vec -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec << s -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s << s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2d_rcontract()
{

    std::string prd_name = "ega2d right contraction";

    // (A >> B) = rwdg(A, rdual(B)) -> in ega2d = rwdg(A,rcmpl(B))
    //
    // with rwdg(A,B) = lcmpl( wdg(rcmpl(A),rcmpl(B)) ) in ega2d

    // create the complement from the input multivector (=forward transformation)
    auto basis = mv2d_basis;

    // lhs
    auto lhs = apply_rules_to_mv(basis, rcmpl_ega2d_rules);
    // fmt::println("lhs:");
    // print_mvec(mv2d_coeff_svps, lhs);
    // fmt::println("");

    // lhs (dualize with rcmpl as inner function before doing the rcmpl for rwdg)
    auto rhs =
        apply_rules_to_mv(apply_rules_to_mv(basis, rcmpl_ega2d_rules), rcmpl_ega2d_rules);
    // fmt::println("rhs:");
    // print_mvec(mv2d_coeff_svps, rhs);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega2d_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_ega2d_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv >> mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e >> mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> ps -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps >> mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> vec -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> mv -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> s -> mv");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> mv -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> ps -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps >> mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> mv_e -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_M_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> mv_e -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_M_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> ps -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> vec -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> s -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec >> vec -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec >> s -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> vec -> 0");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s >> s -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_ega2d_lexpand()
{

    std::string prd_name = "ega2d left expansion";
    auto basis = mv2d_basis;

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto lhs = apply_rules_to_mv(basis, lcmpl_ega2d_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2d_coeff_svps, basis_cmpl_func);
    // fmt::println("");

    auto rhs = basis;

    auto basis_tab =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega2d_rules);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "lexpand(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv,s) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(s,mv) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(mv_e,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(s,mv_e) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(ps,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(s,ps) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(vec,s) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "lexpand(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "lexpand(s,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2d_rexpand()
{

    // ega2d right expansion
    std::string prd_name = "ega2d right expansion";
    auto basis = mv2d_basis;

    auto lhs = basis;

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto rhs = apply_rules_to_mv(basis, rcmpl_ega2d_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2d_coeff_svps, basis_cmpl_func);
    // fmt::println("");

    auto basis_tab =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega2d_rules);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    // fmt::println("{}:", prd_name + space_str + "rexpand(mv,mv) -> mv");
    // auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B);
    // auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,mv) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(s,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,mv_e) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(mv_e,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,ps) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,ps) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(ps,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(s,ps) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,vec) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rexpand(s,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rexpand(s,s) -> ps");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}


void generate_and_print_ega2d_dot()
{

    std::string prd_name = "ega2d inner product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str), dot_ega2d_rules);
    auto basis = mv2d_basis;

    // fmt::println("equal = {}", basis_tab == dot_ega2d_basis_tab);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv,mv) -> s");
    auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "dot(mv_e,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(ps,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2d_rwdg()
{

    // Checking the roundtrip lcmpl(rcmpl(mvec)) == mvec
    //
    // fmt::println("mv_basis:");
    // print_mvec(mv2d_coeff_svps, mv2d_basis);
    // fmt::println("");
    // auto double_cmpl_func = apply_rules_to_mv(mv2d_basis, rcmpl_ega2d_rules);
    // fmt::println("double_cmpl 1st:");
    // print_mvec(mv2d_coeff_svps, double_cmpl_func);
    // fmt::println("");
    // double_cmpl_func = apply_rules_to_mv(double_cmpl_func, lcmpl_ega2d_rules);
    // fmt::println("double_cmpl 2nd:");
    // print_mvec(mv2d_coeff_svps, double_cmpl_func);
    // fmt::println("");

    std::string prd_name = "ega2d regressive wedge product";
    auto basis = mv2d_basis;

    // fmt::println("mv_basis for rwdg:");
    // print_mvec(mv2d_coeff_svps, basis);
    // fmt::println("");

    // create the complement from the input multivector
    // (=forward transformation of arguments)
    auto basis_cmpl_func = apply_rules_to_mv(basis, rcmpl_ega2d_rules);
    // fmt::println("cmpl:");
    // print_mvec(mv2d_coeff_svps, basis_cmpl_func);
    // fmt::println("");

    // product between complements
    auto basis_tab_with_rules = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str),
        wdg_ega2d_rules);
    // fmt::println("basis_tab_with_rules:");
    // print_prd_tab(basis_tab_with_rules);
    // fmt::println("");

    // create complements of the table entries (=backward transformation)
    // (returns the product table in terms of the unmodified input multivectors)
    auto basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_ega2d_rules);
    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // create specific products
    fmt::println("{}:", prd_name + space_str + "rwdg(mv,mv) -> mv");
    auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv, mv_e) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,ps) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,vec) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv) -> mv");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,ps) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(ps,mv_e) -> mv_e");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,vec) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,mv_e) -> vec");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(mv_e,s) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_A_even, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,mv_e) -> s");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_B_even);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(ps,ps) -> ps");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(ps,vec) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(vec,ps) -> vec");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(ps,s) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "rwdg(s,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "rwdg(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(vec,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    // fmt::println("{}:", prd_name + space_str + "rwdg(s,vec) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps, mv2d_coeff_svps);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    // print_mvec(prd_mv, basis);
    // fmt::println("");


    // fmt::println("{}:", prd_name + space_str + "rwdg(s,s) -> 0");
    // prd_tab = get_prd_tab(basis_tab, mv2d_coeff_svps1, mv2d_coeff_svps2);
    // prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    // print_mvec(prd_mv, basis);
    // fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2d_rotor_gpr()
{

    std::string prd_name = "ega2d sandwich product";
    fmt::println("ega2d sandwich product: rotor * object * rev(rotor):");
    fmt::println("");

    auto basis = mv2d_basis;
    // fmt::println("mv_basis for sandwich product:");
    // print_mvec(mv2d_coeff_svps, basis);
    // fmt::println("");

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str), gpr_ega2d_rules);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    // first product between multivectors in basis_tab (R * v)
    fmt::println("{}:", prd_name + space_str + "mv_e * vec -> vec_tmp");
    auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_R_even, mv2d_coeff_svps);
    // fmt::println("{} - product table with coeffs:", prd_name);
    // print_prd_tab(prd_tab);
    // fmt::println("");

    auto vec_tmp = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec,
                                       brace_switch::use_braces);
    fmt::println("vec_tmp:");
    print_mvec(vec_tmp, basis);
    fmt::println("");

    // second product between multivectors for the product v * rev(R)
    fmt::println("{}:", prd_name + space_str + "vec_tmp * rev(mv_e) -> vec_res");
    prd_tab = get_prd_tab(basis_tab, vec_tmp, mv2d_coeff_R_rev_even);
    // fmt::println("prd_tab:");
    // print_prd_tab(prd_tab);
    // fmt::println("");

    auto vec_res = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(vec_res, basis);
    fmt::println("");
}