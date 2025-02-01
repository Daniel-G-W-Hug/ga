// author: Daniel Hug, 2024

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
    auto prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * ps -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * vec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv * s -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * ps -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * mv_e -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e * s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * ps -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec * ps -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps * s -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * ps -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * vec -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec * s -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s * vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s * s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
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
    auto prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps ^ mv -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ vec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv ^ s -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e ^ mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps ^ mv_e -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e ^ vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ mv_e -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e ^ s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ vec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec ^ ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps ^ s -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ ps -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ vec -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec ^ s -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s ^ vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s ^ s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2d_lcontr()
{

    // ega2d left contraction
    std::string prd_name = "ega2d left contraction";

    auto basis_tab =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, lcontr_str),
                           lcontr_ega2d_rules);
    auto basis = mv2d_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv << mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv << mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv << ps -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps << mv -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv << vec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv << s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e << mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e << ps -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps << mv_e -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e << vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << mv_e -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e << s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << ps -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << vec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec << ps -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps << s -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << ps -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec << vec -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec << s -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s << vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s << s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2d_rcontr()
{

    // ega2d right contraction
    std::string prd_name = "ega2d right contraction";

    auto basis_tab =
        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, rcontr_str),
                           rcontr_ega2d_rules);
    auto basis = mv2d_basis;

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv >> mv -> mv");
    auto prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> mv_e -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "mv_e >> mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> ps -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps >> mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> vec -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> mv -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv >> s -> mv");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> mv -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> ps -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "ps >> mv_e -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> mv_e -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "mv_e >> s -> mv_e");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> mv_e -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> ps -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> vec -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "vec >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "ps >> s -> ps");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> ps -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec >> vec -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "vec >> s -> vec");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "s >> vec -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "s >> s -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}

void generate_and_print_ega2d_dot()
{

    // ega2d scalar product
    std::string prd_name = "ega2d scalar product";

    auto basis_tab = apply_rules_to_tab(
        mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str), dot_ega2d_rules);
    auto basis = mv2d_basis;

    // fmt::println("equal = {}", basis_tab == dot_ega2d_basis_tab);

    fmt::println("{} - basis product table:", prd_name);
    print_prd_tab(basis_tab);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv,mv) -> s");
    auto prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,mv_e) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(mv_e,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(ps,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(vec,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv,s) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(s,mv) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv_e,mv_e) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv_e,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(ps,mv_e) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv_e,vec) -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(vec,mv_e) -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(mv_e,s) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_A_even, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::mv_e, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(s,mv_e) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_B_even);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::mv_e);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(ps,ps) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(ps,vec) -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(vec,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(ps,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::ps, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(s,ps) -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::ps);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(vec,vec) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(vec,s) -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::vec, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("{}:", prd_name + space_str + "dot(s,vec) -> 0");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps, mv2d_coeff_svps);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::vec);
    print_mvec(prd_mv, basis);
    fmt::println("");


    fmt::println("{}:", prd_name + space_str + "dot(s,s) -> s");
    prd_tab = get_prd_tab(basis_tab, basis, mv2d_coeff_svps1, mv2d_coeff_svps2);
    prd_mv = get_mv_from_prd_tab(prd_tab, basis, filter_2d::s, filter_2d::s);
    print_mvec(prd_mv, basis);
    fmt::println("");

    fmt::println("-------------------------------------------------------------------\n");

    return;
}