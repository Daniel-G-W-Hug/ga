// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_generator.hpp"
#include <fmt/core.h>
#include <stdexcept>

// Include mathematical definitions
#include "ga_prdxpr_ega2d.hpp"
#include "ga_prdxpr_ega3d.hpp"
#include "ga_prdxpr_pga2dp.hpp"
#include "ga_prdxpr_pga3dp.hpp"

using namespace configurable;

void ConfigurableGenerator::generate_product_expressions(const AlgebraData& algebra,
                                                         const ProductConfig& config)
{
    // No separator logic here - separators are handled at the end of each product

    // Get the basis table using EXISTING mathematical functions
    auto basis_tab = get_basis_table_for_product(algebra, config.product_name);

    // Print basis table if enabled (matches reference format)
    if (config.show_basis_table) {
        print_basis_table(algebra, config, basis_tab);
        fmt::println("");
    }

    // Generate each configured case
    if (config.is_sandwich_product) {
        // For sandwich products, generate the complete two-step process regardless of
        // individual case settings
        OutputCase dummy_case("", "", "", "", "");
        generate_sandwich_case(algebra, config, dummy_case, basis_tab);
    }
    else {
        for (const auto& case_def : config.cases) {
            if (case_def.is_two_step) {
                generate_sandwich_case(algebra, config, case_def, basis_tab);
            }
            else {
                generate_single_case(algebra, config, case_def, basis_tab);
            }
        }
    }
}

void ConfigurableGenerator::generate_single_case(const AlgebraData& algebra,
                                                 const ProductConfig& config,
                                                 const OutputCase& case_def,
                                                 const prd_table& basis_tab)
{
    // Get coefficients from algebra data (maps to existing coefficient objects)
    auto left_coeff_it = algebra.coefficients.find(case_def.left_coeff_name);
    auto right_coeff_it = algebra.coefficients.find(case_def.right_coeff_name);

    if (left_coeff_it == algebra.coefficients.end()) {
        throw std::invalid_argument("Unknown left coefficient: " +
                                    case_def.left_coeff_name);
    }
    if (right_coeff_it == algebra.coefficients.end()) {
        throw std::invalid_argument("Unknown right coefficient: " +
                                    case_def.right_coeff_name);
    }

    const auto& left_coeff = left_coeff_it->second;
    const auto& right_coeff = right_coeff_it->second;

    // Get filters based on algebra dimension
    mvec_coeff prd_mv;
    auto prd_tab = get_prd_tab(basis_tab, left_coeff, right_coeff);

    if (algebra.dimension == 2) {
        auto left_filter = get_filter_2d(algebra, case_def.left_filter_name);
        auto right_filter = get_filter_2d(algebra, case_def.right_filter_name);
        prd_mv = get_mv_from_prd_tab(prd_tab, algebra.basis, left_filter, right_filter);
    }
    else if (algebra.dimension == 3) {
        auto left_filter = get_filter_3d(algebra, case_def.left_filter_name);
        auto right_filter = get_filter_3d(algebra, case_def.right_filter_name);
        prd_mv = get_mv_from_prd_tab(prd_tab, algebra.basis, left_filter, right_filter);
    }
    else if (algebra.dimension == 4) {
        auto left_filter = get_filter_4d(algebra, case_def.left_filter_name);
        auto right_filter = get_filter_4d(algebra, case_def.right_filter_name);
        prd_mv = get_mv_from_prd_tab(prd_tab, algebra.basis, left_filter, right_filter);
    }
    else {
        throw std::invalid_argument("Unsupported algebra dimension: " +
                                    std::to_string(algebra.dimension));
    }

    // Format output to match reference implementation exactly
    print_case_header(algebra, config, case_def.case_name);
    print_case_result(prd_mv, algebra.basis);
    fmt::println("");
}

void ConfigurableGenerator::generate_sandwich_case(const AlgebraData& algebra,
                                                   const ProductConfig& config,
                                                   const OutputCase& /* case_def */,
                                                   const prd_table& basis_tab)
{
    std::string prd_name = algebra.name + " " + config.display_name;

    if (algebra.name == "ega2d") {

        // First product between multivectors in basis_tab (R * v)
        fmt::println("{}:", prd_name + space_str + "mv_e * vec -> vec_tmp");
        auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_R_even, mv2d_coeff_svps);

        auto vec_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_2d::mv_e,
                                           filter_2d::vec, brace_switch::use_braces);
        fmt::println("vec_tmp:");
        print_mvec(vec_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product v * rev(R)
        fmt::println("{}:", prd_name + space_str + "vec_tmp * rev(mv_e) -> vec_res");
        prd_tab = get_prd_tab(basis_tab, vec_tmp, mv2d_coeff_R_rev_even);

        auto vec_res =
            get_mv_from_prd_tab(prd_tab, algebra.basis, filter_2d::vec, filter_2d::mv_e);
        print_mvec(vec_res, algebra.basis);
        fmt::println("");
    }

    else if (algebra.name == "ega3d") {

        // First product between multivectors in basis_tab (R * v)
        fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u_tmp");
        auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_R_even, mv3d_coeff_svBps);

        auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_3d::mv_e,
                                            filter_3d::vec, brace_switch::use_braces);
        fmt::println("mv_u_tmp:");
        print_mvec(mv_u_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product mv_u_tmp * rev(R)
        fmt::println("{}:", prd_name + space_str + "mv_u_tmp * rev(mv_e) -> mv_u_res");
        auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv3d_coeff_R_rev_even);

        auto mv_u_res_v = get_mv_from_prd_tab(prd_tab_v, algebra.basis, filter_3d::mv_u,
                                              filter_3d::mv_e);
        print_mvec(mv_u_res_v, algebra.basis);
        fmt::println("");

        //// Bivector case

        // First product between multivectors in basis_tab (R * B)
        fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e_tmp");
        auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_3d::mv_e,
                                            filter_3d::bivec, brace_switch::use_braces);
        fmt::println("mv_e_tmp:");
        print_mvec(mv_e_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product mv_e_tmp * rev(R)
        fmt::println("{}:", prd_name + space_str + "mv_e_tmp * rev(mv_e) -> mv_e_res");
        auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv3d_coeff_R_rev_even);

        auto mv_e_res_B = get_mv_from_prd_tab(prd_tab_B, algebra.basis, filter_3d::mv_e,
                                              filter_3d::mv_e);
        print_mvec(mv_e_res_B, algebra.basis);
        fmt::println("");
    }

    else if (algebra.name == "pga2dp") {

        // Vector case: first product between multivectors in basis_tab (R * v)
        fmt::println("{}:", prd_name + space_str + "rgpr(mv_u, vec) -> mv_u_tmp");
        auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_R_odd, mv2dp_coeff_svBps);

        auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_3d::mv_u,
                                            filter_3d::vec, brace_switch::use_braces);
        fmt::println("mv_u_tmp:");
        print_mvec(mv_u_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product mv_u_tmp * rrev(R)
        fmt::println("{}:", prd_name + space_str + "mv_u_tmp * rrev(mv_u) -> mv_u_res");
        auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv2dp_coeff_R_rrev_odd);

        auto mv_u_res_v = get_mv_from_prd_tab(prd_tab_v, algebra.basis, filter_3d::mv_u,
                                              filter_3d::mv_u);
        print_mvec(mv_u_res_v, algebra.basis);
        fmt::println("");

        //// Bivector case

        // First product between multivectors in basis_tab (R * B)
        fmt::println("{}:", prd_name + space_str + "mv_u * bivec -> mv_e_tmp");
        auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_3d::mv_u,
                                            filter_3d::bivec, brace_switch::use_braces);
        fmt::println("mv_e_tmp:");
        print_mvec(mv_e_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product mv_e_tmp * rrev(R)
        fmt::println("{}:", prd_name + space_str + "mv_e_tmp * rrev(mv_u) -> mv_e_res");
        auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv2dp_coeff_R_rrev_odd);

        auto mv_e_res_B = get_mv_from_prd_tab(prd_tab_B, algebra.basis, filter_3d::mv_e,
                                              filter_3d::mv_u);
        print_mvec(mv_e_res_B, algebra.basis);
        fmt::println("");
    }

    else if (algebra.name == "pga3dp") {

        // Vector case: first product between multivectors in basis_tab (R * v)
        fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u_tmp");
        auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_R_even, mv3dp_coeff_svBtps);

        auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_4d::mv_e,
                                            filter_4d::vec, brace_switch::use_braces);
        fmt::println("mv_u_tmp:");
        print_mvec(mv_u_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product mv_u_tmp * rrev(R)
        fmt::println("{}:", prd_name + space_str + "mv_u_tmp * rrev(mv_e) -> mv_u_res");
        auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv3dp_coeff_R_rrev_even);

        auto mv_u_res_v = get_mv_from_prd_tab(prd_tab_v, algebra.basis, filter_4d::mv_u,
                                              filter_4d::mv_e);
        print_mvec(mv_u_res_v, algebra.basis);
        fmt::println("");

        //// Bivector case

        // First product between multivectors in basis_tab (R * B)
        fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e_tmp");
        auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_4d::mv_e,
                                            filter_4d::bivec, brace_switch::use_braces);
        fmt::println("mv_e_tmp:");
        print_mvec(mv_e_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product mv_e_tmp * rrev(R)
        fmt::println("{}:", prd_name + space_str + "mv_e_tmp * rrev(mv_e) -> mv_e_res");
        auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv3dp_coeff_R_rrev_even);

        auto mv_e_res_B = get_mv_from_prd_tab(prd_tab_B, algebra.basis, filter_4d::mv_e,
                                              filter_4d::mv_e);
        print_mvec(mv_e_res_B, algebra.basis);
        fmt::println("");

        //// Trivector case

        // First product between multivectors in basis_tab (R * T)
        fmt::println("{}:", prd_name + space_str + "mv_e * trivec -> mv_u_tmp_t");
        auto mv_u_tmp_t =
            get_mv_from_prd_tab(prd_tab, algebra.basis, filter_4d::mv_e,
                                filter_4d::trivec, brace_switch::use_braces);
        fmt::println("mv_u_tmp_t:");
        print_mvec(mv_u_tmp_t, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product mv_u_tmp_t * rrev(R)
        fmt::println("{}:",
                     prd_name + space_str + "mv_u_tmp_t * rrev(mv_e) -> mv_u_res_t");
        auto prd_tab_t = get_prd_tab(basis_tab, mv_u_tmp_t, mv3dp_coeff_R_rrev_even);

        auto mv_u_res_t = get_mv_from_prd_tab(prd_tab_t, algebra.basis, filter_4d::mv_u,
                                              filter_4d::mv_e);
        print_mvec(mv_u_res_t, algebra.basis);
        fmt::println("");
    }
}

prd_table
ConfigurableGenerator::get_basis_table_for_product(const AlgebraData& algebra,
                                                   const std::string& product_name)
{
    // This function maps product names to the EXISTING mathematical basis table
    // generation from the reference implementation - no changes to mathematical logic

    if (algebra.name == "ega2d") {

        if (product_name == "gpr") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str),
                gpr_ega2d_rules);
        }

        else if (product_name == "cmt") {
            // Commutator product (=asymmetric part of the geometric product)
            auto basis_tab = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str),
                gpr_ega2d_rules);
            return get_prd_tab_asym(basis_tab); // use the asymmetric part only
        }

        else if (product_name == "wdg") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, wdg_str),
                wdg_ega2d_rules);
        }

        else if (product_name == "dot") {
            // Inner product (=dot product)
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str),
                dot_ega2d_rules);
        }

        else if (product_name == "dot (alternative)") {
            // Inner product (=dot product) - defined by outer product
            // (see Grassmann Algebra, John Browne, p. 35, p. 373)
            // remember: A ^ rcmpl(A) = I_2d
            // inner product = rwdg(A,rcmpl(A))
            //               = lcmpl( wdg( rcmpl(A), rcmpl( rcmpl(A) ) ) )
            // => only valid, if both operands are of the same grade, i.e. only the
            // main diagonal is valid

            auto cmpl_func = apply_rules_to_mv(mv2d_basis, rcmpl_ega2d_rules);
            auto dbl_cmpl_func = apply_rules_to_mv(cmpl_func, rcmpl_ega2d_rules);
            // fmt::println("");
            // fmt::println("mv2d_basis = {}", mv2d_basis);
            // fmt::println("cmpl_func = {}", cmpl_func);
            // fmt::println("dbl_cmpl_func = {}", dbl_cmpl_func);
            // fmt::println("");
            auto tmp_tab = mv_coeff_to_coeff_prd_tab(cmpl_func, dbl_cmpl_func, wdg_str);
            // fmt::println("tmp_tab = {}", tmp_tab);
            // fmt::println("");
            auto basis_tab_with_rules = apply_rules_to_tab(tmp_tab, wdg_ega2d_rules);
            // fmt::println("basis_tab_with_rules = {}", basis_tab_with_rules);
            // fmt::println("");
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_ega2d_rules);
        }

        else if (product_name == "rwdg") {
            // Regressive wedge: rwdg(A,B) = lcmpl(wdg(rcmpl(A), rcmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2d_basis, rcmpl_ega2d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str),
                wdg_ega2d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_ega2d_rules);
        }

        else if (product_name == "lcontract") {
            // A << B = rwdg(lcmpl(A), B)
            //        = lcmpl( wdg( rcmpl(lcmpl(A)), rcmpl(B) ) )
            auto lhs = apply_rules_to_mv(apply_rules_to_mv(mv2d_basis, lcmpl_ega2d_rules),
                                         rcmpl_ega2d_rules);
            auto rhs = apply_rules_to_mv(mv2d_basis, rcmpl_ega2d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega2d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_ega2d_rules);
        }

        else if (product_name == "rcontract") {
            // A >> B = rwdg(A,rcmpl(B))
            //        = lcmpl( wdg( rcmpl(A),rcmpl(rcmpl(B)) )  )
            auto lhs = apply_rules_to_mv(mv2d_basis, rcmpl_ega2d_rules);
            auto rhs = apply_rules_to_mv(apply_rules_to_mv(mv2d_basis, rcmpl_ega2d_rules),
                                         rcmpl_ega2d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega2d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_ega2d_rules);
        }

        else if (product_name == "lexpand") {
            // Left expansion: lexpand(A,B) = wdg(lcmpl(A), B)
            auto lhs = apply_rules_to_mv(mv2d_basis, lcmpl_ega2d_rules);
            auto rhs = mv2d_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_ega2d_rules);
        }

        else if (product_name == "rexpand") {
            // Right expansion: rexpand(A,B) = wdg(A, rcmpl(B))
            auto lhs = mv2d_basis;
            auto rhs = apply_rules_to_mv(mv2d_basis, rcmpl_ega2d_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_ega2d_rules);
        }

        else if (product_name == "sandwich_gpr") {
            // Rotor sandwich product basis table
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str),
                gpr_ega2d_rules);
        }
    }

    else if (algebra.name == "ega3d") {

        if (product_name == "gpr") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str),
                gpr_ega3d_rules);
        }

        else if (product_name == "cmt") {
            // Commutator product (=asymmetric part of the geometric product)
            auto basis_tab = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str),
                gpr_ega3d_rules);
            return get_prd_tab_asym(basis_tab); // use the asymmetric part only
        }

        else if (product_name == "wdg") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, wdg_str),
                wdg_ega3d_rules);
        }

        else if (product_name == "dot") {
            // Inner product (=dot product)
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str),
                dot_ega3d_rules);
        }

        else if (product_name == "rwdg") {
            // Regressive wedge: rwdg(A,B) = lcmpl(wdg(rcmpl(A), rcmpl(B)))
            // For EGA3D: lcmpl = rcmpl = cmpl
            auto basis_cmpl_func = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str),
                wdg_ega3d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
        }

        else if (product_name == "lcontract") {
            // Left contraction: A << B = rwdg(lcmpl(A), B)
            // For EGA3D: lcmpl = cmpl
            //                          = cmpl( wdg( cmpl(cmpl(A)), cmpl(B) ) )
            auto lhs = apply_rules_to_mv(apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules),
                                         cmpl_ega3d_rules);
            auto rhs = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega3d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
        }

        else if (product_name == "rcontract") {
            // Right contraction: A >> B = rwdg(A, rcmpl(B)))
            // For EGA3D: lcmpl = rcmpl = cmpl
            //                           = cmpl( wdg( cmpl(A),cmpl(cmpl(B)) )  )
            auto lhs = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto rhs = apply_rules_to_mv(apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules),
                                         cmpl_ega3d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_ega3d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
        }

        else if (product_name == "lexpand") {
            // Left expansion: lexpand(A,B) = wdg(lcmpl(A), B)
            // For EGA3D: lcmpl = cmpl
            auto lhs = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto rhs = mv3d_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_ega3d_rules);
        }

        else if (product_name == "rexpand") {
            // Right expansion: rexpand(A,B) = wdg(A, rcmpl(B))
            // For EGA3D: rcmpl = cmpl
            auto lhs = mv3d_basis;
            auto rhs = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_ega3d_rules);
        }

        else if (product_name == "sandwich_gpr") {
            // Rotor sandwich product basis table
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str),
                gpr_ega3d_rules);
        }
    }

    else if (algebra.name == "pga2dp") {

        if (product_name == "gpr") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str),
                gpr_pga2dp_rules);
        }

        else if (product_name == "cmt") {
            // Commutator product (=asymmetric part of the geometric product)
            //                   cmt(A,B) = asym(gpr(A,B))
            auto basis_tab = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str),
                gpr_pga2dp_rules);
            return get_prd_tab_asym(basis_tab); // use the asymmetric part only
        }

        else if (product_name == "wdg") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, wdg_str),
                wdg_pga2dp_rules);
        }

        else if (product_name == "dot") {
            // Inner product (=dot product)
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str),
                dot_pga2dp_rules);
        }

        else if (product_name == "rwdg") {
            // Regressive wedge: rwdg(A,B) = cmpl(wdg(cmpl(A), cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str),
                wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "rdot") {
            // Regressive inner product: rdot(A,B) = cmpl(dot(cmpl(A), cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
                dot_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "rgpr") {
            // Regressive geometric product: rgpr(A,B) = cmpl(gpr(cmpl(A), cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
                gpr_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "rcmt") {
            // Commutator product: cmt(A,B) = asym(gpr(A,B))
            // Regressive commutator product:
            // rcmt(A,B) = asym(cmpl(gpr(cmpl(A),cmpl(B))))
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
                gpr_pga2dp_rules);
            auto full_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
            return get_prd_tab_asym(full_tab);
        }

        else if (product_name == "right_bulk_contract") {
            // Right bulk contraction: A >> B = cmpl(wdg(cmpl(A), cmpl(bulk_dual(B))))
            auto lhs = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto rhs =
                apply_rules_to_mv(apply_rules_to_mv(mv2dp_basis, bulk_dual_pga2dp_rules),
                                  cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "right_weight_contract") {
            // Right weight contraction: A >> B = cmpl(wdg(cmpl(A), cmpl(weight_dual(B))))
            auto lhs = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto rhs = apply_rules_to_mv(
                apply_rules_to_mv(mv2dp_basis, weight_dual_pga2dp_rules),
                cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "left_bulk_contract") {
            // Left bulk contraction: A << B = cmpl(wdg(cmpl(bulk_dual(A)), cmpl(B)))
            auto lhs =
                apply_rules_to_mv(apply_rules_to_mv(mv2dp_basis, bulk_dual_pga2dp_rules),
                                  cmpl_pga2dp_rules);
            auto rhs = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "left_weight_contract") {
            // Left weight contraction: A << B = cmpl(wdg(cmpl(weight_dual(A)), cmpl(B)))
            auto lhs = apply_rules_to_mv(
                apply_rules_to_mv(mv2dp_basis, weight_dual_pga2dp_rules),
                cmpl_pga2dp_rules);
            auto rhs = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "right_bulk_expand") {
            // Right bulk expansion: A >< B = wdg(A, bulk_dual(B))
            auto lhs = mv2dp_basis;
            auto rhs = apply_rules_to_mv(mv2dp_basis, bulk_dual_pga2dp_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_pga2dp_rules);
        }

        else if (product_name == "right_weight_expand") {
            // Right weight expansion: A >< B = wdg(A, weight_dual(B))
            auto lhs = mv2dp_basis;
            auto rhs = apply_rules_to_mv(mv2dp_basis, weight_dual_pga2dp_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_pga2dp_rules);
        }

        else if (product_name == "left_bulk_expand") {
            // Left bulk expansion: A <> B = wdg(bulk_dual(A), B)
            auto lhs = apply_rules_to_mv(mv2dp_basis, bulk_dual_pga2dp_rules);
            auto rhs = mv2dp_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_pga2dp_rules);
        }

        else if (product_name == "left_weight_expand") {
            // Left weight expansion: A <> B = wdg(weight_dual(A), B)
            auto lhs = apply_rules_to_mv(mv2dp_basis, weight_dual_pga2dp_rules);
            auto rhs = mv2dp_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_pga2dp_rules);
        }

        else if (product_name == "sandwich_rgpr") {
            // Regressive sandwich product basis table - same as rgpr for PGA2DP
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
                gpr_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }
    }

    else if (algebra.name == "pga3dp") {

        if (product_name == "gpr") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, mul_str),
                gpr_pga3dp_rules);
        }

        else if (product_name == "cmt") {
            // Commutator product (=asymmetric part of the geometric product)
            //                   cmt(A,B) = asym(gpr(A,B))
            auto basis_tab = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, mul_str),
                gpr_pga3dp_rules);
            return get_prd_tab_asym(basis_tab); // use the asymmetric part only
        }

        else if (product_name == "wdg") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, wdg_str),
                wdg_pga3dp_rules);
        }

        else if (product_name == "dot") {
            // Inner product (=dot product)
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, mul_str),
                dot_pga3dp_rules);
        }

        else if (product_name == "rwdg") {
            // Regressive wedge: rwdg(A,B) = lcmpl(wdg(rcmpl(A), rcmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str),
                wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
        }

        else if (product_name == "rdot") {
            // Regressive inner: rdot(A,B) = lcmpl(dot(rcmpl(A), rcmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
                dot_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
        }

        else if (product_name == "rgpr") {
            // Regressive geometric: rgpr(A,B) = lcmpl(gpr(rcmpl(A), rcmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
                gpr_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
        }

        else if (product_name == "rcmt") {
            // Commutator product: cmt(A,B) = asym(gpr(A,B))
            // Regressive commutator product:
            // rcmt(A,B) = asym(lcmpl(gpr(rcmpl(A),rcmpl(B))))
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
                gpr_pga3dp_rules);
            auto full_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
            return get_prd_tab_asym(full_tab);
        }

        else if (product_name == "right_bulk_contract") {
            // Right bulk contraction:
            // right_bulk_contract(A,B) = lcmpl(wdg(rcmpl(A), rcmpl(right_bulk_dual(B))))
            auto lhs = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
            auto rhs = apply_rules_to_mv(
                apply_rules_to_mv(mv3dp_basis, right_bulk_dual_pga3dp_rules),
                rcmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
        }

        else if (product_name == "right_weight_contract") {
            // Right weight contraction: right_weight_contract(A,B) =
            // lcmpl(wdg(rcmpl(A), rcmpl(right_weight_dual(B))))
            auto lhs = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
            auto rhs = apply_rules_to_mv(
                apply_rules_to_mv(mv3dp_basis, right_weight_dual_pga3dp_rules),
                rcmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
        }

        else if (product_name == "left_bulk_contract") {
            // Left bulk contraction: left_bulk_contract(A,B) =
            // lcmpl(wdg(rcmpl(left_bulk_dual(A)), rcmpl(B)))
            auto lhs = apply_rules_to_mv(
                apply_rules_to_mv(mv3dp_basis, left_bulk_dual_pga3dp_rules),
                rcmpl_pga3dp_rules);
            auto rhs = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
        }

        else if (product_name == "left_weight_contract") {
            // Left weight contraction: left_weight_contract(A,B) =
            // lcmpl(wdg(rcmpl(left_weight_dual(A)), rcmpl(B)))
            auto lhs = apply_rules_to_mv(
                apply_rules_to_mv(mv3dp_basis, left_weight_dual_pga3dp_rules),
                rcmpl_pga3dp_rules);
            auto rhs = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
        }

        else if (product_name == "right_bulk_expand") {
            // Right bulk expansion: right_bulk_expand(A,B) = wdg(A, right_bulk_dual(B))
            auto lhs = mv3dp_basis;
            auto rhs = apply_rules_to_mv(mv3dp_basis, right_bulk_dual_pga3dp_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_pga3dp_rules);
        }

        else if (product_name == "right_weight_expand") {
            // Right weight expansion: right_weight_expand(A,B) = wdg(A,
            // right_weight_dual(B))
            auto lhs = mv3dp_basis;
            auto rhs = apply_rules_to_mv(mv3dp_basis, right_weight_dual_pga3dp_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_pga3dp_rules);
        }

        else if (product_name == "left_bulk_expand") {
            // Left bulk expansion: left_bulk_expand(A,B) = wdg(left_bulk_dual(A), B)
            auto lhs = apply_rules_to_mv(mv3dp_basis, left_bulk_dual_pga3dp_rules);
            auto rhs = mv3dp_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_pga3dp_rules);
        }

        else if (product_name == "left_weight_expand") {
            // Left weight expansion:
            // left_weight_expand(A,B) = wdg(left_weight_dual(A), B)
            auto lhs = apply_rules_to_mv(mv3dp_basis, left_weight_dual_pga3dp_rules);
            auto rhs = mv3dp_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str),
                                      wdg_pga3dp_rules);
        }

        else if (product_name == "sandwich_rgpr") {
            // Regressive sandwich product basis table - same as rgpr for PGA3DP
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, rcmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str),
                gpr_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, lcmpl_pga3dp_rules);
        }
    }


    // Add other algebras as needed...
    // For now, throw an error for unimplemented combinations
    throw std::invalid_argument("Unsupported product: " + algebra.name +
                                "::" + product_name);
}

filter_2d ConfigurableGenerator::get_filter_2d(const AlgebraData& algebra,
                                               const std::string& filter_name)
{
    auto it = algebra.filters_2d.find(filter_name);
    if (it == algebra.filters_2d.end()) {
        throw std::invalid_argument("Unknown 2D filter: " + filter_name);
    }
    return it->second;
}

filter_3d ConfigurableGenerator::get_filter_3d(const AlgebraData& algebra,
                                               const std::string& filter_name)
{
    auto it = algebra.filters_3d.find(filter_name);
    if (it == algebra.filters_3d.end()) {
        throw std::invalid_argument("Unknown 3D filter: " + filter_name);
    }
    return it->second;
}

filter_4d ConfigurableGenerator::get_filter_4d(const AlgebraData& algebra,
                                               const std::string& filter_name)
{
    auto it = algebra.filters_4d.find(filter_name);
    if (it == algebra.filters_4d.end()) {
        throw std::invalid_argument("Unknown 4D filter: " + filter_name);
    }
    return it->second;
}

void ConfigurableGenerator::print_product_header(const AlgebraData& algebra,
                                                 const ProductConfig& config)
{
    fmt::println("{} {} - basis product table:", algebra.name, config.display_name);
}

void ConfigurableGenerator::print_basis_table(const AlgebraData& algebra,
                                              const ProductConfig& config,
                                              const prd_table& basis_tab)
{
    print_product_header(algebra, config);
    // Use EXISTING print function from reference implementation
    print_prd_tab(basis_tab);
    fmt::println("\nsymmetric part:");
    print_prd_tab(get_prd_tab_sym(basis_tab));
    fmt::println("\nasymmetric part:");
    print_prd_tab(get_prd_tab_asym(basis_tab));
}

void ConfigurableGenerator::print_case_header(const AlgebraData& algebra,
                                              const ProductConfig& config,
                                              const std::string& case_name)
{
    // reference format, e.g. "ega2d geometric product mv * mv -> mv:"
    fmt::println("{} {} {}:", algebra.name, config.display_name, case_name);
}

void ConfigurableGenerator::print_case_result(const mvec_coeff& result,
                                              const mvec_coeff& basis)
{
    // Use EXISTING print function from reference implementation
    print_mvec(result, basis);
}