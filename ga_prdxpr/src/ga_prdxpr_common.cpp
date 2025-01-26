// author: Daniel Hug, 2024

#include "ga_prdxpr_common.hpp"

#include <exception>

////////////////////////////////////////////////////////////////////////////////
// user related functions
////////////////////////////////////////////////////////////////////////////////

prd_table mv_coeff_to_coeff_prd_tab(mvec_coeff const& lcoeff, mvec_coeff const& rcoeff)
{
    if (lcoeff.size() != rcoeff.size()) {
        throw std::runtime_error("Multivector sizes must match.");
    }

    // reserve a full table of the required size
    prd_table prd_coeff_tab(lcoeff.size(), mvec_coeff(rcoeff.size()));

    // fill-in all products
    for (size_t i = 0; i < lcoeff.size(); ++i) {
        for (size_t j = 0; j < rcoeff.size(); ++j) {
            prd_coeff_tab[i][j] = lcoeff[i] + space_str + mul_str + space_str + rcoeff[j];
        }
    }

    return prd_coeff_tab;
}


prd_table combine_coeff_and_basis_prd_tabs(prd_table const& coeff_tab,
                                           prd_table const& basis_tab)
{

    // make sure sizes match as required
    if (coeff_tab.size() != basis_tab.size()) {
        throw std::runtime_error("Multivector sizes must match.");
    }
    for (size_t i = 0; i < coeff_tab.size(); ++i) {
        if (coeff_tab[i].size() != basis_tab[i].size()) {
            throw std::runtime_error("Multivector sizes must match.");
        }
    }


    // reserve a full table of the required size
    prd_table prd_tab(coeff_tab.size(), mvec_coeff(coeff_tab.size()));


    for (size_t i = 0; i < coeff_tab.size(); ++i) {
        for (size_t j = 0; j < coeff_tab.size(); ++j) {

            if (basis_tab[i][j] == zero_str) {
                prd_tab[i][j] = zero_str;
            }
            else if (basis_tab[i][j].starts_with(minus_str)) {
                prd_tab[i][j] = minus_str + coeff_tab[i][j] + space_str +
                                basis_tab[i][j].substr(1, basis_tab[i][j].size());
            }
            else {
                prd_tab[i][j] = coeff_tab[i][j] + space_str + basis_tab[i][j];
            }
        }
    }

    return prd_tab;
}


prd_table get_prd_tab(prd_table const& basis_tab, mvec_coeff const& mv_basis,
                      mvec_coeff const& mv_lcoeff, mvec_coeff const& mv_rcoeff)
{
    auto mv_coeff_tab = mv_coeff_to_coeff_prd_tab(mv_lcoeff, mv_rcoeff);
    // fmt::println("{} - coeff product table:", prd_name);
    // print_prd_tab(mv_coeff_tab);
    // fmt::println("");

    auto prd_tab = combine_coeff_and_basis_prd_tabs(mv_coeff_tab, basis_tab);
    // fmt::println("{}: product table:", prd_name);
    // print_prd_tab(prd_tab);
    // fmt::println("");

    return prd_tab;
}


mvec_coeff get_mv_from_prd_tab(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                               filter_2d lfilter, filter_2d rfilter)
{

    // make sure sizes match as required
    if (prd_tab.size() != mv_basis.size()) {
        throw std::runtime_error(
            "Multivector size of product table and multivector basis size must match.");
    }

    return extractor(prd_tab, mv_basis, get_coeff_filter(lfilter),
                     get_coeff_filter(rfilter));
}

mvec_coeff get_mv_from_prd_tab(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                               filter_3d lfilter, filter_3d rfilter)
{

    // make sure sizes match as required
    if (prd_tab.size() != mv_basis.size()) {
        throw std::runtime_error(
            "Multivector size of product table and multivector basis size must match.");
    }

    return extractor(prd_tab, mv_basis, get_coeff_filter(lfilter),
                     get_coeff_filter(rfilter));
}


mvec_coeff extractor(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                     mvec_coeff_filter const& lcoeff_filter,
                     mvec_coeff_filter const& rcoeff_filter)
{

    // make sure sizes match as required
    if ((prd_tab.size() != mv_basis.size()) ||
        (mv_basis.size() != lcoeff_filter.size()) ||
        (mv_basis.size() != rcoeff_filter.size())) {
        throw std::runtime_error(
            "Multivector size of product table and multivector basis size must match.");
    }

    // fmt::println("");
    // fmt::println("lcoeff_filter: {}", fmt::join(lcoeff_filter, ", "));
    // fmt::println("rcoeff_filter: {}", fmt::join(rcoeff_filter, ", "));
    // fmt::println("");

    // reserve a full multivector for output
    mvec_coeff mv_prd(mv_basis.size());

    for (size_t b = 0; b < mv_basis.size(); ++b) {
        // for each basis element sum-up all contributions to that basis element

        auto basis_element = mv_basis[b];

        for (size_t i = 0; i < prd_tab.size(); ++i) {
            if (lcoeff_filter[i] == 0) continue; // skip filtered elements on lhs

            for (size_t j = 0; j < prd_tab.size(); ++j) {
                if (rcoeff_filter[j] == 0) continue; // skip filtered elements on rhs

                if (prd_tab[i][j].ends_with(space_str + basis_element)) {
                    // found a contribution
                    if (mv_prd[b] == empty_str) {
                        mv_prd[b] = prd_tab[i][j].substr(0, prd_tab[i][j].size() -
                                                                mv_basis[b].size());
                    }
                    else {

                        if (prd_tab[i][j].starts_with(minus_str)) {
                            mv_prd[b] +=
                                minus_str + space_str +
                                prd_tab[i][j].substr(1, prd_tab[i][j].size() -
                                                            mv_basis[b].size() - 1);
                        }
                        else {
                            mv_prd[b] += plus_str + space_str +
                                         prd_tab[i][j].substr(0, prd_tab[i][j].size() -
                                                                     mv_basis[b].size());
                        }
                    }
                }
            }
        }
    }
    return mv_prd;
}


////////////////////////////////////////////////////////////////////////////////
// printing
////////////////////////////////////////////////////////////////////////////////

void print_mvec(mvec_coeff const& mv, mvec_coeff const& mv_basis)
{

    // make sure sizes match as required
    if (mv.size() != mv_basis.size()) {
        throw std::runtime_error("Size of multivector and basis must match.");
    }

    size_t max_width = 0;
    for (const auto& e : mv) {
        if (e.size() > max_width) max_width = e.size();
    }
    // fmt::println("max_width = {}", max_width);

    // print elements right aligned fitting to max_width
    for (size_t i = 0; i < mv.size(); ++i) {
        const auto& e = mv[i];
        fmt::print("[ ");
        fmt::print("{:>{w}}", e, fmt::arg("w", max_width));
        if (&e != &mv.back()) {
            fmt::println("] {},", mv_basis[i]);
        }
        else {
            fmt::println("] {}", mv_basis[i]);
        }
    }
}

void print_prd_tab(prd_table const& tab)
{
    size_t max_width = 0;
    for (const auto& e : tab) {
        for (const auto& sub : e) {
            if (sub.size() > max_width) max_width = sub.size();
        }
    }
    // fmt::println("max_width = {}", max_width);

    // print elements right aligned fitting to max_width
    for (const auto& e : tab) {
        fmt::print("[ ");
        fmt::print("{:>{w}}", fmt::join(e, ", "), fmt::arg("w", max_width));
        if (&e != &tab.back()) {
            fmt::println("], ");
        }
        else {
            fmt::println("] ");
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
// helper functions (not directly intended for user)
////////////////////////////////////////////////////////////////////////////////
mvec_coeff_filter get_coeff_filter(filter_2d filter)
{

    mvec_coeff_filter filter_vec(4); // 2d multivector has 4 basis elements

    switch (filter) {
        case filter_2d::mv:
            filter_vec = {1, 1, 1, 1}; // all elements
            break;
        case filter_2d::mv_e:
            filter_vec = {1, 0, 0, 1}; // even grade elements
            break;
        case filter_2d::s:
            filter_vec = {1, 0, 0, 0}; // scalar element
            break;
        case filter_2d::vec:
            filter_vec = {0, 1, 1, 0}; // vector elements
            break;
        case filter_2d::ps:
            filter_vec = {0, 0, 0, 1}; // pseudoscalar element
            break;
        default:
            std::unreachable();
    }
    return filter_vec;
}

mvec_coeff_filter get_coeff_filter(filter_3d filter)
{

    mvec_coeff_filter filter_vec(8); // 3d multivector has 8 basis elements

    switch (filter) {
        case filter_3d::mv:
            filter_vec = {1, 1, 1, 1, 1, 1, 1, 1}; // all elements
            break;
        case filter_3d::mv_e:
            filter_vec = {1, 0, 0, 0, 1, 1, 1, 0}; // even grade elements
            break;
        case filter_3d::mv_u:
            filter_vec = {0, 1, 1, 1, 0, 0, 0, 1}; // uneven grade elements
            break;
        case filter_3d::s:
            filter_vec = {1, 0, 0, 0, 0, 0, 0, 0}; // scalar element
            break;
        case filter_3d::vec:
            filter_vec = {0, 1, 1, 1, 0, 0, 0, 0}; // vector elements
            break;
        case filter_3d::bivec:
            filter_vec = {0, 0, 0, 0, 1, 1, 1, 0}; // bivector elements
            break;
        case filter_3d::ps:
            filter_vec = {0, 0, 0, 0, 0, 0, 0, 1}; // pseudoscalar element
            break;
        default:
            std::unreachable();
    }
    return filter_vec;
}