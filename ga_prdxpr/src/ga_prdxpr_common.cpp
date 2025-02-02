// author: Daniel Hug, 2024

#include "ga_prdxpr_common.hpp"

#include <exception>

////////////////////////////////////////////////////////////////////////////////
// user related functions
////////////////////////////////////////////////////////////////////////////////

prd_table mv_coeff_to_coeff_prd_tab(mvec_coeff const& lcoeff, mvec_coeff const& rcoeff,
                                    std::string const& operator_str)
{
    if (lcoeff.size() != rcoeff.size()) {
        throw std::runtime_error("Multivector sizes must match.");
    }

    // check rules: there must be not leading or trailing space_str
    // (check for potential user error, when manually entering the coefficients)
    for (auto const& e : lcoeff) {
        if (e.starts_with(space_str) || e.ends_with(space_str)) {
            fmt::println("mvec: {}", lcoeff);
            fmt::println("mvec coefficient: '{}'", e);
            throw std::runtime_error("Product coefficients on left hand side must not "
                                     "start or end with space_str.");
        }
    }
    for (auto const& e : rcoeff) {
        if (e.starts_with(space_str) || e.ends_with(space_str)) {
            fmt::println("mvec: {}", rcoeff);
            fmt::println("mvec coefficient: '{}'", e);
            throw std::runtime_error("Product coefficients on right hand side must not "
                                     "start or end with space_str.");
        }
    }

    // reserve a full table of the required size
    prd_table prd_coeff_tab(lcoeff.size(), mvec_coeff(rcoeff.size()));

    // fill-in all products
    for (size_t i = 0; i < lcoeff.size(); ++i) {
        for (size_t j = 0; j < rcoeff.size(); ++j) {

            bool is_negative{false}; // must be "false" for toggle logic below to work
            std::string lhs{lcoeff[i]};
            std::string rhs{rcoeff[j]};

            if (lcoeff[i].starts_with(minus_str)) {
                toggle_bool(is_negative);
                lhs = lcoeff[i].substr(1, lcoeff[i].size() - 1);
            }

            if (rcoeff[j].starts_with(minus_str)) {
                toggle_bool(is_negative);
                rhs = rcoeff[j].substr(1, rcoeff[j].size() - 1);
            }

            if (is_negative) {
                prd_coeff_tab[i][j] =
                    minus_str + lhs + space_str + operator_str + space_str + rhs;
            }
            else {
                prd_coeff_tab[i][j] = lhs + space_str + operator_str + space_str + rhs;
            }
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

mvec_coeff apply_rules_to_mv(mvec_coeff const& coeff, mvec_rules const& rules)
{

    // check rules: there must be not leading or trailing space_str
    // (check for potential user error, when defininig the rules manually entered)
    for (auto const& r : rules) {
        if (r.first.starts_with(space_str) || r.first.ends_with(space_str) ||
            r.second.starts_with(space_str) || r.second.ends_with(space_str)) {
            fmt::println("Rule key: '{}', rule value: '{}'", r.first, r.second);
            throw std::runtime_error("Rules must not start or end with space_str.");
        }
    }

    mvec_coeff mvec(coeff.size());

    for (size_t i = 0; i < coeff.size(); ++i) {

        bool is_negative{false}; // must be "false" for toggle logic below to work
        std::string value{coeff[i]};

        // remove the minus-sign, if it is present in the input
        if (coeff[i].starts_with(minus_str)) {
            toggle_bool(is_negative);
            value = coeff[i].substr(1, coeff[i].size() - 1);
        }

        // fmt::println("i: {}, value: '{}', is_negative: {}", i, value, is_negative);

        // apply the rule to the input, after the input is without intial minus_str
        // zero_str value will be left unchanged, i.e. implicit rule "0" -> "0"
        if (value != zero_str) {
            value = rules.at(value);
        }

        // fmt::println("i: {}, value: '{}', is_negative: {}", i, value, is_negative);

        // in case the replaced value does already have a minus sign on its own
        if (value.starts_with(minus_str)) {
            toggle_bool(is_negative);
            value = value.substr(1, value.size() - 1);
        }

        // fmt::println("i: {}, value: '{}', is_negative: {}", i, value, is_negative);

        if (is_negative && (value != zero_str)) { // zero doesn't need a minus sign
            mvec[i] = minus_str + value;
        }
        else {
            mvec[i] = value;
        }
    }

    return mvec;
}

prd_table apply_rules_to_tab(prd_table const& tab, prd_rules const& rules)
{

    for (size_t i = 0; i < tab.size(); ++i) {
        if (tab[i].size() != tab.size()) {
            throw std::runtime_error("Product tables must be square matrices. Sizes of "
                                     "rows and columns must match.");
        }
    }

    // check rules: there must be not leading or trailing space_str
    // (check for potential user error, when defininig the rules manually entered)
    for (auto const& r : rules) {
        if (r.first.starts_with(space_str) || r.first.ends_with(space_str) ||
            r.second.starts_with(space_str) || r.second.ends_with(space_str)) {
            fmt::println("Rule key: '{}', rule value: '{}'", r.first, r.second);
            throw std::runtime_error("Rules must not start or end with empty_str.");
        }
    }

    // reserve a full table of the required size
    prd_table prd_tab(tab.size(), mvec_coeff(tab.size()));

    for (size_t i = 0; i < tab.size(); ++i) {
        for (size_t j = 0; j < tab.size(); ++j) {

            bool is_negative{false}; // must be "false" for toggle logic below to work
            std::string value{tab[i][j]};

            // remove the minus-sign, if it is present in the input
            if (tab[i][j].starts_with(minus_str)) {
                toggle_bool(is_negative);
                value = tab[i][j].substr(1, tab[i][j].size() - 1);
            }

            // apply the rule to the input, after the input is without intial minus_str
            // zero_str value will be left unchanged, i.e. implicit rule "0" -> "0"
            if (value != zero_str) {
                value = rules.at(value);
            }

            // fmt::println("i: {}, j: {}, value: '{}', is_negative: {}", i, j, value,
            //              is_negative);

            // in case the replaced value does already have a minus sign on its own
            if (value.starts_with(minus_str)) {
                toggle_bool(is_negative);
                value = value.substr(1, value.size() - 1);
            }

            // fmt::println("i: {}, j: {}, value: '{}', is_negative: {}", i, j, value,
            //              is_negative);

            if (is_negative && (value != zero_str)) { // zero doesn't need a minus sign
                prd_tab[i][j] = minus_str + value;
            }
            else {
                prd_tab[i][j] = value;
            }
        }
    }

    return prd_tab;
}


prd_table get_prd_tab(prd_table const& basis_tab, mvec_coeff const& mv_lcoeff,
                      mvec_coeff const& mv_rcoeff)
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
        throw std::runtime_error("Multivector size of product table and multivector "
                                 "basis size must match.");
    }

    return extractor(prd_tab, mv_basis, get_coeff_filter(lfilter),
                     get_coeff_filter(rfilter));
}

mvec_coeff get_mv_from_prd_tab(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                               filter_3d lfilter, filter_3d rfilter)
{

    // make sure sizes match as required
    if (prd_tab.size() != mv_basis.size()) {
        throw std::runtime_error("Multivector size of product table and multivector "
                                 "basis size must match.");
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
        throw std::runtime_error("Multivector size of product table and multivector "
                                 "basis size must match.");
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
                        // add the first entry (including a minus-sign, if present)
                        // remove empty_str and basis element from prd_tab entry
                        mv_prd[b] = prd_tab[i][j].substr(0, prd_tab[i][j].size() -
                                                                mv_basis[b].size() - 1);
                    }
                    else {
                        // add the 2nd and following entries
                        if (prd_tab[i][j].starts_with(minus_str)) {
                            // remove empty_str and basis element and minus sign from
                            // prd_tab entry
                            mv_prd[b] +=
                                space_str + minus_str + space_str +
                                prd_tab[i][j].substr(1, prd_tab[i][j].size() -
                                                            mv_basis[b].size() - 2);
                        }
                        else {
                            // remove empty_str and basis element from prd_tab entry
                            mv_prd[b] +=
                                space_str + plus_str + space_str +
                                prd_tab[i][j].substr(0, prd_tab[i][j].size() -
                                                            mv_basis[b].size() - 1);
                            ;
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
            fmt::println(" ] {},", mv_basis[i]);
        }
        else {
            fmt::println(" ] {}", mv_basis[i]);
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
            fmt::println(" ], ");
        }
        else {
            fmt::println(" ] ");
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


void toggle_bool(bool& truth_value)
{
    if (truth_value) {
        truth_value = false;
    }
    else {
        truth_value = true;
    }
    return;
}