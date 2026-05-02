// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "core/ga_prdxpr_common.hpp"

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

    // check rules: there must be not leading or trailing space_str()
    // (check for potential user error, when manually entering the coefficients)
    for (auto const& e : lcoeff) {
        if (e.starts_with(space_str()) || e.ends_with(space_str())) {
            fmt::println("mvec: {}", lcoeff);
            fmt::println("mvec coefficient: '{}'", e);
            throw std::runtime_error("Product coefficients on left hand side must not "
                                     "start or end with space_str().");
        }
    }
    for (auto const& e : rcoeff) {
        if (e.starts_with(space_str()) || e.ends_with(space_str())) {
            fmt::println("mvec: {}", rcoeff);
            fmt::println("mvec coefficient: '{}'", e);
            throw std::runtime_error("Product coefficients on right hand side must not "
                                     "start or end with space_str().");
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

            if (lcoeff[i].starts_with(minus_str())) {
                toggle_bool(is_negative);
                lhs = lcoeff[i].substr(1, lcoeff[i].size() - 1);
            }

            if (rcoeff[j].starts_with(minus_str())) {
                toggle_bool(is_negative);
                rhs = rcoeff[j].substr(1, rcoeff[j].size() - 1);
            }

            if (lhs == zero_str() || rhs == zero_str()) {
                // product result is zero, thus coefficient becomes zero
                prd_coeff_tab[i][j] = zero_str();
            }
            else if (is_negative) {
                prd_coeff_tab[i][j] =
                    minus_str() + lhs + space_str() + operator_str + space_str() + rhs;
            }
            else {
                prd_coeff_tab[i][j] =
                    lhs + space_str() + operator_str + space_str() + rhs;
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

            bool is_negative{false}; // must be "false" for toggle logic below to work
            std::string lhs{coeff_tab[i][j]};
            std::string rhs{basis_tab[i][j]};

            if (coeff_tab[i][j].starts_with(minus_str())) {
                toggle_bool(is_negative);
                lhs = coeff_tab[i][j].substr(1, coeff_tab[i][j].size() - 1);
            }

            if (basis_tab[i][j].starts_with(minus_str())) {
                toggle_bool(is_negative);
                rhs = basis_tab[i][j].substr(1, basis_tab[i][j].size() - 1);
            }

            if (basis_tab[i][j] == zero_str() || coeff_tab[i][j] == zero_str()) {
                // if either coefficient or basis entry is zero the product must be zero
                prd_tab[i][j] = zero_str();
            }
            else if (is_negative) {
                prd_tab[i][j] = minus_str() + lhs + space_str() + rhs;
            }
            else {
                prd_tab[i][j] = lhs + space_str() + rhs;
            }
        }
    }

    return prd_tab;
}


mvec_coeff apply_rules_to_mv(mvec_coeff const& coeff, mvec_rules const& rules)
{

    // check rules: there must be not leading or trailing space_str()
    // (check for potential user error, when defininig the rules manually entered)
    for (auto const& r : rules) {
        if (r.first.starts_with(space_str()) || r.first.ends_with(space_str()) ||
            r.second.starts_with(space_str()) || r.second.ends_with(space_str())) {
            fmt::println("Rule key: '{}', rule value: '{}'", r.first, r.second);
            throw std::runtime_error("Rules must not start or end with space_str().");
        }
    }

    mvec_coeff mvec(coeff.size());

    for (size_t i = 0; i < coeff.size(); ++i) {

        bool is_negative{false}; // must be "false" for toggle logic below to work
        std::string value{coeff[i]};

        // remove the minus-sign, if it is present in the input
        if (coeff[i].starts_with(minus_str())) {
            toggle_bool(is_negative);
            value = coeff[i].substr(1, coeff[i].size() - 1);
        }

        // fmt::println("i: {}, value: '{}', is_negative: {}", i, value, is_negative);

        // apply the rule to the input, after the input is without intial minus_str()
        // zero_str() value will be left unchanged, i.e. implicit rule "0" -> "0"
        if (value != zero_str()) {
            value = rules.at(value);
        }

        // fmt::println("i: {}, value: '{}', is_negative: {}", i, value, is_negative);

        // in case the replaced value does already have a minus sign on its own
        if (value.starts_with(minus_str())) {
            toggle_bool(is_negative);
            value = value.substr(1, value.size() - 1);
        }

        // fmt::println("i: {}, value: '{}', is_negative: {}", i, value, is_negative);

        if (is_negative && (value != zero_str())) { // zero doesn't need a minus sign
            mvec[i] = minus_str() + value;
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

    // check rules: there must be not leading or trailing space_str()
    // (check for potential user error, when defininig the rules manually entered)
    for (auto const& r : rules) {
        if (r.first.starts_with(space_str()) || r.first.ends_with(space_str()) ||
            r.second.starts_with(space_str()) || r.second.ends_with(space_str())) {
            fmt::println("Rule key: '{}', rule value: '{}'", r.first, r.second);
            throw std::runtime_error("Rules must not start or end with empty_str().");
        }
    }

    // reserve a full table of the required size
    prd_table prd_tab(tab.size(), mvec_coeff(tab.size()));

    for (size_t i = 0; i < tab.size(); ++i) {
        for (size_t j = 0; j < tab.size(); ++j) {

            bool is_negative{false}; // must be "false" for toggle logic below to work
            std::string value{tab[i][j]};

            // remove the minus-sign, if it is present in the input
            if (tab[i][j].starts_with(minus_str())) {
                toggle_bool(is_negative);
                value = tab[i][j].substr(1, tab[i][j].size() - 1);
            }

            // fmt::println("i: {}, j: {}, value: '{}', is_negative: {}", i, j, value,
            //              is_negative);

            // apply the rule to the input, after the input is without intial
            // minus_str() zero_str() value will be left unchanged, i.e. implicit rule "0"
            // -> "0"
            if (value != zero_str()) {
                value = rules.at(value);
            }

            // fmt::println("i: {}, j: {}, value: '{}', is_negative: {}", i, j, value,
            //              is_negative);

            // in case the replaced value does already have a minus sign on its own
            if (value.starts_with(minus_str())) {
                toggle_bool(is_negative);
                value = value.substr(1, value.size() - 1);
            }

            // fmt::println("i: {}, j: {}, value: '{}', is_negative: {}", i, j, value,
            //              is_negative);

            if (is_negative && (value != zero_str())) { // zero doesn't need a minus sign
                prd_tab[i][j] = minus_str() + value;
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
    // fmt::println("get_prd_tab - coeff product table:");
    // print_prd_tab(mv_coeff_tab);
    // fmt::println("");

    auto prd_tab = combine_coeff_and_basis_prd_tabs(mv_coeff_tab, basis_tab);
    // fmt::println("{}: product table:", prd_name);
    // print_prd_tab(prd_tab);
    // fmt::println("");

    return prd_tab;
}

prd_table get_prd_tab_sym(prd_table const& tab)
{
    //
    // a*b = 0.5*(a*b + b*a) + 0.5*(a*b - b*a) = (a*b)_sym + (a*b)_asym
    //
    // provide the symmetric part of the product table (a*b)_sym = 0.5*(a*b + b*a)
    //

    // make sure the table has quadradic shape
    for (size_t i = 0; i < tab.size(); ++i) {
        if (tab[i].size() != tab.size()) {
            throw std::runtime_error("Product tables must be square matrices. Sizes of "
                                     "rows and columns must match.");
        }
    }

    prd_table sym{tab}; // create a table of the right size by copying

    for (size_t i = 0; i < tab.size(); ++i) {
        for (size_t j = 0; j < tab[i].size(); ++j) {

            // sym[i][j] = tab[i][j] is already contained from initialization
            if (sym[i][j] == tab[j][i]) continue; // symmetric part already contained

            if ((sym[i][j].starts_with(minus_str()) &&
                 sym[i][j].substr(1, sym[i][j].size()) == tab[j][i]) ||
                (tab[j][i].starts_with(minus_str()) &&
                 tab[j][i].substr(1, tab[j][i].size()) == sym[i][j])) {
                sym[i][j] = zero_str(); // asymmetric part, so set element to zero_str()
            }
            else if (sym[i][j].starts_with(minus_str()) && tab[j][i] == zero_str()) {
                sym[i][j] = "-0.5 * " + sym[i][j].substr(1, sym[i][j].size());
            }
            else if (sym[i][j] == zero_str() && tab[j][i].starts_with(minus_str())) {
                sym[i][j] = "-0.5 * " + tab[j][i].substr(1, tab[j][i].size());
            }
            else if (!sym[i][j].starts_with(minus_str()) && tab[j][i] == zero_str()) {
                sym[i][j] = "0.5 * " + sym[i][j];
            }
            else if (sym[i][j] == zero_str() && !tab[j][i].starts_with(minus_str())) {
                sym[i][j] = "0.5 * " + tab[j][i];
            }
            else {
                // this occurs only when tab does not contain expected base
                // vector products
                sym[i][j] = "0.5 * (" + tab[i][j] + space_str() + plus_str() +
                            space_str() + tab[j][i] + ")";
            }
        }
    }

    // fmt::println("");
    // print_prd_tab(sym);
    // fmt::println("");

    return sym;
}
prd_table get_prd_tab_asym(prd_table const& tab)
{
    //
    // a*b = 0.5*(a*b + b*a) + 0.5*(a*b - b*a) = (a*b)_sym + (a*b)_asym
    //
    // provide the asymmetric part of the product table (a*b)_asym = 0.5*(a*b - b*a)
    //

    // make sure the table has quadradic shape
    for (size_t i = 0; i < tab.size(); ++i) {
        if (tab[i].size() != tab.size()) {
            throw std::runtime_error("Product tables must be square matrices. Sizes of "
                                     "rows and columns must match.");
        }
    }

    prd_table asym{tab}; // create a table of the right size by copying

    for (size_t i = 0; i < tab.size(); ++i) {
        for (size_t j = 0; j < tab[i].size(); ++j) {

            // asym[i][j] = tab[i][j] is already contained from initialization
            if ((asym[i][j].starts_with(minus_str()) &&
                 asym[i][j].substr(1, asym[i][j].size()) == tab[j][i]) ||
                (tab[j][i].starts_with(minus_str()) &&
                 tab[j][i].substr(1, tab[j][i].size()) == asym[i][j])) {
                continue; // asymmetric part already contained
            }

            if (asym[i][j] == tab[j][i]) {
                asym[i][j] = zero_str(); // symmetric part, so set element to zero_str()
            }
            else if (asym[i][j].starts_with(minus_str()) && tab[j][i] == zero_str()) {
                asym[i][j] = "-0.5 * " + asym[i][j].substr(1, asym[i][j].size());
            }
            else if (asym[i][j] == zero_str() && tab[j][i].starts_with(minus_str())) {
                asym[i][j] = "0.5 * " + tab[j][i].substr(1, tab[j][i].size());
            }
            else if (!asym[i][j].starts_with(minus_str()) && tab[j][i] == zero_str()) {
                asym[i][j] = "0.5 * " + asym[i][j];
            }
            else if (asym[i][j] == zero_str() && !tab[j][i].starts_with(minus_str())) {
                asym[i][j] = "-0.5 * " + tab[j][i];
            }
            else {
                // this occurs only when tab does not contain expected base
                // vector products
                asym[i][j] = "0.5 * (" + tab[i][j] + space_str() + minus_str() +
                             space_str() + tab[j][i] + ")";
            }
        }
    }

    // fmt::println("");
    // print_prd_tab(asym);
    // fmt::println("");

    return asym;
}

prd_table add_prd_tab(prd_table const& ta, prd_table const& tb)
{
    // add two product tables

    // make sure the table a has quadradic shape
    for (size_t i = 0; i < ta.size(); ++i) {
        if (ta[i].size() != ta.size()) {
            throw std::runtime_error("Product tables must be square matrices. Sizes of "
                                     "rows and columns must match.");
        }
    }
    // make sure the table b has quadradic shape
    for (size_t i = 0; i < tb.size(); ++i) {
        if (tb[i].size() != tb.size()) {
            throw std::runtime_error("Product tables must be square matrices. Sizes of "
                                     "rows and columns must match.");
        }
    }
    // make sure both matrices are of same size
    if (ta.size() != tb.size()) {
        throw std::runtime_error("Product table sizes must match to add tables.");
    }

    prd_table tres;
    tres.resize(ta.size());
    for (size_t i = 0; i < ta.size(); ++i) {
        tres[i].resize(ta.size());
    }

    for (size_t i = 0; i < ta.size(); ++i) {
        for (size_t j = 0; j < ta[i].size(); ++j) {

            tres[i][j] = ta[i][j];

            if (tres[i][j] != zero_str() && tb[i][j] != zero_str()) {

                tres[i][j] += space_str() + plus_str() + space_str() + tb[i][j];
            }

            if (tres[i][j] == zero_str() && tb[i][j] != zero_str()) {

                tres[i][j] = tb[i][j];
            }
        }
    }

    return tres;
}

prd_table negate_prd_tab(prd_table const& t)
{
    // negate the product table (change only entries != "0")

    prd_table tres;
    tres.resize(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        tres[i].resize(t.size());
    }

    for (size_t i = 0; i < t.size(); ++i) {
        for (size_t j = 0; j < t[i].size(); ++j) {

            tres[i][j] = t[i][j];

            if (tres[i][j] != zero_str()) {

                if (tres[i][j].starts_with(minus_str())) {
                    tres[i][j] = tres[i][j].substr(1);
                }
                else {
                    tres[i][j] = minus_str() + tres[i][j];
                }
            }
        }
    }

    return tres;
}

mvec_coeff get_mv_from_prd_tab(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                               filter_2d lfilter, filter_2d rfilter, brace_switch brsw)
{

    // make sure sizes match as required
    if (prd_tab.size() != mv_basis.size()) {
        throw std::runtime_error("Multivector size of product table and multivector "
                                 "basis size must match.");
    }

    return extractor(prd_tab, mv_basis, get_coeff_filter(lfilter),
                     get_coeff_filter(rfilter), brsw);
}

mvec_coeff get_mv_from_prd_tab(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                               filter_3d lfilter, filter_3d rfilter, brace_switch brsw)
{

    // make sure sizes match as required
    if (prd_tab.size() != mv_basis.size()) {
        throw std::runtime_error("Multivector size of product table and multivector "
                                 "basis size must match.");
    }

    return extractor(prd_tab, mv_basis, get_coeff_filter(lfilter),
                     get_coeff_filter(rfilter), brsw);
}

mvec_coeff get_mv_from_prd_tab(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                               filter_4d lfilter, filter_4d rfilter, brace_switch brsw)
{

    // make sure sizes match as required
    if (prd_tab.size() != mv_basis.size()) {
        throw std::runtime_error("Multivector size of product table and multivector "
                                 "basis size must match.");
    }

    return extractor(prd_tab, mv_basis, get_coeff_filter(lfilter),
                     get_coeff_filter(rfilter), brsw);
}


mvec_coeff extractor(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                     mvec_coeff_filter const& lcoeff_filter,
                     mvec_coeff_filter const& rcoeff_filter, brace_switch brsw)
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

                if (prd_tab[i][j].ends_with(space_str() + basis_element)) {
                    // found a contribution
                    if (mv_prd[b] == empty_str()) {
                        // add the first entry (including a minus-sign, if present)
                        // remove empty_str() and basis element from prd_tab entry
                        mv_prd[b] += prd_tab[i][j].substr(0, prd_tab[i][j].size() -
                                                                 mv_basis[b].size() - 1);
                    }
                    else {
                        // add the 2nd and following entries
                        if (prd_tab[i][j].starts_with(minus_str())) {
                            // remove empty_str() and basis element and minus sign from
                            // prd_tab entry
                            mv_prd[b] +=
                                space_str() + minus_str() + space_str() +
                                prd_tab[i][j].substr(1, prd_tab[i][j].size() -
                                                            mv_basis[b].size() - 2);
                        }
                        else {
                            // remove empty_str() and basis element from prd_tab entry
                            mv_prd[b] +=
                                space_str() + plus_str() + space_str() +
                                prd_tab[i][j].substr(0, prd_tab[i][j].size() -
                                                            mv_basis[b].size() - 1);
                        }
                    }
                }
            }
        }
        // replace remaining empty elements by zero_str()
        if (mv_prd[b] == empty_str()) {
            mv_prd[b] = zero_str();
        }
        // add braces for non-empty elements, if requested
        if (brsw == brace_switch::use_braces && mv_prd[b] != zero_str()) {
            mv_prd[b] = brace_open_str() + mv_prd[b] + brace_close_str();
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

    // print elements right aligned fitting to max_width with vertical separator like
    // tables
    for (size_t i = 0; i < mv.size(); ++i) {
        const auto& e = mv[i];
        fmt::print("[ ");
        fmt::print("{:>{w}}", e, fmt::arg("w", max_width));
        fmt::println(" ] | {}", mv_basis[i]);
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

void print_prd_tab_with_headers(prd_table const& tab, mvec_coeff const& basis)
{
    if (tab.empty() || basis.empty()) {
        fmt::println("Error: Empty table or basis");
        return;
    }

    // Calculate maximum width for table cells and headers
    size_t max_cell_width = 0;
    for (const auto& row : tab) {
        for (const auto& cell : row) {
            max_cell_width = std::max(max_cell_width, cell.size());
        }
    }

    size_t max_header_width = 0;
    for (const auto& header : basis) {
        max_header_width = std::max(max_header_width, header.size());
    }

    size_t cell_width = std::max(max_cell_width, max_header_width);

    // Calculate max width for row headers
    size_t max_row_header_width = 0;
    for (const auto& header : basis) {
        max_row_header_width = std::max(max_row_header_width, header.size());
    }

    // Print column headers (multivector basis elements)
    // Headers align with table entries but without brackets: "  entry, entry, entry |
    // header"
    fmt::print("  "); // Small indent to align with table content inside brackets
    for (size_t col = 0; col < basis.size() && col < tab[0].size(); ++col) {
        fmt::print("{:>{w}}", basis[col], fmt::arg("w", cell_width));
        if (col < tab[0].size() - 1) {
            fmt::print(", ");
        }
    }
    fmt::print("   |");
    fmt::println("");

    // Calculate total width for horizontal separator
    size_t table_width = 2; // "[ "
    for (size_t col = 0; col < basis.size() && col < tab[0].size(); ++col) {
        table_width += cell_width;
        if (col < tab[0].size() - 1) {
            table_width += 2; // ", "
        }
    }
    table_width += 5; // " ] | "
    table_width += max_row_header_width;

    // Print horizontal separator line
    fmt::print("{:-<{}}", "", table_width);
    fmt::println("");

    // Print table rows with row headers on the right
    for (size_t row = 0; row < tab.size() && row < basis.size(); ++row) {
        // Print table row
        fmt::print("[ ");
        for (size_t col = 0; col < tab[row].size(); ++col) {
            fmt::print("{:>{w}}", tab[row][col], fmt::arg("w", cell_width));
            if (col < tab[row].size() - 1) {
                fmt::print(", ");
            }
        }

        // Print vertical separator and row header (basis element for this row) on the
        // right, left-aligned
        fmt::print(" ] | {}", basis[row]);
        fmt::println("");
    }
}

void print_prd_rules(prd_rules const& rules, std::string const& title,
                     mvec_coeff const& basis_order, std::string const& operator_symbol)
{
    fmt::println("\n=== {} ===", title);

    // Print in grade order by iterating through basis in order
    for (const auto& a : basis_order) {
        for (const auto& b : basis_order) {
            std::string key = a + " " + operator_symbol + " " + b;

            auto it = rules.find(key);
            if (it != rules.end()) {
                fmt::println("    {{\"{}\", \"{}\"}},", key, it->second);
            }
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
// helper functions
////////////////////////////////////////////////////////////////////////////////
mvec_coeff_filter get_coeff_filter(filter_2d filter)
{

    mvec_coeff_filter filter_vec(4); // 2d multivector has 4 basis components
                                     // 1 scalar, 2 vector, and 1 pseudoscalar components
    switch (filter) {
        case filter_2d::mv:
            filter_vec = {1, 1, 1, 1}; // all components
            break;
        case filter_2d::mv_e:
            filter_vec = {1, 0, 0, 1}; // even grade components
            break;
        case filter_2d::s:
            filter_vec = {1, 0, 0, 0}; // scalar element
            break;
        case filter_2d::vec:
            filter_vec = {0, 1, 1, 0}; // vector components
            break;
        case filter_2d::ps:
            filter_vec = {0, 0, 0, 1}; // pseudoscalar element (=bivector in 2d)
            break;
        default:
            std::unreachable();
    }
    return filter_vec;
}

mvec_coeff_filter get_coeff_filter(filter_3d filter)
{

    mvec_coeff_filter filter_vec(8); // 3d multivector has 8 basis components:
                                     // 1 scalar, 3 vector, 3 bivector,
                                     // and 1 pseudoscalar components

    switch (filter) {
        case filter_3d::mv:
            filter_vec = {1, 1, 1, 1, 1, 1, 1, 1}; // all components
            break;
        case filter_3d::mv_e:
            filter_vec = {1, 0, 0, 0, 1, 1, 1, 0}; // even grade components
            break;
        case filter_3d::mv_u:
            filter_vec = {0, 1, 1, 1, 0, 0, 0, 1}; // odd grade components
            break;
        case filter_3d::s:
            filter_vec = {1, 0, 0, 0, 0, 0, 0, 0}; // scalar element
            break;
        case filter_3d::vec:
            filter_vec = {0, 1, 1, 1, 0, 0, 0, 0}; // vector components
            break;
        case filter_3d::bivec:
            filter_vec = {0, 0, 0, 0, 1, 1, 1, 0}; // bivector components
            break;
        case filter_3d::ps:
            filter_vec = {0, 0, 0, 0, 0, 0, 0, 1}; // pseudoscalar element
                                                   // (=trivector in 3d)
            break;
        default:
            std::unreachable();
    }
    return filter_vec;
}

mvec_coeff_filter get_coeff_filter(filter_4d filter)
{

    mvec_coeff_filter filter_vec(16); // 4d multivector has 16 basis components:
                                      // 1 scalar, 4 vector, 6 bivector, 4 trivector,
                                      // and 1 pseudoscalar components

    switch (filter) {
        case filter_4d::mv:
            filter_vec = {1, 1, 1, 1, 1, 1, 1, 1,
                          1, 1, 1, 1, 1, 1, 1, 1}; // all components
            break;
        case filter_4d::mv_e:
            filter_vec = {1, 0, 0, 0, 0, 1, 1, 1,
                          1, 1, 1, 0, 0, 0, 0, 1}; // even grade components
            break;
        case filter_4d::mv_u:
            filter_vec = {0, 1, 1, 1, 1, 0, 0, 0,
                          0, 0, 0, 1, 1, 1, 1, 0}; // odd grade components
            break;
        case filter_4d::s:
            filter_vec = {1, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0}; // scalar element
            break;
        case filter_4d::vec:
            filter_vec = {0, 1, 1, 1, 1, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 0}; // vector components
            break;
        case filter_4d::bivec:
            filter_vec = {0, 0, 0, 0, 0, 1, 1, 1,
                          1, 1, 1, 0, 0, 0, 0, 0}; // bivector components
            break;
        case filter_4d::trivec:
            filter_vec = {0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 1, 1, 1, 1, 0}; // trivector components
            break;
        case filter_4d::ps:
            filter_vec = {0, 0, 0, 0, 0, 0, 0, 0,
                          0, 0, 0, 0, 0, 0, 0, 1}; // pseudoscalar element
                                                   // (=quadvector in 4d)
            break;
        default:
            std::unreachable();
    }
    return filter_vec;
}

////////////////////////////////////////////////////////////////////////////////
// Basis prefix validation and extraction
////////////////////////////////////////////////////////////////////////////////

std::string extract_basis_prefix(mvec_coeff const& vector_basis)
{
    if (vector_basis.empty()) {
        throw std::runtime_error(
            "extract_basis_prefix: vector_basis is empty, cannot extract prefix");
    }

    std::string const& first_vector = vector_basis[0];
    if (first_vector.empty()) {
        throw std::runtime_error(
            "extract_basis_prefix: first vector element is empty, cannot extract prefix");
    }

    // Extract all leading non-digit characters as the prefix
    // Examples: "e1" -> "e", "e12" -> "e", "g0" -> "g", "g123" -> "g"
    std::string prefix;
    for (char c : first_vector) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            break;
        }
        prefix += c;
    }

    if (prefix.empty()) {
        throw std::runtime_error(fmt::format("extract_basis_prefix: could not extract "
                                             "prefix from first vector element '{}'",
                                             first_vector));
    }

    return prefix;
}

void validate_basis_consistency(mvec_coeff const& multivector_basis,
                                std::vector<mvec_coeff> const& basis_kvec,
                                std::string const& expected_prefix,
                                std::string const& scalar_name)
{
    // Validate that basis_kvec is not empty
    if (basis_kvec.empty()) {
        throw std::runtime_error(
            "validate_basis_consistency: basis_kvec is empty, cannot validate");
    }

    // Validate that multivector_basis is not empty
    if (multivector_basis.empty()) {
        throw std::runtime_error(
            "validate_basis_consistency: multivector_basis is empty, cannot validate");
    }

    // Validate that scalar element is present in multivector_basis
    if (multivector_basis[0] != scalar_name) {
        throw std::runtime_error(fmt::format(
            "validate_basis_consistency: first element of multivector_basis must be "
            "scalar '{}', but found '{}'",
            scalar_name, multivector_basis[0]));
    }

    // Validate that scalar element is present in basis_kvec[0]
    if (basis_kvec[0].empty()) {
        throw std::runtime_error(
            "validate_basis_consistency: basis_kvec[0] (scalar grade) is empty");
    }

    if (basis_kvec[0].size() != 1) {
        throw std::runtime_error(fmt::format(
            "validate_basis_consistency: basis_kvec[0] (scalar grade) must contain "
            "exactly 1 element, but contains {} elements",
            basis_kvec[0].size()));
    }

    if (basis_kvec[0][0] != scalar_name) {
        throw std::runtime_error(fmt::format(
            "validate_basis_consistency: basis_kvec[0][0] must be scalar '{}', but "
            "found '{}'",
            scalar_name, basis_kvec[0][0]));
    }

    // Check consistency in multivector_basis (all non-scalar elements)
    for (size_t i = 1; i < multivector_basis.size(); ++i) {
        std::string const& element = multivector_basis[i];

        // Check if element starts with expected_prefix
        if (!element.starts_with(expected_prefix)) {
            throw std::runtime_error(fmt::format(
                "validate_basis_consistency: multivector_basis[{}] element '{}' does "
                "not start with expected prefix '{}'",
                i, element, expected_prefix));
        }
    }

    // Check consistency across all grades in basis_kvec (skip grade 0, already checked)
    for (size_t grade = 1; grade < basis_kvec.size(); ++grade) {
        mvec_coeff const& grade_basis = basis_kvec[grade];

        for (std::string const& element : grade_basis) {
            // Check if element starts with expected_prefix
            if (!element.starts_with(expected_prefix)) {
                throw std::runtime_error(fmt::format(
                    "validate_basis_consistency: basis_kvec[{}] element '{}' does "
                    "not start with expected prefix '{}'",
                    grade, element, expected_prefix));
            }
        }
    }

    // Validate that multivector_basis matches flattened basis_kvec
    // Flatten basis_kvec into a single vector
    mvec_coeff flattened_basis_kvec;
    for (mvec_coeff const& grade_basis : basis_kvec) {
        for (std::string const& element : grade_basis) {
            flattened_basis_kvec.push_back(element);
        }
    }

    // Check that sizes match
    if (multivector_basis.size() != flattened_basis_kvec.size()) {
        throw std::runtime_error(fmt::format(
            "validate_basis_consistency: multivector_basis size ({}) does not match "
            "flattened basis_kvec size ({})",
            multivector_basis.size(), flattened_basis_kvec.size()));
    }

    // Check that all elements match in order
    for (size_t i = 0; i < multivector_basis.size(); ++i) {
        if (multivector_basis[i] != flattened_basis_kvec[i]) {
            throw std::runtime_error(fmt::format(
                "validate_basis_consistency: multivector_basis[{}] = '{}' does not "
                "match flattened basis_kvec[{}] = '{}'",
                i, multivector_basis[i], i, flattened_basis_kvec[i]));
        }
    }
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

////////////////////////////////////////////////////////////////////////////////
// transwedge product helper functions for alternative geometric product calculation
////////////////////////////////////////////////////////////////////////////////

mvec_coeff init_zero_multivector(size_t size)
{
    mvec_coeff c0(size);
    for (auto& e : c0) {
        e = zero_str();
    }
    return c0;
}

prd_table init_zero_product_table(size_t size)
{
    auto c0 = init_zero_multivector(size);
    prd_table tab;
    tab.resize(size);
    for (size_t i = 0; i < size; ++i) {
        tab[i] = c0;
    }
    return tab;
}

void print_product_tables_by_grade(std::vector<prd_table> const& tables,
                                   std::vector<mvec_coeff> const& basis_kvec)
{
    if (tables.empty()) return;

    auto print_res0 = init_zero_product_table(tables[0].size());

    // flatten basis_kvec to multivector basis (for printing)
    mvec_coeff basis;
    for (size_t k = 0; k < basis_kvec.size(); ++k) {
        for (size_t cnt = 0; cnt < basis_kvec[k].size(); ++cnt) {
            basis.push_back(basis_kvec[k][cnt]);
        }
    }

    size_t idx = 0;
    for (size_t k = 0; k < basis_kvec.size(); ++k) {
        auto p_tab = print_res0; // empty table initialized with "0" elements

        for (size_t cnt = 0; cnt < basis_kvec[k].size(); ++cnt) {
            if (idx < tables.size()) {
                p_tab = add_prd_tab(p_tab, tables[idx]);
            }
            ++idx;
        }
        fmt::println("order k = {}:", k);
        print_prd_tab_with_headers(p_tab, basis);
        fmt::println("");
    }
}