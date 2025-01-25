#pragma once

// author: Daniel Hug, 2024

#include <string>
#include <vector>

// multivector coefficients (scalars)
using mvec_coeff = std::vector<std::string>;

// multivector product table - contains products of scalars and basis elements
using prd_table = std::vector<mvec_coeff>;

// combine all products
using namespace std::literals::string_literals;
constexpr std::string empty_str{""s};
constexpr std::string zero_str{"0"s};
constexpr std::string one_str{"1"s};
constexpr std::string plus_str{"+"s};
constexpr std::string minus_str{"-"s};
constexpr std::string mul_str{"*"s};
constexpr std::string space_str{" "s};

prd_table mv_coeff_to_coeff_prd_tab(mvec_coeff const& lcoeff, mvec_coeff const& rcoeff);
prd_table combine_coeff_and_basis_prd_tabs(prd_table const& coeff_tab,
                                           prd_table const& basis_tab);
mvec_coeff extract_mv_from_prd_tab(prd_table const& prd_tab, mvec_coeff const& mv_basis);

void generate_and_print_product(std::string const& prd_name, prd_table const& basis_tab,
                                mvec_coeff const& mv_basis, mvec_coeff const& mv_lcoeff,
                                mvec_coeff const& mv_rcoeff);

void print_mvec(mvec_coeff const& mv, mvec_coeff const& mv_basis);
void print_prd_tab(prd_table const& tab);