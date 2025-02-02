#pragma once

// author: Daniel Hug, 2025

#include <map>
#include <string>
#include <vector>

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// multivector coefficients (scalars)
using mvec_coeff = std::vector<std::string>;

// multivector coefficient filter (1: coeff is used, 0: coeff is not used)
using mvec_coeff_filter = std::vector<int>;

// rules to allow for multivector function maps
using mvec_rules = std::map<std::string, std::string>;

// rules to simplify product mappings in tables
using prd_rules = std::map<std::string, std::string>;

// multivector: mv, even grade multivector: mv_e, scalar: s, vector: vec, pseudoscalar: ps
enum class filter_2d { mv, mv_e, s, vec, ps };

// multivector: mv, even grade multivector: mv_e, uneven grade multivector: mv_u,
// scalar: s, vector: vec, bivector: bivec, pseudoscalar: ps
enum class filter_3d { mv, mv_e, mv_u, s, vec, bivec, ps };

// multivector product table - contains products of scalars and basis elements
using prd_table = std::vector<mvec_coeff>;

// required literals for string handling
using namespace std::literals::string_literals;
constexpr inline std::string empty_str{""s};
constexpr inline std::string zero_str{"0"s};
constexpr inline std::string one_str{"1"s};
constexpr inline std::string plus_str{"+"s};
constexpr inline std::string minus_str{"-"s};
constexpr inline std::string mul_str{"*"s};
constexpr inline std::string wdg_str{"^"s};
constexpr inline std::string lcontr_str{"<<"s};
constexpr inline std::string rcontr_str{">>"s};
constexpr inline std::string space_str{" "s};

////////////////////////////////////////////////////////////////////////////////
// user related functions
////////////////////////////////////////////////////////////////////////////////

mvec_coeff get_filtered_mv(mvec_coeff const& mv, filter_2d filter = filter_2d::mv);
mvec_coeff get_filtered_mv(mvec_coeff const& mv, filter_3d filter = filter_3d::mv);

prd_table mv_coeff_to_coeff_prd_tab(mvec_coeff const& lcoeff, mvec_coeff const& rcoeff,
                                    std::string const& operator_str = mul_str);

prd_table combine_coeff_and_basis_prd_tabs(prd_table const& coeff_tab,
                                           prd_table const& basis_tab);

mvec_coeff apply_rules_to_mv(mvec_coeff const& coeff, mvec_rules const& rules);
prd_table apply_rules_to_tab(prd_table const& tab, prd_rules const& rules);

prd_table get_prd_tab(prd_table const& basis_tab, mvec_coeff const& mv_lcoeff,
                      mvec_coeff const& mv_rcoeff);

mvec_coeff get_mv_from_prd_tab(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                               filter_2d lfilter, filter_2d rfilter);
mvec_coeff get_mv_from_prd_tab(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                               filter_3d lfilter, filter_3d rfilter);

mvec_coeff extractor(prd_table const& prd_tab, mvec_coeff const& mv_basis,
                     mvec_coeff_filter const& lcoeff_filter,
                     mvec_coeff_filter const& rcoeff_filter);


////////////////////////////////////////////////////////////////////////////////
// printing
////////////////////////////////////////////////////////////////////////////////

void print_mvec(mvec_coeff const& mv, mvec_coeff const& mv_basis);
void print_prd_tab(prd_table const& tab);

////////////////////////////////////////////////////////////////////////////////
// helper functions (not directly intended for user)
////////////////////////////////////////////////////////////////////////////////
mvec_coeff_filter get_coeff_filter(filter_2d filter = filter_2d::mv);
mvec_coeff_filter get_coeff_filter(filter_3d filter = filter_3d::mv);

void toggle_bool(bool& value);