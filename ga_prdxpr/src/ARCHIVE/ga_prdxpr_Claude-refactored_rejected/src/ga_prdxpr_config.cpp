// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_config.hpp"
#include "ga_prdxpr_config_ega2d.hpp"
#include "ga_prdxpr_config_ega3d.hpp"
#include "ga_prdxpr_config_pga2dp.hpp"
#include "ga_prdxpr_config_pga3dp.hpp"

namespace configs {

////////////////////////////////////////////////////////////////////////////////
// Helper function implementations
////////////////////////////////////////////////////////////////////////////////

filter_2d get_filter_2d_from_string(const std::string& filter_name) {
    if (filter_name == "s") return filter_2d::s;
    if (filter_name == "vec") return filter_2d::vec;
    if (filter_name == "ps") return filter_2d::ps;
    if (filter_name == "mv_e") return filter_2d::mv_e;
    if (filter_name == "mv") return filter_2d::mv;
    throw std::invalid_argument("Unknown filter_2d: " + filter_name);
}

filter_3d get_filter_3d_from_string(const std::string& filter_name) {
    if (filter_name == "s") return filter_3d::s;
    if (filter_name == "vec") return filter_3d::vec;
    if (filter_name == "bivec") return filter_3d::bivec;
    if (filter_name == "ps") return filter_3d::ps;
    if (filter_name == "mv_e") return filter_3d::mv_e;
    if (filter_name == "mv_u") return filter_3d::mv_u;
    if (filter_name == "mv") return filter_3d::mv;
    throw std::invalid_argument("Unknown filter_3d: " + filter_name);
}

filter_4d get_filter_4d_from_string(const std::string& filter_name) {
    if (filter_name == "s") return filter_4d::s;
    if (filter_name == "vec") return filter_4d::vec;
    if (filter_name == "bivec") return filter_4d::bivec;
    if (filter_name == "trivec") return filter_4d::trivec;
    if (filter_name == "ps") return filter_4d::ps;
    if (filter_name == "mv_e") return filter_4d::mv_e;
    if (filter_name == "mv_u") return filter_4d::mv_u;
    if (filter_name == "mv") return filter_4d::mv;
    throw std::invalid_argument("Unknown filter_4d: " + filter_name);
}

bool validate_algebra_config(const AlgebraConfig& config) {
    // Basic validation - can be extended
    return !config.name.empty() && 
           !config.basis.empty() && 
           !config.product_rules.empty() &&
           !config.coefficients.empty();
}

} // namespace configs