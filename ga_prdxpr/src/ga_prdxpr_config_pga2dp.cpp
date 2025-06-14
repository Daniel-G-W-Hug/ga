// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_config_pga2dp.hpp"
#include "ga_prdxpr_pga2dp.hpp"

namespace configs {

////////////////////////////////////////////////////////////////////////////////
// PGA2DP Configuration
////////////////////////////////////////////////////////////////////////////////

const AlgebraConfig pga2dp_config = {
    .name = "pga2dp",
    .description = "Projective Geometric Algebra 2D+ - G(2,0,1)",
    .basis = mv2dp_basis,
    .product_rules = {
        {"gpr", gpr_pga2dp_rules},
        {"wdg", wdg_pga2dp_rules},
        {"dot", dot_pga2dp_rules},
        {"rgpr", gpr_pga2dp_rules}  // regressive geometric uses gpr with complement transformations
    },
    .coefficients = {
        {"A", mv2dp_coeff_A},
        {"B", mv2dp_coeff_B},
        {"M", mv2dp_coeff_M},
        {"M1", mv2dp_coeff_M1},
        {"M2", mv2dp_coeff_M2},
        {"M_even", mv2dp_coeff_M_even},
        {"M_odd", mv2dp_coeff_M_odd},
        {"A_even", mv2dp_coeff_A_even},
        {"B_even", mv2dp_coeff_B_even},
        {"A_odd", mv2dp_coeff_A_odd},
        {"B_odd", mv2dp_coeff_B_odd},
        {"R_even", mv2dp_coeff_R_even},
        {"R_odd", mv2dp_coeff_R_odd},
        {"R_rev_even", mv2dp_coeff_R_rev_even},
        {"R_rev_odd", mv2dp_coeff_R_rev_odd},
        {"R_rrev_even", mv2dp_coeff_R_rrev_even},
        {"R_rrev_odd", mv2dp_coeff_R_rrev_odd},
        {"svBps", mv2dp_coeff_svBps},
        {"svBps1", mv2dp_coeff_svBps1},
        {"svBps2", mv2dp_coeff_svBps2}
    },
    .complement_rules = {
        {"cmpl", cmpl_pga2dp_rules},
        {"lcmpl", cmpl_pga2dp_rules},  // left complement = main complement in PGA2DP
        {"rcmpl", cmpl_pga2dp_rules},  // right complement = main complement in PGA2DP
        {"bulk_dual", bulk_dual_pga2dp_rules},
        {"weight_dual", weight_dual_pga2dp_rules}
    },
    .filter_names = {"s", "vec", "bivec", "ps", "mv_e", "mv_u", "mv"},
    .dimension = 3  // 2d projective = 3d homogeneous
};

// PGA2DP Product-specific case definitions
const std::vector<ProductCase> pga2dp_geometric_cases = {
    // 36 cases from reference - PGA2DP geometric product cases
    {"mv * mv -> mv", "A", "B", "mv", "mv", true},
    {"mv * mv_e -> mv", "A", "B_even", "mv", "mv_e", true},
    {"mv_e * mv -> mv", "A_even", "B", "mv_e", "mv", true},
    {"mv * mv_u -> mv", "A", "B_odd", "mv", "mv_u", true},
    {"mv_u * mv -> mv", "A_odd", "B", "mv_u", "mv", true},
    {"mv * ps -> mv", "A", "svBps", "mv", "ps", true},
    {"ps * mv -> mv", "svBps", "B", "ps", "mv", true},
    {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e", true},
    {"mv_u * mv_u -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u", true},
    {"mv_e * mv_u -> mv_u", "A_even", "B_odd", "mv_e", "mv_u", true},
    {"mv_u * mv_e -> mv_u", "A_odd", "B_even", "mv_u", "mv_e", true},
    {"mv_e * ps -> mv_u", "A_even", "svBps", "mv_e", "ps", true},
    {"ps * mv_e -> mv_u", "svBps", "B_even", "ps", "mv_e", true},
    {"mv_u * ps -> mv_e", "A_odd", "svBps", "mv_u", "ps", true},
    {"ps * mv_u -> mv_e", "svBps", "B_odd", "ps", "mv_u", true},
    {"mv_e * bivec -> mv_e", "M_even", "svBps", "mv_e", "bivec", true},
    {"bivec * mv_e -> mv_e", "svBps", "M_even", "bivec", "mv_e", true},
    {"mv_u * bivec -> mv_u", "M_odd", "svBps", "mv_u", "bivec", true},
    {"bivec * mv_u -> mv_u", "svBps", "M_odd", "bivec", "mv_u", true},
    {"mv_e * vec -> mv_u", "A_even", "svBps", "mv_e", "vec", true},
    {"vec * mv_e -> mv_u", "svBps", "B_even", "vec", "mv_e", true},
    {"ps * ps -> 0", "svBps1", "svBps2", "ps", "ps", true},
    {"ps * bivec -> vec", "svBps", "svBps", "ps", "bivec", true},
    {"bivec * ps -> vec", "svBps", "svBps", "bivec", "ps", true},
    {"ps * vec -> bivec", "svBps", "svBps", "ps", "vec", true},
    {"vec * ps -> bivec", "svBps", "svBps", "vec", "ps", true},
    {"ps * s -> ps", "svBps", "svBps", "ps", "s", true},
    {"s * ps -> ps", "svBps", "svBps", "s", "ps", true},
    {"bivec * bivec -> mv_e", "svBps1", "svBps2", "bivec", "bivec", true},
    {"bivec * vec -> mv_u", "svBps", "svBps", "bivec", "vec", true},
    {"vec * bivec -> mv_u", "svBps", "svBps", "vec", "bivec", true},
    {"bivec * s -> bivec", "svBps", "svBps", "bivec", "s", true},
    {"s * bivec -> bivec", "svBps", "svBps", "s", "bivec", true},
    {"vec * vec -> mv_e", "svBps1", "svBps2", "vec", "vec", true},
    {"vec * s -> vec", "svBps", "svBps", "vec", "s", true},
    {"s * vec -> vec", "svBps", "svBps", "s", "vec", true},
    {"s * s -> s", "svBps1", "svBps2", "s", "s", true}
};

const std::vector<ProductCase> pga2dp_commutator_cases = {
    // 5 cases from reference - PGA2DP commutator product cases  
    {"cmt(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"cmt(bivec,bivec) -> bivec", "svBps1", "svBps2", "bivec", "bivec", true},
    {"cmt(bivec,vec) -> vec", "svBps", "svBps", "bivec", "vec", true},
    {"cmt(vec,bivec) -> vec", "svBps", "svBps", "vec", "bivec", true},
    {"cmt(vec,vec) -> bivec", "svBps1", "svBps2", "vec", "vec", true}
};

const std::vector<ProductCase> pga2dp_wedge_cases = {
    // 17 cases from reference - PGA2DP wedge product cases
    {"mv ^ mv -> mv", "A", "B", "mv", "mv", true},
    {"ps ^ ps -> 0", "svBps1", "svBps2", "ps", "ps", true},
    {"ps ^ bivec -> 0", "svBps", "svBps", "ps", "bivec", true},
    {"bivec ^ ps -> 0", "svBps", "svBps", "bivec", "ps", true},
    {"ps ^ vec -> 0", "svBps", "svBps", "ps", "vec", true},
    {"vec ^ ps -> 0", "svBps", "svBps", "vec", "ps", true},
    {"ps ^ s -> ps", "svBps", "svBps", "ps", "s", true},
    {"s ^ ps -> ps", "svBps", "svBps", "s", "ps", true},
    {"bivec ^ bivec -> 0", "svBps1", "svBps2", "bivec", "bivec", true},
    {"bivec ^ vec -> ps", "svBps", "svBps", "bivec", "vec", true},
    {"vec ^ bivec -> ps", "svBps", "svBps", "vec", "bivec", true},
    {"bivec ^ s -> bivec", "svBps", "svBps", "bivec", "s", true},
    {"s ^ bivec -> bivec", "svBps", "svBps", "s", "bivec", true},
    {"vec ^ vec -> bivec", "svBps1", "svBps2", "vec", "vec", true},
    {"vec ^ s -> vec", "svBps", "svBps", "vec", "s", true},
    {"s ^ vec -> vec", "svBps", "svBps", "s", "vec", true},
    {"s ^ s -> s", "svBps1", "svBps2", "s", "s", true}
};

const std::vector<ProductCase> pga2dp_inner_cases = {
    // 5 cases from reference - PGA2DP inner product cases
    {"dot(mv,mv) -> s", "A", "B", "mv", "mv", true},
    {"dot(ps,ps) -> 0", "svBps1", "svBps2", "ps", "ps", true},
    {"dot(bivec,bivec) -> s", "svBps1", "svBps2", "bivec", "bivec", true},
    {"dot(vec,vec) -> s", "svBps1", "svBps2", "vec", "vec", true},
    {"dot(s,s) -> s", "svBps1", "svBps2", "s", "s", true}
};

const std::vector<ProductCase> pga2dp_regressive_wedge_cases = {
    // 6 cases from reference - PGA2DP regressive wedge product cases
    {"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"rwdg(ps,bivec) -> bivec", "svBps", "svBps", "ps", "bivec", true},
    {"rwdg(bivec,ps) -> bivec", "svBps", "svBps", "bivec", "ps", true},
    {"rwdg(bivec,bivec) -> vec", "svBps1", "svBps2", "bivec", "bivec", true},
    {"rwdg(bivec,vec) -> s", "svBps", "svBps", "bivec", "vec", true},
    {"rwdg(vec,bivec) -> s", "svBps", "svBps", "vec", "bivec", true}
};

const std::vector<ProductCase> pga2dp_regressive_inner_cases = {
    // 4 cases from reference - PGA2DP regressive inner product cases
    {"rdot(ps,ps) -> ps", "svBps1", "svBps2", "ps", "ps", true},
    {"rdot(bivec,bivec) -> ps", "svBps1", "svBps2", "bivec", "bivec", true},
    {"rdot(vec,vec) -> ps", "svBps1", "svBps2", "vec", "vec", true},
    {"rdot(s,s) -> 0", "svBps1", "svBps2", "s", "s", true}
};

const std::vector<ProductCase> pga2dp_regressive_geometric_cases = {
    // 12 cases from reference - PGA2DP regressive geometric product cases
    {"rgpr(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"rgpr(mv_e,mv_e) -> mv_u", "A_even", "B_even", "mv_e", "mv_e", true},
    {"rgpr(mv_u,mv_u) -> mv_u", "A_odd", "B_odd", "mv_u", "mv_u", true},
    {"rgpr(mv_e,mv_u) -> mv_e", "A_even", "B_odd", "mv_e", "mv_u", true},
    {"rgpr(mv_u,mv_e) -> mv_e", "A_odd", "B_even", "mv_u", "mv_e", true},
    {"rgpr(mv_e,bivec) -> mv_u", "M_even", "svBps", "mv_e", "bivec", true},
    {"rgpr(mv_u,bivec) -> mv_e", "M_odd", "svBps", "mv_u", "bivec", true},
    {"rgpr(mv_e,vec) -> mv_e", "M_even", "svBps", "mv_e", "vec", true},
    {"rgpr(mv_u,vec) -> mv_u", "M_odd", "svBps", "mv_u", "vec", true},
    {"rgpr(bivec,bivec) -> mv_u", "svBps1", "svBps2", "bivec", "bivec", true},
    {"rgpr(bivec,vec) -> mv_e", "svBps", "svBps", "bivec", "vec", true},
    {"rgpr(vec,vec) -> mv_u", "svBps1", "svBps2", "vec", "vec", true}
};

const std::vector<ProductCase> pga2dp_right_bulk_contraction_cases = {
    // 16 cases from reference - PGA2DP right bulk contraction cases
    {"rbulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"rbulk_contract(ps,ps) -> 0", "svBps1", "svBps2", "ps", "ps", true},
    {"rbulk_contract(ps,bivec) -> vec", "svBps", "svBps", "ps", "bivec", true},
    {"rbulk_contract(bivec,ps) -> 0", "svBps", "svBps", "bivec", "ps", true},
    {"rbulk_contract(vec,ps) -> 0", "svBps", "svBps", "vec", "ps", true},
    {"rbulk_contract(ps,s) -> ps", "svBps", "svBps", "ps", "s", true},
    {"rbulk_contract(s,ps) -> 0", "svBps", "svBps", "s", "ps", true},
    {"rbulk_contract(bivec,bivec) -> s", "svBps1", "svBps2", "bivec", "bivec", true},
    {"rbulk_contract(bivec,vec) -> vec", "svBps", "svBps", "bivec", "vec", true},
    {"rbulk_contract(vec,bivec) -> 0", "svBps", "svBps", "vec", "bivec", true},
    {"rbulk_contract(bivec,s) -> bivec", "svBps", "svBps", "bivec", "s", true},
    {"rbulk_contract(s,bivec) -> 0", "svBps", "svBps", "s", "bivec", true},
    {"rbulk_contract(vec,vec) -> s", "svBps1", "svBps2", "vec", "vec", true},
    {"rbulk_contract(vec,s) -> vec", "svBps", "svBps", "vec", "s", true},
    {"rbulk_contract(s,vec) -> 0", "svBps", "svBps", "s", "vec", true},
    {"rbulk_contract(s,s) -> s", "svBps1", "svBps2", "s", "s", true}
};

// Right weight contraction cases (based on reference - appears to only show basis table)
const std::vector<ProductCase> pga2dp_right_weight_contraction_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductCase> pga2dp_left_bulk_contraction_cases = {
    // 16 cases from reference - PGA2DP left bulk contraction cases
    {"lbulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"lbulk_contract(ps,ps) -> 0", "svBps1", "svBps2", "ps", "ps", true},
    {"lbulk_contract(ps,bivec) -> 0", "svBps", "svBps", "ps", "bivec", true},
    {"lbulk_contract(bivec,ps) -> vec", "svBps", "svBps", "bivec", "ps", true},
    {"lbulk_contract(ps,vec) -> 0", "svBps", "svBps", "ps", "vec", true},
    {"lbulk_contract(vec,ps) -> bivec", "svBps", "svBps", "vec", "ps", true},
    {"lbulk_contract(ps,s) -> 0", "svBps", "svBps", "ps", "s", true},
    {"lbulk_contract(s,ps) -> ps", "svBps", "svBps", "s", "ps", true},
    {"lbulk_contract(bivec,bivec) -> s", "svBps1", "svBps2", "bivec", "bivec", true},
    {"lbulk_contract(bivec,vec) -> 0", "svBps", "svBps", "bivec", "vec", true},
    {"lbulk_contract(vec,bivec) -> vec", "svBps", "svBps", "vec", "bivec", true},
    {"lbulk_contract(bivec,s) -> 0", "svBps", "svBps", "bivec", "s", true},
    {"lbulk_contract(s,bivec) -> bivec", "svBps", "svBps", "s", "bivec", true},
    {"lbulk_contract(vec,vec) -> s", "svBps1", "svBps2", "vec", "vec", true},
    {"lbulk_contract(vec,s) -> 0", "svBps", "svBps", "vec", "s", true},
    {"lbulk_contract(s,vec) -> vec", "svBps", "svBps", "s", "vec", true},
    {"lbulk_contract(s,s) -> s", "svBps1", "svBps2", "s", "s", true}
};

// Left weight contraction cases (based on reference - appears to only show basis table)
const std::vector<ProductCase> pga2dp_left_weight_contraction_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

// Expansion cases (based on reference - appears to only show basis tables)
const std::vector<ProductCase> pga2dp_right_bulk_expansion_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductCase> pga2dp_right_weight_expansion_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductCase> pga2dp_left_bulk_expansion_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductCase> pga2dp_left_weight_expansion_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductDefinition> pga2dp_products = {
    {"gpr", "geometric product", "gpr", mul_str, ProductType::geometric_product, false, true, pga2dp_geometric_cases},
    {"cmt", "commutator product", "gpr", mul_str, ProductType::commutator_product, true, true, pga2dp_commutator_cases},
    {"wdg", "wedge product", "wdg", wdg_str, ProductType::wedge_product, false, true, pga2dp_wedge_cases},
    {"dot", "inner product", "dot", mul_str, ProductType::inner_product, false, true, pga2dp_inner_cases},
    {"rwdg", "regressive wedge product", "wdg", wdg_str, ProductType::regressive_wedge, false, true, pga2dp_regressive_wedge_cases},
    {"rdot", "regressive inner product", "dot", mul_str, ProductType::regressive_inner, false, true, pga2dp_regressive_inner_cases},
    {"rgpr", "regressive geometric product", "gpr", mul_str, ProductType::regressive_geometric, false, true, pga2dp_regressive_geometric_cases},
    {"rbulk_contract", "right bulk contraction", "wdg", ">>", ProductType::right_bulk_contraction, false, true, pga2dp_right_bulk_contraction_cases},
    {"rweight_contract", "right weight contraction", "wdg", ">>", ProductType::right_weight_contraction, false, true, pga2dp_right_weight_contraction_cases},
    {"lbulk_contract", "left bulk contraction", "wdg", "<<", ProductType::left_bulk_contraction, false, true, pga2dp_left_bulk_contraction_cases},
    {"lweight_contract", "left weight contraction", "wdg", "<<", ProductType::left_weight_contraction, false, true, pga2dp_left_weight_contraction_cases},
    {"rbulk_expand", "right bulk expansion", "wdg", "><", ProductType::right_bulk_expansion, false, true, pga2dp_right_bulk_expansion_cases},
    {"rweight_expand", "right weight expansion", "wdg", "><", ProductType::right_weight_expansion, false, true, pga2dp_right_weight_expansion_cases},
    {"lbulk_expand", "left bulk expansion", "wdg", "<>", ProductType::left_bulk_expansion, false, true, pga2dp_left_bulk_expansion_cases},
    {"lweight_expand", "left weight expansion", "wdg", "<>", ProductType::left_weight_expansion, false, true, pga2dp_left_weight_expansion_cases},
    
    // Regressive sandwich products - motor * object * rrev(motor)
    {"regressive_sandwich", "regressive sandwich product", "rgpr", mul_str, ProductType::regressive_sandwich, false, true, {}}
};

} // namespace configs