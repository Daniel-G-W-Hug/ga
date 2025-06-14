// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_config_pga3dp.hpp"
#include "ga_prdxpr_pga3dp.hpp"

namespace configs {

////////////////////////////////////////////////////////////////////////////////
// PGA3DP Configuration
////////////////////////////////////////////////////////////////////////////////

const AlgebraConfig pga3dp_config = {
    .name = "pga3dp",
    .description = "Projective Geometric Algebra 3D+ - G(3,0,1)",
    .basis = mv3dp_basis,
    .product_rules = {{"gpr", gpr_pga3dp_rules},
                      {"wdg", wdg_pga3dp_rules},
                      {"dot", dot_pga3dp_rules},
                      {"rgpr", gpr_pga3dp_rules}},  // regressive geometric uses gpr with complement transformations
    .coefficients = {{"A", mv3dp_coeff_A},
                     {"B", mv3dp_coeff_B},
                     {"M", mv3dp_coeff_M},
                     {"M1", mv3dp_coeff_M1},
                     {"M2", mv3dp_coeff_M2},
                     {"M_even", mv3dp_coeff_M_even},
                     {"M_odd", mv3dp_coeff_M_odd},
                     {"A_even", mv3dp_coeff_A_even},
                     {"B_even", mv3dp_coeff_B_even},
                     {"A_odd", mv3dp_coeff_A_odd},
                     {"B_odd", mv3dp_coeff_B_odd},
                     {"R_even", mv3dp_coeff_R_even},
                     {"R_odd", mv3dp_coeff_R_odd},
                     {"R_rev_even", mv3dp_coeff_R_rev_even},
                     {"R_rev_odd", mv3dp_coeff_R_rev_odd},
                     {"R_rrev_even", mv3dp_coeff_R_rrev_even},
                     {"R_rrev_odd", mv3dp_coeff_R_rrev_odd},
                     {"svBtps", mv3dp_coeff_svBtps},
                     {"svBtps1", mv3dp_coeff_svBtps1},
                     {"svBtps2", mv3dp_coeff_svBtps2}},
    .complement_rules = {{"lcmpl", lcmpl_pga3dp_rules},
                         {"rcmpl", rcmpl_pga3dp_rules},
                         {"cmpl",
                          lcmpl_pga3dp_rules}, // Use lcmpl for regressive products
                         {"bulk_dual", right_bulk_dual_pga3dp_rules},
                         {"left_bulk_dual", left_bulk_dual_pga3dp_rules},
                         {"weight_dual", right_weight_dual_pga3dp_rules},
                         {"left_weight_dual", left_weight_dual_pga3dp_rules}},
    .filter_names = {"s", "vec", "bivec", "trivec", "ps", "mv_e", "mv_u", "mv"},
    .dimension = 4 // 3d projective = 4d homogeneous
};

// PGA3DP Product-specific case definitions
const std::vector<ProductCase> pga3dp_geometric_cases = {
    // 41 cases from reference - PGA3DP geometric product cases
    {"mv * mv -> mv", "A", "B", "mv", "mv", true},
    {"mv * mv_e -> mv", "A", "B_even", "mv", "mv_e", true},
    {"mv_e * mv -> mv", "A_even", "B", "mv_e", "mv", true},
    {"mv * mv_u -> mv", "A", "B_odd", "mv", "mv_u", true},
    {"mv_u * mv -> mv", "A_odd", "B", "mv_u", "mv", true},
    {"mv * ps -> mv", "A", "svBtps", "mv", "ps", true},
    {"ps * mv -> mv", "svBtps", "B", "ps", "mv", true},
    {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e", true},
    {"mv_u * mv_u -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u", true},
    {"mv_e * mv_u -> mv_u", "A_even", "B_odd", "mv_e", "mv_u", true},
    {"mv_u * mv_e -> mv_u", "A_odd", "B_even", "mv_u", "mv_e", true},
    {"mv_e * ps -> mv_e", "A_even", "svBtps", "mv_e", "ps", true},
    {"ps * mv_e -> mv_e", "svBtps", "B_even", "ps", "mv_e", true},
    {"mv_u * ps -> mv_u", "A_odd", "svBtps", "mv_u", "ps", true},
    {"ps * mv_u -> mv_u", "svBtps", "B_odd", "ps", "mv_u", true},
    {"mv_e * trivec -> mv_u", "M_even", "svBtps", "mv_e", "trivec", true},
    {"trivec * mv_e -> mv_u", "svBtps", "M_even", "trivec", "mv_e", true},
    {"mv_e * bivec -> mv_e", "M_even", "svBtps", "mv_e", "bivec", true},
    {"bivec * mv_e -> mv_e", "svBtps", "M_even", "bivec", "mv_e", true},
    {"mv_u * bivec -> mv_u", "M_odd", "svBtps", "mv_u", "bivec", true},
    {"bivec * mv_u -> mv_u", "svBtps", "M_odd", "bivec", "mv_u", true},
    {"mv_e * vec -> mv_u", "A_even", "svBtps", "mv_e", "vec", true},
    {"vec * mv_e -> mv_u", "svBtps", "B_even", "vec", "mv_e", true},
    {"ps * s -> ps", "svBtps", "svBtps", "ps", "s", true},
    {"s * ps -> ps", "svBtps", "svBtps", "s", "ps", true},
    {"trivec * trivec -> mv_e", "svBtps1", "svBtps2", "trivec", "trivec", true},
    {"trivec * bivec -> mv_u", "svBtps", "svBtps", "trivec", "bivec", true},
    {"bivec * trivec -> mv_u", "svBtps", "svBtps", "bivec", "trivec", true},
    {"trivec * vec -> mv_e", "svBtps", "svBtps", "trivec", "vec", true},
    {"vec * trivec -> mv_e", "svBtps", "svBtps", "vec", "trivec", true},
    {"trivec * s -> trivec", "svBtps", "svBtps", "trivec", "s", true},
    {"s * trivec -> trivec", "svBtps", "svBtps", "s", "trivec", true},
    {"bivec * bivec -> mv_e", "svBtps1", "svBtps2", "bivec", "bivec", true},
    {"bivec * vec -> mv_u", "svBtps", "svBtps", "bivec", "vec", true},
    {"vec * bivec -> mv_u", "svBtps", "svBtps", "vec", "bivec", true},
    {"bivec * s -> bivec", "svBtps", "svBtps", "bivec", "s", true},
    {"s * bivec -> bivec", "svBtps", "svBtps", "s", "bivec", true},
    {"vec * vec -> mv_e", "svBtps1", "svBtps2", "vec", "vec", true},
    {"vec * s -> vec", "svBtps", "svBtps", "vec", "s", true},
    {"s * vec -> vec", "svBtps", "svBtps", "s", "vec", true},
    {"s * s -> s", "svBtps1", "svBtps2", "s", "s", true}};

const std::vector<ProductCase> pga3dp_commutator_cases = {
    // 7 cases from reference - PGA3DP commutator product cases
    {"cmt(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"cmt(trivec,trivec) -> bivec", "svBtps1", "svBtps2", "trivec", "trivec", true},
    {"cmt(trivec,bivec) -> trivec", "svBtps", "svBtps", "trivec", "bivec", true},
    {"cmt(bivec,trivec) -> trivec", "svBtps", "svBtps", "bivec", "trivec", true},
    {"cmt(bivec,bivec) -> bivec", "svBtps1", "svBtps2", "bivec", "bivec", true},
    {"cmt(bivec,vec) -> vec", "svBtps", "svBtps", "bivec", "vec", true},
    {"cmt(vec,bivec) -> vec", "svBtps", "svBtps", "vec", "bivec", true},
    {"cmt(vec,vec) -> bivec", "svBtps1", "svBtps2", "vec", "vec", true}};

const std::vector<ProductCase> pga3dp_wedge_cases = {
    // 16 cases from reference - PGA3DP wedge product cases
    {"mv ^ mv -> mv", "A", "B", "mv", "mv", true},
    {"ps ^ s -> ps", "svBtps", "svBtps", "ps", "s", true},
    {"s ^ ps -> ps", "svBtps", "svBtps", "s", "ps", true},
    {"trivec ^ vec -> ps", "svBtps", "svBtps", "trivec", "vec", true},
    {"vec ^ trivec -> ps", "svBtps", "svBtps", "vec", "trivec", true},
    {"trivec ^ s -> trivec", "svBtps", "svBtps", "trivec", "s", true},
    {"s ^ trivec -> trivec", "svBtps", "svBtps", "s", "trivec", true},
    {"bivec ^ bivec -> ps", "svBtps1", "svBtps2", "bivec", "bivec", true},
    {"bivec ^ vec -> trivec", "svBtps", "svBtps", "bivec", "vec", true},
    {"vec ^ bivec -> trivec", "svBtps", "svBtps", "vec", "bivec", true},
    {"bivec ^ s -> bivec", "svBtps", "svBtps", "bivec", "s", true},
    {"s ^ bivec -> bivec", "svBtps", "svBtps", "s", "bivec", true},
    {"vec ^ vec -> bivec", "svBtps1", "svBtps2", "vec", "vec", true},
    {"vec ^ s -> vec", "svBtps", "svBtps", "vec", "s", true},
    {"s ^ vec -> vec", "svBtps", "svBtps", "s", "vec", true},
    {"s ^ s -> s", "svBtps1", "svBtps2", "s", "s", true}};

const std::vector<ProductCase> pga3dp_inner_cases = {
    // 6 cases from reference - PGA3DP inner product cases
    {"dot(mv,mv) -> s", "A", "B", "mv", "mv", true},
    {"dot(ps,ps) -> s", "svBtps1", "svBtps2", "ps", "ps", true},
    {"dot(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec", "trivec", true},
    {"dot(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec", true},
    {"dot(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec", true},
    {"dot(s,s) -> s", "svBtps1", "svBtps2", "s", "s", true}};

const std::vector<ProductCase> pga3dp_regressive_wedge_cases = {
    // 9 cases from reference - PGA3DP regressive wedge product cases
    {"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"rwdg(ps,trivec) -> trivec", "svBtps", "svBtps", "ps", "trivec", true},
    {"rwdg(trivec,ps) -> trivec", "svBtps", "svBtps", "trivec", "ps", true},
    {"rwdg(trivec,trivec) -> bivec", "svBtps1", "svBtps2", "trivec", "trivec", true},
    {"rwdg(trivec,bivec) -> vec", "svBtps", "svBtps", "trivec", "bivec", true},
    {"rwdg(bivec,trivec) -> vec", "svBtps", "svBtps", "bivec", "trivec", true},
    {"rwdg(trivec,vec) -> s", "svBtps", "svBtps", "trivec", "vec", true},
    {"rwdg(vec,trivec) -> s", "svBtps", "svBtps", "vec", "trivec", true},
    {"rwdg(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec", true}};

const std::vector<ProductCase> pga3dp_regressive_inner_cases = {
    // 6 cases from reference - PGA3DP regressive inner product cases
    {"rdot(mv,mv) -> ps", "A", "B", "mv", "mv", true},
    {"rdot(ps,ps) -> ps", "svBtps1", "svBtps2", "ps", "ps", true},
    {"rdot(trivec,trivec) -> ps", "svBtps1", "svBtps2", "trivec", "trivec", true},
    {"rdot(bivec,bivec) -> ps", "svBtps1", "svBtps2", "bivec", "bivec", true},
    {"rdot(vec,vec) -> ps", "svBtps1", "svBtps2", "vec", "vec", true},
    {"rdot(s,s) -> ps", "svBtps1", "svBtps2", "s", "s", true}};

const std::vector<ProductCase> pga3dp_regressive_geometric_cases = {
    // 13 cases from reference - PGA3DP regressive geometric product cases
    {"rgpr(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"rgpr(mv_e,mv_e) -> mv_e", "A_even", "B_even", "mv_e", "mv_e", true},
    {"rgpr(mv_u,mv_e) -> mv_u", "A_odd", "B_even", "mv_u", "mv_e", true},
    {"rgpr(mv_e,trivec) -> mv_u", "M_even", "svBtps", "mv_e", "trivec", true},
    {"rgpr(mv_u,trivec) -> mv_e", "M_odd", "svBtps", "mv_u", "trivec", true},
    {"rgpr(mv_e,bivec) -> mv_e", "M_even", "svBtps", "mv_e", "bivec", true},
    {"rgpr(mv_u,bivec) -> mv_u", "M_odd", "svBtps", "mv_u", "bivec", true},
    {"rgpr(mv_e,vec) -> mv_u", "M_even", "svBtps", "mv_e", "vec", true},
    {"rgpr(trivec,trivec) -> mv_e", "svBtps1", "svBtps2", "trivec", "trivec", true},
    {"rgpr(trivec,bivec) -> mv_u", "svBtps", "svBtps", "trivec", "bivec", true},
    {"rgpr(trivec,vec) -> mv_e", "svBtps", "svBtps", "trivec", "vec", true},
    {"rgpr(bivec,bivec) -> mv_e", "svBtps1", "svBtps2", "bivec", "bivec", true},
    {"rgpr(vec,vec) -> mv_e", "svBtps1", "svBtps2", "vec", "vec", true}};

const std::vector<ProductCase> pga3dp_right_bulk_contraction_cases = {
    // Cases from reference - PGA3DP right bulk contraction cases
    {"rbulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"rbulk_contract(ps,ps) -> 0", "svBtps1", "svBtps2", "ps", "ps", true},
    {"rbulk_contract(ps,trivec) -> vec", "svBtps", "svBtps", "ps", "trivec", true},
    {"rbulk_contract(trivec,ps) -> 0", "svBtps", "svBtps", "trivec", "ps", true},
    {"rbulk_contract(ps,bivec) -> bivec", "svBtps", "svBtps", "ps", "bivec", true},
    {"rbulk_contract(bivec,ps) -> 0", "svBtps", "svBtps", "bivec", "ps", true},
    {"rbulk_contract(ps,vec) -> trivec", "svBtps", "svBtps", "ps", "vec", true},
    {"rbulk_contract(vec,ps) -> 0", "svBtps", "svBtps", "vec", "ps", true},
    {"rbulk_contract(ps,s) -> ps", "svBtps", "svBtps", "ps", "s", true},
    {"rbulk_contract(s,ps) -> 0", "svBtps", "svBtps", "s", "ps", true},
    {"rbulk_contract(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec", "trivec", true},
    {"rbulk_contract(trivec,bivec) -> vec", "svBtps", "svBtps", "trivec", "bivec", true},
    {"rbulk_contract(bivec,trivec) -> 0", "svBtps", "svBtps", "bivec", "trivec", true},
    {"rbulk_contract(trivec,vec) -> bivec", "svBtps", "svBtps", "trivec", "vec", true},
    {"rbulk_contract(vec,trivec) -> 0", "svBtps", "svBtps", "vec", "trivec", true},
    {"rbulk_contract(trivec,s) -> trivec", "svBtps", "svBtps", "trivec", "s", true},
    {"rbulk_contract(s,trivec) -> 0", "svBtps", "svBtps", "s", "trivec", true},
    {"rbulk_contract(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec", true},
    {"rbulk_contract(bivec,vec) -> vec", "svBtps", "svBtps", "bivec", "vec", true},
    {"rbulk_contract(vec,bivec) -> 0", "svBtps", "svBtps", "vec", "bivec", true},
    {"rbulk_contract(bivec,s) -> bivec", "svBtps", "svBtps", "bivec", "s", true},
    {"rbulk_contract(s,bivec) -> 0", "svBtps", "svBtps", "s", "bivec", true},
    {"rbulk_contract(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec", true},
    {"rbulk_contract(vec,s) -> vec", "svBtps", "svBtps", "vec", "s", true},
    {"rbulk_contract(s,vec) -> 0", "svBtps", "svBtps", "s", "vec", true},
    {"rbulk_contract(s,s) -> s", "svBtps1", "svBtps2", "s", "s", true}
};

const std::vector<ProductCase> pga3dp_right_weight_contraction_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductCase> pga3dp_left_bulk_contraction_cases = {
    // Cases from reference - PGA3DP left bulk contraction cases
    {"lbulk_contract(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"lbulk_contract(ps,ps) -> 0", "svBtps1", "svBtps2", "ps", "ps", true},
    {"lbulk_contract(ps,trivec) -> 0", "svBtps", "svBtps", "ps", "trivec", true},
    {"lbulk_contract(trivec,ps) -> vec", "svBtps", "svBtps", "trivec", "ps", true},
    {"lbulk_contract(ps,bivec) -> 0", "svBtps", "svBtps", "ps", "bivec", true},
    {"lbulk_contract(bivec,ps) -> bivec", "svBtps", "svBtps", "bivec", "ps", true},
    {"lbulk_contract(ps,vec) -> 0", "svBtps", "svBtps", "ps", "vec", true},
    {"lbulk_contract(vec,ps) -> trivec", "svBtps", "svBtps", "vec", "ps", true},
    {"lbulk_contract(ps,s) -> 0", "svBtps", "svBtps", "ps", "s", true},
    {"lbulk_contract(s,ps) -> ps", "svBtps", "svBtps", "s", "ps", true},
    {"lbulk_contract(trivec,trivec) -> s", "svBtps1", "svBtps2", "trivec", "trivec", true},
    {"lbulk_contract(trivec,bivec) -> 0", "svBtps", "svBtps", "trivec", "bivec", true},
    {"lbulk_contract(bivec,trivec) -> vec", "svBtps", "svBtps", "bivec", "trivec", true},
    {"lbulk_contract(trivec,vec) -> 0", "svBtps", "svBtps", "trivec", "vec", true},
    {"lbulk_contract(vec,trivec) -> bivec", "svBtps", "svBtps", "vec", "trivec", true},
    {"lbulk_contract(trivec,s) -> 0", "svBtps", "svBtps", "trivec", "s", true},
    {"lbulk_contract(s,trivec) -> trivec", "svBtps", "svBtps", "s", "trivec", true},
    {"lbulk_contract(bivec,bivec) -> s", "svBtps1", "svBtps2", "bivec", "bivec", true},
    {"lbulk_contract(bivec,vec) -> 0", "svBtps", "svBtps", "bivec", "vec", true},
    {"lbulk_contract(vec,bivec) -> vec", "svBtps", "svBtps", "vec", "bivec", true},
    {"lbulk_contract(bivec,s) -> 0", "svBtps", "svBtps", "bivec", "s", true},
    {"lbulk_contract(s,bivec) -> bivec", "svBtps", "svBtps", "s", "bivec", true},
    {"lbulk_contract(vec,vec) -> s", "svBtps1", "svBtps2", "vec", "vec", true},
    {"lbulk_contract(vec,s) -> 0", "svBtps", "svBtps", "vec", "s", true},
    {"lbulk_contract(s,vec) -> vec", "svBtps", "svBtps", "s", "vec", true},
    {"lbulk_contract(s,s) -> s", "svBtps1", "svBtps2", "s", "s", true}
};

const std::vector<ProductCase> pga3dp_left_weight_contraction_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductCase> pga3dp_right_bulk_expansion_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductCase> pga3dp_right_weight_expansion_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductCase> pga3dp_left_bulk_expansion_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductCase> pga3dp_left_weight_expansion_cases = {
    // Only basis product table shown in reference - no specific product cases listed
};

const std::vector<ProductDefinition> pga3dp_products = {
    {"gpr", "geometric product", "gpr", mul_str, ProductType::geometric_product, false,
     true, pga3dp_geometric_cases},
    {"cmt", "commutator product", "gpr", mul_str, ProductType::commutator_product, true,
     true, pga3dp_commutator_cases},
    {"wdg", "wedge product", "wdg", wdg_str, ProductType::wedge_product, false, true,
     pga3dp_wedge_cases},
    {"dot", "inner product", "dot", mul_str, ProductType::inner_product, false, true,
     pga3dp_inner_cases},
    {"rwdg", "regressive wedge product", "wdg", wdg_str, ProductType::regressive_wedge,
     false, true, pga3dp_regressive_wedge_cases},
    {"rdot", "regressive inner product", "dot", mul_str, ProductType::regressive_inner,
     false, true, pga3dp_regressive_inner_cases},
    {"rgpr", "regressive geometric product", "gpr", mul_str,
     ProductType::regressive_geometric, false, true, pga3dp_regressive_geometric_cases},
    {"rbulk_contract", "right bulk contraction", "wdg", ">>", ProductType::right_bulk_contraction, false, true, pga3dp_right_bulk_contraction_cases},
    {"rweight_contract", "right weight contraction", "wdg", ">>", ProductType::right_weight_contraction, false, true, pga3dp_right_weight_contraction_cases},
    {"lbulk_contract", "left bulk contraction", "wdg", "<<", ProductType::left_bulk_contraction, false, true, pga3dp_left_bulk_contraction_cases},
    {"lweight_contract", "left weight contraction", "wdg", "<<", ProductType::left_weight_contraction, false, true, pga3dp_left_weight_contraction_cases},
    {"rbulk_expand", "right bulk expansion", "wdg", "><", ProductType::right_bulk_expansion, false, true, pga3dp_right_bulk_expansion_cases},
    {"rweight_expand", "right weight expansion", "wdg", "><", ProductType::right_weight_expansion, false, true, pga3dp_right_weight_expansion_cases},
    {"lbulk_expand", "left bulk expansion", "wdg", "<>", ProductType::left_bulk_expansion, false, true, pga3dp_left_bulk_expansion_cases},
    {"lweight_expand", "left weight expansion", "wdg", "<>", ProductType::left_weight_expansion, false, true, pga3dp_left_weight_expansion_cases},
    
    // Regressive sandwich products - motor * object * rrev(motor)
    {"regressive_sandwich", "regressive sandwich product", "rgpr", mul_str, ProductType::regressive_sandwich, false, true, {}}
};

} // namespace configs