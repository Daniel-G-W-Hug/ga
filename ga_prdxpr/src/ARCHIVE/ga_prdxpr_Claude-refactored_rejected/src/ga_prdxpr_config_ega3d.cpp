// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_config_ega3d.hpp"
#include "ga_prdxpr_ega3d.hpp"

namespace configs {

////////////////////////////////////////////////////////////////////////////////
// EGA3D Configuration
////////////////////////////////////////////////////////////////////////////////

const AlgebraConfig ega3d_config = {
    .name = "ega3d",
    .description = "Euclidean Geometric Algebra 3D - G(3,0,0)",
    .basis = mv3d_basis,
    .product_rules = {
        {"gpr", gpr_ega3d_rules},
        {"wdg", wdg_ega3d_rules},
        {"dot", dot_ega3d_rules}
    },
    .coefficients = {
        {"A", mv3d_coeff_A},
        {"B", mv3d_coeff_B},
        {"M", mv3d_coeff_M},
        {"M1", mv3d_coeff_M1},
        {"M2", mv3d_coeff_M2},
        {"A_even", mv3d_coeff_A_even},
        {"B_even", mv3d_coeff_B_even},
        {"A_odd", mv3d_coeff_A_odd},
        {"B_odd", mv3d_coeff_B_odd},
        {"M_even", mv3d_coeff_M_even},
        {"M_odd", mv3d_coeff_M_odd},
        {"R_even", mv3d_coeff_R_even},
        {"R_rev_even", mv3d_coeff_R_rev_even},
        {"svBps", mv3d_coeff_svBps},
        {"svBps1", mv3d_coeff_svBps1},
        {"svBps2", mv3d_coeff_svBps2}
    },
    .complement_rules = {
        {"cmpl", cmpl_ega3d_rules},
        {"lcmpl", cmpl_ega3d_rules},  // EGA3D uses same complement for both
        {"rcmpl", cmpl_ega3d_rules}
    },
    .filter_names = {"s", "vec", "bivec", "ps", "mv_e", "mv_u", "mv"},
    .dimension = 3
};

// EGA3D Product-specific case definitions
const std::vector<ProductCase> ega3d_geometric_cases = {
    // Following EGA2D pattern: geometric product always uses A/B coefficients
    {"mv * mv -> mv", "A", "B", "mv", "mv", true},
    {"mv * mv_e -> mv", "A", "B_even", "mv", "mv_e", true},
    {"mv_e * mv -> mv", "A_even", "B", "mv_e", "mv", true},
    {"mv * mv_u -> mv", "A", "B_odd", "mv", "mv_u", true},
    {"mv_u * mv -> mv", "A_odd", "B", "mv_u", "mv", true},
    {"mv * ps -> mv", "A", "svBps", "mv", "ps", true},         // ps uses A 
    {"ps * mv -> mv", "svBps", "B", "ps", "mv", true},         // ps uses B
    {"mv * bivec -> mv", "M", "svBps", "mv", "bivec", true},   // bivec uses M
    {"bivec * mv -> mv", "svBps", "M", "bivec", "mv", true},   // bivec uses M  
    {"mv * vec -> mv", "M", "svBps", "mv", "vec", true},       // vec uses M
    {"vec * mv -> mv", "svBps", "M", "vec", "mv", true},       // vec uses M
    {"mv * s -> mv", "M", "svBps", "mv", "s", true},           // s uses M
    {"s * mv -> mv", "svBps", "M", "s", "mv", true},           // s uses M
    {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e", true},
    {"mv_u * mv_u -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u", true},
    {"mv_e * mv_u -> mv_u", "A_even", "B_odd", "mv_e", "mv_u", true},
    {"mv_u * mv_e -> mv_u", "A_odd", "B_even", "mv_u", "mv_e", true},
    {"mv_e * ps -> mv_u", "M_even", "svBps", "mv_e", "ps", true},     // M_even for mv_e×ps
    {"ps * mv_e -> mv_u", "svBps", "M_even", "ps", "mv_e", true},     // M_even for ps×mv_e  
    {"mv_u * ps -> mv_e", "M_odd", "svBps", "mv_u", "ps", true},      // M_odd for mv_u×ps
    {"ps * mv_u -> mv_e", "svBps", "M_odd", "ps", "mv_u", true},      // M_odd for ps×mv_u
    {"mv_e * bivec -> mv_e", "M_even", "svBps", "mv_e", "bivec", true},  // M_even for mv_e×bivec -> mv_e
    {"bivec * mv_e -> mv_e", "svBps", "M_even", "bivec", "mv_e", true},  // Fixed: bivec×mv_e -> mv_e
    {"mv_u * bivec -> mv_u", "M_odd", "svBps", "mv_u", "bivec", true},   // M_odd for mv_u×bivec -> mv_u  
    {"bivec * mv_u -> mv_u", "svBps", "M_odd", "bivec", "mv_u", true},   // Fixed: bivec×mv_u -> mv_u
    {"mv_e * vec -> mv_u", "M_even", "svBps", "mv_e", "vec", true},      // M_even for mv_e×vec -> mv_u
    {"vec * mv_e -> mv_u", "svBps", "M_even", "vec", "mv_e", true},      // Fixed: vec×mv_e -> mv_u
    {"mv_u * vec -> mv_e", "M_odd", "svBps", "mv_u", "vec", true},       // M_odd for mv_u×vec -> mv_e
    {"vec * mv_u -> mv_e", "svBps", "M_odd", "vec", "mv_u", true},       // Fixed: vec×mv_u -> mv_e
    {"mv_e * s -> mv_e", "M_even", "svBps", "mv_e", "s", true},          // M_even for mv_e×s
    {"s * mv_e -> mv_e", "svBps", "M_even", "s", "mv_e", true},          // M_even for s×mv_e
    {"mv_u * s -> mv_u", "M_odd", "svBps", "mv_u", "s", true},           // M_odd for mv_u×s
    {"s * mv_u -> mv_u", "svBps", "M_odd", "s", "mv_u", true},           // M_odd for s×mv_u
    {"ps * ps -> s", "svBps1", "svBps2", "ps", "ps", true},
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

const std::vector<ProductCase> ega3d_commutator_cases = {
    // 5 cases from reference: cmt(mv,mv), cmt(bivec,bivec), cmt(bivec,vec), cmt(vec,bivec), cmt(vec,vec)
    {"cmt(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"cmt(bivec,bivec) -> bivec", "svBps1", "svBps2", "bivec", "bivec", true},
    {"cmt(bivec,vec) -> vec", "svBps", "svBps", "bivec", "vec", true},
    {"cmt(vec,bivec) -> vec", "svBps", "svBps", "vec", "bivec", true},
    {"cmt(vec,vec) -> bivec", "svBps1", "svBps2", "vec", "vec", true}
};

const std::vector<ProductCase> ega3d_wedge_cases = {
    // 23 cases from reference - wedge product uses A/B for mv×mv, M for mv×bivec/vec/s cases
    {"mv ^ mv -> mv", "A", "B", "mv", "mv", true},
    {"mv ^ bivec -> mv", "M", "svBps", "mv", "bivec", true},
    {"bivec ^ mv -> mv", "svBps", "M", "bivec", "mv", true}, 
    {"mv ^ vec -> mv", "M", "svBps", "mv", "vec", true},
    {"vec ^ mv -> mv", "svBps", "M", "vec", "mv", true},
    {"mv ^ s -> mv", "M", "svBps", "mv", "s", true},
    {"s ^ mv -> mv", "svBps", "M", "s", "mv", true},
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

const std::vector<ProductCase> ega3d_inner_cases = {
    // 5 cases from reference: dot(mv,mv), dot(ps,ps), dot(bivec,bivec), dot(vec,vec), dot(s,s)
    {"dot(mv,mv) -> s", "A", "B", "mv", "mv", true},
    {"dot(ps,ps) -> s", "svBps1", "svBps2", "ps", "ps", true},
    {"dot(bivec,bivec) -> s", "svBps1", "svBps2", "bivec", "bivec", true},
    {"dot(vec,vec) -> s", "svBps1", "svBps2", "vec", "vec", true},
    {"dot(s,s) -> s", "svBps1", "svBps2", "s", "s", true}
};

const std::vector<ProductCase> ega3d_regressive_wedge_cases = {
    // 11 cases from reference - regressive wedge product 
    {"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"rwdg(ps,ps) -> ps", "svBps1", "svBps2", "ps", "ps", true},
    {"rwdg(ps,bivec) -> bivec", "svBps", "svBps", "ps", "bivec", true},
    {"rwdg(bivec,ps) -> bivec", "svBps", "svBps", "bivec", "ps", true},
    {"rwdg(ps,vec) -> vec", "svBps", "svBps", "ps", "vec", true},
    {"rwdg(vec,ps) -> vec", "svBps", "svBps", "vec", "ps", true},
    {"rwdg(ps,s) -> s", "svBps", "svBps", "ps", "s", true},
    {"rwdg(s,ps) -> s", "svBps", "svBps", "s", "ps", true},
    {"rwdg(bivec,bivec) -> vec", "svBps1", "svBps2", "bivec", "bivec", true},
    {"rwdg(bivec,vec) -> s", "svBps", "svBps", "bivec", "vec", true},
    {"rwdg(vec,bivec) -> s", "svBps", "svBps", "vec", "bivec", true}
};

const std::vector<ProductCase> ega3d_left_contraction_cases = {
    // 22 cases from reference - left contraction product
    {"mv << mv -> mv", "A", "B", "mv", "mv", true},
    {"bivec << mv -> mv", "svBps", "M", "bivec", "mv", true},
    {"vec << mv -> mv", "svBps", "M", "vec", "mv", true},
    {"s << mv -> mv", "svBps", "M", "s", "mv", true},
    {"s << mv_e -> mv_e", "svBps", "M_even", "s", "mv_e", true},
    {"s << mv_u -> mv_u", "svBps", "M_odd", "s", "mv_u", true},
    {"ps << ps -> s", "svBps1", "svBps2", "ps", "ps", true},
    {"ps << bivec -> 0", "svBps", "svBps", "ps", "bivec", true},
    {"bivec << ps -> vec", "svBps", "svBps", "bivec", "ps", true},
    {"ps << vec -> 0", "svBps", "svBps", "ps", "vec", true},
    {"vec << ps -> bivec", "svBps", "svBps", "vec", "ps", true},
    {"ps << s -> 0", "svBps", "svBps", "ps", "s", true},
    {"s << ps -> ps", "svBps", "svBps", "s", "ps", true},
    {"bivec << bivec -> s", "svBps1", "svBps2", "bivec", "bivec", true},
    {"bivec << vec -> 0", "svBps", "svBps", "bivec", "vec", true},
    {"vec << bivec -> vec", "svBps", "svBps", "vec", "bivec", true},
    {"bivec << s -> 0", "svBps", "svBps", "bivec", "s", true},
    {"s << bivec -> bivec", "svBps", "svBps", "s", "bivec", true},
    {"vec << vec -> s", "svBps1", "svBps2", "vec", "vec", true},
    {"vec << s -> 0", "svBps", "svBps", "vec", "s", true},
    {"s << vec -> vec", "svBps", "svBps", "s", "vec", true},
    {"s << s -> s", "svBps1", "svBps2", "s", "s", true}
};

const std::vector<ProductCase> ega3d_right_contraction_cases = {
    // 22 cases from reference - right contraction product
    {"mv >> mv -> mv", "A", "B", "mv", "mv", true},
    {"mv >> bivec -> mv", "M", "svBps", "mv", "bivec", true},
    {"mv >> vec -> mv", "M", "svBps", "mv", "vec", true},
    {"mv >> s -> mv", "M", "svBps", "mv", "s", true},
    {"mv_e >> s -> mv_e", "M_even", "svBps", "mv_e", "s", true},
    {"mv_u >> s -> mv_u", "M_odd", "svBps", "mv_u", "s", true},
    {"ps >> ps -> s", "svBps1", "svBps2", "ps", "ps", true},
    {"ps >> bivec -> vec", "svBps", "svBps", "ps", "bivec", true},
    {"bivec >> ps -> 0", "svBps", "svBps", "bivec", "ps", true},
    {"ps >> vec -> bivec", "svBps", "svBps", "ps", "vec", true},
    {"vec >> ps -> 0", "svBps", "svBps", "vec", "ps", true},
    {"ps >> s -> ps", "svBps", "svBps", "ps", "s", true},
    {"s >> ps -> 0", "svBps", "svBps", "s", "ps", true},
    {"bivec >> bivec -> s", "svBps1", "svBps2", "bivec", "bivec", true},
    {"bivec >> vec -> vec", "svBps", "svBps", "bivec", "vec", true},
    {"vec >> bivec -> 0", "svBps", "svBps", "vec", "bivec", true},
    {"bivec >> s -> bivec", "svBps", "svBps", "bivec", "s", true},
    {"s >> bivec -> 0", "svBps", "svBps", "s", "bivec", true},
    {"vec >> vec -> s", "svBps1", "svBps2", "vec", "vec", true},
    {"vec >> s -> vec", "svBps", "svBps", "vec", "s", true},
    {"s >> vec -> 0", "svBps", "svBps", "s", "vec", true},
    {"s >> s -> s", "svBps1", "svBps2", "s", "s", true}
};

const std::vector<ProductCase> ega3d_left_expansion_cases = {
    // 2 cases from reference - left expansion product
    {"lexpand(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"lexpand(bivec,vec) -> bivec", "svBps", "svBps", "bivec", "vec", true}
};

const std::vector<ProductCase> ega3d_right_expansion_cases = {
    // 2 cases from reference - right expansion product
    {"rexpand(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"rexpand(vec,bivec) -> bivec", "svBps", "svBps", "vec", "bivec", true}
};

const std::vector<ProductDefinition> ega3d_products = {
    {"gpr", "geometric product", "gpr", mul_str, ProductType::geometric_product, false, true, ega3d_geometric_cases},
    {"cmt", "commutator product", "gpr", mul_str, ProductType::commutator_product, false, true, ega3d_commutator_cases},
    {"wdg", "wedge product", "wdg", wdg_str, ProductType::wedge_product, false, true, ega3d_wedge_cases},
    {"dot", "inner product", "dot", "dot", ProductType::inner_product, false, true, ega3d_inner_cases},
    {"rwdg", "regressive wedge product", "wdg", "rwdg", ProductType::regressive_wedge, false, true, ega3d_regressive_wedge_cases},
    {"lcontract", "left contraction", "wdg", lcontr_str, ProductType::left_contraction, false, true, ega3d_left_contraction_cases},
    {"rcontract", "right contraction", "wdg", rcontr_str, ProductType::right_contraction, false, true, ega3d_right_contraction_cases},
    {"lexpand", "left expansion", "wdg", lcontr_str, ProductType::left_expansion, false, true, ega3d_left_expansion_cases},
    {"rexpand", "right expansion", "wdg", rcontr_str, ProductType::right_expansion, false, true, ega3d_right_expansion_cases},
    
    // Sandwich products - rotor * object * rev(rotor)
    {"sandwich", "sandwich product", "gpr", mul_str, ProductType::sandwich, false, true, {}}
};

} // namespace configs