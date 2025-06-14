// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_config_ega2d.hpp"
#include "ga_prdxpr_ega2d.hpp"

namespace configs {

////////////////////////////////////////////////////////////////////////////////
// EGA2D Configuration
////////////////////////////////////////////////////////////////////////////////

const AlgebraConfig ega2d_config = {
    .name = "ega2d",
    .description = "Euclidean Geometric Algebra 2D - G(2,0,0)",
    .basis = mv2d_basis,
    .product_rules = {
        {"gpr", gpr_ega2d_rules},
        {"wdg", wdg_ega2d_rules},
        {"dot", dot_ega2d_rules}
        // Note: regressive and contraction rules will be computed from these
    },
    .coefficients = {
        {"A", mv2d_coeff_A},
        {"B", mv2d_coeff_B},
        {"M", mv2d_coeff_M},
        {"M1", mv2d_coeff_M1},
        {"M2", mv2d_coeff_M2},
        {"A_even", mv2d_coeff_A_even},
        {"B_even", mv2d_coeff_B_even},
        {"M_even", mv2d_coeff_M_even},
        {"R_even", mv2d_coeff_R_even},
        {"R_rev_even", mv2d_coeff_R_rev_even},
        {"svps", mv2d_coeff_svps},
        {"svps1", mv2d_coeff_svps1},
        {"svps2", mv2d_coeff_svps2}
    },
    .complement_rules = {
        {"lcmpl", lcmpl_ega2d_rules},
        {"rcmpl", rcmpl_ega2d_rules}
    },
    .filter_names = {"s", "vec", "ps", "mv_e", "mv"},
    .dimension = 2
};

// EGA2D Product-specific case definitions
const std::vector<ProductCase> ega2d_geometric_cases = {
    // These 26 cases match the current generate_and_print_ega2d_gpr() function calls
    {"mv * mv -> mv", "A", "B", "mv", "mv", true},
    {"mv * mv_e -> mv", "A", "B_even", "mv", "mv_e", true}, 
    {"mv_e * mv -> mv", "A_even", "B", "mv_e", "mv", true},
    {"mv * ps -> mv", "A", "svps", "mv", "ps", true},
    {"ps * mv -> mv", "svps", "B", "ps", "mv", true},
    {"mv * vec -> mv", "A", "svps", "mv", "vec", true},
    {"vec * mv -> mv", "svps", "B", "vec", "mv", true},
    {"mv * s -> mv", "A", "svps", "mv", "s", true},
    {"s * mv -> mv", "svps", "B", "s", "mv", true},
    {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e", true},
    {"mv_e * ps -> mv_e", "A_even", "svps", "mv_e", "ps", true},
    {"ps * mv_e -> mv_e", "svps", "B_even", "ps", "mv_e", true},
    {"mv_e * vec -> vec", "A_even", "svps", "mv_e", "vec", true},
    {"vec * mv_e -> vec", "svps", "B_even", "vec", "mv_e", true},
    {"mv_e * s -> mv_e", "A_even", "svps", "mv_e", "s", true},
    {"s * mv_e -> mv_e", "svps", "B_even", "s", "mv_e", true},
    {"ps * ps -> s", "svps1", "svps2", "ps", "ps", true},
    {"ps * vec -> vec", "svps", "svps", "ps", "vec", true},
    {"vec * ps -> vec", "svps", "svps", "vec", "ps", true},
    {"ps * s -> ps", "svps", "svps", "ps", "s", true},
    {"s * ps -> ps", "svps", "svps", "s", "ps", true},
    {"vec * vec -> mv_e", "svps1", "svps2", "vec", "vec", true},
    {"vec * s -> vec", "svps", "svps", "vec", "s", true},
    {"s * vec -> vec", "svps", "svps", "s", "vec", true},
    {"s * s -> s", "svps1", "svps2", "s", "s", true}
};

const std::vector<ProductCase> ega2d_commutator_cases = {
    // 4 cases from reference: cmt(mv,mv), cmt(ps,vec), cmt(vec,ps), cmt(vec,vec)
    {"cmt(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"cmt(ps,vec) -> vec", "svps", "svps", "ps", "vec", true},
    {"cmt(vec,ps) -> vec", "svps", "svps", "vec", "ps", true},
    {"cmt(vec,vec) -> ps", "svps1", "svps2", "vec", "vec", true}    // Fixed: svps1, svps2 for v1.x * v2.y pattern
};

const std::vector<ProductCase> ega2d_wedge_cases = {
    // Wedge product cases matching exact coefficient usage from generate_and_print_ega2d_wdg()
    {"mv ^ mv -> mv", "A", "B", "mv", "mv", true},
    {"mv ^ mv_e -> mv", "A", "B_even", "mv", "mv_e", true}, 
    {"mv_e ^ mv -> mv", "A_even", "B", "mv_e", "mv", true},
    {"mv ^ ps -> ps", "M", "svps", "mv", "ps", true},        // M for mv with ps/vec/s
    {"ps ^ mv -> ps", "svps", "M", "ps", "mv", true},        // M for mv with ps/vec/s
    {"mv ^ vec -> mv", "M", "svps", "mv", "vec", true},      // M for mv with ps/vec/s
    {"vec ^ mv -> mv", "svps", "M", "vec", "mv", true},      // M for mv with ps/vec/s
    {"mv ^ s -> mv", "M", "svps", "mv", "s", true},          // M for mv with ps/vec/s
    {"s ^ mv -> mv", "svps", "M", "s", "mv", true},          // M for mv with ps/vec/s
    {"mv_e ^ mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e", true},
    {"mv_e ^ ps -> ps", "M_even", "svps", "mv_e", "ps", true},     // M_even for mv_e with ps/vec/s
    {"ps ^ mv_e -> ps", "svps", "M_even", "ps", "mv_e", true},     // M_even for mv_e with ps/vec/s
    {"mv_e ^ vec -> vec", "M_even", "svps", "mv_e", "vec", true},  // M_even for mv_e with ps/vec/s
    {"vec ^ mv_e -> vec", "svps", "M_even", "vec", "mv_e", true},  // M_even for mv_e with ps/vec/s
    {"mv_e ^ s -> mv_e", "M_even", "svps", "mv_e", "s", true},     // M_even for mv_e with ps/vec/s
    {"s ^ mv_e -> mv_e", "svps", "M_even", "s", "mv_e", true},     // M_even for mv_e with ps/vec/s
    {"ps ^ ps -> 0", "svps1", "svps2", "ps", "ps", true},
    {"ps ^ vec -> 0", "svps", "svps", "ps", "vec", true},
    {"vec ^ ps -> 0", "svps", "svps", "vec", "ps", true},
    {"ps ^ s -> ps", "svps", "svps", "ps", "s", true},
    {"s ^ ps -> ps", "svps", "svps", "s", "ps", true},
    {"vec ^ vec -> ps", "svps1", "svps2", "vec", "vec", true},
    {"vec ^ s -> vec", "svps", "svps", "vec", "s", true},
    {"s ^ vec -> vec", "svps", "svps", "s", "vec", true},
    {"s ^ s -> s", "svps1", "svps2", "s", "s", true}
};

const std::vector<ProductCase> ega2d_inner_cases = {
    // 4 cases from reference: dot(mv,mv), dot(ps,ps), dot(vec,vec), dot(s,s)
    {"dot(mv,mv) -> s", "A", "B", "mv", "mv", true},
    {"dot(ps,ps) -> s", "svps1", "svps2", "ps", "ps", true},
    {"dot(vec,vec) -> s", "svps1", "svps2", "vec", "vec", true},
    {"dot(s,s) -> s", "svps1", "svps2", "s", "s", true}
};

const std::vector<ProductCase> ega2d_regressive_wedge_cases = {
    // 7 cases from reference: rwdg(mv,mv), rwdg(ps,ps), rwdg(ps,vec), rwdg(vec,ps), rwdg(ps,s), rwdg(s,ps), rwdg(vec,vec)
    {"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv", true},
    {"rwdg(ps,ps) -> ps", "svps1", "svps2", "ps", "ps", true},
    {"rwdg(ps,vec) -> vec", "svps", "svps", "ps", "vec", true},
    {"rwdg(vec,ps) -> vec", "svps", "svps", "vec", "ps", true},
    {"rwdg(ps,s) -> s", "svps", "svps", "ps", "s", true},
    {"rwdg(s,ps) -> s", "svps", "svps", "s", "ps", true},
    {"rwdg(vec,vec) -> s", "svps1", "svps2", "vec", "vec", true}
};

const std::vector<ProductCase> ega2d_left_contraction_cases = {
    // Left contraction cases matching exact coefficient usage from generate_and_print_ega2d_lcontract()
    {"mv << mv -> mv", "A", "B", "mv", "mv", true},
    {"mv << mv_e -> mv", "A", "B_even", "mv", "mv_e", true}, 
    {"mv_e << mv -> mv", "A_even", "B", "mv_e", "mv", true},
    {"mv << ps -> mv", "M", "svps", "mv", "ps", true},        // M for mv with ps/vec/s
    {"ps << mv -> s", "svps", "M", "ps", "mv", true},         // M for mv with ps/vec/s + output type s
    {"mv << vec -> mv", "M", "svps", "mv", "vec", true},      // M for mv with ps/vec/s
    {"vec << mv -> mv", "svps", "M", "vec", "mv", true},      // M for mv with ps/vec/s
    {"mv << s -> s", "M", "svps", "mv", "s", true},           // M for mv with ps/vec/s + output type s
    {"s << mv -> mv", "svps", "M", "s", "mv", true},          // M for mv with ps/vec/s
    {"mv_e << mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e", true},
    {"mv_e << ps -> mv_e", "M_even", "svps", "mv_e", "ps", true},     // M_even for mv_e with ps/vec/s
    {"ps << mv_e -> s", "svps", "M_even", "ps", "mv_e", true},        // M_even for mv_e with ps/vec/s + output type s
    {"mv_e << vec -> vec", "M_even", "svps", "mv_e", "vec", true},    // M_even for mv_e with ps/vec/s
    {"vec << mv_e -> vec", "svps", "M_even", "vec", "mv_e", true},    // M_even for mv_e with ps/vec/s
    {"mv_e << s -> s", "M_even", "svps", "mv_e", "s", true},          // M_even for mv_e with ps/vec/s + output type s
    {"s << mv_e -> mv_e", "svps", "M_even", "s", "mv_e", true},       // M_even for mv_e with ps/vec/s
    {"ps << ps -> s", "svps1", "svps2", "ps", "ps", true},
    {"ps << vec -> 0", "svps", "svps", "ps", "vec", true},
    {"vec << ps -> vec", "svps", "svps", "vec", "ps", true},
    {"ps << s -> 0", "svps", "svps", "ps", "s", true},
    {"s << ps -> ps", "svps", "svps", "s", "ps", true},
    {"vec << vec -> s", "svps1", "svps2", "vec", "vec", true},
    {"vec << s -> 0", "svps", "svps", "vec", "s", true},
    {"s << vec -> vec", "svps", "svps", "s", "vec", true},
    {"s << s -> s", "svps1", "svps2", "s", "s", true}
};

const std::vector<ProductCase> ega2d_right_contraction_cases = {
    // Right contraction cases matching exact coefficient usage from generate_and_print_ega2d_rcontract()
    {"mv >> mv -> mv", "A", "B", "mv", "mv", true},
    {"mv >> mv_e -> mv", "A", "B_even", "mv", "mv_e", true}, 
    {"mv_e >> mv -> mv", "A_even", "B", "mv_e", "mv", true},
    {"mv >> ps -> s", "M", "svps", "mv", "ps", true},         // M for mv with ps/vec/s + output type s
    {"ps >> mv -> mv", "svps", "M", "ps", "mv", true},        // M for mv with ps/vec/s
    {"mv >> vec -> mv", "M", "svps", "mv", "vec", true},      // M for mv with ps/vec/s
    {"vec >> mv -> mv", "svps", "M", "vec", "mv", true},      // M for mv with ps/vec/s
    {"mv >> s -> mv", "M", "svps", "mv", "s", true},          // M for mv with ps/vec/s
    {"s >> mv -> s", "svps", "M", "s", "mv", true},           // M for mv with ps/vec/s + output type s
    {"mv_e >> mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e", true},
    {"mv_e >> ps -> s", "M_even", "svps", "mv_e", "ps", true},        // M_even for mv_e with ps/vec/s + output type s
    {"ps >> mv_e -> mv_e", "svps", "M_even", "ps", "mv_e", true},     // M_even for mv_e with ps/vec/s
    {"mv_e >> vec -> vec", "M_even", "svps", "mv_e", "vec", true},    // M_even for mv_e with ps/vec/s
    {"vec >> mv_e -> vec", "svps", "M_even", "vec", "mv_e", true},    // M_even for mv_e with ps/vec/s
    {"mv_e >> s -> mv_e", "M_even", "svps", "mv_e", "s", true},       // M_even for mv_e with ps/vec/s
    {"s >> mv_e -> s", "svps", "M_even", "s", "mv_e", true},          // M_even for mv_e with ps/vec/s + output type s
    {"ps >> ps -> s", "svps1", "svps2", "ps", "ps", true},
    {"ps >> vec -> vec", "svps", "svps", "ps", "vec", true},
    {"vec >> ps -> 0", "svps", "svps", "vec", "ps", true},
    {"ps >> s -> ps", "svps", "svps", "ps", "s", true},
    {"s >> ps -> 0", "svps", "svps", "s", "ps", true},
    {"vec >> vec -> s", "svps1", "svps2", "vec", "vec", true},
    {"vec >> s -> vec", "svps", "svps", "vec", "s", true},
    {"s >> vec -> 0", "svps", "svps", "s", "vec", true},
    {"s >> s -> s", "svps1", "svps2", "s", "s", true}
};

// EGA2D Products with their specific cases
const std::vector<ProductDefinition> ega2d_products = {
    {"gpr", "geometric product", "gpr", mul_str, ProductType::geometric_product, false, true, ega2d_geometric_cases},
    {"cmt", "commutator product", "gpr", mul_str, ProductType::commutator_product, false, true, ega2d_commutator_cases},
    {"wdg", "wedge product", "wdg", wdg_str, ProductType::wedge_product, false, true, ega2d_wedge_cases},
    {"dot", "inner product", "dot", "dot", ProductType::inner_product, false, true, ega2d_inner_cases},
    {"rwdg", "regressive wedge product", "wdg", "rwdg", ProductType::regressive_wedge, false, true, ega2d_regressive_wedge_cases},
    {"lcontract", "left contraction", "wdg", lcontr_str, ProductType::left_contraction, false, true, ega2d_left_contraction_cases},
    {"rcontract", "right contraction", "wdg", rcontr_str, ProductType::right_contraction, false, true, ega2d_right_contraction_cases},
    
    // Expansions show basis product tables but no specific cases for ega2d
    {"lexpand", "left expansion", "wdg", lcontr_str, ProductType::left_expansion, false, true, {}},
    {"rexpand", "right expansion", "wdg", rcontr_str, ProductType::right_expansion, false, true, {}},
    
    // Sandwich products - rotor * object * rev(rotor)
    {"sandwich", "sandwich product", "gpr", mul_str, ProductType::sandwich, false, true, {}}
};

} // namespace configs