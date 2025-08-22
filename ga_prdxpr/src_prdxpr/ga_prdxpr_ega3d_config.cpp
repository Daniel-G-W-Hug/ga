// Copyright 2024-2025, Daniel Hug. All rights reserved.

/*
 * GA PRODUCT RULE GENERATION SYSTEM GUIDE
 * =======================================
 *
 * This configuration file defines two types of product rule generation:
 *
 * 1. AUTOMATIC GENERATION (AlgebraConfig): Complement rules generated mathematically
 * 2. MANUAL CASES (ProductConfig): Expression cases configured manually
 *
 * AUTOMATIC RULE GENERATION:
 * --------------------------
 * The get_ega3d_algebra_config() function provides mathematical algebra parameters:
 * - basis_vectors: Basis vector names (e.g., {"e1", "e2", "e3"})
 * - metric_signature: Quadratic form values (e.g., {+1, +1, +1})
 * - multivector_basis: Complete basis element names (e.g., {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e123"})
 * - scalar_name: Name for scalar element (typically "1")
 * - basis_prefix: Prefix for basis elements (typically "e")
 *
 * From this configuration, the system automatically generates:
 * - Geometric product rules (geometric multiplication with metric)
 * - Wedge product rules (antisymmetric exterior product)
 * - Dot product rules (symmetric contraction with extended metric)
 * - Complement rules (computed from wedge product tables)
 *
 * COMPLEMENT RULE GENERATION:
 * The complement rules are now AUTOMATICALLY GENERATED using mathematical algorithms:
 * - Even algebras (EGA2D, PGA3DP): Generate left_complement and right_complement
 * - Odd algebras (EGA3D, PGA2DP): Generate single complement
 * - Algorithm: Search wedge product table for pseudoscalar relationships
 * - Result: Mathematically verified complement rules with zero transcription errors
 *
 * MANUAL CASE CONFIGURATION:
 * --------------------------
 * Each ProductConfig contains a .cases array with mathematical operation descriptions.
 * Cases define which coefficient combinations and type filters to use for code
 * generation.
 *
 * STANDARD FORMAT (5 parameters):
 * {"case_name", "left_coeff", "right_coeff", "left_filter", "right_filter"}
 *
 * SANDWICH FORMAT (7 parameters):
 * {"case_name", "left_coeff", "right_coeff", "left_filter", "right_filter", is_two_step,
 * "intermediate_name"}
 *
 * PARAMETER DESCRIPTIONS:
 * ----------------------
 * 1. case_name:         Mathematical operation description
 *    Examples: "mv * mv -> mv", "vec ^ s -> vec", "dot(bivec,bivec) -> s"
 *    Format: "left_type OPERATOR right_type -> result_type"
 *    Special: Use function notation for non-infix operators: "dot(A,B) -> result"
 *
 * 2. left_coeff_name:   Coefficient name for left operand (see COEFFICIENT SYSTEM below)
 * 3. right_coeff_name:  Coefficient name for right operand (see COEFFICIENT SYSTEM below)
 * 4. left_filter_name:  Type filter for left operand (see AVAILABLE FILTERS below)
 * 5. right_filter_name: Type filter for right operand (see AVAILABLE FILTERS below)
 *
 * SANDWICH PRODUCTS (7 parameters):
 * 6. is_two_step:       Always true for sandwich products
 * 7. intermediate_name: Name for intermediate result (e.g., "vec_tmp")
 *
 * COEFFICIENT SYSTEM LINKAGE:
 * ---------------------------
 * Coefficient names in config files map to mvec_coeff definitions in corresponding header
 * files. The mapping is established in create_[algebra]_algebra_data() via
 * AlgebraData.coefficients:
 *
 * Config Name -> Header Definition:
 * "A"         -> mv2d_coeff_A / mv2dp_coeff_A / mv3d_coeff_A / mv3dp_coeff_A
 * "B"         -> mv2d_coeff_B / mv2dp_coeff_B / mv3d_coeff_B / mv3dp_coeff_B
 * "A_even"    -> mv2d_coeff_A_even / mv2dp_coeff_A_even / etc.
 * "R_even"    -> mv2d_coeff_R_even / mv2dp_coeff_R_even / etc.
 * "svps"      -> mv2d_coeff_svps (EGA2D/3D naming)
 * "svBps"     -> mv2dp_coeff_svBps (PGA2DP/3DP naming)
 *
 * Each mvec_coeff contains component strings matching the algebra's basis:
 * - EGA2D: {"1", "e1", "e2", "e12"} -> {"A.c0", "A.c1", "A.c2", "A.c3"}
 * - EGA3D: {"1", "e1", "e2", "e3", "e23", "e31", "e12", "e321"} -> {"A.c0"..."A.c7"}
 *
 * TO EXTEND THE SYSTEM:
 * 1. Define new mvec_coeff in the algebra header file (e.g., ga_prdxpr_ega3d.hpp)
 * 2. Add mapping in create_[algebra]_algebra_data() coefficients map
 * 3. Use the coefficient name in .cases arrays
 *
 * Example extension for new coefficient "C":
 * Header: const mvec_coeff mv3d_coeff_C = {"C.c0", "C.c1", "C.c2", "C.c3", "C.c4",
 * "C.c5", "C.c6", "C.c7"}; Config: ega3d.coefficients = {..."C", mv3d_coeff_C}... Usage:
 * {"mv * mv -> mv", "C", "B", "mv", "mv"}
 *
 * AVAILABLE COEFFICIENTS (EGA3D):
 * - General: A, B, M, M1, M2 (full multivectors)
 * - Even/Odd: A_even, B_even, A_odd, B_odd, M_even, M_odd (grade-filtered)
 * - Rotors: R_even, R_rev_even (even elements for rotation)
 * - Symmetric: svBps, svBps1, svBps2 (scalar, vector, bivector, pseudoscalar patterns)
 *
 * AVAILABLE FILTERS (EGA3D):
 * - s (scalar), vec (vector), bivec (bivector), ps (pseudoscalar)
 * - mv_e (even), mv_u (odd), mv (multivector)
 *
 * NAMING CONVENTIONS:
 * - Use lowercase for all type names: s, vec, bivec, ps, mv, mv_e, mv_u
 * - Result types should match geometric algebra conventions
 * - Zero results: Use "0" as result type for operations that yield zero
 */

#include "ga_prdxpr_ega3d_config.hpp"

// Automatic rule generation configuration for EGA3D
AlgebraConfig get_ega3d_algebra_config()
{
    return {.basis_vectors = {"e1", "e2", "e3"},
            .metric_signature = {+1, +1, +1},
            .multivector_basis = {"1", "e1", "e2", "e3", "e23", "e31", "e12",
                                  "e123"}, // Copy from mv3d_basis
            .scalar_name = one_str(),
            .basis_prefix = "e"};
}

////////////////////////////////////////////////////////////////////////////////
// ALGEBRA CONFIGURATION - EGA3D
//
// This file contains pure configuration for the EGA3D geometric algebra.
// Product rule generation is handled in ga_prdxpr_ega3d.cpp.
//
// Configuration includes:
// - Algebra parameters (basis vectors, metric signature)
// - Coefficient definitions
// - Product case specifications
////////////////////////////////////////////////////////////////////////////////

namespace configurable {

AlgebraData create_ega3d_algebra_data()
{
    AlgebraData ega3d;
    ega3d.name = "ega3d";
    ega3d.basis = mv3d_basis; // Use existing basis
    ega3d.dimension = 3;

    // Map coefficient names to existing coefficient objects from header
    ega3d.coefficients = {{"A", mv3d_coeff_A},
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
                          {"svBps2", mv3d_coeff_svBps2}};

    // Map filter names to existing filter enums
    ega3d.filters_3d = {{"s", filter_3d::s},         {"vec", filter_3d::vec},
                        {"bivec", filter_3d::bivec}, {"ps", filter_3d::ps},
                        {"mv_e", filter_3d::mv_e},   {"mv_u", filter_3d::mv_u},
                        {"mv", filter_3d::mv}};

    return ega3d;
}

ProductConfig get_ega3d_gpr_config()
{
    return {.product_name = "gpr",
            .description = "ega3d geometric product",
            .display_name = "geometric product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv * mv -> mv", "A", "B", "mv", "mv"},
                      {"mv * mv_e -> mv", "A", "B_even", "mv", "mv_e"},
                      {"mv_e * mv -> mv", "A_even", "B", "mv_e", "mv"},
                      {"mv * mv_u -> mv", "A", "B_odd", "mv", "mv_u"},
                      {"mv_u * mv -> mv", "A_odd", "B", "mv_u", "mv"},
                      {"mv * ps -> mv", "A", "svBps", "mv", "ps"},
                      {"ps * mv -> mv", "svBps", "B", "ps", "mv"},
                      {"mv * bivec -> mv", "M", "svBps", "mv", "bivec"},
                      {"bivec * mv -> mv", "svBps", "M", "bivec", "mv"},
                      {"mv * vec -> mv", "M", "svBps", "mv", "vec"},
                      {"vec * mv -> mv", "svBps", "M", "vec", "mv"},
                      {"mv * s -> mv", "M", "svBps", "mv", "s"},
                      {"s * mv -> mv", "svBps", "M", "s", "mv"},
                      {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e"},
                      {"mv_u * mv_u -> mv_e", "A_odd", "B_odd", "mv_u", "mv_u"},
                      {"mv_e * mv_u -> mv_u", "A_even", "B_odd", "mv_e", "mv_u"},
                      {"mv_u * mv_e -> mv_u", "A_odd", "B_even", "mv_u", "mv_e"},
                      {"mv_e * ps -> mv_u", "M_even", "svBps", "mv_e", "ps"},
                      {"ps * mv_e -> mv_u", "svBps", "M_even", "ps", "mv_e"},
                      {"mv_u * ps -> mv_e", "M_odd", "svBps", "mv_u", "ps"},
                      {"ps * mv_u -> mv_e", "svBps", "M_odd", "ps", "mv_u"},
                      {"mv_e * bivec -> mv_e", "M_even", "svBps", "mv_e", "bivec"},
                      {"bivec * mv_e -> mv_e", "svBps", "M_even", "bivec", "mv_e"},
                      {"mv_u * bivec -> mv_u", "M_odd", "svBps", "mv_u", "bivec"},
                      {"bivec * mv_u -> mv_u", "svBps", "M_odd", "bivec", "mv_u"},
                      {"mv_e * vec -> mv_u", "M_even", "svBps", "mv_e", "vec"},
                      {"vec * mv_e -> mv_u", "svBps", "M_even", "vec", "mv_e"},
                      {"mv_u * vec -> mv_e", "M_odd", "svBps", "mv_u", "vec"},
                      {"vec * mv_u -> mv_e", "svBps", "M_odd", "vec", "mv_u"},
                      {"mv_e * s -> mv_e", "M_even", "svBps", "mv_e", "s"},
                      {"s * mv_e -> mv_e", "svBps", "M_even", "s", "mv_e"},
                      {"mv_u * s -> mv_u", "M_odd", "svBps", "mv_u", "s"},
                      {"s * mv_u -> mv_u", "svBps", "M_odd", "s", "mv_u"},
                      {"ps * ps -> s", "svBps1", "svBps2", "ps", "ps"},
                      {"ps * bivec -> vec", "svBps", "svBps", "ps", "bivec"},
                      {"bivec * ps -> vec", "svBps", "svBps", "bivec", "ps"},
                      {"ps * vec -> bivec", "svBps", "svBps", "ps", "vec"},
                      {"vec * ps -> bivec", "svBps", "svBps", "vec", "ps"},
                      {"ps * s -> ps", "svBps", "svBps", "ps", "s"},
                      {"s * ps -> ps", "svBps", "svBps", "s", "ps"},
                      {"bivec * bivec -> mv_e", "svBps1", "svBps2", "bivec", "bivec"},
                      {"bivec * vec -> mv_u", "svBps", "svBps", "bivec", "vec"},
                      {"vec * bivec -> mv_u", "svBps", "svBps", "vec", "bivec"},
                      {"bivec * s -> bivec", "svBps", "svBps", "bivec", "s"},
                      {"s * bivec -> bivec", "svBps", "svBps", "s", "bivec"},
                      {"vec * vec -> mv_e", "svBps1", "svBps2", "vec", "vec"},
                      {"vec * s -> vec", "svBps", "svBps", "vec", "s"},
                      {"s * vec -> vec", "svBps", "svBps", "s", "vec"},
                      {"s * s -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_cmt_config()
{
    return {.product_name = "cmt",
            .description = "ega3d commutator product",
            .display_name = "commutator product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"cmt(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"cmt(bivec,bivec) -> bivec", "svBps1", "svBps2", "bivec", "bivec"},
                      {"cmt(bivec,vec) -> vec", "svBps", "svBps", "bivec", "vec"},
                      {"cmt(vec,bivec) -> vec", "svBps", "svBps", "vec", "bivec"},
                      {"cmt(vec,vec) -> bivec", "svBps1", "svBps2", "vec", "vec"}}};
}

ProductConfig get_ega3d_wdg_config()
{
    return {.product_name = "wdg",
            .description = "ega3d wedge product",
            .display_name = "wedge product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv ^ mv -> mv", "A", "B", "mv", "mv"},
                      {"mv ^ bivec -> mv", "M", "svBps", "mv", "bivec"},
                      {"bivec ^ mv -> mv", "svBps", "M", "bivec", "mv"},
                      {"mv ^ vec -> mv", "M", "svBps", "mv", "vec"},
                      {"vec ^ mv -> mv", "svBps", "M", "vec", "mv"},
                      {"mv ^ s -> mv", "M", "svBps", "mv", "s"},
                      {"s ^ mv -> mv", "svBps", "M", "s", "mv"},
                      {"ps ^ ps -> 0", "svBps1", "svBps2", "ps", "ps"},
                      {"ps ^ bivec -> 0", "svBps", "svBps", "ps", "bivec"},
                      {"bivec ^ ps -> 0", "svBps", "svBps", "bivec", "ps"},
                      {"ps ^ vec -> 0", "svBps", "svBps", "ps", "vec"},
                      {"vec ^ ps -> 0", "svBps", "svBps", "vec", "ps"},
                      {"ps ^ s -> ps", "svBps", "svBps", "ps", "s"},
                      {"s ^ ps -> ps", "svBps", "svBps", "s", "ps"},
                      {"bivec ^ bivec -> 0", "svBps1", "svBps2", "bivec", "bivec"},
                      {"bivec ^ vec -> ps", "svBps", "svBps", "bivec", "vec"},
                      {"vec ^ bivec -> ps", "svBps", "svBps", "vec", "bivec"},
                      {"bivec ^ s -> bivec", "svBps", "svBps", "bivec", "s"},
                      {"s ^ bivec -> bivec", "svBps", "svBps", "s", "bivec"},
                      {"vec ^ vec -> bivec", "svBps1", "svBps2", "vec", "vec"},
                      {"vec ^ s -> vec", "svBps", "svBps", "vec", "s"},
                      {"s ^ vec -> vec", "svBps", "svBps", "s", "vec"},
                      {"s ^ s -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_dot_config()
{
    return {.product_name = "dot",
            .description = "ega3d inner product",
            .display_name = "inner product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"dot(mv,mv) -> s", "A", "B", "mv", "mv"},
                      {"dot(ps,ps) -> s", "svBps1", "svBps2", "ps", "ps"},
                      {"dot(bivec,bivec) -> s", "svBps1", "svBps2", "bivec", "bivec"},
                      {"dot(vec,vec) -> s", "svBps1", "svBps2", "vec", "vec"},
                      {"dot(s,s) -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_rwdg_config()
{
    return {.product_name = "rwdg",
            .description = "ega3d regressive wedge product",
            .display_name = "regressive wedge product",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"rwdg(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"rwdg(ps,ps) -> ps", "svBps1", "svBps2", "ps", "ps"},
                      {"rwdg(ps,bivec) -> bivec", "svBps", "svBps", "ps", "bivec"},
                      {"rwdg(bivec,ps) -> bivec", "svBps", "svBps", "bivec", "ps"},
                      {"rwdg(ps,vec) -> vec", "svBps", "svBps", "ps", "vec"},
                      {"rwdg(vec,ps) -> vec", "svBps", "svBps", "vec", "ps"},
                      {"rwdg(ps,s) -> s", "svBps", "svBps", "ps", "s"},
                      {"rwdg(s,ps) -> s", "svBps", "svBps", "s", "ps"},
                      {"rwdg(bivec,bivec) -> vec", "svBps1", "svBps2", "bivec", "bivec"},
                      {"rwdg(bivec,vec) -> s", "svBps", "svBps", "bivec", "vec"},
                      {"rwdg(vec,bivec) -> s", "svBps", "svBps", "vec", "bivec"},
                      {"rwdg(bivec,s) -> 0", "svBps", "svBps", "bivec", "s"},
                      {"rwdg(s,bivec) -> 0", "svBps", "svBps", "s", "bivec"},
                      {"rwdg(vec,vec) -> 0", "svBps1", "svBps2", "vec", "vec"},
                      {"rwdg(vec,s) -> 0", "svBps", "svBps", "vec", "s"},
                      {"rwdg(s,vec) -> 0", "svBps", "svBps", "s", "vec"},
                      {"rwdg(s,s) -> 0", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_lcontract_config()
{
    return {.product_name = "lcontract",
            .description = "ega3d left contraction",
            .display_name = "left contraction",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv << mv -> mv", "A", "B", "mv", "mv"},
                      {"bivec << mv -> mv", "svBps", "M", "bivec", "mv"},
                      {"vec << mv -> mv", "svBps", "M", "vec", "mv"},
                      {"s << mv -> mv", "svBps", "M", "s", "mv"},
                      {"s << mv_e -> mv_e", "svBps", "M_even", "s", "mv_e"},
                      {"s << mv_u -> mv_u", "svBps", "M_odd", "s", "mv_u"},
                      {"ps << ps -> s", "svBps1", "svBps2", "ps", "ps"},
                      {"ps << bivec -> 0", "svBps", "svBps", "ps", "bivec"},
                      {"bivec << ps -> vec", "svBps", "svBps", "bivec", "ps"},
                      {"ps << vec -> 0", "svBps", "svBps", "ps", "vec"},
                      {"vec << ps -> bivec", "svBps", "svBps", "vec", "ps"},
                      {"ps << s -> 0", "svBps", "svBps", "ps", "s"},
                      {"s << ps -> ps", "svBps", "svBps", "s", "ps"},
                      {"bivec << bivec -> s", "svBps1", "svBps2", "bivec", "bivec"},
                      {"bivec << vec -> 0", "svBps", "svBps", "bivec", "vec"},
                      {"vec << bivec -> vec", "svBps", "svBps", "vec", "bivec"},
                      {"bivec << s -> 0", "svBps", "svBps", "bivec", "s"},
                      {"s << bivec -> bivec", "svBps", "svBps", "s", "bivec"},
                      {"vec << vec -> s", "svBps1", "svBps2", "vec", "vec"},
                      {"vec << s -> 0", "svBps", "svBps", "vec", "s"},
                      {"s << vec -> vec", "svBps", "svBps", "s", "vec"},
                      {"s << s -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_rcontract_config()
{
    return {.product_name = "rcontract",
            .description = "ega3d right contraction",
            .display_name = "right contraction",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv >> mv -> mv", "A", "B", "mv", "mv"},
                      {"mv >> bivec -> mv", "M", "svBps", "mv", "bivec"},
                      {"mv >> vec -> mv", "M", "svBps", "mv", "vec"},
                      {"mv >> s -> mv", "M", "svBps", "mv", "s"},
                      {"mv_e >> s -> mv_e", "M_even", "svBps", "mv_e", "s"},
                      {"mv_u >> s -> mv_u", "M_odd", "svBps", "mv_u", "s"},
                      {"ps >> ps -> s", "svBps1", "svBps2", "ps", "ps"},
                      {"ps >> bivec -> vec", "svBps", "svBps", "ps", "bivec"},
                      {"bivec >> ps -> 0", "svBps", "svBps", "bivec", "ps"},
                      {"ps >> vec -> bivec", "svBps", "svBps", "ps", "vec"},
                      {"vec >> ps -> 0", "svBps", "svBps", "vec", "ps"},
                      {"ps >> s -> ps", "svBps", "svBps", "ps", "s"},
                      {"s >> ps -> 0", "svBps", "svBps", "s", "ps"},
                      {"bivec >> bivec -> s", "svBps1", "svBps2", "bivec", "bivec"},
                      {"bivec >> vec -> vec", "svBps", "svBps", "bivec", "vec"},
                      {"vec >> bivec -> 0", "svBps", "svBps", "vec", "bivec"},
                      {"bivec >> s -> bivec", "svBps", "svBps", "bivec", "s"},
                      {"s >> bivec -> 0", "svBps", "svBps", "s", "bivec"},
                      {"vec >> vec -> s", "svBps1", "svBps2", "vec", "vec"},
                      {"vec >> s -> vec", "svBps", "svBps", "vec", "s"},
                      {"s >> vec -> 0", "svBps", "svBps", "s", "vec"},
                      {"s >> s -> s", "svBps1", "svBps2", "s", "s"}}};
}

ProductConfig get_ega3d_lexpand_config()
{
    return {.product_name = "lexpand",
            .description = "ega3d left expansion",
            .display_name = "left expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"lexpand(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"lexpand(bivec,vec) -> bivec", "svBps", "svBps", "bivec", "vec"}}};
}

ProductConfig get_ega3d_rexpand_config()
{
    return {.product_name = "rexpand",
            .description = "ega3d right expansion",
            .display_name = "right expansion",
            // Format: {"operation(A,B) -> result", "left_coeff", "right_coeff",
            // "left_filter", "right_filter"}
            .cases = {{"rexpand(mv,mv) -> mv", "A", "B", "mv", "mv"},
                      {"rexpand(vec,bivec) -> bivec", "svBps", "svBps", "vec", "bivec"}}};
}

ProductConfig get_ega3d_sandwich_gpr_config()
{
    return {.product_name = "sandwich_gpr",
            .description = "ega3d rotor geometric product",
            .display_name = "sandwich product",
            // Format: {"case_name", "left_coeff", "right_coeff", "left_filter",
            // "right_filter"}
            .cases = {{"mv_e * vec -> mv_u", "R_even", "svBps", "mv_e", "vec"},
                      {"mv_e * bivec -> mv_e", "R_even", "svBps", "mv_e", "bivec"}},
            .is_sandwich_product = true};
}

} // namespace configurable