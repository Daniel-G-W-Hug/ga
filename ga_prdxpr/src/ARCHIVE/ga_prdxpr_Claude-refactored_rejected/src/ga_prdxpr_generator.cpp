// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_generator.hpp"
#include "ga_prdxpr_common.hpp"
#include "ga_prdxpr_config.hpp"

////////////////////////////////////////////////////////////////////////////////
// ProductExpressionGenerator implementation
////////////////////////////////////////////////////////////////////////////////

void ProductExpressionGenerator::generate_all_algebras()
{
    generate_ega2d_expressions();
    generate_ega3d_expressions();
    generate_pga2dp_expressions();
    generate_pga3dp_expressions();
}

void ProductExpressionGenerator::generate_algebra(
    const AlgebraConfig& config, const std::vector<ProductDefinition>& products)
{
    // TODO: Re-enable structural elements later when output matches reference 100%
    // fmt::println("===============================================================================");
    // fmt::println("Generating expressions for: {} - {}", config.name,
    // config.description);
    // fmt::println("===============================================================================\n");

    // Generate each enabled product type
    for (const auto& product : products) {
        if (product.enabled) {
            generate_product_type(config, product);
        }
    }

    // TODO: Re-enable structural elements later
    // fmt::println("===============================================================================");
    // fmt::println("Completed: {}", config.name);
    // fmt::println("===============================================================================\n\n");
}

void ProductExpressionGenerator::generate_product_type(
    const AlgebraConfig& config, const ProductDefinition& product)
{
    try {
        // Get the basis product table for this product type
        auto basis_tab = get_basis_product_table(config, product);

        // Print header for this product type
        std::string prd_name = config.name + " " + product.display_name;
        fmt::println("{} - basis product table:", prd_name);
        print_prd_tab(basis_tab);
        fmt::println("");

        // Generate each enabled case for this product using product's own cases
        for (const auto& case_def : product.cases) {
            if (case_def.enabled) {
                generate_product_case(config, product, case_def, basis_tab);
            }
        }

        // Generate sandwich products after all cases
        if (product.type == ProductType::sandwich || product.type == ProductType::regressive_sandwich) {
            generate_sandwich_products(config, product, basis_tab);
        }

        fmt::println(
            "-------------------------------------------------------------------\n");
    }
    catch (const std::exception& e) {
        fmt::println("Error generating {} {}: {}", config.name, product.display_name,
                     e.what());
    }
}

void ProductExpressionGenerator::generate_product_case(const AlgebraConfig& config,
                                                       const ProductDefinition& product,
                                                       const ProductCase& case_def,
                                                       const prd_table& basis_tab)
{
    // Dispatch to dimension-specific implementation
    switch (config.dimension) {
        case 2:
            generate_case_2d(config, product, case_def, basis_tab);
            break;
        case 3:
            generate_case_3d(config, product, case_def, basis_tab);
            break;
        case 4:
            generate_case_4d(config, product, case_def, basis_tab);
            break;
        default:
            throw std::invalid_argument("Unsupported dimension: " +
                                        std::to_string(config.dimension));
    }
}

void ProductExpressionGenerator::generate_case_2d(const AlgebraConfig& config,
                                                  const ProductDefinition& product,
                                                  const ProductCase& case_def,
                                                  const prd_table& basis_tab)
{
    // Look up coefficients
    auto left_coeff_it = config.coefficients.find(case_def.left_coeff_key);
    auto right_coeff_it = config.coefficients.find(case_def.right_coeff_key);

    if (left_coeff_it == config.coefficients.end()) {
        throw std::invalid_argument("Unknown left coefficient: " +
                                    case_def.left_coeff_key);
    }
    if (right_coeff_it == config.coefficients.end()) {
        throw std::invalid_argument("Unknown right coefficient: " +
                                    case_def.right_coeff_key);
    }

    const auto& left_coeff = left_coeff_it->second;
    const auto& right_coeff = right_coeff_it->second;

    // Convert filter strings to enums
    auto left_filter = configs::get_filter_2d_from_string(case_def.left_filter);
    auto right_filter = configs::get_filter_2d_from_string(case_def.right_filter);

    // Generate the product
    fmt::println("{}:",
                 config.name + " " + product.display_name + " " + case_def.description);
    auto prd_tab = get_prd_tab(basis_tab, left_coeff, right_coeff);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, config.basis, left_filter, right_filter);
    print_mvec(prd_mv, config.basis);
    fmt::println("");

    // Print comment if available
    if (!case_def.comment.empty()) {
        fmt::println("  // {}", case_def.comment);
    }
}

void ProductExpressionGenerator::generate_case_3d(const AlgebraConfig& config,
                                                  const ProductDefinition& product,
                                                  const ProductCase& case_def,
                                                  const prd_table& basis_tab)
{
    // Look up coefficients
    auto left_coeff_it = config.coefficients.find(case_def.left_coeff_key);
    auto right_coeff_it = config.coefficients.find(case_def.right_coeff_key);

    if (left_coeff_it == config.coefficients.end()) {
        throw std::invalid_argument("Unknown left coefficient: " +
                                    case_def.left_coeff_key);
    }
    if (right_coeff_it == config.coefficients.end()) {
        throw std::invalid_argument("Unknown right coefficient: " +
                                    case_def.right_coeff_key);
    }

    const auto& left_coeff = left_coeff_it->second;
    const auto& right_coeff = right_coeff_it->second;

    // Convert filter strings to enums
    auto left_filter = configs::get_filter_3d_from_string(case_def.left_filter);
    auto right_filter = configs::get_filter_3d_from_string(case_def.right_filter);

    // Generate the product
    fmt::println("{}:",
                 config.name + " " + product.display_name + " " + case_def.description);
    auto prd_tab = get_prd_tab(basis_tab, left_coeff, right_coeff);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, config.basis, left_filter, right_filter);
    print_mvec(prd_mv, config.basis);
    fmt::println("");

    // Print comment if available
    if (!case_def.comment.empty()) {
        fmt::println("  // {}", case_def.comment);
    }
}

void ProductExpressionGenerator::generate_case_4d(const AlgebraConfig& config,
                                                  const ProductDefinition& product,
                                                  const ProductCase& case_def,
                                                  const prd_table& basis_tab)
{
    // Look up coefficients
    auto left_coeff_it = config.coefficients.find(case_def.left_coeff_key);
    auto right_coeff_it = config.coefficients.find(case_def.right_coeff_key);

    if (left_coeff_it == config.coefficients.end()) {
        throw std::invalid_argument("Unknown left coefficient: " +
                                    case_def.left_coeff_key);
    }
    if (right_coeff_it == config.coefficients.end()) {
        throw std::invalid_argument("Unknown right coefficient: " +
                                    case_def.right_coeff_key);
    }

    const auto& left_coeff = left_coeff_it->second;
    const auto& right_coeff = right_coeff_it->second;

    // Convert filter strings to enums
    auto left_filter = configs::get_filter_4d_from_string(case_def.left_filter);
    auto right_filter = configs::get_filter_4d_from_string(case_def.right_filter);

    // Generate the product
    fmt::println("{}:",
                 config.name + " " + product.display_name + " " + case_def.description);
    auto prd_tab = get_prd_tab(basis_tab, left_coeff, right_coeff);
    auto prd_mv = get_mv_from_prd_tab(prd_tab, config.basis, left_filter, right_filter);
    print_mvec(prd_mv, config.basis);
    fmt::println("");

    // Print comment if available
    if (!case_def.comment.empty()) {
        fmt::println("  // {}", case_def.comment);
    }
}

prd_table
ProductExpressionGenerator::get_basis_product_table(const AlgebraConfig& config,
                                                    const ProductDefinition& product)
{
    // Find the underlying product rules (usually wedge for computed products)
    auto rules_it = config.product_rules.find(product.rule_key);
    if (rules_it == config.product_rules.end()) {
        throw std::invalid_argument("Unknown product rule: " + product.rule_key);
    }

    const auto& rules = rules_it->second;

    // Create basis product table based on product type
    prd_table basis_tab;

    switch (product.type) {
        case ProductType::geometric_product: {
            // Geometric product: apply gpr rules to (* operator)
            basis_tab =
                apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(config.basis, config.basis, mul_str),
                                   rules);
            break;
        }
        
        case ProductType::commutator_product: {
            // Commutator product: geometric product basis table + asymmetric part
            basis_tab =
                apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(config.basis, config.basis, mul_str),
                                   rules);
            basis_tab = get_prd_tab_asym(basis_tab);  // use the asymmetric part only
            break;
        }
        
        case ProductType::wedge_product: {
            // Wedge product: apply wdg rules to (^ operator)
            basis_tab =
                apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(config.basis, config.basis, wdg_str),
                                   rules);
            break;
        }
        
        case ProductType::inner_product: {
            // Inner product: apply dot rules to (* operator)  
            basis_tab =
                apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(config.basis, config.basis, mul_str),
                                   rules);
            break;
        }

        case ProductType::regressive_wedge: {
            // rwdg(A,B) = lcmpl(wdg(rcmpl(A), rcmpl(B)))
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");

            if (rcmpl_rules_it != config.complement_rules.end() &&
                lcmpl_rules_it != config.complement_rules.end()) {

                // Transform basis elements: rcmpl(A) and rcmpl(B)
                auto rcmpl_left = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);
                auto rcmpl_right =
                    apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // Compute wedge product of complements
                auto transformed_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(rcmpl_left, rcmpl_right, wdg_str), rules);

                // Apply left complement to result: lcmpl(result)
                basis_tab = apply_rules_to_tab(transformed_tab, lcmpl_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Complement rules not found for regressive wedge");
            }
            break;
        }

        case ProductType::left_contraction: {
            // A << B = rwdg(lcmpl(A), B)
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");

            if (rcmpl_rules_it != config.complement_rules.end() &&
                lcmpl_rules_it != config.complement_rules.end()) {

                // lhs: apply lcmpl then rcmpl to A: rcmpl(lcmpl(A))
                auto lhs = apply_rules_to_mv(
                    apply_rules_to_mv(config.basis, lcmpl_rules_it->second),
                    rcmpl_rules_it->second);

                // rhs: apply rcmpl to B: rcmpl(B)
                auto rhs = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // Compute wedge product
                auto transformed_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), rules);

                // Apply left complement to result
                basis_tab = apply_rules_to_tab(transformed_tab, lcmpl_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Complement rules not found for left contraction");
            }
            break;
        }

        case ProductType::right_contraction: {
            // A >> B = rwdg(A, rcmpl(B))
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");

            if (rcmpl_rules_it != config.complement_rules.end() &&
                lcmpl_rules_it != config.complement_rules.end()) {

                // lhs: apply rcmpl to A: rcmpl(A)
                auto lhs = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // rhs: apply rcmpl to rcmpl(B): rcmpl(rcmpl(B))
                auto rhs = apply_rules_to_mv(
                    apply_rules_to_mv(config.basis, rcmpl_rules_it->second),
                    rcmpl_rules_it->second);

                // Compute wedge product
                auto transformed_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), rules);

                // Apply left complement to result
                basis_tab = apply_rules_to_tab(transformed_tab, lcmpl_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Complement rules not found for right contraction");
            }
            break;
        }

        case ProductType::left_expansion: {
            // lexpand(A,B) = wdg(lcmpl(A), B)
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");

            if (lcmpl_rules_it != config.complement_rules.end()) {

                // lhs: apply lcmpl to A: lcmpl(A)
                auto lhs = apply_rules_to_mv(config.basis, lcmpl_rules_it->second);

                // rhs: B unchanged
                auto rhs = config.basis;

                // Compute wedge product directly (no final complement transform)
                basis_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), rules);
            }
            else {
                throw std::invalid_argument(
                    "Complement rules not found for left expansion");
            }
            break;
        }

        case ProductType::right_expansion: {
            // rexpand(A,B) = wdg(A, rcmpl(B))
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");

            if (rcmpl_rules_it != config.complement_rules.end()) {

                // lhs: A unchanged
                auto lhs = config.basis;

                // rhs: apply rcmpl to B: rcmpl(B)
                auto rhs = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // Compute wedge product directly (no final complement transform)
                basis_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str), rules);
            }
            else {
                throw std::invalid_argument(
                    "Complement rules not found for right expansion");
            }
            break;
        }

        case ProductType::regressive_inner: {
            // rdot(A,B) = lcmpl(dot(rcmpl(A), rcmpl(B)))
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");

            if (rcmpl_rules_it != config.complement_rules.end() &&
                lcmpl_rules_it != config.complement_rules.end()) {
                // Transform basis elements: rcmpl(A) and rcmpl(B)
                auto rcmpl_left = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);
                auto rcmpl_right = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // Compute inner product of complements
                auto transformed_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(rcmpl_left, rcmpl_right, mul_str), rules);

                // Apply left complement to result: lcmpl(result)
                basis_tab = apply_rules_to_tab(transformed_tab, lcmpl_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Complement rules not found for regressive inner (need both rcmpl and lcmpl)");
            }
            break;
        }

        case ProductType::regressive_geometric: {
            // rgpr(A,B) = lcmpl(gpr(rcmpl(A), rcmpl(B)))
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");

            if (rcmpl_rules_it != config.complement_rules.end() &&
                lcmpl_rules_it != config.complement_rules.end()) {
                // Transform basis elements: rcmpl(A) and rcmpl(B)
                auto rcmpl_left = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);
                auto rcmpl_right = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // Compute geometric product of complements
                auto transformed_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(rcmpl_left, rcmpl_right, mul_str), rules);

                // Apply left complement to result: lcmpl(result)
                basis_tab = apply_rules_to_tab(transformed_tab, lcmpl_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Complement rules not found for regressive geometric (need both rcmpl and lcmpl)");
            }
            break;
        }

        case ProductType::right_bulk_contraction: {
            // A >> B = rwdg(A, bulk_dual(B)) = lcmpl(wdg(rcmpl(A), rcmpl(bulk_dual(B))))
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");
            auto bulk_dual_rules_it = config.complement_rules.find("bulk_dual");
            auto wdg_rules_it = config.product_rules.find("wdg");

            if (rcmpl_rules_it != config.complement_rules.end() &&
                lcmpl_rules_it != config.complement_rules.end() &&
                bulk_dual_rules_it != config.complement_rules.end() &&
                wdg_rules_it != config.product_rules.end()) {

                // Transform left operand: rcmpl(A)
                auto lhs_rcmpl = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // Transform right operand: rcmpl(bulk_dual(B))
                auto rhs_bulk_dual = apply_rules_to_mv(config.basis, bulk_dual_rules_it->second);
                auto rhs_bulk_dual_rcmpl = apply_rules_to_mv(rhs_bulk_dual, rcmpl_rules_it->second);

                // Compute wedge product of transformed operands
                auto transformed_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_bulk_dual_rcmpl, wdg_str),
                    wdg_rules_it->second);

                // Apply left complement to result: lcmpl(result)
                basis_tab = apply_rules_to_tab(transformed_tab, lcmpl_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Required rules not found for right bulk contraction");
            }
            break;
        }

        case ProductType::right_weight_contraction: {
            // A >> B = rwdg(A, weight_dual(B)) = lcmpl(wdg(rcmpl(A), rcmpl(weight_dual(B))))
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");
            auto weight_dual_rules_it = config.complement_rules.find("weight_dual");
            auto wdg_rules_it = config.product_rules.find("wdg");

            if (rcmpl_rules_it != config.complement_rules.end() &&
                lcmpl_rules_it != config.complement_rules.end() &&
                weight_dual_rules_it != config.complement_rules.end() &&
                wdg_rules_it != config.product_rules.end()) {

                // Transform left operand: rcmpl(A)
                auto lhs_rcmpl = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // Transform right operand: rcmpl(weight_dual(B))
                auto rhs_weight_dual = apply_rules_to_mv(config.basis, weight_dual_rules_it->second);
                auto rhs_weight_dual_rcmpl = apply_rules_to_mv(rhs_weight_dual, rcmpl_rules_it->second);

                // Compute wedge product of transformed operands
                auto transformed_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_weight_dual_rcmpl, wdg_str),
                    wdg_rules_it->second);

                // Apply left complement to result: lcmpl(result)
                basis_tab = apply_rules_to_tab(transformed_tab, lcmpl_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Required rules not found for right weight contraction");
            }
            break;
        }

        case ProductType::left_bulk_contraction: {
            // A << B = rwdg(left_bulk_dual(A), B) = lcmpl(wdg(rcmpl(left_bulk_dual(A)), rcmpl(B)))
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");
            auto left_bulk_dual_rules_it = config.complement_rules.find("left_bulk_dual");
            auto bulk_dual_rules_it = config.complement_rules.find("bulk_dual"); // fallback
            auto wdg_rules_it = config.product_rules.find("wdg");

            // Use left_bulk_dual if available, otherwise fallback to bulk_dual
            auto dual_rules_it = (left_bulk_dual_rules_it != config.complement_rules.end()) 
                ? left_bulk_dual_rules_it : bulk_dual_rules_it;

            if (rcmpl_rules_it != config.complement_rules.end() &&
                lcmpl_rules_it != config.complement_rules.end() &&
                dual_rules_it != config.complement_rules.end() &&
                wdg_rules_it != config.product_rules.end()) {

                // Transform left operand: rcmpl(left_bulk_dual(A))
                auto lhs_bulk_dual = apply_rules_to_mv(config.basis, dual_rules_it->second);
                auto lhs_bulk_dual_rcmpl = apply_rules_to_mv(lhs_bulk_dual, rcmpl_rules_it->second);

                // Transform right operand: rcmpl(B)
                auto rhs_rcmpl = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // Compute wedge product of transformed operands
                auto transformed_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_bulk_dual_rcmpl, rhs_rcmpl, wdg_str),
                    wdg_rules_it->second);

                // Apply left complement to result: lcmpl(result)
                basis_tab = apply_rules_to_tab(transformed_tab, lcmpl_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Required rules not found for left bulk contraction");
            }
            break;
        }

        case ProductType::left_weight_contraction: {
            // A << B = rwdg(left_weight_dual(A), B) = lcmpl(wdg(rcmpl(left_weight_dual(A)), rcmpl(B)))
            auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
            auto lcmpl_rules_it = config.complement_rules.find("lcmpl");
            auto left_weight_dual_rules_it = config.complement_rules.find("left_weight_dual");
            auto weight_dual_rules_it = config.complement_rules.find("weight_dual"); // fallback
            auto wdg_rules_it = config.product_rules.find("wdg");

            // Use left_weight_dual if available, otherwise fallback to weight_dual
            auto dual_rules_it = (left_weight_dual_rules_it != config.complement_rules.end()) 
                ? left_weight_dual_rules_it : weight_dual_rules_it;

            if (rcmpl_rules_it != config.complement_rules.end() &&
                lcmpl_rules_it != config.complement_rules.end() &&
                dual_rules_it != config.complement_rules.end() &&
                wdg_rules_it != config.product_rules.end()) {

                // Transform left operand: rcmpl(left_weight_dual(A))
                auto lhs_weight_dual = apply_rules_to_mv(config.basis, dual_rules_it->second);
                auto lhs_weight_dual_rcmpl = apply_rules_to_mv(lhs_weight_dual, rcmpl_rules_it->second);

                // Transform right operand: rcmpl(B)
                auto rhs_rcmpl = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                // Compute wedge product of transformed operands
                auto transformed_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_weight_dual_rcmpl, rhs_rcmpl, wdg_str),
                    wdg_rules_it->second);

                // Apply left complement to result: lcmpl(result)
                basis_tab = apply_rules_to_tab(transformed_tab, lcmpl_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Required rules not found for left weight contraction");
            }
            break;
        }

        case ProductType::right_bulk_expansion: {
            // A >< B = wdg(A, bulk_dual(B))
            auto bulk_dual_rules_it = config.complement_rules.find("bulk_dual");
            auto wdg_rules_it = config.product_rules.find("wdg");

            if (bulk_dual_rules_it != config.complement_rules.end() &&
                wdg_rules_it != config.product_rules.end()) {

                // Transform right operand: bulk_dual(B)
                auto rhs_bulk_dual = apply_rules_to_mv(config.basis, bulk_dual_rules_it->second);

                // Compute wedge product: wdg(A, bulk_dual(B))
                basis_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(config.basis, rhs_bulk_dual, wdg_str),
                    wdg_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Required rules not found for right bulk expansion");
            }
            break;
        }

        case ProductType::right_weight_expansion: {
            // A >< B = wdg(A, weight_dual(B))
            auto weight_dual_rules_it = config.complement_rules.find("weight_dual");
            auto wdg_rules_it = config.product_rules.find("wdg");

            if (weight_dual_rules_it != config.complement_rules.end() &&
                wdg_rules_it != config.product_rules.end()) {

                // Transform right operand: weight_dual(B)
                auto rhs_weight_dual = apply_rules_to_mv(config.basis, weight_dual_rules_it->second);

                // Compute wedge product: wdg(A, weight_dual(B))
                basis_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(config.basis, rhs_weight_dual, wdg_str),
                    wdg_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Required rules not found for right weight expansion");
            }
            break;
        }

        case ProductType::left_bulk_expansion: {
            // A >< B = wdg(left_bulk_dual(A), B)
            auto left_bulk_dual_rules_it = config.complement_rules.find("left_bulk_dual");
            auto bulk_dual_rules_it = config.complement_rules.find("bulk_dual"); // fallback
            auto wdg_rules_it = config.product_rules.find("wdg");

            // Use left_bulk_dual if available, otherwise fallback to bulk_dual
            auto dual_rules_it = (left_bulk_dual_rules_it != config.complement_rules.end()) 
                ? left_bulk_dual_rules_it : bulk_dual_rules_it;

            if (dual_rules_it != config.complement_rules.end() &&
                wdg_rules_it != config.product_rules.end()) {

                // Transform left operand: left_bulk_dual(A)
                auto lhs_bulk_dual = apply_rules_to_mv(config.basis, dual_rules_it->second);

                // Compute wedge product: wdg(bulk_dual(A), B)
                basis_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_bulk_dual, config.basis, wdg_str),
                    wdg_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Required rules not found for left bulk expansion");
            }
            break;
        }

        case ProductType::left_weight_expansion: {
            // A >< B = wdg(left_weight_dual(A), B)
            auto left_weight_dual_rules_it = config.complement_rules.find("left_weight_dual");
            auto weight_dual_rules_it = config.complement_rules.find("weight_dual"); // fallback
            auto wdg_rules_it = config.product_rules.find("wdg");

            // Use left_weight_dual if available, otherwise fallback to weight_dual
            auto dual_rules_it = (left_weight_dual_rules_it != config.complement_rules.end()) 
                ? left_weight_dual_rules_it : weight_dual_rules_it;

            if (dual_rules_it != config.complement_rules.end() &&
                wdg_rules_it != config.product_rules.end()) {

                // Transform left operand: left_weight_dual(A)
                auto lhs_weight_dual = apply_rules_to_mv(config.basis, dual_rules_it->second);

                // Compute wedge product: wdg(weight_dual(A), B)
                basis_tab = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_weight_dual, config.basis, wdg_str),
                    wdg_rules_it->second);
            }
            else {
                throw std::invalid_argument(
                    "Required rules not found for left weight expansion");
            }
            break;
        }

        // Sandwich products
        case ProductType::sandwich: {
            // Regular sandwich products: rotor * object * rev(rotor) - use geometric product
            basis_tab =
                apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(config.basis, config.basis,
                                                             product.operator_str),
                                   rules);
            break;
        }
        
        case ProductType::regressive_sandwich: {
            // Regressive sandwich: motor * object * rrev(motor) - use regressive geometric product
            // PGA3DP uses rcmpl(A) -> gpr -> lcmpl sequence, others use cmpl -> gpr -> cmpl
            if (config.name == "pga3dp") {
                // PGA3DP specific sequence: rcmpl(gpr(rcmpl(A), rcmpl(B)))
                auto rcmpl_rules_it = config.complement_rules.find("rcmpl");
                auto lcmpl_rules_it = config.complement_rules.find("lcmpl");
                
                if (rcmpl_rules_it != config.complement_rules.end() &&
                    lcmpl_rules_it != config.complement_rules.end()) {
                    
                    // Forward transformation: rcmpl(A), rcmpl(B)
                    auto lhs_rcmpl = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);
                    auto rhs_rcmpl = apply_rules_to_mv(config.basis, rcmpl_rules_it->second);

                    // Apply geometric product to transformed operands
                    auto basis_tab_with_rules = apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_rcmpl, product.operator_str),
                        rules);

                    // Backward transformation: lcmpl(result)
                    basis_tab = apply_rules_to_tab(basis_tab_with_rules, lcmpl_rules_it->second);
                } else {
                    // Fallback to regular product
                    basis_tab =
                        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(config.basis, config.basis,
                                                                     product.operator_str),
                                           rules);
                }
            } else {
                // Other PGA algebras: cmpl(gpr(cmpl(A), cmpl(B)))
                auto cmpl_rules_it = config.complement_rules.find("cmpl");
                if (cmpl_rules_it != config.complement_rules.end()) {
                    // Transform both operands: cmpl(A), cmpl(B)
                    auto lhs_cmpl = apply_rules_to_mv(config.basis, cmpl_rules_it->second);
                    auto rhs_cmpl = apply_rules_to_mv(config.basis, cmpl_rules_it->second);

                    // Apply geometric product to transformed operands
                    auto basis_tab_with_rules = apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_cmpl, rhs_cmpl, product.operator_str),
                        rules);

                    // Transform result back: cmpl(result)
                    basis_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_rules_it->second);
                } else {
                    // Fallback to regular product if complement rules not found
                    basis_tab =
                        apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(config.basis, config.basis,
                                                                     product.operator_str),
                                           rules);
                }
            }
            break;
        }

        // Other PGA products - placeholder implementations
        default: {
            // For now, use basic product table - full implementation would use
            // appropriate transformations
            basis_tab =
                apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(config.basis, config.basis,
                                                             product.operator_str),
                                   rules);
            break;
        }
    }

    // Note: Asymmetric transformation is handled directly in ProductType::commutator_product

    return basis_tab;
}

void ProductExpressionGenerator::generate_sandwich_products(const AlgebraConfig& config,
                                                          const ProductDefinition& product,
                                                          const prd_table& basis_tab)
{
    std::string prd_name = config.name + " " + product.display_name;
    
    // EGA2D sandwich: mv_e * vec -> vec_tmp, then vec_tmp * rev(mv_e) -> vec_res
    if (config.name == "ega2d") {
        // Step 1: mv_e * vec -> vec_tmp
        fmt::println("{}:", prd_name + space_str + "mv_e * vec -> vec_tmp");
        
        auto coeff_R_even_it = config.coefficients.find("R_even");
        auto coeff_svps_it = config.coefficients.find("svps");
        
        if (coeff_R_even_it != config.coefficients.end() && 
            coeff_svps_it != config.coefficients.end()) {
            
            auto prd_tab = get_prd_tab(basis_tab, coeff_R_even_it->second, coeff_svps_it->second);
            auto vec_tmp = get_mv_from_prd_tab(prd_tab, config.basis, filter_2d::mv_e, filter_2d::vec,
                                               brace_switch::use_braces);
            
            fmt::println("vec_tmp:");
            print_mvec(vec_tmp, config.basis);
            fmt::println("");
            
            // Step 2: vec_tmp * rev(mv_e) -> vec_res
            fmt::println("{}:", prd_name + space_str + "vec_tmp * rev(mv_e) -> vec_res");
            
            auto coeff_R_rev_even_it = config.coefficients.find("R_rev_even");
            if (coeff_R_rev_even_it != config.coefficients.end()) {
                // Generate product table for vec_tmp (which is the result from step 1) * R_rev_even
                auto prd_tab2 = get_prd_tab(basis_tab, vec_tmp, coeff_R_rev_even_it->second);
                auto vec_res = get_mv_from_prd_tab(prd_tab2, config.basis, filter_2d::vec, filter_2d::mv_e);
                
                print_mvec(vec_res, config.basis);
                fmt::println("");
            }
        }
    }
    
    // EGA3D sandwich: Two products: mv_e * vec and mv_e * bivec
    else if (config.name == "ega3d") {
        auto coeff_R_even_it = config.coefficients.find("R_even");
        auto coeff_svBps_it = config.coefficients.find("svBps");
        auto coeff_R_rev_even_it = config.coefficients.find("R_rev_even");
        
        if (coeff_R_even_it != config.coefficients.end() && 
            coeff_svBps_it != config.coefficients.end() &&
            coeff_R_rev_even_it != config.coefficients.end()) {
            
            // First sandwich: mv_e * vec -> mv_u_tmp, then mv_u_tmp * rev(mv_e) -> mv_u_res
            fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u_tmp");
            
            auto prd_tab = get_prd_tab(basis_tab, coeff_R_even_it->second, coeff_svBps_it->second);
            auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, config.basis, filter_3d::mv_e, filter_3d::vec,
                                               brace_switch::use_braces);
            
            fmt::println("mv_u_tmp:");
            print_mvec(mv_u_tmp, config.basis);
            fmt::println("");
            
            // Step 2: mv_u_tmp * rev(mv_e) -> mv_u_res
            fmt::println("{}:", prd_name + space_str + "mv_u_tmp * rev(mv_e) -> mv_u_res");
            
            auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, coeff_R_rev_even_it->second);
            auto mv_u_res = get_mv_from_prd_tab(prd_tab_v, config.basis, filter_3d::mv_u, filter_3d::mv_e);
            
            print_mvec(mv_u_res, config.basis);
            fmt::println("");
            
            // Second sandwich: mv_e * bivec -> mv_e_tmp, then mv_e_tmp * rev(mv_e) -> mv_e_res
            fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e_tmp");
            
            auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, config.basis, filter_3d::mv_e, filter_3d::bivec,
                                               brace_switch::use_braces);
            fmt::println("mv_e_tmp:");
            print_mvec(mv_e_tmp, config.basis);
            fmt::println("");
            
            // Step 2: mv_e_tmp * rev(mv_e) -> mv_e_res
            fmt::println("{}:", prd_name + space_str + "mv_e_tmp * rev(mv_e) -> mv_e_res");
            
            auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, coeff_R_rev_even_it->second);
            auto mv_e_res = get_mv_from_prd_tab(prd_tab_B, config.basis, filter_3d::mv_e, filter_3d::mv_e);
            
            print_mvec(mv_e_res, config.basis);
            fmt::println("");
        }
    }
    
    // PGA2DP regressive sandwich: Two products: mv_u * vec and mv_u * bivec
    else if (config.name == "pga2dp") {
        auto coeff_R_odd_it = config.coefficients.find("R_odd");
        auto coeff_svBps_it = config.coefficients.find("svBps");
        auto coeff_R_rrev_odd_it = config.coefficients.find("R_rrev_odd");
        
        if (coeff_R_odd_it != config.coefficients.end() && 
            coeff_svBps_it != config.coefficients.end() &&
            coeff_R_rrev_odd_it != config.coefficients.end()) {
            
            // First sandwich: rgpr(mv_u, vec) -> mv_u_tmp, then mv_u_tmp * rrev(mv_u) -> mv_u_res
            fmt::println("{}:", prd_name + space_str + "rgpr(mv_u, vec) -> mv_u_tmp");
            
            auto prd_tab = get_prd_tab(basis_tab, coeff_R_odd_it->second, coeff_svBps_it->second);
            auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, config.basis, filter_2dp::mv_u, filter_2dp::vec,
                                               brace_switch::use_braces);
            
            fmt::println("mv_u_tmp:");
            print_mvec(mv_u_tmp, config.basis);
            fmt::println("");
            
            // Step 2: mv_u_tmp * rrev(mv_u) -> mv_u_res
            fmt::println("{}:", prd_name + space_str + "mv_u_tmp * rrev(mv_u) -> mv_u_res");
            
            auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, coeff_R_rrev_odd_it->second);
            auto mv_u_res = get_mv_from_prd_tab(prd_tab_v, config.basis, filter_2dp::mv_u, filter_2dp::mv_u);
            
            print_mvec(mv_u_res, config.basis);
            fmt::println("");
            
            // Second sandwich: mv_u * bivec -> mv_e_tmp, then mv_e_tmp * rrev(mv_u) -> mv_e_res
            fmt::println("{}:", prd_name + space_str + "mv_u * bivec -> mv_e_tmp");
            
            auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, config.basis, filter_2dp::mv_u, filter_2dp::bivec,
                                               brace_switch::use_braces);
            fmt::println("mv_e_tmp:");
            print_mvec(mv_e_tmp, config.basis);
            fmt::println("");
            
            // Step 2: mv_e_tmp * rrev(mv_u) -> mv_e_res
            fmt::println("{}:", prd_name + space_str + "mv_e_tmp * rrev(mv_u) -> mv_e_res");
            
            auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, coeff_R_rrev_odd_it->second);
            auto mv_e_res = get_mv_from_prd_tab(prd_tab_B, config.basis, filter_2dp::mv_e, filter_2dp::mv_u);
            
            print_mvec(mv_e_res, config.basis);
            fmt::println("");
        }
    }
    
    // PGA3DP regressive sandwich: Three products: mv_e * vec, mv_e * bivec, mv_e * trivec
    else if (config.name == "pga3dp") {
        auto coeff_R_even_it = config.coefficients.find("R_even");
        auto coeff_svBtps_it = config.coefficients.find("svBtps");
        auto coeff_R_rrev_even_it = config.coefficients.find("R_rrev_even");
        
        if (coeff_R_even_it != config.coefficients.end() && 
            coeff_svBtps_it != config.coefficients.end() &&
            coeff_R_rrev_even_it != config.coefficients.end()) {
            
            // Create shared product table for all three sandwich products
            auto prd_tab = get_prd_tab(basis_tab, coeff_R_even_it->second, coeff_svBtps_it->second);
            
            // First sandwich: mv_e * vec -> mv_u_tmp, then mv_u_tmp * rrev(mv_e) -> mv_u_res
            fmt::println("{}:", prd_name + space_str + "mv_e * vec -> mv_u_tmp");
            
            auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, config.basis, filter_3dp::mv_e, filter_3dp::vec,
                                               brace_switch::use_braces);
            
            fmt::println("mv_u_tmp:");
            print_mvec(mv_u_tmp, config.basis);
            fmt::println("");
            
            // Step 2: mv_u_tmp * rrev(mv_e) -> mv_u_res
            fmt::println("{}:", prd_name + space_str + "mv_u_tmp * rrev(mv_e) -> mv_u_res");
            
            auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, coeff_R_rrev_even_it->second);
            auto mv_u_res = get_mv_from_prd_tab(prd_tab_v, config.basis, filter_3dp::mv_u, filter_3dp::mv_e);
            
            print_mvec(mv_u_res, config.basis);
            fmt::println("");
            
            // Second sandwich: mv_e * bivec -> mv_e_tmp, then mv_e_tmp * rrev(mv_e) -> mv_e_res
            fmt::println("{}:", prd_name + space_str + "mv_e * bivec -> mv_e_tmp");
            
            auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, config.basis, filter_3dp::mv_e, filter_3dp::bivec,
                                               brace_switch::use_braces);
            fmt::println("mv_e_tmp:");
            print_mvec(mv_e_tmp, config.basis);
            fmt::println("");
            
            // Step 2: mv_e_tmp * rrev(mv_e) -> mv_e_res
            fmt::println("{}:", prd_name + space_str + "mv_e_tmp * rrev(mv_e) -> mv_e_res");
            
            auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, coeff_R_rrev_even_it->second);
            auto mv_e_res = get_mv_from_prd_tab(prd_tab_B, config.basis, filter_3dp::mv_e, filter_3dp::mv_e);
            
            print_mvec(mv_e_res, config.basis);
            fmt::println("");
            
            // Third sandwich: mv_e * trivec -> mv_u_tmp_t, then mv_u_tmp_t * rrev(mv_e) -> mv_u_res_t
            fmt::println("{}:", prd_name + space_str + "mv_e * trivec -> mv_u_tmp_t");
            
            auto mv_u_tmp_t = get_mv_from_prd_tab(prd_tab, config.basis, filter_3dp::mv_e, filter_3dp::trivec,
                                                 brace_switch::use_braces);
            fmt::println("mv_u_tmp_t:");
            print_mvec(mv_u_tmp_t, config.basis);
            fmt::println("");
            
            // Step 2: mv_u_tmp_t * rrev(mv_e) -> mv_u_res_t
            fmt::println("{}:", prd_name + space_str + "mv_u_tmp_t * rrev(mv_e) -> mv_u_res_t");
            
            auto prd_tab_t = get_prd_tab(basis_tab, mv_u_tmp_t, coeff_R_rrev_even_it->second);
            auto mv_u_res_t = get_mv_from_prd_tab(prd_tab_t, config.basis, filter_3dp::mv_u, filter_3dp::mv_e);
            
            print_mvec(mv_u_res_t, config.basis);
            fmt::println("");
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Convenience functions for generating specific algebras
////////////////////////////////////////////////////////////////////////////////

void generate_ega2d_expressions()
{
    ProductExpressionGenerator generator;
    generator.generate_algebra(configs::ega2d_config, configs::ega2d_products);
}

void generate_ega3d_expressions()
{
    ProductExpressionGenerator generator;
    generator.generate_algebra(configs::ega3d_config, configs::ega3d_products);
}

void generate_pga2dp_expressions()
{
    ProductExpressionGenerator generator;
    generator.generate_algebra(configs::pga2dp_config, configs::pga2dp_products);
}

void generate_pga3dp_expressions()
{
    ProductExpressionGenerator generator;
    generator.generate_algebra(configs::pga3dp_config, configs::pga3dp_products);
}