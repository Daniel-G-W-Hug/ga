// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_trafo_nary_expression.hpp"
#include "ga_prdxpr_trafo_expr_simplifier.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>

#include "fmt/format.h"

///////////////////////////////////////////////////////////////////////////////
// NAryTerm implementation
///////////////////////////////////////////////////////////////////////////////

NAryTerm NAryTerm::operator*(const NAryTerm& other) const
{
    NAryTerm result;
    result.coefficient = this->coefficient * other.coefficient;

    // Combine variable maps
    result.variables = this->variables;
    for (const auto& [var, power] : other.variables) {
        result.variables[var] += power;
    }

    return result;
}

bool NAryTerm::canCombineWith(const NAryTerm& other) const
{
    // Terms can combine if they have the same variables with same powers
    return this->variables == other.variables;
}

NAryTerm NAryTerm::operator+(const NAryTerm& other) const
{
    if (!canCombineWith(other)) {
        throw std::runtime_error("Cannot add terms with different variable factors");
    }

    NAryTerm result;
    result.variables = this->variables;
    result.coefficient = this->coefficient + other.coefficient;

    return result;
}

bool NAryTerm::isZero() const { return std::abs(coefficient) < 1e-10; }

void NAryTerm::applyCommutativity()
{
    // Apply canonical ordering to variables using centralized function
    // Note: This recreates the map, which loses ordering, but ensures consistent internal
    // state
    // Use default patterns for backward compatibility
    variables =
        GAAlgebraRules::reorderCommutativeFactors(variables, GeometricVariablePatterns{});
}

void NAryTerm::applyCommutativity(const GeometricVariablePatterns& patterns)
{
    // Apply canonical ordering to variables using custom patterns
    variables = GAAlgebraRules::reorderCommutativeFactors(variables, patterns);
}

std::string NAryTerm::toString() const
{
    // Use default patterns for backward compatibility
    return toString(GeometricVariablePatterns{});
}

std::string NAryTerm::toString(const GeometricVariablePatterns& patterns) const
{
    if (isZero()) return "0";

    std::string result;

    // Handle coefficient
    if (std::abs(coefficient - 1.0) > 1e-10) {
        if (coefficient == -1.0 && !variables.empty()) {
            result = "-";
        }
        else {
            // Use floating point format for consistency with manual reference
            if (coefficient == floor(coefficient)) {
                result = fmt::format("{:.1f}", coefficient);
            }
            else {
                result = fmt::format("{}", coefficient);
            }
        }
    }
    else if (coefficient == -1.0 && variables.empty()) {
        result = "-1.0";
    }
    else if (variables.empty()) {
        result = "1.0";
    }

    // Add variables using pattern-based canonical ordering
    auto var_pairs = GAAlgebraRules::getSortedVariablePairs(variables, patterns);

    bool first_var = (result.empty() || result == "-");
    for (const auto& [var, power] : var_pairs) {
        if (!first_var) {
            result += " * ";
        }

        // Handle special FACTORED variables by removing the prefix
        if (var.starts_with("FACTORED:")) {
            result += var.substr(9); // Remove "FACTORED:" prefix
        }
        else {
            result += var;
        }

        if (power > 1) {
            result += fmt::format("^{}", power);
        }
        first_var = false;
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////
// NAryExpression implementation
///////////////////////////////////////////////////////////////////////////////

void NAryExpression::expandDistributiveProducts()
{
    // Step 1: Resolve braces by multiplication
    // In our current representation, terms are already expanded from the binary AST
    // This step is mainly handled during the conversion from binary AST

    // However, we need to handle cases where we have parenthesized sub-expressions
    // For now, we assume the conversion process handles the distributive expansion
}

void NAryExpression::normalizeSignsAndCommutativity()
{
    // Step 2: Sign normalization and apply commutativity (default patterns)
    for (auto& term : terms) {
        term.applyCommutativity();
    }
}

void NAryExpression::normalizeSignsAndCommutativity(
    const GeometricVariablePatterns& patterns)
{
    // Step 2: Sign normalization and apply commutativity with custom patterns
    for (auto& term : terms) {
        term.applyCommutativity(patterns);
    }
}

void NAryExpression::combineTermsAndRegroup()
{
    // Step 3: Combine like terms AND factor by common geometric variables
    std::vector<NAryTerm> result;


    // First pass: combine exact matches
    for (const auto& term : terms) {
        if (term.isZero()) continue;

        bool combined = false;

        // Try to combine with existing terms (exact variable match)
        for (auto& existing : result) {
            if (existing.canCombineWith(term)) {
                existing = existing + term;
                combined = true;
                break;
            }
        }

        if (!combined) {
            result.push_back(term);
        }
    }

    // Second pass: factor by common geometric variables (v.x, v.y, etc.)
    std::map<std::string, std::vector<size_t>> groups_by_geom_var;

    // Group terms by their geometric variable
    for (size_t i = 0; i < result.size(); ++i) {
        std::string geom_var = "";

        // Find the geometric variable (v.x, v.y, B.x, B.vx, B.mx, t.x, etc.)
        for (const auto& [var, power] : result[i].variables) {
            if (var.starts_with("v.") || var.starts_with("B.") || var.starts_with("t.")) {
                geom_var = var;
                break;
            }
        }

        if (!geom_var.empty()) {
            groups_by_geom_var[geom_var].push_back(i);
        }
    }


    // Create factored terms for groups with multiple terms
    std::vector<NAryTerm> final_result;
    std::set<size_t> processed;

    // Always factor by geometric variable - this is the correct behavior
    for (const auto& [geom_var, indices] : groups_by_geom_var) {
        if (indices.size() > 1) {
            // Multiple terms with this geometric variable - create factored expression

            // Create a factored term: (sum of rotor expressions) * geom_var
            NAryTerm factored_term;
            factored_term.coefficient = 1.0;
            factored_term.variables[geom_var] = 1;

            // Sort indices to match manual ordering: by rotor coefficient order first
            // (R.c0, R.c1, R.c2, R.c3, ...)
            auto sorted_indices = indices;
            std::sort(sorted_indices.begin(), sorted_indices.end(),
                      [&result](size_t a, size_t b) {
                          // Extract coefficient pattern recognition for sorting
                          auto get_first_coeff = [](const NAryTerm& term) -> std::string {
                              for (const auto& [var, power] : term.variables) {
                                  // Support common coefficient patterns: R.c, M.c, etc.
                                  if (var.contains(".c") &&
                                      (var.starts_with("R.") || var.starts_with("M."))) {
                                      return var;
                                  }
                              }
                              return "";
                          };

                          // Sort by coefficient order (R.c0 < R.c1 < R.c2 < R.c3 or M.c0
                          // < M.c1 < M.c2 < M.c3)
                          return get_first_coeff(result[a]) < get_first_coeff(result[b]);
                      });

            // Check for common factors to extract (like 2 in "-2*A + 2*B" -> "2*(-A +
            // B)")
            double common_factor = 0.0;
            bool can_factor = false;

            if (sorted_indices.size() > 1) {
                // Find the greatest common factor of all coefficients
                double first_abs_coeff = std::abs(result[sorted_indices[0]].coefficient);
                common_factor = first_abs_coeff;

                for (size_t i = 1; i < sorted_indices.size(); ++i) {
                    double curr_abs_coeff =
                        std::abs(result[sorted_indices[i]].coefficient);
                    // Find GCD of coefficients
                    if (curr_abs_coeff > 0 && first_abs_coeff > 0) {
                        double a = std::max(common_factor, curr_abs_coeff);
                        double b = std::min(common_factor, curr_abs_coeff);
                        while (b > 1e-10) {
                            double temp = b;
                            b = std::fmod(a, b);
                            a = temp;
                        }
                        common_factor = a;
                    }
                }

                // Only factor out if the common factor is greater than 1
                can_factor = (common_factor > 1.0 + 1e-10);
            }

            // Build coefficient expression string
            std::string coeff_expr = "";
            bool first = true;

            for (size_t idx : sorted_indices) {
                const auto& term = result[idx];

                // Build rotor coefficient string
                std::string rotor_str = "";
                double coeff =
                    can_factor ? (term.coefficient / common_factor) : term.coefficient;

                // SKIP ZERO COEFFICIENTS - they should not appear in the final expression
                if (std::abs(coeff) < 1e-10) {
                    processed.insert(idx);
                    continue;
                }

                // Build the rotor part (everything except geom_var)
                auto rotor_vars = term.variables;
                rotor_vars.erase(geom_var);

                if (!rotor_vars.empty()) {
                    for (const auto& [var, power] : rotor_vars) {
                        // Expand powers using multiplication (R.c0^2 -> R.c0 * R.c0)
                        for (int p = 0; p < power; ++p) {
                            if (!rotor_str.empty()) rotor_str += " * ";
                            rotor_str += var;
                        }
                    }

                    if (std::abs(coeff - 1.0) > 1e-10) {
                        if (std::abs(coeff + 1.0) < 1e-10) {
                            rotor_str = "-" + rotor_str; // -1 * R.c0 -> -R.c0
                        }
                        else {
                            // Use floating point format for consistency
                            if (coeff == floor(coeff)) {
                                rotor_str = fmt::format("{:.1f} * {}", coeff, rotor_str);
                            }
                            else {
                                rotor_str = fmt::format("{} * {}", coeff, rotor_str);
                            }
                        }
                    }
                    else if (coeff < 0) {
                        rotor_str = "-" + rotor_str;
                    }
                }
                else {
                    // Use floating point format for consistency
                    if (coeff == floor(coeff)) {
                        rotor_str = fmt::format("{:.1f}", coeff);
                    }
                    else {
                        rotor_str = fmt::format("{}", coeff);
                    }
                }

                // Add to coefficient expression
                if (first) {
                    coeff_expr = rotor_str;
                    first = false;
                }
                else {
                    if (coeff >= 0 && !rotor_str.starts_with("-")) {
                        coeff_expr += " + " + rotor_str;
                    }
                    else {
                        // For negative coefficients, create proper subtraction
                        std::string positive_rotor =
                            rotor_str.starts_with("-") ? rotor_str.substr(1) : rotor_str;
                        if (positive_rotor.starts_with("1 * ")) {
                            positive_rotor = positive_rotor.substr(4); // Remove "1 * "
                        }
                        coeff_expr += " - " + positive_rotor;
                    }
                }

                processed.insert(idx);
            }

            // If all coefficients were zero, skip this factored term entirely
            if (coeff_expr.empty()) {
                continue;
            }

            // Apply common factor if factored
            if (can_factor) {
                // Use floating point format for consistency
                if (common_factor == floor(common_factor)) {
                    coeff_expr = fmt::format("{:.1f} * ({})", common_factor, coeff_expr);
                }
                else {
                    coeff_expr = fmt::format("{} * ({})", common_factor, coeff_expr);
                }
            }
            else {
                // Even without common factor, we need parentheses for multiple terms
                if (sorted_indices.size() > 1) {
                    coeff_expr = fmt::format("({})", coeff_expr);
                }
            }

            // Store the factored expression as a variable name
            factored_term.variables.clear();
            factored_term.variables[coeff_expr] = 1;
            factored_term.variables[geom_var] = 1;

            final_result.push_back(factored_term);
        }
        else if (indices.size() == 1) {
            // Single term with this geometric variable
            final_result.push_back(result[indices[0]]);
            processed.insert(indices[0]);
        }
    }

    // Add terms without geometric variables
    for (size_t i = 0; i < result.size(); ++i) {
        if (processed.find(i) == processed.end()) {
            final_result.push_back(result[i]);
        }
    }

    // Sort terms by geometric variable to match manual output (user requirement)
    // Manual output is sorted according to the first variable in the term
    std::sort(final_result.begin(), final_result.end(),
              [](const NAryTerm& a, const NAryTerm& b) {
                  // Extract geometric variables for comparison
                  std::string geom_var_a = "";
                  std::string geom_var_b = "";

                  for (const auto& [var, power] : a.variables) {
                      if (var.starts_with("v.") || var.starts_with("B.") ||
                          var.starts_with("t.")) {
                          geom_var_a = var;
                          break;
                      }
                  }

                  for (const auto& [var, power] : b.variables) {
                      if (var.starts_with("v.") || var.starts_with("B.") ||
                          var.starts_with("t.")) {
                          geom_var_b = var;
                          break;
                      }
                  }

                  // Sort by geometric variable in dimensional order
                  // Vectors: v.x, v.y, v.z, v.w
                  // Bivectors: B.vx, B.vy, B.vz, B.mx, B.my, B.mz
                  // Trivectors: t.x, t.y, t.z, t.w
                  auto getDimensionalOrder = [](const std::string& var) -> int {
                      // Vector ordering
                      if (var == "v.x" || var == "B.x") return 1;
                      if (var == "v.y" || var == "B.y") return 2;
                      if (var == "v.z" || var == "B.z") return 3;
                      if (var == "v.w" || var == "B.w") return 4;

                      // Bivector ordering: B.vx, B.vy, B.vz, B.mx, B.my, B.mz
                      if (var == "B.vx") return 5;
                      if (var == "B.vy") return 6;
                      if (var == "B.vz") return 7;
                      if (var == "B.mx") return 8;
                      if (var == "B.my") return 9;
                      if (var == "B.mz") return 10;

                      // Trivector ordering: t.x, t.y, t.z, t.w
                      if (var == "t.x") return 11;
                      if (var == "t.y") return 12;
                      if (var == "t.z") return 13;
                      if (var == "t.w") return 14;

                      return 0; // fallback for other variables
                  };

                  int order_a = getDimensionalOrder(geom_var_a);
                  int order_b = getDimensionalOrder(geom_var_b);

                  if (order_a != 0 && order_b != 0) {
                      return order_a < order_b;
                  }

                  // Fallback to alphabetical for other cases
                  return geom_var_a < geom_var_b;
              });

    terms = final_result;
}

void NAryExpression::factorCommonVariables()
{
    // Step 4: Skip factorization for now - the issue is in AST conversion, not in N-ary
    // logic
    // TODO: Implement proper factorization later for term sorting and grouping
}

void NAryExpression::removeZeroTerms()
{
    terms.erase(std::remove_if(terms.begin(), terms.end(),
                               [](const NAryTerm& t) { return t.isZero(); }),
                terms.end());
}

std::string NAryExpression::toString() const
{
    if (terms.empty()) return "0";

    std::string result;
    for (size_t i = 0; i < terms.size(); ++i) {
        if (i > 0) {
            if (terms[i].coefficient >= 0) {
                result += " + ";
            }
            else {
                result += " - ";
                // Create a positive version of the term for display
                NAryTerm positive_term = terms[i];
                positive_term.coefficient = std::abs(positive_term.coefficient);
                result += positive_term.toString();
                continue;
            }
        }
        result += terms[i].toString();
    }

    return result;
}

std::string NAryExpression::toString(const GeometricVariablePatterns& patterns) const
{
    if (terms.empty()) return "0";

    std::string result;
    for (size_t i = 0; i < terms.size(); ++i) {
        if (i > 0) {
            if (terms[i].coefficient >= 0) {
                result += " + ";
            }
            else {
                result += " - ";
                // Create a positive version of the term for display
                NAryTerm positive_term = terms[i];
                positive_term.coefficient = std::abs(positive_term.coefficient);
                result += positive_term.toString(patterns); // Use patterns here
                continue;
            }
        }
        result += terms[i].toString(patterns); // Use patterns here
    }

    return result;
}

void NAryExpression::printDebug(const std::string& step_name) const
{
    fmt::println("=== {} ===", step_name);
    fmt::println("Terms ({}): {}", terms.size(), toString());
    for (size_t i = 0; i < terms.size(); ++i) {
        fmt::println("  Term {}: coeff={}, vars={}", i, terms[i].coefficient,
                     terms[i].variables.size());
        for (const auto& [var, power] : terms[i].variables) {
            fmt::println("    {} ^ {}", var, power);
        }
    }
    fmt::println("");
}

///////////////////////////////////////////////////////////////////////////////
// NAryConverter implementation
///////////////////////////////////////////////////////////////////////////////

NAryExpression NAryConverter::fromBinaryAST(std::shared_ptr<ast_node> ast)
{
    std::vector<NAryTerm> terms;

    if (auto expr = std::dynamic_pointer_cast<Expression>(ast)) {
        extractTermsFromExpression(expr, terms);
    }
    else if (auto term = std::dynamic_pointer_cast<Term>(ast)) {
        extractTermsFromTerm(term, terms);
    }
    else if (auto factor = std::dynamic_pointer_cast<Factor>(ast)) {
        NAryTerm single_term;
        extractFactorsFromTerm(factor, single_term);
        terms.push_back(single_term);
    }
    else if (auto primary = std::dynamic_pointer_cast<Primary>(ast)) {
        NAryTerm single_term;
        if (primary->type == Primary::Primary_t::NUMBER) {
            single_term.coefficient = primary->num_value;
        }
        else if (primary->type == Primary::Primary_t::VARIABLE) {
            single_term.variables[primary->str_value] = 1;
        }
        terms.push_back(single_term);
    }

    return NAryExpression(terms);
}

void NAryConverter::extractTermsFromExpression(std::shared_ptr<Expression> expr,
                                               std::vector<NAryTerm>& terms, double sign)
{
    // Handle left side
    if (expr->left) {
        if (auto left_expr = std::dynamic_pointer_cast<Expression>(expr->left)) {
            extractTermsFromExpression(left_expr, terms, sign);
        }
        else if (auto left_term = std::dynamic_pointer_cast<Term>(expr->left)) {
            extractTermsFromTerm(left_term, terms, sign);
        }
    }

    // Handle right side with correct sign
    if (expr->right) {
        double right_sign = (expr->op == '+') ? sign : -sign;
        if (auto right_expr = std::dynamic_pointer_cast<Expression>(expr->right)) {
            extractTermsFromExpression(right_expr, terms, right_sign);
        }
        else if (auto right_term = std::dynamic_pointer_cast<Term>(expr->right)) {
            extractTermsFromTerm(right_term, terms, right_sign);
        }
    }
}

void NAryConverter::extractTermsFromTerm(std::shared_ptr<Term> term,
                                         std::vector<NAryTerm>& terms, double sign)
{
    // For multiplication terms, we need to handle distributive property
    // This is complex - for now, let's handle simple cases

    NAryTerm current_term;
    current_term.coefficient = sign;

    // Extract factors from left and right
    if (term->left) {
        extractFactorsFromTerm(term->left, current_term);
    }
    if (term->right) {
        extractFactorsFromTerm(term->right, current_term);
    }

    terms.push_back(current_term);
}

void NAryConverter::extractFactorsFromTerm(std::shared_ptr<ast_node> node,
                                           NAryTerm& current_term, double sign)
{
    if (auto factor = std::dynamic_pointer_cast<Factor>(node)) {
        double factor_sign = (factor->sign == '-') ? -1.0 : 1.0;
        current_term.coefficient *= sign * factor_sign;

        if (auto primary = std::dynamic_pointer_cast<Primary>(factor->prim_val)) {
            if (primary->type == Primary::Primary_t::NUMBER) {
                current_term.coefficient *= primary->num_value;
            }
            else if (primary->type == Primary::Primary_t::VARIABLE) {
                current_term.variables[primary->str_value]++;
            }
            else if (primary->type == Primary::Primary_t::EXPRESSION) {
                // This is the critical case: we have a parenthesized expression as a
                // factor We need to flag this for distributive expansion, but we can't
                // handle it here because we're building a single term. This needs to be
                // handled at a higher level.
                fmt::println("WARNING: Nested expressions need distributive handling - "
                             "this requires rework");
            }
        }
    }
    else if (auto term = std::dynamic_pointer_cast<Term>(node)) {
        // Recursively handle term multiplication
        extractFactorsFromTerm(term->left, current_term, sign);
        extractFactorsFromTerm(term->right, current_term, sign);
    }
}

std::shared_ptr<ast_node> NAryConverter::toBinaryAST(const NAryExpression& expr)
{
    if (expr.terms.empty()) {
        auto primary = std::make_shared<Primary>();
        primary->type = Primary::Primary_t::NUMBER;
        primary->str_value = "0";
        primary->num_value = 0.0;
        return primary;
    }

    if (expr.terms.size() == 1) {
        return buildTermAST(expr.terms[0]);
    }

    // Build expression tree for multiple terms
    return buildExpressionAST(expr.terms);
}

std::shared_ptr<ast_node> NAryConverter::buildTermAST(const NAryTerm& term)
{
    std::shared_ptr<ast_node> result;

    // Handle coefficient
    if (std::abs(term.coefficient - 1.0) > 1e-10) {
        auto coeff_primary = std::make_shared<Primary>();
        coeff_primary->type = Primary::Primary_t::NUMBER;
        // Use floating point format for consistency with manual reference
        double abs_coeff = std::abs(term.coefficient);
        if (abs_coeff == floor(abs_coeff)) {
            coeff_primary->str_value = fmt::format("{:.1f}", abs_coeff);
        }
        else {
            coeff_primary->str_value = fmt::format("{}", abs_coeff);
        }
        coeff_primary->num_value = abs_coeff;

        auto coeff_factor = std::make_shared<Factor>('\0', coeff_primary);
        if (term.coefficient < 0) {
            coeff_factor = std::make_shared<Factor>('-', coeff_primary);
        }

        result = coeff_factor;
    }

    // Add variables in the correct order: FACTORED expressions first, then geometric
    // variables This ensures the format: coefficient * factored_expression *
    // geometric_variable

    // First pass: Process FACTORED expressions
    for (const auto& [var, power] : term.variables) {
        if (var.starts_with("FACTORED:")) {
            // Handle factored expressions specially
            std::string factored_expr = var.substr(9); // Remove "FACTORED:" prefix

            // Parse the factored expression
            Parser factored_parser(factored_expr);
            auto factored_ast = factored_parser.parse();

            // Check if this is a simple coefficient * (...) expression that doesn't need
            // extra parentheses
            bool needs_parentheses = true;
            if (factored_expr.find(" * (") != std::string::npos &&
                (factored_expr.starts_with("2.0 *") ||
                 factored_expr.starts_with("-2.0 *") ||
                 factored_expr.starts_with("1.0 *") ||
                 factored_expr.starts_with("-1.0 *"))) {
                needs_parentheses = false;
            }

            std::shared_ptr<ast_node> factored_node;
            if (needs_parentheses) {
                // Wrap in parentheses for complex expressions
                auto paren_primary = std::make_shared<Primary>();
                paren_primary->type = Primary::Primary_t::EXPRESSION;
                paren_primary->expr = factored_ast;
                factored_node = std::make_shared<Factor>('\0', paren_primary);
            }
            else {
                // Use the expression directly for simple coefficient * (...) cases
                factored_node = factored_ast;
            }

            if (!result) {
                result = factored_node;
            }
            else {
                result = std::make_shared<Term>(result, '*', factored_node);
            }
        }
    }

    // Second pass: Process geometric variables (v.x, B.y, etc.)
    for (const auto& [var, power] : term.variables) {
        if (!var.starts_with("FACTORED:")) {
            // Handle normal variables
            auto var_primary = std::make_shared<Primary>();
            var_primary->type = Primary::Primary_t::VARIABLE;
            var_primary->str_value = var;

            auto var_factor = std::make_shared<Factor>('\0', var_primary);

            // Handle powers > 1 by repeated multiplication
            for (int p = 0; p < power; ++p) {
                if (!result) {
                    result = var_factor;
                }
                else {
                    result = std::make_shared<Term>(result, '*', var_factor);
                }
            }
        }
    }

    // Handle case where we only have coefficient 1 with no variables
    if (!result) {
        auto one_primary = std::make_shared<Primary>();
        one_primary->type = Primary::Primary_t::NUMBER;
        // Use floating point format for consistency
        one_primary->str_value = (term.coefficient >= 0) ? "1.0" : "-1.0";
        one_primary->num_value = term.coefficient;

        char sign = (term.coefficient >= 0) ? '\0' : '-';
        return std::make_shared<Factor>(sign, one_primary);
    }

    return result;
}

std::shared_ptr<ast_node>
NAryConverter::buildExpressionAST(const std::vector<NAryTerm>& terms)
{
    std::shared_ptr<ast_node> result = buildTermAST(terms[0]);

    for (size_t i = 1; i < terms.size(); ++i) {
        char op = (terms[i].coefficient >= 0) ? '+' : '-';

        // Create positive version of term for negative coefficients
        NAryTerm term_for_ast = terms[i];
        if (terms[i].coefficient < 0) {
            term_for_ast.coefficient = std::abs(terms[i].coefficient);
        }

        auto right_term = buildTermAST(term_for_ast);
        result = std::make_shared<Expression>(result, op, right_term);
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////
// NAryExpressionSimplifier implementation
///////////////////////////////////////////////////////////////////////////////

std::shared_ptr<ast_node>
NAryExpressionSimplifier::simplify(std::shared_ptr<ast_node> ast)
{
    // HYBRID APPROACH: Use existing distributive expansion first, then convert to N-ary

    // Step 1: Use the working distributive expansion from the original system
    auto expanded_ast = ExpressionSimplifier::expandProducts(ast);

    // Step 2: Convert the already-expanded AST to N-ary (this should be simpler now)
    auto nary_expr = NAryConverter::fromBinaryAST(expanded_ast);

    // Step 3: Apply remaining transformations in N-ary space (no conversions!)
    nary_expr.normalizeSignsAndCommutativity();
    nary_expr.combineTermsAndRegroup();
    nary_expr.removeZeroTerms();

    // Step 4: Convert back to binary AST (once)
    auto result = NAryConverter::toBinaryAST(nary_expr);
    return result;
}

NAryExpression NAryExpressionSimplifier::simplifyToNAry(std::shared_ptr<ast_node> ast)
{
    auto nary_expr = NAryConverter::fromBinaryAST(ast);
    nary_expr.simplify();
    return nary_expr;
}