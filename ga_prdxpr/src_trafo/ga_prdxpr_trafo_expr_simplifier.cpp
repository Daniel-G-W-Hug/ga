// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_trafo_expr_simplifier.hpp"

#include <algorithm>
#include <cmath>
#include <regex>
#include <sstream>

#include "fmt/format.h"

///////////////////////////////////////////////////////////////////////////////
// SimplifiedTerm implementation
///////////////////////////////////////////////////////////////////////////////

SimplifiedTerm::SimplifiedTerm() : numeric_coefficient(1.0), sign(1) {}

SimplifiedTerm::SimplifiedTerm(double coeff, const std::map<std::string, int>& vars,
                               int s) :
    numeric_coefficient(std::abs(coeff)), factors(vars), sign(s * (coeff >= 0 ? 1 : -1))
{
    normalize();
}

SimplifiedTerm SimplifiedTerm::fromString(const std::string& term_str)
{
    SimplifiedTerm term;

    // Parse string like "2 * R.c0 * R.c1 * v.x" or "-R.c2 * v.z"
    std::string clean_str = term_str;

    // Handle sign
    if (clean_str.starts_with("-")) {
        term.sign = -1;
        clean_str = clean_str.substr(1);
    }
    else if (clean_str.starts_with("+")) {
        clean_str = clean_str.substr(1);
    }

    // Split by * and parse factors
    std::regex factor_regex(R"(\s*\*?\s*([^\*\s]+)\s*)");
    std::sregex_iterator iter(clean_str.begin(), clean_str.end(), factor_regex);
    std::sregex_iterator end;

    for (; iter != end; ++iter) {
        std::string factor = (*iter)[1].str();

        // Try to parse as number
        try {
            double num = std::stod(factor);
            term.numeric_coefficient *= num;
        }
        catch (...) {
            // It's a variable - add to factors map
            if (!factor.empty()) {
                term.factors[factor]++;
            }
        }
    }

    term.normalize();
    return term;
}

SimplifiedTerm SimplifiedTerm::operator*(const SimplifiedTerm& other) const
{
    SimplifiedTerm result;

    result.numeric_coefficient = this->numeric_coefficient * other.numeric_coefficient;
    result.sign = this->sign * other.sign;

    // Combine variable factors
    result.factors = this->factors;
    for (const auto& [var, power] : other.factors) {
        result.factors[var] += power;
    }

    result.normalize();
    return result;
}

SimplifiedTerm SimplifiedTerm::operator+(const SimplifiedTerm& other) const
{
    if (!canCombineWith(other)) {
        throw std::runtime_error("Cannot add terms with different variable factors");
    }

    SimplifiedTerm result;
    // Use commutatively reordered factors for consistency
    result.factors = GAAlgebraRules::reorderCommutativeFactors(this->factors);

    double this_coeff = this->sign * this->numeric_coefficient;
    double other_coeff = other.sign * other.numeric_coefficient;
    double combined_coeff = this_coeff + other_coeff;

    result.numeric_coefficient = std::abs(combined_coeff);
    result.sign = (combined_coeff >= 0) ? 1 : -1;

    result.normalize();
    return result;
}

bool SimplifiedTerm::canCombineWith(const SimplifiedTerm& other) const
{
    // Check if terms can be combined after applying commutativity
    auto this_reordered = GAAlgebraRules::reorderCommutativeFactors(this->factors);
    auto other_reordered = GAAlgebraRules::reorderCommutativeFactors(other.factors);
    return this_reordered == other_reordered;
}

bool SimplifiedTerm::isZero() const { return std::abs(numeric_coefficient) < 1e-10; }

bool SimplifiedTerm::isEquivalentTo(const SimplifiedTerm& other) const
{
    // Check if terms are equivalent after applying commutativity
    auto this_reordered = GAAlgebraRules::reorderCommutativeFactors(this->factors);
    auto other_reordered = GAAlgebraRules::reorderCommutativeFactors(other.factors);

    return this_reordered == other_reordered &&
           std::abs(this->numeric_coefficient - other.numeric_coefficient) < 1e-10 &&
           this->sign == other.sign;
}

std::string SimplifiedTerm::toString() const
{
    if (isZero()) return "0";

    std::string result;

    // Add sign
    if (sign < 0) result += "-";

    // Add numeric coefficient if not 1
    if (std::abs(numeric_coefficient - 1.0) > 1e-10) {
        result += fmt::format("{}", numeric_coefficient);
    }

    // Add variable factors using centralized canonical ordering
    auto factor_pairs = GAAlgebraRules::getSortedVariablePairs(factors);

    bool has_coefficient = (std::abs(numeric_coefficient - 1.0) > 1e-10);
    bool first_factor = true;
    for (const auto& [var, power] : factor_pairs) {
        if (!first_factor || has_coefficient) {
            result += " * ";
        }

        result += var;
        if (power > 1) {
            result += fmt::format("^{}", power);
        }
        first_factor = false;
    }

    // Handle case where we only have a numeric coefficient
    if (factors.empty() && std::abs(numeric_coefficient - 1.0) <= 1e-10) {
        result = (sign > 0) ? "1" : "-1";
    }

    return result;
}

void SimplifiedTerm::normalize()
{
    // Remove zero-power factors
    for (auto it = factors.begin(); it != factors.end();) {
        if (it->second == 0) {
            it = factors.erase(it);
        }
        else {
            ++it;
        }
    }

    // Handle zero coefficient
    if (isZero()) {
        factors.clear();
        sign = 1;
        numeric_coefficient = 0.0;
    }
}

///////////////////////////////////////////////////////////////////////////////
// ExpressionSimplifier implementation
///////////////////////////////////////////////////////////////////////////////

std::shared_ptr<ast_node> ExpressionSimplifier::simplify(std::shared_ptr<ast_node> ast)
{
    // Multi-pass simplification pipeline
    auto step1 = expandProducts(ast);
    auto step2 = applyCommutativity(step1);
    auto step3 = combineTerms(step2);
    auto step4 = cancelZeros(step3);

    return step4;
}

std::shared_ptr<ast_node>
ExpressionSimplifier::expandProducts(std::shared_ptr<ast_node> ast)
{
    // WORKAROUND: Use string-based expansion to avoid AST traversal bugs
    // This ensures we capture all terms correctly

    std::string ast_str = ast->to_string();

    // Apply repeated re-parsing until we get a stable result
    // This handles cases where nested expressions don't expand correctly on first pass
    std::string prev_str = "";
    std::string curr_str = ast_str;
    int iterations = 0;
    const int max_iterations = 5;

    while (prev_str != curr_str && iterations < max_iterations) {
        prev_str = curr_str;

        // Parse and extract terms
        Parser parser(curr_str);
        auto parsed_ast = parser.parse();
        auto terms = astToTerms(parsed_ast);

        // Convert back to string
        auto new_ast = termsToAst(terms);
        curr_str = new_ast->to_string();

        iterations++;
    }

    // Final parse of the stabilized string
    Parser final_parser(curr_str);
    return final_parser.parse();
}

std::shared_ptr<ast_node>
ExpressionSimplifier::combineTerms(std::shared_ptr<ast_node> ast)
{
    auto terms = astToTerms(ast);
    auto combined = combineTermsVector(terms);
    return termsToAst(combined);
}

std::shared_ptr<ast_node>
ExpressionSimplifier::applyCommutativity(std::shared_ptr<ast_node> ast)
{
    auto terms = astToTerms(ast);

    // Apply commutative reordering to each term
    for (auto& term : terms) {
        term = applyCommutativeReordering(term);
    }

    return termsToAst(terms);
}


std::shared_ptr<ast_node> ExpressionSimplifier::cancelZeros(std::shared_ptr<ast_node> ast)
{
    auto terms = astToTerms(ast);

    // Remove zero terms
    terms.erase(std::remove_if(terms.begin(), terms.end(),
                               [](const SimplifiedTerm& t) { return t.isZero(); }),
                terms.end());

    return termsToAst(terms);
}

std::vector<SimplifiedTerm>
ExpressionSimplifier::astToTerms(std::shared_ptr<ast_node> ast)
{
    std::vector<SimplifiedTerm> terms;

    if (auto expr = std::dynamic_pointer_cast<Expression>(ast)) {
        extractTermsFromExpression(expr, terms);
    }
    else if (auto term = std::dynamic_pointer_cast<Term>(ast)) {
        extractTermsFromTerm(term, terms);
    }
    else if (auto factor = std::dynamic_pointer_cast<Factor>(ast)) {
        SimplifiedTerm simple_term;
        extractTermsFromFactor(factor, simple_term);
        terms.push_back(simple_term);
    }
    else if (auto primary = std::dynamic_pointer_cast<Primary>(ast)) {
        SimplifiedTerm simple_term;
        if (primary->type == Primary::Primary_t::NUMBER) {
            simple_term.numeric_coefficient = primary->num_value;
        }
        else if (primary->type == Primary::Primary_t::VARIABLE) {
            simple_term.factors[primary->str_value] = 1;
        }
        terms.push_back(simple_term);
    }

    return terms;
}

void ExpressionSimplifier::extractTermsFromExpression(std::shared_ptr<Expression> expr,
                                                      std::vector<SimplifiedTerm>& terms,
                                                      int sign)
{
    if (expr->left) {
        if (auto left_expr = std::dynamic_pointer_cast<Expression>(expr->left)) {
            extractTermsFromExpression(left_expr, terms, sign);
        }
        else if (auto left_term = std::dynamic_pointer_cast<Term>(expr->left)) {
            extractTermsFromTerm(left_term, terms, sign);
        }
    }

    if (expr->right) {
        int right_sign = (expr->op == '+') ? sign : -sign;
        if (auto right_expr = std::dynamic_pointer_cast<Expression>(expr->right)) {
            extractTermsFromExpression(right_expr, terms, right_sign);
        }
        else if (auto right_term = std::dynamic_pointer_cast<Term>(expr->right)) {
            extractTermsFromTerm(right_term, terms, right_sign);
        }
    }
}

void ExpressionSimplifier::extractTermsFromTerm(std::shared_ptr<Term> term,
                                                std::vector<SimplifiedTerm>& terms,
                                                int sign)
{
    // Check if this term contains parenthesized expressions that need distribution
    std::vector<std::vector<SimplifiedTerm>> factor_groups;

    // Process left side
    if (auto left_term = std::dynamic_pointer_cast<Term>(term->left)) {
        std::vector<SimplifiedTerm> left_terms;
        extractTermsFromTerm(left_term, left_terms, 1);
        factor_groups.push_back(left_terms);
    }
    else if (auto left_factor = std::dynamic_pointer_cast<Factor>(term->left)) {
        std::vector<SimplifiedTerm> left_terms = extractFactorTerms(left_factor, 1);
        factor_groups.push_back(left_terms);
    }

    // Process right side
    if (auto right_factor = std::dynamic_pointer_cast<Factor>(term->right)) {
        std::vector<SimplifiedTerm> right_terms = extractFactorTerms(right_factor, 1);
        factor_groups.push_back(right_terms);
    }

    // Apply distributive property: multiply all combinations
    if (factor_groups.empty()) {
        return;
    }

    std::vector<SimplifiedTerm> result_terms = factor_groups[0];

    for (size_t i = 1; i < factor_groups.size(); ++i) {
        std::vector<SimplifiedTerm> new_result;
        for (const auto& left_term : result_terms) {
            for (const auto& right_term : factor_groups[i]) {
                SimplifiedTerm combined = left_term * right_term;
                combined.sign *= sign;
                new_result.push_back(combined);
            }
        }
        result_terms = new_result;
    }

    // Add all resulting terms
    for (const auto& result_term : result_terms) {
        terms.push_back(result_term);
    }
}

void ExpressionSimplifier::extractTermsFromFactor(std::shared_ptr<Factor> factor,
                                                  SimplifiedTerm& current_term, int sign)
{
    int factor_sign = (factor->sign == '-') ? -1 : 1;
    current_term.sign *= sign * factor_sign;

    if (auto primary = std::dynamic_pointer_cast<Primary>(factor->prim_val)) {
        if (primary->type == Primary::Primary_t::NUMBER) {
            current_term.numeric_coefficient *= primary->num_value;
        }
        else if (primary->type == Primary::Primary_t::VARIABLE) {
            current_term.factors[primary->str_value]++;
        }
        else if (primary->type == Primary::Primary_t::EXPRESSION) {
            // Handle parenthesized expressions
            auto sub_terms = astToTerms(primary->expr);
            // For now, treat as single complex factor - more sophisticated handling
            // needed
            if (sub_terms.size() == 1) {
                current_term = current_term * sub_terms[0];
            }
        }
    }
}

std::vector<SimplifiedTerm>
ExpressionSimplifier::extractFactorTerms(std::shared_ptr<Factor> factor, int sign)
{
    int factor_sign = (factor->sign == '-') ? -1 : 1;
    int total_sign = sign * factor_sign;

    if (auto primary = std::dynamic_pointer_cast<Primary>(factor->prim_val)) {
        if (primary->type == Primary::Primary_t::NUMBER) {
            SimplifiedTerm term;
            term.numeric_coefficient = primary->num_value;
            term.sign = total_sign;
            return {term};
        }
        else if (primary->type == Primary::Primary_t::VARIABLE) {
            SimplifiedTerm term;
            term.factors[primary->str_value] = 1;
            term.sign = total_sign;
            return {term};
        }
        else if (primary->type == Primary::Primary_t::EXPRESSION) {
            // Handle parenthesized expressions - this is the key fix!
            auto sub_terms = astToTerms(primary->expr);
            // Apply the sign to all sub-terms
            for (auto& term : sub_terms) {
                term.sign *= total_sign;
            }
            return sub_terms;
        }
    }

    // Default case - return empty term
    SimplifiedTerm empty_term;
    empty_term.sign = total_sign;
    return {empty_term};
}

std::vector<SimplifiedTerm>
ExpressionSimplifier::combineTermsVector(const std::vector<SimplifiedTerm>& terms)
{
    std::vector<SimplifiedTerm> result;

    for (auto term : terms) {
        // Apply commutative reordering first
        term = applyCommutativeReordering(term);

        bool combined = false;

        // Try to combine with existing terms
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

    return result;
}

SimplifiedTerm
ExpressionSimplifier::applyCommutativeReordering(const SimplifiedTerm& term)
{
    SimplifiedTerm result = term;
    result.factors = GAAlgebraRules::reorderCommutativeFactors(term.factors);
    return result;
}

std::shared_ptr<ast_node>
ExpressionSimplifier::termsToAst(const std::vector<SimplifiedTerm>& terms)
{
    if (terms.empty()) {
        // Return zero
        auto primary = std::make_shared<Primary>();
        primary->type = Primary::Primary_t::NUMBER;
        primary->str_value = "0";
        primary->num_value = 0.0;
        return primary;
    }

    if (terms.size() == 1) {
        // Single term - convert directly
        return convertSingleTermToAst(terms[0]);
    }

    // Multiple terms - build expression tree
    std::shared_ptr<ast_node> result = convertSingleTermToAst(terms[0]);

    for (size_t i = 1; i < terms.size(); ++i) {
        char op = (terms[i].sign > 0) ? '+' : '-';
        auto right_term = convertSingleTermToAst(terms[i]);

        // Make right term positive since sign is handled by operator
        if (terms[i].sign < 0) {
            // Need to make the term positive since we're using '-' operator
            auto modified_term = terms[i];
            modified_term.sign = 1;
            right_term = convertSingleTermToAst(modified_term);
        }

        result = std::make_shared<Expression>(result, op, right_term);
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////
// Helper function implementation
///////////////////////////////////////////////////////////////////////////////

std::shared_ptr<ast_node> convertSingleTermToAst(const SimplifiedTerm& term)
{
    if (term.isZero()) {
        auto primary = std::make_shared<Primary>();
        primary->type = Primary::Primary_t::NUMBER;
        primary->str_value = "0";
        primary->num_value = 0.0;
        return primary;
    }

    // Build AST for: [sign] * [numeric_coeff] * [var1] * [var2] * ...
    std::shared_ptr<ast_node> result;

    // Start with numeric coefficient if not 1
    if (std::abs(term.numeric_coefficient - 1.0) > 1e-10) {
        auto coeff_primary = std::make_shared<Primary>();
        coeff_primary->type = Primary::Primary_t::NUMBER;
        coeff_primary->str_value = fmt::format("{}", term.numeric_coefficient);
        coeff_primary->num_value = term.numeric_coefficient;
        result = coeff_primary;
    }

    // Add variable factors
    for (const auto& [var, power] : term.factors) {
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

    // Handle case where we only have coefficient 1 with no variables
    if (!result) {
        auto one_primary = std::make_shared<Primary>();
        one_primary->type = Primary::Primary_t::NUMBER;
        one_primary->str_value = "1";
        one_primary->num_value = 1.0;
        result = one_primary;
    }

    // Apply sign
    if (term.sign < 0) {
        result = std::make_shared<Factor>('-', result);
    }

    return result;
}

///////////////////////////////////////////////////////////////////////////////
// GAAlgebraRules implementation
///////////////////////////////////////////////////////////////////////////////

int GAAlgebraRules::getCanonicalOrderPriority(const std::string& var)
{
    // Comprehensive canonical ordering for all GA algebras (EGA2D, EGA3D, PGA2DP, PGA3DP)
    // Lower numbers = higher priority (appear first in expressions)

    // GROUP 0: Numeric coefficients (0-4) - Always first
    try {
        return 1; // All numeric coefficients get priority 1
    }
    catch (...) {
        // Not a numeric coefficient, continue to other checks
    }

    // GROUP 1: Parenthesized rotor coefficient expressions (5-9) - After numeric
    // coefficients
    if (var.starts_with("(") && var.find("R.c") != std::string::npos) {
        return 5; // All parenthesized rotor expressions get priority 5
    }

    // GROUP 2: Simple rotor/motor coefficients (10-99) - After parenthesized expressions
    if (var.starts_with("R.c")) {
        try {
            int coeff_num = std::stoi(var.substr(3));
            return 10 + coeff_num; // R.c0=10, R.c1=11, R.c2=12, ..., R.c7=17
        }
        catch (...) {
            return 99; // fallback for malformed R.c variables
        }
    }

    // GROUP 3: Bivector components (1000-1099) - Before vectors in output expressions
    if (var.starts_with("B.")) {
        // Standard bivector components (for all algebras)
        if (var == "B.x") return 1000; // EGA3D/PGA bivector x-component
        if (var == "B.y") return 1001; // EGA3D/PGA bivector y-component
        if (var == "B.z") return 1002; // EGA3D/PGA bivector z-component
        if (var == "B.w") return 1003; // PGA bivector w-component

        // PGA3DP bivector components (moment + velocity)
        if (var == "B.vx") return 1010; // PGA3DP bivector velocity x
        if (var == "B.vy") return 1011; // PGA3DP bivector velocity y
        if (var == "B.vz") return 1012; // PGA3DP bivector velocity z
        if (var == "B.mx") return 1020; // PGA3DP bivector moment x
        if (var == "B.my") return 1021; // PGA3DP bivector moment y
        if (var == "B.mz") return 1022; // PGA3DP bivector moment z

        return 1099; // fallback for other B. variables
    }

    // GROUP 3: Vector components (2000-2099)
    if (var.starts_with("v.")) {
        if (var == "v.x") return 2000; // Vector x-component (all algebras)
        if (var == "v.y") return 2001; // Vector y-component (all algebras)
        if (var == "v.z") return 2002; // Vector z-component (EGA3D, PGA)
        if (var == "v.w") return 2003; // Vector w-component (PGA3DP)
        return 2099;                   // fallback for other v. variables
    }

    // GROUP 4: Trivector components (3000-3099) - PGA3DP only
    if (var.starts_with("t.")) {
        if (var == "t.x") return 3000; // Trivector x-component
        if (var == "t.y") return 3001; // Trivector y-component
        if (var == "t.z") return 3002; // Trivector z-component
        if (var == "t.w") return 3003; // Trivector w-component
        return 3099;                   // fallback for other t. variables
    }

    // GROUP 5: Everything else (5000+) - Alphabetical ordering
    return 5000 + static_cast<int>(var[0]);
}

std::vector<std::pair<std::string, int>>
GAAlgebraRules::getSortedVariablePairs(const std::map<std::string, int>& factors)
{

    // Convert map to vector for custom sorting
    std::vector<std::pair<std::string, int>> var_pairs(factors.begin(), factors.end());

    // Sort by canonical order priority, then alphabetically for ties
    std::sort(var_pairs.begin(), var_pairs.end(), [](const auto& a, const auto& b) {
        int order_a = getCanonicalOrderPriority(a.first);
        int order_b = getCanonicalOrderPriority(b.first);
        if (order_a != order_b) {
            return order_a < order_b;
        }
        return a.first < b.first; // alphabetical for ties
    });

    return var_pairs;
}

std::map<std::string, int>
GAAlgebraRules::reorderCommutativeFactors(const std::map<std::string, int>& factors)
{

    // Use the new canonical ordering system
    auto sorted_pairs = getSortedVariablePairs(factors);

    // Rebuild map (NOTE: this loses the ordering due to std::map's nature)
    // This function is kept for backward compatibility but toString methods
    // should use getSortedVariablePairs directly
    std::map<std::string, int> result;
    for (const auto& [var, power] : sorted_pairs) {
        result[var] = power;
    }

    return result;
}

bool GAAlgebraRules::isRotorCoefficient(const std::string& var)
{
    return var.starts_with("R.c");
}

bool GAAlgebraRules::isVectorComponent(const std::string& var)
{
    return var.starts_with("v.") || var.starts_with("B.") || var.starts_with("t.");
}

///////////////////////////////////////////////////////////////////////////////
// Test implementation for validation
///////////////////////////////////////////////////////////////////////////////

void SimplificationTests::testEGA2DExpansion()
{
    fmt::println("=== Testing EGA2D Expression Expansion ===");

    std::string input =
        "(R.c0 * v.x + R.c1 * v.y) * R.c0 + (R.c0 * v.y - R.c1 * v.x) * R.c1";
    fmt::println("Input: {}", input);

    Parser parser(input);
    auto ast = parser.parse();

    auto simplified = ExpressionSimplifier::simplify(ast);

    fmt::println("Simplified: {}", simplified->to_string());
    fmt::println("");
}