// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_dual_calc.hpp"
#include "ga_prdxpr_metric_calc.hpp"
#include <mdspan>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
// Helper Functions for Dual Expression Construction
////////////////////////////////////////////////////////////////////////////////

// Build a dual expression string from matrix-vector multiplication result
// Constructs: coefficient * basis_element (with sign handling)
static std::string build_dual_expression(
    std::vector<int> const& coefficients,
    mvec_coeff const& basis,
    prd_rules const& complement_rules)
{
    std::vector<std::string> terms;

    for (size_t j = 0; j < coefficients.size(); ++j) {
        int coeff = coefficients[j];

        if (coeff == 0) continue;  // Skip zero terms

        // Get complement of basis[j]
        auto it = complement_rules.find(basis[j]);
        if (it == complement_rules.end()) continue;

        std::string const& complement_value = it->second;

        // Handle zero complement
        if (complement_value == zero_str() || complement_value == "0") {
            continue;
        }

        // Build term: coeff * complement(basis[j])
        std::string term;

        if (coeff == 1) {
            term = complement_value;
        } else if (coeff == -1) {
            // Handle sign: if complement already has minus, remove it; else add minus
            if (complement_value.find(minus_str()) == 0) {
                term = complement_value.substr(minus_str().length());
            } else {
                term = minus_str() + complement_value;
            }
        } else {
            // Non-unit coefficient (shouldn't happen for orthogonal bases)
            std::ostringstream oss;
            oss << coeff << "*" << complement_value;
            term = oss.str();
        }

        terms.push_back(term);
    }

    // Combine terms
    if (terms.empty()) {
        return zero_str();
    }

    if (terms.size() == 1) {
        return terms[0];
    }

    // Multiple terms: join with " + " (signs already in terms)
    std::ostringstream result;
    for (size_t i = 0; i < terms.size(); ++i) {
        if (i > 0) {
            // Add separator between terms
            if (terms[i].find(minus_str()) == 0) {
                result << " " << terms[i];  // Minus sign already present
            } else {
                result << " + " << terms[i];
            }
        } else {
            result << terms[i];
        }
    }

    return result.str();
}

////////////////////////////////////////////////////////////////////////////////
// Matrix-Vector Multiplication for Dual Calculation
////////////////////////////////////////////////////////////////////////////////

// Compute matrix-vector product: result[i] = Σⱼ G[i,j] * v[j]
// Where v[j] represents the complement transformation
static std::vector<int> matrix_vector_multiply(
    std::vector<int> const& G_data,
    std::vector<int> const& v_data,
    size_t n)
{
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>>
        G{G_data.data(), n, n};

    std::vector<int> result(n, 0);

    for (size_t i = 0; i < n; ++i) {
        int sum = 0;
        for (size_t j = 0; j < n; ++j) {
            sum += G[i, j] * v_data[j];
        }
        result[i] = sum;
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
// Dual Rule Calculation
////////////////////////////////////////////////////////////////////////////////

prd_rules calculate_dual_rules(
    AlgebraConfig const& config,
    std::vector<int> const& G_data,
    prd_rules const& complement_rules)
{
    auto const& basis = config.multivector_basis;
    size_t const n = basis.size();

    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>>
        G{G_data.data(), n, n};

    prd_rules dual_rules;

    // For each basis element i: dual(basis[i]) = Σⱼ G[i,j] * complement(basis[j])
    for (size_t i = 0; i < n; ++i) {
        std::string const& elem = basis[i];

        // Extract row i from G matrix
        std::vector<int> row_coeffs(n);
        for (size_t j = 0; j < n; ++j) {
            row_coeffs[j] = G[i, j];
        }

        // Build dual expression from coefficients and complements
        std::string dual_expr = build_dual_expression(row_coeffs, basis, complement_rules);

        dual_rules[elem] = dual_expr;
    }

    return dual_rules;
}

prd_rules calculate_regressive_dual_rules(
    AlgebraConfig const& config,
    std::vector<int> const& G_reg_data,
    prd_rules const& complement_rules)
{
    auto const& basis = config.multivector_basis;
    size_t const n = basis.size();

    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>>
        G_reg{G_reg_data.data(), n, n};

    prd_rules regressive_dual_rules;

    // For each basis element i: regressive_dual(basis[i]) = Σⱼ Ḡ[i,j] * complement(basis[j])
    for (size_t i = 0; i < n; ++i) {
        std::string const& elem = basis[i];

        // Extract row i from Ḡ matrix
        std::vector<int> row_coeffs(n);
        for (size_t j = 0; j < n; ++j) {
            row_coeffs[j] = G_reg[i, j];
        }

        // Build regressive dual expression from coefficients and complements
        std::string dual_expr = build_dual_expression(row_coeffs, basis, complement_rules);

        regressive_dual_rules[elem] = dual_expr;
    }

    return regressive_dual_rules;
}
