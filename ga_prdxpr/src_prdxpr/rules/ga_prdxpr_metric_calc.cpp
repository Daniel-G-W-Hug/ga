// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "rules/ga_prdxpr_metric_calc.hpp"
#include "rules/ga_prdxpr_rule_generator.hpp" // For calculate_extended_metric()
#include <algorithm>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
// Integer Determinant Computation for Gram Matrices
////////////////////////////////////////////////////////////////////////////////

int compute_integer_determinant_2x2(
    std::mdspan<int const,
                std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M)
{
    return M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0];
}

int compute_integer_determinant_3x3(
    std::mdspan<int const,
                std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M)
{
    return M[0, 0] * (M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]) -
           M[0, 1] * (M[1, 0] * M[2, 2] - M[1, 2] * M[2, 0]) +
           M[0, 2] * (M[1, 0] * M[2, 1] - M[1, 1] * M[2, 0]);
}

int compute_integer_determinant_4x4(
    std::mdspan<int const,
                std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M)
{
    // Cofactor expansion along first row
    int det = 0;

    // Build 3×3 minors for each column
    for (int col = 0; col < 4; ++col) {
        // Build 3×3 minor by excluding row 0 and column col
        std::vector<int> minor_data(9, 0);
        std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>>
            minor{minor_data.data(), 3, 3};

        int minor_row = 0;
        for (int i = 1; i < 4; ++i) { // Skip row 0
            int minor_col = 0;
            for (int j = 0; j < 4; ++j) {
                if (j != col) { // Skip column col
                    minor[minor_row, minor_col] = M[i, j];
                    ++minor_col;
                }
            }
            ++minor_row;
        }

        // Compute minor determinant
        int minor_det = compute_integer_determinant_3x3(minor);

        // Add cofactor: (-1)^(0+col) * M[0,col] * minor_det
        int sign = (col % 2 == 0) ? 1 : -1;
        det += sign * M[0, col] * minor_det;
    }

    return det;
}

int compute_integer_determinant_5x5(
    std::mdspan<int const,
                std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M)
{
    // Cofactor expansion along first row
    int det = 0;

    // Build 4×4 minors for each column
    for (int col = 0; col < 5; ++col) {
        // Build 4×4 minor by excluding row 0 and column col
        std::vector<int> minor_data(16, 0);
        std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>>
            minor{minor_data.data(), 4, 4};

        int minor_row = 0;
        for (int i = 1; i < 5; ++i) { // Skip row 0
            int minor_col = 0;
            for (int j = 0; j < 5; ++j) {
                if (j != col) { // Skip column col
                    minor[minor_row, minor_col] = M[i, j];
                    ++minor_col;
                }
            }
            ++minor_row;
        }

        // Compute minor determinant
        int minor_det = compute_integer_determinant_4x4(minor);

        // Add cofactor: (-1)^(0+col) * M[0,col] * minor_det
        int sign = (col % 2 == 0) ? 1 : -1;
        det += sign * M[0, col] * minor_det;
    }

    return det;
}

int compute_integer_determinant(
    std::mdspan<int const,
                std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M,
    int n)
{
    if (n == 1) return M[0, 0];
    if (n == 2) return compute_integer_determinant_2x2(M);
    if (n == 3) return compute_integer_determinant_3x3(M);
    if (n == 4) return compute_integer_determinant_4x4(M);
    if (n == 5) return compute_integer_determinant_5x5(M);

    // Should not reach here for standard GA (max grade 5 for CGA3D)
    throw std::runtime_error("Integer determinant not implemented for n > 5");
}

////////////////////////////////////////////////////////////////////////////////
// Vector Metric and Gram Matrix Computation
////////////////////////////////////////////////////////////////////////////////

int get_vector_metric(int index_i, int index_j, AlgebraConfig const& config)
{
    // Detect 1-based vs 0-based indexing by inspecting basis vectors
    bool is_one_based = true;
    for (auto const& bv : config.basis_vectors) {
        auto bv_indices = parse_indices(bv, config.basis_prefix);
        if (!bv_indices.empty() && bv_indices[0] == 0) {
            is_one_based = false;
            break;
        }
    }

    // Convert to metric array indices
    int metric_idx_i = is_one_based ? (index_i - 1) : index_i;
    int metric_idx_j = is_one_based ? (index_j - 1) : index_j;

    // Same vector: use diagonal metric
    if (metric_idx_i == metric_idx_j) {
        if (metric_idx_i >= 0 &&
            metric_idx_i < static_cast<int>(config.metric_signature.size())) {
            return config.metric_signature[metric_idx_i];
        }
        return 0; // Degenerate direction
    }

    // Different vectors: check if full metric matrix provided (future CGA support)
    // For now, all supported algebras have orthogonal bases → off-diagonal = 0
    // TODO: When AlgebraConfig gets metric_matrix field, use:
    // if (config.has_full_metric()) {
    //     return (*config.metric_matrix)[metric_idx_i][metric_idx_j];
    // }

    // Orthogonal basis (current algebras): off-diagonal = 0
    return 0;
}

int compute_gram_determinant(std::vector<int> const& indices_a,
                             std::vector<int> const& indices_b,
                             AlgebraConfig const& config)
{
    int k = static_cast<int>(indices_a.size());

    if (k != static_cast<int>(indices_b.size())) {
        return 0; // Different grades → zero inner product
    }

    if (k == 0) {
        return 1; // Scalars: ⟨1, 1⟩ = 1
    }

    if (k == 1) {
        // Vectors: use vector metric directly
        return get_vector_metric(indices_a[0], indices_b[0], config);
    }

    // Build k×k Gram matrix: Gram[p,q] = metric(a[p], b[q])
    std::vector<int> gram_data(k * k, 0);
    std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> Gram{
        gram_data.data(), static_cast<size_t>(k), static_cast<size_t>(k)};

    for (int p = 0; p < k; ++p) {
        for (int q = 0; q < k; ++q) {
            Gram[p, q] = get_vector_metric(indices_a[p], indices_b[q], config);
        }
    }

    // Return determinant using integer arithmetic
    return compute_integer_determinant(Gram, k);
}

int compute_blade_inner_product(std::string const& blade_a, std::string const& blade_b,
                                AlgebraConfig const& config)
{
    // Special case: scalars
    if (blade_a == config.scalar_name && blade_b == config.scalar_name) {
        return 1; // ⟨1, 1⟩ = 1
    }
    if (blade_a == config.scalar_name || blade_b == config.scalar_name) {
        return 0; // ⟨scalar, vector⟩ = 0
    }

    // Parse indices from blade names (e.g., "e12" → {1, 2})
    auto indices_a = parse_indices(blade_a, config.basis_prefix);
    auto indices_b = parse_indices(blade_b, config.basis_prefix);

    // Different grades → zero inner product
    if (indices_a.size() != indices_b.size()) {
        return 0;
    }

    // Same grade: compute Gram matrix determinant
    return compute_gram_determinant(indices_a, indices_b, config);
}

////////////////////////////////////////////////////////////////////////////////
// Full Extended Metric Matrix Computation
////////////////////////////////////////////////////////////////////////////////

std::vector<int> calculate_extended_metric_matrix_full(AlgebraConfig const& config)
{
    auto const& basis = config.multivector_basis;
    size_t const n = basis.size();

    // Flattened row-major storage: matrix[i,j] = data[i*n + j]
    std::vector<int> matrix_data(n * n, 0);

    // Create mdspan view with square bracket access: G[i,j]
    std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> G{
        matrix_data.data(), n, n};

    // Get the correct diagonal values from calculate_extended_metric()
    // This uses the STA4D special case logic and is known to be correct
    auto diagonal_values = calculate_extended_metric(config);

    // Compute ALL elements G[i,j] = ⟨basis[i], basis[j]⟩
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) { // Upper triangular (symmetric matrix)

            int metric_value;

            if (i == j) {
                // Use diagonal calculation as source of truth
                // This handles STA4D special cases correctly
                metric_value = diagonal_values[i];
            }
            else {
                // Off-diagonal: compute inner product of different basis blades
                std::string const& blade_i = basis[i];
                std::string const& blade_j = basis[j];
                metric_value = compute_blade_inner_product(blade_i, blade_j, config);
            }

            G[i, j] = metric_value;
            G[j, i] = metric_value; // Symmetric: G[j,i] = G[i,j]
        }
    }

    return matrix_data;
}

////////////////////////////////////////////////////////////////////////////////
// Regressive Extended Metric Matrix Computation (for PGA)
////////////////////////////////////////////////////////////////////////////////

std::pair<int, int> parse_signed_basis_element(std::string const& signed_element,
                                               mvec_coeff const& basis)
{
    // Handle zero case
    if (signed_element == zero_str() || signed_element == "0") {
        return {-1, 0};
    }

    // Check for minus sign
    bool has_minus = (signed_element.find(minus_str()) == 0);
    std::string element =
        has_minus ? signed_element.substr(minus_str().length()) : signed_element;

    // Find element in basis
    auto it = std::find(basis.begin(), basis.end(), element);
    if (it == basis.end()) {
        return {-1, 0}; // Not found
    }

    int basis_index = static_cast<int>(std::distance(basis.begin(), it));
    int sign = has_minus ? -1 : +1;

    return {basis_index, sign};
}

std::vector<int> build_cmpl_matrix(mvec_coeff const& basis,
                                   prd_rules const& complement_rules)
{
    size_t const n = basis.size();

    // Initialize n×n matrix with zeros
    std::vector<int> C_data(n * n, 0);
    std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> C{
        C_data.data(), n, n};

    // Build matrix representation: C[i,j] = coeff of basis[j] in complement(basis[i])
    for (size_t i = 0; i < n; ++i) {
        std::string const& elem = basis[i];

        auto it = complement_rules.find(elem);
        if (it != complement_rules.end()) {
            std::string const& complement_value = it->second;

            // Parse the complement result to find which basis element and sign
            auto [basis_idx, sign] = parse_signed_basis_element(complement_value, basis);

            if (basis_idx >= 0) {
                C[i, basis_idx] = sign;
            }
            // If basis_idx == -1, complement is zero (C row stays all zeros)
        }
    }

    return C_data;
}

std::vector<int> matrix_multiply(std::vector<int> const& A_data,
                                 std::vector<int> const& B_data, size_t n)
{
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>>
        A{A_data.data(), n, n};
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>>
        B{B_data.data(), n, n};

    // Result matrix
    std::vector<int> R_data(n * n, 0);
    std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> R{
        R_data.data(), n, n};

    // Standard matrix multiplication: R[i,j] = Σ A[i,k] * B[k,j]
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            int sum = 0;
            for (size_t k = 0; k < n; ++k) {
                sum += A[i, k] * B[k, j];
            }
            R[i, j] = sum;
        }
    }

    return R_data;
}

std::vector<int> matrix_transpose(std::vector<int> const& A_data, size_t n)
{
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>>
        A{A_data.data(), n, n};

    // Result matrix
    std::vector<int> R_data(n * n, 0);
    std::mdspan<int, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> R{
        R_data.data(), n, n};

    // Transpose: R[i,j] = A[j,i]
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            R[i, j] = A[j, i];
        }
    }

    return R_data;
}

std::vector<int> matrix_triple_product(std::vector<int> const& C_data,
                                       std::vector<int> const& G_data, size_t n)
{
    // Compute: R = C · G · C^T

    // Step 1: Compute C^T
    auto CT_data = matrix_transpose(C_data, n);

    // Step 2: Compute T = G · C^T
    auto T_data = matrix_multiply(G_data, CT_data, n);

    // Step 3: Compute R = C · T = C · (G · C^T)
    auto R_data = matrix_multiply(C_data, T_data, n);

    return R_data;
}

std::vector<int>
calculate_regressive_extended_metric_matrix_full(AlgebraConfig const& config,
                                                 prd_rules const& complement_rules)
{
    auto const& basis = config.multivector_basis;
    size_t const n = basis.size();

    // Step 1: Calculate standard extended metric G
    auto G_data = calculate_extended_metric_matrix_full(config);

    // Step 2: Build complement transformation matrix C
    auto C_data = build_cmpl_matrix(basis, complement_rules);

    // Step 3: Compute regressive metric Ḡ = C · G · C^T
    auto G_reg_data = matrix_triple_product(C_data, G_data, n);

    return G_reg_data;
}
