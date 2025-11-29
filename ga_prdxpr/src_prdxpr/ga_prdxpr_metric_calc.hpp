// Copyright 2024-2025, Daniel Hug. All rights reserved.

#ifndef GA_PRDXPR_METRIC_CALC_HPP
#define GA_PRDXPR_METRIC_CALC_HPP

#include "ga_prdxpr_rule_generator.hpp"
#include <mdspan>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// Integer Determinant Computation for Gram Matrices
////////////////////////////////////////////////////////////////////////////////

// Compute determinant of 2×2 integer matrix using direct formula
// Returns: det(M) = M[0,0]*M[1,1] - M[0,1]*M[1,0]
int compute_integer_determinant_2x2(
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M);

// Compute determinant of 3×3 integer matrix using direct formula (Sarrus rule)
// Returns: det(M) using cofactor expansion along first row
int compute_integer_determinant_3x3(
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M);

// Compute determinant of 4×4 integer matrix using cofactor expansion
// Returns: det(M) using expansion along first row
int compute_integer_determinant_4x4(
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M);

// Compute determinant of 5×5 integer matrix using cofactor expansion
// Returns: det(M) using expansion along first row
int compute_integer_determinant_5x5(
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M);

// Generic integer determinant dispatcher
// Calls appropriate specialized function based on matrix size
// Supports n = 1, 2, 3, 4, 5 (sufficient for CGA3D)
int compute_integer_determinant(
    std::mdspan<int const, std::extents<size_t, std::dynamic_extent, std::dynamic_extent>> const& M,
    int n);

////////////////////////////////////////////////////////////////////////////////
// Vector Metric and Gram Matrix Computation
////////////////////////////////////////////////////////////////////////////////

// Get metric value for two basis vector indices
// Handles both 1-based and 0-based indexing automatically
// For diagonal metric: returns config.metric_signature[i] when i == j, else 0
// For full metric (future CGA): returns config.metric_matrix[i][j]
int get_vector_metric(int index_i, int index_j, AlgebraConfig const& config);

// Compute determinant of Gram matrix for two blades of the same grade
// Gram[p,q] = metric(indices_a[p], indices_b[q])
// Returns: det(Gram) using integer arithmetic
int compute_gram_determinant(std::vector<int> const& indices_a,
                             std::vector<int> const& indices_b,
                             AlgebraConfig const& config);

// Compute inner product of two basis blades using Gram matrix determinant
// Returns: ⟨blade_a, blade_b⟩ = det(Gram) for same grade, 0 for different grades
int compute_blade_inner_product(std::string const& blade_a,
                                std::string const& blade_b,
                                AlgebraConfig const& config);

////////////////////////////////////////////////////////////////////////////////
// Full Extended Metric Matrix Computation
////////////////////////////////////////////////////////////////////////////////

// Calculate full extended metric matrix (including off-diagonal elements)
// Returns: n×n matrix in flattened row-major storage where n = basis.size()
// Mathematical definition: G[i,j] = ⟨basis[i], basis[j]⟩
// Properties:
//   - Symmetric: G[i,j] = G[j,i]
//   - Block-diagonal by grade (different grades → 0)
//   - Diagonal for orthogonal bases (EGA, PGA, STA)
//   - Non-diagonal for null bases (CGA with n₀, n∞)
std::vector<int> calculate_extended_metric_matrix_full(AlgebraConfig const& config);

////////////////////////////////////////////////////////////////////////////////
// Regressive Extended Metric Matrix Computation (for PGA)
////////////////////////////////////////////////////////////////////////////////

// Parse a signed basis element string to extract the basis element and sign
// Returns: pair of (basis_index, sign) where basis_index is position in basis vector
// Example: "-e23" → (index_of_e23, -1), "e1" → (index_of_e1, +1), "0" → (-1, 0)
std::pair<int, int> parse_signed_basis_element(std::string const& signed_element,
                                                 mvec_coeff const& basis);

// Build complement transformation matrix C from complement rules
// C[i,j] = coefficient of basis[j] in complement(basis[i])
// For orthogonal bases, C is a permutation matrix with ±1 entries
// Returns: n×n matrix in flattened row-major storage
std::vector<int> build_complement_matrix(mvec_coeff const& basis,
                                          prd_rules const& complement_rules);

// Compute matrix multiplication: R = A · B
// All matrices in flattened row-major storage
// Returns: n×n result matrix
std::vector<int> matrix_multiply(std::vector<int> const& A_data,
                                 std::vector<int> const& B_data,
                                 size_t n);

// Compute matrix transpose: R = A^T
// Returns: n×n transposed matrix in flattened row-major storage
std::vector<int> matrix_transpose(std::vector<int> const& A_data,
                                  size_t n);

// Compute triple matrix product: R = C · G · C^T
// Used for regressive metric: Ḡ = C · G · C^T
// Returns: n×n result matrix in flattened row-major storage
std::vector<int> matrix_triple_product(std::vector<int> const& C_data,
                                       std::vector<int> const& G_data,
                                       size_t n);

// Calculate regressive extended metric matrix Ḡ via complement transformation
// Mathematical definition: Ḡ = C · G · C^T
// where C is the complement transformation matrix
// Returns: n×n matrix in flattened row-major storage
//
// For diagonal PGA metrics, this reduces to: Ḡ[i,i] = 1 - G[i,i]
// But the general formula works for all cases including future non-diagonal
std::vector<int> calculate_regressive_extended_metric_matrix_full(
    AlgebraConfig const& config,
    prd_rules const& complement_rules);

#endif // GA_PRDXPR_METRIC_CALC_HPP
