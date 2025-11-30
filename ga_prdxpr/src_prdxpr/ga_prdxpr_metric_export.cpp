// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_metric_export.hpp"
#include "ga_prdxpr_metric_calc.hpp"
#include <iostream>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
// Helper Functions
////////////////////////////////////////////////////////////////////////////////

// Calculate the maximum row data length for alignment
static size_t calculate_max_row_length(std::vector<int> const& matrix_data, size_t n)
{
    size_t max_length = 0;

    for (size_t row = 0; row < n; ++row) {
        std::ostringstream temp;
        temp << "    ";

        for (size_t col = 0; col < n; ++col) {
            size_t idx = row * n + col;
            temp << matrix_data[idx];

            if (row < n - 1 || col < n - 1) {
                temp << ", ";
            }
        }

        max_length = std::max(max_length, temp.str().length());
    }

    return max_length;
}

// Format matrix row with basis element comment
static std::string format_matrix_row(std::vector<int> const& matrix_data, size_t row,
                                     size_t n, std::string const& basis_element,
                                     bool is_last_row, size_t max_row_length)
{
    std::ostringstream oss;

    // Row data with trailing comment
    oss << "    ";

    for (size_t col = 0; col < n; ++col) {
        size_t idx = row * n + col;
        oss << matrix_data[idx];

        if (row < n - 1 || col < n - 1) { // Not the last element
            oss << ", ";
        }
    }

    // Calculate padding needed for alignment
    size_t current_length = oss.str().length();
    size_t padding =
        max_row_length > current_length ? max_row_length - current_length : 0;

    // Add padding and basis element comment
    oss << std::string(padding, ' ') << "  // " << basis_element;

    if (!is_last_row) {
        oss << "\n";
    }

    return oss.str();
}

////////////////////////////////////////////////////////////////////////////////
// Export Functions
////////////////////////////////////////////////////////////////////////////////

std::string export_extended_metric_to_code(std::string const& algebra_name,
                                           mvec_coeff const& basis,
                                           std::vector<int> const& matrix_data)
{
    size_t const n = basis.size();
    size_t const total_elements = n * n;

    std::ostringstream oss;

    // Header comment
    oss << "// extended metric for " << algebra_name << " (" << n << "×" << n
        << " matrix, " << total_elements << " elements)\n";

    // Basis comment
    oss << "// Basis: ";
    for (size_t i = 0; i < n; ++i) {
        oss << basis[i];
        if (i < n - 1) oss << ", ";
    }
    oss << "\n";

    // Array declaration
    oss << "inline constexpr std::array<int, " << total_elements << "> " << algebra_name
        << "_metric = {\n";

    // Calculate maximum row length for comment alignment
    size_t max_row_length = calculate_max_row_length(matrix_data, n);

    // Matrix data with row comments
    for (size_t row = 0; row < n; ++row) {
        bool is_last_row = (row == n - 1);
        oss << format_matrix_row(matrix_data, row, n, basis[row], is_last_row,
                                 max_row_length);
    }

    oss << "\n};\n";

    // Accessor function
    oss << "\n";
    oss << "// Accessor function for convenient 2D indexing: G[i,j]\n";
    oss << "inline constexpr auto " << algebra_name << "_metric_view() {\n";
    oss << "    return std::mdspan<int const, std::extents<size_t, " << n << ", " << n
        << ">>{\n";
    oss << "        " << algebra_name << "_metric.data()\n";
    oss << "    };\n";
    oss << "}\n";

    return oss.str();
}

std::string export_regressive_metric_to_code(std::string const& algebra_name,
                                             mvec_coeff const& basis,
                                             std::vector<int> const& matrix_data)
{
    size_t const n = basis.size();
    size_t const total_elements = n * n;

    std::ostringstream oss;

    // Header comment
    oss << "// Regressive extended metric for " << algebra_name << " (" << n << "×" << n
        << " matrix, " << total_elements << " elements)\n";

    // Basis comment
    oss << "// Basis: ";
    for (size_t i = 0; i < n; ++i) {
        oss << basis[i];
        if (i < n - 1) oss << ", ";
    }
    oss << "\n";

    // Array declaration
    oss << "inline constexpr std::array<int, " << total_elements << "> " << algebra_name
        << "_rmetric = {\n";

    // Calculate maximum row length for comment alignment
    size_t max_row_length = calculate_max_row_length(matrix_data, n);

    // Matrix data with row comments
    for (size_t row = 0; row < n; ++row) {
        bool is_last_row = (row == n - 1);
        oss << format_matrix_row(matrix_data, row, n, basis[row], is_last_row,
                                 max_row_length);
    }

    oss << "\n};\n";

    // Accessor function
    oss << "\n";
    oss << "// Accessor function for convenient 2D indexing: Ḡ[i,j]\n";
    oss << "inline constexpr auto " << algebra_name << "_rmetric_view() {\n";
    oss << "    return std::mdspan<int const, std::extents<size_t, " << n << ", " << n
        << ">>{\n";
    oss << "        " << algebra_name << "_rmetric.data()\n";
    oss << "    };\n";
    oss << "}\n";

    return oss.str();
}

std::string export_pga_metrics_to_code(std::string const& algebra_name,
                                       mvec_coeff const& basis,
                                       std::vector<int> const& extended_metric_data,
                                       std::vector<int> const& regressive_metric_data)
{
    std::ostringstream oss;

    oss << "/////////////////////////////////////////////////////////////////////////////"
           "///\n";
    oss << "// " << algebra_name << " - Projective Geometric Algebra Metrics\n";
    oss << "/////////////////////////////////////////////////////////////////////////////"
           "///\n";
    oss << "\n";

    // Extended metric
    oss << export_extended_metric_to_code(algebra_name, basis, extended_metric_data);
    oss << "\n";

    // Regressive metric
    oss << export_regressive_metric_to_code(algebra_name, basis, regressive_metric_data);

    return oss.str();
}

// Generate algebra name from basis prefix (e.g., "e" -> "algebra", "g" -> "sta")
static std::string generate_algebra_name(AlgebraConfig const& config)
{
    // Simple heuristic based on basis prefix and size
    if (config.basis_prefix == "g") {
        return "sta" + std::to_string(config.basis_vectors.size()) + "d";
    }
    else if (config.metric_signature.back() == 0) {
        // PGA (has null vector)
        int euclidean_dims = static_cast<int>(config.basis_vectors.size()) - 1;
        return "pga" + std::to_string(euclidean_dims) + "dp";
    }
    else {
        // EGA
        return "ega" + std::to_string(config.basis_vectors.size()) + "d";
    }
}

void print_metrics_for_algebra(AlgebraConfig const& config)
{
    auto const& basis = config.multivector_basis;
    std::string algebra_name = generate_algebra_name(config);

    std::cout << "///////////////////////////////////////////////////////////////////////"
                 "/////////\n";
    std::cout << "// extended metric for " << algebra_name << "\n";
    std::cout << "///////////////////////////////////////////////////////////////////////"
                 "/////////\n";
    std::cout << "\n";

    // Calculate extended metric
    auto G_data = calculate_extended_metric_matrix_full(config);

    // Export extended metric
    std::string extended_code =
        export_extended_metric_to_code(algebra_name, basis, G_data);
    std::cout << extended_code << "\n";

    // For PGA algebras, also export regressive metric
    bool is_pga =
        (config.metric_signature.back() == 0); // Last basis vector is degenerate

    if (is_pga) {
        // Generate complement rules for regressive metric calculation
        auto wedge_rules = generate_wedge_product_rules(config);
        auto complement_rules = generate_complement_rules(config, wedge_rules);

        // Calculate regressive metric
        auto G_reg_data =
            calculate_regressive_extended_metric_matrix_full(config, complement_rules);

        // Export regressive metric
        std::string regressive_code =
            export_regressive_metric_to_code(algebra_name, basis, G_reg_data);
        std::cout << regressive_code << "\n";
    }
}
