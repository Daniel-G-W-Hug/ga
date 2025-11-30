// Copyright 2024-2025, Daniel Hug. All rights reserved.

#ifndef GA_PRDXPR_METRIC_EXPORT_HPP
#define GA_PRDXPR_METRIC_EXPORT_HPP

#include "ga_prdxpr_rule_generator.hpp"
#include <mdspan>
#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////////////////
// Metric Export for ga_usr_consts.hpp Integration
////////////////////////////////////////////////////////////////////////////////

// Export extended metric matrix as C++ constexpr array code
// Generates code ready for copy-paste into ga_usr_consts.hpp
//
// Output format:
//   // Extended metric for ega3d (8×8 matrix, 64 elements)
//   // Basis: 1, e1, e2, e3, e23, e31, e12, e123
//   inline constexpr std::array<int, 64> ega3d_extended_metric = {
//       // Row 0 (1):
//       1, 0, 0, 0, 0, 0, 0, 0,
//       // Row 1 (e1):
//       0, 1, 0, 0, 0, 0, 0, 0,
//       ...
//   };
std::string export_extended_metric_to_code(std::string const& algebra_name,
                                           mvec_coeff const& basis,
                                           std::vector<int> const& matrix_data);

// Export regressive extended metric matrix as C++ constexpr array code
// Used for PGA algebras with bulk/weight dual operations
//
// Output format similar to extended metric export
std::string export_regressive_metric_to_code(std::string const& algebra_name,
                                             mvec_coeff const& basis,
                                             std::vector<int> const& matrix_data);

// Export both extended and regressive metrics for a PGA algebra
// Returns combined code block with both metrics
std::string export_pga_metrics_to_code(std::string const& algebra_name,
                                       mvec_coeff const& basis,
                                       std::vector<int> const& extended_metric_data,
                                       std::vector<int> const& regressive_metric_data);

// Print all metrics for an algebra to stdout (ready for copy-paste)
void print_metrics_for_algebra(AlgebraConfig const& config);

#endif // GA_PRDXPR_METRIC_EXPORT_HPP
