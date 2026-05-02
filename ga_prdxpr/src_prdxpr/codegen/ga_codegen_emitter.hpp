#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "../ga_prdxpr_common.hpp"
#include "../ga_prdxpr_config_types.hpp"
#include "ga_codegen_types.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace codegen {

// Map a ProductConfig.product_name to the C++ function name that should be emitted.
// Returns nullopt for products that have no per-overload C++ implementation
// (wrappers, display-only "alternative" computations, sandwich placeholders, etc.).
std::optional<std::string>
product_to_cpp_function(std::string const& product_name);

// Extract the result-type filter name from a ProductCase.case_name string.
// Examples:
//   "mv * mv -> mv"            -> "mv"
//   "vec ^ bivec -> ps"        -> "ps"
//   "dot(bivec, bivec) -> s"   -> "s"
//   "bivec << ps -> vec"       -> "vec"
//   "vec ^ vec -> 0"           -> "0"   (zero-result marker, default scalar type)
//   "rdot(s, s) -> 0 ps"       -> "0"   (zero-result, explicit type = "ps")
// Use parse_explicit_zero_override() to recover the explicit type (if any).
// Throws std::invalid_argument if the "->" arrow is missing.
std::string parse_result_filter(std::string const& case_name);

// For case_names of the form "... -> 0 X" (zero result with explicit type X),
// returns X. Returns std::nullopt if the case_name does not have that form.
std::optional<std::string>
parse_explicit_zero_override(std::string const& case_name);

// Strip a single component suffix from an access string.
//   "A.c0"  -> "A"
//   "v.x"   -> "v"
//   "B.vx"  -> "B"
//   "s"     -> "s"   (no suffix to strip)
//   "0"     -> "0"   (zero literal stays as-is)
std::string strip_component(std::string const& access);

// Determine the variable name an operand uses, given:
//   - `coeff`         : the per-basis access strings (size = full algebra basis)
//   - `basis_indices` : positions within the basis that belong to the operand's filter
//
// Reads non-empty / non-"0" entries at the indicated positions, strips the component
// suffix, and asserts they all agree on a single name. Returns that name.
//
// Examples:
//   coeff = {"A.c0","A.c1","A.c2","A.c3"},  indices = {0,1,2,3}        -> "A"
//   coeff = {"s","v.x","v.y","ps"},         indices = {1,2}            -> "v"
//   coeff = {"R.c0","0","0","R.c1"},        indices = {0,3}            -> "R"
//   coeff = {"s","v.x","v.y","ps"},         indices = {0}              -> "s"
//
// Throws std::runtime_error if no usable entries exist or names disagree.
std::string param_name_from_coeff(mvec_coeff const& coeff,
                                  std::vector<std::size_t> const& basis_indices);

// Emit a complete C++ function for one (non-sandwich) OutputCase.
//
// Returns nullopt for cases that cannot or should not be emitted at this stage:
//  - non-primitive product (wrapper, alternative, ...): see product_to_cpp_function
//  - zero-result cases ("-> 0"): deferred to Phase 6
//  - missing/unknown filter or coefficient: skipped with a reason
// In skip cases, *skip_reason (if non-null) receives a one-line diagnostic.
//
// Otherwise returns the rendered function body, including a trailing newline.
std::optional<std::string>
emit_function(configurable::AlgebraData const& algebra,
              configurable::ProductConfig const& product,
              configurable::OutputCase const& case_def, mvec_coeff const& prd_mv,
              TypeRegistry const& registry, std::string* skip_reason = nullptr);

} // namespace codegen
