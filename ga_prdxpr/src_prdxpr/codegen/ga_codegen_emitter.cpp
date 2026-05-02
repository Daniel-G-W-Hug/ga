// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "codegen/ga_codegen_emitter.hpp"

#include <algorithm>
#include <cctype>
#include <fmt/core.h>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace codegen {

std::optional<std::string>
product_to_cpp_function(std::string const& product_name)
{
    static std::unordered_map<std::string, std::string> const m = {
        {"gpr", "operator*"},
        {"wdg", "wdg"},
        {"dot", "dot"},
        {"cmt", "cmt"},
        {"rgpr", "rgpr"},
        {"rwdg", "rwdg"},
        {"rdot", "rdot"},
        {"rcmt", "rcmt"},
        {"twdg1", "twdg1"},
        {"rtwdg1", "rtwdg1"},
        {"l_contract", "operator<<"},
        {"r_contract", "operator>>"},
    };
    auto it = m.find(product_name);
    if (it == m.end()) return std::nullopt;
    return it->second;
}

std::string parse_result_filter(std::string const& case_name)
{
    auto const arrow_pos = case_name.find("->");
    if (arrow_pos == std::string::npos) {
        throw std::invalid_argument("parse_result_filter: missing '->' in case_name: " +
                                    case_name);
    }
    // Skip past "->" and any whitespace.
    std::size_t i = arrow_pos + 2;
    while (i < case_name.size() && std::isspace(static_cast<unsigned char>(case_name[i]))) {
        ++i;
    }
    // Read the result-filter token: alnum + underscore.
    std::size_t j = i;
    while (j < case_name.size() &&
           (std::isalnum(static_cast<unsigned char>(case_name[j])) ||
            case_name[j] == '_')) {
        ++j;
    }
    if (j == i) {
        throw std::invalid_argument(
            "parse_result_filter: empty result token in case_name: " + case_name);
    }
    return case_name.substr(i, j - i);
}

std::optional<std::string>
parse_explicit_zero_override(std::string const& case_name)
{
    auto const arrow_pos = case_name.find("->");
    if (arrow_pos == std::string::npos) return std::nullopt;
    std::size_t i = arrow_pos + 2;
    auto skip_ws = [&]() {
        while (i < case_name.size() &&
               std::isspace(static_cast<unsigned char>(case_name[i]))) {
            ++i;
        }
    };
    skip_ws();
    // First token must be exactly "0".
    if (i >= case_name.size() || case_name[i] != '0') return std::nullopt;
    if (i + 1 < case_name.size() &&
        (std::isalnum(static_cast<unsigned char>(case_name[i + 1])) ||
         case_name[i + 1] == '_')) {
        return std::nullopt; // "0..." token is something other than just "0"
    }
    ++i;
    skip_ws();
    // Optional second token = explicit override type.
    std::size_t j = i;
    while (j < case_name.size() &&
           (std::isalnum(static_cast<unsigned char>(case_name[j])) ||
            case_name[j] == '_')) {
        ++j;
    }
    if (j == i) return std::nullopt; // bare "-> 0" form, no explicit type
    return case_name.substr(i, j - i);
}

std::string strip_component(std::string const& access)
{
    auto const dot_pos = access.find('.');
    if (dot_pos == std::string::npos) return access;
    return access.substr(0, dot_pos);
}

std::string param_name_from_coeff(mvec_coeff const& coeff,
                                  std::vector<std::size_t> const& basis_indices)
{
    std::string found;
    for (auto idx : basis_indices) {
        if (idx >= coeff.size()) {
            throw std::out_of_range(
                "param_name_from_coeff: basis index out of range: " +
                std::to_string(idx));
        }
        std::string const& v = coeff[idx];
        if (v.empty() || v == "0") continue;
        std::string name = strip_component(v);
        if (found.empty()) {
            found = name;
        }
        else if (found != name) {
            throw std::runtime_error(
                "param_name_from_coeff: inconsistent names within filter region: '" +
                found + "' vs '" + name + "'");
        }
    }
    if (found.empty()) {
        throw std::runtime_error(
            "param_name_from_coeff: no non-zero coefficient entries in filter region");
    }
    return found;
}

namespace {

// Return the expression to use for a given basis index in the result.
// `prd_mv` stores "0" for slots with no contribution; we convert that to "0.0"
// so it's a valid C++ literal.
std::string expr_at(mvec_coeff const& prd_mv, std::size_t idx)
{
    if (idx >= prd_mv.size()) return std::string("0.0");
    std::string const& s = prd_mv[idx];
    if (s.empty() || s == "0") return std::string("0.0");
    return s;
}

// Render a parameter declaration ("MVec3d<T> const& A" or "Scalar2d<T> s").
std::string render_param(TypeInfo const& info, std::string const& tparam,
                         std::string const& name)
{
    std::ostringstream os;
    os << info.cpp_type << "<" << tparam << ">";
    if (info.pass_by_value) {
        os << " " << name;
    }
    else {
        os << " const& " << name;
    }
    return os.str();
}

// Render the inline arguments for a Named/Indexed/SingleValue (sub-)type by reading
// `prd_mv` at the given basis indices (in order) and joining them with ", ".
std::string render_args(mvec_coeff const& prd_mv,
                        std::vector<std::size_t> const& basis_indices)
{
    std::string out;
    for (std::size_t k = 0; k < basis_indices.size(); ++k) {
        if (k > 0) out += ", ";
        out += expr_at(prd_mv, basis_indices[k]);
    }
    return out;
}

// Build the body's return statement for an Indexed result, using temp-vars
// (ctype const cN = ...; ... return Type<ctype>(c0, ..., cN-1);).
std::string build_indexed_body(TypeInfo const& result_info, mvec_coeff const& prd_mv)
{
    std::ostringstream os;
    auto const& idx = result_info.basis_indices;
    for (std::size_t k = 0; k < idx.size(); ++k) {
        os << "    ctype const c" << k << " = " << expr_at(prd_mv, idx[k]) << ";\n";
    }
    os << "    return " << result_info.cpp_type << "<ctype>(";
    for (std::size_t k = 0; k < idx.size(); ++k) {
        if (k > 0) os << ", ";
        os << "c" << k;
    }
    os << ");\n";
    return os.str();
}

// Build the body's return statement for a Composite result, using nested constructors.
std::string build_composite_body(TypeInfo const& result_info,
                                 TypeRegistry const& registry,
                                 mvec_coeff const& prd_mv)
{
    std::ostringstream os;
    os << "    return " << result_info.cpp_type << "<ctype>(";
    for (std::size_t k = 0; k < result_info.sub_filter_names.size(); ++k) {
        if (k > 0) os << ", ";
        auto const& sub = registry.get(result_info.sub_filter_names[k]);
        os << sub.cpp_type << "<ctype>(" << render_args(prd_mv, sub.basis_indices)
           << ")";
    }
    os << ");\n";
    return os.str();
}

// Build the body's return statement for Named or SingleValue results
// (no temp-vars; expressions inline in the constructor call).
std::string build_inline_body(TypeInfo const& result_info, mvec_coeff const& prd_mv)
{
    std::ostringstream os;
    os << "    return " << result_info.cpp_type << "<ctype>("
       << render_args(prd_mv, result_info.basis_indices) << ");\n";
    return os.str();
}

// True if every basis-index slot of `result_info` in `prd_mv` is empty or "0".
// Drives [[maybe_unused]] zero-result emission.
bool is_zero_result(TypeInfo const& result_info, mvec_coeff const& prd_mv)
{
    auto walk = [&](std::size_t idx) {
        if (idx >= prd_mv.size()) return true;
        std::string const& s = prd_mv[idx];
        return s.empty() || s == "0";
    };
    if (result_info.storage == StorageKind::Composite) {
        // Composite spans the union of sub-types; the flat basis_indices already
        // contains every position we care about.
        for (auto idx : result_info.basis_indices) {
            if (!walk(idx)) return false;
        }
        return true;
    }
    for (auto idx : result_info.basis_indices) {
        if (!walk(idx)) return false;
    }
    return true;
}

// Render an unnamed [[maybe_unused]] parameter declaration for a zero-result
// function. Examples:
//   PScalar2d<T>      ->  "[[maybe_unused]] PScalar2d<T>"
//   Vec2d<T> const&   ->  "[[maybe_unused]] Vec2d<T> const&"
std::string render_param_anon(TypeInfo const& info, std::string const& tparam)
{
    std::ostringstream os;
    os << "[[maybe_unused]] " << info.cpp_type << "<" << tparam << ">";
    if (!info.pass_by_value) os << " const&";
    return os.str();
}

// Wrap bare scalar/pscalar identifier(s) in expressions with `ctype(name)` so
// that strong-typed Scalar/PScalar values participate correctly in arithmetic
// alongside underlying `ctype` values. Without this, code like
//
//     ctype const c0 = s * M.c0;
//
// fails to compile because `s` (Scalar2d<T>) doesn't implicitly convert to
// `ctype`. The conventional fix is `ctype(s) * M.c0`.
//
// `bare_names` is the set of bare identifiers that need wrapping. Walks each
// expression in `prd_mv` and replaces \b<name>\b with `ctype(<name>)` exactly
// when not already wrapped.
mvec_coeff wrap_bare_scalars(mvec_coeff const& prd_mv,
                             std::vector<std::string> const& bare_names)
{
    if (bare_names.empty()) return prd_mv;
    mvec_coeff out = prd_mv;
    for (auto& expr : out) {
        for (auto const& name : bare_names) {
            // Replace each whole-word occurrence of `name` with `ctype(name)`,
            // but skip any that are already inside a `ctype(...)` wrap.
            std::string updated;
            std::size_t pos = 0;
            while (pos < expr.size()) {
                auto found = expr.find(name, pos);
                if (found == std::string::npos) {
                    updated.append(expr, pos, std::string::npos);
                    break;
                }
                // Word boundary checks.
                auto is_word_char = [](char c) {
                    return std::isalnum(static_cast<unsigned char>(c)) || c == '_';
                };
                bool left_ok =
                    (found == 0) || !is_word_char(expr[found - 1]);
                bool right_ok = (found + name.size() == expr.size()) ||
                                !is_word_char(expr[found + name.size()]);
                if (!left_ok || !right_ok) {
                    updated.append(expr, pos, found - pos + 1);
                    pos = found + 1;
                    continue;
                }
                // Skip if already inside a "ctype(" wrap.
                if (found >= 6 && expr.compare(found - 6, 6, "ctype(") == 0) {
                    updated.append(expr, pos, found - pos + name.size());
                    pos = found + name.size();
                    continue;
                }
                updated.append(expr, pos, found - pos);
                updated += "ctype(";
                updated += name;
                updated += ")";
                pos = found + name.size();
            }
            expr = std::move(updated);
        }
    }
    return out;
}

// Detect bare-scalar parameter names by inspecting an operand's coefficient
// definition restricted to the operand's filter region. SingleValue operands
// have exactly one bare identifier (no '.') at their basis index.
std::vector<std::string>
collect_bare_scalars(mvec_coeff const& lhs_coeff, TypeInfo const& lhs_info,
                     mvec_coeff const& rhs_coeff, TypeInfo const& rhs_info)
{
    std::vector<std::string> out;
    auto add_if_bare = [&](mvec_coeff const& coeff,
                           std::vector<std::size_t> const& indices, bool is_single) {
        if (!is_single) return;
        if (indices.size() != 1) return;
        auto idx = indices[0];
        if (idx >= coeff.size()) return;
        std::string const& s = coeff[idx];
        if (s.empty() || s == "0") return;
        if (s.find('.') != std::string::npos) return; // not a bare scalar
        out.push_back(s);
    };
    add_if_bare(lhs_coeff, lhs_info.basis_indices,
                lhs_info.storage == StorageKind::SingleValue);
    add_if_bare(rhs_coeff, rhs_info.basis_indices,
                rhs_info.storage == StorageKind::SingleValue);
    return out;
}

} // namespace

std::optional<std::string>
emit_function(configurable::AlgebraData const& algebra,
              configurable::ProductConfig const& product,
              configurable::OutputCase const& case_def, mvec_coeff const& prd_mv,
              TypeRegistry const& registry, std::string* skip_reason)
{
    auto set_skip = [&](std::string const& msg) {
        if (skip_reason) *skip_reason = msg;
    };

    // 1. Function name.
    auto func_name = product_to_cpp_function(product.product_name);
    if (!func_name) {
        set_skip("non-primitive product '" + product.product_name + "'");
        return std::nullopt;
    }

    // 2. Result filter. A case_name like "vec ^ ps -> 0" marks an explicit zero
    //    result; the convention is to return Scalar of the algebra in that case.
    //    A case_name like "rdot(s, s) -> 0 ps" specifies the zero-result type
    //    explicitly (here: PScalar).
    std::string const result_filter_raw = parse_result_filter(case_def.case_name);
    std::string result_filter = result_filter_raw;
    if (result_filter_raw == "0") {
        auto override_type = parse_explicit_zero_override(case_def.case_name);
        result_filter = override_type.value_or("s");
    }

    // 3. TypeInfo lookups.
    if (!registry.has(case_def.left_filter_name)) {
        set_skip("unknown left filter '" + case_def.left_filter_name + "'");
        return std::nullopt;
    }
    if (!registry.has(case_def.right_filter_name)) {
        set_skip("unknown right filter '" + case_def.right_filter_name + "'");
        return std::nullopt;
    }
    if (!registry.has(result_filter)) {
        set_skip("unknown result filter '" + result_filter + "'");
        return std::nullopt;
    }
    auto const& lhs_info = registry.get(case_def.left_filter_name);
    auto const& rhs_info = registry.get(case_def.right_filter_name);
    auto const& result_info = registry.get(result_filter);

    // 4. Coefficient lookups.
    auto lit = algebra.coefficients.find(case_def.left_coeff_name);
    auto rit = algebra.coefficients.find(case_def.right_coeff_name);
    if (lit == algebra.coefficients.end()) {
        set_skip("unknown left coeff '" + case_def.left_coeff_name + "'");
        return std::nullopt;
    }
    if (rit == algebra.coefficients.end()) {
        set_skip("unknown right coeff '" + case_def.right_coeff_name + "'");
        return std::nullopt;
    }

    // 5. Wrap bare Scalar/PScalar operand names in `ctype(...)` so the
    //    strong-typed values participate correctly in arithmetic (cf. comment on
    //    wrap_bare_scalars).
    auto const bare_names =
        collect_bare_scalars(lit->second, lhs_info, rit->second, rhs_info);
    mvec_coeff const prd_mv_wrapped = wrap_bare_scalars(prd_mv, bare_names);

    // 6. Detect zero-result and emit the [[maybe_unused]] form. This catches both
    //    explicit "-> 0" markers AND cases where the case_name advertises a
    //    non-zero result type but the computed coefficients are all zero (degenerate
    //    products).
    bool const all_zero =
        (result_filter_raw == "0") || is_zero_result(result_info, prd_mv_wrapped);

    std::ostringstream os;
    os << "template <typename T, typename U>\n";
    os << "    requires(numeric_type<T> && numeric_type<U>)\n";

    if (all_zero) {
        os << "constexpr " << result_info.cpp_type << "<std::common_type_t<T, U>> "
           << *func_name << "(" << render_param_anon(lhs_info, "T") << ", "
           << render_param_anon(rhs_info, "U") << ")\n";
        os << "{\n";
        os << "    using ctype = std::common_type_t<T, U>;\n";
        os << "    return " << result_info.cpp_type << "<ctype>(0.0);\n";
        os << "}\n";
        return os.str();
    }

    // 7. Non-zero path: extract parameter names from coefficient definitions.
    std::string lhs_param;
    std::string rhs_param;
    try {
        lhs_param = param_name_from_coeff(lit->second, lhs_info.basis_indices);
        rhs_param = param_name_from_coeff(rit->second, rhs_info.basis_indices);
    }
    catch (std::exception const& e) {
        set_skip(std::string("param-name extraction failed: ") + e.what());
        return std::nullopt;
    }

    os << "constexpr " << result_info.cpp_type << "<std::common_type_t<T, U>> "
       << *func_name << "(" << render_param(lhs_info, "T", lhs_param) << ", "
       << render_param(rhs_info, "U", rhs_param) << ")\n";
    os << "{\n";
    os << "    using ctype = std::common_type_t<T, U>;\n";

    switch (result_info.storage) {
        case StorageKind::Indexed:
        case StorageKind::Named:
            // Both flat-constructor cases emit temp-vars + Type(c0, c1, ...) for
            // uniform readability across the codebase. The compiler optimizer
            // collapses the temporaries.
            os << build_indexed_body(result_info, prd_mv_wrapped);
            break;
        case StorageKind::Composite:
            os << build_composite_body(result_info, registry, prd_mv_wrapped);
            break;
        case StorageKind::SingleValue:
            os << build_inline_body(result_info, prd_mv_wrapped);
            break;
    }

    os << "}\n";
    return os.str();
}

} // namespace codegen
