// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "codegen/ga_codegen_emitter.hpp"
#include "codegen/ga_codegen_types.hpp"
#include "generator/ga_prdxpr_generator.hpp"
#include "sandwich/ga_prdxpr_sandwich_simplifier.hpp"
#include <cctype>
#include <fmt/core.h>
#include <limits>
#include <map>
#include <optional>
#include <regex>
#include <set>
#include <stdexcept>
#include <string_view>

// Include mathematical definitions
#include "algebras/ga_prdxpr_ega2d.hpp"
#include "algebras/ga_prdxpr_ega3d.hpp"
#include "algebras/ga_prdxpr_pga2dp.hpp"
#include "algebras/ga_prdxpr_pga3dp.hpp"
#include "algebras/ga_prdxpr_sta4d.hpp"

using namespace configurable;

////////////////////////////////////////////////////////////////////////////////
// User-config validation: warn (don't fail) on inconsistent OutputCase entries.
//
// Catches three kinds of mistakes in algebra-config tuples like
//   {"mv_e * mv_e -> mv_e", "A_even", "B_even", "mv_e", "mv_e"}:
//
//   A. case_name's LHS/RHS tokens don't match left_filter_name/right_filter_name.
//   B. case_name's result type is not a known filter in the algebra (typo, or
//      filter from another algebra).
//   C. The product's actually-computed components fall outside the declared
//      result type's basis support (e.g. "mv_e * mv_e -> vec" produces a
//      scalar+bivec, both of which lie outside `vec`'s basis indices).
//
// Warnings go to stderr so they don't pollute stdout consumed by callers
// (--output=code, validation pipelines). Generation continues regardless.
////////////////////////////////////////////////////////////////////////////////
namespace {

std::string trim_ws(std::string_view sv)
{
    auto const b = sv.find_first_not_of(" \t");
    if (b == std::string_view::npos) return "";
    auto const e = sv.find_last_not_of(" \t");
    return std::string(sv.substr(b, e - b + 1));
}

bool is_zero_expr(std::string const& s)
{
    auto const t = trim_ws(s);
    return t.empty() || t == "0" || t == "0.0";
}

struct CaseTokens {
    std::string lhs;
    std::string rhs;
    std::string result; // bare token after "->"; "0" for zero-result cases
};

// Parse "LHS OP RHS -> RESULT" (OP in *, ^, <<, >>) or
// "func(LHS,RHS) -> RESULT" / "func(LHS, RHS) -> RESULT".
// Returns nullopt on parse failure (validation only warns, never throws).
std::optional<CaseTokens> parse_case_tokens(std::string const& case_name)
{
    auto const arrow = case_name.find("->");
    if (arrow == std::string::npos) return std::nullopt;

    // Result-side token: alnum + underscore, after optional whitespace.
    std::size_t i = arrow + 2;
    while (i < case_name.size() &&
           std::isspace(static_cast<unsigned char>(case_name[i]))) {
        ++i;
    }
    std::size_t j = i;
    while (
        j < case_name.size() &&
        (std::isalnum(static_cast<unsigned char>(case_name[j])) || case_name[j] == '_')) {
        ++j;
    }
    if (j == i) return std::nullopt;

    CaseTokens out;
    out.result = case_name.substr(i, j - i);

    std::string const lhs_part = trim_ws(std::string_view(case_name).substr(0, arrow));

    // Functional form: "func(LHS,RHS)".
    auto const lparen = lhs_part.find('(');
    auto const rparen = lhs_part.rfind(')');
    if (lparen != std::string::npos && rparen != std::string::npos && rparen > lparen) {
        std::string const inside = lhs_part.substr(lparen + 1, rparen - lparen - 1);
        auto const comma = inside.find(',');
        if (comma == std::string::npos) return std::nullopt;
        out.lhs = trim_ws(std::string_view(inside).substr(0, comma));
        out.rhs = trim_ws(std::string_view(inside).substr(comma + 1));
        return out;
    }

    // Infix form: split on operator surrounded by spaces. Order matters for the
    // multi-char ops so they don't get partially matched.
    static constexpr std::string_view ops[] = {" << ", " >> ", " * ", " ^ "};
    for (auto const op : ops) {
        auto const pos = lhs_part.find(op);
        if (pos != std::string::npos) {
            out.lhs = trim_ws(std::string_view(lhs_part).substr(0, pos));
            out.rhs = trim_ws(std::string_view(lhs_part).substr(pos + op.size()));
            return out;
        }
    }
    return std::nullopt;
}

// Look up the basis-index mask for a filter name in the algebra. Returns
// nullopt if the name is unknown (caller decides whether to warn).
std::optional<mvec_coeff_filter> get_filter_mask(AlgebraData const& algebra,
                                                 std::string const& filter_name)
{
    if (algebra.dimension == 2) {
        auto const it = algebra.filters_2d.find(filter_name);
        if (it == algebra.filters_2d.end()) return std::nullopt;
        return get_coeff_filter(it->second);
    }
    if (algebra.dimension == 3) {
        auto const it = algebra.filters_3d.find(filter_name);
        if (it == algebra.filters_3d.end()) return std::nullopt;
        return get_coeff_filter(it->second);
    }
    if (algebra.dimension == 4) {
        auto const it = algebra.filters_4d.find(filter_name);
        if (it == algebra.filters_4d.end()) return std::nullopt;
        return get_coeff_filter(it->second);
    }
    return std::nullopt;
}

bool filter_name_known(AlgebraData const& algebra, std::string const& filter_name)
{
    return get_filter_mask(algebra, filter_name).has_value();
}

// Among the algebra's known filters, find the one with the smallest basis
// support that still covers every non-zero component of `prd_mv`. The full
// multivector (`mv`) always covers, so a suggestion is always returned when
// `prd_mv` has any non-zero entry. Ties broken alphabetically for determinism.
std::optional<std::string> suggest_minimal_result_type(AlgebraData const& algebra,
                                                       mvec_coeff const& prd_mv)
{
    auto for_each_filter = [&](auto&& fn) {
        if (algebra.dimension == 2) {
            for (auto const& [name, f] : algebra.filters_2d)
                fn(name, get_coeff_filter(f));
        }
        else if (algebra.dimension == 3) {
            for (auto const& [name, f] : algebra.filters_3d)
                fn(name, get_coeff_filter(f));
        }
        else if (algebra.dimension == 4) {
            for (auto const& [name, f] : algebra.filters_4d)
                fn(name, get_coeff_filter(f));
        }
    };

    std::optional<std::string> best_name;
    std::size_t best_size = std::numeric_limits<std::size_t>::max();

    for_each_filter([&](std::string const& name, mvec_coeff_filter const& mask) {
        if (mask.size() != prd_mv.size()) return;
        for (std::size_t k = 0; k < prd_mv.size(); ++k) {
            if (!is_zero_expr(prd_mv[k]) && mask[k] == 0) return; // doesn't cover
        }
        std::size_t sz = 0;
        for (auto const v : mask) {
            if (v) ++sz;
        }
        if (sz < best_size || (sz == best_size && best_name && name < *best_name)) {
            best_size = sz;
            best_name = name;
        }
    });
    return best_name;
}

// Buffer of all warnings emitted across the run, so the end-of-run summary
// can reproduce them without forcing the user to scroll up. Single-threaded
// generator, so no synchronization needed.
std::vector<std::string>& warning_buffer()
{
    static std::vector<std::string> buf;
    return buf;
}

void warn_case(AlgebraData const& algebra, ProductConfig const& config,
               OutputCase const& case_def, std::string const& msg)
{
    auto const line = fmt::format("WARNING [{} :: {} :: {}]: {}", algebra.name,
                                  config.product_name, case_def.case_name, msg);
    fmt::print(stderr, "{}\n", line);
    warning_buffer().push_back(line);
}

// Run all checks (A, B, C) against a single OutputCase. `prd_mv` may be empty;
// if so, check C is skipped. `prd_mv` (when provided) must align with
// `algebra.basis` index-by-index.
void validate_case(AlgebraData const& algebra, ProductConfig const& config,
                   OutputCase const& case_def, mvec_coeff const& prd_mv)
{
    // Pre-check: filter_name fields must reference filters that actually exist
    // in the algebra. If they don't, the existing generation paths throw with
    // an opaque message — surfacing it here gives better context.
    bool const left_filter_known = filter_name_known(algebra, case_def.left_filter_name);
    bool const right_filter_known =
        filter_name_known(algebra, case_def.right_filter_name);
    if (!left_filter_known) {
        warn_case(algebra, config, case_def,
                  fmt::format("left_filter '{}' is not a known filter for "
                              "algebra '{}'",
                              case_def.left_filter_name, algebra.name));
    }
    if (!right_filter_known) {
        warn_case(algebra, config, case_def,
                  fmt::format("right_filter '{}' is not a known filter for "
                              "algebra '{}'",
                              case_def.right_filter_name, algebra.name));
    }

    auto const tokens = parse_case_tokens(case_def.case_name);
    if (!tokens) {
        warn_case(algebra, config, case_def,
                  "could not parse case_name (expected 'LHS OP RHS -> RESULT' "
                  "or 'func(LHS,RHS) -> RESULT')");
        return;
    }

    // Check D: when case_name uses functional form, the function name should
    // match ProductConfig.product_name (or its first whitespace-delimited word,
    // so suffixes like " (alternative)" still match). Catches misplaced cases
    // like a "wdg(...)" entry inside a gpr ProductConfig.
    auto const lparen = case_def.case_name.find('(');
    if (lparen != std::string::npos) {
        std::size_t a = 0;
        while (a < lparen &&
               std::isspace(static_cast<unsigned char>(case_def.case_name[a]))) {
            ++a;
        }
        std::size_t b = a;
        while (b < lparen &&
               (std::isalnum(static_cast<unsigned char>(case_def.case_name[b])) ||
                case_def.case_name[b] == '_')) {
            ++b;
        }
        if (b == lparen && b > a) {
            std::string const func_name = case_def.case_name.substr(a, b - a);
            std::size_t pa = 0;
            while (pa < config.product_name.size() &&
                   std::isspace(static_cast<unsigned char>(config.product_name[pa]))) {
                ++pa;
            }
            std::size_t pb = pa;
            while (pb < config.product_name.size() &&
                   (std::isalnum(static_cast<unsigned char>(config.product_name[pb])) ||
                    config.product_name[pb] == '_')) {
                ++pb;
            }
            std::string const product_word = config.product_name.substr(pa, pb - pa);
            if (!product_word.empty() && func_name != product_word) {
                warn_case(algebra, config, case_def,
                          fmt::format("case_name function '{}' does not match "
                                      "product '{}' -- expected '{}(...)'",
                                      func_name, config.product_name, product_word));
            }
        }
    }

    // Check A: case_name input tokens match the filter_name fields. When the
    // case_name token is itself a known filter, suggest aligning the filter to
    // it (case_name is descriptive — typos there are rarer than in the filter).
    if (tokens->lhs != case_def.left_filter_name) {
        std::string msg =
            fmt::format("declared LHS type '{}' does not match left_filter '{}'",
                        tokens->lhs, case_def.left_filter_name);
        if (filter_name_known(algebra, tokens->lhs)) {
            msg += fmt::format(" -- suggested left_filter: '{}' (to match "
                               "case_name LHS)",
                               tokens->lhs);
        }
        warn_case(algebra, config, case_def, msg);
    }
    if (tokens->rhs != case_def.right_filter_name) {
        std::string msg =
            fmt::format("declared RHS type '{}' does not match right_filter '{}'",
                        tokens->rhs, case_def.right_filter_name);
        if (filter_name_known(algebra, tokens->rhs)) {
            msg += fmt::format(" -- suggested right_filter: '{}' (to match "
                               "case_name RHS)",
                               tokens->rhs);
        }
        warn_case(algebra, config, case_def, msg);
    }

    // Check F: even/odd coefficient names must pair with the matching grade
    // filter. A coefficient like A_even has zero placeholders at odd basis
    // positions, so pairing it with a non-mv_e filter (e.g. vec) silently
    // produces all zeros at the kept positions. Catches a class of bugs that
    // checks A and E might miss together (operand still parses, result is just
    // numerically zero where the user expected content).
    auto const check_parity = [&](std::string const& side, std::string const& coeff,
                                  std::string const& filter_name) {
        std::string_view const sv = coeff;
        if (sv.ends_with("_even") && filter_name != "mv_e") {
            warn_case(algebra, config, case_def,
                      fmt::format("{}_coeff '{}' has '_even' suffix but "
                                  "{}_filter is '{}' (expected 'mv_e')",
                                  side, coeff, side, filter_name));
        }
        else if (sv.ends_with("_odd") && filter_name != "mv_u") {
            warn_case(algebra, config, case_def,
                      fmt::format("{}_coeff '{}' has '_odd' suffix but "
                                  "{}_filter is '{}' (expected 'mv_u')",
                                  side, coeff, side, filter_name));
        }
    };
    check_parity("left", case_def.left_coeff_name, case_def.left_filter_name);
    check_parity("right", case_def.right_coeff_name, case_def.right_filter_name);

    // Check B: result token is a known filter (or "0" for zero-result cases).
    if (tokens->result == "0") {
        // "-> 0 X" form: X must also be a known filter if present.
        auto const override_t = codegen::parse_explicit_zero_override(case_def.case_name);
        if (override_t && !filter_name_known(algebra, *override_t)) {
            warn_case(algebra, config, case_def,
                      fmt::format("explicit zero-result override type '{}' is "
                                  "not a known filter for algebra '{}'",
                                  *override_t, algebra.name));
        }
        // Check G: typed-zero convention. Bare '-> 0' is the scalar-zero
        // default; '-> 0 ps' is reserved for cases whose natural return type
        // is pseudoscalar. No other typed-zero forms are used in this
        // codebase — '-> 0 s' is redundant (use bare '-> 0'), and forms
        // like '-> 0 vec' / '-> 0 bivec' break the convention.
        if (override_t && *override_t != "ps") {
            std::string const fix =
                (*override_t == "s")
                    ? "use bare '-> 0' instead (scalar zero is the default)"
                    : "only '-> 0 ps' is allowed as a typed zero — use bare "
                      "'-> 0' for scalar zero";
            warn_case(algebra, config, case_def,
                      fmt::format("typed zero '-> 0 {}' violates convention "
                                  "-- {}",
                                  *override_t, fix));
        }
        return; // No check C on zero results.
    }
    auto const result_mask = get_filter_mask(algebra, tokens->result);
    if (!result_mask) {
        warn_case(algebra, config, case_def,
                  fmt::format("declared result type '{}' is not a known filter "
                              "for algebra '{}'",
                              tokens->result, algebra.name));
        return; // Cannot run check C without a mask.
    }

    // Check C: actually-computed non-zero components must fit inside the
    // declared result type's basis support. (We allow the reverse — a sparser
    // result than the declared type — because legitimate cases produce sparse
    // multivectors of a wider declared type.)
    if (prd_mv.empty()) return;
    if (prd_mv.size() != result_mask->size() || prd_mv.size() != algebra.basis.size()) {
        return; // Shape mismatch — let the existing pipeline surface it.
    }
    std::vector<std::string> stray;
    bool all_zero = true;
    for (std::size_t k = 0; k < prd_mv.size(); ++k) {
        if (!is_zero_expr(prd_mv[k])) {
            all_zero = false;
            if ((*result_mask)[k] == 0) {
                stray.push_back(
                    fmt::format("{}='{}'", algebra.basis[k], trim_ws(prd_mv[k])));
            }
        }
    }
    if (!stray.empty()) {
        std::string joined;
        for (std::size_t k = 0; k < stray.size(); ++k) {
            if (k) joined += ", ";
            joined += stray[k];
        }
        std::string msg = fmt::format("declared result type '{}' cannot hold computed "
                                      "components: {}",
                                      tokens->result, joined);
        auto const suggestion = suggest_minimal_result_type(algebra, prd_mv);
        if (suggestion && *suggestion != tokens->result) {
            msg += fmt::format(" -- suggested minimal result type: '{}'", *suggestion);
        }
        warn_case(algebra, config, case_def, msg);
    }
    else if (all_zero) {
        // Check E: declared a non-zero result type, but the actual computation
        // is identically zero across all basis positions. Catches mistakes like
        // declaring "-> s" for a product whose grade arithmetic forbids any
        // non-zero outcome (e.g. rwdg(bivec,vec) in 4D, where 2+1<4 makes it 0).
        // Skipped when stray components fired check C — that's the more
        // informative warning.
        //
        // Convention: bare '-> 0' means a scalar zero (the codegen default);
        // '-> 0 ps' is reserved for cases whose natural type is pseudoscalar
        // and so needs the typed-zero form to keep PScalar in the C++ return
        // type. No other typed zeros are used.
        std::string const suggestion = (tokens->result == "ps")
                                           ? "'-> 0 ps' (typed pseudoscalar zero)"
                                           : "'-> 0' (default scalar zero)";
        warn_case(algebra, config, case_def,
                  fmt::format("declared result type '{}' but the computation "
                              "is identically zero -- use {}",
                              tokens->result, suggestion));
    }
}

// Compute prd_mv for one case using the same path as generate_single_case.
// Returns an empty vector if any prerequisite is missing — validate_case will
// then skip check C, and the existing generation path will surface a clearer
// error when it actually tries to use these inputs.
mvec_coeff compute_prd_mv_for_validation(AlgebraData const& algebra,
                                         OutputCase const& case_def,
                                         prd_table const& basis_tab)
{
    auto const left_it = algebra.coefficients.find(case_def.left_coeff_name);
    auto const right_it = algebra.coefficients.find(case_def.right_coeff_name);
    if (left_it == algebra.coefficients.end() || right_it == algebra.coefficients.end()) {
        return {};
    }
    if (!filter_name_known(algebra, case_def.left_filter_name) ||
        !filter_name_known(algebra, case_def.right_filter_name)) {
        return {};
    }
    auto const prd_tab = get_prd_tab(basis_tab, left_it->second, right_it->second);
    if (algebra.dimension == 2) {
        auto const lf = algebra.filters_2d.at(case_def.left_filter_name);
        auto const rf = algebra.filters_2d.at(case_def.right_filter_name);
        return get_mv_from_prd_tab(prd_tab, algebra.basis, lf, rf);
    }
    if (algebra.dimension == 3) {
        auto const lf = algebra.filters_3d.at(case_def.left_filter_name);
        auto const rf = algebra.filters_3d.at(case_def.right_filter_name);
        return get_mv_from_prd_tab(prd_tab, algebra.basis, lf, rf);
    }
    if (algebra.dimension == 4) {
        auto const lf = algebra.filters_4d.at(case_def.left_filter_name);
        auto const rf = algebra.filters_4d.at(case_def.right_filter_name);
        return get_mv_from_prd_tab(prd_tab, algebra.basis, lf, rf);
    }
    return {};
}

} // namespace

void configurable::print_validation_summary()
{
    auto& buf = warning_buffer();
    if (buf.empty()) return;

    constexpr char sep[] =
        "==========================================================================";
    fmt::print(stderr, "{}\n", sep);
    fmt::print(stderr, "Validation summary: {} warning{}\n", buf.size(),
               buf.size() == 1 ? "" : "s");
    fmt::print(stderr, "{}\n", sep);
    for (auto const& line : buf) {
        fmt::print(stderr, "  {}\n", line);
    }
    fmt::print(stderr, "{}\n", sep);
    buf.clear();
}

////////////////////////////////////////////////////////////////////////////////
// Template function for (regressive) transwedge geometric product calculation
////////////////////////////////////////////////////////////////////////////////

template <typename GetLhsFunc, typename GetRhsFunc>
std::pair<prd_table, std::vector<prd_table>> calculate_transwedge_geometric_product(
    mvec_coeff const& basis, std::vector<mvec_coeff> const& basis_kvec,
    prd_rules const& wedge_rules, GetLhsFunc const& get_lhs, GetRhsFunc const& get_rhs)
{
    std::vector<prd_table> tab; // result tables for each index
    prd_table tab_res = init_zero_product_table(basis.size());

    // calculate the transwedge product and the geometric product resulting
    // from summing up its components
    size_t idx = 0; // index into multivector
    for (size_t k = 0; k < basis_kvec.size(); ++k) {
        // for each order k
        for (size_t cnt = 0; cnt < basis_kvec[k].size(); ++cnt) {
            // for each element within each order k
            auto c = init_zero_multivector(basis.size());
            c[idx] = basis[idx]; // c is a multivector with one entry != "0"

            auto lhs = get_lhs(c, idx); // lhs of transwedge product of order k
            auto rhs = get_rhs(c, idx); // rhs of transwedge product of order k

            tab.emplace_back(apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wedge_rules));

            if ((k * (k - 1) / 2) % 2 != 0) {
                // negate tab for odd values of k*(k-1)/2
                // sign = (-1)^[k*(k-1)/2]
                tab.back() = negate_prd_tab(tab.back());
            }
            tab_res = add_prd_tab(tab_res, tab.back());

            ++idx; // advance to next multivector entry
        }
    } // tab_res contains geometric product table after loops over k & cnt
      // tab[k] contains (-1)^[(k * (k - 1) / 2] * trw[k]

    return {tab_res, tab};
}

template <typename GetLhsFunc, typename GetRhsFunc>
std::pair<prd_table, std::vector<prd_table>>
calculate_regressive_transwedge_geometric_product(
    mvec_coeff const& basis, std::vector<mvec_coeff> const& basis_kvec,
    prd_rules const& wedge_rules, prd_rules const& r_cmpl_rules,
    GetLhsFunc const& get_lhs, GetRhsFunc const& get_rhs)
{
    std::vector<prd_table> tab; // result tables for each index
    prd_table tab_res = init_zero_product_table(basis.size());

    // calculate the transwedge product and the geometric product resulting
    // from summing up its components
    size_t idx = 0; // index into multivector
    for (size_t k = 0; k < basis_kvec.size(); ++k) {
        // for each order k
        for (size_t cnt = 0; cnt < basis_kvec[k].size(); ++cnt) {
            // for each element within each order k
            auto c = init_zero_multivector(basis.size());
            c[idx] = basis[idx]; // c is a multivector with one entry != "0"

            auto lhs = get_lhs(c, idx); // lhs of transwedge product of order k
            auto rhs = get_rhs(c, idx); // rhs of transwedge product of order k

            tab.emplace_back(apply_rules_to_tab(
                apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                   wedge_rules),
                r_cmpl_rules));

            if ((k * (k - 1) / 2) % 2 != 0) {
                // negate tab for odd values of k*(k-1)/2
                // sign = (-1)^[k*(k-1)/2]
                tab.back() = negate_prd_tab(tab.back());
            }
            tab_res = add_prd_tab(tab_res, tab.back());

            ++idx; // advance to next multivector entry
        }
    } // tab_res contains geometric product table after loops over k & cnt
      // tab[k] contains (-1)^[(k * (k - 1) / 2] * rtrw[k]

    return {tab_res, tab};
}

void ConfigurableGenerator::generate_product_expressions(AlgebraData const& algebra,
                                                         ProductConfig const& config,
                                                         GeneratorOptions const& options)
{
    // No separator logic here - separators are handled at the end of each product

    // Skip "alternative" product variants when only code output is requested. They
    // are display-only computations (the canonical product is computed by another
    // path) and `get_basis_table_for_product` emits intermediate basis tables for
    // them as a side-effect, which would pollute --output=code with non-C++ text.
    bool const code_only = options.should_show_code() && !options.should_show_coeffs() &&
                           !options.should_show_tables() &&
                           !options.should_show_metrics();
    if (code_only && config.product_name.find("alternative") != std::string::npos) {
        return;
    }

    // Get the basis table using EXISTING mathematical functions
    auto basis_tab = get_basis_table_for_product(algebra, config.product_name);

    // Validate user-supplied OutputCase tuples and warn (to stderr) on
    // inconsistencies. Sandwich and two-step cases use a different config shape
    // and are skipped here.
    if (!config.is_sandwich_product) {
        for (auto const& case_def : config.cases) {
            if (case_def.is_two_step) continue;
            auto const prd_mv =
                compute_prd_mv_for_validation(algebra, case_def, basis_tab);
            validate_case(algebra, config, case_def, prd_mv);
        }
    }

    // Print basis table if enabled (matches reference format) and if tables are requested
    if (config.show_basis_table && options.should_show_tables()) {
        print_basis_table(algebra, config, basis_tab, options.show_symmetric_asymmetric);
        fmt::println("");
    }

    // Generate each configured case (only if product expressions are requested)
    if (options.should_show_coeffs()) {
        if (config.is_sandwich_product) {
            // For sandwich products, generate the complete two-step process regardless of
            // individual case settings
            OutputCase dummy_case("", "", "", "", "");
            generate_sandwich_case(algebra, config, dummy_case, basis_tab);
        }
        else {
            for (const auto& case_def : config.cases) {
                if (case_def.is_two_step) {
                    generate_sandwich_case(algebra, config, case_def, basis_tab);
                }
                else {
                    generate_single_case(algebra, config, case_def, basis_tab);
                }
            }
        }
    }

    // Emit C++ code (opt-in via --output=code). Sandwich cases are skipped here; their
    // structure differs enough that they're handled by a separate path.
    if (options.should_show_code() && !config.is_sandwich_product) {
        for (const auto& case_def : config.cases) {
            if (case_def.is_two_step) continue;
            emit_single_case_code(algebra, config, case_def, basis_tab);
        }
    }
}

void ConfigurableGenerator::generate_single_case(AlgebraData const& algebra,
                                                 ProductConfig const& config,
                                                 OutputCase const& case_def,
                                                 prd_table const& basis_tab)
{
    // Get coefficients from algebra data (maps to existing coefficient objects)
    auto left_coeff_it = algebra.coefficients.find(case_def.left_coeff_name);
    auto right_coeff_it = algebra.coefficients.find(case_def.right_coeff_name);

    if (left_coeff_it == algebra.coefficients.end()) {
        throw std::invalid_argument("Unknown left coefficient: " +
                                    case_def.left_coeff_name);
    }
    if (right_coeff_it == algebra.coefficients.end()) {
        throw std::invalid_argument("Unknown right coefficient: " +
                                    case_def.right_coeff_name);
    }

    auto const& left_coeff = left_coeff_it->second;
    auto const& right_coeff = right_coeff_it->second;

    // Get filters based on algebra dimension
    mvec_coeff prd_mv;
    auto prd_tab = get_prd_tab(basis_tab, left_coeff, right_coeff);

    if (algebra.dimension == 2) {
        auto left_filter = get_filter_2d(algebra, case_def.left_filter_name);
        auto right_filter = get_filter_2d(algebra, case_def.right_filter_name);
        prd_mv = get_mv_from_prd_tab(prd_tab, algebra.basis, left_filter, right_filter);
    }
    else if (algebra.dimension == 3) {
        auto left_filter = get_filter_3d(algebra, case_def.left_filter_name);
        auto right_filter = get_filter_3d(algebra, case_def.right_filter_name);
        prd_mv = get_mv_from_prd_tab(prd_tab, algebra.basis, left_filter, right_filter);
    }
    else if (algebra.dimension == 4) {
        auto left_filter = get_filter_4d(algebra, case_def.left_filter_name);
        auto right_filter = get_filter_4d(algebra, case_def.right_filter_name);
        prd_mv = get_mv_from_prd_tab(prd_tab, algebra.basis, left_filter, right_filter);
    }
    else {
        throw std::invalid_argument("Unsupported algebra dimension: " +
                                    std::to_string(algebra.dimension));
    }

    // Format output to match reference implementation exactly
    print_case_header(algebra, config, case_def.case_name);
    print_case_result(prd_mv, algebra.basis);
    fmt::println("");
}

void ConfigurableGenerator::emit_single_case_code(AlgebraData const& algebra,
                                                  ProductConfig const& config,
                                                  OutputCase const& case_def,
                                                  prd_table const& basis_tab)
{
    // Lazy per-algebra registry. Codegen is supported for
    // ega2d/ega3d/pga2dp/pga3dp/sta4d; STA4D `.cases` arrays are still empty in
    // algebras/ga_prdxpr_sta4d_config.cpp, so running --output=code --algebra=sta4d will
    // emit nothing until those are populated.
    static std::map<std::string, codegen::TypeRegistry> registries;
    auto reg_it = registries.find(algebra.name);
    if (reg_it == registries.end()) {
        try {
            reg_it = registries.emplace(algebra.name, codegen::TypeRegistry(algebra.name))
                         .first;
        }
        catch (std::exception const& e) {
            fmt::println("// SKIP {} {} :: {} -- {}", algebra.name, config.product_name,
                         case_def.case_name, e.what());
            return;
        }
    }
    auto const& registry = reg_it->second;

    // Compute the per-basis result expressions (same path as generate_single_case).
    auto left_coeff_it = algebra.coefficients.find(case_def.left_coeff_name);
    auto right_coeff_it = algebra.coefficients.find(case_def.right_coeff_name);
    if (left_coeff_it == algebra.coefficients.end() ||
        right_coeff_it == algebra.coefficients.end()) {
        fmt::println("// SKIP {} {} :: {} -- unknown coefficient", algebra.name,
                     config.product_name, case_def.case_name);
        return;
    }

    mvec_coeff prd_mv;
    auto prd_tab = get_prd_tab(basis_tab, left_coeff_it->second, right_coeff_it->second);
    if (algebra.dimension == 2) {
        auto lf = get_filter_2d(algebra, case_def.left_filter_name);
        auto rf = get_filter_2d(algebra, case_def.right_filter_name);
        prd_mv = get_mv_from_prd_tab(prd_tab, algebra.basis, lf, rf);
    }
    else if (algebra.dimension == 3) {
        auto lf = get_filter_3d(algebra, case_def.left_filter_name);
        auto rf = get_filter_3d(algebra, case_def.right_filter_name);
        prd_mv = get_mv_from_prd_tab(prd_tab, algebra.basis, lf, rf);
    }
    else if (algebra.dimension == 4) {
        auto lf = get_filter_4d(algebra, case_def.left_filter_name);
        auto rf = get_filter_4d(algebra, case_def.right_filter_name);
        prd_mv = get_mv_from_prd_tab(prd_tab, algebra.basis, lf, rf);
    }
    else {
        fmt::println("// SKIP {} {} :: {} -- unsupported dimension {}", algebra.name,
                     config.product_name, case_def.case_name, algebra.dimension);
        return;
    }

    std::string skip_reason;
    auto rendered =
        codegen::emit_function(algebra, config, case_def, prd_mv, registry, &skip_reason);
    if (!rendered) {
        fmt::println("// SKIP {} {} :: {} -- {}", algebra.name, config.product_name,
                     case_def.case_name, skip_reason);
        return;
    }
    fmt::println("// {} {} :: {}", algebra.name, config.product_name, case_def.case_name);
    fmt::print("{}", *rendered);
    fmt::println("");
}

void ConfigurableGenerator::generate_sandwich_case(AlgebraData const& algebra,
                                                   ProductConfig const& config,
                                                   OutputCase const& /* case_def */,
                                                   prd_table const& basis_tab)
{
    std::string prd_name = algebra.name + " " + config.display_name;

    if (algebra.name == "ega2d") {

        // First product between multivectors in basis_tab (R * v)
        fmt::println("{}:", prd_name + space_str() + "mv_e * vec -> vec_tmp");
        auto prd_tab = get_prd_tab(basis_tab, mv2d_coeff_R_even, mv2d_coeff_svps);

        auto vec_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_2d::mv_e,
                                           filter_2d::vec, brace_switch::use_braces);
        fmt::println("vec_tmp:");
        print_mvec(vec_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product v * rev(R)
        fmt::println("{}:", prd_name + space_str() + "vec_tmp * rev(mv_e) -> vec_res");
        prd_tab = get_prd_tab(basis_tab, vec_tmp, mv2d_coeff_R_rev_even);

        auto vec_res =
            get_mv_from_prd_tab(prd_tab, algebra.basis, filter_2d::vec, filter_2d::mv_e);
        print_mvec(vec_res, algebra.basis);
        fmt::println("");

        // Add transformation output using consistent format
        print_transformed_result(vec_res, algebra.basis, algebra, config);
    }

    else if (algebra.name == "ega3d") {

        // First product between multivectors in basis_tab (R * v)
        fmt::println("{}:", prd_name + space_str() + "mv_e * vec -> mv_u_tmp");
        auto prd_tab = get_prd_tab(basis_tab, mv3d_coeff_R_even, mv3d_coeff_svBps);

        auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_3d::mv_e,
                                            filter_3d::vec, brace_switch::use_braces);
        fmt::println("mv_u_tmp:");
        print_mvec(mv_u_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product mv_u_tmp * rev(R)
        fmt::println("{}:", prd_name + space_str() + "mv_u_tmp * rev(mv_e) -> mv_u_res");
        auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv3d_coeff_R_rev_even);

        auto mv_u_res_v = get_mv_from_prd_tab(prd_tab_v, algebra.basis, filter_3d::mv_u,
                                              filter_3d::mv_e);
        print_mvec(mv_u_res_v, algebra.basis);
        fmt::println("");

        // Add transformation output for vector case
        print_transformed_result(mv_u_res_v, algebra.basis, algebra, config);

        //// Bivector case

        // First product between multivectors in basis_tab (R * B)
        fmt::println("{}:", prd_name + space_str() + "mv_e * bivec -> mv_e_tmp");
        auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_3d::mv_e,
                                            filter_3d::bivec, brace_switch::use_braces);
        fmt::println("mv_e_tmp:");
        print_mvec(mv_e_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product mv_e_tmp * rev(R)
        fmt::println("{}:", prd_name + space_str() + "mv_e_tmp * rev(mv_e) -> mv_e_res");
        auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv3d_coeff_R_rev_even);

        auto mv_e_res_B = get_mv_from_prd_tab(prd_tab_B, algebra.basis, filter_3d::mv_e,
                                              filter_3d::mv_e);
        print_mvec(mv_e_res_B, algebra.basis);
        fmt::println("");

        // Add transformation output for bivector case
        print_transformed_result(mv_e_res_B, algebra.basis, algebra, config);
    }

    else if (algebra.name == "pga2dp") {

        // Vector case: first product between multivectors in basis_tab rgpr(M,v)
        fmt::println("{}:", prd_name + space_str() + "rgpr(mv_u, vec) -> mv_u_tmp");
        auto prd_tab = get_prd_tab(basis_tab, mv2dp_coeff_R_odd, mv2dp_coeff_svBps);

        auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_3d::mv_u,
                                            filter_3d::vec, brace_switch::use_braces);
        fmt::println("mv_u_tmp:");
        print_mvec(mv_u_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product rgpr(mv_u_tmp, rrev(M))
        fmt::println("{}:",
                     prd_name + space_str() + "rgpr(mv_u_tmp, rrev(mv_u)) -> mv_u_res");
        auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv2dp_coeff_R_rrev_odd);

        auto mv_u_res_v = get_mv_from_prd_tab(prd_tab_v, algebra.basis, filter_3d::mv_u,
                                              filter_3d::mv_u);
        print_mvec(mv_u_res_v, algebra.basis);
        fmt::println("");

        // Add transformation output for vector case
        print_transformed_result(mv_u_res_v, algebra.basis, algebra, config);

        //// Bivector case

        // First product between multivectors in basis_tab rgpr(M,B)
        fmt::println("{}:", prd_name + space_str() + "rgpr(mv_u, bivec) -> mv_e_tmp");
        auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_3d::mv_u,
                                            filter_3d::bivec, brace_switch::use_braces);
        fmt::println("mv_e_tmp:");
        print_mvec(mv_e_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product rgpr(mv_e_tmp,rrev(M))
        fmt::println("{}:",
                     prd_name + space_str() + "rgpr(mv_e_tmp, rrev(mv_u)) -> mv_e_res");
        auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv2dp_coeff_R_rrev_odd);

        auto mv_e_res_B = get_mv_from_prd_tab(prd_tab_B, algebra.basis, filter_3d::mv_e,
                                              filter_3d::mv_u);
        print_mvec(mv_e_res_B, algebra.basis);
        fmt::println("");

        // Add transformation output for bivector case
        print_transformed_result(mv_e_res_B, algebra.basis, algebra, config);
    }

    else if (algebra.name == "pga3dp") {

        // Vector case: first product between multivectors in basis_tab rgpr(M,v)
        fmt::println("{}:", prd_name + space_str() + "rgpr(mv_e, vec) -> mv_u_tmp");
        auto prd_tab = get_prd_tab(basis_tab, mv3dp_coeff_R_even, mv3dp_coeff_svBtps);

        auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_4d::mv_e,
                                            filter_4d::vec, brace_switch::use_braces);
        fmt::println("mv_u_tmp:");
        print_mvec(mv_u_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product rgpr(mv_u_tmp,rrev(M))
        fmt::println("{}:",
                     prd_name + space_str() + "rgpr(mv_u_tmp, rrev(mv_e)) -> mv_u_res");
        auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mv3dp_coeff_R_rrev_even);

        auto mv_u_res_v = get_mv_from_prd_tab(prd_tab_v, algebra.basis, filter_4d::mv_u,
                                              filter_4d::mv_e);
        print_mvec(mv_u_res_v, algebra.basis);
        fmt::println("");

        // Add transformation output for vector case
        print_transformed_result(mv_u_res_v, algebra.basis, algebra, config);

        //// Bivector case

        // First product between multivectors in basis_tab rgpr(M,B)
        fmt::println("{}:", prd_name + space_str() + "rgpr(mv_e, bivec) -> mv_e_tmp");
        auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_4d::mv_e,
                                            filter_4d::bivec, brace_switch::use_braces);
        fmt::println("mv_e_tmp:");
        print_mvec(mv_e_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product rgpr(mv_e_tmp, rrev(M))
        fmt::println("{}:",
                     prd_name + space_str() + "rgpr(mv_e_tmp, rrev(mv_e)) -> mv_e_res");
        auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mv3dp_coeff_R_rrev_even);

        auto mv_e_res_B = get_mv_from_prd_tab(prd_tab_B, algebra.basis, filter_4d::mv_e,
                                              filter_4d::mv_e);
        print_mvec(mv_e_res_B, algebra.basis);
        fmt::println("");

        // Add transformation output for bivector case
        print_transformed_result(mv_e_res_B, algebra.basis, algebra, config);

        //// Trivector case

        // First product between multivectors in basis_tab rgpr(M, T)
        fmt::println("{}:", prd_name + space_str() + "rgpr(mv_e, trivec) -> mv_u_tmp_t");
        auto mv_u_tmp_t =
            get_mv_from_prd_tab(prd_tab, algebra.basis, filter_4d::mv_e,
                                filter_4d::trivec, brace_switch::use_braces);
        fmt::println("mv_u_tmp_t:");
        print_mvec(mv_u_tmp_t, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product rgpr(mv_u_tmp_t,rrev(M))
        fmt::println("{}:", prd_name + space_str() +
                                "rgpr(mv_u_tmp_t, rrev(mv_e)) -> mv_u_res_t");
        auto prd_tab_t = get_prd_tab(basis_tab, mv_u_tmp_t, mv3dp_coeff_R_rrev_even);

        auto mv_u_res_t = get_mv_from_prd_tab(prd_tab_t, algebra.basis, filter_4d::mv_u,
                                              filter_4d::mv_e);
        print_mvec(mv_u_res_t, algebra.basis);
        fmt::println("");

        // Add transformation output for trivector case
        print_transformed_result(mv_u_res_t, algebra.basis, algebra, config);
    }

    else if (algebra.name == "sta4d") {

        // Vector case: first product between multivectors in basis_tab gpr(M,v)
        fmt::println("{}:", prd_name + space_str() + "gpr(mv_e, vec) -> mv_u_tmp");
        auto prd_tab = get_prd_tab(basis_tab, mvsta4d_coeff_R_even, mvsta4d_coeff_svBtps);

        auto mv_u_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_4d::mv_e,
                                            filter_4d::vec, brace_switch::use_braces);
        fmt::println("mv_u_tmp:");
        print_mvec(mv_u_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product gpr(mv_u_tmp,rev(M))
        fmt::println("{}:",
                     prd_name + space_str() + "gpr(mv_u_tmp, rev(mv_e)) -> mv_u_res");
        auto prd_tab_v = get_prd_tab(basis_tab, mv_u_tmp, mvsta4d_coeff_R_rev_even);

        auto mv_u_res_v = get_mv_from_prd_tab(prd_tab_v, algebra.basis, filter_4d::mv_u,
                                              filter_4d::mv_e);
        print_mvec(mv_u_res_v, algebra.basis);
        fmt::println("");

        // Add transformation output for vector case
        print_transformed_result(mv_u_res_v, algebra.basis, algebra, config);

        //// Bivector case

        // First product between multivectors in basis_tab rgpr(M,B)
        fmt::println("{}:", prd_name + space_str() + "gpr(mv_e, bivec) -> mv_e_tmp");
        auto mv_e_tmp = get_mv_from_prd_tab(prd_tab, algebra.basis, filter_4d::mv_e,
                                            filter_4d::bivec, brace_switch::use_braces);
        fmt::println("mv_e_tmp:");
        print_mvec(mv_e_tmp, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product gpr(mv_e_tmp, rrev(M))
        fmt::println("{}:",
                     prd_name + space_str() + "gpr(mv_e_tmp, rev(mv_e)) -> mv_e_res");
        auto prd_tab_B = get_prd_tab(basis_tab, mv_e_tmp, mvsta4d_coeff_R_rev_even);

        auto mv_e_res_B = get_mv_from_prd_tab(prd_tab_B, algebra.basis, filter_4d::mv_e,
                                              filter_4d::mv_e);
        print_mvec(mv_e_res_B, algebra.basis);
        fmt::println("");

        // Add transformation output for bivector case
        print_transformed_result(mv_e_res_B, algebra.basis, algebra, config);

        //// Trivector case

        // First product between multivectors in basis_tab gpr(M, T)
        fmt::println("{}:", prd_name + space_str() + "gpr(mv_e, trivec) -> mv_u_tmp_t");
        auto mv_u_tmp_t =
            get_mv_from_prd_tab(prd_tab, algebra.basis, filter_4d::mv_e,
                                filter_4d::trivec, brace_switch::use_braces);
        fmt::println("mv_u_tmp_t:");
        print_mvec(mv_u_tmp_t, algebra.basis);
        fmt::println("");

        // Second product between multivectors for the product gpr(mv_u_tmp_t,rev(M))
        fmt::println("{}:", prd_name + space_str() +
                                "rgpr(mv_u_tmp_t, rev(mv_e)) -> mv_u_res_t");
        auto prd_tab_t = get_prd_tab(basis_tab, mv_u_tmp_t, mvsta4d_coeff_R_rev_even);

        auto mv_u_res_t = get_mv_from_prd_tab(prd_tab_t, algebra.basis, filter_4d::mv_u,
                                              filter_4d::mv_e);
        print_mvec(mv_u_res_t, algebra.basis);
        fmt::println("");

        // Add transformation output for trivector case
        print_transformed_result(mv_u_res_t, algebra.basis, algebra, config);
    }
}

prd_table
ConfigurableGenerator::get_basis_table_for_product(const AlgebraData& algebra,
                                                   const std::string& product_name)
{
    // This function maps product names to the EXISTING mathematical basis table
    // generation from the reference implementation - no changes to mathematical logic

    if (algebra.name == "ega2d") {

        if (product_name == "gpr") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str()),
                gpr_ega2d_rules);
        }

        else if (product_name == "gpr (alternative)") {

            // calculation of the geometric product based on the "transwedge product"
            // as proposed by Lengyel (https://terathon.com/blog/transwedge-product.html)

            // advantage: geometric product is derived from primitives of Grassmann
            //            algebra exclusively (only requires wdg, rwdg, cmpl and dual).
            //            Shows that geometric product is not more fundamental, since
            //            it can be derived from other primitives.

            // Lambda for calculating left-hand side: rwdg(l_cmpl(c),a) =
            // l_cmpl(wdg(a,r_cmpl(a)))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_rcmpl =
                    coeff; // r_cmpl(l_cmpl(coeff)) is identity transformation
                auto rhs_rcmpl = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_rcmpl, wdg_str()),
                    wdg_ega2d_rules);
                auto lhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side: rwdg(b, r_dual(c)) =
            // l_cmpl(wdg(r_cmpl(b),r_cmpl(r_dual(c)))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_rcmpl = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
                auto rhs_rcmpl = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, r_dual_ega2d_rules), r_cmpl_ega2d_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_rcmpl, wdg_str()),
                    wdg_ega2d_rules);
                auto rhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_transwedge_geometric_product(
                mv2d_basis, mv2d_basis_kvec, wdg_ega2d_rules, get_lhs, get_rhs);

            // now print the resulting product table for each order k
            fmt::println("ega2d geometric product (alternative definition) - "
                         "intermediate results:");
            fmt::println("");
            print_product_tables_by_grade(tab, mv2d_basis_kvec);

            return tab_res;
        }

        else if (product_name == "twdg1") {

            // calculation of the geometric product based on the "transwedge product"
            // as proposed by Lengyel (https://terathon.com/blog/transwedge-product.html)

            // advantage: geometric product is derived from primitives of Grassmann
            //            algebra exclusively (only requires wdg, rwdg, cmpl and dual).
            //            Shows that geometric product is not more fundamental, since
            //            it can be derived from other primitives.

            // Lambda for calculating left-hand side: rwdg(l_cmpl(c),a) =
            // l_cmpl(wdg(a,r_cmpl(a)))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_rcmpl =
                    coeff; // r_cmpl(l_cmpl(coeff)) is identity transformation
                auto rhs_rcmpl = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_rcmpl, wdg_str()),
                    wdg_ega2d_rules);
                auto lhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side: rwdg(b, r_dual(c)) =
            // l_cmpl(wdg(r_cmpl(b),r_cmpl(r_dual(c)))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_rcmpl = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
                auto rhs_rcmpl = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, r_dual_ega2d_rules), r_cmpl_ega2d_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_rcmpl, wdg_str()),
                    wdg_ega2d_rules);
                auto rhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_transwedge_geometric_product(
                mv2d_basis, mv2d_basis_kvec, wdg_ega2d_rules, get_lhs, get_rhs);

            tab_res = add_prd_tab(tab[1], tab[2]);
            return tab_res; // return transwedge product for k=1
        }

        else if (product_name == "cmt") {
            // Commutator product (=asymmetric part of the geometric product)
            auto basis_tab = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str()),
                gpr_ega2d_rules);
            return get_prd_tab_asym(basis_tab); // use the asymmetric part only
        }

        else if (product_name == "wdg") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, wdg_str()),
                wdg_ega2d_rules);
        }

        else if (product_name == "dot") {
            // Inner product (=dot product)
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str()),
                dot_ega2d_rules);
        }

        else if (product_name == "dot (alternative)") {
            // Inner product (=dot product) - defined by outer product
            // (see Grassmann Algebra, John Browne, p. 35, p. 373)
            // remember: A ^ r_cmpl(A) = I_2d
            // inner product = rwdg(A,r_cmpl(A))
            //               = l_cmpl( wdg( r_cmpl(A), r_cmpl( r_cmpl(A) ) ) )
            // => only valid, if both operands are of the same grade, i.e. only the
            // main diagonal is valid

            auto cmpl_func = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
            auto dbl_cmpl_func = apply_rules_to_mv(cmpl_func, r_cmpl_ega2d_rules);
            // fmt::println("");
            // fmt::println("mv2d_basis = {}", mv2d_basis);
            // fmt::println("cmpl_func = {}", cmpl_func);
            // fmt::println("dbl_cmpl_func = {}", dbl_cmpl_func);
            // fmt::println("");
            auto tmp_tab = mv_coeff_to_coeff_prd_tab(cmpl_func, dbl_cmpl_func, wdg_str());
            // fmt::println("tmp_tab = {}", tmp_tab);
            // fmt::println("");
            auto basis_tab_with_rules = apply_rules_to_tab(tmp_tab, wdg_ega2d_rules);
            // fmt::println("basis_tab_with_rules = {}", basis_tab_with_rules);
            // fmt::println("");
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
        }

        else if (product_name == "l_contract") {
            // A << B = rwdg(l_cmpl(A), B)
            //        = l_cmpl( wdg( r_cmpl(l_cmpl(A)), r_cmpl(B) ) )
            //        = l_cmpl( wdg(               A,   r_cmpl(B) ) )
            //
            // auto lhs = apply_rules_to_mv(apply_rules_to_mv(mv2d_basis,
            // l_cmpl_ega2d_rules), r_cmpl_ega2d_rules);
            auto lhs = mv2d_basis;
            auto rhs = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_ega2d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
        }

        else if (product_name == "r_contract") {
            // A >> B = rwdg(A,r_cmpl(B))
            //        = l_cmpl( wdg( r_cmpl(A),r_cmpl(r_cmpl(B)) )  )
            auto lhs = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
            auto rhs = apply_rules_to_mv(
                apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules), r_cmpl_ega2d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_ega2d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
        }

        else if (product_name == "l_expand") {
            // Left expansion: l_expand(A,B) = wdg(l_cmpl(A), B)
            auto lhs = apply_rules_to_mv(mv2d_basis, l_cmpl_ega2d_rules);
            auto rhs = mv2d_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_ega2d_rules);
        }

        else if (product_name == "r_expand") {
            // Right expansion: r_expand(A,B) = wdg(A, r_cmpl(B))
            auto lhs = mv2d_basis;
            auto rhs = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_ega2d_rules);
        }

        else if (product_name == "rgpr") {
            // Regressive geometric: rgpr(A,B) = l_cmpl(gpr(r_cmpl(A), r_cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_ega2d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
        }

        else if (product_name == "rgpr (alternative)") {

            // calculation of the regressive geometric product based on the "anti
            // transwedge product" as proposed by Lengyel
            // (https://terathon.com/blog/transwedge-product.html)

            // advantage: regressive geometric product is derived from primitives of
            //            Grassmann algebra exclusively (only requires wdg, rwdg, cmpl and
            //            dual). Shows that geometric product is not more fundamental,
            //            since it can be derived from other primitives.

            // Lambda for calculating left-hand side: l_cmpl(wdg(c,a))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = coeff;
                auto rhs_arg = mv2d_basis;
                auto lhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_ega2d_rules),
                    l_cmpl_ega2d_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side:
            // l_cmpl(wdg(b,r_cmpl(r_dual(c))))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = mv2d_basis;
                auto rhs_arg = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, r_dual_ega2d_rules), r_cmpl_ega2d_rules);
                auto rhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_ega2d_rules),
                    l_cmpl_ega2d_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_regressive_transwedge_geometric_product(
                mv2d_basis, mv2d_basis_kvec, wdg_ega2d_rules, r_cmpl_ega2d_rules, get_lhs,
                get_rhs);

            // now print the resulting product table for each order k
            fmt::println("ega2d regressive geometric product (alternative definition) - "
                         "intermediate results:");
            fmt::println("");
            print_product_tables_by_grade(tab, mv2d_basis_kvec);

            return tab_res;
        }

        else if (product_name == "rtwdg1") {

            // calculation of the regressive geometric product based on the "anti
            // transwedge product" as proposed by Lengyel
            // (https://terathon.com/blog/transwedge-product.html)

            // advantage: regressive geometric product is derived from primitives of
            //            Grassmann algebra exclusively (only requires wdg, rwdg, cmpl and
            //            dual). Shows that geometric product is not more fundamental,
            //            since it can be derived from other primitives.

            // Lambda for calculating left-hand side: l_cmpl(wdg(c,a))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = coeff;
                auto rhs_arg = mv2d_basis;
                auto lhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_ega2d_rules),
                    l_cmpl_ega2d_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side:
            // l_cmpl(wdg(b,r_cmpl(r_dual(c))))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = mv2d_basis;
                auto rhs_arg = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, r_dual_ega2d_rules), r_cmpl_ega2d_rules);
                auto rhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_ega2d_rules),
                    l_cmpl_ega2d_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_regressive_transwedge_geometric_product(
                mv2d_basis, mv2d_basis_kvec, wdg_ega2d_rules, r_cmpl_ega2d_rules, get_lhs,
                get_rhs);

            tab_res = add_prd_tab(tab[1], tab[2]);
            return tab_res; // return transwedge product for k=1
        }

        else if (product_name == "rcmt") {
            // Commutator product: cmt(A,B) = asym(gpr(A,B))
            // Regressive commutator product:
            //                    rcmt(A,B) = asym(rgpr(A,B))
            //                    rcmt(A,B) = asym(l_cmpl(gpr(r_cmpl(A),r_cmpl(B))))
            auto basis_cmpl_func = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_ega2d_rules);
            auto full_tab = apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
            return get_prd_tab_asym(full_tab);
        }

        else if (product_name == "rwdg") {
            // Regressive wedge: rwdg(A,B) = l_cmpl(wdg(r_cmpl(A), r_cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str()),
                wdg_ega2d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
        }

        else if (product_name == "rdot") {
            // Regressive inner: rdot(A,B) = l_cmpl(dot(r_cmpl(A), r_cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2d_basis, r_cmpl_ega2d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                dot_ega2d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_ega2d_rules);
        }

        else if (product_name == "sandwich_gpr") {
            // Rotor sandwich product basis table
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2d_basis, mv2d_basis, mul_str()),
                gpr_ega2d_rules);
        }
    }

    else if (algebra.name == "ega3d") {

        if (product_name == "gpr") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str()),
                gpr_ega3d_rules);
        }

        else if (product_name == "gpr (alternative)") {

            // calculation of the geometric product based on the "transwedge product"
            // as proposed by Lengyel (https://terathon.com/blog/transwedge-product.html)

            // advantage: geometric product is derived from primitives of Grassmann
            //            algebra exclusively (only requires wdg, rwdg, cmpl and dual).
            //            Shows that geometric product is not more fundamental, since
            //            it can be derived from other primitives.

            // Lambda for calculating left-hand side: rwdg(cmpl(c),a) =
            // cmpl(wdg(a,cmpl(a)))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_cmpl = coeff; // cmpl(cmpl(coeff)) is identity transformation
                auto rhs_cmpl = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_cmpl, rhs_cmpl, wdg_str()),
                    wdg_ega3d_rules);
                auto lhs_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side: rwdg(b, dual(c)) =
            // cmpl(wdg(cmpl(b),cmpl(dual(c)))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_cmpl = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
                auto rhs_cmpl = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, dual_ega3d_rules), cmpl_ega3d_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_cmpl, rhs_cmpl, wdg_str()),
                    wdg_ega3d_rules);
                auto rhs_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_transwedge_geometric_product(
                mv3d_basis, mv3d_basis_kvec, wdg_ega3d_rules, get_lhs, get_rhs);

            // now print the resulting product table for each order k
            fmt::println("ega3d geometric product (alternative definition) - "
                         "intermediate results:");
            fmt::println("");
            print_product_tables_by_grade(tab, mv3d_basis_kvec);

            return tab_res; // return full geometric product
        }

        else if (product_name == "twdg1") {

            // calculation of the geometric product based on the "transwedge product"
            // as proposed by Lengyel (https://terathon.com/blog/transwedge-product.html)

            // advantage: geometric product is derived from primitives of Grassmann
            //            algebra exclusively (only requires wdg, rwdg, cmpl and dual).
            //            Shows that geometric product is not more fundamental, since
            //            it can be derived from other primitives.

            // Lambda for calculating left-hand side: rwdg(cmpl(c),a) =
            // cmpl(wdg(a,cmpl(a)))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_cmpl = coeff; // cmpl(cmpl(coeff)) is identity transformation
                auto rhs_cmpl = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_cmpl, rhs_cmpl, wdg_str()),
                    wdg_ega3d_rules);
                auto lhs_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side: rwdg(b, dual(c)) =
            // cmpl(wdg(cmpl(b),cmpl(dual(c)))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_cmpl = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
                auto rhs_cmpl = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, dual_ega3d_rules), cmpl_ega3d_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_cmpl, rhs_cmpl, wdg_str()),
                    wdg_ega3d_rules);
                auto rhs_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_transwedge_geometric_product(
                mv3d_basis, mv3d_basis_kvec, wdg_ega3d_rules, get_lhs, get_rhs);

            tab_res = add_prd_tab(tab[1], tab[2]);
            tab_res = add_prd_tab(tab_res, tab[3]);
            return tab_res; // return transwedge product for k=1
        }

        else if (product_name == "cmt") {
            // Commutator product (=asymmetric part of the geometric product)
            auto basis_tab = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str()),
                gpr_ega3d_rules);
            return get_prd_tab_asym(basis_tab); // use the asymmetric part only
        }

        else if (product_name == "wdg") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, wdg_str()),
                wdg_ega3d_rules);
        }

        else if (product_name == "dot") {
            // Inner product (=dot product)
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str()),
                dot_ega3d_rules);
        }

        else if (product_name == "l_contract") {
            // Left contraction: A << B = rwdg(l_cmpl(A), B)
            // For EGA3D: lcmpl = cmpl
            //                          = cmpl( wdg( cmpl(cmpl(A)), cmpl(B) ) )
            auto lhs = apply_rules_to_mv(apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules),
                                         cmpl_ega3d_rules);
            auto rhs = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_ega3d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
        }

        else if (product_name == "r_contract") {
            // Right contraction: A >> B = rwdg(A, r_cmpl(B)))
            // For EGA3D: lcmpl = rcmpl = cmpl
            //                           = cmpl( wdg( cmpl(A),cmpl(cmpl(B)) )  )
            auto lhs = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto rhs = apply_rules_to_mv(apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules),
                                         cmpl_ega3d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_ega3d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
        }

        else if (product_name == "l_expand") {
            // Left expansion: l_expand(A,B) = wdg(l_cmpl(A), B)
            // For EGA3D: lcmpl = cmpl
            auto lhs = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto rhs = mv3d_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_ega3d_rules);
        }

        else if (product_name == "r_expand") {
            // Right expansion: r_expand(A,B) = wdg(A, r_cmpl(B))
            // For EGA3D: rcmpl = cmpl
            auto lhs = mv3d_basis;
            auto rhs = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_ega3d_rules);
        }

        else if (product_name == "rgpr") {
            // Regressive geometric: rgpr(A,B) = l_cmpl(gpr(r_cmpl(A), r_cmpl(B)))
            //                                 =  cmpl(gpr( cmpl(A),  cmpl(B))) in ega3d
            auto basis_cmpl_func = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_ega3d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
        }

        else if (product_name == "rgpr (alternative)") {

            // calculation of the regressive geometric product based on the "anti
            // transwedge product" as proposed by Lengyel
            // (https://terathon.com/blog/transwedge-product.html)

            // advantage: regressive geometric product is derived from primitives of
            //            Grassmann algebra exclusively (only requires wdg, rwdg, cmpl and
            //            dual). Shows that geometric product is not more fundamental,
            //            since it can be derived from other primitives.

            // Lambda for calculating left-hand side: cmpl(wdg(c,a))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = coeff;
                auto rhs_arg = mv3d_basis;
                auto lhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_ega3d_rules),
                    cmpl_ega3d_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side:
            // cmpl(wdg(b,cmpl(dual(c))))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = mv3d_basis;
                auto rhs_arg = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, dual_ega3d_rules), cmpl_ega3d_rules);
                auto rhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_ega3d_rules),
                    cmpl_ega3d_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_regressive_transwedge_geometric_product(
                mv3d_basis, mv3d_basis_kvec, wdg_ega3d_rules, cmpl_ega3d_rules, get_lhs,
                get_rhs);

            // now print the resulting product table for each order k
            fmt::println("ega3d regressive geometric product (alternative definition) - "
                         "intermediate results:");
            fmt::println("");
            print_product_tables_by_grade(tab, mv3d_basis_kvec);

            return tab_res;
        }

        else if (product_name == "rtwdg1") {

            // calculation of the regressive geometric product based on the "anti
            // transwedge product" as proposed by Lengyel
            // (https://terathon.com/blog/transwedge-product.html)

            // advantage: regressive geometric product is derived from primitives of
            //            Grassmann algebra exclusively (only requires wdg, rwdg, cmpl and
            //            dual). Shows that geometric product is not more fundamental,
            //            since it can be derived from other primitives.

            // Lambda for calculating left-hand side: cmpl(wdg(c,a))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = coeff;
                auto rhs_arg = mv3d_basis;
                auto lhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_ega3d_rules),
                    cmpl_ega3d_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side:
            // cmpl(wdg(b,cmpl(dual(c))))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = mv3d_basis;
                auto rhs_arg = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, dual_ega3d_rules), cmpl_ega3d_rules);
                auto rhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_ega3d_rules),
                    cmpl_ega3d_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_regressive_transwedge_geometric_product(
                mv3d_basis, mv3d_basis_kvec, wdg_ega3d_rules, cmpl_ega3d_rules, get_lhs,
                get_rhs);

            tab_res = add_prd_tab(tab[1], tab[2]);
            tab_res = add_prd_tab(tab_res, tab[3]);
            return tab_res; // return regressive transwedge product for k=1
        }

        else if (product_name == "rcmt") {
            // Commutator product: cmt(A,B) = asym(gpr(A,B))
            // Regressive commutator product:
            //                    rcmt(A,B) = asym(rgpr(A,B))
            //                    rcmt(A,B) = asym(cmpl(gpr(cmpl(A),cmpl(B))))
            auto basis_cmpl_func = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_ega3d_rules);
            auto full_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
            return get_prd_tab_asym(full_tab);
        }

        else if (product_name == "rwdg") {
            // Regressive wedge: rwdg(A,B) = l_cmpl(wdg(r_cmpl(A), r_cmpl(B)))
            //                             =   cmpl(wdg(  cmpl(A),   cmpl(B))) in ega3d
            auto basis_cmpl_func = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str()),
                wdg_ega3d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
        }

        else if (product_name == "rdot") {
            // Regressive inner: rdot(A,B) = l_cmpl(dot(r_cmpl(A), r_cmpl(B)))
            //                             =   cmpl(dot(  cmpl(A),   cmpl(B))) in ega3d
            auto basis_cmpl_func = apply_rules_to_mv(mv3d_basis, cmpl_ega3d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                dot_ega3d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_ega3d_rules);
        }

        else if (product_name == "sandwich_gpr") {
            // Rotor sandwich product basis table
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3d_basis, mv3d_basis, mul_str()),
                gpr_ega3d_rules);
        }
    }

    else if (algebra.name == "pga2dp") {

        if (product_name == "gpr") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str()),
                gpr_pga2dp_rules);
        }

        else if (product_name == "gpr (alternative)") {

            // calculation of the geometric product based on the "transwedge product"
            // as proposed by Lengyel (https://terathon.com/blog/transwedge-product.html)

            // advantage: geometric product is derived from primitives of Grassmann
            //            algebra exclusively (only requires wdg, rwdg, cmpl and dual).
            //            Shows that geometric product is not more fundamental, since
            //            it can be derived from other primitives.

            // Lambda for calculating left-hand side: rwdg(cmpl(c),a) =
            // cmpl(wdg(a,cmpl(a)))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_cmpl = coeff; // cmpl(cmpl(coeff)) is identity transformation
                auto rhs_cmpl = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_cmpl, rhs_cmpl, wdg_str()),
                    wdg_pga2dp_rules);
                auto lhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side: rwdg(b, bulk_dual(c)) =
            // cmpl(wdg(cmpl(b),cmpl(bulk_dual(c)))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_cmpl = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
                auto rhs_cmpl = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, bulk_dual_pga2dp_rules), cmpl_pga2dp_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_cmpl, rhs_cmpl, wdg_str()),
                    wdg_pga2dp_rules);
                auto rhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_transwedge_geometric_product(
                mv2dp_basis, mv2dp_basis_kvec, wdg_pga2dp_rules, get_lhs, get_rhs);

            // now print the resulting product table for each order k
            fmt::println("pga2dp geometric product (alternative definition) - "
                         "intermediate results:");
            fmt::println("");
            print_product_tables_by_grade(tab, mv2dp_basis_kvec);

            return tab_res;
        }

        else if (product_name == "twdg1") {

            // calculation of the geometric product based on the "transwedge product"
            // as proposed by Lengyel (https://terathon.com/blog/transwedge-product.html)

            // advantage: geometric product is derived from primitives of Grassmann
            //            algebra exclusively (only requires wdg, rwdg, cmpl and dual).
            //            Shows that geometric product is not more fundamental, since
            //            it can be derived from other primitives.

            // Lambda for calculating left-hand side: rwdg(cmpl(c),a) =
            // cmpl(wdg(a,cmpl(a)))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_cmpl = coeff; // cmpl(cmpl(coeff)) is identity transformation
                auto rhs_cmpl = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_cmpl, rhs_cmpl, wdg_str()),
                    wdg_pga2dp_rules);
                auto lhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side: rwdg(b, bulk_dual(c)) =
            // cmpl(wdg(cmpl(b),cmpl(bulk_dual(c)))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_cmpl = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
                auto rhs_cmpl = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, bulk_dual_pga2dp_rules), cmpl_pga2dp_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_cmpl, rhs_cmpl, wdg_str()),
                    wdg_pga2dp_rules);
                auto rhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_transwedge_geometric_product(
                mv2dp_basis, mv2dp_basis_kvec, wdg_pga2dp_rules, get_lhs, get_rhs);

            tab_res = add_prd_tab(tab[1], tab[2]);
            tab_res = add_prd_tab(tab_res, tab[3]);
            return tab_res; // return transwedge product for k=1
        }

        else if (product_name == "cmt") {
            // Commutator product (=asymmetric part of the geometric product)
            //                   cmt(A,B) = asym(gpr(A,B))
            auto basis_tab = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str()),
                gpr_pga2dp_rules);
            return get_prd_tab_asym(basis_tab); // use the asymmetric part only
        }

        else if (product_name == "wdg") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, wdg_str()),
                wdg_pga2dp_rules);
        }

        else if (product_name == "dot") {
            // Inner product (=dot product)
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv2dp_basis, mv2dp_basis, mul_str()),
                dot_pga2dp_rules);
        }

        else if (product_name == "l_bulk_contract") {
            // Left bulk contraction: A << B = cmpl(wdg(cmpl(bulk_dual(A)), cmpl(B)))
            auto lhs =
                apply_rules_to_mv(apply_rules_to_mv(mv2dp_basis, bulk_dual_pga2dp_rules),
                                  cmpl_pga2dp_rules);
            auto rhs = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "r_bulk_contract") {
            // Right bulk contraction: A >> B = cmpl(wdg(cmpl(A), cmpl(bulk_dual(B))))
            auto lhs = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto rhs =
                apply_rules_to_mv(apply_rules_to_mv(mv2dp_basis, bulk_dual_pga2dp_rules),
                                  cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }


        else if (product_name == "l_weight_contract") {
            // Left weight contraction: A << B = cmpl(wdg(cmpl(weight_dual(A)),
            // cmpl(B)))
            auto lhs = apply_rules_to_mv(
                apply_rules_to_mv(mv2dp_basis, weight_dual_pga2dp_rules),
                cmpl_pga2dp_rules);
            auto rhs = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "r_weight_contract") {
            // Right weight contraction: A >> B = cmpl(wdg(cmpl(A),
            // cmpl(weight_dual(B))))
            auto lhs = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto rhs = apply_rules_to_mv(
                apply_rules_to_mv(mv2dp_basis, weight_dual_pga2dp_rules),
                cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "l_bulk_expand") {
            // Left bulk expansion: A <> B = wdg(bulk_dual(A), B)
            auto lhs = apply_rules_to_mv(mv2dp_basis, bulk_dual_pga2dp_rules);
            auto rhs = mv2dp_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_pga2dp_rules);
        }

        else if (product_name == "r_bulk_expand") {
            // Right bulk expansion: A >< B = wdg(A, bulk_dual(B))
            auto lhs = mv2dp_basis;
            auto rhs = apply_rules_to_mv(mv2dp_basis, bulk_dual_pga2dp_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_pga2dp_rules);
        }

        else if (product_name == "l_weight_expand") {
            // Left weight expansion: A <> B = wdg(weight_dual(A), B)
            auto lhs = apply_rules_to_mv(mv2dp_basis, weight_dual_pga2dp_rules);
            auto rhs = mv2dp_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_pga2dp_rules);
        }

        else if (product_name == "r_weight_expand") {
            // Right weight expansion: A >< B = wdg(A, weight_dual(B))
            auto lhs = mv2dp_basis;
            auto rhs = apply_rules_to_mv(mv2dp_basis, weight_dual_pga2dp_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_pga2dp_rules);
        }

        else if (product_name == "rgpr") {
            // Regressive geometric product: rgpr(A,B) = cmpl(gpr(cmpl(A), cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "rgpr (alternative)") {

            // calculation of the regressive geometric product based on the "anti
            // transwedge product" as proposed by Lengyel
            // (https://terathon.com/blog/transwedge-product.html)

            // advantage: regressive geometric product is derived from primitives of
            //            Grassmann algebra exclusively (only requires wdg, rwdg, cmpl and
            //            dual). Shows that geometric product is not more fundamental,
            //            since it can be derived from other primitives.

            // Lambda for calculating left-hand side: cmpl(wdg(c,a))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = coeff;
                auto rhs_arg = mv2dp_basis;
                auto lhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_pga2dp_rules),
                    cmpl_pga2dp_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side:
            // cmpl(wdg(b,cmpl(dual(c))))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = mv2dp_basis;
                auto rhs_arg = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, bulk_dual_pga2dp_rules), cmpl_pga2dp_rules);
                auto rhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_pga2dp_rules),
                    cmpl_pga2dp_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_regressive_transwedge_geometric_product(
                mv2dp_basis, mv2dp_basis_kvec, wdg_pga2dp_rules, cmpl_pga2dp_rules,
                get_lhs, get_rhs);

            // now print the resulting product table for each order k
            fmt::println("pga2dp regressive geometric product (alternative definition) - "
                         "intermediate results:");
            fmt::println("");
            print_product_tables_by_grade(tab, mv2dp_basis_kvec);

            return tab_res;
        }

        else if (product_name == "rtwdg1") {

            // calculation of the regressive geometric product based on the "anti
            // transwedge product" as proposed by Lengyel
            // (https://terathon.com/blog/transwedge-product.html)

            // advantage: regressive geometric product is derived from primitives of
            //            Grassmann algebra exclusively (only requires wdg, rwdg, cmpl and
            //            dual). Shows that geometric product is not more fundamental,
            //            since it can be derived from other primitives.

            // Lambda for calculating left-hand side: cmpl(wdg(c,a))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = coeff;
                auto rhs_arg = mv2dp_basis;
                auto lhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_pga2dp_rules),
                    cmpl_pga2dp_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side:
            // cmpl(wdg(b,cmpl(dual(c))))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = mv2dp_basis;
                auto rhs_arg = apply_rules_to_mv(
                    apply_rules_to_mv(coeff, bulk_dual_pga2dp_rules), cmpl_pga2dp_rules);
                auto rhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_pga2dp_rules),
                    cmpl_pga2dp_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_regressive_transwedge_geometric_product(
                mv2dp_basis, mv2dp_basis_kvec, wdg_pga2dp_rules, cmpl_pga2dp_rules,
                get_lhs, get_rhs);

            tab_res = add_prd_tab(tab[1], tab[2]);
            tab_res = add_prd_tab(tab_res, tab[3]);
            return tab_res; // return regressive transwedge product for k=1
        }

        else if (product_name == "rcmt") {
            // Commutator product: cmt(A,B) = asym(gpr(A,B))
            // Regressive commutator product:
            //                    rcmt(A,B) = asym(rgpr(A,B))
            //                    rcmt(A,B) = asym(cmpl(gpr(cmpl(A),cmpl(B))))
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_pga2dp_rules);
            auto full_tab = apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
            return get_prd_tab_asym(full_tab);
        }

        else if (product_name == "rwdg") {
            // Regressive wedge: rwdg(A,B) = cmpl(wdg(cmpl(A), cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str()),
                wdg_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "rdot") {
            // Regressive inner product: rdot(A,B) = cmpl(dot(cmpl(A), cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                dot_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }

        else if (product_name == "sandwich_rgpr") {
            // Regressive sandwich product basis table - same as rgpr for PGA2DP
            auto basis_cmpl_func = apply_rules_to_mv(mv2dp_basis, cmpl_pga2dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_pga2dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, cmpl_pga2dp_rules);
        }
    }

    else if (algebra.name == "pga3dp") {

        if (product_name == "gpr") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, mul_str()),
                gpr_pga3dp_rules);
        }

        else if (product_name == "gpr (alternative)") {

            // calculation of the geometric product based on the "transwedge product"
            // as proposed by Lengyel (https://terathon.com/blog/transwedge-product.html)

            // advantage: geometric product is derived from primitives of Grassmann
            //            algebra exclusively (only requires wdg, rwdg, cmpl and dual).
            //            Shows that geometric product is not more fundamental, since
            //            it can be derived from other primitives.

            // Lambda for calculating left-hand side: rwdg(l_cmpl(c),a) =
            // l_cmpl(wdg(a,r_cmpl(a)))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_rcmpl =
                    coeff; // r_cmpl(l_cmpl(coeff)) is identity transformation
                auto rhs_rcmpl = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_rcmpl, wdg_str()),
                    wdg_pga3dp_rules);
                auto lhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side: rwdg(b, r_bulk_dual(c)) =
            // l_cmpl(wdg(r_cmpl(b),r_cmpl(r_bulk_dual(c)))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_rcmpl = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
                auto rhs_rcmpl =
                    apply_rules_to_mv(apply_rules_to_mv(coeff, r_bulk_dual_pga3dp_rules),
                                      r_cmpl_pga3dp_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_rcmpl, wdg_str()),
                    wdg_pga3dp_rules);
                auto rhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_transwedge_geometric_product(
                mv3dp_basis, mv3dp_basis_kvec, wdg_pga3dp_rules, get_lhs, get_rhs);

            // now print the resulting product table for each order k
            fmt::println("pga3dp geometric product (alternative definition) - "
                         "intermediate results:");
            fmt::println("");
            print_product_tables_by_grade(tab, mv3dp_basis_kvec);

            return tab_res;
        }

        else if (product_name == "twdg1") {

            // calculation of the geometric product based on the "transwedge product"
            // as proposed by Lengyel (https://terathon.com/blog/transwedge-product.html)

            // advantage: geometric product is derived from primitives of Grassmann
            //            algebra exclusively (only requires wdg, rwdg, cmpl and dual).
            //            Shows that geometric product is not more fundamental, since
            //            it can be derived from other primitives.

            // Lambda for calculating left-hand side: rwdg(l_cmpl(c),a) =
            // l_cmpl(wdg(a,r_cmpl(a)))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_rcmpl =
                    coeff; // r_cmpl(l_cmpl(coeff)) is identity transformation
                auto rhs_rcmpl = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_rcmpl, wdg_str()),
                    wdg_pga3dp_rules);
                auto lhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side: rwdg(b, r_bulk_dual(c)) =
            // l_cmpl(wdg(r_cmpl(b),r_cmpl(r_bulk_dual(c)))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_rcmpl = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
                auto rhs_rcmpl =
                    apply_rules_to_mv(apply_rules_to_mv(coeff, r_bulk_dual_pga3dp_rules),
                                      r_cmpl_pga3dp_rules);
                auto basis_tab_with_rules = apply_rules_to_tab(
                    mv_coeff_to_coeff_prd_tab(lhs_rcmpl, rhs_rcmpl, wdg_str()),
                    wdg_pga3dp_rules);
                auto rhs_tab =
                    apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_transwedge_geometric_product(
                mv3dp_basis, mv3dp_basis_kvec, wdg_pga3dp_rules, get_lhs, get_rhs);

            tab_res = add_prd_tab(tab[1], tab[2]);
            tab_res = add_prd_tab(tab_res, tab[3]);
            tab_res = add_prd_tab(tab_res, tab[4]);
            return tab_res; // return transwedge product for k=1
        }

        else if (product_name == "cmt") {
            // Commutator product (=asymmetric part of the geometric product)
            //                   cmt(A,B) = asym(gpr(A,B))
            auto basis_tab = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, mul_str()),
                gpr_pga3dp_rules);
            return get_prd_tab_asym(basis_tab); // use the asymmetric part only
        }

        else if (product_name == "wdg") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, wdg_str()),
                wdg_pga3dp_rules);
        }

        else if (product_name == "dot") {
            // Inner product (=dot product)
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mv3dp_basis, mv3dp_basis, mul_str()),
                dot_pga3dp_rules);
        }

        else if (product_name == "l_bulk_contract") {
            // Left bulk contraction: l_bulk_contract(A,B) =
            // l_cmpl(wdg(r_cmpl(l_bulk_dual(A)), r_cmpl(B)))
            auto lhs = apply_rules_to_mv(
                apply_rules_to_mv(mv3dp_basis, l_bulk_dual_pga3dp_rules),
                r_cmpl_pga3dp_rules);
            auto rhs = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
        }

        else if (product_name == "r_bulk_contract") {
            // Right bulk contraction:
            // r_bulk_contract(A,B) = l_cmpl(wdg(r_cmpl(A),
            // r_cmpl(r_bulk_dual(B))))
            auto lhs = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
            auto rhs = apply_rules_to_mv(
                apply_rules_to_mv(mv3dp_basis, r_bulk_dual_pga3dp_rules),
                r_cmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
        }

        else if (product_name == "l_weight_contract") {
            // Left weight contraction: l_weight_contract(A,B) =
            // l_cmpl(wdg(r_cmpl(l_weight_dual(A)), r_cmpl(B)))
            auto lhs = apply_rules_to_mv(
                apply_rules_to_mv(mv3dp_basis, l_weight_dual_pga3dp_rules),
                r_cmpl_pga3dp_rules);
            auto rhs = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
        }

        else if (product_name == "r_weight_contract") {
            // Right weight contraction: r_weight_contract(A,B) =
            // l_cmpl(wdg(r_cmpl(A), r_cmpl(r_weight_dual(B))))
            auto lhs = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
            auto rhs = apply_rules_to_mv(
                apply_rules_to_mv(mv3dp_basis, r_weight_dual_pga3dp_rules),
                r_cmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
        }

        else if (product_name == "l_bulk_expand") {
            // Left bulk expansion: l_bulk_expand(A,B) = wdg(l_bulk_dual(A), B)
            auto lhs = apply_rules_to_mv(mv3dp_basis, l_bulk_dual_pga3dp_rules);
            auto rhs = mv3dp_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_pga3dp_rules);
        }

        else if (product_name == "r_bulk_expand") {
            // Right bulk expansion: r_bulk_expand(A,B) = wdg(A,
            // r_bulk_dual(B))
            auto lhs = mv3dp_basis;
            auto rhs = apply_rules_to_mv(mv3dp_basis, r_bulk_dual_pga3dp_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_pga3dp_rules);
        }

        else if (product_name == "l_weight_expand") {
            // Left weight expansion:
            // l_weight_expand(A,B) = wdg(l_weight_dual(A), B)
            auto lhs = apply_rules_to_mv(mv3dp_basis, l_weight_dual_pga3dp_rules);
            auto rhs = mv3dp_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_pga3dp_rules);
        }

        else if (product_name == "r_weight_expand") {
            // Right weight expansion: r_weight_expand(A,B) = wdg(A,
            // r_weight_dual(B))
            auto lhs = mv3dp_basis;
            auto rhs = apply_rules_to_mv(mv3dp_basis, r_weight_dual_pga3dp_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_pga3dp_rules);
        }

        else if (product_name == "rgpr") {
            // Regressive geometric: rgpr(A,B) = l_cmpl(gpr(r_cmpl(A), r_cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
        }

        else if (product_name == "rgpr (alternative)") {

            // calculation of the regressive geometric product based on the "anti
            // transwedge product" as proposed by Lengyel
            // (https://terathon.com/blog/transwedge-product.html)

            // advantage: regressive geometric product is derived from primitives of
            //            Grassmann algebra exclusively (only requires wdg, rwdg, cmpl and
            //            dual). Shows that geometric product is not more fundamental,
            //            since it can be derived from other primitives.

            // Lambda for calculating left-hand side: l_cmpl(wdg(c,a))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = coeff;
                auto rhs_arg = mv3dp_basis;
                auto lhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_pga3dp_rules),
                    l_cmpl_pga3dp_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side:
            // l_cmpl(wdg(b,r_cmpl(r_bulk_dual(c))))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = mv3dp_basis;
                auto rhs_arg =
                    apply_rules_to_mv(apply_rules_to_mv(coeff, r_bulk_dual_pga3dp_rules),
                                      r_cmpl_pga3dp_rules);
                auto rhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_pga3dp_rules),
                    l_cmpl_pga3dp_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_regressive_transwedge_geometric_product(
                mv3dp_basis, mv3dp_basis_kvec, wdg_pga3dp_rules, r_cmpl_pga3dp_rules,
                get_lhs, get_rhs);

            // now print the resulting product table for each order k
            fmt::println("pga3dp regressive geometric product (alternative definition) - "
                         "intermediate results:");
            fmt::println("");
            print_product_tables_by_grade(tab, mv3dp_basis_kvec);

            return tab_res;
        }

        else if (product_name == "rtwdg1") {

            // calculation of the regressive geometric product based on the "anti
            // transwedge product" as proposed by Lengyel
            // (https://terathon.com/blog/transwedge-product.html)

            // advantage: regressive geometric product is derived from primitives of
            //            Grassmann algebra exclusively (only requires wdg, rwdg, cmpl and
            //            dual). Shows that geometric product is not more fundamental,
            //            since it can be derived from other primitives.

            // Lambda for calculating left-hand side: l_cmpl(wdg(c,a))
            auto get_lhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = coeff;
                auto rhs_arg = mv3dp_basis;
                auto lhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_pga3dp_rules),
                    l_cmpl_pga3dp_rules);
                return lhs_tab[index];
            };

            // Lambda for calculating right-hand side:
            // l_cmpl(wdg(b,r_cmpl(r_bulk_dual(c))))
            auto get_rhs = [&](mvec_coeff const& coeff, size_t index) {
                auto lhs_arg = mv3dp_basis;
                auto rhs_arg =
                    apply_rules_to_mv(apply_rules_to_mv(coeff, r_bulk_dual_pga3dp_rules),
                                      r_cmpl_pga3dp_rules);
                auto rhs_tab = apply_rules_to_tab(
                    apply_rules_to_tab(
                        mv_coeff_to_coeff_prd_tab(lhs_arg, rhs_arg, wdg_str()),
                        wdg_pga3dp_rules),
                    l_cmpl_pga3dp_rules);
                mvec_coeff rhs = rhs_tab[index];
                for (size_t i = 0; i < coeff.size(); ++i) {
                    rhs[i] = rhs_tab[i][index];
                }
                return rhs;
            };

            auto [tab_res, tab] = calculate_regressive_transwedge_geometric_product(
                mv3dp_basis, mv3dp_basis_kvec, wdg_pga3dp_rules, r_cmpl_pga3dp_rules,
                get_lhs, get_rhs);

            tab_res = add_prd_tab(tab[1], tab[2]);
            tab_res = add_prd_tab(tab_res, tab[3]);
            tab_res = add_prd_tab(tab_res, tab[4]);
            return tab_res; // return regressive transwedge product for k=1
        }

        else if (product_name == "rcmt") {
            // Commutator product: cmt(A,B) = asym(gpr(A,B))
            // Regressive commutator product:
            //                    rcmt(A,B) = asym(rgpr(A,B))
            //                    rcmt(A,B) = asym(l_cmpl(gpr(r_cmpl(A),r_cmpl(B))))
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_pga3dp_rules);
            auto full_tab = apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
            return get_prd_tab_asym(full_tab);
        }

        else if (product_name == "rwdg") {
            // Regressive wedge: rwdg(A,B) = l_cmpl(wdg(r_cmpl(A), r_cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str()),
                wdg_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
        }

        else if (product_name == "rdot") {
            // Regressive inner: rdot(A,B) = l_cmpl(dot(r_cmpl(A), r_cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                dot_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
        }

        else if (product_name == "sandwich_rgpr") {
            // Regressive sandwich product basis table - same as rgpr for PGA3DP
            auto basis_cmpl_func = apply_rules_to_mv(mv3dp_basis, r_cmpl_pga3dp_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_pga3dp_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_pga3dp_rules);
        }
    }

    else if (algebra.name == "sta4d") {

        if (product_name == "gpr") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mvsta4d_basis, mvsta4d_basis, mul_str()),
                gpr_sta4d_rules);
        }

        else if (product_name == "cmt") {
            // Commutator product (=asymmetric part of the geometric product)
            //                   cmt(A,B) = asym(gpr(A,B))
            auto basis_tab = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mvsta4d_basis, mvsta4d_basis, mul_str()),
                gpr_sta4d_rules);
            return get_prd_tab_asym(basis_tab); // use the asymmetric part only
        }

        else if (product_name == "wdg") {
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mvsta4d_basis, mvsta4d_basis, wdg_str()),
                wdg_sta4d_rules);
        }

        else if (product_name == "dot") {
            // Inner product (=dot product)
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mvsta4d_basis, mvsta4d_basis, mul_str()),
                dot_sta4d_rules);
        }

        else if (product_name == "l_contract") {
            // Left contraction: l_contract(A,B) =
            // l_cmpl(wdg(r_cmpl(l_dual(A)), r_cmpl(B)))
            auto lhs = apply_rules_to_mv(
                apply_rules_to_mv(mvsta4d_basis, l_dual_sta4d_rules), r_cmpl_sta4d_rules);
            auto rhs = apply_rules_to_mv(mvsta4d_basis, r_cmpl_sta4d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_sta4d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_sta4d_rules);
        }

        else if (product_name == "r_contract") {
            // Right contraction:
            // r_contract(A,B) = l_cmpl(wdg(r_cmpl(A),
            // r_cmpl(r_dual(B))))
            auto lhs = apply_rules_to_mv(mvsta4d_basis, r_cmpl_sta4d_rules);
            auto rhs = apply_rules_to_mv(
                apply_rules_to_mv(mvsta4d_basis, r_dual_sta4d_rules), r_cmpl_sta4d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()), wdg_sta4d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_sta4d_rules);
        }

        else if (product_name == "l_expand") {
            // Left expansion: l_expand(A,B) = wdg(l_dual(A), B)
            auto lhs = apply_rules_to_mv(mvsta4d_basis, l_dual_sta4d_rules);
            auto rhs = mvsta4d_basis;
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_sta4d_rules);
        }

        else if (product_name == "r_expand") {
            // Right expansion: r_expand(A,B) = wdg(A, r_dual(B))
            auto lhs = mvsta4d_basis;
            auto rhs = apply_rules_to_mv(mvsta4d_basis, r_dual_sta4d_rules);
            return apply_rules_to_tab(mv_coeff_to_coeff_prd_tab(lhs, rhs, wdg_str()),
                                      wdg_sta4d_rules);
        }

        else if (product_name == "rgpr") {
            // Regressive geometric: rgpr(A,B) = l_cmpl(gpr(r_cmpl(A), r_cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mvsta4d_basis, r_cmpl_sta4d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_sta4d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_sta4d_rules);
        }

        else if (product_name == "rcmt") {
            // Commutator product: cmt(A,B) = asym(gpr(A,B))
            // Regressive commutator product:
            //                    rcmt(A,B) = asym(rgpr(A,B))
            //                    rcmt(A,B) = asym(l_cmpl(gpr(r_cmpl(A),r_cmpl(B))))
            auto basis_cmpl_func = apply_rules_to_mv(mvsta4d_basis, r_cmpl_sta4d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                gpr_sta4d_rules);
            auto full_tab = apply_rules_to_tab(basis_tab_with_rules, l_cmpl_sta4d_rules);
            return get_prd_tab_asym(full_tab);
        }

        else if (product_name == "rwdg") {
            // Regressive wedge: rwdg(A,B) = l_cmpl(wdg(r_cmpl(A), r_cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mvsta4d_basis, r_cmpl_sta4d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, wdg_str()),
                wdg_sta4d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_sta4d_rules);
        }

        else if (product_name == "rdot") {
            // Regressive inner: rdot(A,B) = l_cmpl(dot(r_cmpl(A), r_cmpl(B)))
            auto basis_cmpl_func = apply_rules_to_mv(mvsta4d_basis, r_cmpl_sta4d_rules);
            auto basis_tab_with_rules = apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(basis_cmpl_func, basis_cmpl_func, mul_str()),
                dot_sta4d_rules);
            return apply_rules_to_tab(basis_tab_with_rules, l_cmpl_sta4d_rules);
        }

        else if (product_name == "sandwich_gpr") {
            // Sandwich product basis table - same as gpr for sta4d
            return apply_rules_to_tab(
                mv_coeff_to_coeff_prd_tab(mvsta4d_basis, mvsta4d_basis, mul_str()),
                gpr_sta4d_rules);
        }
    }

    // Add other algebras as needed...
    // For now, throw an error for unimplemented combinations
    throw std::invalid_argument("Unsupported product: " + algebra.name +
                                "::" + product_name);
}

filter_2d ConfigurableGenerator::get_filter_2d(AlgebraData const& algebra,
                                               std::string const& filter_name)
{
    auto it = algebra.filters_2d.find(filter_name);
    if (it == algebra.filters_2d.end()) {
        throw std::invalid_argument("Unknown 2D filter: " + filter_name);
    }
    return it->second;
}

filter_3d ConfigurableGenerator::get_filter_3d(AlgebraData const& algebra,
                                               std::string const& filter_name)
{
    auto it = algebra.filters_3d.find(filter_name);
    if (it == algebra.filters_3d.end()) {
        throw std::invalid_argument("Unknown 3D filter: " + filter_name);
    }
    return it->second;
}

filter_4d ConfigurableGenerator::get_filter_4d(AlgebraData const& algebra,
                                               std::string const& filter_name)
{
    auto it = algebra.filters_4d.find(filter_name);
    if (it == algebra.filters_4d.end()) {
        throw std::invalid_argument("Unknown 4D filter: " + filter_name);
    }
    return it->second;
}

void ConfigurableGenerator::print_product_header(AlgebraData const& algebra,
                                                 ProductConfig const& config)
{
    fmt::println("{} {} - basis product table:", algebra.name, config.display_name);
}

void ConfigurableGenerator::print_basis_table(AlgebraData const& algebra,
                                              ProductConfig const& config,
                                              prd_table const& basis_tab,
                                              bool show_sym_asym)
{
    print_product_header(algebra, config);
    print_prd_tab_with_headers(basis_tab, algebra.basis);

    if (show_sym_asym) {
        fmt::println("\nsymmetric part:");
        print_prd_tab_with_headers(get_prd_tab_sym(basis_tab), algebra.basis);
        fmt::println("\nasymmetric part:");
        print_prd_tab_with_headers(get_prd_tab_asym(basis_tab), algebra.basis);
    }
}

void ConfigurableGenerator::print_case_header(AlgebraData const& algebra,
                                              ProductConfig const& config,
                                              std::string const& case_name)
{
    // reference format, e.g. "ega2d geometric product mv * mv -> mv:"
    fmt::println("{} {} {}:", algebra.name, config.display_name, case_name);
}

void ConfigurableGenerator::print_case_result(mvec_coeff const& result,
                                              mvec_coeff const& basis)
{
    // Use EXISTING print function from reference implementation
    print_mvec(result, basis);
}

void ConfigurableGenerator::print_transformed_result(mvec_coeff const& result,
                                                     mvec_coeff const& basis,
                                                     AlgebraData const& algebra,
                                                     ProductConfig const& config)
{
    // Create algebra-specific header using correct display name
    fmt::println("{} {} (after transformation):", algebra.name, config.display_name);

    try {
        // Convert each component to string for transformation
        std::vector<std::string> component_expressions;
        for (size_t i = 0; i < result.size() && i < basis.size(); ++i) {
            if (!result[i].empty()) {
                component_expressions.push_back(result[i]);
            }
        }

        // Transform using the sandwich transformer
        auto transformed_expressions = SandwichTransformer::transformSandwichMultivector(
            component_expressions, algebra.name);

        // Convert back to mvec_coeff format with explicit "0" for empty results
        mvec_coeff transformed_result(basis.size());
        for (size_t i = 0; i < basis.size(); ++i) {
            if (i < transformed_expressions.size() &&
                !transformed_expressions[i].empty() &&
                transformed_expressions[i] != "0") {
                transformed_result[i] = transformed_expressions[i];
            }
            else {
                // Use explicit "0" for empty/zero components to make copy-paste ready
                transformed_result[i] = "0";
            }
        }

        // Apply coefficient alignment before printing
        apply_coefficient_alignment(transformed_result, algebra.name);

        // Use the same print function as original results for consistent
        // right-alignment
        print_mvec(transformed_result, basis);
    }
    catch (const std::exception& e) {
        fmt::println("  Transformation error: {}", e.what());
    }
    fmt::println("");
}

void ConfigurableGenerator::apply_coefficient_alignment(mvec_coeff& expressions,
                                                        std::string const& algebra_name)
{
    if (expressions.empty()) return;

    // Dynamically find geometric variables (avoid rotor coefficients like R.c)
    std::regex var_regex(R"([a-zA-Z]+\.[a-zA-Z]+)");
    std::set<std::string> all_variables;

    for (const auto& expr : expressions) {
        if (expr.empty() || expr == "0") continue;

        std::sregex_iterator iter(expr.begin(), expr.end(), var_regex);
        std::sregex_iterator end;
        for (; iter != end; ++iter) {
            std::string var = iter->str();
            // Only include variables that represent geometric objects (v.*, B.*, t.*)
            // Skip rotor/motor coefficients (R.*)
            if (!var.starts_with("R.")) {
                all_variables.insert(var);
            }
        }
    }

    if (all_variables.empty()) return;

    // Get appropriate geometric patterns based on algebra
    GeometricVariablePatterns patterns;
    if (algebra_name == "pga3dp") {
        patterns = GeometricVariablePatterns::createPGA3DPPatterns();
    }
    else if (algebra_name == "pga2dp") {
        patterns = GeometricVariablePatterns::createPGA2DPPatterns();
    }
    else if (algebra_name == "ega3d") {
        patterns = GeometricVariablePatterns::createEGA3DPatterns();
    }
    else if (algebra_name == "ega2d") {
        patterns = GeometricVariablePatterns::createEGA2DPatterns();
    }
    else if (algebra_name == "sta4d") {
        patterns = GeometricVariablePatterns::createSTA4DPatterns();
    }
    else {
        // Default to PGA3DP patterns as fallback
        patterns = GeometricVariablePatterns::createPGA3DPPatterns();
        fmt::println("WARNING: Default pattern used. Provide algebra specific pattern.");
    }

    // Create factors map from variables (all with power 1)
    std::map<std::string, int> factors_map;
    for (const auto& var : all_variables) {
        factors_map[var] = 1;
    }

    // Use canonical ordering
    auto sorted_pairs = GAAlgebraRules::getSortedVariablePairs(factors_map, patterns);

    // Extract just the variable names in canonical order
    std::vector<std::string> sorted_variables;
    for (const auto& pair : sorted_pairs) {
        sorted_variables.push_back(pair.first);
    }

    if (sorted_variables.empty()) return;

    // Step 1: Calculate maximum coefficient space needed for each variable across ALL
    // expressions
    std::map<std::string, size_t> max_coeff_widths;

    for (const auto& var : sorted_variables) {
        size_t max_width = 0;
        bool found_in_any_expr = false;

        for (const auto& expr : expressions) {
            if (expr.empty() || expr == "0") continue;

            std::string pattern = " * " + var;
            size_t pos = expr.find(pattern);

            if (pos != std::string::npos) {
                found_in_any_expr = true;

                // Find coefficient start by walking backwards respecting parentheses
                size_t coeff_start = pos;
                int paren_depth = 0;

                while (coeff_start > 0) {
                    char c = expr[coeff_start - 1];
                    if (c == ')') {
                        paren_depth++;
                    }
                    else if (c == '(') {
                        paren_depth--;
                    }
                    else if (paren_depth == 0) {
                        if (c == '+' || (c == '-' && coeff_start > 1 &&
                                         expr[coeff_start - 2] == ' ')) {
                            if (c == '-') {
                                coeff_start--;
                            }
                            break;
                        }
                        else if (c == '[') {
                            break;
                        }
                    }
                    coeff_start--;
                }

                // Skip leading spaces for calculation
                while (coeff_start < pos && expr[coeff_start] == ' ') {
                    coeff_start++;
                }

                if (coeff_start < pos) {
                    size_t coeff_width = pos - coeff_start;
                    max_width = std::max(max_width, coeff_width);
                }
            }
        }

        if (found_in_any_expr) {
            max_coeff_widths[var] = max_width;
        }
    }

    // Step 2: For each expression, rebuild it with proper column alignment
    for (auto& expr : expressions) {
        if (expr.empty() || expr == "0") {
            expr = " 0";
            continue;
        }

        // Preserve the prefix (like "[0] = ")
        std::string prefix = "";
        size_t bracket_end = expr.find("] = ");
        if (bracket_end != std::string::npos) {
            prefix = expr.substr(0, bracket_end + 4);
        }

        // Build new expression with proper column alignment
        std::string result = " "; // Start with leading space

        for (const auto& var : sorted_variables) {
            if (max_coeff_widths.find(var) == max_coeff_widths.end()) continue;

            std::string pattern = " * " + var;
            size_t pos = expr.find(pattern);

            // Also try to match patterns for cases like ") * var"
            if (pos == std::string::npos) {
                // Try pattern with closing parenthesis
                std::string paren_pattern = ") * " + var;
                size_t paren_pos = expr.find(paren_pattern);
                if (paren_pos != std::string::npos) {
                    pos = paren_pos + 1; // Position at the space before *
                }
            }


            if (pos != std::string::npos) {
                // Variable is present - extract and format coefficient
                size_t coeff_start = pos;
                int paren_depth = 0;

                while (coeff_start > 0) {
                    char c = expr[coeff_start - 1];
                    if (c == ')') {
                        paren_depth++;
                    }
                    else if (c == '(') {
                        paren_depth--;
                    }
                    else if (paren_depth == 0) {
                        if (c == '+' || (c == '-' && coeff_start > 1 &&
                                         expr[coeff_start - 2] == ' ')) {
                            if (c == '-') {
                                coeff_start--;
                            }
                            break;
                        }
                        else if (c == '[') {
                            break;
                        }
                    }
                    coeff_start--;
                }

                // Skip leading spaces
                while (coeff_start < pos && expr[coeff_start] == ' ') {
                    coeff_start++;
                }

                // Extract coefficient
                std::string coefficient = expr.substr(coeff_start, pos - coeff_start);


                // Pad coefficient to maximum width for this variable
                size_t target_width = max_coeff_widths[var];
                size_t padding_needed = (coefficient.length() < target_width)
                                            ? target_width - coefficient.length()
                                            : 0;


                // Check if any subsequent variables exist in this expression
                bool has_next_var = false;
                auto current_var_it =
                    std::find(sorted_variables.begin(), sorted_variables.end(), var);
                for (auto next_it = current_var_it + 1; next_it != sorted_variables.end();
                     ++next_it) {
                    if (max_coeff_widths.find(*next_it) != max_coeff_widths.end()) {
                        std::string next_pattern = " * " + *next_it;
                        if (expr.find(next_pattern) != std::string::npos) {
                            has_next_var = true;
                            break;
                        }
                    }
                }

                // Add properly aligned variable term
                if (has_next_var) {
                    // Check if the next coefficient starts with a minus sign
                    bool next_coeff_has_minus = false;
                    for (auto next_it = current_var_it + 1;
                         next_it != sorted_variables.end(); ++next_it) {
                        if (max_coeff_widths.find(*next_it) != max_coeff_widths.end()) {
                            std::string next_pattern = " * " + *next_it;
                            size_t next_pos = expr.find(next_pattern);
                            if (next_pos != std::string::npos) {
                                // Find the next coefficient start
                                size_t next_coeff_start = next_pos;
                                int next_paren_depth = 0;
                                while (next_coeff_start > 0) {
                                    char c = expr[next_coeff_start - 1];
                                    if (c == ')') {
                                        next_paren_depth++;
                                    }
                                    else if (c == '(') {
                                        next_paren_depth--;
                                    }
                                    else if (next_paren_depth == 0) {
                                        if (c == '+' ||
                                            (c == '-' && next_coeff_start > 1 &&
                                             expr[next_coeff_start - 2] == ' ')) {
                                            if (c == '-') {
                                                next_coeff_start--;
                                            }
                                            break;
                                        }
                                        else if (c == '[') {
                                            break;
                                        }
                                    }
                                    next_coeff_start--;
                                }

                                // Skip leading spaces and check for minus
                                while (next_coeff_start < next_pos &&
                                       expr[next_coeff_start] == ' ') {
                                    next_coeff_start++;
                                }
                                if (next_coeff_start < next_pos &&
                                    expr[next_coeff_start] == '-') {
                                    next_coeff_has_minus = true;
                                }
                                break;
                            }
                        }
                    }

                    if (next_coeff_has_minus) {
                        // Add 2 extra spaces when the next coefficient has a minus to
                        // match target format
                        result += std::string(padding_needed + 2, ' ') + coefficient +
                                  " * " + var + "   ";
                    }
                    else {
                        result += std::string(padding_needed, ' ') + coefficient + " * " +
                                  var + " + ";
                    }
                }
                else {
                    result +=
                        std::string(padding_needed, ' ') + coefficient + " * " + var;
                }
            }
            else {
                // Variable is missing - check if any subsequent variables exist in
                // this expression
                bool has_next_var = false;
                auto current_var_it =
                    std::find(sorted_variables.begin(), sorted_variables.end(), var);
                for (auto next_it = current_var_it + 1; next_it != sorted_variables.end();
                     ++next_it) {
                    if (max_coeff_widths.find(*next_it) != max_coeff_widths.end()) {
                        std::string next_pattern = " * " + *next_it;
                        if (expr.find(next_pattern) != std::string::npos) {
                            has_next_var = true;
                            break;
                        }
                    }
                }

                // Reserve space for missing variable
                size_t target_coeff_width = max_coeff_widths[var];
                size_t total_column_width;
                if (has_next_var) {
                    // Variable has successors: " + coeff * var + " (with trailing " +
                    // ")
                    total_column_width = 3 + target_coeff_width + 3 + var.length() +
                                         3; // " + " + coeff + " * " + var + " + "
                }
                else {
                    // Variable is last: " + coeff * var" (no trailing " + ")
                    total_column_width = 3 + target_coeff_width + 3 +
                                         var.length(); // " + " + coeff + " * " + var
                }
                result += std::string(total_column_width, ' ');
            }
        }

        // Remove trailing " + "
        if (result.length() >= 3 && result.substr(result.length() - 3) == " + ") {
            result = result.substr(0, result.length() - 3);
        }

        // Apply prefix if it exists
        expr = prefix + result;
    }

    // Step 3: Clean up excessive leading whitespace
    for (auto& expr : expressions) {
        if (!expr.empty() && expr != "0") {
            size_t first_non_space = expr.find_first_not_of(' ');
            if (first_non_space > 1) {
                expr = " " + expr.substr(first_non_space);
            }
        }
    }
}