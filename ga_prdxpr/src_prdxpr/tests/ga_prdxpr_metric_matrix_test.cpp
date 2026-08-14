// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Tests for the non-orthogonal metric (metric_matrix) machinery:
//
//   - config validation (Gram matrix, diagonalization, inverse round-trip)
//   - extended metric matrix via Gram determinants (signed permutation)
//   - dot rules with off-diagonal (partner) entries
//   - multi-term geometric product rules: diagonal-basis detour vs direct
//     non-orthogonal rewrite (both routes must agree — always-on gate), plus
//     hand-derived expected products and algebraic identities
//   - complement + dual generation (single-term, non-degenerate)
//   - multi-term basis table build + fused coefficient extraction
//
// The test config is a conformal-style null-vector basis over 2D Euclidean
// space (4D algebra, 16 basis elements): e1, e2 Euclidean; e3, e4 a null pair
// with e3*e3 = e4*e4 = 0 and e3.e4 = -1. Internal orthogonal basis:
// (e1, e2, e-, e+) with signature (+1, +1, -1, +1) and
//
//     e3 = 0.5*(e- - e+),   e4 = e- + e+
//     e- = e3 + 0.5*e4,     e+ = -e3 + 0.5*e4

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "rules/ga_prdxpr_metric_calc.hpp"
#include "rules/ga_prdxpr_rule_generator.hpp"

namespace {

AlgebraConfig make_null_basis_config()
{
    AlgebraConfig config;
    config.basis_vectors = {"e1", "e2", "e3", "e4"};
    config.metric_signature = {1, 1, 0, 0}; // diagonal of metric_matrix
    config.multivector_basis = {"1",    "e1",   "e2",   "e3",   "e4",  "e31",
                                "e32",  "e12",  "e14",  "e24",  "e34", "e314",
                                "e324", "e124", "e321", "e1234"};
    config.scalar_name = one_str();
    config.basis_prefix = "e";

    config.metric_matrix = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, -1}, {0, 0, -1, 0}};
    config.basis_change = {{1.0, 0.0, 0.0, 0.0},
                           {0.0, 1.0, 0.0, 0.0},
                           {0.0, 0.0, 0.5, -0.5},
                           {0.0, 0.0, 1.0, 1.0}};
    config.basis_change_inv = {{1.0, 0.0, 0.0, 0.0},
                               {0.0, 1.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.5},
                               {0.0, 0.0, -1.0, 0.5}};
    config.diag_signature = {1, 1, -1, 1};
    return config;
}

// find a term's coefficient in a prd_terms (0 if absent)
int coeff_of(prd_terms const& terms, std::string const& blade)
{
    for (auto const& t : terms) {
        if (t.blade == blade) return t.coeff;
    }
    return 0;
}

} // namespace

TEST_SUITE("non-orthogonal metric machinery")
{

    TEST_CASE("config validation accepts the null-basis config")
    {
        auto const config = make_null_basis_config();
        CHECK_NOTHROW(validate_metric_matrix_config(config));
    }

    TEST_CASE("config validation rejects inconsistent inputs")
    {
        SUBCASE("asymmetric metric matrix")
        {
            auto config = make_null_basis_config();
            config.metric_matrix[0][1] = 1; // breaks symmetry
            CHECK_THROWS(validate_metric_matrix_config(config));
        }
        SUBCASE("signature not equal to matrix diagonal")
        {
            auto config = make_null_basis_config();
            config.metric_signature = {1, 1, 1, 0};
            CHECK_THROWS(validate_metric_matrix_config(config));
        }
        SUBCASE("basis change does not reproduce the metric")
        {
            auto config = make_null_basis_config();
            config.diag_signature = {1, 1, 1, 1}; // wrong internal signature
            CHECK_THROWS(validate_metric_matrix_config(config));
        }
        SUBCASE("degenerate metric matrix is rejected (exact determinant)")
        {
            auto config = make_null_basis_config();
            // remove the null-pair coupling: diagonal still matches the
            // signature, matrix stays symmetric, but det(G) becomes 0 — the
            // signature-based PGA classification is skipped for matrix
            // configs, so degeneracy must be rejected up front
            config.metric_matrix[2][3] = 0;
            config.metric_matrix[3][2] = 0;
            CHECK_THROWS(validate_metric_matrix_config(config));
        }
    }

    TEST_CASE("extended metric matrix is the expected signed permutation")
    {
        auto const config = make_null_basis_config();
        auto const& basis = config.multivector_basis;
        size_t const n = basis.size();
        auto const G = calculate_extended_metric_matrix_full(config);

        // every row holds exactly one non-zero entry, and it is +1 or -1
        for (size_t i = 0; i < n; ++i) {
            int nonzero = 0;
            for (size_t j = 0; j < n; ++j) {
                int const v = G[i * n + j];
                if (v != 0) {
                    ++nonzero;
                    CHECK((v == 1 || v == -1));
                }
            }
            CHECK(nonzero == 1);
        }

        // spot checks against hand-computed Gram determinants
        auto entry = [&](std::string const& r, std::string const& c) {
            auto const ir = static_cast<size_t>(std::find(basis.begin(), basis.end(), r) -
                                                basis.begin());
            auto const ic = static_cast<size_t>(std::find(basis.begin(), basis.end(), c) -
                                                basis.begin());
            return G[ir * n + ic];
        };
        CHECK(entry("1", "1") == 1);
        CHECK(entry("e1", "e1") == 1);
        CHECK(entry("e3", "e3") == 0);
        CHECK(entry("e3", "e4") == -1);       // the null pair
        CHECK(entry("e34", "e34") == -1);     // det[[0,-1],[-1,0]] = -1
        CHECK(entry("e31", "e14") == 1);      // w-blade pairs with u-blade
        CHECK(entry("e1234", "e1234") == -1); // det(G) = -1
    }

    TEST_CASE("dot rules carry the off-diagonal partner entries")
    {
        auto const config = make_null_basis_config();
        auto const dot_rules = generate_dot_product_rules(config);

        CHECK(dot_rules.at("e1 * e1") == "1");
        CHECK(dot_rules.at("e3 * e3") == "0");
        CHECK(dot_rules.at("e3 * e4") == "-1");
        CHECK(dot_rules.at("e4 * e3") == "-1");
        CHECK(dot_rules.at("e34 * e34") == "-1");
        CHECK(dot_rules.at("e31 * e14") == "1");
        CHECK(dot_rules.at("e1234 * e1234") == "-1");
        CHECK(dot_rules.at("e1 * e2") == "0");
    }

    TEST_CASE("multi-term gpr rules: two independent routes agree pairwise")
    {
        auto const config = make_null_basis_config();
        for (auto const& a : config.multivector_basis) {
            for (auto const& b : config.multivector_basis) {
                auto const detour = multiply_basis_elements_via_diagonal(a, b, config);
                auto const direct = multiply_basis_elements_general(a, b, config);
                CHECK_MESSAGE(detour == direct, "mismatch for ", a, " * ", b, ": '",
                              prd_terms_to_string(detour), "' vs '",
                              prd_terms_to_string(direct), "'");
            }
        }
    }

    TEST_CASE("multi-term gpr rules: hand-derived products")
    {
        auto const config = make_null_basis_config();
        auto const rules = generate_geometric_product_rules_mt(config);

        // e3 * e4 = e3.e4 + e3^e4 = -1 + e34
        CHECK(prd_terms_to_string(rules.at("e3 * e4")) == "-1 + e34");
        // e4 * e3 = -1 - e34
        CHECK(prd_terms_to_string(rules.at("e4 * e3")) == "-1 - e34");
        // null vectors square to zero
        CHECK(rules.at("e3 * e3").empty());
        CHECK(rules.at("e4 * e4").empty());
        // the null-pair bivector squares to +1
        CHECK(prd_terms_to_string(rules.at("e34 * e34")) == "1");
        // Euclidean part behaves as usual
        CHECK(prd_terms_to_string(rules.at("e1 * e1")) == "1");
        CHECK(prd_terms_to_string(rules.at("e1 * e2")) == "e12");
        CHECK(prd_terms_to_string(rules.at("e12 * e12")) == "-1");
        // scalar row/column
        CHECK(prd_terms_to_string(rules.at("1 * e34")) == "e34");
        CHECK(prd_terms_to_string(rules.at("e34 * 1")) == "e34");
    }

    TEST_CASE("multi-term gpr rules: algebraic identities")
    {
        auto const config = make_null_basis_config();
        auto const rules = generate_geometric_product_rules_mt(config);
        auto const& basis = config.multivector_basis;
        size_t const n = basis.size();
        auto const G = calculate_extended_metric_matrix_full(config);

        auto pos_of = [&](std::string const& e) {
            return static_cast<size_t>(std::find(basis.begin(), basis.end(), e) -
                                       basis.begin());
        };

        SUBCASE("vector anticommutators reproduce the metric: "
                "e_i e_j + e_j e_i = 2*G[i,j]")
        {
            for (auto const& vi : config.basis_vectors) {
                for (auto const& vj : config.basis_vectors) {
                    auto const ij = rules.at(vi + " * " + vj);
                    auto const ji = rules.at(vj + " * " + vi);
                    int const scalar_sum = coeff_of(ij, config.scalar_name) +
                                           coeff_of(ji, config.scalar_name);
                    auto const si = pos_of(vi) - 1; // vector slots follow the scalar
                    auto const sj = pos_of(vj) - 1;
                    CHECK(scalar_sum == 2 * config.metric_matrix[si][sj]);
                }
            }
        }

        SUBCASE("scalar part of e_S * e_T equals the reversed extended metric entry")
        {
            // gr0(e_S ⟑ e_T) = <rev(e_S), e_T> = (-1)^(k(k-1)/2) * G_ext[S,T]
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    auto const& terms = rules.at(basis[i] + " * " + basis[j]);
                    size_t const k = parse_indices(basis[i], config.basis_prefix).size();
                    int const rev_sign = ((k * (k - 1) / 2) % 2 == 0) ? 1 : -1;
                    CHECK(coeff_of(terms, config.scalar_name) == rev_sign * G[i * n + j]);
                }
            }
        }

        SUBCASE("top-grade part of e_S * e_T equals the wedge rules")
        {
            auto const wdg_rules = generate_wedge_product_rules(config);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    auto const key_w = basis[i] + " ^ " + basis[j];
                    auto const& terms = rules.at(basis[i] + " * " + basis[j]);
                    size_t const gi = parse_indices(basis[i], config.basis_prefix).size();
                    size_t const gj = parse_indices(basis[j], config.basis_prefix).size();

                    // collect the grade-(gi+gj) part of the gpr terms
                    std::string top = zero_str();
                    for (auto const& t : terms) {
                        if (parse_indices(t.blade, config.basis_prefix).size() ==
                            gi + gj) {
                            REQUIRE((t.coeff == 1 || t.coeff == -1));
                            top = (t.coeff == 1) ? t.blade : ("-" + t.blade);
                        }
                    }
                    CHECK(top == wdg_rules.at(key_w));
                }
            }
        }

        SUBCASE("associativity on all basis triples")
        {
            // multiply a prd_terms combination by a basis element via the rules
            auto mul_right = [&](prd_terms const& terms, std::string const& e) {
                std::map<std::string, int> acc;
                for (auto const& t : terms) {
                    for (auto const& r : rules.at(t.blade + " * " + e)) {
                        acc[r.blade] += t.coeff * r.coeff;
                    }
                }
                return acc;
            };
            auto mul_left = [&](std::string const& e, prd_terms const& terms) {
                std::map<std::string, int> acc;
                for (auto const& t : terms) {
                    for (auto const& r : rules.at(e + " * " + t.blade)) {
                        acc[r.blade] += t.coeff * r.coeff;
                    }
                }
                return acc;
            };

            for (auto const& a : basis) {
                for (auto const& b : basis) {
                    auto const ab = rules.at(a + " * " + b);
                    for (auto const& c : basis) {
                        auto lhs = mul_right(ab, c);                     // (a*b)*c
                        auto rhs = mul_left(a, rules.at(b + " * " + c)); // a*(b*c)
                        std::erase_if(lhs, [](auto const& kv) { return kv.second == 0; });
                        std::erase_if(rhs, [](auto const& kv) { return kv.second == 0; });
                        CHECK_MESSAGE(lhs == rhs, "associativity failed for (", a, " * ",
                                      b, ") * ", c);
                    }
                }
            }
        }
    }

    TEST_CASE("full rule generation: complements and duals are single-term")
    {
        auto const config = make_null_basis_config();
        auto const rules = generate_algebra_rules(config);

        // gpr lives in the mt slot; the single-term slot stays empty
        CHECK(rules.geometric_product.empty());
        CHECK(rules.geometric_product_mt.size() ==
              config.multivector_basis.size() * config.multivector_basis.size());

        // even-dimensional algebra: l/r complements; non-degenerate: l/r duals
        CHECK(rules.l_cmpl.size() == config.multivector_basis.size());
        CHECK(rules.r_cmpl.size() == config.multivector_basis.size());
        CHECK(rules.l_dual.size() == config.multivector_basis.size());
        CHECK(rules.r_dual.size() == config.multivector_basis.size());
        CHECK(rules.bulk_dual.empty());   // not a degenerate (PGA) algebra
        CHECK(rules.weight_dual.empty()); // not a degenerate (PGA) algebra

        // non-degenerate metric: every dual is a single signed basis element
        for (auto const& rule_set : {rules.l_dual, rules.r_dual}) {
            for (auto const& [key, value] : rule_set) {
                CHECK(value != zero_str());
                CHECK(value.find('+') == std::string::npos);
                CHECK(value.find(' ') == std::string::npos);
            }
        }
    }

    TEST_CASE("mt basis table build and fused coefficient extraction")
    {
        auto const config = make_null_basis_config();
        auto const rules = generate_geometric_product_rules_mt(config);
        auto const& basis = config.multivector_basis;

        auto const tab = build_mt_basis_table(basis, basis, rules, mul_str());
        REQUIRE(tab.size() == basis.size());

        // the e3/e4 cell carries the two-term product
        auto pos_of = [&](std::string const& e) {
            return static_cast<size_t>(std::find(basis.begin(), basis.end(), e) -
                                       basis.begin());
        };
        CHECK(prd_terms_to_string(tab[pos_of("e3")][pos_of("e4")]) == "-1 + e34");

        SUBCASE("single-entry coefficients: a3*e3 times b4*e4")
        {
            mvec_coeff lcoeff(basis.size(), zero_str());
            mvec_coeff rcoeff(basis.size(), zero_str());
            lcoeff[pos_of("e3")] = "a3";
            rcoeff[pos_of("e4")] = "b4";
            mvec_coeff_filter const all(basis.size(), 1);

            auto const mv = get_mv_from_mt_tab(tab, lcoeff, rcoeff, basis, all, all);
            CHECK(mv[0] == "-a3 * b4");            // scalar part: -1
            CHECK(mv[pos_of("e34")] == "a3 * b4"); // bivector part: +e34
            for (size_t k = 0; k < basis.size(); ++k) {
                if (k != 0 && k != pos_of("e34")) {
                    CHECK(mv[k] == zero_str());
                }
            }
        }

        SUBCASE("sign folding with negative coefficient strings")
        {
            mvec_coeff lcoeff(basis.size(), zero_str());
            mvec_coeff rcoeff(basis.size(), zero_str());
            lcoeff[pos_of("e3")] = "-a3";
            rcoeff[pos_of("e4")] = "b4";
            mvec_coeff_filter const all(basis.size(), 1);

            auto const mv = get_mv_from_mt_tab(tab, lcoeff, rcoeff, basis, all, all);
            CHECK(mv[0] == "a3 * b4");
            CHECK(mv[pos_of("e34")] == "-a3 * b4");
        }

        SUBCASE("coefficient magnitude > 1 is emitted as a numeric prefactor")
        {
            // artificial 2-element system exercising the "2.0 * " formatting
            mvec_coeff const mini_basis{"1", "x"};
            prd_rules_mt mini_rules;
            mini_rules["1 * 1"] = {prd_term{1, "1"}};
            mini_rules["1 * x"] = {prd_term{1, "x"}};
            mini_rules["x * 1"] = {prd_term{1, "x"}};
            mini_rules["x * x"] = {prd_term{2, "x"}};

            auto const mini_tab =
                build_mt_basis_table(mini_basis, mini_basis, mini_rules, mul_str());
            mvec_coeff const lcoeff{"0", "a"};
            mvec_coeff const rcoeff{"0", "b"};
            mvec_coeff_filter const all(2, 1);

            auto const mv =
                get_mv_from_mt_tab(mini_tab, lcoeff, rcoeff, mini_basis, all, all);
            CHECK(mv[0] == "0");
            CHECK(mv[1] == "2.0 * a * b");
        }
    }
}
