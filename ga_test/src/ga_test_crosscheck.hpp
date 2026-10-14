#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Support for the generated cross-checks (ga_<alg>_crosscheck_test.hpp, written by
// ga_test/utilities/gen_crosscheck_tests.py).
//
// Every typed product overload -- wdg(BiVec3dp, Vec3dp), operator*(Vec2d, MVec2d_E),
// ... -- is a specialization of one product on the full multivector, with the operand
// types fixing which coefficients can be non-zero and the declared result type fixing
// which can survive. So each typed overload must give EXACTLY what the full product
// gives on the same operands lifted into the full multivector: a dropped term, a wrong
// sign, a swapped coefficient or a result type too narrow for the true result all show
// up as a disagreement. The reference is the library's own full-multivector product, so
// this is a consistency check between the two generated layers, not an independent
// derivation of the algebra -- the identity tests in the per-algebra suites are that.
//
// The same holds for the unary operations of the *_ops_basics.hpp headers (involutions,
// complements, duals, bulk/weight, norms, normalizations): each typed overload is the
// full-multivector operation restricted to its type.
//
// Operands are drawn from a seeded generator, so a failure reproduces exactly.

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <string>
#include <type_traits>
#include <utility>

#include "doctest/doctest.h"
#include "ga/ga_value_t.hpp"

namespace xcheck {

using hd::ga::value_t;

struct rng {
    std::mt19937_64 gen;
    std::uniform_real_distribution<value_t> dist{value_t(-1), value_t(1)};

    explicit rng(std::uint64_t seed) : gen(seed) {}

    value_t operator()() { return dist(gen); }
};

// every value type stores exactly sizeof(X) / sizeof(value_t) coefficients
template <typename X> inline constexpr size_t n_coeff = sizeof(X) / sizeof(value_t);

// a random value of type X. A type that cannot be built from its coefficients directly
// (a geometric object only constructible from its underlying blade) gets an explicit
// specialization in the generated file that needs it.
template <typename X> X make(rng& r)
{
    std::array<value_t, n_coeff<X>> v{};
    for (auto& c : v)
        c = r(); // drawn in order, so the stream is the same on every compiler
    return [&]<size_t... I>(std::index_sequence<I...>) {
        return X(v[I]...);
    }(std::make_index_sequence<n_coeff<X>>{});
}

template <typename M> std::array<value_t, n_coeff<M>> coeffs(M const& m)
{
    static_assert(std::is_trivially_copyable_v<M>);
    return std::bit_cast<std::array<value_t, n_coeff<M>>>(m);
}

// the two full multivectors agree to a relative tolerance: the typed and the full
// product may sum the same terms in a different order, so exact equality is not owed
template <typename M> bool agree(M const& a, M const& b, value_t rel_tol)
{
    auto const ca = coeffs(a);
    auto const cb = coeffs(b);
    value_t scale = 1;
    for (size_t i = 0; i < ca.size(); ++i)
        scale = std::max({scale, std::abs(ca[i]), std::abs(cb[i])});
    for (size_t i = 0; i < ca.size(); ++i)
        if (!(std::abs(ca[i] - cb[i]) <= rel_tol * scale)) return false;
    return true;
}

inline constexpr value_t tol = 1.0e4 * std::numeric_limits<value_t>::epsilon();
inline constexpr int draws = 3;

// f(a, b) on the typed operands, lifted, against ref(lift(a), lift(b)), lifted
template <typename A, typename B, typename F, typename Ref, typename Lift>
void binary(rng& r, F f, Ref ref, Lift lift, char const* sig)
{
    for (int k = 0; k < draws; ++k) {
        A const a = make<A>(r);
        B const b = make<B>(r);
        auto const typed = lift(f(a, b));
        auto const full = lift(ref(lift(a), lift(b)));
        INFO(std::string(sig)); // a char const* would be logged as an address
        CHECK(agree(typed, full, tol));
    }
}

// u (x) inverse(u) == one == inverse(u) (x) u for the product (x) the inverse belongs to
template <typename A, typename Inv, typename Prd, typename Lift, typename One>
void inverse(rng& r, Inv inverse_of, Prd prd, Lift lift, One one, value_t rel_tol,
             char const* sig)
{
    for (int k = 0; k < draws; ++k) {
        A const a = make<A>(r);
        auto const ai = inverse_of(a);
        INFO(std::string(sig)); // a char const* would be logged as an address
        CHECK(agree(lift(prd(lift(a), lift(ai))), lift(one), rel_tol));
        CHECK(agree(lift(prd(lift(ai), lift(a))), lift(one), rel_tol));
    }
}

// f(a) on the typed operand against f on the lifted operand. A result that has a full
// multivector form (a blade, a versor) is compared lifted; one that has not (a norm, a
// dual number) is compared as the type both overloads return.
template <typename A, typename F, typename Lift>
void unary(rng& r, F f, Lift lift, char const* sig)
{
    for (int k = 0; k < draws; ++k) {
        A const a = make<A>(r);
        auto const typed = f(a);
        auto const full = f(lift(a));
        INFO(std::string(sig)); // a char const* would be logged as an address
        if constexpr (std::is_same_v<decltype(typed), decltype(full)>)
            CHECK(agree(typed, full, tol));
        else CHECK(agree(lift(typed), lift(full), tol));
    }
}

} // namespace xcheck
