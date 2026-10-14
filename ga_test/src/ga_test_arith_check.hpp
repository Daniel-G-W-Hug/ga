#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// The arithmetic contract every value type of the library keeps, asked of one type at a
// time (the storage templates in ga/detail/type_t/ and their per-algebra
// specializations):
//
//   - each compound assignment is its binary operator: a += b gives a + b, and so on
//   - a scalar multiplies from either side, and division by s is multiplication by 1/s
//   - -a is the additive inverse, and a - b is a + (-b)
//   - != is the negation of ==
//   - coeff_sq is the plain sum of the squared coefficients (no metric)
//   - a float value converts to value_t without loss (the converting constructor)
//   - in an extended-test build, division by zero is refused, not carried out
//
// operator== compares relatively (see "operator==" in the library notes), which absorbs
// the rounding difference between a / s and a * (1 / s).

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "doctest/doctest.h"
#include "ga/ga_value_t.hpp"

namespace arith_check {

using hd::ga::value_t;

// the same storage template with another coefficient type
template <typename X> struct rebind;
template <template <typename, typename> class S, typename T, typename Tag>
struct rebind<S<T, Tag>> {
    template <typename U> using to = S<U, Tag>;
};

template <typename X> inline constexpr size_t n_comp = sizeof(X) / sizeof(value_t);

// float-representable coefficients, so the float -> value_t conversion is exact
template <typename C> constexpr C coeff(size_t i, int shift)
{
    return C(0.5 * double(i % 7) - 1.25 + 0.25 * shift + (i % 2 == 0 ? 0.0 : 0.125));
}

// N coefficients of type C (N is the value_t type's count: a float type is half the size)
template <typename X, typename C = value_t, size_t N = n_comp<X>> X make(int shift)
{
    return [&]<size_t... I>(std::index_sequence<I...>) {
        return X(coeff<C>(I, shift)...);
    }(std::make_index_sequence<N>{});
}

template <typename X> struct is_scalar_t : std::false_type {};
template <typename T, typename Tag>
struct is_scalar_t<hd::ga::Scalar_t<T, Tag>> : std::true_type {};

template <typename X> void check(std::string_view name)
{
    INFO(std::string(name));
    X const a = make<X>(0);
    X const b = make<X>(3);
    value_t const s = 2.5;

    {
        X c = a;
        c += b;
        CHECK(c == a + b);
    }
    {
        X c = a;
        c -= b;
        CHECK(c == a - b);
        CHECK(c == a + (-b));
    }
    {
        X c = a;
        c *= s;
        CHECK(c == a * s);
        CHECK(c == s * a);
    }
    {
        X c = a;
        c /= s;
        CHECK(c == a / s);
        CHECK(c == a * (value_t(1) / s));
    }
    CHECK(a + (-a) == X{});
    CHECK(-(-a) == a);
    CHECK(a != b);
    CHECK_FALSE(a != a);

    value_t sum = 0;
    for (size_t i = 0; i < n_comp<X>; ++i)
        sum += coeff<value_t>(i, 0) * coeff<value_t>(i, 0);
    CHECK(hd::ga::detail::coeff_sq(a) == doctest::Approx(sum));

    using Xf = typename rebind<X>::template to<float>;
    Xf const af = make<Xf, float, n_comp<X>>(0);
    if constexpr (is_scalar_t<X>::value)
        CHECK(X(value_t(float(af))) == a); // a Scalar_t converts through its value
    else CHECK(X(af) == a); // every other storage template has a converting ctor

#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    CHECK_THROWS_AS((void)(a / value_t(0)), std::runtime_error);
    X c = a;
    CHECK_THROWS_AS(c /= value_t(0), std::runtime_error);
#endif
}

} // namespace arith_check
