#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// The contract every value-type formatter of the library keeps, asked of one type at a
// time:
//
//   - "{}" prints the type name and every component, in storage order, comma-separated
//   - a numeric spec ("{:.2f}") reaches every component, not the object as a whole
//   - a spec given at RUN time (fmt::runtime) gives the same text: this is the path that
//     executes the formatter's parse() at run time -- a compile-time checked format
//     string runs parse() only during compilation
//   - an invalid spec is rejected with fmt::format_error, not silently ignored
//   - a negative zero prints as "0"
//   - operator<< on a std::ostream, where the type has one, prints the same components
//     without the type name
//
// The components are distinct and exactly representable, so a swapped or a dropped
// component shows as a wrong string, and "{}" prints them identically through fmt and
// through an ostream.

#include <cstddef>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

#include "fmt/format.h"

#include "doctest/doctest.h"
#include "ga/ga_value_t.hpp"

namespace fmt_check {

using hd::ga::value_t;

inline constexpr value_t vals[32] = {1,    -2,   0.5,  3,     -4.25, 5,     6.5,  -7,
                                     8,    9.75, -10,  11,    12.5,  -13,   14,   15.25,
                                     -16,  17,   18.5, -19,   20,    21.25, -22,  23,
                                     24.5, -25,  26,   27.75, -28,   29,    30.5, -31};

// the number of components: every value type stores exactly that many value_t
template <typename X> inline constexpr size_t n_comp = sizeof(X) / sizeof(value_t);

template <typename X, size_t... I> X make_seq(std::index_sequence<I...>)
{
    return X(vals[I]...);
}

template <typename X, size_t... I> X make_negzero_seq(std::index_sequence<I...>)
{
    return X(((void)I, value_t(-0.0))...);
}

// a value with the components vals[0..n), and one with every component -0.0
template <typename X> X make()
{
    return make_seq<X>(std::make_index_sequence<n_comp<X>>{});
}

template <typename X> X make_negzero()
{
    return make_negzero_seq<X>(std::make_index_sequence<n_comp<X>>{});
}

// "v0,v1,...", each component through the given (run-time) spec
inline std::string joined(size_t n, std::string_view spec, bool zeros = false)
{
    std::string s;
    for (size_t i = 0; i < n; ++i) {
        if (i > 0) s += ",";
        s += fmt::format(fmt::runtime(spec), zeros ? value_t(0) : vals[i]);
    }
    return s;
}

template <typename X>
void check_value(X const& x, X const& negzero, std::string_view name)
{
    size_t const n = n_comp<X>;
    std::string const plain = fmt::format("{}({})", name, joined(n, "{}"));
    std::string const fixed = fmt::format("{}({})", name, joined(n, "{:.2f}"));

    // compile-time checked format strings
    CHECK(fmt::format("{}", x) == plain);
    CHECK(fmt::format("{:.2f}", x) == fixed);

    // the same specs at run time
    CHECK(fmt::format(fmt::runtime("{}"), x) == plain);
    CHECK(fmt::format(fmt::runtime("{:.2f}"), x) == fixed);

    // an invalid spec is an error, not ignored
    CHECK_THROWS_AS((void)fmt::format(fmt::runtime("{:q}"), x), fmt::format_error);

    // negative zero is suppressed in every component
    CHECK(fmt::format("{}", negzero) ==
          fmt::format("{}({})", name, joined(n, "{}", true)));

    // iostream: same components, no type name
    if constexpr (requires(std::ostream& os) { os << x; }) {
        std::ostringstream os;
        os << x;
        CHECK(os.str() == "(" + joined(n, "{}") + ")");
    }
}

template <typename X> void check(std::string_view name)
{
    check_value(make<X>(), make_negzero<X>(), name);
}

} // namespace fmt_check
