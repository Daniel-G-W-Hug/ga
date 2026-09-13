#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <array>    // std::array (componentwise comparison)
#include <concepts> // numeric_type<T>
#include <ostream>  // std::ostream

#include "../ga_error_handling.hpp"
#include "ga_type_tags.hpp"

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// MVec8_t<T, Tag> definition of a multivector with 8 components c0, ..., c7
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(numeric_type<T>)
struct MVec8_t {

    // ctors

    // (all grades = 0)
    MVec8_t() = default;

    // assign all components directly
    constexpr MVec8_t(T s, T x, T y, T z, T yz, T zx, T xy, T ps) :
        c0(s), c1(x), c2(y), c3(z), c4(yz), c5(zx), c6(xy), c7(ps)
    {
    }

    // cp & mv ctor
    MVec8_t(MVec8_t const&) = default;
    MVec8_t(MVec8_t&&) noexcept = default;

    // cp & mv assignment
    MVec8_t& operator=(MVec8_t const&) = default;
    MVec8_t& operator=(MVec8_t&&) noexcept = default;

    // dtor
    ~MVec8_t() = default;

    // floating point type conversion
    template <typename U>
        requires(numeric_type<U>)
    constexpr MVec8_t(MVec8_t<U, Tag> const& v) :
        MVec8_t(v.c0, v.c1, v.c2, v.c3, v.c4, v.c5, v.c6, v.c7)
    {
    }

    friend void swap(MVec8_t& lhs, MVec8_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.c0, rhs.c0);
        swap(lhs.c1, rhs.c1);
        swap(lhs.c2, rhs.c2);
        swap(lhs.c3, rhs.c3);
        swap(lhs.c4, rhs.c4);
        swap(lhs.c5, rhs.c5);
        swap(lhs.c6, rhs.c6);
        swap(lhs.c7, rhs.c7);
    }

    ////////////////////////////////////////////////////////////////////////////
    // component definition
    ////////////////////////////////////////////////////////////////////////////

    T c0{}; // scalar
    T c1{}; // vector 3d, 1st component   (x)  - maps to basis vector    e1
    T c2{}; // vector 3d, 2nd component   (y)  - maps to basis vector    e2
    T c3{}; // vector 3d, 3rd component   (z)  - maps to basis vector    e3
    T c4{}; // bivector 3d, 1st component (yz) - maps to basis bivector  e2^e3
    T c5{}; // bivector 3d, 2nd component (zx) - maps to basis bivector  e3^e1
    T c6{}; // bivector 3d, 3rd component (xy) - maps to basis bivector  e1^e2
    T c7{}; // trivector 3d = 3d pseudoscalar  - maps to basis trivector e1^e2^e3

    template <typename U>
        requires(numeric_type<U>)
    MVec8_t& operator+=(MVec8_t<U, Tag> const& v) noexcept
    {
        c0 += v.c0;
        c1 += v.c1;
        c2 += v.c2;
        c3 += v.c3;
        c4 += v.c4;
        c5 += v.c5;
        c6 += v.c6;
        c7 += v.c7;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec8_t& operator-=(MVec8_t<U, Tag> const& v) noexcept
    {
        c0 -= v.c0;
        c1 -= v.c1;
        c2 -= v.c2;
        c3 -= v.c3;
        c4 -= v.c4;
        c5 -= v.c5;
        c6 -= v.c6;
        c7 -= v.c7;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec8_t& operator*=(U s) noexcept
    {
        c0 *= s;
        c1 *= s;
        c2 *= s;
        c3 *= s;
        c4 *= s;
        c5 *= s;
        c6 *= s;
        c7 *= s;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec8_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 8 comp.");
        c0 /= s;
        c1 /= s;
        c2 /= s;
        c3 /= s;
        c4 /= s;
        c5 /= s;
        c6 /= s;
        c7 /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec8_t<T, Tag> core operations
////////////////////////////////////////////////////////////////////////////////

// equality - allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator==(MVec8_t<T, Tag> const& lhs, MVec8_t<U, Tag> const& rhs)
{
    // componentwise comparison against the library's one comparison rule
    // (detail::coeffs_equal, in ga_error_handling.hpp -- see the note there)
    return detail::coeffs_equal(
        std::array{lhs.c0, lhs.c1, lhs.c2, lhs.c3, lhs.c4, lhs.c5, lhs.c6, lhs.c7},
        std::array{rhs.c0, rhs.c1, rhs.c2, rhs.c3, rhs.c4, rhs.c5, rhs.c6, rhs.c7});
}

// inequality - allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator!=(MVec8_t<T, Tag> const& lhs, MVec8_t<U, Tag> const& rhs)
{
    return !(lhs == rhs);
}

// unary minus
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr MVec8_t<T, Tag> operator-(MVec8_t<T, Tag> const& v)
{
    return MVec8_t<T, Tag>(-v.c0, -v.c1, -v.c2, -v.c3, -v.c4, -v.c5, -v.c6, -v.c7);
}

// add multivectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec8_t<std::common_type_t<T, U>, Tag> operator+(MVec8_t<T, Tag> const& v1,
                                                           MVec8_t<U, Tag> const& v2)
{
    return MVec8_t<std::common_type_t<T, U>, Tag>(
        v1.c0 + v2.c0, v1.c1 + v2.c1, v1.c2 + v2.c2, v1.c3 + v2.c3, v1.c4 + v2.c4,
        v1.c5 + v2.c5, v1.c6 + v2.c6, v1.c7 + v2.c7);
}

// substract multivectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec8_t<std::common_type_t<T, U>, Tag> operator-(MVec8_t<T, Tag> const& v1,
                                                           MVec8_t<U, Tag> const& v2)
{
    return MVec8_t<std::common_type_t<T, U>, Tag>(
        v1.c0 - v2.c0, v1.c1 - v2.c1, v1.c2 - v2.c2, v1.c3 - v2.c3, v1.c4 - v2.c4,
        v1.c5 - v2.c5, v1.c6 - v2.c6, v1.c7 - v2.c7);
}

// multiply a multivector with a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec8_t<std::common_type_t<T, U>, Tag> operator*(MVec8_t<T, Tag> const& v, U s)
{
    return MVec8_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s,
                                                  v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s);
}

template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec8_t<std::common_type_t<T, U>, Tag> operator*(T s, MVec8_t<U, Tag> const& v)
{
    return MVec8_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s, v.c2 * s, v.c3 * s,
                                                  v.c4 * s, v.c5 * s, v.c6 * s, v.c7 * s);
}

// devide a multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
inline MVec8_t<std::common_type_t<T, U>, Tag> operator/(MVec8_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "multivector division");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec8_t<ctype, Tag>(v.c0 * inv, v.c1 * inv, v.c2 * inv, v.c3 * inv, v.c4 * inv,
                               v.c5 * inv, v.c6 * inv, v.c7 * inv);
}

////////////////////////////////////////////////////////////////////////////////
// MVec8_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(numeric_type<T>)
std::ostream& operator<<(std::ostream& os, MVec8_t<T, Tag> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << "," << v.c2 << "," << v.c3 << "," << v.c4 << ","
       << v.c5 << "," << v.c6 << "," << v.c7 << ")";
    return os;
}


////////////////////////////////////////////////////////////////////////////////
// coeff_sq: the sum of the squared coefficients -- a METRIC-FREE gauge of "how big is
// this object", used to judge degeneracy relative to the object itself (see
// detail::check_invertible in ga_error_handling.hpp). It is NOT a norm: no metric enters,
// which is the point -- nrm_sq vanishes on a null blade in a Lorentzian or conformal
// metric and on an ideal one in a degenerate metric, exactly where a gauge is needed.
////////////////////////////////////////////////////////////////////////////////

namespace detail {

template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T coeff_sq(MVec8_t<T, Tag> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1 + v.c2 * v.c2 + v.c3 * v.c3 + v.c4 * v.c4 +
           v.c5 * v.c5 + v.c6 * v.c6 + v.c7 * v.c7;
}

} // namespace detail

} // namespace hd::ga