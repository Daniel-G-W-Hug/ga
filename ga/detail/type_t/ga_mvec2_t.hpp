#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <array>    // std::array (componentwise comparison)
#include <cmath>    // std::abs, std::sqrt
#include <concepts> // numeric_type<T>
#include <iostream> // std::ostream

#include "../ga_error_handling.hpp"
#include "ga_scalar_t.hpp" // Scalar_t<T, Tag> (graded ctor)
#include "ga_type_tags.hpp"

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T, Tag> definition of a multivector with 2 components c0, c1
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires(numeric_type<T>)
struct MVec2_t {

    // ctors

    // (all grades = 0)
    MVec2_t() = default;

    // assign all components (works with raw floating-point types and implicit
    // conversions)
    constexpr MVec2_t(T s, T ps) : c0(s), c1(ps) {}

    template <typename Tag_S, typename Tag_PS>
        requires(numeric_type<T> && same_base_class<Tag_S, Tag_PS> &&
                 !std::same_as<Tag_S, Tag_PS>)
    constexpr MVec2_t(Scalar_t<T, Tag_S> const& s, Scalar_t<T, Tag_PS> const& ps) :
        c0(T(s)), c1(T(ps))
    {
    }

    // cp & mv ctor
    MVec2_t(MVec2_t const&) = default;
    MVec2_t(MVec2_t&&) noexcept = default;

    // cp & mv assignment
    MVec2_t& operator=(MVec2_t const&) = default;
    MVec2_t& operator=(MVec2_t&&) noexcept = default;

    // dtor
    ~MVec2_t() = default;

    // floating point type conversion
    template <typename U>
        requires(numeric_type<U>)
    constexpr MVec2_t(MVec2_t<U, Tag> const& v) : MVec2_t(v.c0, v.c1)
    {
    }

    friend void swap(MVec2_t& lhs, MVec2_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.c0, rhs.c0);
        swap(lhs.c1, rhs.c1);
    }

    ////////////////////////////////////////////////////////////////////////////
    // component definition
    ////////////////////////////////////////////////////////////////////////////

    T c0{}; // scalar component
    T c1{}; // bivector component (2d pseudoscalar)

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator+=(MVec2_t<U, Tag> const& v) noexcept
    {
        c0 += v.c0;
        c1 += v.c1;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator-=(MVec2_t<U, Tag> const& v) noexcept
    {
        c0 -= v.c0;
        c1 -= v.c1;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator*=(U s) noexcept
    {
        c0 *= s;
        c1 *= s;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 2 comp.");
        c0 /= s;
        c1 /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// equality - allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator==(MVec2_t<T, Tag> const& lhs, MVec2_t<U, Tag> const& rhs)
{
    // componentwise comparison against the library's one comparison rule
    // (detail::coeffs_equal, in ga_error_handling.hpp -- see the note there)
    return detail::coeffs_equal(std::array{lhs.c0, lhs.c1}, std::array{rhs.c0, rhs.c1});
}

// inequality - allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator!=(MVec2_t<T, Tag> const& lhs, MVec2_t<U, Tag> const& rhs)
{
    return !(lhs == rhs);
}

// unary minus for multivectors from the even subalgebra
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr MVec2_t<T, Tag> operator-(MVec2_t<T, Tag> const& v)
{
    return MVec2_t<T, Tag>(-v.c0, -v.c1);
}

// add multivectors from the even subalgebra
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2_t<std::common_type_t<T, U>, Tag> operator+(MVec2_t<T, Tag> const& v1,
                                                           MVec2_t<U, Tag> const& v2)
{
    return MVec2_t<std::common_type_t<T, U>, Tag>(v1.c0 + v2.c0, v1.c1 + v2.c1);
}

// substract multivectors multivectors from the even subalgebra
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2_t<std::common_type_t<T, U>, Tag> operator-(MVec2_t<T, Tag> const& v1,
                                                           MVec2_t<U, Tag> const& v2)
{
    return MVec2_t<std::common_type_t<T, U>, Tag>(v1.c0 - v2.c0, v1.c1 - v2.c1);
}


// multiply a multivector multivectors from the even subalgebra with a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2_t<std::common_type_t<T, U>, Tag> operator*(MVec2_t<T, Tag> const& v, U s)
{
    return MVec2_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s);
}

template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2_t<std::common_type_t<T, U>, Tag> operator*(T s, MVec2_t<U, Tag> const& v)
{
    return MVec2_t<std::common_type_t<T, U>, Tag>(v.c0 * s, v.c1 * s);
}

// devide an even multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
inline MVec2_t<std::common_type_t<T, U>, Tag> operator/(MVec2_t<T, Tag> const& v, U s)
{
    detail::check_division_by_zero<T, U>(s, "multivector division");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return MVec2_t<ctype, Tag>(v.c0 * inv, v.c1 * inv);
}

////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(numeric_type<T>)
std::ostream& operator<<(std::ostream& os, MVec2_t<T, Tag> const& v)
{
    os << "(" << v.c0 << "," << v.c1 << ")";
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
constexpr T coeff_sq(MVec2_t<T, Tag> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1;
}

} // namespace detail

} // namespace hd::ga