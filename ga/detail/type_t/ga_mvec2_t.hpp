#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <cmath>    // std::abs, std::sqrt
#include <concepts> // numeric_type<T>
#include <iostream> // std::ostream

#include "../ga_error_handling.hpp"
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
    // componentwise comparison
    // equality implies same magnitude and direction
    // comparison is not exact, but accepts epsilon deviations
    auto abs_delta_c0 = std::abs(lhs.c0 - rhs.c0);
    auto abs_delta_c1 = std::abs(lhs.c1 - rhs.c1);
    auto constexpr delta_eps = detail::safe_epsilon<T, U>();
    return (abs_delta_c0 < delta_eps && abs_delta_c1 < delta_eps);
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

// magnitude of a k-vector
// return squared magnitude (e.g. of complex number)
// |Z|^2 = Z rev(Z) = c0^2 + c1^2
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr T nrm_sq(MVec2_t<T, Tag> const& v)
{
    return v.c0 * v.c0 + v.c1 * v.c1;
}

template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr T nrm(MVec2_t<T, Tag> const& v)
{
    return sqrt(nrm_sq(v));
}

// return a multivector M normalized to nrm(M) == 1.0
template <typename T, typename Tag>
    requires(numeric_type<T>)
inline MVec2_t<T, Tag> normalize(MVec2_t<T, Tag> const& M)
{
    T m = nrm(M);
    detail::check_normalization<T>(m, "multivector");
    T inv = T(1.0) / m; // for multiplication with inverse of norm
    return MVec2_t<T, Tag>(M.c0 * inv, M.c1 * inv);
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

} // namespace hd::ga