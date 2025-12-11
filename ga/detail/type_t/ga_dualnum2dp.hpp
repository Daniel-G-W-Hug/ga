#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type2dp.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// DualNum2dp<T> M = c0 1 + c1 I_2dp (dual number with s 1 + ps I_2dp,
//                                    where I_2dp = e321, e321^2 = 0)
/////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec2_t including its ctors and add specific ctors for MVec2_t<T, Tag>
// by using partial template specialization for the Tag=dual_number2dp_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct MVec2_t<T, dual_number2dp_tag> : public MVec2_t<T, default_tag> {

    using MVec2_t<T, default_tag>::MVec2_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec2_t(Scalar2dp<T> s) : MVec2_t(T(s), T(0.0)) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec2_t(PScalar2dp<T> ps) : MVec2_t(T(0.0), T(ps)) {}

    // assign from a scalar and a pseudoscalar part
    MVec2_t(Scalar2dp<T> s, PScalar2dp<T> ps) : MVec2_t(T(s), T(ps)) {}

    // Override compound assignment operators to return correct derived type
    // This ensures GCC+doctest can properly deduce the tag type without needing cross-tag
    // comparisons
    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator+=(MVec2_t<U, dual_number2dp_tag> const& v) noexcept
    {
        this->c0 += v.c0;
        this->c1 += v.c1;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator-=(MVec2_t<U, dual_number2dp_tag> const& v) noexcept
    {
        this->c0 -= v.c0;
        this->c1 -= v.c1;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator*=(U s) noexcept
    {
        this->c0 *= s;
        this->c1 *= s;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 2 comp.");
        this->c0 /= s;
        this->c1 /= s;
        return *this;
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization
// DualNum2dp<T> = MVec2_t<T, dual_number2dp_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector from the even subalgebra
//
// grade 0: gr0() - scalar
// grade 3: gr3() - pseudoscalar in pga2dp

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dp<T> gr0(DualNum2dp<T> const& M)
{
    return Scalar2dp<T>(M.c0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dp<T> gr3(DualNum2dp<T> const& M)
{
    return PScalar2dp<T>(M.c1);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine dual numbers
////////////////////////////////////////////////////////////////////////////////

// scalar + pseudoscalar => dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator+(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(s, ps);
}

// pseudoscalar + scalar => dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator+(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(s, ps);
}

// scalar + dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator+(Scalar2dp<T> s,
                                                         DualNum2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(T(s) + M.c0, M.c1);
}

// dual number + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator+(DualNum2dp<T> const& M,
                                                         Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(U(s) + M.c0, M.c1);
}

// pseudosscalar + dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator+(PScalar2dp<T> ps,
                                                         DualNum2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(M.c0, T(ps) + M.c1);
}

// dual number + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator+(DualNum2dp<T> const& M,
                                                         PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(M.c0, U(ps) + M.c1);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars and dual numbers
////////////////////////////////////////////////////////////////////////////////

// scalar - pseudoscalar => dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator-(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(s, -ps);
}

// pseudoscalar - scalar => dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator-(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(-s, ps);
}

// scalar - dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator-(Scalar2dp<T> s,
                                                         DualNum2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(T(s) - M.c0, -M.c1);
}

// dual number - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator-(DualNum2dp<T> const& M,
                                                         Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(M.c0 - U(s), M.c1);
}

// pseudosscalar - dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator-(PScalar2dp<T> ps,
                                                         DualNum2dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(-M.c0, T(ps) - M.c1);
}

// dual number - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum2dp<std::common_type_t<T, U>> operator-(DualNum2dp<T> const& M,
                                                         PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum2dp<ctype>(M.c0, M.c1 - U(ps));
}

} // namespace hd::ga