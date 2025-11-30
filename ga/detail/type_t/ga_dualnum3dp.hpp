#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type3dp.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// DualNum3dp<T> M = c0 1 + c1 I_3dp (dual number with s 1 + ps I_3dp,
//                                    where I_3dp = e1234, e1234^2 = 0)
/////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec2_t including its ctors and add specific ctors for MVec2_t<T, Tag>
// by using partial template specialization for the Tag=dual_number3dp_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct MVec2_t<T, dual_number3dp_tag> : public MVec2_t<T, default_tag> {

    using MVec2_t<T, default_tag>::MVec2_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec2_t(Scalar3dp<T> s) : MVec2_t(T(s), T(0.0)) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec2_t(PScalar3dp<T> ps) : MVec2_t(T(0.0), T(ps)) {}

    // assign from a scalar and a pseudoscalar part
    MVec2_t(Scalar3dp<T> s, PScalar3dp<T> ps) : MVec2_t(T(s), T(ps)) {}

    // Override compound assignment operators to return correct derived type
    // This ensures GCC+doctest can properly deduce the tag type without needing cross-tag
    // comparisons
    template <typename U>
        requires(std::floating_point<U>)
    MVec2_t& operator+=(MVec2_t<U, dual_number3dp_tag> const& v) noexcept
    {
        this->c0 += v.c0;
        this->c1 += v.c1;
        return *this;
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec2_t& operator-=(MVec2_t<U, dual_number3dp_tag> const& v) noexcept
    {
        this->c0 -= v.c0;
        this->c1 -= v.c1;
        return *this;
    }

    template <typename U>
        requires(std::floating_point<U>)
    MVec2_t& operator*=(U s) noexcept
    {
        this->c0 *= s;
        this->c1 *= s;
        return *this;
    }

    template <typename U>
        requires(std::floating_point<U>)
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
// DualNum3dp<T> = MVec2_t<T, dual_number3dp_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector from the even subalgebra
//
// grade 0: gr0() - scalar
// grade 4: gr4() - pseudoscalar in pga3dp

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> gr0(DualNum3dp<T> const& M)
{
    return Scalar3dp<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> gr4(DualNum3dp<T> const& M)
{
    return PScalar3dp<T>(M.c1);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine dual numbers
////////////////////////////////////////////////////////////////////////////////

// scalar + pseudoscalar => dual number
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator+(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(s, ps);
}

// pseudoscalar + scalar => dual number
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator+(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(s, ps);
}

// scalar + dual number
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator+(Scalar3dp<T> s,
                                                         DualNum3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(T(s) + M.c0, M.c1);
}

// dual number + scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator+(DualNum3dp<T> const& M,
                                                         Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(U(s) + M.c0, M.c1);
}

// pseudosscalar + dual number
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator+(PScalar3dp<T> ps,
                                                         DualNum3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(M.c0, T(ps) + M.c1);
}

// dual number + pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator+(DualNum3dp<T> const& M,
                                                         PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(M.c0, U(ps) + M.c1);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars and dual numbers
////////////////////////////////////////////////////////////////////////////////

// scalar - pseudoscalar => dual number
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator-(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(s, -ps);
}

// pseudoscalar - scalar => dual number
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator-(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(-s, ps);
}

// scalar - dual number
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator-(Scalar3dp<T> s,
                                                         DualNum3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(T(s) - M.c0, -M.c1);
}

// dual number - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator-(DualNum3dp<T> const& M,
                                                         Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(M.c0 - U(s), M.c1);
}

// pseudosscalar - dual number
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator-(PScalar3dp<T> ps,
                                                         DualNum3dp<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(-M.c0, T(ps) - M.c1);
}

// dual number - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr DualNum3dp<std::common_type_t<T, U>> operator-(DualNum3dp<T> const& M,
                                                         PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum3dp<ctype>(M.c0, M.c1 - U(ps));
}

} // namespace hd::ga