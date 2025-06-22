#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type2d.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// MVec2d_E<T> M = c0 + c1 I (with I being the bivector of the plane e1^e2)
//
// This is the definition of a multivector in the even subalgebra of G<2,0,0>,
// i.e. it models only multivectors with even grades 0 and 2 in the plane e1^e2
// (e.g. complex numbers).
// This subalgebra is closed under addition and multiplication.
/////////////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of a fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec2_t including its ctors and add specific ctors for MVec2_t<T, Tag>
// by using partial template specialization for the Tag=mvec2d_e_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec2_t<T, mvec2d_e_tag> : public MVec2_t<T, default_tag> {

    using MVec2_t<T, default_tag>::MVec2_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec2_t(Scalar2d<T> s) : MVec2_t(s, T(0.0)) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec2_t(PScalar2d<T> ps) : MVec2_t(T(0.0), ps) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to full geometric product)
    MVec2_t(Scalar2d<T> s, PScalar2d<T> ps) : MVec2_t(T(s), T(ps)) {}
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec2_t<T, mvec2d_e_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector from the even subalgebra
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector (= pseudoscalar in 2d)

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> gr0(MVec2d_E<T> const& M)
{
    return Scalar2d<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> gr2(MVec2d_E<T> const& M)
{
    return PScalar2d<T>(M.c1);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar + pseudoscalar (=bivector) => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator+(Scalar2d<T> s,
                                                              PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(s, ps);
}

// pseudoscalar (=bivector) + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator+(PScalar2d<T> ps,
                                                              Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(s, ps);
}

// scalar + even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator+(Scalar2d<T> s,
                                                              MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(s + M.c0, M.c1);
}

// even grade multivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator+(MVec2d_E<T> const& M,
                                                              Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(s + M.c0, M.c1);
}

// pseudosscalar + even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator+(PScalar2d<T> ps,
                                                              MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(M.c0, ps + M.c1);
}

// even grade multivector + pseudoscalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator+(MVec2d_E<T> const& M,
                                                              PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(M.c0, ps + M.c1);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar - pseudoscalar (=bivector) => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator-(Scalar2d<T> s,
                                                              PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(s, -ps);
}

// pseudoscalar (=bivector) - scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator-(PScalar2d<T> ps,
                                                              Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(-s, ps);
}

// scalar - even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator-(Scalar2d<T> s,
                                                              MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(s - M.c0, -M.c1);
}

// even grade multivector - scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator-(MVec2d_E<T> const& M,
                                                              Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(M.c0 - s, M.c1);
}

// pseudosscalar (=bivector) - even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator-(PScalar2d<T> ps,
                                                              MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(-M.c0, ps - M.c1);
}

// even grade multivector - pseudoscalar (=bivector) => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2d_E<std::common_type_t<T, U>> operator-(MVec2d_E<T> const& M,
                                                              PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(M.c0, M.c1 - ps);
}

} // namespace hd::ga