#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type2d.hpp" // type aliases for scalars, vector, pseudoscalar, etc.

#include "ga_mvec2d_e.hpp" // even-grade multivector

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec4_t including its ctors and add specific ctors for MVec4_t<T, Tag>
// by using partial template specialization for the Tag=mvec2d_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
struct MVec4_t<T, mvec2d_tag> : public MVec4_t<T, default_tag> {

    using MVec4_t<T, default_tag>::MVec4_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    constexpr MVec4_t(Scalar2d<T> s) : MVec4_t(T(s), T(0.0), T(0.0), T(0.0)) {}

    // assign a vector part exclusively (other grades = 0)
    constexpr MVec4_t(Vec2d<T> const& v) : MVec4_t(T(0.0), v.x, v.y, T(0.0)) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    constexpr MVec4_t(PScalar2d<T> ps) : MVec4_t(T(0.0), T(0.0), T(0.0), T(ps)) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to full geometric product)
    constexpr MVec4_t(Scalar2d<T> s, PScalar2d<T> ps) :
        MVec4_t(T(s), T(0.0), T(0.0), T(ps))
    {
    }

    // assign from a complex number, i.e. from the even subalgebra
    constexpr MVec4_t(MVec2d_E<T> const& M) : MVec4_t(M.c0, T(0.0), T(0.0), M.c1) {}

    // assign a full multivector
    constexpr MVec4_t(Scalar2d<T> s, Vec2d<T> const& v, PScalar2d<T> ps) :
        MVec4_t(T(s), v.x, v.y, T(ps))
    {
    }

    // Override compound assignment operators to return correct derived type
    // This ensures GCC+doctest can properly deduce the tag type without needing cross-tag
    // comparisons The base class members c0, c1, c2, c3 are directly accessible
    template <typename U>
        requires(numeric_type<U>)
    MVec4_t& operator+=(MVec4_t<U, mvec2d_tag> const& v) noexcept
    {
        this->c0 += v.c0;
        this->c1 += v.c1;
        this->c2 += v.c2;
        this->c3 += v.c3;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec4_t& operator-=(MVec4_t<U, mvec2d_tag> const& v) noexcept
    {
        this->c0 -= v.c0;
        this->c1 -= v.c1;
        this->c2 -= v.c2;
        this->c3 -= v.c3;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec4_t& operator*=(U s) noexcept
    {
        this->c0 *= s;
        this->c1 *= s;
        this->c2 *= s;
        this->c3 *= s;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec4_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 4 comp.");
        this->c0 /= s;
        this->c1 /= s;
        this->c2 /= s;
        this->c3 /= s;
        return *this;
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec4_t<T, mvec2d_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector (= pseudoscalar in 2d)

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2d<T> gr0(MVec2d<T> const& M)
{
    return Scalar2d<T>(M.c0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2d<T> gr1(MVec2d<T> const& M)
{
    return Vec2d<T>(M.c1, M.c2);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2d<T> gr2(MVec2d<T> const& M)
{
    return PScalar2d<T>(M.c3);
}


////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec2d_e.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(s, v, PScalar2d<ctype>{});
}

// vector + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(Vec2d<T> const& v, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(s, v, PScalar2d<ctype>{});
}

// pseudoscalar + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(Scalar2d<ctype>{}, v, ps);
}

// vector + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(Scalar2d<ctype>{}, v, ps);
}

// multivector + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M) + s, gr1(M), gr2(M));
}

// multivector + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(MVec2d<T> const& M,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M), gr1(M) + v, gr2(M));
}

// multivector + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(MVec2d<T> const& M, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M), gr1(M), gr2(M) + ps);
}

////////////////////////////////////////////////////////////////////////////////
// Addition operations completing the type matrix: the remaining pairs whose grade
// union fits no subalgebra, hence returning the full multivector
////////////////////////////////////////////////////////////////////////////////

// scalar + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(s + gr0(M), gr1(M), gr2(M));
}

// vector + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(Vec2d<T> const& v,
                                                     MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M), v + gr1(M), gr2(M));
}

// pseudoscalar + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(PScalar2d<T> ps, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M), gr1(M), ps + gr2(M));
}

// even grade multivector + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(MVec2d_E<T> const& E,
                                                     MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(E) + gr0(M), gr1(M), gr2(E) + gr2(M));
}

// multivector + even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(MVec2d<T> const& M,
                                                     MVec2d_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M) + gr0(E), gr1(M), gr2(M) + gr2(E));
}

// vector + even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(Vec2d<T> const& v,
                                                     MVec2d_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(E), v, gr2(E));
}

// even grade multivector + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator+(MVec2d_E<T> const& E,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(E), v, gr2(E));
}


////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec2d_e.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(s, -v, PScalar2d<ctype>{});
}

// vector - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(Vec2d<T> const& v, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(-s, v, PScalar2d<ctype>{});
}

// pseudoscalar - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(Scalar2d<ctype>{}, -v, ps);
}

// vector - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(Scalar2d<ctype>{}, v, -ps);
}

// multivector - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M) - s, gr1(M), gr2(M));
}

// multivector - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(MVec2d<T> const& M,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M), gr1(M) - v, gr2(M));
}

// multivector - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(MVec2d<T> const& M, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M), gr1(M), gr2(M) - ps);
}

////////////////////////////////////////////////////////////////////////////////
// Subtraction operations completing the type matrix: the remaining pairs whose grade
// union fits no subalgebra, hence returning the full multivector
////////////////////////////////////////////////////////////////////////////////

// scalar - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(s - gr0(M), -gr1(M), -gr2(M));
}

// vector - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(Vec2d<T> const& v,
                                                     MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(-gr0(M), v - gr1(M), -gr2(M));
}

// pseudoscalar - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(PScalar2d<T> ps, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(-gr0(M), -gr1(M), ps - gr2(M));
}

// even grade multivector - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(MVec2d_E<T> const& E,
                                                     MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(E) - gr0(M), -gr1(M), gr2(E) - gr2(M));
}

// multivector - even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(MVec2d<T> const& M,
                                                     MVec2d_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(M) - gr0(E), gr1(M), gr2(M) - gr2(E));
}

// vector - even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(Vec2d<T> const& v,
                                                     MVec2d_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(-gr0(E), v, -gr2(E));
}

// even grade multivector - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator-(MVec2d_E<T> const& E,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(gr0(E), -v, gr2(E));
}

} // namespace hd::ga