#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type2dc.hpp" // type aliases for scalars, vector, pseudoscalar, etc.

#include "ga_mvec2dc_e.hpp" // even-grade multivector
#include "ga_mvec2dc_u.hpp" // uneven- or odd-grade multivector


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec16_t including its ctors and add specific ctors for MVec16_t<T, Tag>
// by using partial template specialization for the Tag=mvec2dc_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
struct MVec16_t<T, mvec2dc_tag> : public MVec16_t<T, default_tag> {

    using MVec16_t<T, default_tag>::MVec16_t; // inherit base class ctors

    ////////////////////////////////////////////////////////////////////////////
    // ctors for MVec2dc<T>
    ////////////////////////////////////////////////////////////////////////////

    // assign a scalar part exclusively (other grades = 0)
    constexpr MVec16_t(Scalar2dc<T> s) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a vector part exclusively (other grades = 0)
    constexpr MVec16_t(Vec2dc<T> const& v) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    constexpr MVec16_t(BiVec2dc<T> const& B) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx, B.my,
                 B.mz, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a trivector part exclusively (other grades = 0)
    constexpr MVec16_t(TriVec2dc<T> const& t) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    constexpr MVec16_t(PScalar2dc<T> ps) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even-grade multivector from scalar, bivector, pseudoscalar
    constexpr MVec16_t(Scalar2dc<T> s, BiVec2dc<T> const& B, PScalar2dc<T> ps) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz,
                 T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even-grade multivector
    constexpr MVec16_t(MVec2dc_E<T> const& M) :
        MVec16_t(M.c0, T(0.0), T(0.0), T(0.0), T(0.0), M.c1, M.c2, M.c3, M.c4, M.c5, M.c6,
                 T(0.0), T(0.0), T(0.0), T(0.0), M.c7)
    {
    }

    // assign components of an odd-grade multivector from vector and trivector
    constexpr MVec16_t(Vec2dc<T> const& v, TriVec2dc<T> const& t) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign components of an odd-grade multivector
    constexpr MVec16_t(MVec2dc_U<T> const& M) :
        MVec16_t(T(0.0), M.c0, M.c1, M.c2, M.c3, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), M.c4, M.c5, M.c6, M.c7, T(0.0))
    {
    }

    // assign a full multivector
    // assign components of an even-grade multivector
    constexpr MVec16_t(Scalar2dc<T> s, Vec2dc<T> const& v, BiVec2dc<T> const& B,
                       TriVec2dc<T> const& t, PScalar2dc<T> ps) :
        MVec16_t(T(s), v.x, v.y, v.z, v.w, B.vx, B.vy, B.vz, B.mx, B.my, B.mz, t.x, t.y,
                 t.z, t.w, T(ps))
    {
    }

    // Override compound assignment operators to return correct derived type
    // This ensures GCC+doctest can properly deduce the tag type without needing cross-tag
    // comparisons
    template <typename U>
        requires(numeric_type<U>)
    MVec16_t& operator+=(MVec16_t<U, mvec2dc_tag> const& v) noexcept
    {
        this->c0 += v.c0;
        this->c1 += v.c1;
        this->c2 += v.c2;
        this->c3 += v.c3;
        this->c4 += v.c4;
        this->c5 += v.c5;
        this->c6 += v.c6;
        this->c7 += v.c7;
        this->c8 += v.c8;
        this->c9 += v.c9;
        this->c10 += v.c10;
        this->c11 += v.c11;
        this->c12 += v.c12;
        this->c13 += v.c13;
        this->c14 += v.c14;
        this->c15 += v.c15;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec16_t& operator-=(MVec16_t<U, mvec2dc_tag> const& v) noexcept
    {
        this->c0 -= v.c0;
        this->c1 -= v.c1;
        this->c2 -= v.c2;
        this->c3 -= v.c3;
        this->c4 -= v.c4;
        this->c5 -= v.c5;
        this->c6 -= v.c6;
        this->c7 -= v.c7;
        this->c8 -= v.c8;
        this->c9 -= v.c9;
        this->c10 -= v.c10;
        this->c11 -= v.c11;
        this->c12 -= v.c12;
        this->c13 -= v.c13;
        this->c14 -= v.c14;
        this->c15 -= v.c15;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec16_t& operator*=(U s) noexcept
    {
        this->c0 *= s;
        this->c1 *= s;
        this->c2 *= s;
        this->c3 *= s;
        this->c4 *= s;
        this->c5 *= s;
        this->c6 *= s;
        this->c7 *= s;
        this->c8 *= s;
        this->c9 *= s;
        this->c10 *= s;
        this->c11 *= s;
        this->c12 *= s;
        this->c13 *= s;
        this->c14 *= s;
        this->c15 *= s;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec16_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 16 comp.");
        this->c0 /= s;
        this->c1 /= s;
        this->c2 /= s;
        this->c3 /= s;
        this->c4 /= s;
        this->c5 /= s;
        this->c6 /= s;
        this->c7 /= s;
        this->c8 /= s;
        this->c9 /= s;
        this->c10 /= s;
        this->c11 /= s;
        this->c12 /= s;
        this->c13 /= s;
        this->c14 /= s;
        this->c15 /= s;
        return *this;
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec16_t<T, mvec2dc_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector
// grade 4: gr4() - quadvector (= pseudoscalar in 2dc)

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2dc<T> gr0(MVec2dc<T> const& M)
{
    return Scalar2dc<T>(M.c0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2dc<T> gr1(MVec2dc<T> const& M)
{
    return Vec2dc<T>(M.c1, M.c2, M.c3, M.c4);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec2dc<T> gr2(MVec2dc<T> const& M)
{
    return BiVec2dc<T>(M.c5, M.c6, M.c7, M.c8, M.c9, M.c10);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec2dc<T> gr3(MVec2dc<T> const& M)
{
    return TriVec2dc<T>(M.c11, M.c12, M.c13, M.c14);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2dc<T> gr4(MVec2dc<T> const& M)
{
    return PScalar2dc<T>(M.c15);
}


////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec2dc_e.hpp and mvec2dc_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(Scalar2dc<T> s, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s, v, BiVec2dc<ctype>{}, TriVec2dc<ctype>{},
                          PScalar2dc<ctype>{});
}

// vector + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(Vec2dc<T> const& v, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s, v, BiVec2dc<ctype>{}, TriVec2dc<ctype>{},
                          PScalar2dc<ctype>{});
}

// scalar + trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(Scalar2dc<T> s,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s, Vec2dc<ctype>{}, BiVec2dc<ctype>{}, t, PScalar2dc<ctype>{});
}

// trivector + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(TriVec2dc<T> const& t,
                                                      Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s, Vec2dc<ctype>{}, BiVec2dc<ctype>{}, t, PScalar2dc<ctype>{});
}

// (scalar + pseudoscalar) and (pseudoscalar + scalar) are BOTH even in 2dc (the
// pseudoscalar is grade 4), so they are defined in ga_mvec2dc_e.hpp with the
// narrower MVec2dc_E result.

// vector + bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(Vec2dc<T> const& v,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, v, B, TriVec2dc<ctype>{},
                          PScalar2dc<ctype>{});
}

// bivector + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(BiVec2dc<T> const& B,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, v, B, TriVec2dc<ctype>{},
                          PScalar2dc<ctype>{});
}

// vector + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(Vec2dc<T> const& v,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, v, BiVec2dc<ctype>{}, TriVec2dc<ctype>{},
                          ps);
}

// pseudoscalar + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(PScalar2dc<T> ps,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, v, BiVec2dc<ctype>{}, TriVec2dc<ctype>{},
                          ps);
}

// bivector + trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(BiVec2dc<T> const& B,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, Vec2dc<ctype>{}, B, t, PScalar2dc<ctype>{});
}

// trivector + bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(TriVec2dc<T> const& t,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, Vec2dc<ctype>{}, B, t, PScalar2dc<ctype>{});
}

// pseudoscalar + trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(PScalar2dc<T> ps,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, Vec2dc<ctype>{}, BiVec2dc<ctype>{}, t, ps);
}

// trivector + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(TriVec2dc<T> const& t,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, Vec2dc<ctype>{}, BiVec2dc<ctype>{}, t, ps);
}

// multivector + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc<T> const& M, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M) + s, gr1(M), gr2(M), gr3(M), gr4(M));
}

// multivector + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc<T> const& M,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M) + v, gr2(M), gr3(M), gr4(M));
}

// multivector + bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc<T> const& M,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M), gr2(M) + B, gr3(M), gr4(M));
}

// multivector + trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc<T> const& M,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) + t, gr4(M));
}

// multivector + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc<T> const& M,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M), gr2(M), gr3(M), gr4(M) + ps);
}

////////////////////////////////////////////////////////////////////////////////
// Addition operations completing the type matrix: the remaining pairs whose grade
// union fits no subalgebra, hence returning the full multivector
////////////////////////////////////////////////////////////////////////////////

// scalar + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(Scalar2dc<T> s, MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s + gr0(M), gr1(M), gr2(M), gr3(M), gr4(M));
}

// vector + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(Vec2dc<T> const& v,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), v + gr1(M), gr2(M), gr3(M), gr4(M));
}

// bivector + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(BiVec2dc<T> const& B,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M), B + gr2(M), gr3(M), gr4(M));
}

// trivector + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(TriVec2dc<T> const& t,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M), gr2(M), t + gr3(M), gr4(M));
}

// pseudoscalar + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(PScalar2dc<T> ps,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M), gr2(M), gr3(M), ps + gr4(M));
}

// even grade multivector + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc_E<T> const& E,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E) + gr0(M), gr1(M), gr2(E) + gr2(M), gr3(M),
                          gr4(E) + gr4(M));
}

// odd grade multivector + multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc_U<T> const& O,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(O) + gr1(M), gr2(M), gr3(O) + gr3(M), gr4(M));
}

// multivector + even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc<T> const& M,
                                                      MVec2dc_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M) + gr0(E), gr1(M), gr2(M) + gr2(E), gr3(M),
                          gr4(M) + gr4(E));
}

// multivector + odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc<T> const& M,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M) + gr1(O), gr2(M), gr3(M) + gr3(O), gr4(M));
}

// even grade multivector + odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc_E<T> const& E,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E), gr1(O), gr2(E), gr3(O), gr4(E));
}

// odd grade multivector + even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc_U<T> const& O,
                                                      MVec2dc_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E), gr1(O), gr2(E), gr3(O), gr4(E));
}

// vector + even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(Vec2dc<T> const& v,
                                                      MVec2dc_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E), v, gr2(E), TriVec2dc<ctype>{}, gr4(E));
}

// even grade multivector + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc_E<T> const& E,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E), v, gr2(E), TriVec2dc<ctype>{}, gr4(E));
}

// trivector + even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(TriVec2dc<T> const& t,
                                                      MVec2dc_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E), Vec2dc<ctype>{}, gr2(E), t, gr4(E));
}

// even grade multivector + trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc_E<T> const& E,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E), Vec2dc<ctype>{}, gr2(E), t, gr4(E));
}

// scalar + odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(Scalar2dc<T> s,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s, gr1(O), BiVec2dc<ctype>{}, gr3(O), PScalar2dc<ctype>{});
}

// odd grade multivector + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc_U<T> const& O,
                                                      Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s, gr1(O), BiVec2dc<ctype>{}, gr3(O), PScalar2dc<ctype>{});
}

// bivector + odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(BiVec2dc<T> const& B,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, gr1(O), B, gr3(O), PScalar2dc<ctype>{});
}

// odd grade multivector + bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc_U<T> const& O,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, gr1(O), B, gr3(O), PScalar2dc<ctype>{});
}

// pseudoscalar + odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(PScalar2dc<T> ps,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, gr1(O), BiVec2dc<ctype>{}, gr3(O), ps);
}

// odd grade multivector + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator+(MVec2dc_U<T> const& O,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, gr1(O), BiVec2dc<ctype>{}, gr3(O), ps);
}


////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec2dc_e.hpp and mvec2dc_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(Scalar2dc<T> s, Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s, -v, BiVec2dc<ctype>{}, TriVec2dc<ctype>{},
                          PScalar2dc<ctype>{});
}

// vector - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(Vec2dc<T> const& v, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-s, v, BiVec2dc<ctype>{}, TriVec2dc<ctype>{},
                          PScalar2dc<ctype>{});
}

// scalar - trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(Scalar2dc<T> s,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s, Vec2dc<ctype>{}, BiVec2dc<ctype>{}, -t, PScalar2dc<ctype>{});
}

// trivector - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(TriVec2dc<T> const& t,
                                                      Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-s, Vec2dc<ctype>{}, BiVec2dc<ctype>{}, t, PScalar2dc<ctype>{});
}

// (scalar - pseudoscalar) and (pseudoscalar - scalar) are BOTH even in 2dc (the
// pseudoscalar is grade 4), so they are defined in ga_mvec2dc_e.hpp with the
// narrower MVec2dc_E result.

// vector - bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(Vec2dc<T> const& v,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, v, -B, TriVec2dc<ctype>{},
                          PScalar2dc<ctype>{});
}

// bivector - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(BiVec2dc<T> const& B,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, -v, B, TriVec2dc<ctype>{},
                          PScalar2dc<ctype>{});
}

// vector - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(Vec2dc<T> const& v,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, v, BiVec2dc<ctype>{}, TriVec2dc<ctype>{},
                          -ps);
}

// pseudoscalar - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(PScalar2dc<T> ps,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, -v, BiVec2dc<ctype>{}, TriVec2dc<ctype>{},
                          ps);
}

// bivector - trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(BiVec2dc<T> const& B,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, Vec2dc<ctype>{}, B, -t,
                          PScalar2dc<ctype>{});
}

// trivector - bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(TriVec2dc<T> const& t,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, Vec2dc<ctype>{}, -B, t,
                          PScalar2dc<ctype>{});
}

// pseudoscalar - trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(PScalar2dc<T> ps,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, Vec2dc<ctype>{}, BiVec2dc<ctype>{}, -t, ps);
}

// trivector - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(TriVec2dc<T> const& t,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, Vec2dc<ctype>{}, BiVec2dc<ctype>{}, t, -ps);
}

// multivector - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc<T> const& M, Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M) - s, gr1(M), gr2(M), gr3(M), gr4(M));
}

// multivector - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc<T> const& M,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M) - v, gr2(M), gr3(M), gr4(M));
}

// multivector - bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc<T> const& M,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M), gr2(M) - B, gr3(M), gr4(M));
}

// multivector - trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc<T> const& M,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) - t, gr4(M));
}

// multivector - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc<T> const& M,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M), gr2(M), gr3(M), gr4(M) - ps);
}

////////////////////////////////////////////////////////////////////////////////
// Subtraction operations completing the type matrix: the remaining pairs whose grade
// union fits no subalgebra, hence returning the full multivector
////////////////////////////////////////////////////////////////////////////////

// scalar - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(Scalar2dc<T> s, MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s - gr0(M), -gr1(M), -gr2(M), -gr3(M), -gr4(M));
}

// vector - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(Vec2dc<T> const& v,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-gr0(M), v - gr1(M), -gr2(M), -gr3(M), -gr4(M));
}

// bivector - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(BiVec2dc<T> const& B,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-gr0(M), -gr1(M), B - gr2(M), -gr3(M), -gr4(M));
}

// trivector - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(TriVec2dc<T> const& t,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-gr0(M), -gr1(M), -gr2(M), t - gr3(M), -gr4(M));
}

// pseudoscalar - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(PScalar2dc<T> ps,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-gr0(M), -gr1(M), -gr2(M), -gr3(M), ps - gr4(M));
}

// even grade multivector - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc_E<T> const& E,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E) - gr0(M), -gr1(M), gr2(E) - gr2(M), -gr3(M),
                          gr4(E) - gr4(M));
}

// odd grade multivector - multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc_U<T> const& O,
                                                      MVec2dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-gr0(M), gr1(O) - gr1(M), -gr2(M), gr3(O) - gr3(M), -gr4(M));
}

// multivector - even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc<T> const& M,
                                                      MVec2dc_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M) - gr0(E), gr1(M), gr2(M) - gr2(E), gr3(M),
                          gr4(M) - gr4(E));
}

// multivector - odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc<T> const& M,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(M), gr1(M) - gr1(O), gr2(M), gr3(M) - gr3(O), gr4(M));
}

// even grade multivector - odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc_E<T> const& E,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E), -gr1(O), gr2(E), -gr3(O), gr4(E));
}

// odd grade multivector - even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc_U<T> const& O,
                                                      MVec2dc_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-gr0(E), gr1(O), -gr2(E), gr3(O), -gr4(E));
}

// vector - even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(Vec2dc<T> const& v,
                                                      MVec2dc_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-gr0(E), v, -gr2(E), TriVec2dc<ctype>{}, -gr4(E));
}

// even grade multivector - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc_E<T> const& E,
                                                      Vec2dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E), -v, gr2(E), TriVec2dc<ctype>{}, gr4(E));
}

// trivector - even grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(TriVec2dc<T> const& t,
                                                      MVec2dc_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-gr0(E), Vec2dc<ctype>{}, -gr2(E), t, -gr4(E));
}

// even grade multivector - trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc_E<T> const& E,
                                                      TriVec2dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(gr0(E), Vec2dc<ctype>{}, gr2(E), -t, gr4(E));
}

// scalar - odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(Scalar2dc<T> s,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(s, -gr1(O), BiVec2dc<ctype>{}, -gr3(O), PScalar2dc<ctype>{});
}

// odd grade multivector - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc_U<T> const& O,
                                                      Scalar2dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(-s, gr1(O), BiVec2dc<ctype>{}, gr3(O), PScalar2dc<ctype>{});
}

// bivector - odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(BiVec2dc<T> const& B,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, -gr1(O), B, -gr3(O), PScalar2dc<ctype>{});
}

// odd grade multivector - bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc_U<T> const& O,
                                                      BiVec2dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, gr1(O), -B, gr3(O), PScalar2dc<ctype>{});
}

// pseudoscalar - odd grade multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(PScalar2dc<T> ps,
                                                      MVec2dc_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, -gr1(O), BiVec2dc<ctype>{}, -gr3(O), ps);
}

// odd grade multivector - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dc<std::common_type_t<T, U>> operator-(MVec2dc_U<T> const& O,
                                                      PScalar2dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dc<ctype>(Scalar2dc<ctype>{}, gr1(O), BiVec2dc<ctype>{}, gr3(O), -ps);
}

} // namespace hd::ga