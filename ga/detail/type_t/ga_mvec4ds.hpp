#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type4ds.hpp" // type aliases for scalars, vector, pseudoscalar, etc.

#include "ga_mvec4ds_e.hpp" // even grade multivector
#include "ga_mvec4ds_u.hpp" // uneven (or odd) grade multivector


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec16_t including its ctors and add specific ctors for MVec16_t<T, Tag>
// by using partial template specialization for the Tag=mvec4ds_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
struct MVec16_t<T, mvec4ds_tag> : public MVec16_t<T, default_tag> {

    using MVec16_t<T, default_tag>::MVec16_t; // inherit base class ctors

    ////////////////////////////////////////////////////////////////////////////
    // ctors for MVec4ds<T>
    ////////////////////////////////////////////////////////////////////////////

    // assign a scalar part exclusively (other grades = 0)
    constexpr MVec16_t(Scalar4ds<T> s) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a vector part exclusively (other grades = 0)
    constexpr MVec16_t(Vec4ds<T> const& v) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    constexpr MVec16_t(BiVec4ds<T> const& B) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx, B.my,
                 B.mz, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a trivector part exclusively (other grades = 0)
    constexpr MVec16_t(TriVec4ds<T> const& t) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    constexpr MVec16_t(PScalar4ds<T> ps) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even grade multivector
    constexpr MVec16_t(Scalar4ds<T> s, BiVec4ds<T> const& B, PScalar4ds<T> ps) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz,
                 T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign components of an even grade multivector
    constexpr MVec16_t(MVec4ds_E<T> const& M) :
        MVec16_t(M.c0, T(0.0), T(0.0), T(0.0), T(0.0), M.c1, M.c2, M.c3, M.c4, M.c5, M.c6,
                 T(0.0), T(0.0), T(0.0), T(0.0), M.c7)
    {
    }

    // assign components of an odd grade multivector
    constexpr MVec16_t(Vec4ds<T> const& v, TriVec4ds<T> const& t) :
        MVec16_t(T(0.0), v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), t.x, t.y, t.z, t.w, T(0.0))
    {
    }

    // assign components of an odd grade multivector
    constexpr MVec16_t(MVec4ds_U<T> const& M) :
        MVec16_t(T(0.0), M.c0, M.c1, M.c2, M.c3, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), M.c4, M.c5, M.c6, M.c7, T(0.0))
    {
    }

    // assign a full multivector
    // assign components of an even grade multivector
    constexpr MVec16_t(Scalar4ds<T> s, Vec4ds<T> const& v, BiVec4ds<T> const& B,
                       TriVec4ds<T> const& t, PScalar4ds<T> ps) :
        MVec16_t(T(s), v.x, v.y, v.z, v.w, B.vx, B.vy, B.vz, B.mx, B.my, B.mz, t.x, t.y,
                 t.z, t.w, T(ps))
    {
    }

    // Override compound assignment operators to return correct derived type
    // This ensures GCC+doctest can properly deduce the tag type without needing cross-tag
    // comparisons
    template <typename U>
        requires(numeric_type<U>)
    MVec16_t& operator+=(MVec16_t<U, mvec4ds_tag> const& v) noexcept
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
    MVec16_t& operator-=(MVec16_t<U, mvec4ds_tag> const& v) noexcept
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
// define grade operations for partial specialization MVec16_t<T, mvec4ds_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector
// grade 4: gr4() - quadvector (= pseudoscalar in 4ds)

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> gr0(MVec4ds<T> const& M)
{
    return Scalar4ds<T>(M.c0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> gr1(MVec4ds<T> const& M)
{
    return Vec4ds<T>(M.c1, M.c2, M.c3, M.c4);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> gr2(MVec4ds<T> const& M)
{
    return BiVec4ds<T>(M.c5, M.c6, M.c7, M.c8, M.c9, M.c10);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> gr3(MVec4ds<T> const& M)
{
    return TriVec4ds<T>(M.c11, M.c12, M.c13, M.c14);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> gr4(MVec4ds<T> const& M)
{
    return PScalar4ds<T>(M.c15);
}


////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec4ds_e.hpp and mvec4ds_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(Scalar4ds<T> s, Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(s, v, BiVec4ds<ctype>{}, TriVec4ds<ctype>{},
                          PScalar4ds<ctype>{});
}

// vector + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(Vec4ds<T> const& v, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(s, v, BiVec4ds<ctype>{}, TriVec4ds<ctype>{},
                          PScalar4ds<ctype>{});
}

// scalar + trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(Scalar4ds<T> s,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(s, Vec4ds<ctype>{}, BiVec4ds<ctype>{}, t, PScalar4ds<ctype>{});
}

// trivector + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(TriVec4ds<T> const& t,
                                                      Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(s, Vec4ds<ctype>{}, BiVec4ds<ctype>{}, t, PScalar4ds<ctype>{});
}

// (scalar + pseudoscalar) and (pseudoscalar + scalar) intentionally NOT
// defined here — see ga_dualnum4ds.hpp for the canonical DualNum4ds result.

// vector + bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(Vec4ds<T> const& v,
                                                      BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, v, B, TriVec4ds<ctype>{},
                          PScalar4ds<ctype>{});
}

// bivector + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(BiVec4ds<T> const& B,
                                                      Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, v, B, TriVec4ds<ctype>{},
                          PScalar4ds<ctype>{});
}

// vector + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(Vec4ds<T> const& v,
                                                      PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, v, BiVec4ds<ctype>{}, TriVec4ds<ctype>{},
                          ps);
}

// pseudoscalar + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(PScalar4ds<T> ps,
                                                      Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, v, BiVec4ds<ctype>{}, TriVec4ds<ctype>{},
                          ps);
}

// bivector + trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(BiVec4ds<T> const& B,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, Vec4ds<ctype>{}, B, t, PScalar4ds<ctype>{});
}

// trivector + bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(TriVec4ds<T> const& t,
                                                      BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, Vec4ds<ctype>{}, B, t, PScalar4ds<ctype>{});
}

// pseudoscalar + trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(PScalar4ds<T> ps,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, Vec4ds<ctype>{}, BiVec4ds<ctype>{}, t, ps);
}

// trivector + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(TriVec4ds<T> const& t,
                                                      PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, Vec4ds<ctype>{}, BiVec4ds<ctype>{}, t, ps);
}

// multivector + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(MVec4ds<T> const& M, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M) + s, gr1(M), gr2(M), gr3(M), gr4(M));
}

// multivector + vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(MVec4ds<T> const& M,
                                                      Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M), gr1(M) + v, gr2(M), gr3(M), gr4(M));
}

// multivector + bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(MVec4ds<T> const& M,
                                                      BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M), gr1(M), gr2(M) + B, gr3(M), gr4(M));
}

// multivector + trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(MVec4ds<T> const& M,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) + t, gr4(M));
}

// multivector + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(MVec4ds<T> const& M,
                                                      PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M), gr1(M), gr2(M), gr3(M), gr4(M) + ps);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec4ds_e.hpp and mvec4ds_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(Scalar4ds<T> s, Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(s, -v, BiVec4ds<ctype>{}, TriVec4ds<ctype>{},
                          PScalar4ds<ctype>{});
}

// vector - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(Vec4ds<T> const& v, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(-s, v, BiVec4ds<ctype>{}, TriVec4ds<ctype>{},
                          PScalar4ds<ctype>{});
}

// scalar - trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(Scalar4ds<T> s,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(s, Vec4ds<ctype>{}, BiVec4ds<ctype>{}, -t, PScalar4ds<ctype>{});
}

// trivector - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(TriVec4ds<T> const& t,
                                                      Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(-s, Vec4ds<ctype>{}, BiVec4ds<ctype>{}, t, PScalar4ds<ctype>{});
}

// (scalar - pseudoscalar) and (pseudoscalar - scalar) intentionally NOT
// defined here — see ga_dualnum4ds.hpp for the canonical DualNum4ds result.

// vector - bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(Vec4ds<T> const& v,
                                                      BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, v, -B, TriVec4ds<ctype>{},
                          PScalar4ds<ctype>{});
}

// bivector - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(BiVec4ds<T> const& B,
                                                      Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, -v, B, TriVec4ds<ctype>{},
                          PScalar4ds<ctype>{});
}

// vector - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(Vec4ds<T> const& v,
                                                      PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, v, BiVec4ds<ctype>{}, TriVec4ds<ctype>{},
                          -ps);
}

// pseudoscalar - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(PScalar4ds<T> ps,
                                                      Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, -v, BiVec4ds<ctype>{}, TriVec4ds<ctype>{},
                          ps);
}

// bivector - trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(BiVec4ds<T> const& B,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, Vec4ds<ctype>{}, B, -t,
                          PScalar4ds<ctype>{});
}

// trivector - bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(TriVec4ds<T> const& t,
                                                      BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, Vec4ds<ctype>{}, -B, t,
                          PScalar4ds<ctype>{});
}

// pseudoscalar - trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(PScalar4ds<T> ps,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, Vec4ds<ctype>{}, BiVec4ds<ctype>{}, -t, ps);
}

// trivector - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(TriVec4ds<T> const& t,
                                                      PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>{}, Vec4ds<ctype>{}, BiVec4ds<ctype>{}, t, -ps);
}

// multivector - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(MVec4ds<T> const& M, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M) - s, gr1(M), gr2(M), gr3(M), gr4(M));
}

// multivector - vector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(MVec4ds<T> const& M,
                                                      Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M), gr1(M) - v, gr2(M), gr3(M), gr4(M));
}

// multivector - bivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(MVec4ds<T> const& M,
                                                      BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M), gr1(M), gr2(M) - B, gr3(M), gr4(M));
}

// multivector - trivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(MVec4ds<T> const& M,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) - t, gr4(M));
}

// multivector - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(MVec4ds<T> const& M,
                                                      PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(gr0(M), gr1(M), gr2(M), gr3(M), gr4(M) - ps);
}

} // namespace hd::ga