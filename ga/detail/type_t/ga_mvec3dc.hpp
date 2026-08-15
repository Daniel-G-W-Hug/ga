#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_mvec3dc_e.hpp"
#include "ga_mvec3dc_u.hpp"
#include "ga_type3dc.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec32_t including its ctors and add specific ctors for MVec32_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dc_tag
// (cga3dc layout: c0 = scalar, c1..c5 = vector, c6..c15 = bivector,
//  c16..c25 = trivector, c26..c30 = quadvector, c31 = pseudoscalar)
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
struct MVec32_t<T, mvec3dc_tag> : public MVec32_t<T, default_tag> {

    using MVec32_t<T, default_tag>::MVec32_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    constexpr MVec32_t(Scalar3dc<T> s) :
        MVec32_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a vector part exclusively (other grades = 0)
    constexpr MVec32_t(Vec3dc<T> const& v) :
        MVec32_t(T(0.0), v.x, v.y, v.z, v.w, v.u, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    constexpr MVec32_t(BiVec3dc<T> const& B) :
        MVec32_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx,
                 B.my, B.mz, B.px, B.py, B.pz, B.pw, T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a trivector part exclusively (other grades = 0)
    constexpr MVec32_t(TriVec3dc<T> const& t) :
        MVec32_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), t.vx, t.vy, t.vz,
                 t.mx, t.my, t.mz, t.px, t.py, t.pz, t.pw, T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0))
    {
    }

    // assign a quadvector part exclusively (other grades = 0)
    constexpr MVec32_t(QuadVec3dc<T> const& Q) :
        MVec32_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), Q.x, Q.y,
                 Q.z, Q.w, Q.u, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    constexpr MVec32_t(PScalar3dc<T> ps) :
        MVec32_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign the even grade parts (odd grades = 0)
    constexpr MVec32_t(Scalar3dc<T> s, BiVec3dc<T> const& B, QuadVec3dc<T> const& Q) :
        MVec32_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), B.vx, B.vy, B.vz, B.mx,
                 B.my, B.mz, B.px, B.py, B.pz, B.pw, T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), Q.x, Q.y, Q.z, Q.w, Q.u,
                 T(0.0))
    {
    }

    // assign the odd grade parts (even grades = 0)
    constexpr MVec32_t(Vec3dc<T> const& v, TriVec3dc<T> const& t, PScalar3dc<T> ps) :
        MVec32_t(T(0.0), v.x, v.y, v.z, v.w, v.u, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), t.vx, t.vy, t.vz, t.mx, t.my,
                 t.mz, t.px, t.py, t.pz, t.pw, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(ps))
    {
    }

    // assign all six grade parts explicitly
    constexpr MVec32_t(Scalar3dc<T> s, Vec3dc<T> const& v, BiVec3dc<T> const& B,
                       TriVec3dc<T> const& t, QuadVec3dc<T> const& Q, PScalar3dc<T> ps) :
        MVec32_t(T(s), v.x, v.y, v.z, v.w, v.u, B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px,
                 B.py, B.pz, B.pw, t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py, t.pz,
                 t.pw, Q.x, Q.y, Q.z, Q.w, Q.u, T(ps))
    {
    }

    // assign from an even-grade multivector
    constexpr MVec32_t(MVec3dc_E<T> const& M) :
        MVec32_t(M.c0, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), M.c1, M.c2, M.c3, M.c4,
                 M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), M.c11, M.c12, M.c13,
                 M.c14, M.c15, T(0.0))
    {
    }

    // assign from an uneven (odd) grade multivector
    constexpr MVec32_t(MVec3dc_U<T> const& M) :
        MVec32_t(T(0.0), M.c0, M.c1, M.c2, M.c3, M.c4, T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), M.c5, M.c6, M.c7, M.c8,
                 M.c9, M.c10, M.c11, M.c12, M.c13, M.c14, T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), M.c15)
    {
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec32_t& operator+=(MVec32_t<U, mvec3dc_tag> const& v) noexcept
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
        this->c16 += v.c16;
        this->c17 += v.c17;
        this->c18 += v.c18;
        this->c19 += v.c19;
        this->c20 += v.c20;
        this->c21 += v.c21;
        this->c22 += v.c22;
        this->c23 += v.c23;
        this->c24 += v.c24;
        this->c25 += v.c25;
        this->c26 += v.c26;
        this->c27 += v.c27;
        this->c28 += v.c28;
        this->c29 += v.c29;
        this->c30 += v.c30;
        this->c31 += v.c31;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec32_t& operator-=(MVec32_t<U, mvec3dc_tag> const& v) noexcept
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
        this->c16 -= v.c16;
        this->c17 -= v.c17;
        this->c18 -= v.c18;
        this->c19 -= v.c19;
        this->c20 -= v.c20;
        this->c21 -= v.c21;
        this->c22 -= v.c22;
        this->c23 -= v.c23;
        this->c24 -= v.c24;
        this->c25 -= v.c25;
        this->c26 -= v.c26;
        this->c27 -= v.c27;
        this->c28 -= v.c28;
        this->c29 -= v.c29;
        this->c30 -= v.c30;
        this->c31 -= v.c31;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec32_t& operator*=(U s) noexcept
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
        this->c16 *= s;
        this->c17 *= s;
        this->c18 *= s;
        this->c19 *= s;
        this->c20 *= s;
        this->c21 *= s;
        this->c22 *= s;
        this->c23 *= s;
        this->c24 *= s;
        this->c25 *= s;
        this->c26 *= s;
        this->c27 *= s;
        this->c28 *= s;
        this->c29 *= s;
        this->c30 *= s;
        this->c31 *= s;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec32_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 32 comp.");
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
        this->c16 /= s;
        this->c17 /= s;
        this->c18 /= s;
        this->c19 /= s;
        this->c20 /= s;
        this->c21 /= s;
        this->c22 /= s;
        this->c23 /= s;
        this->c24 /= s;
        this->c25 /= s;
        this->c26 /= s;
        this->c27 /= s;
        this->c28 /= s;
        this->c29 /= s;
        this->c30 /= s;
        this->c31 /= s;
        return *this;
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec32_t<T, mvec3dc_tag>
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> gr0(MVec3dc<T> const& M)
{
    return Scalar3dc<T>(M.c0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> gr1(MVec3dc<T> const& M)
{
    return Vec3dc<T>(M.c1, M.c2, M.c3, M.c4, M.c5);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> gr2(MVec3dc<T> const& M)
{
    return BiVec3dc<T>(M.c6, M.c7, M.c8, M.c9, M.c10, M.c11, M.c12, M.c13, M.c14, M.c15);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> gr3(MVec3dc<T> const& M)
{
    return TriVec3dc<T>(M.c16, M.c17, M.c18, M.c19, M.c20, M.c21, M.c22, M.c23, M.c24,
                        M.c25);
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> gr4(MVec3dc<T> const& M)
{
    return QuadVec3dc<T>(M.c26, M.c27, M.c28, M.c29, M.c30);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> gr5(MVec3dc<T> const& M)
{
    return PScalar3dc<T>(M.c31);
}

////////////////////////////////////////////////////////////////////////////////
// addition/subtraction operations for mixed-parity grade combinations
// (same-parity combinations live in ga_mvec3dc_e.hpp / ga_mvec3dc_u.hpp)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(Scalar3dc<T> s, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s), v.x, v.y, v.z, v.w, v.u, ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(Scalar3dc<T> s,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                          t.pz, t.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(Scalar3dc<T> s, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(Vec3dc<T> const& v, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s), v.x, v.y, v.z, v.w, v.u, ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(Vec3dc<T> const& v,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), v.x, v.y, v.z, v.w, v.u, B.vx, B.vy, B.vz, B.mx,
                          B.my, B.mz, B.px, B.py, B.pz, B.pw, ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(Vec3dc<T> const& v,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), v.x, v.y, v.z, v.w, v.u, ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), Q.x, Q.y, Q.z, Q.w, Q.u,
                          ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(BiVec3dc<T> const& B,
                                                      Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), v.x, v.y, v.z, v.w, v.u, B.vx, B.vy, B.vz, B.mx,
                          B.my, B.mz, B.px, B.py, B.pz, B.pw, ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(BiVec3dc<T> const& B,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                          B.pz, B.pw, t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                          t.pz, t.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(BiVec3dc<T> const& B,
                                                      PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                          B.pz, B.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(TriVec3dc<T> const& t,
                                                      Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                          t.pz, t.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(TriVec3dc<T> const& t,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                          B.pz, B.pw, t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                          t.pz, t.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(TriVec3dc<T> const& t,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                          t.pz, t.pw, Q.x, Q.y, Q.z, Q.w, Q.u, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(QuadVec3dc<T> const& Q,
                                                      Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), v.x, v.y, v.z, v.w, v.u, ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), Q.x, Q.y, Q.z, Q.w, Q.u,
                          ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(QuadVec3dc<T> const& Q,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                          t.pz, t.pw, Q.x, Q.y, Q.z, Q.w, Q.u, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(QuadVec3dc<T> const& Q,
                                                      PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), Q.x, Q.y, Q.z, Q.w, Q.u, ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(PScalar3dc<T> ps, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(PScalar3dc<T> ps,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                          B.pz, B.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(PScalar3dc<T> ps,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), Q.x, Q.y, Q.z, Q.w, Q.u, ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(Scalar3dc<T> s, Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s), -v.x, -v.y, -v.z, -v.w, -v.u, ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(Scalar3dc<T> s,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), -t.vx, -t.vy, -t.vz, -t.mx, -t.my, -t.mz, -t.px,
                          -t.py, -t.pz, -t.pw, ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(Scalar3dc<T> s, PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), -ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(Vec3dc<T> const& v, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(-ctype(s), v.x, v.y, v.z, v.w, v.u, ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(Vec3dc<T> const& v,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), v.x, v.y, v.z, v.w, v.u, -B.vx, -B.vy, -B.vz, -B.mx,
                          -B.my, -B.mz, -B.px, -B.py, -B.pz, -B.pw, ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(Vec3dc<T> const& v,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), v.x, v.y, v.z, v.w, v.u, ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), -Q.x, -Q.y, -Q.z, -Q.w,
                          -Q.u, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(BiVec3dc<T> const& B,
                                                      Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), -v.x, -v.y, -v.z, -v.w, -v.u, B.vx, B.vy, B.vz,
                          B.mx, B.my, B.mz, B.px, B.py, B.pz, B.pw, ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(BiVec3dc<T> const& B,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                          B.pz, B.pw, -t.vx, -t.vy, -t.vz, -t.mx, -t.my, -t.mz, -t.px,
                          -t.py, -t.pz, -t.pw, ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(BiVec3dc<T> const& B,
                                                      PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                          B.pz, B.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), -ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(TriVec3dc<T> const& t,
                                                      Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(-ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                          t.pz, t.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(TriVec3dc<T> const& t,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), -B.vx, -B.vy, -B.vz, -B.mx, -B.my, -B.mz, -B.px,
                          -B.py, -B.pz, -B.pw, t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px,
                          t.py, t.pz, t.pw, ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(TriVec3dc<T> const& t,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                          t.pz, t.pw, -Q.x, -Q.y, -Q.z, -Q.w, -Q.u, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(QuadVec3dc<T> const& Q,
                                                      Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), -v.x, -v.y, -v.z, -v.w, -v.u, ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), Q.x, Q.y, Q.z,
                          Q.w, Q.u, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(QuadVec3dc<T> const& Q,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), -t.vx, -t.vy, -t.vz, -t.mx, -t.my, -t.mz, -t.px,
                          -t.py, -t.pz, -t.pw, Q.x, Q.y, Q.z, Q.w, Q.u, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(QuadVec3dc<T> const& Q,
                                                      PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), Q.x, Q.y, Q.z, Q.w, Q.u, -ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(PScalar3dc<T> ps, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(-ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(PScalar3dc<T> ps,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), -B.vx, -B.vy, -B.vz, -B.mx, -B.my, -B.mz, -B.px,
                          -B.py, -B.pz, -B.pw, ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(PScalar3dc<T> ps,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                          ctype(0.0), -Q.x, -Q.y, -Q.z, -Q.w, -Q.u, ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(Scalar3dc<T> s, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s) + M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8,
                          M.c9, M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17,
                          M.c18, M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25, M.c26,
                          M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(MVec3dc<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0 + ctype(s), M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8,
                          M.c9, M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17,
                          M.c18, M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25, M.c26,
                          M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(Vec3dc<T> const& v,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, v.x + M.c1, v.y + M.c2, v.z + M.c3, v.w + M.c4,
                          v.u + M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11, M.c12, M.c13,
                          M.c14, M.c15, M.c16, M.c17, M.c18, M.c19, M.c20, M.c21, M.c22,
                          M.c23, M.c24, M.c25, M.c26, M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(MVec3dc<T> const& M,
                                                      Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1 + v.x, M.c2 + v.y, M.c3 + v.z, M.c4 + v.w,
                          M.c5 + v.u, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11, M.c12, M.c13,
                          M.c14, M.c15, M.c16, M.c17, M.c18, M.c19, M.c20, M.c21, M.c22,
                          M.c23, M.c24, M.c25, M.c26, M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(BiVec3dc<T> const& B,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, B.vx + M.c6, B.vy + M.c7,
                          B.vz + M.c8, B.mx + M.c9, B.my + M.c10, B.mz + M.c11,
                          B.px + M.c12, B.py + M.c13, B.pz + M.c14, B.pw + M.c15, M.c16,
                          M.c17, M.c18, M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25,
                          M.c26, M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(MVec3dc<T> const& M,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6 + B.vx, M.c7 + B.vy,
                          M.c8 + B.vz, M.c9 + B.mx, M.c10 + B.my, M.c11 + B.mz,
                          M.c12 + B.px, M.c13 + B.py, M.c14 + B.pz, M.c15 + B.pw, M.c16,
                          M.c17, M.c18, M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25,
                          M.c26, M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(TriVec3dc<T> const& t,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                          M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, t.vx + M.c16,
                          t.vy + M.c17, t.vz + M.c18, t.mx + M.c19, t.my + M.c20,
                          t.mz + M.c21, t.px + M.c22, t.py + M.c23, t.pz + M.c24,
                          t.pw + M.c25, M.c26, M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(MVec3dc<T> const& M,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                          M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16 + t.vx,
                          M.c17 + t.vy, M.c18 + t.vz, M.c19 + t.mx, M.c20 + t.my,
                          M.c21 + t.mz, M.c22 + t.px, M.c23 + t.py, M.c24 + t.pz,
                          M.c25 + t.pw, M.c26, M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(QuadVec3dc<T> const& Q,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                          M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17, M.c18,
                          M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25, Q.x + M.c26,
                          Q.y + M.c27, Q.z + M.c28, Q.w + M.c29, Q.u + M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(MVec3dc<T> const& M,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                          M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17, M.c18,
                          M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25, M.c26 + Q.x,
                          M.c27 + Q.y, M.c28 + Q.z, M.c29 + Q.w, M.c30 + Q.u, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(PScalar3dc<T> ps,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                          M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17, M.c18,
                          M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25, M.c26, M.c27,
                          M.c28, M.c29, M.c30, ctype(ps) + M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator+(MVec3dc<T> const& M,
                                                      PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                          M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17, M.c18,
                          M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25, M.c26, M.c27,
                          M.c28, M.c29, M.c30, M.c31 + ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(Scalar3dc<T> s, MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(ctype(s) - M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6,
                          -M.c7, -M.c8, -M.c9, -M.c10, -M.c11, -M.c12, -M.c13, -M.c14,
                          -M.c15, -M.c16, -M.c17, -M.c18, -M.c19, -M.c20, -M.c21, -M.c22,
                          -M.c23, -M.c24, -M.c25, -M.c26, -M.c27, -M.c28, -M.c29, -M.c30,
                          -M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(MVec3dc<T> const& M, Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0 - ctype(s), M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8,
                          M.c9, M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17,
                          M.c18, M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25, M.c26,
                          M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(Vec3dc<T> const& v,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(-M.c0, v.x - M.c1, v.y - M.c2, v.z - M.c3, v.w - M.c4,
                          v.u - M.c5, -M.c6, -M.c7, -M.c8, -M.c9, -M.c10, -M.c11, -M.c12,
                          -M.c13, -M.c14, -M.c15, -M.c16, -M.c17, -M.c18, -M.c19, -M.c20,
                          -M.c21, -M.c22, -M.c23, -M.c24, -M.c25, -M.c26, -M.c27, -M.c28,
                          -M.c29, -M.c30, -M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(MVec3dc<T> const& M,
                                                      Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1 - v.x, M.c2 - v.y, M.c3 - v.z, M.c4 - v.w,
                          M.c5 - v.u, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11, M.c12, M.c13,
                          M.c14, M.c15, M.c16, M.c17, M.c18, M.c19, M.c20, M.c21, M.c22,
                          M.c23, M.c24, M.c25, M.c26, M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(BiVec3dc<T> const& B,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(
        -M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, B.vx - M.c6, B.vy - M.c7, B.vz - M.c8,
        B.mx - M.c9, B.my - M.c10, B.mz - M.c11, B.px - M.c12, B.py - M.c13, B.pz - M.c14,
        B.pw - M.c15, -M.c16, -M.c17, -M.c18, -M.c19, -M.c20, -M.c21, -M.c22, -M.c23,
        -M.c24, -M.c25, -M.c26, -M.c27, -M.c28, -M.c29, -M.c30, -M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(MVec3dc<T> const& M,
                                                      BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6 - B.vx, M.c7 - B.vy,
                          M.c8 - B.vz, M.c9 - B.mx, M.c10 - B.my, M.c11 - B.mz,
                          M.c12 - B.px, M.c13 - B.py, M.c14 - B.pz, M.c15 - B.pw, M.c16,
                          M.c17, M.c18, M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25,
                          M.c26, M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(TriVec3dc<T> const& t,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(
        -M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, -M.c7, -M.c8, -M.c9, -M.c10,
        -M.c11, -M.c12, -M.c13, -M.c14, -M.c15, t.vx - M.c16, t.vy - M.c17, t.vz - M.c18,
        t.mx - M.c19, t.my - M.c20, t.mz - M.c21, t.px - M.c22, t.py - M.c23,
        t.pz - M.c24, t.pw - M.c25, -M.c26, -M.c27, -M.c28, -M.c29, -M.c30, -M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(MVec3dc<T> const& M,
                                                      TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                          M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16 - t.vx,
                          M.c17 - t.vy, M.c18 - t.vz, M.c19 - t.mx, M.c20 - t.my,
                          M.c21 - t.mz, M.c22 - t.px, M.c23 - t.py, M.c24 - t.pz,
                          M.c25 - t.pw, M.c26, M.c27, M.c28, M.c29, M.c30, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(QuadVec3dc<T> const& Q,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(-M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, -M.c7, -M.c8,
                          -M.c9, -M.c10, -M.c11, -M.c12, -M.c13, -M.c14, -M.c15, -M.c16,
                          -M.c17, -M.c18, -M.c19, -M.c20, -M.c21, -M.c22, -M.c23, -M.c24,
                          -M.c25, Q.x - M.c26, Q.y - M.c27, Q.z - M.c28, Q.w - M.c29,
                          Q.u - M.c30, -M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(MVec3dc<T> const& M,
                                                      QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                          M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17, M.c18,
                          M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25, M.c26 - Q.x,
                          M.c27 - Q.y, M.c28 - Q.z, M.c29 - Q.w, M.c30 - Q.u, M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(PScalar3dc<T> ps,
                                                      MVec3dc<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(-M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, -M.c7, -M.c8,
                          -M.c9, -M.c10, -M.c11, -M.c12, -M.c13, -M.c14, -M.c15, -M.c16,
                          -M.c17, -M.c18, -M.c19, -M.c20, -M.c21, -M.c22, -M.c23, -M.c24,
                          -M.c25, -M.c26, -M.c27, -M.c28, -M.c29, -M.c30,
                          ctype(ps) - M.c31);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc<std::common_type_t<T, U>> operator-(MVec3dc<T> const& M,
                                                      PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                          M.c10, M.c11, M.c12, M.c13, M.c14, M.c15, M.c16, M.c17, M.c18,
                          M.c19, M.c20, M.c21, M.c22, M.c23, M.c24, M.c25, M.c26, M.c27,
                          M.c28, M.c29, M.c30, M.c31 - ctype(ps));
}

} // namespace hd::ga
