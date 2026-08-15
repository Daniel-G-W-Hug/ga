#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type3dc.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec16_t including its ctors and add specific ctors for MVec16_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dc_e_tag
// (cga3dc EVEN layout: c0 = scalar, c1..c10 = bivector, c11..c15 = quadvector)
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
struct MVec16_t<T, mvec3dc_e_tag> : public MVec16_t<T, default_tag> {

    using MVec16_t<T, default_tag>::MVec16_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    constexpr MVec16_t(Scalar3dc<T> s) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    constexpr MVec16_t(BiVec3dc<T> const& B) :
        MVec16_t(T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py, B.pz, B.pw,
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a quadvector part exclusively (other grades = 0)
    constexpr MVec16_t(QuadVec3dc<T> const& Q) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), Q.x, Q.y, Q.z, Q.w, Q.u)
    {
    }

    // assign scalar and bivector part (quadvector part = 0)
    constexpr MVec16_t(Scalar3dc<T> s, BiVec3dc<T> const& B) :
        MVec16_t(T(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py, B.pz, B.pw, T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign bivector and quadvector part (scalar part = 0)
    constexpr MVec16_t(BiVec3dc<T> const& B, QuadVec3dc<T> const& Q) :
        MVec16_t(T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py, B.pz, B.pw, Q.x,
                 Q.y, Q.z, Q.w, Q.u)
    {
    }

    // assign scalar and quadvector part (bivector part = 0)
    constexpr MVec16_t(Scalar3dc<T> s, QuadVec3dc<T> const& Q) :
        MVec16_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), Q.x, Q.y, Q.z, Q.w, Q.u)
    {
    }

    // assign all three grade parts explicitly
    constexpr MVec16_t(Scalar3dc<T> s, BiVec3dc<T> const& B, QuadVec3dc<T> const& Q) :
        MVec16_t(T(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py, B.pz, B.pw, Q.x,
                 Q.y, Q.z, Q.w, Q.u)
    {
    }

    // Override compound assignment operators to return correct derived type
    template <typename U>
        requires(numeric_type<U>)
    MVec16_t& operator+=(MVec16_t<U, mvec3dc_e_tag> const& v) noexcept
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
    MVec16_t& operator-=(MVec16_t<U, mvec3dc_e_tag> const& v) noexcept
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
// define grade operations for partial specialization MVec16_t<T, mvec3dc_e_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of an even 3dc multivector
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector
// grade 4: gr4() - quadvector

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar3dc<T> gr0(MVec3dc_E<T> const& M)
{
    return Scalar3dc<T>(M.c0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec3dc<T> gr2(MVec3dc_E<T> const& M)
{
    return BiVec3dc<T>(M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9, M.c10);
}

template <typename T>
    requires(numeric_type<T>)
constexpr QuadVec3dc<T> gr4(MVec3dc_E<T> const& M)
{
    return QuadVec3dc<T>(M.c11, M.c12, M.c13, M.c14, M.c15);
}

////////////////////////////////////////////////////////////////////////////////
// addition/subtraction operations to combine even grade parts
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(Scalar3dc<T> s,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                            B.pz, B.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(BiVec3dc<T> const& B,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                            B.pz, B.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(Scalar3dc<T> s,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), Q.x, Q.y, Q.z, Q.w, Q.u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(QuadVec3dc<T> const& Q,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), Q.x, Q.y, Q.z, Q.w, Q.u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(BiVec3dc<T> const& B,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                            B.pz, B.pw, Q.x, Q.y, Q.z, Q.w, Q.u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(QuadVec3dc<T> const& Q,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                            B.pz, B.pw, Q.x, Q.y, Q.z, Q.w, Q.u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(Scalar3dc<T> s,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(s) + M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7,
                            M.c8, M.c9, M.c10, M.c11, M.c12, M.c13, M.c14, M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(MVec3dc_E<T> const& M,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(M.c0 + ctype(s), M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7,
                            M.c8, M.c9, M.c10, M.c11, M.c12, M.c13, M.c14, M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(BiVec3dc<T> const& B,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(M.c0, B.vx + M.c1, B.vy + M.c2, B.vz + M.c3, B.mx + M.c4,
                            B.my + M.c5, B.mz + M.c6, B.px + M.c7, B.py + M.c8,
                            B.pz + M.c9, B.pw + M.c10, M.c11, M.c12, M.c13, M.c14, M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(MVec3dc_E<T> const& M,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(M.c0, M.c1 + B.vx, M.c2 + B.vy, M.c3 + B.vz, M.c4 + B.mx,
                            M.c5 + B.my, M.c6 + B.mz, M.c7 + B.px, M.c8 + B.py,
                            M.c9 + B.pz, M.c10 + B.pw, M.c11, M.c12, M.c13, M.c14, M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(QuadVec3dc<T> const& Q,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                            M.c10, Q.x + M.c11, Q.y + M.c12, Q.z + M.c13, Q.w + M.c14,
                            Q.u + M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator+(MVec3dc_E<T> const& M,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                            M.c10, M.c11 + Q.x, M.c12 + Q.y, M.c13 + Q.z, M.c14 + Q.w,
                            M.c15 + Q.u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(Scalar3dc<T> s,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(s), -B.vx, -B.vy, -B.vz, -B.mx, -B.my, -B.mz, -B.px,
                            -B.py, -B.pz, -B.pw, ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(BiVec3dc<T> const& B,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(-ctype(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                            B.pz, B.pw, ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(Scalar3dc<T> s,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), -Q.x, -Q.y, -Q.z, -Q.w, -Q.u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(QuadVec3dc<T> const& Q,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(-ctype(s), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), Q.x, Q.y, Q.z, Q.w, Q.u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(BiVec3dc<T> const& B,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, B.px, B.py,
                            B.pz, B.pw, -Q.x, -Q.y, -Q.z, -Q.w, -Q.u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(QuadVec3dc<T> const& Q,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(0.0), -B.vx, -B.vy, -B.vz, -B.mx, -B.my, -B.mz, -B.px,
                            -B.py, -B.pz, -B.pw, Q.x, Q.y, Q.z, Q.w, Q.u);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(Scalar3dc<T> s,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(ctype(s) - M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6,
                            -M.c7, -M.c8, -M.c9, -M.c10, -M.c11, -M.c12, -M.c13, -M.c14,
                            -M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(MVec3dc_E<T> const& M,
                                                        Scalar3dc<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(M.c0 - ctype(s), M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7,
                            M.c8, M.c9, M.c10, M.c11, M.c12, M.c13, M.c14, M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(BiVec3dc<T> const& B,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(-M.c0, B.vx - M.c1, B.vy - M.c2, B.vz - M.c3, B.mx - M.c4,
                            B.my - M.c5, B.mz - M.c6, B.px - M.c7, B.py - M.c8,
                            B.pz - M.c9, B.pw - M.c10, -M.c11, -M.c12, -M.c13, -M.c14,
                            -M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(MVec3dc_E<T> const& M,
                                                        BiVec3dc<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(M.c0, M.c1 - B.vx, M.c2 - B.vy, M.c3 - B.vz, M.c4 - B.mx,
                            M.c5 - B.my, M.c6 - B.mz, M.c7 - B.px, M.c8 - B.py,
                            M.c9 - B.pz, M.c10 - B.pw, M.c11, M.c12, M.c13, M.c14, M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(QuadVec3dc<T> const& Q,
                                                        MVec3dc_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(-M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, -M.c7, -M.c8,
                            -M.c9, -M.c10, Q.x - M.c11, Q.y - M.c12, Q.z - M.c13,
                            Q.w - M.c14, Q.u - M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_E<std::common_type_t<T, U>> operator-(MVec3dc_E<T> const& M,
                                                        QuadVec3dc<U> const& Q)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_E<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                            M.c10, M.c11 - Q.x, M.c12 - Q.y, M.c13 - Q.z, M.c14 - Q.w,
                            M.c15 - Q.u);
}

} // namespace hd::ga
