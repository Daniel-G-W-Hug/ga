#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type3dc.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec16_t including its ctors and add specific ctors for MVec16_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dc_u_tag
// (cga3dc ODD layout: c0..c4 = vector, c5..c14 = trivector, c15 = pseudoscalar)
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
struct MVec16_t<T, mvec3dc_u_tag> : public MVec16_t<T, default_tag> {

    using MVec16_t<T, default_tag>::MVec16_t; // inherit base class ctors

    // assign a vector part exclusively (other grades = 0)
    constexpr MVec16_t(Vec3dc<T> const& v) :
        MVec16_t(v.x, v.y, v.z, v.w, v.u, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a trivector part exclusively (other grades = 0)
    constexpr MVec16_t(TriVec3dc<T> const& t) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), t.vx, t.vy, t.vz, t.mx, t.my,
                 t.mz, t.px, t.py, t.pz, t.pw, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    constexpr MVec16_t(PScalar3dc<T> ps) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign vector and trivector part (pseudoscalar part = 0)
    constexpr MVec16_t(Vec3dc<T> const& v, TriVec3dc<T> const& t) :
        MVec16_t(v.x, v.y, v.z, v.w, v.u, t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                 t.pz, t.pw, T(0.0))
    {
    }

    // assign trivector and pseudoscalar part (vector part = 0)
    constexpr MVec16_t(TriVec3dc<T> const& t, PScalar3dc<T> ps) :
        MVec16_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), t.vx, t.vy, t.vz, t.mx, t.my,
                 t.mz, t.px, t.py, t.pz, t.pw, T(ps))
    {
    }

    // assign vector and pseudoscalar part (trivector part = 0)
    constexpr MVec16_t(Vec3dc<T> const& v, PScalar3dc<T> ps) :
        MVec16_t(v.x, v.y, v.z, v.w, v.u, T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0),
                 T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign all three grade parts explicitly
    constexpr MVec16_t(Vec3dc<T> const& v, TriVec3dc<T> const& t, PScalar3dc<T> ps) :
        MVec16_t(v.x, v.y, v.z, v.w, v.u, t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py,
                 t.pz, t.pw, T(ps))
    {
    }

    // Override compound assignment operators to return correct derived type
    template <typename U>
        requires(numeric_type<U>)
    MVec16_t& operator+=(MVec16_t<U, mvec3dc_u_tag> const& v) noexcept
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
    MVec16_t& operator-=(MVec16_t<U, mvec3dc_u_tag> const& v) noexcept
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
// define grade operations for partial specialization MVec16_t<T, mvec3dc_u_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of an uneven (odd) 3dc multivector
//
// grade 1: gr1() - vector
// grade 3: gr3() - trivector
// grade 5: gr5() - pseudoscalar

template <typename T>
    requires(numeric_type<T>)
constexpr Vec3dc<T> gr1(MVec3dc_U<T> const& M)
{
    return Vec3dc<T>(M.c0, M.c1, M.c2, M.c3, M.c4);
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec3dc<T> gr3(MVec3dc_U<T> const& M)
{
    return TriVec3dc<T>(M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11, M.c12, M.c13, M.c14);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar3dc<T> gr5(MVec3dc_U<T> const& M)
{
    return PScalar3dc<T>(M.c15);
}

////////////////////////////////////////////////////////////////////////////////
// addition/subtraction operations to combine odd grade parts
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(Vec3dc<T> const& v,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(v.x, v.y, v.z, v.w, v.u, t.vx, t.vy, t.vz, t.mx, t.my, t.mz,
                            t.px, t.py, t.pz, t.pw, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(TriVec3dc<T> const& t,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(v.x, v.y, v.z, v.w, v.u, t.vx, t.vy, t.vz, t.mx, t.my, t.mz,
                            t.px, t.py, t.pz, t.pw, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(Vec3dc<T> const& v,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(v.x, v.y, v.z, v.w, v.u, ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(PScalar3dc<T> ps,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(v.x, v.y, v.z, v.w, v.u, ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(TriVec3dc<T> const& t,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py, t.pz, t.pw,
                            ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(PScalar3dc<T> ps,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py, t.pz, t.pw,
                            ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(Vec3dc<T> const& v,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(v.x + M.c0, v.y + M.c1, v.z + M.c2, v.w + M.c3, v.u + M.c4,
                            M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11, M.c12, M.c13,
                            M.c14, M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(MVec3dc_U<T> const& M,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(M.c0 + v.x, M.c1 + v.y, M.c2 + v.z, M.c3 + v.w, M.c4 + v.u,
                            M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11, M.c12, M.c13,
                            M.c14, M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(TriVec3dc<T> const& t,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, t.vx + M.c5, t.vy + M.c6,
                            t.vz + M.c7, t.mx + M.c8, t.my + M.c9, t.mz + M.c10,
                            t.px + M.c11, t.py + M.c12, t.pz + M.c13, t.pw + M.c14,
                            M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(MVec3dc_U<T> const& M,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5 + t.vx, M.c6 + t.vy,
                            M.c7 + t.vz, M.c8 + t.mx, M.c9 + t.my, M.c10 + t.mz,
                            M.c11 + t.px, M.c12 + t.py, M.c13 + t.pz, M.c14 + t.pw,
                            M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(PScalar3dc<T> ps,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                            M.c10, M.c11, M.c12, M.c13, M.c14, ctype(ps) + M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator+(MVec3dc_U<T> const& M,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                            M.c10, M.c11, M.c12, M.c13, M.c14, M.c15 + ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(Vec3dc<T> const& v,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(v.x, v.y, v.z, v.w, v.u, -t.vx, -t.vy, -t.vz, -t.mx, -t.my,
                            -t.mz, -t.px, -t.py, -t.pz, -t.pw, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(TriVec3dc<T> const& t,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(-v.x, -v.y, -v.z, -v.w, -v.u, t.vx, t.vy, t.vz, t.mx, t.my,
                            t.mz, t.px, t.py, t.pz, t.pw, ctype(0.0));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(Vec3dc<T> const& v,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(v.x, v.y, v.z, v.w, v.u, ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), -ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(PScalar3dc<T> ps,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(-v.x, -v.y, -v.z, -v.w, -v.u, ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            ctype(0.0), ctype(0.0), ctype(0.0), ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(TriVec3dc<T> const& t,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            t.vx, t.vy, t.vz, t.mx, t.my, t.mz, t.px, t.py, t.pz, t.pw,
                            -ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(PScalar3dc<T> ps,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0), ctype(0.0),
                            -t.vx, -t.vy, -t.vz, -t.mx, -t.my, -t.mz, -t.px, -t.py, -t.pz,
                            -t.pw, ctype(ps));
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(Vec3dc<T> const& v,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(v.x - M.c0, v.y - M.c1, v.z - M.c2, v.w - M.c3, v.u - M.c4,
                            -M.c5, -M.c6, -M.c7, -M.c8, -M.c9, -M.c10, -M.c11, -M.c12,
                            -M.c13, -M.c14, -M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(MVec3dc_U<T> const& M,
                                                        Vec3dc<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(M.c0 - v.x, M.c1 - v.y, M.c2 - v.z, M.c3 - v.w, M.c4 - v.u,
                            M.c5, M.c6, M.c7, M.c8, M.c9, M.c10, M.c11, M.c12, M.c13,
                            M.c14, M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(TriVec3dc<T> const& t,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(-M.c0, -M.c1, -M.c2, -M.c3, -M.c4, t.vx - M.c5, t.vy - M.c6,
                            t.vz - M.c7, t.mx - M.c8, t.my - M.c9, t.mz - M.c10,
                            t.px - M.c11, t.py - M.c12, t.pz - M.c13, t.pw - M.c14,
                            -M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(MVec3dc_U<T> const& M,
                                                        TriVec3dc<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5 - t.vx, M.c6 - t.vy,
                            M.c7 - t.vz, M.c8 - t.mx, M.c9 - t.my, M.c10 - t.mz,
                            M.c11 - t.px, M.c12 - t.py, M.c13 - t.pz, M.c14 - t.pw,
                            M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(PScalar3dc<T> ps,
                                                        MVec3dc_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(-M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, -M.c7, -M.c8,
                            -M.c9, -M.c10, -M.c11, -M.c12, -M.c13, -M.c14,
                            ctype(ps) - M.c15);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec3dc_U<std::common_type_t<T, U>> operator-(MVec3dc_U<T> const& M,
                                                        PScalar3dc<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dc_U<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7, M.c8, M.c9,
                            M.c10, M.c11, M.c12, M.c13, M.c14, M.c15 - ctype(ps));
}

} // namespace hd::ga
