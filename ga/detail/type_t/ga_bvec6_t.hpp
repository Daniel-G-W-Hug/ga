#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <algorithm> // std::max
#include <cmath>     // std::abs, std::sqrt
#include <concepts>  // numeric_type<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "../ga_error_handling.hpp"
#include "ga_type_tags.hpp"
#include "ga_vec3_t.hpp" // Required for Vec3_t<T> in constructor

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T, Tag> definition of a bivector with 6 components
// vx, vy, vz, mx, my, mz
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T, typename Tag = default_tag>
    requires numeric_type<T>
struct BVec6_t {

    // ctors

    // (all grades = 0)
    BVec6_t() = default;

    // assign all components directly
    constexpr BVec6_t(T vx_in, T vy_in, T vz_in, T mx_in, T my_in, T mz_in) :
        vx(vx_in), vy(vy_in), vz(vz_in), mx(mx_in), my(my_in), mz(mz_in)
    {
    }

    // assign the vector parts separately to the bivector
    // (both vectors must have the same tag to ensure type safety)
    template <typename VecTag>
    constexpr BVec6_t(Vec3_t<T, VecTag> const& lv, Vec3_t<T, VecTag> const& lm) :
        vx(lv.x), vy(lv.y), vz(lv.z), mx(lm.x), my(lm.y), mz(lm.z)
    {
    }

    // assign the vector parts from two Vec3_t with different tags
    // (for special cases where semantically different types are combined,
    //  e.g., Line3d: direction vector + moment bivector)
    template <typename Vec1Tag, typename Vec2Tag>
    constexpr BVec6_t(Vec3_t<T, Vec1Tag> const& lv, Vec3_t<T, Vec2Tag> const& lm) :
        vx(lv.x), vy(lv.y), vz(lv.z), mx(lm.x), my(lm.y), mz(lm.z)
    {
    }

    // cp & mv ctor
    BVec6_t(BVec6_t const&) = default;
    BVec6_t(BVec6_t&&) noexcept = default;

    // cp & mv assignment
    BVec6_t& operator=(BVec6_t const&) = default;
    BVec6_t& operator=(BVec6_t&&) noexcept = default;

    // dtor
    ~BVec6_t() = default;

    // floating point type conversion
    template <typename U>
        requires(numeric_type<U>)
    constexpr BVec6_t(BVec6_t<U, Tag> const& B) :
        vx(B.vx), vy(B.vy), vz(B.vz), mx(B.mx), my(B.my), mz(B.mz)
    {
    }

    friend void swap(BVec6_t& lhs, BVec6_t& rhs) noexcept
    {
        using std::swap;
        swap(lhs.vx, rhs.vx);
        swap(lhs.vy, rhs.vy);
        swap(lhs.vz, rhs.vz);
        swap(lhs.mx, rhs.mx);
        swap(lhs.my, rhs.my);
        swap(lhs.mz, rhs.mz);
    }

    T vx{}; // as BiVec3dp<T> maps to basis bivector e4^e1 - as Line3dp<T> to vx
    T vy{}; // as BiVec3dp<T> maps to basis bivector e4^e2 - as Line3dp<T> to vy
    T vz{}; // as BiVec3dp<T> maps to basis bivector e4^e3 - as Line3dp<T> to vz
    T mx{}; // as BiVec3dp<T> maps to basis bivector e2^e3 - as Line3dp<T> to mx
    T my{}; // as BiVec3dp<T> maps to basis bivector e3^e1 - as Line3dp<T> to my
    T mz{}; // as BiVec3dp<T> maps to basis bivector e1^e2 - as Line3dp<T> to mz

    // with l = l(lv,lm) return lv(l) and lm(l) each as vec3d
    //      - lm is a bivec3d, but a vec3d is sometimes more practical,
    //        e.g. for simple use in the dot product between lv and lm to check
    //        for proper lines satisfying the geometric constraint
    //      - if lm is needed as bivec3d use lm_as_bivec = hd::ga::ega::cmpl(lm)
    constexpr Vec3_t<T, vec3d_tag> lv() const { return Vec3_t<T, vec3d_tag>(vx, vy, vz); }
    constexpr Vec3_t<T, vec3d_tag> lm() const { return Vec3_t<T, vec3d_tag>(mx, my, mz); }

    template <typename U>
        requires(numeric_type<U>)
    BVec6_t& operator+=(BVec6_t<U, Tag> const& B) noexcept
    {
        vx += B.vx;
        vy += B.vy;
        vz += B.vz;
        mx += B.mx;
        my += B.my;
        mz += B.mz;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    BVec6_t& operator-=(BVec6_t<U, Tag> const& B) noexcept
    {
        vx -= B.vx;
        vy -= B.vy;
        vz -= B.vz;
        mx -= B.mx;
        my -= B.my;
        mz -= B.mz;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    BVec6_t& operator*=(U s) noexcept
    {
        vx *= s;
        vy *= s;
        vz *= s;
        mx *= s;
        my *= s;
        mz *= s;
        return (*this);
    }

    template <typename U>
        requires(numeric_type<U>)
    BVec6_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "bivector_division 6 comp.");
        vx /= s;
        vy /= s;
        vz /= s;
        mx /= s;
        my /= s;
        mz /= s;
        return (*this);
    }
};

////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T> core operations
////////////////////////////////////////////////////////////////////////////////

// equality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator==(BVec6_t<T, Tag> const& lhs, BVec6_t<U, Tag> const& rhs)
{
    // componentwise comparison
    // equality implies same magnitude and direction
    // comparison is not exact, but accepts epsilon deviations
    auto abs_delta_vx = std::abs(lhs.vx - rhs.vx);
    auto abs_delta_vy = std::abs(lhs.vy - rhs.vy);
    auto abs_delta_vz = std::abs(lhs.vz - rhs.vz);
    auto abs_delta_mx = std::abs(lhs.mx - rhs.mx);
    auto abs_delta_my = std::abs(lhs.my - rhs.my);
    auto abs_delta_mz = std::abs(lhs.mz - rhs.mz);
    auto constexpr delta_eps = detail::safe_epsilon<T, U>();
    return (abs_delta_vx < delta_eps && abs_delta_vy < delta_eps &&
            abs_delta_vz < delta_eps && abs_delta_mx < delta_eps &&
            abs_delta_my < delta_eps && abs_delta_mz < delta_eps);
}

// inequality - only allows comparison between same tag types
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
bool operator!=(BVec6_t<T, Tag> const& lhs, BVec6_t<U, Tag> const& rhs)
{
    return !(lhs == rhs);
}

// unary minus
template <typename T, typename Tag>
    requires(numeric_type<T>)
constexpr BVec6_t<T, Tag> operator-(BVec6_t<T, Tag> const& B)
{
    return BVec6_t<T, Tag>(-B.vx, -B.vy, -B.vz, -B.mx, -B.my, -B.mz);
}

// add multivectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BVec6_t<std::common_type_t<T, U>, Tag> operator+(BVec6_t<T, Tag> const& v1,
                                                           BVec6_t<U, Tag> const& v2)
{
    return BVec6_t<std::common_type_t<T, U>, Tag>(v1.vx + v2.vx, v1.vy + v2.vy,
                                                  v1.vz + v2.vz, v1.mx + v2.mx,
                                                  v1.my + v2.my, v1.mz + v2.mz);
}

// substract multivectors
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BVec6_t<std::common_type_t<T, U>, Tag> operator-(BVec6_t<T, Tag> const& v1,
                                                           BVec6_t<U, Tag> const& v2)
{
    return BVec6_t<std::common_type_t<T, U>, Tag>(v1.vx - v2.vx, v1.vy - v2.vy,
                                                  v1.vz - v2.vz, v1.mx - v2.mx,
                                                  v1.my - v2.my, v1.mz - v2.mz);
}

// multiply a multivector with a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BVec6_t<std::common_type_t<T, U>, Tag> operator*(BVec6_t<T, Tag> const& B, U s)
{
    return BVec6_t<std::common_type_t<T, U>, Tag>(B.vx * s, B.vy * s, B.vz * s, B.mx * s,
                                                  B.my * s, B.mz * s);
}

template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
constexpr BVec6_t<std::common_type_t<T, U>, Tag> operator*(T s, BVec6_t<U, Tag> const& B)
{
    return BVec6_t<std::common_type_t<T, U>, Tag>(B.vx * s, B.vy * s, B.vz * s, B.mx * s,
                                                  B.my * s, B.mz * s);
}

// devide a multivector by a scalar
template <typename T, typename U, typename Tag>
    requires(numeric_type<T> && numeric_type<U>)
inline BVec6_t<std::common_type_t<T, U>, Tag> operator/(BVec6_t<T, Tag> const& B, U s)
{
    detail::check_division_by_zero<T, U>(s, "bivector division");
    using ctype = std::common_type_t<T, U>;
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return BVec6_t<ctype, Tag>(B.vx * inv, B.vy * inv, B.vz * inv, B.mx * inv, B.my * inv,
                               B.mz * inv);
}

////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(numeric_type<T>)
std::ostream& operator<<(std::ostream& os, BVec6_t<T, Tag> const& B)
{
    os << "(" << B.vx << "," << B.vy << "," << B.vz << "," << B.mx << "," << B.my << ","
       << B.mz << ")";
    return os;
}

} // namespace hd::ga