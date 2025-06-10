#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type_tags.hpp"

#include "ga_scalar_t.hpp" // scalar, pseudoscalar
#include "ga_vec3_t.hpp"   // vector, bivector

#include "ga_mvec4_t.hpp" // even and odd multivectors
#include "ga_mvec8_t.hpp" // full multivector

#include "ga_mvec2_t.hpp" // for DualNum2dp
#include "ga_type2d.hpp"  // for Vec2d<T>

/////////////////////////////////////////////////////////////////////////////////////////
//
// consistent type definitions for pga2dp (in namespace hd::ga)
//
// convenience type aliases (in namespace hd::ga::pga; needed there for name resolution)
//     - Vector2d
//     - Point2d
//     - Point2dp
//     - Line2d
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar2dp = Scalar_t<T, scalar2dp_tag>;
template <typename T> using Vec2dp = Vec3_t<T, vec2dp_tag>;
template <typename T> using BiVec2dp = Vec3_t<T, bivec2dp_tag>;
template <typename T> using PScalar2dp = Scalar_t<T, pscalar2dp_tag>;
template <typename T> using MVec2dp_E = MVec4_t<T, mvec2dp_e_tag>;
template <typename T> using MVec2dp_U = MVec4_t<T, mvec2dp_u_tag>;
template <typename T> using MVec2dp = MVec8_t<T, mvec2dp_tag>;

// for the geometric norm as homogeneous quantity,
// which has a scalar part (c0 component) and a pseudoscalar part (c1 component)
template <typename T> using DualNum2dp = MVec2_t<T, dual_number2dp_tag>;


// Vec2dp: represents as direction when z=0.0 (a point at infinity in that direction)
//         the direction is x*e1 + y*e2
//
// the e3 basis vector represents the origin of the projective plane
//
// the e12 basis bivector represents the horizon of the 2d plane
// containing all points at infinity
//
// Vec2dp: represents a point in 2d when z=1.0 with position at x*e1 + y*e2 in the
//         2d plane. The weight is z*e3.
//
// join(p,q) = wdg(p,q)
// meet(p,q) = rwdg(p,q) = !wdg(!a,!b) with !a as the complement of a, where a^!a = I_2dp
//
//  u:     1 |  e1 |  e2 |  e3 | e23 | e31 | e12 | e123
// !u:  e123 | e23 | e31 | e12 |  e1 |  e2 |  e3 |    1
//
//  right complement (a^(!a)=e123) and left complement ((!a)^a=e123) are identical in 2dp
//
// BiVec2dp: represents a line in 2d when created as join from two points p and q as
//           wdg(p, q) = p^q
//              l(p,q) = (p.y * q.z - p.z * q.y) * e23
//                     + (p.z * q.x - p.x * q.z) * e31
//                     + (p.x * q.y - p.y * q.x) * e12
//                     = l.x * e23 + l.y * e31 + l.z * e12
//
//           The line is created where the 2d projective plane through the origin e3 (!e3)
//           and the plane represented by the bivector (which goes through the origin
//           (0,0,0) of the representational space) intersect:
//
//           meet(l,!e3) = !wdg(!l,e3) = !wdg( (-l.x*e1 - l.y*e2 - l.z*e3), e3)
//                       = !( (-l.y)*e23 + (l.x)*e31 + (0.0)*e12 )
//                       = (l.y) * e1 + (-l.x) * e2
//
//           l = x*e23 + y*e31 + z*e12
//               (normal: x*e23 + y*e31)
//               (position: z*e12)

} // namespace hd::ga


namespace hd::ga::pga {

////////////////////////////////////////////////////////////////////////////////
// convenient pga specific type aliases
////////////////////////////////////////////////////////////////////////////////

// Vector2d: 2d vector of projective algebra storing only components x, y explicitly.
// The z component is assumed to be z = 0
template <typename T>
    requires(std::floating_point<T>)
struct Vector2d : public Vec2d<T> {

    using Vec2d<T>::Vec2d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec2d<T>::x;
    using Vec2d<T>::y;

    Vector2d(Vec2d<T> const& v) : Vec2d<T>(v) {};
};

// Point2d: 2d point of projective algebra storing only components x, y explicitly.
// The z component is assumed to by z=1
template <typename T>
    requires(std::floating_point<T>)
struct Point2d : public Vec2d<T> {

    using Vec2d<T>::Vec2d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec2d<T>::x;
    using Vec2d<T>::y;

    Point2d(Vec2d<T> const& v) : Vec2d<T>(v) {};
    Point2d(T x, T y) : Vec2d<T>(x, y) {};
};

// Point2dp: 2d point of projective algebra storing all three components x, y, z
//           explicitly
//
// a Point2dp is a Vec2dp, thus all operations defined for Vec2dp
// work directly for Point2dp - only deviations will be specified
template <typename T>
    requires(std::floating_point<T>)
struct Point2dp : public Vec2dp<T> {

    using Vec2dp<T>::Vec2dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec2dp<T>::x;
    using Vec2dp<T>::y;
    using Vec2dp<T>::z;

    Point2dp(Vec2dp<T> const& v) : Vec2dp<T>(v) {};
    Point2dp(Point2d<T> const& p) : Vec2dp<T>(p.x, p.y, T(1.0)) {};
    Point2dp(Vec2d<T> const& v) : Vec2dp<T>(v.x, v.y, T(1.0)) {};
};

// Line2d: 2d line of projective algebra storing all three components
//          explicitly as components x, y, z of a BiVec2dp
//
// a Line2d is a BiVec2dp, thus all operations defined for BiVec2dp
// work directly for Line2d - only deviations will be specified
template <typename T>
    requires(std::floating_point<T>)
struct Line2d : public BiVec2dp<T> {

    using BiVec2dp<T>::BiVec2dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using BiVec2dp<T>::x;
    using BiVec2dp<T>::y;
    using BiVec2dp<T>::z;

    Line2d() = default;
    Line2d(T x, T y, T z) : BiVec2dp<T>(x, y, z) {};
    Line2d(BiVec2dp<T> const& b) : BiVec2dp<T>(b) {};
    Line2d(Point2d<T> const& p, Point2d<T> const& q) :
        // Line2d = wdg(p,q), but wdg() cannot be used here to avoid circular dependency
        BiVec2dp<T>(p.y - q.y, q.x - p.x, p.x * q.y - p.y * q.x) {};
    Line2d(Point2d<T> const& p, Vec2d<T> const& v) :
        // Line2d constructed from a point and a direction vector
        BiVec2dp<T>(-v.y, v.x, p.x * v.y - p.y * v.x) {};
};

} // namespace hd::ga::pga