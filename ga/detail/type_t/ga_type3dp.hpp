#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type_tags.hpp"

#include "ga_bvec6_t.hpp"  // bivector
#include "ga_scalar_t.hpp" // scalar, pseudoscalar
#include "ga_vec4_t.hpp"   // vector

#include "ga_mvec16_t.hpp" // full multivector
#include "ga_mvec8_t.hpp"  // even and odd multivectors

#include "ga_mvec2_t.hpp" // for DualNum3dp
#include "ga_type3d.hpp"  // for Vec3d<T>

/////////////////////////////////////////////////////////////////////////////////////////
// consistent type definitions in namespace hd::ga for pga3dp
// convenience type aliases in namespace hd::ga::pga (needed there for name resolution)
//     - Vector3d
//     - Point3d
//     - Point3dp
//     - Line3d
//     - Plane3d
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar3dp = Scalar_t<T, scalar3dp_tag>;
template <typename T> using Vec3dp = Vec4_t<T, vec3dp_tag>;
template <typename T> using BiVec3dp = BVec6_t<T, bivec3dp_tag>;
template <typename T> using TriVec3dp = Vec4_t<T, trivec3dp_tag>;
template <typename T> using PScalar3dp = Scalar_t<T, pscalar3dp_tag>;
template <typename T> using MVec3dp_E = MVec8_t<T, mvec3dp_e_tag>;
template <typename T> using MVec3dp_U = MVec8_t<T, mvec3dp_u_tag>;
template <typename T> using MVec3dp = MVec16_t<T, mvec3dp_tag>;

// for the geometric norm as homogeneous quantity,
// which has a scalar part (c0 component) and a pseudoscalar part (c1 component)
template <typename T> using DualNum3dp = MVec2_t<T, dual_number3dp_tag>;


} // namespace hd::ga


namespace hd::ga::pga {
////////////////////////////////////////////////////////////////////////////////
// convenient pga specific type aliases
////////////////////////////////////////////////////////////////////////////////

// Vector3d: 3d vector of projective algebra storing only components x, y, z explicitly.
// The w component is assumed to be w = 0.0
template <typename T>
    requires(std::floating_point<T>)
struct Vector3d : public Vec3d<T> {

    using Vec3d<T>::Vec3d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec3d<T>::x;
    using Vec3d<T>::y;
    using Vec3d<T>::z;

    Vector3d(Vec3d<T> const& v) : Vec3d<T>(v) {};
};

// Point3d: 3d point of projective algebra storing only components x, y, z explicitly.
// The w component is assumed to by w = 1.0
template <typename T>
    requires(std::floating_point<T>)
struct Point3d : public Vec3d<T> {

    using Vec3d<T>::Vec3d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec3d<T>::x;
    using Vec3d<T>::y;
    using Vec3d<T>::z;

    Point3d(Vec3d<T> const& v) : Vec3d<T>(v) {};
    Point3d(T x, T y, T z) : Vec3d<T>(x, y, z) {};
};

// Point3dp: 3d point of projective algebra storing all four components x, y, z, w
//           explicitly
//
// a Point3dp is a Vec3dp, thus all operations defined for Vec3dp
// work directly for Point3dp - only deviations will be specified
template <typename T>
    requires(std::floating_point<T>)
struct Point3dp : public Vec3dp<T> {

    using Vec3dp<T>::Vec3dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec3dp<T>::x;
    using Vec3dp<T>::y;
    using Vec3dp<T>::z;
    using Vec3dp<T>::w;

    Point3dp(Vec3dp<T> const& v) : Vec3dp<T>(v) {};
    Point3dp(Point3d<T> const& p) : Vec3dp<T>(p.x, p.y, p.z, T(1.0)) {};
    Point3dp(Vec3d<T> const& v) : Vec3dp<T>(v.x, v.y, v.z, T(1.0)) {};
};

// Line3d: 3d line of projective algebra storing all six components
//          explicitly as components vx, vy, vz, mx, my, mz of a BiVec3dp
//
// a Line3d is a BiVec3dp, thus all operations defined for BiVec3dp
// work directly for Line3d - only deviations will be specified
template <typename T>
    requires(std::floating_point<T>)
struct Line3d : public BiVec3dp<T> {

    using BiVec3dp<T>::BiVec3dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using BiVec3dp<T>::vx;
    using BiVec3dp<T>::vy;
    using BiVec3dp<T>::vz;
    using BiVec3dp<T>::mx;
    using BiVec3dp<T>::my;
    using BiVec3dp<T>::mz;

    Line3d() = default;
    Line3d(T vx, T vy, T vz, T mx, T my, T mz) : BiVec3dp<T>(vx, vy, vz, mx, my, mz) {};
    Line3d(BiVec3dp<T> const& b) : BiVec3dp<T>(b) {};
    Line3d(Vec3d<T> const& dir, BiVec3d<T> const& mom) :
        // direction vector dir and moment bivector mom must match, i.e. be perpendicular
        // to each other (direction << moment) == 0
        BiVec3dp<T>(dir.x, dir.y, dir.z, mom.x, mom.y, mom.z) {};
};

// Plane3d: 3d plane of projective algebra (a trivector in the modeling 4d space)
template <typename T>
    requires(std::floating_point<T>)
struct Plane3d : public TriVec3dp<T> {

    using TriVec3dp<T>::TriVec3dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using TriVec3dp<T>::x;
    using TriVec3dp<T>::y;
    using TriVec3dp<T>::z;
    using TriVec3dp<T>::w;

    Plane3d(TriVec3dp<T> const& t) : TriVec3dp<T>(t) {};
    // join functions must be accessible at template instantiation time
    Plane3d(Line3d<T> const& l, Point3d<T> const& p) : TriVec3dp<T>(join(l, p)) {};
    Plane3d(Point3d<T> const& p1, Point3d<T> const& p2, Point3d<T> const& p3) :
        TriVec3dp<T>(join(join(p1, p2), p3)) {};
};

} // namespace hd::ga::pga