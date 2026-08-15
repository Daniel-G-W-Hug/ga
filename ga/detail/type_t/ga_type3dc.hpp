#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type_tags.hpp"

#include "ga_bvec10_t.hpp" // bivector and trivector (10 components each)
#include "ga_scalar_t.hpp" // scalar, pseudoscalar
#include "ga_vec5_t.hpp"   // vector and quadvector (5 components each)

#include "ga_mvec16_t.hpp" // even and uneven (odd) multivectors
#include "ga_mvec32_t.hpp" // full multivector

#include "ga_mvec2_t.hpp" // for DualNum3dc

/////////////////////////////////////////////////////////////////////////////////////////
// consistent type and grade definitions in namespace hd::ga for cga3dc
//
// cga3dc: conformal geometric algebra of 3D Euclidean space (32 components).
// Basis vectors: e1, e2, e3 Euclidean; e4 projective (w, same slot as pga3dp);
// e5 round (u). The null pair satisfies e4*e4 = e5*e5 = 0 and e4.e5 = -1
// (non-orthogonal metric; G(4,1,0) after diagonalization).
//
// Field-name reuse (positional):
// - Vec3dc = Vec5_t carries (x, y, z, w, u) directly.
// - QuadVec3dc = Vec5_t stores the five basis quadvectors positionally in the
//   label order given by the omitted basis vector:
//   (x, y, z, w, u) -> (e4235, e4315, e4125, e3215, e1234).
// - BiVec3dc = BVec10_t with fields (vx, vy, vz, mx, my, mz, px, py, pz, pw)
//   on the blades (e41, e42, e43, e23, e31, e12, e15, e25, e35, e45): the
//   v/m fields carry a lifted pga3dp line verbatim, the p fields the round
//   position part, pw the blade e45.
// - TriVec3dc = BVec10_t on the blades (e415, e425, e435, e235, e315, e125,
//   e423, e431, e412, e321): a conformal LINE lands positionally in vx..mz,
//   the p fields carry the flat-point part, pw the pure-Euclidean blade e321.
//
// The conformal geometric object aliases (round point, dipole, circle, sphere
// and the flat counterparts) are added together with the cga3dc operations
// layer.
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar3dc = Scalar_t<T, scalar3dc_tag>;
template <typename T> using Vec3dc = Vec5_t<T, vec3dc_tag>;
template <typename T> using BiVec3dc = BVec10_t<T, bivec3dc_tag>;
template <typename T> using TriVec3dc = BVec10_t<T, trivec3dc_tag>;
template <typename T> using QuadVec3dc = Vec5_t<T, quadvec3dc_tag>;
template <typename T> using PScalar3dc = Scalar_t<T, pscalar3dc_tag>;
template <typename T> using MVec3dc_E = MVec16_t<T, mvec3dc_e_tag>;
template <typename T> using MVec3dc_U = MVec16_t<T, mvec3dc_u_tag>;
template <typename T> using MVec3dc = MVec32_t<T, mvec3dc_tag>;

// for the geometric norm as homogeneous quantity,
// which has a scalar part (c0 component) and a pseudoscalar part (c1 component)
template <typename T> using DualNum3dc = MVec2_t<T, dual_number3dc_tag>;

// return the grades of the basic types

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Scalar3dc<T>)
{
    return 0;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Vec3dc<T> const&)
{
    return 1;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] BiVec3dc<T> const&)
{
    return 2;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] TriVec3dc<T> const&)
{
    return 3;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] QuadVec3dc<T> const&)
{
    return 4;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] PScalar3dc<T>)
{
    return 5;
}

// return the regressive grades of the basic types: rgr(arg) = n - gr(arg)
// with n being the dimension of the space
// => gr(arg) + rgr(arg) = n
//
// (Lengyel calls that anti-grade of an object, i.e. the number of dimensions
//  of the space NOT populated by the object)

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] Scalar3dc<T>)
{
    return 5; // 5 - 0
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] Vec3dc<T> const&)
{
    return 4; // 5 - 1
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] BiVec3dc<T> const&)
{
    return 3; // 5 - 2
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] TriVec3dc<T> const&)
{
    return 2; // 5 - 3
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] QuadVec3dc<T> const&)
{
    return 1; // 5 - 4
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] PScalar3dc<T>)
{
    return 0; // 5 - 5
}

} // namespace hd::ga
