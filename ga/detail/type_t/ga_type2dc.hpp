#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type_tags.hpp"

#include "ga_bvec6_t.hpp"  // bivector
#include "ga_scalar_t.hpp" // scalar, pseudoscalar
#include "ga_vec4_t.hpp"   // vector

#include "ga_mvec16_t.hpp" // full multivector
#include "ga_mvec8_t.hpp"  // even and odd multivectors

#include "ga_mvec2_t.hpp" // for DualNum2dc

/////////////////////////////////////////////////////////////////////////////////////////
// consistent type and grade definitions in namespace hd::ga for cga2dc
//
// cga2dc: conformal geometric algebra of 2D Euclidean space (16 components).
// Basis vectors: e1, e2 Euclidean; e3 projective (w); e4 round (u). The null
// pair satisfies e3*e3 = e4*e4 = 0 and e3.e4 = -1 (non-orthogonal metric;
// G(3,1,0) after diagonalization).
//
// Field-name reuse (positional): Vec2dc = Vec4_t stores the (x, y, w, u)
// components in the fields named x, y, z, w — same pattern as pga2dp, where
// Vec2dp = Vec3_t stores the homogeneous coordinate in field z. TriVec2dc
// covers the basis trivectors (e314, e324, e124, e321) in fields (x, y, z, w);
// its first three components carry a lifted pga2dp line (pga2dp bivectors ∧
// e4), the w field the pure-Euclidean blade e321 (the pga2dp pseudoscalar).
//
// The conformal geometric object aliases (round point, point pair, circle and
// the flat counterparts) are added together with the cga operations layer.
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar2dc = Scalar_t<T, scalar2dc_tag>;
template <typename T> using Vec2dc = Vec4_t<T, vec2dc_tag>;
template <typename T> using BiVec2dc = BVec6_t<T, bivec2dc_tag>;
template <typename T> using TriVec2dc = Vec4_t<T, trivec2dc_tag>;
template <typename T> using PScalar2dc = Scalar_t<T, pscalar2dc_tag>;
template <typename T> using MVec2dc_E = MVec8_t<T, mvec2dc_e_tag>;
template <typename T> using MVec2dc_U = MVec8_t<T, mvec2dc_u_tag>;
template <typename T> using MVec2dc = MVec16_t<T, mvec2dc_tag>;

// for the geometric norm as homogeneous quantity,
// which has a scalar part (c0 component) and a pseudoscalar part (c1 component)
template <typename T> using DualNum2dc = MVec2_t<T, dual_number2dc_tag>;

// return the grades of the basic types

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Scalar2dc<T>)
{
    return 0;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Vec2dc<T> const&)
{
    return 1;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] BiVec2dc<T> const&)
{
    return 2;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] TriVec2dc<T> const&)
{
    return 3;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] PScalar2dc<T>)
{
    return 4;
}

// return the regressive grades of the basic types: rgr(arg) = n - gr(arg)
// with n being the dimension of the space
// => gr(arg) + rgr(arg) = n
//
// (Lengyel calls that anti-grade of an object, i.e. the number of dimensions
//  of the space NOT populated by the object)

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] Scalar2dc<T>)
{
    return 4; // 4 - 0
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] Vec2dc<T> const&)
{
    return 3; // 4 - 1
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] BiVec2dc<T> const&)
{
    return 2; // 4 - 2
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] TriVec2dc<T> const&)
{
    return 1; // 4 - 3
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] PScalar2dc<T>)
{
    return 0; // 4 - 4
}

} // namespace hd::ga
