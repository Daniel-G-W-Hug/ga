#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type_tags.hpp"

#include "ga_bvec6_t.hpp"  // bivector
#include "ga_scalar_t.hpp" // scalar, pseudoscalar
#include "ga_vec4_t.hpp"   // vector

#include "ga_mvec16_t.hpp" // full multivector
#include "ga_mvec8_t.hpp"  // even and odd multivectors

#include "ga_mvec2_t.hpp" // for DualNum4ds

/////////////////////////////////////////////////////////////////////////////////////////
// consistent type and grade definitions in namespace hd::ga for pga4ds
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar4ds = Scalar_t<T, scalar4ds_tag>;
template <typename T> using Vec4ds = Vec4_t<T, vec4ds_tag>;
template <typename T> using BiVec4ds = BVec6_t<T, bivec4ds_tag>;
template <typename T> using TriVec4ds = Vec4_t<T, trivec4ds_tag>;
template <typename T> using PScalar4ds = Scalar_t<T, pscalar4ds_tag>;
template <typename T> using MVec4ds_E = MVec8_t<T, mvec4ds_e_tag>;
template <typename T> using MVec4ds_U = MVec8_t<T, mvec4ds_u_tag>;
template <typename T> using MVec4ds = MVec16_t<T, mvec4ds_tag>;

// for the geometric norm as homogeneous quantity,
// which has a scalar part (c0 component) and a pseudoscalar part (c1 component)
template <typename T> using DualNum4ds = MVec2_t<T, dual_number4ds_tag>;

// return the grades of the basic types

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Scalar4ds<T>)
{
    return 0;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Vec4ds<T> const&)
{
    return 1;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] BiVec4ds<T> const&)
{
    return 2;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] TriVec4ds<T> const&)
{
    return 3;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] PScalar4ds<T>)
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
constexpr size_t rgr([[maybe_unused]] Scalar4ds<T>)
{
    return 4; // 4 - 0
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] Vec4ds<T> const&)
{
    return 3; // 4 - 1
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] BiVec4ds<T> const&)
{
    return 2; // 2 - 2
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] TriVec4ds<T> const&)
{
    return 1; // 4 - 3
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t rgr([[maybe_unused]] PScalar4ds<T>)
{
    return 0; // 4 - 4
}

} // namespace hd::ga
