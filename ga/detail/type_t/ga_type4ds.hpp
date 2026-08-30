#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type_tags.hpp"

#include "ga_bvec6_t.hpp"  // bivector
#include "ga_scalar_t.hpp" // scalar, pseudoscalar
#include "ga_vec4_t.hpp"   // vector

#include "ga_mvec16_t.hpp" // full multivector
#include "ga_mvec8_t.hpp"  // even and odd multivectors

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

// No rgr() here, deliberately. The anti-grade is the grade bookkeeping of the
// ANTI-PRODUCT family (rgpr, rdot, rcmt), which sta4ds does not have -- it carries
// only rwdg. ega2d/ega3d are in exactly the same position and define no rgr either;
// pga2dp/3dp and cga2dc/3dc, which do have the full family, define it.

} // namespace hd::ga
