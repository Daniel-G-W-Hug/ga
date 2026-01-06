#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type_tags.hpp"

#include "ga_scalar_t.hpp"
#include "ga_vec3_t.hpp"

#include "ga_mvec4_t.hpp"
#include "ga_mvec8_t.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// consistent type and grade definitions (ega3d)
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar3d = Scalar_t<T, scalar3d_tag>;
template <typename T> using Vec3d = Vec3_t<T, vec3d_tag>;
template <typename T> using BiVec3d = Vec3_t<T, bivec3d_tag>;
template <typename T> using PScalar3d = Scalar_t<T, pscalar3d_tag>;
template <typename T> using MVec3d_E = MVec4_t<T, mvec3d_e_tag>;
template <typename T> using MVec3d_U = MVec4_t<T, mvec3d_u_tag>;
template <typename T> using MVec3d = MVec8_t<T, mvec3d_tag>;

// return the grades of the basic types

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Scalar3d<T>)
{
    return 0;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Vec3d<T> const&)
{
    return 1;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] BiVec3d<T> const&)
{
    return 2;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] PScalar3d<T>)
{
    return 3;
}

} // namespace hd::ga