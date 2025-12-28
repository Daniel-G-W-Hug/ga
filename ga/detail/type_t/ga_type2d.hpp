#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type_tags.hpp"

#include "ga_scalar_t.hpp"
#include "ga_vec2_t.hpp"

#include "ga_mvec2_t.hpp"
#include "ga_mvec4_t.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// consistent type and grade definitions (ega2d)
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar2d = Scalar_t<T, scalar2d_tag>;
template <typename T> using Vec2d = Vec2_t<T, vec2d_tag>;
template <typename T> using PScalar2d = Scalar_t<T, pscalar2d_tag>;
template <typename T> using MVec2d_E = MVec2_t<T, mvec2d_e_tag>;
template <typename T> using MVec2d = MVec4_t<T, mvec2d_tag>;

// return the grades of the basic types

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Scalar2d<T>)
{
    return 0;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] Vec2d<T> const&)
{
    return 1;
}

template <typename T>
    requires(numeric_type<T>)
constexpr size_t gr([[maybe_unused]] PScalar2d<T>)
{
    return 2;
}

} // namespace hd::ga