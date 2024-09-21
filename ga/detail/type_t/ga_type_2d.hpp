#pragma once

// author: Daniel Hug, 2024

#include "ga_type_tags.hpp"

#include "ga_scalar_t.hpp"
#include "ga_vec2_t.hpp"

#include "ga_mvec2_t.hpp"
#include "ga_mvec4_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar2d = Scalar_t<T, scalar2d_tag>;
template <typename T> using Vec2d = Vec2_t<T, vec2d_tag>;
template <typename T> using PScalar2d = Scalar_t<T, pscalar2d_tag>;
template <typename T> using MVec2d_E = MVec2_t<T, mvec2d_e_tag>;
template <typename T> using MVec2d = MVec4_t<T, mvec2d_tag>;

} // namespace hd::ga