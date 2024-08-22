#pragma once

// author: Daniel Hug, 2024

#include "ga_scalar_t.hpp"
#include "ga_vec3_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions for easy use
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

struct vec2dp_tag {};
struct bivec2dp_tag {};
struct pscalar2dp_tag {};

template <typename T> using Vec2dp = Vec3_t<T, vec2dp_tag>;
template <typename T> using BiVec2dp = Vec3_t<T, bivec2dp_tag>;
template <typename T> using PScalar2dp = Scalar_t<T, pscalar2dp_tag>;

} // namespace hd::ga