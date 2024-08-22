#pragma once

// author: Daniel Hug, 2024

#include "ga_scalar_t.hpp"
#include "ga_vec2_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions for easy use
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

struct vec2d_tag {};
struct pscalar2d_tag {};

template <typename T> using Vec2d = Vec2_t<T, vec2d_tag>;
template <typename T> using PScalar2d = Scalar_t<T, pscalar2d_tag>;

} // namespace hd::ga