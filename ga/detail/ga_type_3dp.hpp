#pragma once

// author: Daniel Hug, 2024

#include "ga_type_tags.hpp"

#include "ga_mvec6_t.hpp"
#include "ga_scalar_t.hpp"
#include "ga_vec4_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Vec3dp = Vec4_t<T, vec3dp_tag>;
template <typename T> using BiVec3dp = BVec6_t<T, bivec3dp_tag>;
template <typename T> using TriVec3dp = Vec4_t<T, trivec3dp_tag>;
template <typename T> using PScalar3dp = Scalar_t<T, pscalar3dp_tag>;

} // namespace hd::ga