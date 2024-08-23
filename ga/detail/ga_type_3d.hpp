#pragma once

// author: Daniel Hug, 2024

#include "ga_type_tags.hpp"

#include "ga_scalar_t.hpp"
#include "ga_vec3_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Vec3d = Vec3_t<T, vec3d_tag>;
template <typename T> using BiVec3d = Vec3_t<T, bivec3d_tag>;
template <typename T> using PScalar3d = Scalar_t<T, pscalar3d_tag>;

} // namespace hd::ga