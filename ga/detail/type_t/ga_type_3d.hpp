#pragma once

// author: Daniel Hug, 2024

#include "ga_type_tags.hpp"

#include "ga_scalar_t.hpp"
#include "ga_vec3_t.hpp"

#include "ga_mvec4_t.hpp"
#include "ga_mvec8_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar3d = Scalar_t<T, scalar3d_tag>;
template <typename T> using Vec3d = Vec3_t<T, vec3d_tag>;
template <typename T> using BiVec3d = Vec3_t<T, bivec3d_tag>;
template <typename T> using PScalar3d = Scalar_t<T, pscalar3d_tag>;
template <typename T> using MVec3d_E = MVec4_t<T, mvec3d_e_tag>;
template <typename T> using MVec3d_U = MVec4_t<T, mvec3d_u_tag>;
template <typename T> using MVec3d = MVec8_t<T, mvec3d_tag>;

} // namespace hd::ga