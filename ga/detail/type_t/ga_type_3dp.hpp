#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type_tags.hpp"

#include "ga_bvec6_t.hpp"
#include "ga_scalar_t.hpp"
#include "ga_vec4_t.hpp"

#include "ga_mvec16_t.hpp"
#include "ga_mvec8_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar3dp = Scalar_t<T, scalar3dp_tag>;
template <typename T> using Vec3dp = Vec4_t<T, vec3dp_tag>;
template <typename T> using BiVec3dp = BVec6_t<T, bivec3dp_tag>;
template <typename T> using TriVec3dp = Vec4_t<T, trivec3dp_tag>;
template <typename T> using PScalar3dp = Scalar_t<T, pscalar3dp_tag>;
template <typename T> using MVec3dp_E = MVec8_t<T, mvec3dp_e_tag>;
template <typename T> using MVec3dp_U = MVec8_t<T, mvec3dp_u_tag>;
template <typename T> using MVec3dp = MVec16_t<T, mvec3dp_tag>;

// for the geometric norm as homogeneous quantity,
// which has a scalar part (c0 component) and a pseudoscalar part (c1 component)
template <typename T> using DualNum3dp = MVec2_t<T, dual_number3dp_tag>;


} // namespace hd::ga