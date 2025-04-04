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

template <typename T> using Scalar4d = Scalar_t<T, scalar4d_tag>;
template <typename T> using Vec4d = Vec4_t<T, vec4d_tag>;
template <typename T> using BiVec4d = BVec6_t<T, bivec4d_tag>;
template <typename T> using TriVec4d = Vec4_t<T, trivec4d_tag>;
template <typename T> using PScalar4d = Scalar_t<T, pscalar4d_tag>;
template <typename T> using MVec4d_E = MVec8_t<T, mvec4d_e_tag>;
template <typename T> using MVec4d_U = MVec8_t<T, mvec4d_u_tag>;
template <typename T> using MVec4d = MVec16_t<T, mvec4d_tag>;

} // namespace hd::ga