#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type_tags.hpp"

#include "ga_scalar_t.hpp"
#include "ga_vec3_t.hpp"

#include "ga_mvec4_t.hpp"
#include "ga_mvec8_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

template <typename T> using Scalar2dp = Scalar_t<T, scalar2dp_tag>;
template <typename T> using Vec2dp = Vec3_t<T, vec2dp_tag>;
template <typename T> using BiVec2dp = Vec3_t<T, bivec2dp_tag>;
template <typename T> using PScalar2dp = Scalar_t<T, pscalar2dp_tag>;
template <typename T> using MVec2dp_E = MVec4_t<T, mvec2dp_e_tag>;
template <typename T> using MVec2dp_U = MVec4_t<T, mvec2dp_u_tag>;
template <typename T> using MVec2dp = MVec8_t<T, mvec2dp_tag>;

// for the geometric norm as homogeneous quantity,
// which has a scalar part (c0 component) and a pseudoscalar part (c1 component)
template <typename T> using DualNum2dp = MVec2_t<T, dual_number2dp_tag>;


// no support for specific interpretations or optimization yet
//
// Vec2dp: represents as direction when z=0.0 (or a point at infinity in that direction)
//         the direction is x*e1 + y*e2
//
// the e3 basis vector represents the origin of the projective plane
//
// the e12 basis bivector represents the horizon of the 2d plane
// containing all points at infinity
//
// Vec2dp: represents a point in 2d when z=1.0 with position at x*e1 + y*e2 in the
//         2d plane. The weight is z*e3.
//
// join(p,q) = wdg(p,q)
// meet(p,q) = !wdg(!a,!b) with !a as the complement of a, where a^!a = I_2dp
//
//  u:     1 |  e1 |  e2 |  e3 | e23 | e31 | e12 | e123
// !u:  e123 | e23 | e31 | e12 |  e1 |  e2 |  e3 |    1
//
//  right complement (a^(!a)=e123) and left complement ((!a)^a=e123) are identical in 2dp
//
// BiVec2dp: represents a line in 2d when created as join from two points p and q as
//           wdg(p, q) = p^q
//              l(p,q) = (p.y * q.z - p.z * q.y) * e23
//                     + (p.z * q.x - p.x * q.z) * e31
//                     + (p.x * q.y - p.y * q.x) * e12
//                     = l.x * e23 + l.y * e31 + l.z * e12
//
//           The line is created where the 2d projective plane through the origin e3 (!e3)
//           and the plane represented by the bivector (which goes through the origin
//           (0,0,0) of the representational space) intersect:
//
//           meet(l,!e3) = !wdg(!l,e3) = !wdg( (-l.x*e1 - l.y*e2 - l.z*e3), e3)
//                       = !( (-l.y)*e23 + (l.x)*e31 + (0.0)*e12 )
//                       = (l.y) * e1 + (-l.x) * e2
//
//           l = x*e23 + y*e31 + z*e12
//               (normal: x*e23 + y*e31)
//               (position: z*e12)

} // namespace hd::ga