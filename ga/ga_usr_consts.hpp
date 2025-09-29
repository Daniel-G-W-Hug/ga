#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_value_t.hpp" // defines value_t either as float oder double

#include "detail/type_t/ga_type2d.hpp"
#include "detail/type_t/ga_type3d.hpp"

#include "detail/type_t/ga_type2dp.hpp"
#include "detail/type_t/ga_type3dp.hpp"

#include "detail/type_t/ga_mvec2d.hpp" // MVec2d<T> and subtypes
#include "detail/type_t/ga_mvec3d.hpp" // MVec3d<T> and subtypes

#include "detail/type_t/ga_mvec2dp.hpp" // MVec2dp<T> and subtypes
#include "detail/type_t/ga_mvec3dp.hpp" // MVec3dp<T> and subtypes

// provide constants for ega and pga in their respective namespaces

namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// 2d constants
/////////////////////////////////////////////////////////////////////////////////////////

auto const e1_2d = Vec2d<value_t>{1.0, 0.0};
auto const e2_2d = Vec2d<value_t>{0.0, 1.0};
auto const x_dir_2d = e1_2d;                // x-direction
auto const y_dir_2d = e2_2d;                // y-direction
auto const O_2d = Vec2d<value_t>{0.0, 0.0}; // origin

auto const e1_2d_mv = MVec2d<value_t>{e1_2d}; // e1_2d as multivector
auto const e2_2d_mv = MVec2d<value_t>{e2_2d}; // e2_2d as multivector

auto const I_2d = PScalar2d<value_t>(1.0);
auto const e12_2d = I_2d;
auto const I_2d_mv = MVec2d<value_t>{I_2d};     // I_2d as multivector
auto const I_2d_mv_e = MVec2d_E<value_t>{I_2d}; // I_2d as even grade multivector

/////////////////////////////////////////////////////////////////////////////////////////
// 3d constants
/////////////////////////////////////////////////////////////////////////////////////////

auto const e1_3d = Vec3d<value_t>{1.0, 0.0, 0.0};
auto const e2_3d = Vec3d<value_t>{0.0, 1.0, 0.0};
auto const e3_3d = Vec3d<value_t>{0.0, 0.0, 1.0};
auto const x_dir_3d = e1_3d;                     // x-direction
auto const y_dir_3d = e2_3d;                     // y-direction
auto const z_dir_3d = e3_3d;                     // z-direction
auto const O_3d = Vec3d<value_t>{0.0, 0.0, 0.0}; // origin

auto const e1_3d_mv = MVec3d<value_t>{e1_3d}; // e1_3d as multivector
auto const e2_3d_mv = MVec3d<value_t>{e2_3d}; // e2_3d as multivector
auto const e3_3d_mv = MVec3d<value_t>{e3_3d}; // e3_3d as multivector

auto const e23_3d = BiVec3d<value_t>{1.0, 0.0, 0.0};
auto const e31_3d = BiVec3d<value_t>{0.0, 1.0, 0.0};
auto const e12_3d = BiVec3d<value_t>{0.0, 0.0, 1.0};
auto const yz_3d = e23_3d; // yz- or 23-plane
auto const zx_3d = e31_3d; // zx- or 31-plane
auto const xy_3d = e12_3d; // xy- or 12-plane

auto const e23_3d_mv = MVec3d<value_t>{e23_3d};     // e23_3d as multivector
auto const e31_3d_mv = MVec3d<value_t>{e31_3d};     // e31_3d as multivector
auto const e12_3d_mv = MVec3d<value_t>{e12_3d};     // e12_3d as multivector
auto const e23_3d_mv_e = MVec3d_E<value_t>{e23_3d}; // e23_3d as even grade multivector
auto const e31_3d_mv_e = MVec3d_E<value_t>{e31_3d}; // e31_3d as even grade multivector
auto const e12_3d_mv_e = MVec3d_E<value_t>{e12_3d}; // e12_3d as even grade multivector

auto const I_3d = PScalar3d<value_t>(1.0);
auto const e123_3d = I_3d;
auto const I_3d_mv = MVec3d<value_t>{I_3d};     // I_3d as multivector
auto const I_3d_mv_u = MVec3d_U<value_t>{I_3d}; // I_3d as odd grade multivector

} // namespace hd::ga::ega


namespace hd::ga::pga {

/////////////////////////////////////////////////////////////////////////////////////////
// 2dp constants
/////////////////////////////////////////////////////////////////////////////////////////

auto const e1_2dp = Vec2dp<value_t>{1.0, 0.0, 0.0};
auto const e2_2dp = Vec2dp<value_t>{0.0, 1.0, 0.0};
auto const e3_2dp = Vec2dp<value_t>{0.0, 0.0, 1.0};
auto const O_2dp = e3_2dp; // origin

auto const x_dir_2dp = e1_2dp; // direction of x-axis towards infinity
auto const y_dir_2dp = e2_2dp; // direction of y-axis towards infinity

auto const e1_2dp_mv = MVec2dp<value_t>{e1_2dp}; // e1_2dp as multivector
auto const e2_2dp_mv = MVec2dp<value_t>{e2_2dp}; // e2_2dp as multivector
auto const e3_2dp_mv = MVec2dp<value_t>{e3_2dp}; // e3_2dp as multivector

auto const e23_2dp = BiVec2dp<value_t>{1.0, 0.0, 0.0};
auto const e31_2dp = BiVec2dp<value_t>{0.0, 1.0, 0.0};
auto const e12_2dp = BiVec2dp<value_t>{0.0, 0.0, 1.0};

// x-axis as bivector
//      x-axis = wdg(vec2dp{0, 0, 1}, vec2dp{1, 0, 1}) = bivt2dp(0, 1, 0)
// att(x-axis) = vec2dp(1, 0, 0) = x_dir_2dp = e1_2dp
auto const x_axis_2dp = BiVec2dp<value_t>{0.0, 1.0, 0.0};

// y-axis as bivector
//      y-axis = wdg(vec2dp{0, 0, 1}, vec2dp{0, 1, 1})
// att(y-axis) = vec2dp(0, 1, 0) = y_dir_2dp = e2_2dp
auto const y_axis_2dp = BiVec2dp<value_t>{-1.0, 0.0, 0.0};

auto const H_2dp = -e12_2dp; // horizon == line at infinity == cmpl(origin)

auto const e23_2dp_mv = MVec2dp<value_t>{e23_2dp}; // e23_2dp as multivector
auto const e31_2dp_mv = MVec2dp<value_t>{e31_2dp}; // e31_2dp as multivector
auto const e12_2dp_mv = MVec2dp<value_t>{e12_2dp}; // e12_2dp as multivector

auto const e23_2dp_mv_e = MVec2dp_E<value_t>{e23_2dp}; // e23_2dp as even multivector
auto const e31_2dp_mv_e = MVec2dp_E<value_t>{e31_2dp}; // e31_2dp as even multivector
auto const e12_2dp_mv_e = MVec2dp_E<value_t>{e12_2dp}; // e12_2dp as even multivector

auto const I_2dp = PScalar2dp<value_t>(1.0); // I_2dp = e3^e2^e1 = e321
auto const e321_2dp = I_2dp;
auto const I_2dp_mv = MVec2dp<value_t>{I_2dp};     // I_2dp as multivector
auto const I_2dp_mv_u = MVec2dp_U<value_t>{I_2dp}; // I_2dp as odd grade multivector

/////////////////////////////////////////////////////////////////////////////////////////
// 3dp constants
/////////////////////////////////////////////////////////////////////////////////////////

auto const e1_3dp = Vec3dp<value_t>{1.0, 0.0, 0.0, 0.0};
auto const e2_3dp = Vec3dp<value_t>{0.0, 1.0, 0.0, 0.0};
auto const e3_3dp = Vec3dp<value_t>{0.0, 0.0, 1.0, 0.0};
auto const e4_3dp = Vec3dp<value_t>{0.0, 0.0, 0.0, 1.0};
auto const O_3dp = e4_3dp; // origin

auto const x_dir_3dp = e1_3dp; // direction of x-axis towards infinity
auto const y_dir_3dp = e2_3dp; // direction of y-axis towards infinity
auto const z_dir_3dp = e3_3dp; // direction of z-axis towards infinity

auto const e1_3dp_mv = MVec3dp<value_t>{e1_3dp}; // e1_3dp as multivector
auto const e2_3dp_mv = MVec3dp<value_t>{e2_3dp}; // e2_3dp as multivector
auto const e3_3dp_mv = MVec3dp<value_t>{e3_3dp}; // e3_3dp as multivector
auto const e4_3dp_mv = MVec3dp<value_t>{e4_3dp}; // e4_3dp as multivector

auto const e41_3dp = BiVec3dp<value_t>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e42_3dp = BiVec3dp<value_t>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
auto const e43_3dp = BiVec3dp<value_t>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
auto const e23_3dp = BiVec3dp<value_t>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
auto const e31_3dp = BiVec3dp<value_t>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
auto const e12_3dp = BiVec3dp<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

// x-axis as bivector
//      x-axis = wdg(vec3dp{0, 0, 0, 1}, vec3dp{1, 0, 0, 1}) = bivec3dp(1, 0, 0, 0, 0, 0)
// att(x-axis) = Vec3dp(1, 0, 0, 0) = x_dir_3dp = e1_3dp
auto const x_axis_3dp = e41_3dp;

// y-axis as bivector
//      y-axis = wdg(vec3dp{0, 0, 0, 1}, vec3dp{0, 1, 0, 1}) = bivec3dp(0, 1, 0, 0, 0, 0)
// att(y-axis) = Vec3dp(0, 1, 0, 0) = y_dir_3dp = e2_3dp
auto const y_axis_3dp = e42_3dp;

// z-axis as bivector
//      z-axis = wdg(vec3dp{0, 0, 0, 1}, vec3dp{0, 0, 1, 1}) = bivec3dp(0, 0, 1, 0, 0, 0)
// att(z-axis) = Vec3dp(0, 0, 1, 0) = z_dir_3dp = e3_3dp
auto const z_axis_3dp = e43_3dp;

auto const e41_3dp_mv = MVec3dp<value_t>{e41_3dp}; // e41_3dp as multivector
auto const e42_3dp_mv = MVec3dp<value_t>{e42_3dp}; // e42_3dp as multivector
auto const e43_3dp_mv = MVec3dp<value_t>{e43_3dp}; // e43_3dp as multivector
auto const e23_3dp_mv = MVec3dp<value_t>{e23_3dp}; // e23_3dp as multivector
auto const e31_3dp_mv = MVec3dp<value_t>{e31_3dp}; // e31_3dp as multivector
auto const e12_3dp_mv = MVec3dp<value_t>{e12_3dp}; // e12_3dp as multivector

auto const e23_3dp_mv_e = MVec3dp_E<value_t>{e23_3dp}; // e23_3dp as even multivector
auto const e31_3dp_mv_e = MVec3dp_E<value_t>{e31_3dp}; // e31_3dp as even multivector
auto const e12_3dp_mv_e = MVec3dp_E<value_t>{e12_3dp}; // e12_3dp as even multivector

auto const e423_3dp = TriVec3dp<value_t>{1.0, 0.0, 0.0, 0.0};
auto const e431_3dp = TriVec3dp<value_t>{0.0, 1.0, 0.0, 0.0};
auto const e412_3dp = TriVec3dp<value_t>{0.0, 0.0, 1.0, 0.0};
auto const e321_3dp = TriVec3dp<value_t>{0.0, 0.0, 0.0, 1.0};
auto const yz_3dp = e423_3dp; // yz- or 23-plane (bound bivector: origin ^ e23 = e4 ^ e23)
auto const zx_3dp = e431_3dp; // zx- or 31-plane (bound bivector: origin ^ e31 = e4 ^ e31)
auto const xy_3dp = e412_3dp; // xy- or 12-plane (bound bivector: origin ^ e12 = e4 ^ e12)
auto const H_3dp = e321_3dp;  // horizon == plane at infinity == rcmpl(origin)

auto const e423_3dp_mv = MVec3dp<value_t>{e423_3dp}; // e423_3dp as multivector
auto const e431_3dp_mv = MVec3dp<value_t>{e431_3dp}; // e431_3dp as multivector
auto const e412_3dp_mv = MVec3dp<value_t>{e412_3dp}; // e412_3dp as multivector
auto const e321_3dp_mv = MVec3dp<value_t>{e321_3dp}; // e321_3dp as multivector

auto const I_3dp = PScalar3dp<value_t>(1.0); // I_3dp = e1^e2^e3^e4 = e1234
auto const e1234_3dp = I_3dp;
auto const I_3dp_mv = MVec3dp<value_t>{I_3dp};     // I_3dp as multivector
auto const I_3dp_mv_e = MVec3dp_E<value_t>{I_3dp}; // I_3dp as even grade multivector

} // namespace hd::ga::pga