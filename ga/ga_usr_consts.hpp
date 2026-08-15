#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_value_t.hpp" // defines value_t either as float oder double

#include "detail/type_t/ga_type2d.hpp"
#include "detail/type_t/ga_type3d.hpp"

#include "detail/type_t/ga_type2dp.hpp"
#include "detail/type_t/ga_type3dp.hpp"

#include "detail/type_t/ga_type4ds.hpp"

#include "detail/type_t/ga_mvec2d.hpp" // MVec2d<T> and subtypes
#include "detail/type_t/ga_mvec3d.hpp" // MVec3d<T> and subtypes

#include "detail/type_t/ga_mvec2dp.hpp" // MVec2dp<T> and subtypes
#include "detail/type_t/ga_mvec3dp.hpp" // MVec3dp<T> and subtypes

#include "detail/type_t/ga_mvec2dc.hpp" // MVec2dc<T> and subtypes
#include "detail/type_t/ga_mvec3dc.hpp" // MVec3dc<T> and subtypes
#include "detail/type_t/ga_mvec4ds.hpp" // MVec4ds<T> and subtypes

#include <mdspan>

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

auto const one_2d = Scalar2d<value_t>(1.0);         // Basis element 1
auto const one_2d_mv = MVec2d<value_t>{one_2d};     // 1 as multivector
auto const one_2d_mv_e = MVec2d_E<value_t>{one_2d}; // 1 as even-grade multivector

auto const I_2d = PScalar2d<value_t>(1.0);      // Basis element I_2d
auto const e12_2d = I_2d;                       // (pseudoscalar)
auto const I_2d_mv = MVec2d<value_t>{I_2d};     // I_2d as multivector
auto const I_2d_mv_e = MVec2d_E<value_t>{I_2d}; // I_2d as even-grade multivector

// extended metric for ega2d (4×4 matrix, 16 elements)
// Basis: 1, e1, e2, e12
inline constexpr std::array<int, 16> ega2d_metric = {
    1, 0, 0, 0, // 1
    0, 1, 0, 0, // e1
    0, 0, 1, 0, // e2
    0, 0, 0, 1  // e12
};

// Accessor function for convenient 2D indexing: G[i,j]
inline constexpr auto ega2d_metric_view()
{
    return std::mdspan<int const, std::extents<size_t, 4, 4>>{ega2d_metric.data()};
}

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
auto const e23_3d_mv_e = MVec3d_E<value_t>{e23_3d}; // e23_3d as even-grade multivector
auto const e31_3d_mv_e = MVec3d_E<value_t>{e31_3d}; // e31_3d as even-grade multivector
auto const e12_3d_mv_e = MVec3d_E<value_t>{e12_3d}; // e12_3d as even-grade multivector

auto const one_3d = Scalar3d<value_t>(1.0);         // Basis element 1
auto const one_3d_mv = MVec3d<value_t>{one_3d};     // 1 as multivector
auto const one_3d_mv_e = MVec3d_E<value_t>{one_3d}; // 1 as even-grade multivector

auto const I_3d = PScalar3d<value_t>(1.0);      // Basis element I_3d
auto const e123_3d = I_3d;                      // (pseudoscalar)
auto const I_3d_mv = MVec3d<value_t>{I_3d};     // I_3d as multivector
auto const I_3d_mv_u = MVec3d_U<value_t>{I_3d}; // I_3d as odd-grade multivector

// extended metric for ega3d (8×8 matrix, 64 elements)
// Basis: 1, e1, e2, e3, e23, e31, e12, e123
inline constexpr std::array<int, 64> ega3d_metric = {
    1, 0, 0, 0, 0, 0, 0, 0, // 1
    0, 1, 0, 0, 0, 0, 0, 0, // e1
    0, 0, 1, 0, 0, 0, 0, 0, // e2
    0, 0, 0, 1, 0, 0, 0, 0, // e3
    0, 0, 0, 0, 1, 0, 0, 0, // e23
    0, 0, 0, 0, 0, 1, 0, 0, // e31
    0, 0, 0, 0, 0, 0, 1, 0, // e12
    0, 0, 0, 0, 0, 0, 0, 1  // e123
};

// Accessor function for convenient 2D indexing: G[i,j]
inline constexpr auto ega3d_metric_view()
{
    return std::mdspan<int const, std::extents<size_t, 8, 8>>{ega3d_metric.data()};
}

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

auto const e31_2dp = BiVec2dp<value_t>{1.0, 0.0, 0.0}; // = wdg(e3,e1)
auto const e32_2dp = BiVec2dp<value_t>{0.0, 1.0, 0.0}; // = wdg(e3,e2)
auto const e12_2dp = BiVec2dp<value_t>{0.0, 0.0, 1.0}; // = wdg(e1,e2)

// x-axis as bivector
//      x-axis = wdg(vec2dp{0, 0, 1}, vec2dp{1, 0, 1}) = wdg(e3,e1) = bivt2dp(1, 0, 0)
// att(x-axis) = vec2dp(1, 0, 0) = x_dir_2dp = e1_2dp
auto const x_axis_2dp = BiVec2dp<value_t>{1.0, 0.0, 0.0}; // = e31_3dp

// y-axis as bivector
//      y-axis = wdg(vec2dp{0, 0, 1}, vec2dp{0, 1, 1}) = wdg(e3,e2) = bivt2dp(0, 1, 0)
// att(y-axis) = vec2dp(0, 1, 0) = y_dir_2dp = e2_2dp
auto const y_axis_2dp = BiVec2dp<value_t>{0.0, 1.0, 0.0}; // = e32_3dp

auto const H_2dp = -e12_2dp; // horizon == line at infinity == cmpl(origin)

auto const e31_2dp_mv = MVec2dp<value_t>{e31_2dp}; // e31_2dp as multivector
auto const e32_2dp_mv = MVec2dp<value_t>{e32_2dp}; // e32_2dp as multivector
auto const e12_2dp_mv = MVec2dp<value_t>{e12_2dp}; // e12_2dp as multivector

auto const e31_2dp_mv_e = MVec2dp_E<value_t>{e31_2dp}; // e31_2dp as even multivector
auto const e32_2dp_mv_e = MVec2dp_E<value_t>{e32_2dp}; // e32_2dp as even multivector
auto const e12_2dp_mv_e = MVec2dp_E<value_t>{e12_2dp}; // e12_2dp as even multivector

auto const one_2dp = Scalar2dp<value_t>(1.0);          // Basis element 1
auto const one_2dp_mv = MVec2dp<value_t>{one_2dp};     // 1 as multivector
auto const one_2dp_mv_e = MVec2dp_E<value_t>{one_2dp}; // 1 as even-grade multivector

auto const I_2dp = PScalar2dp<value_t>(1.0);       // I_2dp = e3^e2^e1 = e321
auto const e321_2dp = I_2dp;                       // (pseudoscalar)
auto const I_2dp_mv = MVec2dp<value_t>{I_2dp};     // I_2dp as multivector
auto const I_2dp_mv_u = MVec2dp_U<value_t>{I_2dp}; // I_2dp as odd-grade multivector

// extended metric for pga2dp (8×8 matrix, 64 elements)
// Basis: 1, e1, e2, e3, e31, e32, e12, e321
inline constexpr std::array<int, 64> pga2dp_metric = {
    1, 0, 0, 0, 0, 0, 0, 0, // 1
    0, 1, 0, 0, 0, 0, 0, 0, // e1
    0, 0, 1, 0, 0, 0, 0, 0, // e2
    0, 0, 0, 0, 0, 0, 0, 0, // e3
    0, 0, 0, 0, 0, 0, 0, 0, // e31
    0, 0, 0, 0, 0, 0, 0, 0, // e32
    0, 0, 0, 0, 0, 0, 1, 0, // e12
    0, 0, 0, 0, 0, 0, 0, 0  // e321
};

// Accessor function for convenient 2D indexing: G[i,j]
inline constexpr auto pga2dp_metric_view()
{
    return std::mdspan<int const, std::extents<size_t, 8, 8>>{pga2dp_metric.data()};
}

// Regressive extended metric for pga2dp (8×8 matrix, 64 elements)
// Basis: 1, e1, e2, e3, e31, e32, e12, e321
inline constexpr std::array<int, 64> pga2dp_rmetric = {
    0, 0, 0, 0, 0, 0, 0, 0, // 1
    0, 0, 0, 0, 0, 0, 0, 0, // e1
    0, 0, 0, 0, 0, 0, 0, 0, // e2
    0, 0, 0, 1, 0, 0, 0, 0, // e3
    0, 0, 0, 0, 1, 0, 0, 0, // e31
    0, 0, 0, 0, 0, 1, 0, 0, // e32
    0, 0, 0, 0, 0, 0, 0, 0, // e12
    0, 0, 0, 0, 0, 0, 0, 1  // e321
};

// Accessor function for convenient 2D indexing: Ḡ[i,j]
inline constexpr auto pga2dp_rmetric_view()
{
    return std::mdspan<int const, std::extents<size_t, 8, 8>>{pga2dp_rmetric.data()};
}

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

auto const e41_3dp = BiVec3dp<value_t>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // = wdg(e4,e1)
auto const e42_3dp = BiVec3dp<value_t>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}; // = wdg(e4,e2)
auto const e43_3dp = BiVec3dp<value_t>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0}; // = wdg(e4,e3)
auto const e23_3dp = BiVec3dp<value_t>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}; // = wdg(e2,e3)
auto const e31_3dp = BiVec3dp<value_t>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0}; // = wdg(e3,e1)
auto const e12_3dp = BiVec3dp<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}; // = wdg(e1,e2)

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

auto const e41_3dp_mv_e = MVec3dp_E<value_t>{e41_3dp}; // e41_3dp as even multivector
auto const e42_3dp_mv_e = MVec3dp_E<value_t>{e42_3dp}; // e42_3dp as even multivector
auto const e43_3dp_mv_e = MVec3dp_E<value_t>{e43_3dp}; // e43_3dp as even multivector
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
auto const H_3dp = e321_3dp;  // horizon == plane at infinity == r_cmpl(origin)

auto const e423_3dp_mv = MVec3dp<value_t>{e423_3dp}; // e423_3dp as multivector
auto const e431_3dp_mv = MVec3dp<value_t>{e431_3dp}; // e431_3dp as multivector
auto const e412_3dp_mv = MVec3dp<value_t>{e412_3dp}; // e412_3dp as multivector
auto const e321_3dp_mv = MVec3dp<value_t>{e321_3dp}; // e321_3dp as multivector

auto const one_3dp = Scalar3dp<value_t>(1.0);          // Basis element 1
auto const one_3dp_mv = MVec3dp<value_t>{one_3dp};     // 1 as multivector
auto const one_3dp_mv_e = MVec3dp_E<value_t>{one_3dp}; // 1 as even-grade multivector

auto const I_3dp = PScalar3dp<value_t>(1.0);       // I_3dp = e1^e2^e3^e4 = e1234
auto const e1234_3dp = I_3dp;                      // (pseudoscalar)
auto const I_3dp_mv = MVec3dp<value_t>{I_3dp};     // I_3dp as multivector
auto const I_3dp_mv_e = MVec3dp_E<value_t>{I_3dp}; // I_3dp as even-grade multivector

// extended metric for pga3dp (16×16 matrix, 256 elements)
// Basis: 1, e1, e2, e3, e4, e41, e42, e43, e23, e31, e12, e423, e431, e412, e321, e1234
inline constexpr std::array<int, 256> pga3dp_metric = {
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 1
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e1
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e2
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e3
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e4
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e41
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e42
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e43
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, // e23
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, // e31
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, // e12
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e423
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e431
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e412
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, // e321
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0  // e1234
};

// Accessor function for convenient 2D indexing: G[i,j]
inline constexpr auto pga3dp_metric_view()
{
    return std::mdspan<int const, std::extents<size_t, 16, 16>>{pga3dp_metric.data()};
}

// Regressive extended metric for pga3dp (16×16 matrix, 256 elements)
// Basis: 1, e1, e2, e3, e4, e41, e42, e43, e23, e31, e12, e423, e431, e412, e321, e1234
inline constexpr std::array<int, 256> pga3dp_rmetric = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e2
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e3
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e4
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e41
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e42
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, // e43
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e23
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e31
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e12
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, // e423
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, // e431
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, // e412
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // e321
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1  // e1234
};

// Accessor function for convenient 2D indexing: Ḡ[i,j]
inline constexpr auto pga3dp_rmetric_view()
{
    return std::mdspan<int const, std::extents<size_t, 16, 16>>{pga3dp_rmetric.data()};
}

} // namespace hd::ga::pga


namespace hd::ga::sta {

/////////////////////////////////////////////////////////////////////////////////////////
// 4ds constants
/////////////////////////////////////////////////////////////////////////////////////////

auto const g1_4ds = Vec4ds<value_t>{1.0, 0.0, 0.0, 0.0};
auto const g2_4ds = Vec4ds<value_t>{0.0, 1.0, 0.0, 0.0};
auto const g3_4ds = Vec4ds<value_t>{0.0, 0.0, 1.0, 0.0};
auto const g4_4ds = Vec4ds<value_t>{0.0, 0.0, 0.0, 1.0};

auto const x_dir_4ds = g1_4ds; // direction of x-axis towards infinity (space x)
auto const y_dir_4ds = g2_4ds; // direction of y-axis towards infinity (space y)
auto const z_dir_4ds = g3_4ds; // direction of z-axis towards infinity (space z)
auto const t_dir_4ds = g4_4ds; // direction of t-axis towards infinity (time t on w-comp.)

auto const g1_4ds_mv = MVec4ds<value_t>{g1_4ds}; // g1_4ds as multivector
auto const g2_4ds_mv = MVec4ds<value_t>{g2_4ds}; // g2_4ds as multivector
auto const g3_4ds_mv = MVec4ds<value_t>{g3_4ds}; // g3_4ds as multivector
auto const g4_4ds_mv = MVec4ds<value_t>{g4_4ds}; // g4_4ds as multivector

auto const g14_4ds = BiVec4ds<value_t>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // = wdg(g1,g4)
auto const g24_4ds = BiVec4ds<value_t>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}; // = wdg(g2,g4)
auto const g34_4ds = BiVec4ds<value_t>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0}; // = wdg(g3,g4)
auto const g23_4ds = BiVec4ds<value_t>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}; // = wdg(g2,g3)
auto const g31_4ds = BiVec4ds<value_t>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0}; // = wdg(g3,g1)
auto const g12_4ds = BiVec4ds<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}; // = wdg(g1,g2)

auto const g14_4ds_mv = MVec4ds<value_t>{g14_4ds}; // g14_4ds as multivector
auto const g24_4ds_mv = MVec4ds<value_t>{g24_4ds}; // g24_4ds as multivector
auto const g34_4ds_mv = MVec4ds<value_t>{g34_4ds}; // g34_4ds as multivector
auto const g23_4ds_mv = MVec4ds<value_t>{g23_4ds}; // g23_4ds as multivector
auto const g31_4ds_mv = MVec4ds<value_t>{g31_4ds}; // g31_4ds as multivector
auto const g12_4ds_mv = MVec4ds<value_t>{g12_4ds}; // g12_4ds as multivector

auto const g14_4ds_mv_e = MVec4ds_E<value_t>{g14_4ds}; // g14_4ds as even multivector
auto const g24_4ds_mv_e = MVec4ds_E<value_t>{g24_4ds}; // g24_4ds as even multivector
auto const g34_4ds_mv_e = MVec4ds_E<value_t>{g34_4ds}; // g34_4ds as even multivector
auto const g23_4ds_mv_e = MVec4ds_E<value_t>{g23_4ds}; // g23_4ds as even multivector
auto const g31_4ds_mv_e = MVec4ds_E<value_t>{g31_4ds}; // g31_4ds as even multivector
auto const g12_4ds_mv_e = MVec4ds_E<value_t>{g12_4ds}; // g12_4ds as even multivector

auto const g234_4ds = TriVec4ds<value_t>{1.0, 0.0, 0.0, 0.0};
auto const g314_4ds = TriVec4ds<value_t>{0.0, 1.0, 0.0, 0.0};
auto const g124_4ds = TriVec4ds<value_t>{0.0, 0.0, 1.0, 0.0};
auto const g123_4ds = TriVec4ds<value_t>{0.0, 0.0, 0.0, 1.0};

auto const g234_4ds_mv = MVec4ds<value_t>{g234_4ds}; // g234_4ds as multivector
auto const g314_4ds_mv = MVec4ds<value_t>{g314_4ds}; // g314_4ds as multivector
auto const g124_4ds_mv = MVec4ds<value_t>{g124_4ds}; // g124_4ds as multivector
auto const g123_4ds_mv = MVec4ds<value_t>{g123_4ds}; // g123_4ds as multivector

auto const one_4ds = Scalar4ds<value_t>(1.0);          // Basis element 1
auto const one_4ds_mv = MVec4ds<value_t>{one_4ds};     // 1 as multivector
auto const one_4ds_mv_e = MVec4ds_E<value_t>{one_4ds}; // 1 as even-grade multivector

auto const I_4ds = PScalar4ds<value_t>(1.0);       // I_4ds = g1^g2^g3^g4 = g1234
auto const g1234_4ds = I_4ds;                      // (pseudoscalar)
auto const I_4ds_mv = MVec4ds<value_t>{I_4ds};     // I_4ds as multivector
auto const I_4ds_mv_e = MVec4ds_E<value_t>{I_4ds}; // I_4ds as even grade multivector

// extended metric for sta4ds (16×16 matrix, 256 elements)
// Basis: 1, g1, g2, g3, g4, g14, g24, g34, g23, g31, g12, g234, g314, g124, g123, g1234
inline constexpr std::array<int, 256> sta4ds_metric = {
    //  the diagonal is the metric exomorphism P(e_S) = prod_{i in S} g_i (P-unify):
    //  vectors carry the signature, higher grades the product (NOT the blade square Q)
    1, 0,  0,  0,  0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0,  0, // 1
    0, -1, 0,  0,  0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0,  0, // g1
    0, 0,  -1, 0,  0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0,  0, // g2
    0, 0,  0,  -1, 0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0,  0, // g3
    0, 0,  0,  0,  1, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0,  0, // g4
    0, 0,  0,  0,  0, -1, 0,  0,  0, 0, 0, 0, 0, 0, 0,  0, // g14
    0, 0,  0,  0,  0, 0,  -1, 0,  0, 0, 0, 0, 0, 0, 0,  0, // g24
    0, 0,  0,  0,  0, 0,  0,  -1, 0, 0, 0, 0, 0, 0, 0,  0, // g34
    0, 0,  0,  0,  0, 0,  0,  0,  1, 0, 0, 0, 0, 0, 0,  0, // g23
    0, 0,  0,  0,  0, 0,  0,  0,  0, 1, 0, 0, 0, 0, 0,  0, // g31
    0, 0,  0,  0,  0, 0,  0,  0,  0, 0, 1, 0, 0, 0, 0,  0, // g12
    0, 0,  0,  0,  0, 0,  0,  0,  0, 0, 0, 1, 0, 0, 0,  0, // g234
    0, 0,  0,  0,  0, 0,  0,  0,  0, 0, 0, 0, 1, 0, 0,  0, // g314
    0, 0,  0,  0,  0, 0,  0,  0,  0, 0, 0, 0, 0, 1, 0,  0, // g124
    0, 0,  0,  0,  0, 0,  0,  0,  0, 0, 0, 0, 0, 0, -1, 0, // g123
    0, 0,  0,  0,  0, 0,  0,  0,  0, 0, 0, 0, 0, 0, 0,  -1 // g1234
};

// Accessor function for convenient 2D indexing: G[i,j]
inline constexpr auto sta4ds_metric_view()
{
    return std::mdspan<int const, std::extents<size_t, 16, 16>>{sta4ds_metric.data()};
}

} // namespace hd::ga::sta

namespace hd::ga::cga {

/////////////////////////////////////////////////////////////////////////////////////////
// 2dc constants (conformal algebra of 2d Euclidean space)
//
// basis vectors: e1, e2 Euclidean; e3 projective (w); e4 round (u).
// The null pair satisfies e3*e3 = e4*e4 = 0 and e3.e4 = -1 — a non-orthogonal
// metric (G(3,1,0) after diagonalization). A normalized point embeds as
//
//     q = x*e1 + y*e2 + e3 + 0.5*(x^2 + y^2)*e4        (q is null: q ⟑ q = 0)
//
// Vec2dc stores (x, y, w, u) in the fields named (x, y, z, w) — positional
// field reuse, see detail/type_t/ga_type2dc.hpp.
/////////////////////////////////////////////////////////////////////////////////////////

auto const e1_2dc = Vec2dc<value_t>{1.0, 0.0, 0.0, 0.0};
auto const e2_2dc = Vec2dc<value_t>{0.0, 1.0, 0.0, 0.0};
auto const e3_2dc = Vec2dc<value_t>{0.0, 0.0, 1.0, 0.0}; // projective (w), null vector
auto const e4_2dc = Vec2dc<value_t>{0.0, 0.0, 0.0, 1.0}; // round (u), null vector

auto const e1_2dc_mv = MVec2dc<value_t>{e1_2dc}; // e1_2dc as multivector
auto const e2_2dc_mv = MVec2dc<value_t>{e2_2dc}; // e2_2dc as multivector
auto const e3_2dc_mv = MVec2dc<value_t>{e3_2dc}; // e3_2dc as multivector
auto const e4_2dc_mv = MVec2dc<value_t>{e4_2dc}; // e4_2dc as multivector

// bivector basis: {e31, e32, e12, e14, e24, e34} in fields (vx, vy, vz, mx, my, mz)
auto const e31_2dc = BiVec2dc<value_t>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // = wdg(e3,e1)
auto const e32_2dc = BiVec2dc<value_t>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}; // = wdg(e3,e2)
auto const e12_2dc = BiVec2dc<value_t>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0}; // = wdg(e1,e2)
auto const e14_2dc = BiVec2dc<value_t>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}; // = wdg(e1,e4)
auto const e24_2dc = BiVec2dc<value_t>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0}; // = wdg(e2,e4)
auto const e34_2dc = BiVec2dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}; // = wdg(e3,e4)

auto const e31_2dc_mv = MVec2dc<value_t>{e31_2dc}; // e31_2dc as multivector
auto const e32_2dc_mv = MVec2dc<value_t>{e32_2dc}; // e32_2dc as multivector
auto const e12_2dc_mv = MVec2dc<value_t>{e12_2dc}; // e12_2dc as multivector
auto const e14_2dc_mv = MVec2dc<value_t>{e14_2dc}; // e14_2dc as multivector
auto const e24_2dc_mv = MVec2dc<value_t>{e24_2dc}; // e24_2dc as multivector
auto const e34_2dc_mv = MVec2dc<value_t>{e34_2dc}; // e34_2dc as multivector

auto const e31_2dc_mv_e = MVec2dc_E<value_t>{e31_2dc}; // e31_2dc as even multivector
auto const e32_2dc_mv_e = MVec2dc_E<value_t>{e32_2dc}; // e32_2dc as even multivector
auto const e12_2dc_mv_e = MVec2dc_E<value_t>{e12_2dc}; // e12_2dc as even multivector
auto const e14_2dc_mv_e = MVec2dc_E<value_t>{e14_2dc}; // e14_2dc as even multivector
auto const e24_2dc_mv_e = MVec2dc_E<value_t>{e24_2dc}; // e24_2dc as even multivector
auto const e34_2dc_mv_e = MVec2dc_E<value_t>{e34_2dc}; // e34_2dc as even multivector

// trivector basis: {e314, e324, e124, e321} in fields (x, y, z, w)
auto const e314_2dc = TriVec2dc<value_t>{1.0, 0.0, 0.0, 0.0};
auto const e324_2dc = TriVec2dc<value_t>{0.0, 1.0, 0.0, 0.0};
auto const e124_2dc = TriVec2dc<value_t>{0.0, 0.0, 1.0, 0.0};
auto const e321_2dc = TriVec2dc<value_t>{0.0, 0.0, 0.0, 1.0};

auto const e314_2dc_mv = MVec2dc<value_t>{e314_2dc}; // e314_2dc as multivector
auto const e324_2dc_mv = MVec2dc<value_t>{e324_2dc}; // e324_2dc as multivector
auto const e124_2dc_mv = MVec2dc<value_t>{e124_2dc}; // e124_2dc as multivector
auto const e321_2dc_mv = MVec2dc<value_t>{e321_2dc}; // e321_2dc as multivector

auto const e314_2dc_mv_u = MVec2dc_U<value_t>{e314_2dc}; // e314_2dc as odd multivector
auto const e324_2dc_mv_u = MVec2dc_U<value_t>{e324_2dc}; // e324_2dc as odd multivector
auto const e124_2dc_mv_u = MVec2dc_U<value_t>{e124_2dc}; // e124_2dc as odd multivector
auto const e321_2dc_mv_u = MVec2dc_U<value_t>{e321_2dc}; // e321_2dc as odd multivector

auto const e1234_2dc = PScalar2dc<value_t>(1.0);
auto const I_2dc = e1234_2dc; // the pseudoscalar of cga2dc


// 3dc constants (conformal algebra of 3d Euclidean space)
//
// e1, e2, e3 = Euclidean, e4 = projective/origin (w), e5 = round/infinity (u);
// the null pair satisfies e4*e4 = e5*e5 = 0 and e4.e5 = -1 (see
// detail/type_t/ga_type3dc.hpp for the field layout of each type).

auto const e1_3dc = Vec3dc<value_t>{1.0, 0.0, 0.0, 0.0, 0.0};
auto const e2_3dc = Vec3dc<value_t>{0.0, 1.0, 0.0, 0.0, 0.0};
auto const e3_3dc = Vec3dc<value_t>{0.0, 0.0, 1.0, 0.0, 0.0};
auto const e4_3dc =
    Vec3dc<value_t>{0.0, 0.0, 0.0, 1.0, 0.0}; // projective (w), null vector
auto const e5_3dc = Vec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 1.0}; // round (u), null vector

auto const e1_3dc_mv = MVec3dc<value_t>{e1_3dc}; // e1_3dc as multivector
auto const e2_3dc_mv = MVec3dc<value_t>{e2_3dc}; // e2_3dc as multivector
auto const e3_3dc_mv = MVec3dc<value_t>{e3_3dc}; // e3_3dc as multivector
auto const e4_3dc_mv = MVec3dc<value_t>{e4_3dc}; // e4_3dc as multivector
auto const e5_3dc_mv = MVec3dc<value_t>{e5_3dc}; // e5_3dc as multivector

auto const e41_3dc = BiVec3dc<value_t>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e42_3dc = BiVec3dc<value_t>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e43_3dc = BiVec3dc<value_t>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e23_3dc = BiVec3dc<value_t>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e31_3dc = BiVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e12_3dc = BiVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
auto const e15_3dc = BiVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
auto const e25_3dc = BiVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
auto const e35_3dc = BiVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
auto const e45_3dc = BiVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

auto const e41_3dc_mv = MVec3dc<value_t>{e41_3dc}; // e41_3dc as multivector
auto const e42_3dc_mv = MVec3dc<value_t>{e42_3dc}; // e42_3dc as multivector
auto const e43_3dc_mv = MVec3dc<value_t>{e43_3dc}; // e43_3dc as multivector
auto const e23_3dc_mv = MVec3dc<value_t>{e23_3dc}; // e23_3dc as multivector
auto const e31_3dc_mv = MVec3dc<value_t>{e31_3dc}; // e31_3dc as multivector
auto const e12_3dc_mv = MVec3dc<value_t>{e12_3dc}; // e12_3dc as multivector
auto const e15_3dc_mv = MVec3dc<value_t>{e15_3dc}; // e15_3dc as multivector
auto const e25_3dc_mv = MVec3dc<value_t>{e25_3dc}; // e25_3dc as multivector
auto const e35_3dc_mv = MVec3dc<value_t>{e35_3dc}; // e35_3dc as multivector
auto const e45_3dc_mv = MVec3dc<value_t>{e45_3dc}; // e45_3dc as multivector

auto const e41_3dc_mv_e = MVec3dc_E<value_t>{e41_3dc}; // e41_3dc as even multivector
auto const e42_3dc_mv_e = MVec3dc_E<value_t>{e42_3dc}; // e42_3dc as even multivector
auto const e43_3dc_mv_e = MVec3dc_E<value_t>{e43_3dc}; // e43_3dc as even multivector
auto const e23_3dc_mv_e = MVec3dc_E<value_t>{e23_3dc}; // e23_3dc as even multivector
auto const e31_3dc_mv_e = MVec3dc_E<value_t>{e31_3dc}; // e31_3dc as even multivector
auto const e12_3dc_mv_e = MVec3dc_E<value_t>{e12_3dc}; // e12_3dc as even multivector
auto const e15_3dc_mv_e = MVec3dc_E<value_t>{e15_3dc}; // e15_3dc as even multivector
auto const e25_3dc_mv_e = MVec3dc_E<value_t>{e25_3dc}; // e25_3dc as even multivector
auto const e35_3dc_mv_e = MVec3dc_E<value_t>{e35_3dc}; // e35_3dc as even multivector
auto const e45_3dc_mv_e = MVec3dc_E<value_t>{e45_3dc}; // e45_3dc as even multivector

auto const e415_3dc =
    TriVec3dc<value_t>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e425_3dc =
    TriVec3dc<value_t>{0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e435_3dc =
    TriVec3dc<value_t>{0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e235_3dc =
    TriVec3dc<value_t>{0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e315_3dc =
    TriVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
auto const e125_3dc =
    TriVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
auto const e423_3dc =
    TriVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
auto const e431_3dc =
    TriVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
auto const e412_3dc =
    TriVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
auto const e321_3dc =
    TriVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

auto const e415_3dc_mv = MVec3dc<value_t>{e415_3dc}; // e415_3dc as multivector
auto const e425_3dc_mv = MVec3dc<value_t>{e425_3dc}; // e425_3dc as multivector
auto const e435_3dc_mv = MVec3dc<value_t>{e435_3dc}; // e435_3dc as multivector
auto const e235_3dc_mv = MVec3dc<value_t>{e235_3dc}; // e235_3dc as multivector
auto const e315_3dc_mv = MVec3dc<value_t>{e315_3dc}; // e315_3dc as multivector
auto const e125_3dc_mv = MVec3dc<value_t>{e125_3dc}; // e125_3dc as multivector
auto const e423_3dc_mv = MVec3dc<value_t>{e423_3dc}; // e423_3dc as multivector
auto const e431_3dc_mv = MVec3dc<value_t>{e431_3dc}; // e431_3dc as multivector
auto const e412_3dc_mv = MVec3dc<value_t>{e412_3dc}; // e412_3dc as multivector
auto const e321_3dc_mv = MVec3dc<value_t>{e321_3dc}; // e321_3dc as multivector

auto const e415_3dc_mv_u = MVec3dc_U<value_t>{e415_3dc}; // e415_3dc as odd multivector
auto const e425_3dc_mv_u = MVec3dc_U<value_t>{e425_3dc}; // e425_3dc as odd multivector
auto const e435_3dc_mv_u = MVec3dc_U<value_t>{e435_3dc}; // e435_3dc as odd multivector
auto const e235_3dc_mv_u = MVec3dc_U<value_t>{e235_3dc}; // e235_3dc as odd multivector
auto const e315_3dc_mv_u = MVec3dc_U<value_t>{e315_3dc}; // e315_3dc as odd multivector
auto const e125_3dc_mv_u = MVec3dc_U<value_t>{e125_3dc}; // e125_3dc as odd multivector
auto const e423_3dc_mv_u = MVec3dc_U<value_t>{e423_3dc}; // e423_3dc as odd multivector
auto const e431_3dc_mv_u = MVec3dc_U<value_t>{e431_3dc}; // e431_3dc as odd multivector
auto const e412_3dc_mv_u = MVec3dc_U<value_t>{e412_3dc}; // e412_3dc as odd multivector
auto const e321_3dc_mv_u = MVec3dc_U<value_t>{e321_3dc}; // e321_3dc as odd multivector

auto const e4235_3dc = QuadVec3dc<value_t>{1.0, 0.0, 0.0, 0.0, 0.0};
auto const e4315_3dc = QuadVec3dc<value_t>{0.0, 1.0, 0.0, 0.0, 0.0};
auto const e4125_3dc = QuadVec3dc<value_t>{0.0, 0.0, 1.0, 0.0, 0.0};
auto const e3215_3dc = QuadVec3dc<value_t>{0.0, 0.0, 0.0, 1.0, 0.0};
auto const e1234_3dc = QuadVec3dc<value_t>{0.0, 0.0, 0.0, 0.0, 1.0};

auto const e4235_3dc_mv = MVec3dc<value_t>{e4235_3dc}; // e4235_3dc as multivector
auto const e4315_3dc_mv = MVec3dc<value_t>{e4315_3dc}; // e4315_3dc as multivector
auto const e4125_3dc_mv = MVec3dc<value_t>{e4125_3dc}; // e4125_3dc as multivector
auto const e3215_3dc_mv = MVec3dc<value_t>{e3215_3dc}; // e3215_3dc as multivector
auto const e1234_3dc_mv = MVec3dc<value_t>{e1234_3dc}; // e1234_3dc as multivector

auto const e4235_3dc_mv_e =
    MVec3dc_E<value_t>{e4235_3dc}; // e4235_3dc as even multivector
auto const e4315_3dc_mv_e =
    MVec3dc_E<value_t>{e4315_3dc}; // e4315_3dc as even multivector
auto const e4125_3dc_mv_e =
    MVec3dc_E<value_t>{e4125_3dc}; // e4125_3dc as even multivector
auto const e3215_3dc_mv_e =
    MVec3dc_E<value_t>{e3215_3dc}; // e3215_3dc as even multivector
auto const e1234_3dc_mv_e =
    MVec3dc_E<value_t>{e1234_3dc}; // e1234_3dc as even multivector

auto const e12345_3dc = PScalar3dc<value_t>(1.0);
auto const I_3dc = e12345_3dc; // the pseudoscalar of cga3dc

} // namespace hd::ga::cga