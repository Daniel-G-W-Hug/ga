#pragma once

#include "detail/ga_value_t.hpp" // defines value_t either as float oder double

#include "detail/ga_scalar_t.hpp"

#include "detail/ga_mvec2d.hpp" // MVec2d<T>
#include "detail/ga_vec2d.hpp"  // Vec2d<T>

#include "detail/ga_bivec3d.hpp" // BiVec3d<T>
#include "detail/ga_mvec3d.hpp"  // MVec3d<T>
#include "detail/ga_vec3d.hpp"   // Vec3d<T>

namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// 2D constants
////////////////////////////////////////////////////////////////////////////////

auto const e1_2d = Vec2d<value_t>{1.0, 0.0};
auto const e2_2d = Vec2d<value_t>{0.0, 1.0};
auto const e1m_2d = MVec2d<value_t>{e1_2d}; // e1_2d as multivector
auto const e2m_2d = MVec2d<value_t>{e2_2d}; // e2_2d as multivector

auto const I_2d = PScalar2d<value_t>(1.0);
auto const Im_2d = MVec2d<value_t>{I_2d};     // I_2d as multivector
auto const Im_2d_E = MVec2d_E<value_t>{I_2d}; // I_2d as even grade multivector

////////////////////////////////////////////////////////////////////////////////
// 3D constants
////////////////////////////////////////////////////////////////////////////////

auto const e1_3d = Vec3d<value_t>{1.0, 0.0, 0.0};
auto const e2_3d = Vec3d<value_t>{0.0, 1.0, 0.0};
auto const e3_3d = Vec3d<value_t>{0.0, 0.0, 1.0};
auto const e1m_3d = MVec3d<value_t>{e1_3d}; // e1_3d as multivector
auto const e2m_3d = MVec3d<value_t>{e2_3d}; // e2_3d as multivector
auto const e3m_3d = MVec3d<value_t>{e3_3d}; // e3_3d as multivector

auto const e23_3d = BiVec3d<value_t>{1.0, 0.0, 0.0};
auto const e31_3d = BiVec3d<value_t>{0.0, 1.0, 0.0};
auto const e12_3d = BiVec3d<value_t>{0.0, 0.0, 1.0};
auto const e23m_3d = MVec3d<value_t>{e23_3d};    // e23_3d as multivector
auto const e31m_3d = MVec3d<value_t>{e31_3d};    // e31_3d as multivector
auto const e12m_3d = MVec3d<value_t>{e12_3d};    // e12_3d as multivector
auto const e23me_3d = MVec3d_E<value_t>{e23_3d}; // e23_3d as even grade multivector
auto const e31me_3d = MVec3d_E<value_t>{e31_3d}; // e31_3d as even grade multivector
auto const e12me_3d = MVec3d_E<value_t>{e12_3d}; // e12_3d as even grade multivector

auto const I_3d = PScalar3d<value_t>(1.0);
auto const Im_3d = MVec3d<value_t>{I_3d};     // I_3d as multivector
auto const Im_3d_U = MVec3d_U<value_t>{I_3d}; // I_3d as uneven grade multivector

} // namespace hd::ga