#pragma once

#include "ga_value_t.hpp" // defines value_t either as float oder double
                          // defines Scalar<T>, PScalar2d<T>, PScalar3d<T>
                          // Vec2d<T>, Vec3d<T>, BiVec3d<T>

#include "detail/ga_mvec2d.hpp" // MVec2d<T> and subtypes
#include "detail/ga_mvec3d.hpp" // MVec3d<T> and subtypes

#include "detail/ga_mvec2dp.hpp" // MVec2dp<T> and subtypes
// #include "detail/ga_mvec3dp.hpp" // MVec3dp<T> and subtypes

// provide constants for ega and pga in their respective namespaces

namespace hd::ga::ega {

////////////////////////////////////////////////////////////////////////////////
// 2d constants
////////////////////////////////////////////////////////////////////////////////

auto const e1_2d = Vec2d<value_t>{1.0, 0.0};
auto const e2_2d = Vec2d<value_t>{0.0, 1.0};
auto const e1m_2d = MVec2d<value_t>{e1_2d}; // e1_2d as multivector
auto const e2m_2d = MVec2d<value_t>{e2_2d}; // e2_2d as multivector

auto const I_2d = PScalar2d<value_t>(1.0);
auto const Im_2d = MVec2d<value_t>{I_2d};     // I_2d as multivector
auto const Im_2d_E = MVec2d_E<value_t>{I_2d}; // I_2d as even grade multivector

////////////////////////////////////////////////////////////////////////////////
// 3d constants
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

} // namespace hd::ga::ega


namespace hd::ga::pga {

////////////////////////////////////////////////////////////////////////////////
// 2dp constants
////////////////////////////////////////////////////////////////////////////////

auto const e1_2dp = Vec2dp<value_t>{1.0, 0.0, 0.0};
auto const e2_2dp = Vec2dp<value_t>{0.0, 1.0, 0.0};
auto const e3_2dp = Vec2dp<value_t>{0.0, 0.0, 1.0};
auto const e1m_2dp = MVec2dp<value_t>{e1_2dp}; // e1_2dp as multivector
auto const e2m_2dp = MVec2dp<value_t>{e2_2dp}; // e2_2dp as multivector
auto const e3m_2dp = MVec2dp<value_t>{e3_2dp}; // e3_2dp as multivector

auto const e23_2dp = BiVec2dp<value_t>{1.0, 0.0, 0.0};
auto const e31_2dp = BiVec2dp<value_t>{0.0, 1.0, 0.0};
auto const e12_2dp = BiVec2dp<value_t>{0.0, 0.0, 1.0};
auto const e23m_2dp = MVec2dp<value_t>{e23_2dp};    // e23_2dp as multivector
auto const e31m_2dp = MVec2dp<value_t>{e31_2dp};    // e31_2dp as multivector
auto const e12m_2dp = MVec2dp<value_t>{e12_2dp};    // e12_2dp as multivector
auto const e23me_2dp = MVec2dp_E<value_t>{e23_2dp}; // e23_2dp as even grade multivector
auto const e31me_2dp = MVec2dp_E<value_t>{e31_2dp}; // e31_2dp as even grade multivector
auto const e12me_2dp = MVec2dp_E<value_t>{e12_2dp}; // e12_2dp as even grade multivector

auto const I_2dp = PScalar2dp<value_t>(1.0);
auto const Im_2dp = MVec2dp<value_t>{I_2dp};     // I_2dp as multivector
auto const Im_2dp_U = MVec2dp_U<value_t>{I_2dp}; // I_2dp as uneven grade multivector

} // namespace hd::ga::pga