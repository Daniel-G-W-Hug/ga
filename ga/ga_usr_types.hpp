#pragma once

// author: Daniel Hug, 2024

#include "ga_cfg_bivec3d.hpp" // BiVec3d<T>
#include "ga_cfg_mvec2d.hpp"  // Vec2d<T>
#include "ga_cfg_mvec3d.hpp"  // Vec3d<T>
#include "ga_cfg_strong_t.hpp"
#include "ga_cfg_value_t.hpp" // defines value_t either as float oder double
#include "ga_cfg_vec2d.hpp"   // Vec2d<T>
#include "ga_cfg_vec3d.hpp"   // Vec3d<T>

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions for easy use
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

using scalar = Scalar<value_t>;

// The pseudoscalar is a one-dimensional entity, i.e. it has scalar character.
// The pseudoscalar describes the "unit volume" of the respective space.
// It is a bivector for 2d space, a trivector for 3d space, etc.
// It's scalar component describes a multiple of the unit volume of that space.
//
// Thus the pseudoscalar needs a unique type for each space in order to not mix up
// pseudoscalars from different algebras in the type system.
//
using pscalar2d = PScalar2d<value_t>;
using pscalar3d = PScalar3d<value_t>;

// vector and multivector types
using vec2d = Vec2d<value_t>;
using vec3d = Vec3d<value_t>;
using bivec3d = BiVec3d<value_t>;
using mvec2d = MVec2d<value_t>;     // fully populated 2d multivector
using mvec2d_e = MVec2d_E<value_t>; // multivector 2d of even subalgebra
using mvec3d = MVec3d<value_t>;     // fully populated 3d multivector
using mvec3d_e = MVec3d_E<value_t>; // multivector 3d of even subalgebra
using mvec3d_u = MVec3d_U<value_t>; // multivector 3d of uneven subalgebra

} // namespace hd::ga