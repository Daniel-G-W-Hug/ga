#pragma once

// author: Daniel Hug, 2024 & 2025
#include "ga_value_t.hpp" // defines value_t either as float oder double

#include "detail/ga_mvec2d.hpp" // 2d types for EGA
#include "detail/ga_mvec3d.hpp" // 3d types for EGA
#include "detail/ga_mvec4d.hpp" // 4d types for EGA

#include "detail/ga_mvec2dp.hpp" // 2dp types for PGA
#include "detail/ga_mvec3dp.hpp" // 3dp types for PGA

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions for easy use
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

// The pseudoscalar is a one-dimensional entity, i.e. it has scalar character.
// The pseudoscalar describes the "unit volume" of the respective space.
// It is a bivector for 2d space, a trivector for 3d space, etc.
// It's scalar component describes a multiple of the unit volume of that space.
//
// Thus the pseudoscalar needs a unique type for each space in order to not mix up
// pseudoscalars from different algebras in the type system.
//

// There is actually only one scalar type!
// But due to the statically typed language, overload resolution and unique scalar types
// for 2d, 3d, 2dp, etc. work better for avoiding ambiguous function signatures.

// vector and multivector types EGA 2d

using scalar2d = Scalar2d<value_t>;
using vec2d = Vec2d<value_t>;
using pscalar2d = PScalar2d<value_t>;
using mvec2d_e = MVec2d_E<value_t>; // multivector 2d of even subalgebra
using mvec2d = MVec2d<value_t>;     // fully populated 2d multivector

// vector and multivector types EGA 3d
using scalar3d = Scalar3d<value_t>;
using vec3d = Vec3d<value_t>;
using bivec3d = BiVec3d<value_t>;
using pscalar3d = PScalar3d<value_t>;
using mvec3d_e = MVec3d_E<value_t>; // multivector 3d of even subalgebra
using mvec3d_u = MVec3d_U<value_t>; // multivector 3d of uneven subalgebra
using mvec3d = MVec3d<value_t>;     // fully populated 3d multivector

// vector and multivector types EGA 4d
using scalar4d = Scalar4d<value_t>;
using vec4d = Vec4d<value_t>;
using bivec4d = BiVec4d<value_t>;
using trivec4d = TriVec4d<value_t>;
using pscalar4d = PScalar4d<value_t>;
using mvec4d_e = MVec4d_E<value_t>; // multivector 4d of even subalgebra
using mvec4d_u = MVec4d_U<value_t>; // multivector 4d of uneven subalgebra
using mvec4d = MVec4d<value_t>;     // fully populated 4d multivector

// vector and multivector types PGA 2dp
using scalar2dp = Scalar2dp<value_t>;
using vec2dp = Vec2dp<value_t>;
using bivec2dp = BiVec2dp<value_t>;
using pscalar2dp = PScalar2dp<value_t>;
using mvec2dp_e = MVec2dp_E<value_t>; // multivector 2dp of even subalgebra
using mvec2dp_u = MVec2dp_U<value_t>; // multivector 2dp of uneven subalgebra
using mvec2dp = MVec2dp<value_t>;     // fully populated 2dp multivector

// vector and multivector types PGA 3dp
using scalar3dp = Scalar3dp<value_t>;
using vec3dp = Vec3dp<value_t>;
using bivec3dp = BiVec3dp<value_t>;
using trivec3dp = TriVec3dp<value_t>;
using pscalar3dp = PScalar3dp<value_t>;
using mvec3dp_e = MVec3dp_E<value_t>; // multivector 3dp of even subalgebra
using mvec3dp_u = MVec3dp_U<value_t>; // multivector 3dp of uneven subalgebra
using mvec3dp = MVec3dp<value_t>;     // fully populated 3dp multivector

} // namespace hd::ga