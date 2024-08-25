#pragma once

// author: Daniel Hug, 2024

// this is a convenience header to include all required headers consistently

////////////////////////////////////////////////////////////////////////////////
// provide data types for representing EGA in 2D, 3D and 4D
// e.g. scalar, vector, bivector, etc.
// and corresponding operations
////////////////////////////////////////////////////////////////////////////////

// types used for EGA (provided as templated types):
//
// Scalar<T>
// Vec2d<T>, PScalar2d<T>, MVec2d_E<T>, MVec2d_U<T>, MVec2d<T>
// Vec3d<T>, BiVec3d<T>, PScalar3d<T>, MVec3d_E<T>, MVec3d_U<T>, MVec3d<T>
// Vec4d<T>, BiVec4d<T>, TriVec4d<T>, PScalar4d<T>, MVec4d_E<T>, MVec4d_U<T>, MVec4d<T>

// types used for EGA (provided as user types based on value_t for convenience):
//
// scalar
// vec2d, pscalar2d, mvec2d_e, mvec2d_u, mvec2d
// vec3d, bivec3d, pscalar3d, mvec3d_e, mvec3d_u, mvec3d
// vec4d, bivec4d, trivec4d, pscalar4d, mvec4d_e, mvec4d_u, mvec4d

// commonly used GA stuff is in namespace hd::ga
#include "detail/ga_value_t.hpp" // default type for scalar, vector & bivector components
                                 // switch value_t between float/double here, if needed!

#include "detail/ga_algebra.hpp" // algebra (define basis components)

#include "ga_usr_consts.hpp"    // user constants (for user convenience)
#include "ga_usr_types.hpp"     // consistent user types (for user convenience)
#include "ga_usr_utilities.hpp" // chelper functions (for user convenience)

#include "detail/ga_mvec2d.hpp" // inclusion of multivector imports all component types
#include "detail/ga_mvec3d.hpp" // inclusion of multivector imports all component types

// EGA-specific operations are in namespace hd::ga::ega
#include "detail/ga_ega_2d_ops.hpp" // ega 2d operations
#include "detail/ga_ega_3d_ops.hpp" // ega 3d operations
// TODO: implement operations for ega_4d

// fmt-support is defined outside of other namespaces
#include "detail/ga_fmt_support.hpp" // printing support (fmt library)