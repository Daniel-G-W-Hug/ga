#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

// this is a convenience header to include all required headers consistently

////////////////////////////////////////////////////////////////////////////////
// provide data types for representing EGA in 2D, 3D and 4D
// e.g. scalar, vector, bivector, etc.
// and corresponding operations
////////////////////////////////////////////////////////////////////////////////

// types used for EGA (provided as templated types):
//
// Scalar2d<T>, Vec2dp<T>, PScalar2dp<T>, MVec2d_E<T>, MVec2dp<T>
// Scalar3d<T>, Vec3dp<T>, BiVec3dp<T>, PScalar3dp<T>, MVec3d_E<T>, MVec3d_U<T>,
// MVec3dp<T>
// Scalar4d<T>, Vec4dp<T>, BiVec4dp<T>, TriVec4dp<T>, PScalar4dp<T>, MVec4d_E<T>,
// MVec4d_U<T>, MVec4dp<T>

// types used for EGA (provided as user types based on value_t for convenience):
//
// scalar2d, vec2d, pscalar2d, mvec2d_e, mvec2d
// scalar3d, vec3d, bivec3d, pscalar3d, mvec3d_e, mvec3d_u, mvec3d
// scalar4d, vec4d, bivec4d, trivec4d, pscalar4d, mvec4d_e, mvec4d_u, mvec4d

// commonly used GA stuff is in namespace hd::ga
#include "ga_value_t.hpp" // default type for scalar, vector & bivector components
                          // switch value_t between float/double here, if needed!

#include "ga_algebra.hpp" // algebra (define basis components)

#include "ga_usr_consts.hpp"    // user constants (for user convenience)
#include "ga_usr_types.hpp"     // consistent user types (for user convenience)
#include "ga_usr_utilities.hpp" // chelper functions (for user convenience)

// EGA-specific operations are in namespace hd::ga::ega
#include "ga_ega_2d_ops.hpp" // ega 2d operations
#include "ga_ega_3d_ops.hpp" // ega 3d operations
// TODO: implement operations for ega_4d (better do sta_4d)

// fmt-support is defined outside of other namespaces
#include "detail/ga_fmt_support.hpp" // printing support (fmt library)