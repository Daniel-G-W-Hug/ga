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

// types used for PGA (provided as user types based on value_t for convenience):
//
// scalar2dp, vec2dp, bivec2dp, pscalar2dp, mvec2dp_e, mvec2dp
// scalar3dp, vec3dp, bivec3dp, trivec3dp, pscalar3dp, mvec3dp_e, mvec3dp_u, mvec3dp

// commonly used GA stuff is in namespace hd::ga
#include "ga_value_t.hpp" // default type for scalar, vector & bivector components
                          // switch value_t between float/double here, if needed!

#include "ga_algebra.hpp" // algebra (define basis components)

#include "ga_usr_consts.hpp"    // user constants (for user convenience)
#include "ga_usr_types.hpp"     // consistent user types (for user convenience)
#include "ga_usr_utilities.hpp" // chelper functions (for user convenience)

// PGA-specific operations are in namespace hd::ga::pga
#include "ga_pga_common_ops.hpp" // include all pga operations for 2dp and 3dp

// fmt-support is defined outside of other namespaces
#include "detail/ga_fmt_support.hpp" // printing support (fmt library)
