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
// Vec2dp<T>, PScalar2dp<T>, MVec2d_E<T>, MVec2d_U<T>, MVec2dp<T>
// Vec3dp<T>, BiVec3dp<T>, PScalar3dp<T>, MVec3d_E<T>, MVec3d_U<T>, MVec3dp<T>
// Vec4dp<T>, BiVec4dp<T>, TriVec4dp<T>, PScalar4dp<T>, MVec4d_E<T>, MVec4d_U<T>,
// MVec4dp<T>

// types used for PGA (provided as user types based on value_t for convenience):
//
// scalar
// vec2dp, bivec2dp, pscalar2dp, mvec2dp_e, mvec2dp_u, mvec2dp
// vec3dp, bivec3dp, trivec3dp, pscalar3dp, mvec3dp_e, mvec3dp_u, mvec3dp

// commonly used GA stuff is in namespace hd::ga
#include "detail/ga_value_t.hpp" // default type for scalar, vector & bivector components
                                 // switch value_t between float/double here, if needed!

#include "detail/ga_algebra.hpp" // algebra (define basis components)

#include "ga_usr_consts.hpp"    // user constants (for user convenience)
#include "ga_usr_types.hpp"     // consistent user types (for user convenience)
#include "ga_usr_utilities.hpp" // chelper functions (for user convenience)

#include "detail/ga_mvec2dp.hpp" // inclusion of multivector imports all component types
// #include "detail/ga_mvec3dp.hpp" // inclusion of multivector imports all component
// types

// PGA-specific operations are in namespace hd::ga::ega
#include "detail/ga_pga_2dp_ops.hpp" // pga 2d operations
// #include "detail/ga_pga_3d_ops.hpp" // pga 3d operations

// fmt-support is defined outside of other namespaces
#include "detail/ga_fmt_support.hpp" // printing support (fmt library)