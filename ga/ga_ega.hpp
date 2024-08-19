#pragma once

// author: Daniel Hug, 2024

// this is a convenience header to include all required headers consistently

////////////////////////////////////////////////////////////////////////////////
// provide data types for representing EGA in 2D and 3D
// e.g. scalar, vector, bivector, etc.
// and corresponding operations
////////////////////////////////////////////////////////////////////////////////

// commonly used GA stuff is in namespace hd::ga
#include "detail/ga_value_t.hpp" // default type for scalar, vector & bivector components
                                 // switch value_t between float/double here, if needed!

#include "detail/ga_algebra.hpp" // algebra (define basis components)

#include "ga_usr_consts.hpp"    // user constants (for user convenience)
#include "ga_usr_types.hpp"     // consistent user types (for user convenience)
#include "ga_usr_utilities.hpp" // chelper functions (for user convenience)

#include "detail/ga_strong_t.hpp" // common type for scalars & pseudoscalars
                                  // Scalar<T>, PScalar2d<T>, PScalar3d<T>

#include "detail/ga_vec2d.hpp" // Vec2d<T>
#include "detail/ga_vec3d.hpp" // Vec3d<T>

#include "detail/ga_bivec3d.hpp"  // BiVec3d<T>
#include "detail/ga_mvec2d.hpp"   // MVec2d<T>
#include "detail/ga_mvec2d_e.hpp" // MVec2d_E<T>
#include "detail/ga_mvec3d.hpp"   // MVec3d<T>
#include "detail/ga_mvec3d_e.hpp" // MVec3d_E<T>
#include "detail/ga_mvec3d_u.hpp" // MVec3d_U<T>

// EGA-specific operations are in namespace hd::ga::ega
#include "detail/ga_ega_2d_ops.hpp" // ega 2d operations
#include "detail/ga_ega_3d_ops.hpp" // ega 3d operations

// fmt-support is defined outside of other namespaces
#include "detail/ga_fmt_support.hpp" // printing support (fmt library)