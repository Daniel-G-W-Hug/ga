#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

// this is a convenience header to include all required headers consistently

/////////////////////////////////////////////////////////////////////////////////////////
// provide data types for representing PGA in 2D and 3D
// e.g. scalar, vector, bivector, etc.
// and corresponding operations
/////////////////////////////////////////////////////////////////////////////////////////

// types used for PGA (provided as templated types):
//
// Scalar2dp<T>, Vec2dp<T>, BiVec2dp<T>, PScalar2dp<T>, MVec2dp_E<T>, MVec2dp_U<T>,
// MVec2dp<T> Scalar3dp<T>, Vec3dp<T>, BiVec3dp<T>, TriVec3dp<T>, PScalar3dp<T>,
// MVec3dp_E<T>, MVec3dp_U<T>, MVec3dp<T> DualNum2dp<T>, DualNum3dp<T>

// types used for PGA (provided as user types based on value_t for convenience):
//
// scalar2dp, vec2dp, bivec2dp, pscalar2dp, mvec2dp_e, mvec2dp
// scalar3dp, vec3dp, bivec3dp, trivec3dp, pscalar3dp, mvec3dp_e, mvec3dp_u, mvec3dp

// commonly used GA stuff is in namespace hd::ga
#include "ga_value_t.hpp" // default type for scalar, vector & bivector components
                          // switch value_t between float/double here, if needed!

// for user convenience
#include "ga_usr_consts.hpp"    // user constants
#include "ga_usr_types.hpp"     // consistent user types
#include "ga_usr_utilities.hpp" // helper functions

// PGA-specific operations are in namespace hd::ga::pga
#include "ga_pga2dp_ops.hpp" // include all pga operations for 2dp
#include "ga_pga3dp_ops.hpp" // include all pga operations for 3dp

// fmt-support is defined outside of other namespaces
#include "detail/ga_fmt_support.hpp" // printing support (fmt library)
