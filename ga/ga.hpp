#pragma once

// author: Daniel Hug, 2024

// this is a convenience header to include all required headers consistently

////////////////////////////////////////////////////////////////////////////////
// provide data types for representing GA in 2D and 3D
// e.g. scalar, vector, bivector, etc.
// and corresponding operations
////////////////////////////////////////////////////////////////////////////////

#include "ga_cfg_algebra.hpp" // algebra

#include "ga_cfg_strong_t.hpp" // strong types for scalars & pseudoscalars
#include "ga_cfg_value_t.hpp"  // default type for scalar, vector & bivector components

#include "ga_cfg_pscalar2d.hpp" // PScalar2d<T>
#include "ga_cfg_pscalar3d.hpp" // PScalar3d<T>
#include "ga_cfg_scalar.hpp"    // Scalar<T>
#include "ga_cfg_vec2d.hpp"     // Vec2d<T>
#include "ga_cfg_vec3d.hpp"     // Vec3d<T>

#include "ga_cfg_bivec3d.hpp"  // BiVec3d<T>
#include "ga_cfg_mvec2d.hpp"   // MVec2d<T>
#include "ga_cfg_mvec2d_e.hpp" // MVec2d_E<T>
#include "ga_cfg_mvec3d.hpp"   // MVec3d<T>
#include "ga_cfg_mvec3d_e.hpp" // MVec3d_E<T>
#include "ga_cfg_mvec3d_u.hpp" // MVec3d_U<T>

#include "ga_cfg_2d_ops.hpp" // 2d operations
#include "ga_cfg_3d_ops.hpp" // 3d operations

#include "ga_usr_consts.hpp"    // user constants
#include "ga_usr_types.hpp"     // consistent user types
#include "ga_usr_utilities.hpp" // convenient helper functions
