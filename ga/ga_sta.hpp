#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// this is a convenience header to include all required headers consistently

/////////////////////////////////////////////////////////////////////////////////////////
// provide data types for representing PGA in 2D and 3D
// e.g. scalar, vector, bivector, etc.
// and corresponding operations
/////////////////////////////////////////////////////////////////////////////////////////

// types used for STA (provided as templated types):
//
// Scalar4ds<T>, Vec4ds<T>, BiVec4ds<T>, TriVec4ds<T>, PScalar4ds<T>,
// MVec4ds_E<T>, MVec4ds_U<T>, MVec4ds<T> DualNum2dp<T>, DualNum4ds<T>

// types used for STA (provided as user types based on value_t for convenience):
//
// scalar4ds, vec4ds, bivec4ds, trivec4ds, pscalar4ds, mvec4ds_e, mvec4ds_u, mvec4ds

// commonly used GA stuff is in namespace hd::ga
#include "ga_value_t.hpp" // default type for scalar, vector, bi-, tri-, quadvector comp.
                          // switch value_t between float/double here, if needed!

// for user convenience
#include "ga_usr_consts.hpp"    // user constants
#include "ga_usr_types.hpp"     // consistent user types
#include "ga_usr_utilities.hpp" // helper functions

// STA-specific operations are in namespace hd::ga::pga
#include "ga_sta4ds_ops.hpp" // include all STA operations for 4ds

// fmt-support is defined outside of other namespaces
#include "detail/ga_fmt_support.hpp" // printing support (fmt library)
