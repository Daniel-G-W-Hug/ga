#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// this is a convenience header to include all required headers consistently

/////////////////////////////////////////////////////////////////////////////////////////
// provide data types for representing CGA in 2D (and later 3D)
// e.g. scalar, vector, bivector, etc.
// and corresponding operations
/////////////////////////////////////////////////////////////////////////////////////////

// types used for CGA (provided as templated types):
//
// Scalar2dc<T>, Vec2dc<T>, BiVec2dc<T>, TriVec2dc<T>, PScalar2dc<T>,
// MVec2dc_E<T>, MVec2dc_U<T>, MVec2dc<T>

// types used for CGA (provided as user types based on value_t for convenience):
//
// scalar2dc, vec2dc, bivec2dc, trivec2dc, pscalar2dc, mvec2dc_e, mvec2dc_u, mvec2dc

// commonly used GA stuff is in namespace hd::ga
#include "ga_value_t.hpp" // default type for scalar, vector & bivector components
                          // switch value_t between float/double here, if needed!

// for user convenience
#include "ga_usr_consts.hpp"    // user constants
#include "ga_usr_types.hpp"     // consistent user types
#include "ga_usr_utilities.hpp" // helper functions

// CGA models Euclidean space, so its Euclidean accessors (position, radius,
// direction, normal) hand back EGA vectors -- pull in the EGA operations too,
// or the returned vectors would be inert here (their types alone come with the
// foundation, their operations do not)
#include "ga_ega.hpp"

// CGA-specific operations are in namespace hd::ga::cga
#include "ga_cga2dc_ops_basics.hpp"   // involutions, complements, duals, antiduals
#include "ga_cga2dc_ops_products.hpp" // products (dot, wdg, contractions, gpr, rgpr,
                                      // inv, rinv, ...)
#include "ga_cga2dc_ops.hpp"          // geometric operations (is_congruent, is_close;
                                      // layer under construction)

#include "ga_cga3dc_ops_basics.hpp"   // involutions, complement, dual, antidual
#include "ga_cga3dc_ops_products.hpp" // products (dot, wdg, contractions, gpr, rgpr,
                                      // inv, rinv, ...)
#include "ga_cga3dc_ops.hpp"          // geometric operations (is_congruent, is_close;
                                      // layer under construction)

// fmt-support is defined outside of other namespaces
#include "detail/ga_fmt_support.hpp" // printing support (fmt library)
