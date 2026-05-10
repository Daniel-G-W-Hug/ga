#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// Selective include for 3D Projective Geometric Algebra only
// Provides complete pga4ds functionality with minimal compile overhead
/////////////////////////////////////////////////////////////////////////////////////////

#include "detail/ga_foundation.hpp"

// pga4ds specific type definitions
#include "detail/type_t/ga_type4ds.hpp"

// pga4ds multivector implementations
#include "detail/type_t/ga_mvec4ds.hpp" // Complete 4ds multivectors incl. even/add grade

// user convenience includes
#include "ga_usr_consts.hpp"
#include "ga_usr_types.hpp" // Contains 4ds type aliases like scalar4ds, vec4ds, etc.
#include "ga_usr_utilities.hpp"

// pga4ds operations (includes basics and products)
#include "ga_sta4ds_ops.hpp"

// Formatting support
#include "detail/fmt/ga_fmt_core.hpp"
#include "detail/fmt/ga_fmt_pga.hpp"

////////////////////////////////////////////////////////////////////////////////
// This header provides complete 4ds PGA functionality including:
//
// **Types available:**
// - Scalar4ds<T>, Vec4ds<T>, BiVec4ds<T>, TriVec4ds<T>, PScalar4ds<T>
// - MVec4ds_E<T> (even subalgebra: scalar + bivector + pseudoscalar)
// - MVec4ds_U<T> (odd subalgebra: vector + trivector)
// - MVec4ds<T> (full multivector)
// - DualNum4ds<T> (dual numbers: scalar + pseudoscalar)
//
// **Operations available:**
// - All 4ds STA operations (geometric product, wedge, commutator, ... products)
// - Transformations (translation, rotation, reflection, scaling)
// - Motor algebra for rigid body transformations
//
// **User types (based on value_t):**
// - scalar4ds, vec4ds, bivec4ds, trivec4ds, pscalar4ds
// - mvec4ds_e, mvec4ds_u, mvec4ds, dualnum4ds
//
// **Usage:**
// #include "ga/ga_pga4ds_only.hpp"
// using namespace hd::ga;       // for types
// using namespace hd::ga::sta;  // for operations
/////////////////////////////////////////////////////////////////////////////////////////