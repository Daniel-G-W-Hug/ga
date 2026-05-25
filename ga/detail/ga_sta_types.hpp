#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// STA (Space-Time Algebra) type definitions
// Complete type system 4ds space-time geometric algebra
/////////////////////////////////////////////////////////////////////////////////////////

#include "ga_core_types.hpp"

// STA algebra-specific type definitions
#include "type_t/ga_type4ds.hpp" // pga4ds concrete types (Scalar4ds, Vec4ds, etc.)

// STA multivector types (includes the even and uneven (odd) subtypes)
#include "type_t/ga_mvec4ds.hpp" // pga4ds multivectors incl. even/uneven mvec2dp_e/u


/////////////////////////////////////////////////////////////////////////////////////////
// This header provides complete PGA type system including:
//
// **sta4ds with metric signature G(3,1,0)  ("mostly positive" convention)
//
// - Scalar4ds<T>, Vec4ds<T>, BiVec4ds<T>, TriVec4ds<T>, PScalar4ds<T>
// - MVec4ds_E<T> (even subalgebra), MVec4ds_U<T> (uneven (or odd) subalgebra)
// - MVec4ds<T> (full multivector)
// - DualNum4ds<T> (dual numbers for 4ds space)
//
// Usage: Include this header when you need complete STA functionality.
// For operations files: #include "detail/ga_pga_types.hpp"
/////////////////////////////////////////////////////////////////////////////////////////