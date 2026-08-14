#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// CGA (Conformal Geometric Algebra) type definitions
// Complete type system for the 2dc conformal geometric algebra
// (the 3dc algebra follows once its 32-component type layer exists)
/////////////////////////////////////////////////////////////////////////////////////////

#include "ga_core_types.hpp"

// CGA algebra-specific type definitions
#include "type_t/ga_type2dc.hpp" // cga2dc concrete types (Scalar2dc, Vec2dc, etc.)

// Complete CGA multivector types (includes the even and uneven (odd) subtypes)
#include "type_t/ga_mvec2dc.hpp" // cga2dc multivectors incl. even/uneven mvec2dc_e/u

/////////////////////////////////////////////////////////////////////////////////////////
// This header provides the complete cga2dc type system:
//
// **cga2dc (G(3,1,0) after diagonalization of the null-basis metric):**
//
// - Scalar2dc<T>, Vec2dc<T>, BiVec2dc<T>, TriVec2dc<T>, PScalar2dc<T>
// - MVec2dc_E<T> (even subalgebra), MVec2dc_U<T> (uneven (or odd) subalgebra)
// - MVec2dc<T> (full multivector)
// - DualNum2dc<T> (dual numbers for 2dc space)
//
// The conformal algebra adds TWO projective dimensions to 2D Euclidean space:
// the homogeneous dimension known from PGA (e3, w) and a round dimension (e4,
// u) that embeds points on a null cone — which makes circles and point pairs
// representable as blades alongside the flat objects. The null pair satisfies
// e3*e3 = e4*e4 = 0 and e3.e4 = -1 (non-orthogonal metric).
//
// Usage: Include this header when you need complete CGA type functionality.
// For operations files: #include "detail/ga_cga_types.hpp"
/////////////////////////////////////////////////////////////////////////////////////////
