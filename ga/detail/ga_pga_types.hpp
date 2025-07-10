#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

/////////////////////////////////////////////////////////////////////////////////////////
// PGA (Projective Geometric Algebra) type definitions
// Complete type system for 2dp and 3dp projective geometric algebras
/////////////////////////////////////////////////////////////////////////////////////////

#include "ga_core_types.hpp"

// PGA algebra-specific type definitions
#include "type_t/ga_type2dp.hpp" // pga2dp concrete types (Scalar2dp, Vec2dp, etc.)
#include "type_t/ga_type3dp.hpp" // pga3dp concrete types (Scalar3dp, Vec3dp, etc.)

// Complete PGA multivector types (includes the even and uneven (odd) subtypes)
#include "type_t/ga_mvec2dp.hpp" // pga2dp multivectors incl. even/uneven mvec2dp_e/u
#include "type_t/ga_mvec3dp.hpp" // pga3dp multivectors incl. even/uneven mvec2dp_e/u


/////////////////////////////////////////////////////////////////////////////////////////
// This header provides complete PGA type system including:
//
// **pga2dp (G₂⁺ = G(2,0,1)):**
//
// - Scalar2dp<T>, Vec2dp<T>, BiVec2dp<T>, PScalar2dp<T>
// - MVec2dp_E<T> (even subalgebra), MVec2dp_U<T> (uneven (or odd) subalgebra)
// - MVec2dp<T> (full multivector)
// - DualNum2dp<T> (dual numbers for 2dp space)
// - Geometric objects: Point2dp<T>, Point2d<T>, Vector2d<T>, Line2d<T>
//
//
// **pga3dp (G₃⁺ = G(3,0,1)):**
//
// - Scalar3dp<T>, Vec3dp<T>, BiVec3dp<T>, TriVec3dp<T>, PScalar3dp<T>
// - MVec3dp_E<T> (even subalgebra), MVec3dp_U<T> (uneven (or odd) subalgebra)
// - MVec3dp<T> (full multivector)
// - DualNum3dp<T> (dual numbers for 3dp space)
// - Geometric objects: Point3dp<T>, Point3d<T>, Vector3d<T>, Line3d<T>, Plane3d<T>
//
//
// The projective algebras add one degenerate dimension (squares to 0) which
// enables elegant handling of points at infinity.
//
// Usage: Include this header when you need complete PGA functionality.
// For operations files: #include "detail/ga_pga_types.hpp"
/////////////////////////////////////////////////////////////////////////////////////////