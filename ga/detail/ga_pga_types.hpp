#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

////////////////////////////////////////////////////////////////////////////////
// PGA (Projective Geometric Algebra) type definitions
// Complete type system for 2dp and 3dp projective geometric algebras
////////////////////////////////////////////////////////////////////////////////

#include "ga_core_types.hpp"

// PGA algebra-specific type definitions
#include "type_t/ga_type_2dp.hpp" // 2dp PGA concrete types (Scalar2dp, Vec2dp, etc.)
#include "type_t/ga_type_3dp.hpp" // 3dp PGA concrete types (Scalar3dp, Vec3dp, etc.)

// PGA multivector implementations
#include "ga_mvec2dp_e.hpp" // 2dp even subalgebra multivectors
#include "ga_mvec2dp_u.hpp" // 2dp odd subalgebra multivectors
#include "ga_mvec3dp_e.hpp" // 3dp even subalgebra multivectors
#include "ga_mvec3dp_u.hpp" // 3dp odd subalgebra multivectors

// Complete PGA multivector types
#include "ga_mvec2dp.hpp" // Full 2dp multivectors (includes _e, _u)
#include "ga_mvec3dp.hpp" // Full 3dp multivectors (includes _e, _u)

// PGA geometric objects
#include "ga_pga_2dp_objects.hpp" // 2dp geometric objects (Point2dp, Line2d, etc.)
#include "ga_pga_3dp_objects.hpp" // 3dp geometric objects (Point3dp, Line3d, Plane3d, etc.)

////////////////////////////////////////////////////////////////////////////////
// This header provides complete PGA type system including:
//
// **2dp PGA (G₂⁺ = G(2,0,1)):**
// - Scalar2dp<T>, Vec2dp<T>, BiVec2dp<T>, PScalar2dp<T>
// - MVec2dp_E<T> (even subalgebra), MVec2dp_U<T> (odd subalgebra)
// - MVec2dp<T> (full multivector)
// - DualNum2dp<T> (dual numbers for 2dp space)
// - Geometric objects: Point2dp<T>, Point2d<T>, Vector2d<T>, Line2d<T>
//
// **3dp PGA (G₃⁺ = G(3,0,1)):**
// - Scalar3dp<T>, Vec3dp<T>, BiVec3dp<T>, TriVec3dp<T>, PScalar3dp<T>
// - MVec3dp_E<T> (even subalgebra), MVec3dp_U<T> (odd subalgebra)
// - MVec3dp<T> (full multivector)
// - DualNum3dp<T> (dual numbers for 3dp space)
// - Geometric objects: Point3dp<T>, Point3d<T>, Vector3d<T>, Line3d<T>, Plane3d<T>
//
// The projective algebras add one degenerate dimension (squares to 0) which
// enables elegant handling of points at infinity and parallel operations.
//
// Usage: Include this header when you need complete PGA functionality.
// For operations files: #include "detail/ga_pga_types.hpp"
////////////////////////////////////////////////////////////////////////////////