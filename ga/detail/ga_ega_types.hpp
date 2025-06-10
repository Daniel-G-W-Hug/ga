#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

/////////////////////////////////////////////////////////////////////////////////////////
// EGA (Euclidean Geometric Algebra) type definitions
// Complete type system for 2D, 3D, and 4D Euclidean geometric algebras
/////////////////////////////////////////////////////////////////////////////////////////

#include "ga_core_types.hpp"

// EGA algebra-specific type definitions
#include "type_t/ga_type2d.hpp" // ega2d concrete types (Scalar2d, Vec2d, etc.)
#include "type_t/ga_type3d.hpp" // ega3d concrete types (Scalar3d, Vec3d, etc.)
#include "type_t/ga_type4d.hpp" // ega4d concrete types (Scalar4d, Vec4d, etc.)

// Complete EGA multivector types (includes the even and uneven (odd) subtypes)
#include "ga_mvec2d.hpp" // ega2d multivectors incl. even mvec2d_e
#include "ga_mvec3d.hpp" // ega3d multivectors incl. even mvec3d_e and uneven mvec3d_u
#include "ga_mvec4d.hpp" // ega4d multivectors incl. even mvec4d_e and uneven mvec4d_u

/////////////////////////////////////////////////////////////////////////////////////////
// This header provides complete EGA type system including:
//
// **ega2d (G₂ = G(2,0,0)):**
// - Scalar2d<T>, Vec2d<T>, PScalar2d<T>
// - MVec2d_E<T> (even subalgebra)
// - MVec2d<T> (full multivector)
//
// **ega3d (G₃ = G(3,0,0)):**
// - Scalar3d<T>, Vec3d<T>, BiVec3d<T>, PScalar3d<T>
// - MVec3d_E<T> (even subalgebra), MVec3d_U<T> (uneven (or odd) subalgebra)
// - MVec3d<T> (full multivector)
//
// **ega4d (G₄ = G(4,0,0)):**
// - Scalar4d<T>, Vec4d<T>, BiVec4d<T>, TriVec4d<T>, PScalar4d<T>
// - MVec4d_E<T> (even subalgebra), MVec4d_U<T> (uneven (or odd) subalgebra)
// - MVec4d<T> (full multivector)
//
// Usage: Include this header when you need complete EGA functionality.
// For operations files: #include "detail/ga_ega_types.hpp"
/////////////////////////////////////////////////////////////////////////////////////////