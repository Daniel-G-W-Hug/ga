#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

////////////////////////////////////////////////////////////////////////////////
// EGA (Euclidean Geometric Algebra) type definitions
// Complete type system for 2D, 3D, and 4D Euclidean geometric algebras
////////////////////////////////////////////////////////////////////////////////

#include "ga_core_types.hpp"

// EGA algebra-specific type definitions
#include "type_t/ga_type_2d.hpp"        // 2D EGA concrete types (Scalar2d, Vec2d, etc.)
#include "type_t/ga_type_3d.hpp"        // 3D EGA concrete types (Scalar3d, Vec3d, etc.) 
#include "type_t/ga_type_4d.hpp"        // 4D EGA concrete types (Scalar4d, Vec4d, etc.)

// EGA multivector implementations
#include "ga_mvec2d_e.hpp"              // 2D even subalgebra multivectors
#include "ga_mvec3d_e.hpp"              // 3D even subalgebra multivectors  
#include "ga_mvec3d_u.hpp"              // 3D odd subalgebra multivectors
#include "ga_mvec4d_e.hpp"              // 4D even subalgebra multivectors
#include "ga_mvec4d_u.hpp"              // 4D odd subalgebra multivectors

// Complete EGA multivector types
#include "ga_mvec2d.hpp"                // Full 2D multivectors (includes _e)
#include "ga_mvec3d.hpp"                // Full 3D multivectors (includes _e, _u)
#include "ga_mvec4d.hpp"                // Full 4D multivectors (includes _e, _u)

////////////////////////////////////////////////////////////////////////////////
// This header provides complete EGA type system including:
//
// **2D EGA (G₂ = G(2,0,0)):**
// - Scalar2d<T>, Vec2d<T>, PScalar2d<T>
// - MVec2d_E<T> (even subalgebra)
// - MVec2d<T> (full multivector)
//
// **3D EGA (G₃ = G(3,0,0)):**  
// - Scalar3d<T>, Vec3d<T>, BiVec3d<T>, PScalar3d<T>
// - MVec3d_E<T> (even subalgebra), MVec3d_U<T> (odd subalgebra)
// - MVec3d<T> (full multivector)
//
// **4D EGA (G₄ = G(4,0,0)):**
// - Scalar4d<T>, Vec4d<T>, BiVec4d<T>, TriVec4d<T>, PScalar4d<T>
// - MVec4d_E<T> (even subalgebra), MVec4d_U<T> (odd subalgebra)  
// - MVec4d<T> (full multivector)
//
// Usage: Include this header when you need complete EGA functionality.
// For operations files: #include "detail/ga_ega_types.hpp"
////////////////////////////////////////////////////////////////////////////////