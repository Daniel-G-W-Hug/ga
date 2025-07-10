#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

/////////////////////////////////////////////////////////////////////////////////////////
// Selective include for 3D Euclidean Geometric Algebra only
// Provides complete ega3d functionality with minimal compile overhead
/////////////////////////////////////////////////////////////////////////////////////////

#include "detail/ga_foundation.hpp"

// 3D EGA specific type definitions
#include "detail/type_t/ga_type3d.hpp"

// 3D EGA multivector implementations
#include "detail/type_t/ga_mvec3d.hpp"   // Complete 3D multivectors
#include "detail/type_t/ga_mvec3d_e.hpp" // 3D even subalgebra
#include "detail/type_t/ga_mvec3d_u.hpp" // 3D odd subalgebra

// User convenience includes
#include "ga_usr_consts.hpp"
#include "ga_usr_types.hpp" // Contains 3D type aliases like scalar3d, vec3d, etc.
#include "ga_usr_utilities.hpp"

// ega3d operations (includes basics and products)
#include "ga_ega3d_ops.hpp"

// Optional: Formatting support (can be omitted for faster compilation)
#ifndef _HD_GA_NO_FMT_SUPPORT
#include "detail/fmt/ga_fmt_core.hpp"
#include "detail/fmt/ga_fmt_ega.hpp"
#endif

////////////////////////////////////////////////////////////////////////////////
// This header provides complete 3D EGA functionality including:
//
// **Types available:**
// - Scalar3d<T>, Vec3d<T>, BiVec3d<T>, PScalar3d<T>
// - MVec3d_E<T> (even subalgebra: scalar + bivector)
// - MVec3d_U<T> (odd subalgebra: vector + pseudoscalar)
// - MVec3d<T> (full multivector: scalar + vector + bivector + pseudoscalar)
//
// **Operations available:**
// - All 3D geometric algebra operations (geometric product, wedge, dot, etc.)
// - Normalization, unitization, inverse operations
// - Rotation operations using bivectors and rotors
// - Cross product equivalents via bivector operations
//
// **User types (based on value_t):**
// - scalar3d, vec3d, bivec3d, pscalar3d, mvec3d_e, mvec3d_u, mvec3d
//
// **Compilation benefits:**
// - Excludes 2D and 4D types and operations
// - Excludes all PGA functionality
// - Reduces compilation time by ~25-35% vs full ga_ega.hpp
// - Smaller object files and faster linking
//
// **Usage:**
// #include "ga/ga_ega3d_only.hpp"
// using namespace hd::ga;
// using namespace hd::ga::ega;  // For operations
//
// vec3d v{1.0, 2.0, 3.0};
// bivec3d b{0.0, 0.0, 1.0};  // Rotation bivector around z-axis
// auto rotor = exp(b * 0.5); // Half-angle for rotation
////////////////////////////////////////////////////////////////////////////////