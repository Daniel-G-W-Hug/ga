#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

/////////////////////////////////////////////////////////////////////////////////////////
// Selective include for 2D Euclidean Geometric Algebra only
// Provides complete ega2d functionality with minimal compile overhead
/////////////////////////////////////////////////////////////////////////////////////////

#include "detail/ga_foundation.hpp"

// ega2d specific type definitions
#include "detail/type_t/ga_type2d.hpp"

// ega2d multivector implementations
#include "detail/type_t/ga_mvec2d.hpp"   // Complete 2D multivectors
#include "detail/type_t/ga_mvec2d_e.hpp" // 2D even subalgebra

// User convenience includes
#include "ga_usr_consts.hpp"
#include "ga_usr_types.hpp" // Contains type aliases like scalar2d, vec2d, etc.
#include "ga_usr_utilities.hpp"

// ega2d operations (includes basics and products)
#include "ga_ega2d_ops.hpp"

// Formatting support
#include "detail/fmt/ga_fmt_core.hpp"
#include "detail/fmt/ga_fmt_ega.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// This header provides complete 2D EGA functionality including:
//
// **Types available:**
// - Scalar2d<T>, Vec2d<T>, PScalar2d<T>
// - MVec2d_E<T> (even subalgebra: scalar + pseudoscalar)
// - MVec2d<T> (full multivector: scalar + vector + pseudoscalar)
//
// **Operations available:**
// - All 2D geometric algebra operations (geometric product, wedge, dot, etc.)
// - Normalization, unitization, inverse operations
// - Trigonometric operations for rotations in 2D
//
// **User types (based on value_t):**
// - scalar2d, vec2d, pscalar2d, mvec2d_e, mvec2d
//
// **Compilation benefits:**
// - Excludes 3D and 4D types and operations
// - Excludes all PGA functionality
// - Reduces compilation time by ~30-40% vs full ga_ega.hpp
// - Smaller object files and faster linking
//
// **Usage:**
// #include "ga/ga_ega2d_only.hpp"
// using namespace hd::ga;
// using namespace hd::ga::ega;  // for operations
//
// vec2d v{1.0, 2.0};
// auto length = norm(v);
//
/////////////////////////////////////////////////////////////////////////////////////////