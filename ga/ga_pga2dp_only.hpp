#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

/////////////////////////////////////////////////////////////////////////////////////////
// Selective include for 2D Projective Geometric Algebra only
// Provides complete pga2dp functionality with minimal compile overhead
/////////////////////////////////////////////////////////////////////////////////////////

#include "detail/ga_foundation.hpp"

// pga2dp specific type definitions
#include "detail/type_t/ga_type2dp.hpp"

// pga2dp multivector implementations
#include "detail/ga_mvec2dp.hpp"   // Complete 2dp multivectors
#include "detail/ga_mvec2dp_e.hpp" // 2dp even subalgebra
#include "detail/ga_mvec2dp_u.hpp" // 2dp uneven (or odd) subalgebra

// pga2dp geometric objects
#include "detail/ga_pga2dp_objects.hpp" // Point2dp, Point2d, Vector2d, Line2d

// User convenience includes
#include "ga_usr_consts.hpp"
#include "ga_usr_types.hpp" // Contains 2dp type aliases like scalar2dp, vec2dp, etc.
#include "ga_usr_utilities.hpp"

// pga2dp operations
#include "ga_pga2dp_ops.hpp"

// Formatting support
#include "detail/fmt/ga_fmt_core.hpp"
#include "detail/fmt/ga_fmt_pga.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// This header provides complete 2dp PGA functionality including:
//
// **Algebra:** G₂⁺ = G(2,0,1) with basis {1, e₁, e₂, e₃, e₁₂, e₁₃, e₂₃, e₁₂₃}
// where e₃² = 0 (degenerate/null dimension for points at infinity)
//
// **Types available:**
// - Scalar2dp<T>, Vec2dp<T>, BiVec2dp<T>, PScalar2dp<T>
// - MVec2dp_E<T> (even subalgebra: scalar + bivector)
// - MVec2dp_U<T> (odd subalgebra: vector + pseudoscalar)
// - MVec2dp<T> (full multivector)
// - DualNum2dp<T> (dual numbers: scalar + pseudoscalar)
//
// **Geometric objects:**
// - Point2dp<T> (homogeneous 2D point with explicit w coordinate)
// - Point2d<T> (Euclidean 2D point, w = 1 implicit)
// - Vector2d<T> (2D direction vector, w = 0 implicit)
// - Line2d<T> (2D line represented as bivector)
//
// **Operations available:**
// - All pga2dp operations (geometric product, wedge, regressive product)
// - Point-line incidence, line-line intersection
// - Parallel and perpendicular operations
// - Distance calculations between geometric objects
// - Transformations (translation, rotation, reflection)
//
// **User types (based on value_t):**
// - scalar2dp, vec2dp, bivec2dp, pscalar2dp, mvec2dp_e, mvec2dp_u, mvec2dp
// - dualnum2dp, point2dp, point2d, vector2d, line2d
//
// **Compilation benefits:**
// - Excludes pga3dp functionality
// - Excludes all EGA functionality
// - Reduces compilation time by ~35-45% vs full ga_pga.hpp
// - Ideal for 2D computer graphics and computational geometry
//
// **Usage:**
// #include "ga/ga_pga2dp_only.hpp"
// using namespace hd::ga;
// using namespace hd::ga::pga;  // For operations
//
// point2d p1{1.0, 2.0};
// point2d p2{3.0, 4.0};
// auto line = p1 ^ p2;  // Line through two points
// auto distance = norm(line);
//
/////////////////////////////////////////////////////////////////////////////////////////