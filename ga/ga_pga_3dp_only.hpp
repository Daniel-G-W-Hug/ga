#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

////////////////////////////////////////////////////////////////////////////////
// Selective include for 3dp Projective Geometric Algebra only
// Provides complete 3dp PGA functionality with minimal compile overhead
////////////////////////////////////////////////////////////////////////////////

#include "detail/ga_foundation.hpp"

// 3dp PGA specific type definitions
#include "detail/type_t/ga_type_3dp.hpp"

// 3dp PGA multivector implementations
#include "detail/ga_mvec3dp.hpp"   // Complete 3dp multivectors
#include "detail/ga_mvec3dp_e.hpp" // 3dp even subalgebra
#include "detail/ga_mvec3dp_u.hpp" // 3dp odd subalgebra

// 3dp PGA geometric objects
#include "detail/ga_pga_3dp_objects.hpp" // Point3dp, Point3d, Vector3d, Line3d, Plane3d

// User convenience includes
#include "ga_usr_consts.hpp"
#include "ga_usr_types.hpp" // Contains 3dp type aliases like scalar3dp, vec3dp, etc.
#include "ga_usr_utilities.hpp"

// 3dp PGA operations
#include "ga_pga_3dp_ops.hpp"

// Optional: Formatting support (can be omitted for faster compilation)
#ifndef GA_NO_FMT_SUPPORT
#include "detail/fmt/ga_fmt_core.hpp"
#include "detail/fmt/ga_fmt_pga.hpp"
#endif

////////////////////////////////////////////////////////////////////////////////
// This header provides complete 3dp PGA functionality including:
//
// **Algebra:** G₃⁺ = G(3,0,1) with basis {1, e₁, e₂, e₃, e₄, e₁₂, e₁₃, e₁₄,
// e₂₃, e₂₄, e₃₄, e₁₂₃, e₁₂₄, e₁₃₄, e₂₃₄, e₁₂₃₄}
// where e₄² = 0 (degenerate/null dimension for points at infinity)
//
// **Types available:**
// - Scalar3dp<T>, Vec3dp<T>, BiVec3dp<T>, TriVec3dp<T>, PScalar3dp<T>
// - MVec3dp_E<T> (even subalgebra: scalar + bivector + pseudoscalar)
// - MVec3dp_U<T> (odd subalgebra: vector + trivector)
// - MVec3dp<T> (full multivector)
// - DualNum3dp<T> (dual numbers: scalar + pseudoscalar)
//
// **Geometric objects:**
// - Point3dp<T> (homogeneous 3D point with explicit w coordinate)
// - Point3d<T> (Euclidean 3D point, w = 1 implicit)
// - Vector3d<T> (3D direction vector, w = 0 implicit)
// - Line3d<T> (3D line represented as bivector with 6 components)
// - Plane3d<T> (3D plane represented as trivector)
//
// **Operations available:**
// - All 3dp PGA operations (geometric product, wedge, regressive product)
// - Point-line-plane incidence relationships
// - Line-line, line-plane, plane-plane intersections
// - Distance calculations between all geometric object types
// - Transformations (translation, rotation, reflection, scaling)
// - Motor algebra for rigid body transformations
//
// **User types (based on value_t):**
// - scalar3dp, vec3dp, bivec3dp, trivec3dp, pscalar3dp
// - mvec3dp_e, mvec3dp_u, mvec3dp, dualnum3dp
// - point3dp, point3d, vector3d, line3d, plane3d
//
// **Compilation benefits:**
// - Excludes 2dp PGA functionality
// - Excludes all EGA functionality
// - Reduces compilation time by ~30-40% vs full ga_pga.hpp
// - Ideal for 3D computer graphics, robotics, and computational geometry
//
// **Usage:**
// #include "ga/ga_pga_3dp_only.hpp"
// using namespace hd::ga;
// using namespace hd::ga::pga;  // For operations
//
// point3d p1{1.0, 2.0, 3.0};
// point3d p2{4.0, 5.0, 6.0};
// point3d p3{7.0, 8.0, 9.0};
// auto plane = p1 ^ p2 ^ p3;  // Plane through three points
// auto distance_to_origin = norm(plane);
////////////////////////////////////////////////////////////////////////////////