#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

////////////////////////////////////////////////////////////////////////////////
// Core type templates for GA library
// Common template types used across all geometric algebras (EGA and PGA)
////////////////////////////////////////////////////////////////////////////////

#include "ga_foundation.hpp"

// Type system infrastructure
#include "type_t/ga_type_tags.hpp"      // Type tag definitions for template specialization

// Basic template types (shared across all algebras)
#include "type_t/ga_scalar_t.hpp"       // Scalar_t<T, Tag> template
#include "type_t/ga_vec2_t.hpp"         // Vec2_t<T, Tag> template  
#include "type_t/ga_vec3_t.hpp"         // Vec3_t<T, Tag> template
#include "type_t/ga_vec4_t.hpp"         // Vec4_t<T, Tag> template
#include "type_t/ga_bvec6_t.hpp"        // BVec6_t<T, Tag> template (6-component bivectors)

// Multivector template types (shared across all algebras)
#include "type_t/ga_mvec2_t.hpp"        // MVec2_t<T, Tag> template (2-component)
#include "type_t/ga_mvec4_t.hpp"        // MVec4_t<T, Tag> template (4-component)  
#include "type_t/ga_mvec8_t.hpp"        // MVec8_t<T, Tag> template (8-component)
#include "type_t/ga_mvec16_t.hpp"       // MVec16_t<T, Tag> template (16-component)

////////////////////////////////////////////////////////////////////////////////
// This header provides the common type template infrastructure that both
// EGA and PGA algebras depend on. It includes:
//
// 1. Foundation infrastructure (value_t, error handling, algebra definitions)
// 2. Type tag system for template specialization
// 3. All basic vector and multivector template types
//
// These templates are then specialized in algebra-specific headers with
// concrete type definitions like:
// - EGA: Scalar2d<T>, Vec2d<T>, MVec2d<T>, etc.
// - PGA: Scalar2dp<T>, Vec2dp<T>, MVec2dp<T>, etc.
//
// Usage: Include this in algebra-specific aggregation headers
// (ga_ega_types.hpp, ga_pga_types.hpp) rather than directly in user code.
////////////////////////////////////////////////////////////////////////////////