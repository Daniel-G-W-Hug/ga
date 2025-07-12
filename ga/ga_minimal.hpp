#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

/////////////////////////////////////////////////////////////////////////////////////////
// Minimal GA library include - core types without operations
// Provides all GA types for data storage and basic manipulation without
// the computational overhead of including all operations
/////////////////////////////////////////////////////////////////////////////////////////

#include "detail/ga_core_types.hpp"
#include "detail/ga_foundation.hpp"

// All concrete type definitions
#include "detail/type_t/ga_type2d.hpp" // ega2d types
#include "detail/type_t/ga_type3d.hpp" // ega3d types
#include "detail/type_t/ga_type4d.hpp" // ega4d types

#include "detail/type_t/ga_type2dp.hpp" // pga2dp types
#include "detail/type_t/ga_type3dp.hpp" // pga3dp types

// Complete multivector types (include the even and uneven grade multivectors)
#include "detail/ga_mvec2d.hpp"
#include "detail/ga_mvec2dp.hpp"
#include "detail/ga_mvec3d.hpp"
#include "detail/ga_mvec3dp.hpp"
#include "detail/ga_mvec4d.hpp"

// User convenience types and constants
#include "ga_usr_consts.hpp"
#include "ga_usr_types.hpp"
#include "ga_usr_utilities.hpp"

// Formatting support
#include "detail/fmt/ga_fmt_core.hpp"
#include "detail/fmt/ga_fmt_ega.hpp"
#include "detail/fmt/ga_fmt_pga.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// This header provides the complete GA type system WITHOUT operations.
// Perfect for:
//
// **Data structures and storage:**
// - Storing GA objects in containers
// - Serialization/deserialization
// - Data exchange between modules
// - Type definitions in header files
//
// **Available types:**
//
// **EGA types:** scalar2d, vec2d, mvec2d, scalar3d, vec3d, bivec3d, mvec3d,
//                scalar4d, vec4d, bivec4d, trivec4d, mvec4d
//
// **PGA types:** scalar2dp, vec2dp, bivec2dp, mvec2dp, dualnum2dp,
//                scalar3dp, vec3dp, bivec3dp, trivec3dp, mvec3dp, dualnum3dp
//
// **PGA objects:** point2d, point2dp, vector2d, line2d,
//                  point3d, point3dp, vector3d, line3d, plane3d
//
// **What's excluded:**
// - All arithmetic operations (+, -, *, ^, etc.)
// - All GA-specific operations (geometric product, wedge, contractions)
// - All geometric operations (norm, unitize, exp, log, etc.)
// - All transformation operations
//
// **Compilation benefits:**
// - Fastest possible compilation (~60-70% faster than full headers)
// - Minimal template instantiation
// - Smallest object files
// - Ideal for header-only interface definitions
//
// **Usage scenario:**
// ```cpp
// // header.hpp
// #include "ga/ga_minimal.hpp"
// struct MyClass {
//     hd::ga::vec3d position;
//     hd::ga::bivec3d orientation;
//     hd::ga::point3d target;
// };
//
// // implementation.cpp
// #include "ga/ga_ega3d_only.hpp"  // Now include operations
// #include "ga/ga_pga3dp_only.hpp"
// #include "header.hpp"
// void MyClass::update() {
//     auto direction = target - position;  // Now operations are available
//     // ... rest of implementation
// }
// ```
//
// **When to use this header:**
// - Interface definitions
// - Data structure declarations
// - When you only need type definitions
// - When compilation speed is critical
// - In template metaprogramming contexts
/////////////////////////////////////////////////////////////////////////////////////////