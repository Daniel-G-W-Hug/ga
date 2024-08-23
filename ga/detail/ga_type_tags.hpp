#pragma once

// author: Daniel Hug, 2024

////////////////////////////////////////////////////////////////////////////////
// type tags for consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

struct default_tag {}; // used as reference for inheriting ctors from base types

struct scalar_tag {}; // used for each algebra regardless of dimension

struct vec2d_tag {};
struct pscalar2d_tag {};
struct mvec2d_tag {};
struct mvec2d_e_tag {};

struct vec2dp_tag {};
struct bivec2dp_tag {};
struct pscalar2dp_tag {};

struct vec3d_tag {};
struct bivec3d_tag {};
struct pscalar3d_tag {};
struct mvec3d_tag {};
struct mvec3d_e_tag {};
struct mvec3d_u_tag {};

struct vec3dp_tag {};
struct bivec3dp_tag {};
struct trivec3dp_tag {};
struct pscalar3dp_tag {};

struct vec4d_tag {};
struct bivec4d_tag {};
struct trivec4d_tag {};
struct pscalar4d_tag {};
struct mvec4d_tag {};

} // namespace hd::ga
