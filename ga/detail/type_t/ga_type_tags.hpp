#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

////////////////////////////////////////////////////////////////////////////////
// type tags for consistent type definitions
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

struct default_tag {}; // used as reference for inheriting ctors from base types

struct ega2d_tag {};
struct scalar2d_tag : public ega2d_tag {};
struct vec2d_tag : public ega2d_tag {};
struct pscalar2d_tag : public ega2d_tag {};
struct mvec2d_e_tag : public ega2d_tag {};
struct mvec2d_tag : public ega2d_tag {};

struct ega3d_tag {};
struct scalar3d_tag : public ega3d_tag {};
struct vec3d_tag : public ega3d_tag {};
struct bivec3d_tag : public ega3d_tag {};
struct pscalar3d_tag : public ega3d_tag {};
struct mvec3d_e_tag : public ega3d_tag {};
struct mvec3d_u_tag : public ega3d_tag {};
struct mvec3d_tag : public ega3d_tag {};

struct ega4d_tag {};
struct scalar4d_tag : public ega4d_tag {};
struct vec4d_tag : public ega4d_tag {};
struct bivec4d_tag : public ega4d_tag {};
struct trivec4d_tag : public ega4d_tag {};
struct pscalar4d_tag : public ega4d_tag {};
struct mvec4d_e_tag : public ega4d_tag {};
struct mvec4d_u_tag : public ega4d_tag {};
struct mvec4d_tag : public ega4d_tag {};

struct pga2dp_tag {};
struct scalar2dp_tag : public pga2dp_tag {};
struct vec2dp_tag : public pga2dp_tag {};
struct bivec2dp_tag : public pga2dp_tag {};
struct pscalar2dp_tag : public pga2dp_tag {};
struct mvec2dp_e_tag : public pga2dp_tag {};
struct mvec2dp_u_tag : public pga2dp_tag {};
struct mvec2dp_tag : public pga2dp_tag {};
struct dual_number2dp_tag : public pga2dp_tag {};

struct pga3dp_tag {};
struct scalar3dp_tag : public pga3dp_tag {};
struct vec3dp_tag : public pga3dp_tag {};
struct bivec3dp_tag : public pga3dp_tag {};
struct trivec3dp_tag : public pga3dp_tag {};
struct pscalar3dp_tag : public pga3dp_tag {};
struct mvec3dp_e_tag : public pga3dp_tag {};
struct mvec3dp_u_tag : public pga3dp_tag {};
struct mvec3dp_tag : public pga3dp_tag {};
struct dual_number3dp_tag : public pga3dp_tag {};

} // namespace hd::ga