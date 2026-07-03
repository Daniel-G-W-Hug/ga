#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// User convenience types for GA mechanics (based on value_t)
//
// Companion to ga_usr_types.hpp: collects the value_t-based aliases for the mechanics types
// in ONE findable location. They live here rather than in ga_usr_types.hpp because they
// alias templates (Inertia2dp, ...) that are only defined by the mechanics ops headers,
// which are included AFTER ga_usr_types.hpp. This header must therefore be included after
// the mechanics ops headers (see ga_pga.hpp).
/////////////////////////////////////////////////////////////////////////////////////////

#include "ga_pga2dp_ops_mechanics.hpp" // Inertia2dp<T>, kinematic frame/system types
#include "ga_pga3dp_ops_mechanics.hpp" // Inertia3dp<T>
#include "ga_value_t.hpp"            // value_t

namespace hd::ga {

using inertia2dp = pga::Inertia2dp<value_t>; // 3x3 inertia matrix for 2D rigid body
using inertia3dp = pga::Inertia3dp<value_t>; // 6x6 inertia matrix for 3D rigid body

} // namespace hd::ga
