#pragma once

// author: Daniel Hug, 2024 & 2025

#include <numbers> // math constants like pi

#include "ga_value_t.hpp"

namespace hd::ga {

using std::numbers::pi; // make pi available for users

inline constexpr value_t deg2rad(value_t value) { return value / 180.0 * pi; }
inline constexpr value_t rad2deg(value_t value) { return value * 180.0 / pi; }

inline constexpr value_t sign(value_t value) { return (value >= 0.0) ? 1.0 : -1.0; }

} // namespace hd::ga