#pragma once

// author: Daniel Hug, 2024 & 2025

#include <limits> // std::numeric_limits

////////////////////////////////////////////////////////////////////////////////
// consistent type definitions for easy use
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

// select the floating point type used for scalars, vector and bivector components, etc.
// here, i.e. choose one of the two following definitions:
// using value_t = float;
using value_t = double;

// eps is used internally for equality testing
auto const eps{5.0 * std::numeric_limits<value_t>::epsilon()};


} // namespace hd::ga