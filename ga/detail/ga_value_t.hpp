#pragma once

// author: Daniel Hug, 2024

#include <limits> // std::numeric_limits

#include "ga_strong_t.hpp"

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

// definition of scalar and pscalar located here, because other
// modules like vec2d, vec3d, 2d_ops, 3d_ops etc. just need the type defintion,
// but not the functionality itself (reduce dependencies between files)

struct scalar_tag {};
struct pscalar2d_tag {};
struct pscalar3d_tag {};

template <typename T> using Scalar = Strong_t<T, scalar_tag>;
template <typename T> using PScalar2d = Strong_t<T, pscalar2d_tag>;
template <typename T> using PScalar3d = Strong_t<T, pscalar3d_tag>;

} // namespace hd::ga