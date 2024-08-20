#pragma once

// author: Daniel Hug, 2024

#include <limits> // std::numeric_limits

#include "ga_scalar_t.hpp"
#include "ga_vec_xyz_t.hpp"

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
// modules just need the type defintion, but not the functionality itself
// (reduce dependencies between files)

struct scalar_tag {};
struct pscalar2d_tag {};
struct pscalar3d_tag {};

struct vec3d_tag {};
struct bivec3d_tag {};

template <typename T> using Scalar = Scalar_t<T, scalar_tag>;
template <typename T> using PScalar2d = Scalar_t<T, pscalar2d_tag>;
template <typename T> using PScalar3d = Scalar_t<T, pscalar3d_tag>;

template <typename T> using Vec3d = Vec_xyz_t<T, vec3d_tag>;
template <typename T> using BiVec3d = Vec_xyz_t<T, bivec3d_tag>;


} // namespace hd::ga