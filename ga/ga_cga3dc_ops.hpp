#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_cga3dc_ops_products.hpp"


namespace hd::ga::cga {

////////////////////////////////////////////////////////////////////////////////
// provides cga3dc geometric operations (layer under construction):
//
// - is_congruent()          -> same subspace up to a non-zero scalar factor
// - is_close()              -> same value within a RELATIVE tolerance
//
// exp/log/sqrt, the conformal transformations, the object constructions and
// the bulk/weight/norm story follow with the cga3dc geometry phase (they
// mirror the cga2dc layer one dimension up)
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Scalar3dc<T> a, Scalar3dc<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero scalars are congruent (represent the same 0-dimensional subspace)
    return true;
}

// For vectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Vec3dc<T> const& a, Vec3dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<5>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w), value_t(a.u)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w), value_t(b.u)},
        tolerance);
}

// For bivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(BiVec3dc<T> const& a, BiVec3dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<10>(
        {value_t(a.vx), value_t(a.vy), value_t(a.vz), value_t(a.mx), value_t(a.my),
         value_t(a.mz), value_t(a.px), value_t(a.py), value_t(a.pz), value_t(a.pw)},
        {value_t(b.vx), value_t(b.vy), value_t(b.vz), value_t(b.mx), value_t(b.my),
         value_t(b.mz), value_t(b.px), value_t(b.py), value_t(b.pz), value_t(b.pw)},
        tolerance);
}

// For trivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(TriVec3dc<T> const& a, TriVec3dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<10>(
        {value_t(a.vx), value_t(a.vy), value_t(a.vz), value_t(a.mx), value_t(a.my),
         value_t(a.mz), value_t(a.px), value_t(a.py), value_t(a.pz), value_t(a.pw)},
        {value_t(b.vx), value_t(b.vy), value_t(b.vz), value_t(b.mx), value_t(b.my),
         value_t(b.mz), value_t(b.px), value_t(b.py), value_t(b.pz), value_t(b.pw)},
        tolerance);
}

// For quadvectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(QuadVec3dc<T> const& a, QuadVec3dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<5>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w), value_t(a.u)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w), value_t(b.u)},
        tolerance);
}

// For pseudoscalars: all non-zero pseudoscalars represent the same subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(PScalar3dc<T> a, PScalar3dc<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars are congruent
    return true;
}


////////////////////////////////////////////////////////////////////////////////
// is_close(): equality within a RELATIVE tolerance
//
// operator== compares against an absolute eps and cannot resolve anything once
// coordinates are large; is_close scales the per-component threshold by the
// operands (see detail::coeffs_close). Reach for it whenever the values carry a
// physical scale. Distinct from is_congruent(), which allows an arbitrary scale
// factor between the operands.
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Scalar3dc<T> a, Scalar3dc<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Vec3dc<T> const& a, Vec3dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<5>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w), value_t(a.u)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w), value_t(b.u)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(BiVec3dc<T> const& a, BiVec3dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<10>(
        {value_t(a.vx), value_t(a.vy), value_t(a.vz), value_t(a.mx), value_t(a.my),
         value_t(a.mz), value_t(a.px), value_t(a.py), value_t(a.pz), value_t(a.pw)},
        {value_t(b.vx), value_t(b.vy), value_t(b.vz), value_t(b.mx), value_t(b.my),
         value_t(b.mz), value_t(b.px), value_t(b.py), value_t(b.pz), value_t(b.pw)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(TriVec3dc<T> const& a, TriVec3dc<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<10>(
        {value_t(a.vx), value_t(a.vy), value_t(a.vz), value_t(a.mx), value_t(a.my),
         value_t(a.mz), value_t(a.px), value_t(a.py), value_t(a.pz), value_t(a.pw)},
        {value_t(b.vx), value_t(b.vy), value_t(b.vz), value_t(b.mx), value_t(b.my),
         value_t(b.mz), value_t(b.px), value_t(b.py), value_t(b.pz), value_t(b.pw)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(QuadVec3dc<T> const& a, QuadVec3dc<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<5>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w), value_t(a.u)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w), value_t(b.u)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(PScalar3dc<T> a, PScalar3dc<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec3dc_E<T> const& a, MVec3dc_E<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<16>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7), value_t(a.c8), value_t(a.c9),
         value_t(a.c10), value_t(a.c11), value_t(a.c12), value_t(a.c13), value_t(a.c14),
         value_t(a.c15)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7), value_t(b.c8), value_t(b.c9),
         value_t(b.c10), value_t(b.c11), value_t(b.c12), value_t(b.c13), value_t(b.c14),
         value_t(b.c15)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec3dc_U<T> const& a, MVec3dc_U<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<16>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7), value_t(a.c8), value_t(a.c9),
         value_t(a.c10), value_t(a.c11), value_t(a.c12), value_t(a.c13), value_t(a.c14),
         value_t(a.c15)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7), value_t(b.c8), value_t(b.c9),
         value_t(b.c10), value_t(b.c11), value_t(b.c12), value_t(b.c13), value_t(b.c14),
         value_t(b.c15)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec3dc<T> const& a, MVec3dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<32>(
        {value_t(a.c0),  value_t(a.c1),  value_t(a.c2),  value_t(a.c3),  value_t(a.c4),
         value_t(a.c5),  value_t(a.c6),  value_t(a.c7),  value_t(a.c8),  value_t(a.c9),
         value_t(a.c10), value_t(a.c11), value_t(a.c12), value_t(a.c13), value_t(a.c14),
         value_t(a.c15), value_t(a.c16), value_t(a.c17), value_t(a.c18), value_t(a.c19),
         value_t(a.c20), value_t(a.c21), value_t(a.c22), value_t(a.c23), value_t(a.c24),
         value_t(a.c25), value_t(a.c26), value_t(a.c27), value_t(a.c28), value_t(a.c29),
         value_t(a.c30), value_t(a.c31)},
        {value_t(b.c0),  value_t(b.c1),  value_t(b.c2),  value_t(b.c3),  value_t(b.c4),
         value_t(b.c5),  value_t(b.c6),  value_t(b.c7),  value_t(b.c8),  value_t(b.c9),
         value_t(b.c10), value_t(b.c11), value_t(b.c12), value_t(b.c13), value_t(b.c14),
         value_t(b.c15), value_t(b.c16), value_t(b.c17), value_t(b.c18), value_t(b.c19),
         value_t(b.c20), value_t(b.c21), value_t(b.c22), value_t(b.c23), value_t(b.c24),
         value_t(b.c25), value_t(b.c26), value_t(b.c27), value_t(b.c28), value_t(b.c29),
         value_t(b.c30), value_t(b.c31)},
        rel_tol);
}

} // namespace hd::ga::cga
