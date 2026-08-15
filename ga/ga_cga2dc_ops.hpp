#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_cga2dc_ops_products.hpp"


namespace hd::ga::cga {

/////////////////////////////////////////////////////////////////////////////////////////
// provides cga2dc geometric operations (layer under construction):
//
// - is_congruent()          -> same subspace up to a non-zero scalar factor
// - is_close()              -> same value within a RELATIVE tolerance
//
// Still to come in this layer: bulk/weight and the norms (their design needs
// the conformal norm conventions -- null objects change the norm story vs. the
// flat algebras), up/down point embeddings, the round objects (point pair,
// circle) and flats with their property extraction (center, radius, attitude,
// carrier), meet/join, the conformal versors (translator, rotor, dilator,
// inversion, transversion) with the sandwich verb, and the versor-equality
// test.
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(Scalar2dc<T> a, Scalar2dc<U> b, value_t tolerance = eps)
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
bool is_congruent(Vec2dc<T> const& a, Vec2dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, tolerance);
}

// For bivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(BiVec2dc<T> const& a, BiVec2dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<6>({value_t(a.vx), value_t(a.vy), value_t(a.vz),
                                        value_t(a.mx), value_t(a.my), value_t(a.mz)},
                                       {value_t(b.vx), value_t(b.vy), value_t(b.vz),
                                        value_t(b.mx), value_t(b.my), value_t(b.mz)},
                                       tolerance);
}

// For trivectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(TriVec2dc<T> const& a, TriVec2dc<U> const& b,
                  value_t tolerance = eps_congruent)
{
    return detail::coeffs_congruent<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, tolerance);
}

// For pseudoscalars: all non-zero pseudoscalars represent the same subspace
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_congruent(PScalar2dc<T> a, PScalar2dc<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars are congruent (represent the full 4D space)
    return true;
}


////////////////////////////////////////////////////////////////////////////////
// is_close(): equality within a RELATIVE tolerance
////////////////////////////////////////////////////////////////////////////////
//
// Same question as operator==, but with the threshold scaled by the operands (see
// detail::coeffs_close). Reach for it whenever the values carry a physical scale:
// operator== measures against an absolute eps, which cannot resolve anything once
// coordinates grow large, where a single ulp already exceeds it. In a conformal
// algebra that regime is the default: the round (u) component of an embedded
// point grows QUADRATICALLY with the coordinates.
//
// Distinct from is_congruent(), which allows an arbitrary scale factor between the
// operands and so answers "same subspace", not "same value".

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Scalar2dc<T> a, Scalar2dc<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(Vec2dc<T> const& a, Vec2dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(BiVec2dc<T> const& a, BiVec2dc<U> const& b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<6>({value_t(a.vx), value_t(a.vy), value_t(a.vz),
                                    value_t(a.mx), value_t(a.my), value_t(a.mz)},
                                   {value_t(b.vx), value_t(b.vy), value_t(b.vz),
                                    value_t(b.mx), value_t(b.my), value_t(b.mz)},
                                   rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(TriVec2dc<T> const& a, TriVec2dc<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<4>(
        {value_t(a.x), value_t(a.y), value_t(a.z), value_t(a.w)},
        {value_t(b.x), value_t(b.y), value_t(b.z), value_t(b.w)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(PScalar2dc<T> a, PScalar2dc<U> b, value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<1>({value_t(a)}, {value_t(b)}, rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec2dc_E<T> const& a, MVec2dc_E<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<8>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec2dc_U<T> const& a, MVec2dc_U<U> const& b,
              value_t rel_tol = eps_congruent)
{
    return detail::coeffs_close<8>(
        {value_t(a.c0), value_t(a.c1), value_t(a.c2), value_t(a.c3), value_t(a.c4),
         value_t(a.c5), value_t(a.c6), value_t(a.c7)},
        {value_t(b.c0), value_t(b.c1), value_t(b.c2), value_t(b.c3), value_t(b.c4),
         value_t(b.c5), value_t(b.c6), value_t(b.c7)},
        rel_tol);
}

template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
bool is_close(MVec2dc<T> const& a, MVec2dc<U> const& b, value_t rel_tol = eps_congruent)
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

} // namespace hd::ga::cga
