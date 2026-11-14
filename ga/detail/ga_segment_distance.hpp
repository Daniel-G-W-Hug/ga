#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <algorithm> // std::clamp, std::max
#include <array>     // std::array
#include <cmath>     // std::sqrt
#include <cstddef>   // size_t
#include <limits>    // std::numeric_limits
#include <vector>    // std::vector

/////////////////////////////////////////////////////////////////////////////////////////
// provides the arithmetic two distances of bounded objects share in every dimension:
//
// - detail::segment_distance()   -> the distance of two segments, from their end points'
//                                   coordinates
// - detail::envelope_min()       -> the least value of the upper envelope of straight
//                                   lines over the unit interval
//
// Both are plain numerics on purpose: the clamping to a segment's ends and the walk
// along an envelope are case splits on scalars in any algebra. The typed functions that
// use them are segment_distance() and convex_distance() of pga2dp and pga3dp.
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga::detail {

// THE DISTANCE OF TWO SEGMENTS, a1 b1 and a2 b2, in N dimensions: the closest points of
// the two lines, clamped to the segments (Ericson, Real-Time Collision Detection,
// 5.1.9). Exact, including parallel, coincident and zero-length segments.
//
// A segment counts as a point when its squared length is below machine epsilon squared
// times the largest squared length in the problem, and two lines as parallel when the
// determinant of their directions is lost in the rounding of its two terms -- both
// relative, so the answer does not depend on the unit of length.
template <typename T, size_t N>
T segment_distance(std::array<T, N> const& a1, std::array<T, N> const& b1,
                   std::array<T, N> const& a2, std::array<T, N> const& b2)
{
    std::array<T, N> d1, d2, r;
    for (size_t i = 0; i < N; ++i) {
        d1[i] = b1[i] - a1[i];
        d2[i] = b2[i] - a2[i];
        r[i] = a1[i] - a2[i];
    }
    auto dot = [](std::array<T, N> const& u, std::array<T, N> const& v) {
        T sum = T(0.0);
        for (size_t i = 0; i < N; ++i)
            sum += u[i] * v[i];
        return sum;
    };
    T const eps_m = std::numeric_limits<T>::epsilon();
    T const a = dot(d1, d1), e = dot(d2, d2), f = dot(d2, r);
    T const tiny = eps_m * eps_m * std::max({a, e, dot(r, r)});
    T s = T(0.0), t = T(0.0);
    if (a <= tiny && e <= tiny) {
        // two points
    }
    else if (a <= tiny) t = std::clamp(f / e, T(0.0), T(1.0));
    else {
        T const c = dot(d1, r);
        if (e <= tiny) s = std::clamp(-c / a, T(0.0), T(1.0));
        else {
            // parallel lines (the determinant zero against its terms): any s serves, 0
            T const b = dot(d1, d2), denom = a * e - b * b;
            if (denom > T(4.5e3) * eps_m * a * e)
                s = std::clamp((b * f - c * e) / denom, T(0.0), T(1.0));
            t = (b * s + f) / e;
            if (t < T(0.0)) {
                t = T(0.0);
                s = std::clamp(-c / a, T(0.0), T(1.0));
            }
            else if (t > T(1.0)) {
                t = T(1.0);
                s = std::clamp((b - c) / a, T(0.0), T(1.0));
            }
        }
    }
    T sum = T(0.0);
    for (size_t i = 0; i < N; ++i) {
        T const w = r[i] + s * d1[i] - t * d2[i];
        sum += w * w;
    }
    return std::sqrt(sum);
}

// THE LEAST VALUE OF AN UPPER ENVELOPE: of max_i (a[i] + t s[i]) over t in [0, 1], for n
// straight lines given by their values at t = 0 and their slopes. The envelope is convex
// and piecewise linear; it is walked from t = 0 along the line on top to the crossing
// where a steeper one takes over, until the line on top no longer falls. Exact, nothing
// is sampled. `t_out`: where the least value is; `which_out`: the line on top there.
// Without lines the value is -infinity.
template <typename T>
T envelope_min(T const* a, T const* s, size_t n, T* t_out = nullptr,
               size_t* which_out = nullptr)
{
    if (n == 0) {
        if (t_out) *t_out = T(0.0);
        return -std::numeric_limits<T>::infinity();
    }
    // the line on top at t = 0; of equal ones the steepest, which is on top just after
    size_t top = 0;
    for (size_t i = 1; i < n; ++i)
        if (a[i] > a[top] || (a[i] == a[top] && s[i] > s[top])) top = i;
    T t = T(0.0);
    for (size_t guard = 0; guard <= n; ++guard) {
        if (s[top] >= T(0.0)) break; // the envelope no longer falls: the least is here
        // the first crossing after t at which a steeper line takes over
        T t_next = T(1.0);
        size_t next = top;
        for (size_t j = 0; j < n; ++j) {
            if (s[j] <= s[top]) continue;
            T const tj = (a[top] - a[j]) / (s[j] - s[top]);
            if (tj < t) continue;
            if (tj < t_next || (tj == t_next && next != top && s[j] > s[next])) {
                t_next = tj;
                next = j;
            }
        }
        t = t_next;
        if (next == top) break; // none before the end: the least is at t = 1
        top = next;
    }
    if (t_out) *t_out = t;
    if (which_out) *which_out = top;
    return a[top] + t * s[top];
}

// envelope_min() for lines filled in by `fill(a, s)`: on the stack for up to 16 lines
// (a box has six faces), on the heap beyond.
template <typename T, typename Fill>
T envelope_min_of(size_t n, Fill&& fill, T* t_out = nullptr, size_t* which_out = nullptr)
{
    constexpr size_t on_stack = 16;
    T a_few[on_stack], s_few[on_stack];
    std::vector<T> many;
    T* a = a_few;
    T* s = s_few;
    if (n > on_stack) {
        many.resize(2 * n);
        a = many.data();
        s = a + n;
    }
    fill(a, s);
    return envelope_min(a, s, n, t_out, which_out);
}

} // namespace hd::ga::detail
