// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// What the division guard costs, measured rather than promised.
//
// Built TWICE from this one file (see utilities/CMakeLists.txt):
// `ga_bench_division_guard` with _HD_GA_EXTENDED_TEST_DIV_BY_ZERO as every consumer gets
// it by default, and `ga_bench_division_guard_off` through ga_division_guard(<target>
// OFF). Run both; the ratio of the two ns/op figures is the guard's cost on this compiler
// for a kernel of this shape. Two binaries because the choice is per translation unit --
// the library is header-only, so the check is instantiated where it is used.
//
// The kernel is a pairwise inverse-cube vector sum over two point sets: for every
// source i and target j, r = p_j - q_i, accumulate (m_i . rhat) rhat / |r|^3 -- one
// normalization (r / |r|, a checked division) and one scalar division per pair, inside
// a double loop with no other work. That is the shape of any N-body-like inner loop
// (fields, potentials, contact sums), and the shape in which the guard was first found
// to cost a factor rather than a per-cent: with the throw INLINE, 1.63x on the whole
// loop, byte-identical output. Both binaries print a checksum so a reader can see the
// two computations agree to the bit.

#include "ga/ga_ega.hpp"

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <vector>

using namespace hd::ga;
using namespace hd::ga::ega;

int main()
{
    // deterministic point sets, no RNG: sources on a ring, targets on a lattice below
    constexpr std::size_t NS = 64, NT = 128;
    std::vector<vec3d> q(NS), m(NS), p(NT);
    for (std::size_t i = 0; i < NS; ++i) {
        double const a = 2.0 * pi * double(i) / double(NS);
        q[i] = vec3d(0.065 * std::cos(a), 0.065 * std::sin(a), 0.0);
        m[i] = vec3d(0.0, 0.0, (i % 2 == 0) ? 1.0 : -1.0);
    }
    for (std::size_t j = 0; j < NT; ++j) {
        p[j] = vec3d(-0.16 + 0.32 * double(j % 16) / 15.0,
                     -0.16 + 0.32 * double(j / 16) / 7.0, -0.018);
    }

    // one pass = NS x NT pairs = 8192; 30000 passes ~ 0.5-1 s per run, long enough for
    // turbo/thermal settling not to be the thing measured (400 passes was ~10 ms and read
    // 4.2 -> 2.5 ns across three rounds on the same binary)
    constexpr int REPS = 30000;
    vec3d acc(0.0, 0.0, 0.0);
    auto const t0 = std::chrono::steady_clock::now();
    for (int rep = 0; rep < REPS; ++rep) {
        for (std::size_t j = 0; j < NT; ++j) {
            vec3d f(0.0, 0.0, 0.0);
            for (std::size_t i = 0; i < NS; ++i) {
                vec3d const r = p[j] - q[i];
                double const rn = nrm(r);
                vec3d const rh = r / rn; // the checked division
                double const c = double(dot(m[i], rh)) / (rn * rn * rn);
                f = f + c * rh;
            }
            acc = acc + f;
        }
    }
    auto const t1 = std::chrono::steady_clock::now();
    double const ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
    double const pairs = double(REPS) * double(NS) * double(NT);

#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    char const* const arm = "guard ON ";
#else
    char const* const arm = "guard OFF";
#endif
    std::printf("%s  %.2f ns per pair   checksum %.12e %.12e %.12e\n", arm, ns / pairs,
                acc.x, acc.y, acc.z);
    return 0;
}
