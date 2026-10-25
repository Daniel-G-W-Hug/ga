#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <algorithm> // std::max
#include <cmath>     // std::abs, std::pow, std::sqrt
#include <cstdio>    // std::printf
#include <vector>    // std::vector

#include "ga/ga_pga.hpp" // ga_usr_magnetostatics.hpp

/////////////////////////////////////////////////////////////////////////////////////////
// The 3D magnetostatic primitives of ga_usr_magnetostatics.hpp, the 2D ones lifted.
// Every check is against an INDEPENDENT form written out classically below -- cross
// products, Biot-Savart by quadrature, the dipole-dipole force of Yung et al. -- and
// each case carries a falsifier.
/////////////////////////////////////////////////////////////////////////////////////////

namespace ms3d {

using namespace hd::ga;

inline vec3d cross3(vec3d const& a, vec3d const& b)
{
    return vec3d{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline double dot3(vec3d const& a, vec3d const& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Biot-Savart by Gauss-Legendre quadrature along the segment: mu0 I/(4 pi) int dl x R/R^3
inline vec3d b_biot_savart(vec3d const& r, current_segment const& s, int panels)
{
    auto const& g = hd::ga::detail::gauss_legendre12();
    vec3d const L = s.b - s.a;
    vec3d b{};
    for (int p = 0; p < panels; ++p)
        for (size_t i = 0; i < 12; ++i) {
            double const t = (p + 0.5 * (g.x[i] + 1.0)) / panels;
            vec3d const R = r - (s.a + t * L);
            double const R3 = std::pow(dot3(R, R), 1.5);
            b += (0.5 * g.w[i] / panels / R3) * cross3(L, R);
        }
    return (mu0 * s.I / (4.0 * pi)) * b;
}

// the force on dipole j in the field of dipole i (Yung, Landecker, Villani 1998),
// r = p_j - p_i -- written out
inline vec3d F_dipoles(vec3d const& r, vec3d const& mi, vec3d const& mj)
{
    double const rn = std::sqrt(dot3(r, r));
    vec3d const h = r / rn;
    double const c = 3.0 * mu0 / (4.0 * pi * std::pow(rn, 4));
    return c * (dot3(mi, mj) * h + dot3(h, mj) * mi + dot3(h, mi) * mj -
                5.0 * dot3(h, mi) * dot3(h, mj) * h);
}

inline double rel(vec3d const& a, vec3d const& ref) { return nrm(a - ref) / nrm(ref); }

} // namespace ms3d

TEST_SUITE("magnetostatics 3d")
{
    using namespace ms3d;

    TEST_CASE("magnetostatics 3d: a segment's field is Biot-Savart, and its long limit "
              "the 2D kernel")
    {
        current_segment const s{vec3d{0.1, -0.2, 0.05}, vec3d{0.4, 0.3, -0.25}, 120.0};
        double worst = 0.0;
        for (vec3d const r : {vec3d{0.5, 0.0, 0.0}, vec3d{-0.3, 0.4, 0.2},
                              vec3d{0.25, 0.09, -0.1}, vec3d{1.0, -1.0, 0.5}})
            worst = std::max(worst, rel(b_segment(r, s), b_biot_savart(r, s, 64)));
        std::printf("3D segment against Biot-Savart: %.1e\n", worst);
        CHECK(worst < 1e-12);
        // the 2D picture: a segment of +-1e4 m along e3 seen at z = 0 is the line current
        line_current const c2{vec2d{0.3, -0.2}, pscalar2d(250.0)};
        current_segment const s3{vec3d{0.3, -0.2, -1e4}, vec3d{0.3, -0.2, 1e4}, 250.0};
        for (vec2d const r : {vec2d{1.0, 0.0}, vec2d{-0.4, 0.7}}) {
            vec2d const b2 = b_line_current(r, c2);
            vec3d const b3 = b_segment(vec3d{r.x, r.y, 0.0}, s3);
            CHECK(rel(b3, vec3d{b2.x, b2.y, 0.0}) < 1e-8); // the ends' (d / L)^2
            CHECK(std::abs(b3.z) < 1e-14 * nrm(b3));
        }
        // the falsifier: the segment reversed is the opposite field
        current_segment const back{s.b, s.a, s.I};
        CHECK(rel(b_segment(vec3d{0.5, 0.0, 0.0}, back),
                  b_biot_savart(vec3d{0.5, 0.0, 0.0}, s, 64)) > 1.0);
    }

    TEST_CASE("magnetostatics 3d: b = curl A for a segment")
    {
        current_segment const s{vec3d{0.1, -0.2, 0.05}, vec3d{0.4, 0.3, -0.25}, 120.0};
        auto A = [&](vec3d const& r) { return A_segment(r, s); };
        // a central difference: the error falls as h^2 (measured 2026-10-25: ratio 100
        // from h 1e-3 to 1e-4)
        double worst[2] = {0.0, 0.0};
        for (int k = 0; k < 2; ++k)
            for (vec3d const r : {vec3d{0.5, 0.0, 0.0}, vec3d{-0.3, 0.4, 0.2}}) {
                // nabla ^ A = I curl A
                vec3d const curl = cmpl(ega::nabla_wdg(A, r, k == 0 ? 1e-3 : 1e-4));
                worst[k] = std::max(worst[k], rel(curl, b_segment(r, s)));
            }
        std::printf("3D segment: curl A against b %.1e (h 1e-3), %.1e (h 1e-4)\n",
                    worst[0], worst[1]);
        CHECK(worst[1] < 1e-7);
        CHECK(std::abs(worst[0] / worst[1] - 100.0) < 10.0);
    }

    TEST_CASE("magnetostatics 3d: the stress sandwich is Maxwell's tensor, and the "
              "bivector form's sign")
    {
        double worst = 0.0, worst_biv = 0.0, worst_flip = 0.0;
        for (int k = 0; k < 24; ++k) {
            vec3d const b{0.9 * std::cos(0.7 * k), 1.3 * std::sin(1.1 * k + 0.2),
                          0.4 * std::cos(0.3 * k + 1.0)};
            vec3d n{std::cos(0.37 * k), std::sin(0.37 * k) * std::cos(0.5 * k),
                    std::sin(0.37 * k) * std::sin(0.5 * k)};
            n = n / nrm(n);
            // tensor form ((b.n) b - (1/2) b^2 n) / mu0, written out
            vec3d const ref = (dot3(b, n) / mu0) * b - (0.5 * dot3(b, b) / mu0) * n;
            double const s = std::max(nrm(ref), 1.0);
            worst = std::max(worst, nrm(traction(b, n) - ref) / s);
            // the bivector convention B = I3 b: the same traction is -(1/(2 mu0)) B n B
            bivec3d const B = cmpl(b);
            vec3d const tB = (-1.0 / (2.0 * mu0)) * gr1(B * n * B);
            worst_biv = std::max(worst_biv, nrm(tB - ref) / s);
            worst_flip = std::max(worst_flip, nrm(-1.0 * tB - ref) / s);
        }
        CHECK(worst < 1e-14);
        CHECK(worst_biv < 1e-14);
        CHECK(worst_flip > 0.1);
    }

    TEST_CASE("magnetostatics 3d: the Lorentz force is I dl x b")
    {
        vec3d const b{0.3, -1.1, 0.7}, Idl{2.0, 0.5, -1.5};
        CHECK(nrm(lorentz_force(b, Idl) - cross3(Idl, b)) < 1e-15);
        CHECK(nrm(lorentz_force(b, Idl) + cross3(Idl, b)) > 1.0); // the falsifier
    }

    TEST_CASE(
        "3D: a loop's far field is the dipole sandwich, the octupole falling as r^-2")
    {
        // a loop of radius R, 256 chords, moment I (N/2) R^2 sin(2 pi / N) n; the
        // central symmetry kills the quadrupole, so the error falls as (R / r)^2
        double const R = 0.01, I = 40.0;
        size_t const N = 256;
        vec3d const n = vec3d{0.2, -0.3, 0.9} / nrm(vec3d{0.2, -0.3, 0.9});
        auto const loop = circle_loop(vec3d{}, n, R, N, I);
        vec3d const mu = (I * 0.5 * N * R * R * std::sin(2.0 * pi / N)) * n;
        vec3d const dir = vec3d{0.6, 0.3, -0.5} / nrm(vec3d{0.6, 0.3, -0.5});
        double err[2];
        for (int k = 0; k < 2; ++k) {
            vec3d const r = (k == 0 ? 50.0 * R : 500.0 * R) * dir;
            err[k] = rel(b_segments(r, loop), b_dipole(r, mu));
            // the dipole sandwich against its vector-calculus form
            double const rn = nrm(r);
            vec3d const rh = r / rn;
            vec3d const ref =
                (mu0 / (4.0 * pi * rn * rn * rn)) * (3.0 * dot3(mu, rh) * rh - mu);
            CHECK(rel(b_dipole(r, mu), ref) < 1e-14);
            // the bivector form is the same expression: no sign change
            bivec3d const m = cmpl(mu);
            bivec3d const Bga =
                (mu0 / (8.0 * pi * rn * rn * rn)) * (m + 3.0 * gr2(rh * m * rh));
            CHECK(nrm(cmpl(Bga) - b_dipole(r, mu)) < 1e-14 * nrm(ref));
        }
        std::printf(
            "3D loop against the dipole: %.2e at 50 R, %.2e at 500 R (ratio %.1f)\n",
            err[0], err[1], err[0] / err[1]);
        CHECK(std::abs(err[0] / err[1] - 100.0) < 5.0);
        // the falsifier: the opposite moment
        CHECK(rel(b_segments(50.0 * R * dir, loop), b_dipole(50.0 * R * dir, -1.0 * mu)) >
              1.0);
    }

    TEST_CASE("magnetostatics 3d: the stress wrench of two dipoles is the force line "
              "plus the couple")
    {
        // two point dipoles; the stress over a sphere around dipole j in the total field
        // gives wdg(P_j, F) + couple(mu_j x b_i(P_j)), F from Yung et al.; a SCREW, the
        // couple coming out of the sum of force lines with no couple term written
        vec3d const pi_{0.0, 0.0, 0.0}, pj{0.06, -0.02, 0.03};
        vec3d const mi{0.0, 0.0, 2.0}, mj{0.5, 1.2, -0.4}; // A m^2
        auto b = [&](vec3d const& r) {
            return b_dipole(r - pi_, mi) + b_dipole(r - pj, mj);
        };
        vec3d const F = F_dipoles(pj - pi_, mi, mj);
        vec3d const T = cross3(mj, b_dipole(pj - pi_, mi));
        bivec3dp const W_ref =
            wdg(vec3dp(pj.x, pj.y, pj.z, 1.0), vec3dp(F.x, F.y, F.z, 0.0)) +
            bivec3dp{0.0, 0.0, 0.0, T.x, T.y, T.z};
        auto const S = sphere_surface(pj, 0.02, 8, 96);
        bivec3dp const W = stress_wrench(b, S);
        auto err = [&](bivec3dp const& X, bivec3dp const& Y) {
            double const f = nrm(F), t = std::max(nrm(T), 1e-30);
            return std::max({std::abs(X.vx - Y.vx) / f, std::abs(X.vy - Y.vy) / f,
                             std::abs(X.vz - Y.vz) / f, std::abs(X.mx - Y.mx) / t,
                             std::abs(X.my - Y.my) / t, std::abs(X.mz - Y.mz) / t});
        };
        std::printf(
            "3D two dipoles: the stress wrench against force line + couple %.1e\n",
            err(W, W_ref));
        CHECK(err(W, W_ref) < 1e-8);
        // action = reaction as ONE bivector identity: the wrench on dipole i is the
        // negative, though the forces are not central and the couples differ
        bivec3dp const sum = W + stress_wrench(b, sphere_surface(pi_, 0.02, 8, 96));
        CHECK(err(sum, bivec3dp{}) < 1e-8);
        // the falsifier: the wedge in the wrong order
        bivec3dp Wbad{};
        for (auto const& q : S) {
            vec3d const t = traction(b(q.P), q.n);
            Wbad +=
                wdg(vec3dp(t.x, t.y, t.z, 0.0), vec3dp(q.P.x, q.P.y, q.P.z, 1.0)) * q.dA;
        }
        CHECK(err(Wbad, W_ref) > 1.0);
    }

    TEST_CASE("magnetostatics 3d: a magnetised sphere as surface currents K = m x n")
    {
        // a uniformly magnetised sphere: inside b = (2/3) mu0 m, outside the dipole of
        // moment (4/3) pi R^3 m -- the bridge from a magnetised body to the point dipole
        vec3d const c{}, m{3.0e5, -2.0e5, 6.0e5}; // A/m
        double const R = 0.02;
        CHECK(nrm(surface_current(m, vec3d{0, 0, 1}) - cross3(m, vec3d{0, 0, 1})) < 1e-9);
        auto const cs = magnetised_sphere_currents(c, R, m, 4, 512);
        auto bf = [&](vec3d const& r) { return b_segments(r, cs); };
        double in = 0.0, out = 0.0;
        for (vec3d const r : {vec3d{0.003, -0.004, 0.002}, vec3d{-0.006, 0.002, -0.005}})
            in = std::max(in, rel(bf(r), (2.0 / 3.0) * mu0 * m));
        vec3d const mu = (4.0 / 3.0 * pi * R * R * R) * m;
        for (vec3d const r : {vec3d{0.05, 0.01, -0.02}, vec3d{-0.03, -0.04, 0.03}})
            out = std::max(out, rel(bf(r), b_dipole(r, mu)));
        std::printf("3D magnetised sphere: inside %.1e, outside %.1e\n", in, out);
        // measured 2026-10-25: 5.0e-6 inside, 2.8e-5 outside -- the rings' chords (512
        // per ring) and 48 Gauss-Legendre rings
        CHECK(in < 2e-5);
        CHECK(out < 1e-4);
        // the falsifier: K = n x m reverses the field
        std::vector<current_segment> flipped;
        for (auto const& s : cs)
            flipped.push_back({s.a, s.b, -s.I});
        CHECK(rel(b_segments(vec3d{0.003, -0.004, 0.002}, flipped),
                  (2.0 / 3.0) * mu0 * m) > 1.0);
    }
}
