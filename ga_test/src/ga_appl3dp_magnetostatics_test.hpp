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
                // nabla ^ A = dual(curl A), and the dual twice is the identity in 3D
                vec3d const curl = dual(ega::nabla_wdg(A, r, k == 0 ? 1e-3 : 1e-4));
                worst[k] = std::max(worst[k], rel(curl, b_segment(r, s)));
            }
        std::printf("3D segment: curl A against b %.1e (h 1e-3), %.1e (h 1e-4)\n",
                    worst[0], worst[1]);
        CHECK(worst[1] < 1e-7);
        CHECK(std::abs(worst[0] / worst[1] - 100.0) < 10.0);
    }

    TEST_CASE("magnetostatics 3d: nabla b = mu0 dual(J) -- one equation for both laws")
    {
        // (a) inside a uniformly conducting round wire through c along J, the field
        // b = (mu0/2) J x (r - c) is linear, so the central difference is exact to
        // rounding, and nabla b has ONE part, the bivector mu0 dual(J): div b = 0 and
        // curl b = mu0 J at once
        vec3d const c{0.1, 0.05, -0.02};
        vec3d const J = 3.0e6 * (vec3d{0.2, -0.4, 0.9} / nrm(vec3d{0.2, -0.4, 0.9}));
        auto inside = [&](vec3d const& r) { return 0.5 * mu0 * cross3(J, r - c); };
        bivec3d const ref = mu0 * dual(J);
        double const s = nrm(ref);
        for (vec3d const r : {vec3d{0.12, 0.04, 0.0}, vec3d{0.09, 0.07, -0.05}}) {
            mvec3d const N = ega::nabla(inside, r, 1e-3);
            CHECK(std::abs(value_t(gr0(N))) < 1e-9 * s);
            CHECK(nrm(gr1(N)) < 1e-9 * s);
            CHECK(nrm(gr2(N) - ref) < 1e-9 * s);
            CHECK(std::abs(value_t(gr3(N))) < 1e-9 * s);
            // the falsifier: the opposite current gives the opposite bivector
            auto back = [&](vec3d const& q) { return 0.5 * mu0 * cross3(q - c, J); };
            CHECK(nrm(gr2(ega::nabla(back, r, 1e-3)) - ref) > s);
        }
        // the 2D equation is the 3D one in the plane: dual(J_z e3) = J_z e12
        CHECK(dual(vec3d{0.0, 0.0, 2.5}) == bivec3d{0.0, 0.0, 2.5});

        // (b) outside the currents of a closed loop the derivative vanishes in ALL
        // grades; by finite differences, so by convergence: the residual relative to
        // |b|/|r| falls as h^2, a ratio of ~4 per halving of h
        auto const loop = circle_loop(vec3d{}, vec3d{0.0, 0.0, 1.0}, 0.05, 64, 10.0);
        auto outside = [&](vec3d const& r) { return b_segments(r, loop); };
        vec3d const r{0.08, 0.03, 0.04};
        double const scale = nrm(outside(r)) / nrm(r);
        double prev = 0.0;
        for (double h : {4e-3, 2e-3, 1e-3}) {
            mvec3d const N = ega::nabla(outside, r, h);
            double const res = (std::abs(value_t(gr0(N))) + nrm(gr2(N))) / scale;
            if (prev > 0.0) {
                CHECK(prev / res > 3.8);
                CHECK(prev / res < 4.2);
            }
            prev = res;
        }
    }

    TEST_CASE("magnetostatics 3d: b = -nabla dual(A) in the Coulomb gauge")
    {
        // (a) inside a uniformly conducting wire through c along J, the potential
        // A = -(mu0/4) |r_perp|^2 J is quadratic and divergence-free: -nabla dual(A) is
        // the field (mu0/2) J x (r - c), nabla^2 dual(A) = -mu0 dual(J), both exact to
        // rounding by central differences
        vec3d const c{0.1, 0.05, -0.02};
        vec3d const J = 3.0e6 * (vec3d{0.2, -0.4, 0.9} / nrm(vec3d{0.2, -0.4, 0.9}));
        vec3d const u = J / nrm(J);
        auto A_in = [&](vec3d const& r) {
            vec3d const d = r - c;
            vec3d const perp = d - dot3(d, u) * u;
            return dual((-0.25 * mu0 * dot3(perp, perp)) * J);
        };
        vec3d const r{0.12, 0.04, 0.0};
        vec3d const b_in = 0.5 * mu0 * cross3(J, r - c);
        mvec3d const N = -1.0 * ega::nabla(A_in, r, 1e-3);
        CHECK(nrm(gr1(N) - b_in) < 1e-9 * nrm(b_in));
        CHECK(std::abs(value_t(gr3(N))) < 1e-9 * nrm(b_in)); // div A = 0
        bivec3d const lap = gr2(ega::laplacian(A_in, r, 1e-3));
        CHECK(nrm(lap + mu0 * dual(J)) < 1e-6 * mu0 * nrm(J));
        CHECK(nrm(lap - mu0 * dual(J)) > mu0 * nrm(J)); // the falsifier: the sign

        // (b) the potential of a closed loop of chords (the Coulomb gauge holds for a
        // closed current): -nabla dual(A) converges to its field as h^2, and the
        // trivector part, the divergence, falls with it
        auto const loop = circle_loop(vec3d{}, vec3d{0.0, 0.3, 1.0}, 0.05, 64, 10.0);
        auto A_loop = [&](vec3d const& q) {
            vec3d a{};
            for (auto const& sg : loop)
                a += A_segment(q, sg);
            return dual(a);
        };
        vec3d const rl{0.04, 0.03, 0.05};
        vec3d const bl = b_segments(rl, loop);
        double err[2], div[2];
        for (int k = 0; k < 2; ++k) {
            mvec3d const M = -1.0 * ega::nabla(A_loop, rl, k == 0 ? 2e-3 : 1e-3);
            err[k] = rel(gr1(M), bl);
            div[k] = std::abs(value_t(gr3(M))) / nrm(bl);
        }
        std::printf("3D loop: -nabla dual(A) against b %.1e (h 2e-3), %.1e (h 1e-3); "
                    "divergence %.1e\n",
                    err[0], err[1], div[1]);
        CHECK(std::abs(err[0] / err[1] - 4.0) < 0.2);
        CHECK(div[1] < err[1]);
    }

    TEST_CASE("magnetostatics 3d: the dual identities of the bivector convention")
    {
        // the identities the text uses to move between b and B = dual(b), each against
        // its classical form written out
        double worst = 0.0, worst_flip = 0.0;
        for (int k = 0; k < 24; ++k) {
            vec3d const a{0.9 * std::cos(0.7 * k), 1.3 * std::sin(1.1 * k + 0.2),
                          0.4 * std::cos(0.3 * k + 1.0)};
            vec3d const b{-0.5 * std::sin(0.9 * k + 0.4), 0.8 * std::cos(0.6 * k),
                          1.1 * std::sin(0.2 * k + 0.7)};
            vec3d n{std::cos(0.37 * k), std::sin(0.37 * k) * std::cos(0.5 * k),
                    std::sin(0.37 * k) * std::sin(0.5 * k)};
            n = n / nrm(n);
            // the cross product is the dual of the wedge
            worst = std::max(worst, nrm(dual(wdg(a, b)) - cross3(a, b)));
            // the scalar product carries over the dual; the scalar part of the
            // geometric product of two bivectors is its negative
            worst =
                std::max(worst, std::abs(value_t(dot(dual(a), dual(b))) - dot3(a, b)));
            worst =
                std::max(worst, std::abs(value_t(gr0(dual(a) * dual(b))) + dot3(a, b)));
            // the commutator of two duals is the dual of the cross product
            worst = std::max(worst, nrm(cmt(dual(b), dual(a)) - dual(cross3(a, b))));
            // the field meeting itself changes sign: dual(b) n dual(b) = -(b n b), and
            // b n b has no trivector part
            worst = std::max(worst, nrm(gr1(dual(b) * n * dual(b)) + gr1(b * n * b)));
            worst = std::max(worst, std::abs(value_t(gr3(b * n * b))));
            worst_flip =
                std::max(worst_flip, nrm(gr1(dual(b) * n * dual(b)) - gr1(b * n * b)));
            // the dual of a reflection sandwich is the sandwich of the dual
            worst = std::max(worst, nrm(dual(gr1(n * a * n)) - gr2(n * dual(a) * n)));
        }
        CHECK(worst < 1e-14);
        CHECK(worst_flip > 0.1); // the falsifier: without the sign change
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
            // the bivector convention B = dual(b): the same traction is
            // -(1/(2 mu0)) B n B
            bivec3d const B = dual(b);
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
        // dF = dual(I dl ^ b) = dual(I dl) >> b, the contraction the same in 2D
        CHECK(nrm(lorentz_force(b, Idl) - dual(wdg(Idl, b))) < 1e-15);
        CHECK(nrm(lorentz_force(b, Idl) - (dual(Idl) >> b)) < 1e-15);
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
            // the dipole field in d dimensions at d = 3:
            // mu0 / (2 S_d r^d) ((d - 2) mu + d rh mu rh), S_3 = 4 pi
            vec3d const bd = (mu0 / (2.0 * 4.0 * pi * rn * rn * rn)) *
                             (1.0 * mu + 3.0 * gr1(rh * mu * rh));
            CHECK(rel(bd, ref) < 1e-14);
            // the bivector form is the same expression: no sign change
            bivec3d const m = dual(mu);
            bivec3d const Bga =
                (mu0 / (8.0 * pi * rn * rn * rn)) * (m + 3.0 * gr2(rh * m * rh));
            CHECK(nrm(Bga - dual(b_dipole(r, mu))) < 1e-14 * nrm(ref));
        }
        std::printf(
            "3D loop against the dipole: %.2e at 50 R, %.2e at 500 R (ratio %.1f)\n",
            err[0], err[1], err[0] / err[1]);
        CHECK(std::abs(err[0] / err[1] - 100.0) < 5.0);
        // the falsifier: the opposite moment
        CHECK(rel(b_segments(50.0 * R * dir, loop), b_dipole(50.0 * R * dir, -1.0 * mu)) >
              1.0);
    }

    TEST_CASE("magnetostatics 3d: a dipole's energy, torque and force from one field")
    {
        // U = -m . B, T = cmt(B, m), F = grad_r (m . B) in the bivector convention,
        // against the classical -mu . b, mu x b and the force of Yung et al.
        vec3d const mi{0.0, 0.0, 2.0}, mj{0.5, 1.2, -0.4}; // A m^2
        vec3d const r{0.06, -0.02, 0.03};                  // m_j relative to m_i
        bivec3d const Mi = dual(mi), Mj = dual(mj);
        bivec3d const B = b_dipole(r, Mi);
        vec3d const b = b_dipole(r, mi);
        double const s = nrm(mj) * nrm(b);
        double const coupling = dipole_coupling(r, Mi, Mj);
        CHECK(std::abs(coupling - dot3(mj, b)) < 1e-14 * s);          // U = -mu . b
        CHECK(std::abs(value_t(gr0(Mj * B)) + coupling) < 1e-14 * s); // <m B>_0 = -m . B
        // the torque is the commutator, the dual of mu x b
        CHECK(nrm(cmt(B, Mj) - dual(cross3(mj, b))) < 1e-14 * s);
        CHECK(nrm(cmt(Mj, B) - dual(cross3(mj, b))) > s); // the falsifier: the order
        // the force is the gradient of the coupling (central differences, 2nd order)
        auto grad = [&](double h) {
            vec3d g{};
            for (vec3d const e : {vec3d{1, 0, 0}, vec3d{0, 1, 0}, vec3d{0, 0, 1}})
                g += ((dipole_coupling(r + h * e, Mi, Mj) -
                       dipole_coupling(r - h * e, Mi, Mj)) /
                      (2.0 * h)) *
                     e;
            return g;
        };
        vec3d const F = F_dipoles(r, mi, mj);
        double const e1 = rel(grad(2e-4), F), e2 = rel(grad(1e-4), F);
        std::printf("3D dipole force as grad(m . B): %.1e (h 2e-4), %.1e (h 1e-4)\n", e1,
                    e2);
        CHECK(e2 < 1e-5);
        CHECK(std::abs(e1 / e2 - 4.0) < 0.2);
        CHECK(rel(-1.0 * grad(1e-4), F) > 1.0); // the falsifier: F = -grad(m . B)
        // the anisotropy read off one reflection: axial / equatorial = -2
        vec3d const ax = b_dipole(0.1 * (mi / nrm(mi)), mi);
        vec3d const eq = b_dipole(vec3d{0.1, 0.0, 0.0}, mi);
        CHECK(nrm(ax + 2.0 * eq) < 1e-14 * nrm(ax));
    }

    TEST_CASE("magnetostatics 3d: a loop of N chords converges to the circle as N^-2")
    {
        // on the axis of a circular loop b = mu0 I R^2 / (2 (R^2 + z^2)^(3/2)) e_z
        double const R = 0.05, I = 10.0;
        double err[2] = {0.0, 0.0};
        for (int k = 0; k < 2; ++k) {
            auto const loop =
                circle_loop(vec3d{}, vec3d{0.0, 0.0, 1.0}, R, k == 0 ? 32 : 64, I);
            for (double z : {0.0, 0.03}) {
                vec3d const ref{0.0, 0.0,
                                mu0 * I * R * R / (2.0 * std::pow(R * R + z * z, 1.5))};
                err[k] = std::max(err[k], rel(b_segments(vec3d{0.0, 0.0, z}, loop), ref));
            }
        }
        std::printf("3D loop of chords against the circle: %.2e (N 32), %.2e (N 64)\n",
                    err[0], err[1]);
        CHECK(std::abs(err[0] / err[1] - 4.0) < 0.2);
    }

    TEST_CASE("magnetostatics 3d: the stress wrench of two dipoles is the force line "
              "plus the couple")
    {
        // two point dipoles; the stress over a sphere around dipole j in the total field
        // gives wdg(P_j, F) + couple(mu_j x b_i(P_j)), F from Yung et al.; a SCREW, the
        // couple coming out of the sum of force lines with no couple term written
        vec3d const pi_{}, pj{0.06, -0.02, 0.03};
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
        // the weight carries the force, the bulk the torque about the reference point
        CHECK(std::abs(std::sqrt(value_t(weight_nrm_sq(W))) - nrm(F)) < 1e-8 * nrm(F));
        // a change of reference is one motor sandwich: about Q the torque is the
        // classical T + (p_j - Q) x F, and the force is unchanged
        vec3d const Q{-0.03, 0.05, 0.01};
        bivec3dp const Wq = move3dp(W, rrev(get_motor(vec3dp(Q.x, Q.y, Q.z, 0.0))));
        vec3d const Tq = T + cross3(pj - Q, F);
        CHECK(err(Wq, bivec3dp{F.x, F.y, F.z, Tq.x, Tq.y, Tq.z}) < 1e-8);
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
        // the form of the text, the same in 2D: dual(K) = m ^ n (the dual twice is the
        // identity in 3D)
        vec3d const nn = vec3d{0.2, 0.9, -0.4} / nrm(vec3d{0.2, 0.9, -0.4});
        CHECK(nrm(dual(surface_current(m, nn)) - wdg(m, nn)) < 1e-15 * nrm(m));
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
