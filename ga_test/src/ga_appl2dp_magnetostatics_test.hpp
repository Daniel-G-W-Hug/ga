#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <algorithm> // std::max
#include <cmath>     // std::abs, std::atan2, std::cos, std::hypot, std::sin
#include <cstddef>   // size_t
#include <vector>    // std::vector

#include "ga/ga_pga.hpp" // ga_usr_magnetostatics.hpp

/////////////////////////////////////////////////////////////////////////////////////////
// The 2D magnetostatic primitives of ga_usr_magnetostatics.hpp: the vector derivative
// as one equation for both laws, the flux function, the stress sandwich and the stress
// wrench. Every check is against an INDEPENDENT closed form written out classically
// below, and each case carries a falsifier.
/////////////////////////////////////////////////////////////////////////////////////////

namespace ms2d {

using namespace hd::ga;

// the classical field of a z-directed line current I at p: B = mu0 I/(2 pi rho) phi_hat
inline vec2d b_classical(vec2d const& r, vec2d const& p, double I)
{
    double const dx = r.x - p.x, dy = r.y - p.y, rho2 = dx * dx + dy * dy;
    double const k = mu0 * I / (2.0 * pi * rho2);
    return vec2d{-k * dy, k * dx}; // phi_hat * rho = (-dy, dx)
}

inline double dist(vec2d const& a, vec2d const& b) { return nrm(a - b); }

} // namespace ms2d

TEST_SUITE("magnetostatics 2d")
{
    using namespace ms2d;

    TEST_CASE("magnetostatics 2d: the field of a line current is the vector inverse")
    {
        line_current const c{vec2d{0.3, -0.2}, pscalar2d(250.0)};
        double worst = 0.0, worst_flip = 0.0;
        for (vec2d const r : {vec2d{1.0, 0.0}, vec2d{-0.4, 0.7}, vec2d{0.31, -0.19},
                              vec2d{2.5, 3.1}, vec2d{-1.0, -2.0}}) {
            vec2d const ref = b_classical(r, c.p, 250.0);
            double const scale = nrm(ref);
            worst = std::max(worst, nrm(b_line_current(r, c) - ref) / scale);
            // the falsifier: the opposite current direction is not the same field
            worst_flip =
                std::max(worst_flip,
                         nrm(b_line_current(r, {c.p, pscalar2d(-250.0)}) - ref) / scale);
        }
        CHECK(worst < 1e-14);
        CHECK(worst_flip > 1.0);
        // the form of the text: the right dual of the vector inverse, (mu0 I_z/(2 pi))
        // r_dual((r - p)^-1), and the potential -(mu0/(2 pi)) ln|r - p| I
        for (vec2d const r : {vec2d{1.0, 0.0}, vec2d{-0.4, 0.7}}) {
            vec2d const ref = (mu0 * 250.0 / (2.0 * pi)) * r_dual(inv(r - c.p));
            CHECK(nrm(b_line_current(r, c) - ref) < 1e-15 * nrm(ref));
            double const dx = r.x - c.p.x, dy = r.y - c.p.y;
            double const A = -(mu0 * 250.0 / (2.0 * pi)) * std::log(std::hypot(dx, dy));
            CHECK(std::abs(value_t(A_line_current(r, c)) - A) < 1e-15 * std::abs(A));
        }
    }

    TEST_CASE("magnetostatics 2d: the circulation of b around a wire is mu0 I")
    {
        // the normalisation of the Green's function: Ampere's law in integral form,
        // with the counter-clockwise tangent of the outward normal written out
        line_current const w{vec2d{0.05, -0.02}, pscalar2d(40.0)};
        auto circulation = [&](vec2d const& centre, double radius) {
            double sum = 0.0;
            for (auto const& q : circle_contour(centre, radius, 720)) {
                vec2d const b = b_line_current(q.P, w);
                sum += (-b.x * q.n.y + b.y * q.n.x) * q.ds;
            }
            return sum;
        };
        double const ref = mu0 * 40.0;
        CHECK(std::abs(circulation(vec2d{0.0, 0.0}, 0.2) - ref) < 1e-12 * ref);
        CHECK(std::abs(circulation(vec2d{0.1, 0.05}, 0.3) - ref) < 1e-12 * ref);
        // a contour that does not enclose the wire: no circulation
        CHECK(std::abs(circulation(vec2d{0.6, 0.4}, 0.2)) < 1e-12 * ref);
    }

    TEST_CASE("magnetostatics 2d: nabla b = mu0 J -- one equation for both laws")
    {
        // (a) inside a uniformly conducting round wire, b = (mu0/2) (r - c) J: the field
        // is
        //     linear, so the finite-difference derivative is exact to rounding, and its
        //     only part is the pseudoscalar mu0 J -- div b = 0 and curl b = mu0 J at once
        vec2d const c{0.1, 0.05};
        pscalar2d const J(3.0e6); // A/m^2
        auto inside = [&](vec2d const& r) { return 0.5 * mu0 * ((r - c) * J); };
        for (vec2d const r : {vec2d{0.12, 0.04}, vec2d{0.09, 0.07}}) {
            mvec2d const N = ega::nabla(inside, r);
            double const ref = mu0 * value_t(J);
            CHECK(std::abs(value_t(gr0(N))) < 1e-12 * ref);
            CHECK(nrm(gr1(N)) < 1e-12 * ref);
            CHECK(std::abs(value_t(gr2(N)) - ref) < 1e-9 * ref);
            // the falsifier: without the 1/2 the equation reads 2 mu0 J
            auto wrong = [&](vec2d const& q) { return mu0 * ((q - c) * J); };
            CHECK(std::abs(value_t(gr2(ega::nabla(wrong, r))) - ref) > 0.5 * ref);
        }

        // (b) outside any current the derivative vanishes in ALL grades; it is computed
        // by
        //     finite differences, so it is checked by convergence: the residual relative
        //     to the field's own scale |b|/rho must fall as h^2 (the default central
        //     scheme) -- a ratio of ~4 per halving of h
        line_current const w{vec2d{0.0, 0.0}, pscalar2d(100.0)};
        auto outside = [&](vec2d const& r) { return b_line_current(r, w); };
        vec2d const r{0.3, 0.4};
        double const scale = nrm(outside(r)) / nrm(r);
        double prev = 0.0;
        for (double h : {1e-2, 5e-3, 2.5e-3}) {
            mvec2d const N = ega::nabla(outside, r, h);
            double const res =
                (std::abs(value_t(gr0(N))) + std::abs(value_t(gr2(N)))) / scale;
            if (prev > 0.0) {
                CHECK(prev / res > 3.8);
                CHECK(prev / res < 4.2);
            }
            prev = res;
        }
    }

    TEST_CASE("magnetostatics 2d: b = -nabla A, and flux is a difference of A")
    {
        line_current const w{vec2d{0.05, -0.02}, pscalar2d(40.0)};
        auto A = [&](vec2d const& r) { return A_line_current(r, w); };

        // b = -nabla A at a point (central differences: the error falls 4x per halving)
        vec2d const r{0.4, 0.3};
        vec2d const b = b_line_current(r, w);
        double prev = 0.0;
        for (double h : {2e-2, 1e-2, 5e-3}) {
            double const err = nrm(gr1(-ega::nabla(A, r, h)) - b) / nrm(b);
            if (prev > 0.0) CHECK(prev / err > 3.8); // central differences: 2nd order
            prev = err;
        }

        // inside a uniformly conducting wire, A = -(mu0/4) |r - c|^2 J: its negative
        // derivative is the field (mu0/2) (r - c) J, and nabla^2 A = -mu0 J -- both
        // exact to rounding for a quadratic potential
        vec2d const cw{0.1, 0.05};
        pscalar2d const J(3.0e6);
        auto A_in = [&](vec2d const& q) {
            return pscalar2d(-0.25 * mu0 * value_t(J) * nrm_sq(q - cw));
        };
        vec2d const rin{0.12, 0.04};
        vec2d const b_in = 0.5 * mu0 * ((rin - cw) * J);
        CHECK(nrm(gr1(-ega::nabla(A_in, rin, 1e-3)) - b_in) < 1e-9 * nrm(b_in));
        double const lap = value_t(gr2(ega::laplacian(A_in, rin, 1e-3)));
        CHECK(std::abs(lap + mu0 * value_t(J)) < 1e-6 * mu0 * value_t(J));
        CHECK(std::abs(lap - mu0 * value_t(J)) > mu0 * value_t(J)); // the falsifier

        // the flux per length through the straight path P -> Q, by 20-point Gauss-
        // Legendre quadrature of b . nu with nu = l_dual(t), equals A(Q) - A(P)
        vec2d const P{0.5, -0.4}, Q{-0.2, 0.6};
        double const L = dist(P, Q);
        vec2d const t = (Q - P) / L;
        vec2d const nu = l_dual(t); // t turned by -90 deg
        // nodes and weights on [-1, 1], 20 points (symmetric pairs)
        double const xg[10] = {0.0765265211334973, 0.2277858511416451, 0.3737060887154195,
                               0.5108670019508271, 0.6360536807265150, 0.7463319064601508,
                               0.8391169718222188, 0.9122344282513259, 0.9639719272779138,
                               0.9931285991850949};
        double const wg[10] = {0.1527533871307258, 0.1491729864726037, 0.1420961093183820,
                               0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
                               0.0832767415767048, 0.0626720483341091, 0.0406014298003869,
                               0.0176140071391521};
        // composite rule: 8 panels of 20 points -- the path passes ~0.15 m from the wire,
        // where one 20-point panel is still 1e-3 off (measured: 40 points 1e-7, 80 2e-15)
        double flux = 0.0;
        int const panels = 8;
        double const Lp = L / panels;
        for (int j = 0; j < panels; ++j)
            for (int k = 0; k < 10; ++k)
                for (double s : {-xg[k], xg[k]}) {
                    vec2d const x = P + (Lp * (double(j) + 0.5 * (s + 1.0))) * t;
                    flux += wg[k] * 0.5 * Lp * value_t(dot(b_line_current(x, w), nu));
                }
        double const dA = value_t(A(Q)) - value_t(A(P));
        CHECK(std::abs(flux - dA) < 1e-12 * std::abs(dA));
        // the falsifier: the other normal, r_dual(t) (t turned by +90 deg), gives the
        // opposite flux
        CHECK(nrm(r_dual(t) + nu) < 1e-15);
        CHECK(std::abs(flux + dA) > 1.0 * std::abs(dA));
    }

    TEST_CASE("magnetostatics 2d: the stress sandwich is Maxwell's tensor")
    {
        double worst = 0.0, worst_flip = 0.0;
        for (int k = 0; k < 24; ++k) {
            vec2d const b{0.9 * std::cos(0.7 * k), 1.3 * std::sin(1.1 * k + 0.2)};
            vec2d const n{std::cos(0.37 * k), std::sin(0.37 * k)};
            // tensor form: ((b.n) b - (1/2) b^2 n) / mu0, written out
            double const bn = b.x * n.x + b.y * n.y, b2 = b.x * b.x + b.y * b.y;
            vec2d const ref{(bn * b.x - 0.5 * b2 * n.x) / mu0,
                            (bn * b.y - 0.5 * b2 * n.y) / mu0};
            double const s = std::max(nrm(ref), 1.0);
            worst = std::max(worst, nrm(traction(b, n) - ref) / s);
            worst_flip = std::max(worst_flip, nrm(-1.0 * traction(b, n) - ref) / s);
        }
        CHECK(worst < 1e-14);
        CHECK(worst_flip > 0.1);
    }

    TEST_CASE("magnetostatics 2d: the stress wrench equals the Lorentz wrench")
    {
        // two parallel line currents; the stress over a circle around wire 1 must give
        // the force on wire 1 ON ITS LINE OF ACTION through p1: W = wdg(p1, F1)
        line_current const w1{vec2d{0.10, 0.05}, pscalar2d(300.0)};
        line_current const w2{vec2d{0.45, -0.25}, pscalar2d(-180.0)};
        std::vector<line_current> const both{w1, w2};
        auto b = [&](vec2d const& r) { return b_line_currents(r, both); };

        vec2d const F1 = lorentz_force(b_line_current(w1.p, w2), w1.I);
        bivec2dp const W_ref = wdg(vec2dp(w1.p.x, w1.p.y, 1.0), vec2dp(F1.x, F1.y, 0.0));
        // the classical magnitude mu0 I1 I2 / (2 pi d), repelling for opposite currents
        double const d = dist(w1.p, w2.p);
        CHECK(std::abs(nrm(F1) - mu0 * 300.0 * 180.0 / (2.0 * pi * d)) < 1e-13 * nrm(F1));
        CHECK(dot(F1, w1.p - w2.p) > 0.0);
        // the form of the text, F = I_z r_dual(b), and the classical I e_z x b
        vec2d const b12 = b_line_current(w1.p, w2);
        CHECK(nrm(F1 - 300.0 * r_dual(b12)) < 1e-15 * nrm(F1));
        CHECK(nrm(F1 - (w1.I >> b12)) < 1e-15 * nrm(F1)); // dual(I) >> b, as in 3D
        CHECK(nrm(F1 - vec2d{-300.0 * b12.y, 300.0 * b12.x}) < 1e-15 * nrm(F1));

        auto const scale = std::max(nrm(F1), 1e-30);
        auto rel = [&](bivec2dp const& W) {
            return std::max({std::abs(W.x - W_ref.x), std::abs(W.y - W_ref.y),
                             std::abs(W.z - W_ref.z) / 0.1}) /
                   scale;
        };
        // two different contours around wire 1 (radius and centre): the same wrench
        auto const cA = circle_contour(w1.p, 0.1, 720);
        auto const cB = circle_contour(vec2d{0.12, 0.02}, 0.2, 720);
        CHECK(rel(stress_wrench(b, cA)) < 1e-10);
        CHECK(rel(stress_wrench(b, cB)) < 1e-10);

        // action = reaction as ONE bivector identity: the wrench on wire 2 (a contour
        // around it) is the negative -- the forces are central, so the lines coincide
        auto const c2 = circle_contour(w2.p, 0.1, 720);
        bivec2dp const sum = stress_wrench(b, cA) + stress_wrench(b, c2);
        CHECK(std::max({std::abs(sum.x), std::abs(sum.y), std::abs(sum.z) / 0.1}) <
              1e-10 * scale);

        // the falsifier: the wedge in the wrong order (force first) is the wrong line
        bivec2dp Wbad{};
        for (auto const& q : cA) {
            vec2d const t = traction(b(q.P), q.n);
            Wbad += wdg(vec2dp(t.x, t.y, 0.0), vec2dp(q.P.x, q.P.y, 1.0)) * q.ds;
        }
        CHECK(rel(Wbad) > 1.0);
    }

    TEST_CASE("magnetostatics 2d: a magnetised disc as surface currents K = m ^ n")
    {
        // a uniformly magnetised disc: inside B = (mu0/2) m, outside the 2D dipole field,
        // (mu0 R^2 / (2 r^2)) (2 (m.r_hat) r_hat - m) -- written out classically
        vec2d const c{0.0, 0.0}, m{6.0e5, 3.5e5}; // A/m
        double const R = 0.02;
        // the surface current of the text, K = m ^ n, written out
        vec2d const nn{0.6, 0.8};
        CHECK(std::abs(value_t(surface_current(m, nn)) - (m.x * nn.y - m.y * nn.x)) <
              1e-15 * nrm(m));
        auto const cs = magnetised_disc_currents(c, R, m, 2048);
        auto bfield = [&](vec2d const& r) { return b_line_currents(r, cs); };

        for (vec2d const r : {vec2d{0.003, -0.004}, vec2d{-0.008, 0.002}}) {
            vec2d const ref = 0.5 * mu0 * m;
            CHECK(nrm(bfield(r) - ref) < 1e-9 * nrm(ref));
        }
        for (vec2d const r : {vec2d{0.05, 0.01}, vec2d{-0.03, -0.04}}) {
            double const rr = nrm(r);
            vec2d const rh = r / rr;
            double const mr = m.x * rh.x + m.y * rh.y;
            vec2d const ref = (mu0 * R * R / (2.0 * rr * rr)) *
                              vec2d{2.0 * mr * rh.x - m.x, 2.0 * mr * rh.y - m.y};
            CHECK(nrm(bfield(r) - ref) < 1e-9 * nrm(ref));
            // the line dipole of the text, moment pi R^2 m: one reflection, no identity
            // term, (mu0 / (2 pi r^2)) rh mu rh
            vec2d const mu = (pi * R * R) * m;
            vec2d const sw = (mu0 / (2.0 * pi * rr * rr)) * (rh * mu * rh);
            CHECK(nrm(sw - ref) < 1e-14 * nrm(ref));
            // the dipole field in d dimensions at d = 2: the identity term drops out
            // with d - 2; S_2 = 2 pi
            vec2d const bd =
                (mu0 / (2.0 * 2.0 * pi * rr * rr)) * (0.0 * mu + 2.0 * (rh * mu * rh));
            CHECK(nrm(bd - ref) < 1e-14 * nrm(ref));
        }
        // the falsifier: K = n ^ m (the opposite orientation) reverses the field
        std::vector<line_current> flipped;
        for (auto const& l : cs)
            flipped.push_back({l.p, -l.I});
        vec2d const r_in{0.003, -0.004};
        CHECK(nrm(b_line_currents(r_in, flipped) - 0.5 * mu0 * m) > nrm(0.5 * mu0 * m));
    }
}
