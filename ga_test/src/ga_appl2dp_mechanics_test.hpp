// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <cmath>
#include <iostream>
#include <stdexcept> // std::runtime_error
#include <string>
using namespace std::string_literals; // enable s-suffix for std::string literals

#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

// hd::ga::det / lu_decomp / lu_backsubs come transitively via ga_pga.hpp.

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)


/////////////////////////////////////////////////////////////////////////////////////////
// PGA2DP physics tests preparation - Inertia matrix for rigid body dynamics
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA2DP: physics tests prep")
{

    TEST_CASE("pga2dp: get_point_inertia - single point mass")
    {
        fmt::println("pga2dp: get_point_inertia - single point mass");

        // Point mass at (1, 2) with mass 3.0
        double m = 3.0;
        Vec2dp<double> X{1.0, 2.0, 1.0}; // unitized point

        auto I = get_point_inertia(m, X);
        auto v = I.view();

        // Expected matrix from tex eq. 539-546:
        // m * [  0     Xz^2   -Xy*Xz  ]     [  0    1   -2 ]
        //     [ -Xz^2  0       Xx*Xz  ]  =  [ -1    0    1 ] * 3
        //     [ -Xx*Xz -Xy*Xz  Xx^2+Xy^2]   [ -1   -2    5 ]

        CHECK(v[0, 0] == doctest::Approx(0.0));
        CHECK(v[0, 1] == doctest::Approx(3.0));  // m * Xz^2 = 3 * 1 = 3
        CHECK(v[0, 2] == doctest::Approx(-6.0)); // m * (-Xy*Xz) = 3 * (-2) = -6

        CHECK(v[1, 0] == doctest::Approx(-3.0)); // m * (-Xz^2) = 3 * (-1) = -3
        CHECK(v[1, 1] == doctest::Approx(0.0));
        CHECK(v[1, 2] == doctest::Approx(3.0)); // m * Xx*Xz = 3 * 1 = 3

        CHECK(v[2, 0] == doctest::Approx(-3.0)); // m * (-Xx*Xz) = 3 * (-1) = -3
        CHECK(v[2, 1] == doctest::Approx(-6.0)); // m * (-Xy*Xz) = 3 * (-2) = -6
        CHECK(v[2, 2] == doctest::Approx(15.0)); // m * (Xx^2+Xy^2) = 3 * 5 = 15
    }

    TEST_CASE("pga2dp: Inertia2dp operator+= - accumulate multiple point masses")
    {
        fmt::println("pga2dp: Inertia2dp operator+= - accumulate multiple point masses");

        // Two point masses
        double m1 = 1.0;
        Vec2dp<double> X1{1.0, 0.0, 1.0};

        double m2 = 2.0;
        Vec2dp<double> X2{0.0, 1.0, 1.0};

        auto I1 = get_point_inertia(m1, X1);
        auto I2 = get_point_inertia(m2, X2);

        // Accumulate
        Inertia2dp<double> I_total{};
        I_total += I1;
        I_total += I2;

        // Verify it equals sum of individual matrices
        for (size_t i = 0; i < 9; ++i) {
            CHECK(I_total.data[i] == doctest::Approx(I1.data[i] + I2.data[i]));
        }
    }

    TEST_CASE("pga2dp: Inertia2dp operator() - matrix-vector multiplication")
    {
        fmt::println("pga2dp: Inertia2dp operator() - matrix-vector multiplication");

        // Simple case: point at origin has zero inertia contributions
        // (except diagonal for rotation about origin)
        double m = 2.0;
        Vec2dp<double> X{1.0, 0.0, 1.0}; // point at (1, 0)

        auto I = get_point_inertia(m, X);

        // Apply to rate of change Omega
        Vec2dp<double> Omega{0.0, 0.0, 1.0}; // pure rotation about origin

        BiVec2dp<double> result = I(Omega); // momentum bivector from map I(Omega)

        // Manual calculation:
        // I = m * [  0    1    0 ]     [  0    2    0 ]
        //         [ -1    0    1 ]  =  [ -2    0    2 ]
        //         [ -1    0    1 ]     [ -2    0    2 ]
        // I * [0, 0, 1]^T = [0, 2, 2]^T

        CHECK(result.x == doctest::Approx(0.0));
        CHECK(result.y == doctest::Approx(2.0));
        CHECK(result.z == doctest::Approx(2.0));
    }

    TEST_CASE("pga2dp: get_inertia_inverse - verify I * I_inv = Identity")
    {
        fmt::println("pga2dp: get_inertia_inverse - verify I * I_inv = Identity");

        // Create a non-trivial inertia matrix from multiple points
        Inertia2dp<double> I{};
        I += get_point_inertia(1.0, Vec2dp<double>{1.0, 0.0, 1.0});
        I += get_point_inertia(2.0, Vec2dp<double>{0.0, 1.0, 1.0});
        I += get_point_inertia(1.5, Vec2dp<double>{1.0, 1.0, 1.0});

        auto I_inv = get_inertia_inverse(I);

        // Verify I * I_inv ≈ Identity by checking I(I_inv(e_i)) ≈ e_i
        Vec2dp<double> e1{1.0, 0.0, 0.0};
        Vec2dp<double> e2{0.0, 1.0, 0.0};
        Vec2dp<double> e3{0.0, 0.0, 1.0};

        auto r1 = I(I_inv(e1));
        auto r2 = I(I_inv(e2));
        auto r3 = I(I_inv(e3));

        CHECK(r1.x == doctest::Approx(1.0).epsilon(1e-10));
        CHECK(r1.y == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(r1.z == doctest::Approx(0.0).epsilon(1e-10));

        CHECK(r2.x == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(r2.y == doctest::Approx(1.0).epsilon(1e-10));
        CHECK(r2.z == doctest::Approx(0.0).epsilon(1e-10));

        CHECK(r3.x == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(r3.y == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(r3.z == doctest::Approx(1.0).epsilon(1e-10));
    }

    TEST_CASE("pga2dp: compute_omega_dot - ODE right-hand side")
    {
        fmt::println("pga2dp: compute_omega_dot - ODE right-hand side");

        // Create inertia from a simple point mass
        Inertia2dp<double> I{};
        I += get_point_inertia(1.0, Vec2dp<double>{1.0, 0.0, 1.0});
        I += get_point_inertia(1.0, Vec2dp<double>{-1.0, 0.0, 1.0});

        auto I_inv = get_inertia_inverse(I);

        // Zero force (a bivector F = wdg(P,f)), zero omega (a vector rate of change)
        // -> zero omega_dot
        BiVec2dp<double> F_zero{0.0, 0.0, 0.0};
        Vec2dp<double> Omega_zero{0.0, 0.0, 0.0};

        auto Omega_dot = compute_omega_dot(I_inv, F_zero, Omega_zero, I);

        CHECK(Omega_dot.x == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(Omega_dot.y == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(Omega_dot.z == doctest::Approx(0.0).epsilon(1e-10));

        // Non-zero force, zero omega -> non-zero omega_dot (linear response)
        BiVec2dp<double> F{0.0, 0.0, 1.0}; // torque about origin
        Omega_dot = compute_omega_dot(I_inv, F, Omega_zero, I);

        // Result should be I_inv(F) since rcmt(0, ...) = 0
        auto expected = I_inv(F);
        CHECK(Omega_dot.x == doctest::Approx(expected.x).epsilon(1e-10));
        CHECK(Omega_dot.y == doctest::Approx(expected.y).epsilon(1e-10));
        CHECK(Omega_dot.z == doctest::Approx(expected.z).epsilon(1e-10));
    }

    TEST_CASE("pga2dp for 2D case (force + moment in vec2dp): pre-study linear motion")
    {
        fmt::println(
            "pga2dp for 2D case (force + moment in vec2dp): pre-study linear motion");

        // pre-study: investigate inertia matrices and their inverses
        //            result: -> not invertible for single mass point
        //                    -> invertrible for more than one mass point

        auto X1 = vec2dp{0, 1, 1};  // X1 initial position
        auto X2 = vec2dp{0, -1, 1}; // X2 initial position
        double const m = 1.0;       // mass

        inertia2dp I_tot;

        auto I_X1 = get_point_inertia(m, X1);
        I_tot += I_X1;

        auto I_X2 = get_point_inertia(m, X2);
        I_tot += I_X2;

        auto det_I_X1 = hd::ga::det(I_X1.view());
        auto det_I_X2 = hd::ga::det(I_X2.view());
        auto det_I_tot = hd::ga::det(I_tot.view());

        auto I_tot_inv = get_inertia_inverse(I_tot);

        fmt::println("I_X1       = {:>-7.3f}", I_X1);
        fmt::println("I_X2       = {:>-7.3f}", I_X2);
        fmt::println("I_tot      = {:>-7.3f}", I_tot);
        fmt::println("");
        fmt::println("det(I_X1)  = {:>-7.3f}", det_I_X1);
        fmt::println("det(I_X2)  = {:>-7.3f}", det_I_X2);
        fmt::println("det(I_tot) = {:>-7.3f}", det_I_tot);
        fmt::println("");
        fmt::println("I_tot_inv  = {:>-7.3f}", I_tot_inv);
        fmt::println("");
        fmt::println("I is not invertible for single mass point.");
        fmt::println("I is invertible for several mass points.");
        fmt::println("");

        // pre-study: momentum calculation at different positions
        //            result: -> total momentum is P = X ^ m d/dt(X) is constant

        auto v = vec2dp{1, 0, 0};   // v_x = 1 m/s, v_y = 0 m/s
        auto X11 = vec2dp{1, 1, 1}; // pos. of X1 after 1s
        auto X12 = vec2dp{2, 1, 1}; // pos. of X1 after 2s

        // auto v = vec2dp{1, 1, 0};   // v_x = 1 m/s, v_y = 1 m/s v' = |v|*rsqrt(2))
        // auto X11 = vec2dp{1, 2, 1}; // pos. of X1 after 1s
        // auto X12 = vec2dp{2, 3, 1}; // pos. of X1 after 2s

        auto P_X1 = wdg(X1, m * v);
        auto P_X11 = wdg(X11, m * v);
        auto P_X12 = wdg(X12, m * v);

        fmt::println("v         = {:>-7.3f}", v);
        fmt::println("X11       = {:>-7.3f}", X11);
        fmt::println("X12       = {:>-7.3f}", X12);
        fmt::println("P_X1      = {:>-7.3f}", P_X1);
        fmt::println("P_X11     = {:>-7.3f}", P_X11);
        fmt::println("P_X12     = {:>-7.3f}", P_X12);
        fmt::println("");
        auto I_X11 = get_point_inertia(m, X11);
        auto I_X12 = get_point_inertia(m, X12);
        fmt::println("I_X11     = {:>-7.3f}", I_X11);
        fmt::println("I_X12     = {:>-7.3f}", I_X12);
        fmt::println("");
        fmt::println("total momentum is P = X ^ m d/dt(X) const. for linear motion.");
        fmt::println("inertia changes in global system -> must be considered locally.");
        fmt::println("");

        // pre-study: linear motion expressed by rate of change vector in 2D
        //            linear motion of point in x-direction

        auto omega_tra = att(bulk_dual(v));
        auto B_1s = omega_tra * 1.0;
        auto B_2s = omega_tra * 2.0;
        //
        auto M_1s = get_motor(v * 1.0); // use translation ds = v * dt
        auto M_2s = get_motor(v * 2.0); // use translation ds = v * dt

        auto Q = vec2dp(0, 1, 1);
        auto n = vec2dp(-v.y, v.x, 0) / to_val(bulk_nrm(vec2dp(-v.y, v.x, 0)));
        auto d_1s = v * 1.0 * 0.5; // move line by 0.5 * v at 1s
        auto d_2s = v * 2.0 * 0.5; // move line by 0.5 * v at 2s

        auto B1 = wdg(Q, n);
        auto B2 = wdg(Q + d_1s, n);
        auto B3 = wdg(Q + d_2s, n);
        auto M_1s_alt = get_motor_from_lines(B1, B2);
        auto M_2s_alt = get_motor_from_lines(B1, B3);

        fmt::println("X1                      = {:>-7.3f}", X1);
        fmt::println("v                       = {:>-7.3f}", v);
        fmt::println("omega_tra               = {:>-7.3f}", omega_tra);
        fmt::println("B(1s) = omega_tra * 1s  = {:>-7.3f}", B_1s);
        fmt::println("B(2s) = omega_tra * 2s  = {:>-7.3f}", B_2s);
        fmt::println("M(1s) = rexp(0.5 * B_1s) = {:>-7.3f}", M_1s);
        fmt::println("M(1s)_alt               = {:>-7.3f}", M_1s_alt);
        fmt::println("M(2s) = rexp(0.5 * B_2s) = {:>-7.3f}", M_2s);
        fmt::println("M(2s)_alt               = {:>-7.3f}", M_2s_alt);
        fmt::println("");
        fmt::println("gr1(rgpr(rgpr(M, v), rrev(M)))) = {:>-7.3f}", move2dp(X1, M_1s));
        fmt::println("gr1(rgpr(rgpr(M, v), rrev(M)))) = {:>-7.3f}", move2dp(X1, M_2s));
        fmt::println("");

        CHECK(X11 == move2dp(X1, M_1s));
        CHECK(X12 == move2dp(X1, M_2s));
        CHECK(M_1s == rexp(0.5 * B_1s));
        CHECK(M_2s == rexp(0.5 * B_2s));

        // pre-study: speed defined from constant linear motion
        //            result: -> speed must be const (position independent)

        auto omega_tra_w = att(bulk_dual(v));

        auto v_X1 = rcmt(omega_tra_w, X1);
        auto v_X11 = rcmt(omega_tra_w, X11);
        auto v_X12 = rcmt(omega_tra_w, X12);

        fmt::println("v_X1  = {:>-7.3f}", v_X1);
        fmt::println("v_X11 = {:>-7.3f}", v_X11);
        fmt::println("v_X12 = {:>-7.3f}", v_X12);
        fmt::println("");
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // PGA2DP: discrete vs. continuous inertia
    //////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga2dp: discrete inertia of rectangular plate (25x25 grid)")
    {
        fmt::println("pga2dp: discrete inertia of rectangular plate (25x25 grid)");
        fmt::println("");
        fmt::println("Plate: e1 in [-2, 2], e2 in [-1, 1], 25x25 equidistant points");
        fmt::println("Total mass M = 1.0, uniform distribution");
        fmt::println("");

        // Grid parameters
        int const nx = 25;
        int const ny = 25;
        double const x_min = -2.0, x_max = 2.0;
        double const y_min = -1.0, y_max = 1.0;
        double const dx = (x_max - x_min) / (nx - 1);
        double const dy = (y_max - y_min) / (ny - 1);

        // Total mass = 1.0 distributed uniformly over all grid points
        double const M_total = 1.0;
        double const m_pt = M_total / (nx * ny); // mass per point

        // Accumulate inertia over all grid points
        Inertia2dp<double> I_grid{};
        for (int ix = 0; ix < nx; ++ix) {
            double const x = x_min + ix * dx;
            for (int iy = 0; iy < ny; ++iy) {
                double const y = y_min + iy * dy;
                Vec2dp<double> const X{x, y, 1.0}; // unitized point (Xz = 1)
                I_grid += get_point_inertia(m_pt, X);
            }
        }

        auto const det_I_grid = hd::ga::det(I_grid.view());
        auto const I_grid_inv = get_inertia_inverse(I_grid);

        fmt::println("I_grid     = {:>-10.6f}", I_grid);
        fmt::println("");
        fmt::println("det(I_grid)= {:>-10.6f}", det_I_grid);
        fmt::println("");
        fmt::println("I_grid_inv = {:>-10.6f}", I_grid_inv);
        fmt::println("");

        // Analytical continuous limit for a rectangular plate of width W=4, height H=2,
        // uniform density, centered at origin (Xz=1 for all points):
        //   I[0,1] = M,  I[1,0] = -M  (total mass on off-diagonals)
        //   I[0,2] = -M*mean(y) = 0,  I[1,2] = M*mean(x) = 0  (center at origin)
        //   I[2,2] = M*(W^2+H^2)/12 = M*(16+4)/12 = M*20/12 = 5/3*M ≈ 1.6667
        //
        // Note on discretization error for I[2,2]:
        //   Endpoint-inclusive equidistant grids overestimate mean(x^2) and mean(y^2)
        //   compared to the continuous integral because boundary points carry the same
        //   weight as interior points. Midpoint or trapezoidal sampling converges faster.
        //   For a 25-point grid from a..b the discrete mean(x^2) = (a^2+ab+b^2)/3
        //   minus a small correction, giving ~8% overestimate vs. the integral
        //   (a^2+ab+b^2)/3 for the symmetric case. As n->inf the grid mean converges to
        //   M*(W^2+H^2)/12.
        double const W = x_max - x_min; // 4.0
        double const H = y_max - y_min; // 2.0
        double const I22_continuous = M_total * (W * W + H * H) / 12.0;

        fmt::println("Continuous limit (uniform rectangular plate):");
        fmt::println("  I[0,1] =  M           = {:>-10.6f}  (discrete: {:>-10.6f})",
                     M_total, I_grid.view()[0, 1]);
        fmt::println("  I[1,0] = -M           = {:>-10.6f}  (discrete: {:>-10.6f})",
                     -M_total, I_grid.view()[1, 0]);
        fmt::println("  I[2,2] = M*(W²+H²)/12 = {:>-10.6f}  (discrete: {:>-10.6f})",
                     I22_continuous, I_grid.view()[2, 2]);
        fmt::println(
            "  (discrete overestimates I[2,2] by ~{:.1f}% due to endpoint sampling)",
            100.0 * (I_grid.view()[2, 2] - I22_continuous) / I22_continuous);
        fmt::println("");

        // Exact continuous result from get_plate_inertia() - no discretization error
        auto const I_plate = get_plate_inertia(M_total, W, H);
        auto const I_plate_inv = get_inertia_inverse(I_plate);

        fmt::println("get_plate_inertia() (exact continuous limit):");
        fmt::println("I_plate    = {:>-10.6f}", I_plate);
        fmt::println("");
        fmt::println("I_plate_inv= {:>-10.6f}", I_plate_inv);
        fmt::println("");
        fmt::println("Comparison: I_grid[2,2]={:.6f}  I_plate[2,2]={:.6f}  diff={:.6f}",
                     I_grid.view()[2, 2], I_plate.view()[2, 2],
                     I_grid.view()[2, 2] - I_plate.view()[2, 2]);
        fmt::println("");

        // get_plate_inertia() matches the continuous limit exactly
        CHECK(I_plate.view()[0, 0] == doctest::Approx(0.0).epsilon(1e-15));
        CHECK(I_plate.view()[0, 1] == doctest::Approx(M_total).epsilon(1e-15));
        CHECK(I_plate.view()[0, 2] == doctest::Approx(0.0).epsilon(1e-15));
        CHECK(I_plate.view()[1, 0] == doctest::Approx(-M_total).epsilon(1e-15));
        CHECK(I_plate.view()[1, 1] == doctest::Approx(0.0).epsilon(1e-15));
        CHECK(I_plate.view()[1, 2] == doctest::Approx(0.0).epsilon(1e-15));
        CHECK(I_plate.view()[2, 0] == doctest::Approx(0.0).epsilon(1e-15));
        CHECK(I_plate.view()[2, 1] == doctest::Approx(0.0).epsilon(1e-15));
        CHECK(I_plate.view()[2, 2] == doctest::Approx(I22_continuous).epsilon(1e-15));

        // Discrete grid overestimates I[2,2] vs. the exact continuous plate
        CHECK(I_grid.view()[2, 2] > I_plate.view()[2, 2]);

        // Exact properties of the discrete grid: mass sum and center-of-mass symmetry
        CHECK(I_grid.view()[0, 1] == doctest::Approx(M_total).epsilon(1e-10));
        CHECK(I_grid.view()[1, 0] == doctest::Approx(-M_total).epsilon(1e-10));
        CHECK(I_grid.view()[0, 2] == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(I_grid.view()[1, 2] == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(I_grid.view()[2, 0] == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(I_grid.view()[2, 1] == doctest::Approx(0.0).epsilon(1e-10));
        // I[2,2] converges to continuous limit from above; verify known 25-pt value.
        // Exact discrete value for 25-pt endpoint-inclusive grid:
        //   mean(x^2) = 13/9, mean(y^2) = 13/36  =>  I[2,2] = 65/36 ≈ 1.805556
        CHECK(I_grid.view()[2, 2] == doctest::Approx(65.0 / 36.0).epsilon(1e-6));
    }

    TEST_CASE("pga2dp: get_disc_inertia / make_disc_body")
    {
        fmt::println("pga2dp: get_disc_inertia (polar moment m r^2 / 2)");

        value_t const m = 3.0, r = 2.0;
        auto const I = get_disc_inertia(m, r);

        // mass block + polar moment about the centroid J = m r^2 / 2
        CHECK(I.view()[0, 1] == doctest::Approx(m));
        CHECK(I.view()[1, 0] == doctest::Approx(-m));
        CHECK(I.view()[2, 2] == doctest::Approx(0.5 * m * r * r)); // = 6.0
        // off-diagonal coupling terms are zero
        CHECK(I.view()[0, 2] == doctest::Approx(0.0));
        CHECK(I.view()[1, 2] == doctest::Approx(0.0));
        CHECK(I.view()[2, 0] == doctest::Approx(0.0));
        CHECK(I.view()[2, 1] == doctest::Approx(0.0));

        // I * I_inv == identity
        auto const I_inv = get_inertia_inverse(I);
        auto const P = I.view(), Q = I_inv.view();
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j) {
                value_t s = 0.0;
                for (size_t k = 0; k < 3; ++k)
                    s += P[i, k] * Q[k, j];
                CHECK(s == doctest::Approx(i == j ? 1.0 : 0.0).epsilon(1e-12));
            }

        // Steiner: rotating about a rim point (pivot at distance r) adds m r^2
        auto const I_rim = get_disc_inertia(m, r, vec2dp{r, 0.0, 1.0});
        CHECK(I_rim.view()[2, 2] ==
              doctest::Approx(0.5 * m * r * r + m * r * r)); // 3/2 m r^2

        // make_disc_body wires the inertia + its inverse + mass
        auto const body = make_disc_body(m, r);
        CHECK(body.mass == doctest::Approx(m));
        CHECK(body.I.view()[2, 2] == doctest::Approx(0.5 * m * r * r));
        fmt::println("");
    }

    TEST_CASE("pga2dp: moment_about - torque of a force line about a point")
    {
        fmt::println("pga2dp: moment_about - torque of a force line about a point");

        // same expression as in 3D: bulk(F) - r ^ att(F). In the plane the result's e12
        // component is the scalar torque r x f with the physical sign; the incidence
        // R ^ F of the point with the force line (the planar form in ga_docu) is the same
        // quantity with the orientation of this algebra's pseudoscalar e321, i.e. -M.
        vec2dp const P{1.5, -0.7, 1.0};
        vec2dp const f{0.3, 1.1, 0.0};
        vec2dp const R{-2.0, 0.8, 1.0};
        bivec2dp const F = wdg(P, f);

        bivec2dp const M = moment_about(R, F);
        CHECK(M.z == doctest::Approx((P.x - R.x) * f.y - (P.y - R.y) * f.x)); // r x f
        CHECK(M.z == doctest::Approx(-value_t(wdg(R, F))));                   // -(R ^ F)
        CHECK(std::abs(M.x) + std::abs(M.y) < 1e-14); // free: no force part

        // independent of the point used to build the line, and of R's scale
        bivec2dp const F2 = wdg(vec2dp{P.x + 3.0 * f.x, P.y + 3.0 * f.y, 1.0}, f);
        CHECK(moment_about(R, F2).z == doctest::Approx(M.z));
        CHECK(moment_about(vec2dp{-4.0, 1.6, 2.0}, F).z == doctest::Approx(M.z));

        // zero on the line of action; additive over force lines
        CHECK(std::abs(moment_about(vec2dp{P.x + 0.5 * f.x, P.y + 0.5 * f.y, 1.0}, F).z) <
              1e-14);
        bivec2dp const G = wdg(vec2dp{0.2, 0.4, 1.0}, vec2dp{-0.9, 0.2, 0.0});
        CHECK(moment_about(R, F + G).z ==
              doctest::Approx(moment_about(R, F).z + moment_about(R, G).z));
        fmt::println("");
    }

} // TEST_SUITE("PGA2DP: physics tests prep")

/////////////////////////////////////////////////////////////////////////////////////////
// PGA2DP physics implementation
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA2DP: physics tests implementation")
{

    TEST_CASE("ega3d: linear motion (multi-body system; force+moment in vec3d)")
    {
        fmt::println("ega3d: linear motion (multi-body system; force+moment in vec3d)");

        // hint: rigid body is simulated as n points (integrated separately)
        //       where all see adequate forces to assure movement as rigid body
        //
        // TODO: Supports independent multi-body problem. Thus could be simplified to just
        //       integrate movement of center of gravity and its orientation to solve with
        //       minimal effort.
        //       However, the idea here is just to get a reference solution for
        //       comparision with pga2dp, not to get a solution with minimal effort

        class sim_ode_rigid_ega3d { // model 2nd order ode by a 1st order system
                                    // integrate a -> v -> s,
                                    // i.e. from acceleration via velocity to position

          public:

            sim_ode_rigid_ega3d(size_t npts_in) :
                npts(npts_in), u_mem(npts * 2), uh_mem(2 * npts * 2), rhs_mem(npts * 2)
            {
            }

            void set_initial_values()
            {

                // Create mdspan view for setting initial values
                auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);

                // use vec3d to encode 2d position/velocity in x- and y-component
                // of 3d vector and angle/angular velocity in z-component of 3d vector

                // Set initial state for all n points n = [0,n) and their components:
                // u[n,0] = position (x  ,   y,   phi)
                // u[n,1] = velocity (v_x, v_y, omega)


                size_t n = 0;
                u[n, 0] = vec3d{0, 1, 0}; // initial position (x,   y,   phi)
                u[n, 1] = vec3d{1, 0, 0}; // initial velocity (v_x, v_y, omega)

                n = 1;
                u[n, 0] = vec3d{0, -1, 0}; // initial position
                u[n, 1] = vec3d{1, 0, 0};  // initial velocity
            }

            void calc_rhs()
            {

                // Create mdspan views for rhs calculation
                auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);
                auto rhs = mdspan<vec3d, dextents<size_t, 2>>(rhs_mem.data(), npts, 2);

                for (size_t n = 0; n < npts; ++n) {

                    // get current state
                    vec3d position = u[n, 0]; // position is stored as point in u[n,0]
                    vec3d velocity = u[n, 1]; // velocity is stored as vector in u[n,1]

                    // forces and torques to change linear and angular acceleration
                    // auto force = vec2d{0.0, 0.0}; // no force (x,y)
                    auto force = vec2d{0.0, -m * 9.81}; // gravitation (x=0, y=m*g)
                    auto torque = 0.0;                  // no torque (z)

                    vec3d acceleration;
                    acceleration.x = velocity.z * velocity.y;  // coupling (v_y, omega)
                    acceleration.x += force.x / m;             // impact of force_x
                    acceleration.y = -velocity.z * velocity.x; // coupling (v_x, omega)
                    acceleration.y += force.y / m;             // impact of force_y

                    // J = m * r^2 with r being the distance from rot-axis
                    auto J = m * nrm_sq(position - O_3d);
                    acceleration.z = torque / J; // angular acceleration in (z)
                                                 // (only valid for const inertia J!)

                    // Set right-hand side for ODE system
                    // u[n,0]' = velocity (linear and angular)
                    // u[n,1]' = acceleration (linear and angular)
                    rhs[n, 0] = velocity;
                    rhs[n, 1] = acceleration;
                }
            }

            void calc_rkstep(double dt)
            {
                // Create mdspan views for RK4 integration (for all npts points)
                auto u = mdspan<vec3d, dextents<size_t, 1>>(u_mem.data(), npts * 2);
                auto uh = mdspan<vec3d, dextents<size_t, 2>>(uh_mem.data(), 2, npts * 2);
                auto rhs =
                    mdspan<vec3d const, dextents<size_t, 1>>(rhs_mem.data(), npts * 2);

                // Perform RK4 integration (4 sub-steps)
                for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
                    calc_rhs();
                    rk4_step(u, uh, rhs, dt, rk_step);
                }
            }

            void print_sim(double t)
            {
                // Create mdspan view for printing
                auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);

                if (npts == 1) {

                    fmt::print("t = {:>-7.3f}", t);

                    // get current state
                    vec3d position = u[0, 0]; // position is stored as point in u[n,0]
                    vec3d velocity = u[0, 1]; // velocity is stored as vector in u[n,1]

                    fmt::println(", pos. = {:>-7.3f}, vel. = {:>-7.3f}", position,
                                 velocity);
                }
                else {

                    fmt::println("t = {:>-7.3f}:", t);

                    for (size_t n = 0; n < npts; ++n) {

                        // get current state
                        vec3d position = u[n, 0];
                        vec3d velocity = u[n, 1];

                        fmt::println("    n = {}, pos. = {:>-7.3f}, vel. = {:>-7.3f}", n,
                                     position, velocity);
                    }
                }
            }

          private:

            // number of independent points to solve system for
            size_t npts;

            // RK4 integration state for point n with system order = 2
            // => [n+0: position, n+1: velocity])
            std::vector<vec3d> u_mem;   // [n+0: position, n+1: velocity]
            std::vector<vec3d> uh_mem;  // helper for integration
            std::vector<vec3d> rhs_mem; // right-hand side values

            // mass (= linear inertia)
            double m = 1.0;
        };

        auto constexpr num_points = 2;
        sim_ode_rigid_ega3d sim(num_points);

        // time range
        auto t_rng = discrete_range(0.0, 2.0, 20);

        sim.set_initial_values();
        sim.print_sim(t_rng.min());

        for (size_t n = 1; n <= t_rng.steps(); n++) {

            // integration from t to t + dt
            double t = t_rng.min() + n * t_rng.delta();
            sim.calc_rkstep(t_rng.delta());

            // print sim status at t+dt
            sim.print_sim(t);
        }

        fmt::println("");
    }


    TEST_CASE("pga2dp: linear motion (rigid-body system)")
    {
        fmt::println("pga2dp: linear motion (rigid-body system)");


        class sim_ode_rigid_pga2dp { // model 2nd order ode for rigid body by a 1st order
                                     // system

            // actually only for the integration of the center of mass and its orientation
            // all points constituting the rigid body will be transformed accordingly

          public:

            sim_ode_rigid_pga2dp(size_t npts_in, std::vector<vec2dp> pos_in,
                                 std::vector<value_t> m_in, vec2dp const& cm_pos_in,
                                 vec2dp const& cm_spd_in, value_t cm_phi_in,
                                 value_t cm_omega_in) :
                npts(npts_in), pts(std::move(pos_in)), m(std::move(m_in)),
                cm_w_pos0(cm_pos_in), cm_w_spd0(cm_spd_in), cm_w_phi0(cm_phi_in),
                cm_w_omega0(cm_omega_in), u_mem(2), uh_mem(2 * 2), rhs_mem(2)
            {
                if (npts < 2) {
                    throw std::invalid_argument("sim_ode: rigid body simulation requires"
                                                " at least two discrete points.");
                }
                if ((pts.size() != npts) || (m.size() != npts)) {
                    throw std::invalid_argument(
                        "sim_ode: position and mass vectors require same size.");
                }

                for (size_t n = 0; n < npts; ++n) {
                    cg += pts[n];
                    m_tot += m[n]; // total mass
                }
                cg *= (1.0 / npts); // center of gravity is now available in cg

                if (bulk_nrm(cg - O_2dp) > eps) {
                    throw std::invalid_argument(
                        "sim_ode: center of gravity must be located at O_2dp.");
                }

                // build inertia map from discrete points
                // maps rate of change (a vector in 2D) to the momentum bivector
                for (size_t n = 0; n < npts; ++n) {
                    I += get_point_inertia(m[n], pts[n]); // total inertia map
                }

                // get inverse inertia map
                // maps a momentum bivector to the rate of change (a vector in 2D)
                I_inv = get_inertia_inverse(I);
            }

            void set_initial_values()
            {
                // HINT: currently input variables for initial position NOT used;
                //       for now assume "no initial transformation for position"

                // Create mdspan view for setting initial values
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);

                // initial transformation of "position" encoded in B-vector
                // encoding:
                //  B_rot=(x0_fix, y0_phi, 1) * phi0
                //  B_tra=(-y0_trans, x0_trans, 0)
                //  B    = B_rot + B_tra
                //
                //  CAUTION (generators add only at the velocity level):
                //  B = B_rot + B_tra is the additive split of ONE twist. rexp(0.5*B) of
                //  that single B is ONE finite screw (3D) / rotation about a shifted
                //  pivot (2D) -- NOT "translate, then rotate". To compose two SEPARATE
                //  finite motions, MULTIPLY their motors: M = rgpr(M2, M1), which is
                //  order-dependent (non-commutative). Adding generators == multiplying
                //  motors only when they commute: for translation+rotation that is
                //  solely the 3D screw case, never in 2D.
                //
                //  typical starting value is "no initial transformation":
                //  B0=(0,0,0) => M0 = rexp(0.5 * B0) = pscalar2dp(1)
                //                (=identity transformation at t=0)
                //
                u[0] = vec2dp(0.0, 0.0, 0.0);

                // initial rate of change transformation of "velocity" dB(t0)/dt = Omega0
                // encoding:
                // rotation:    Omega0_rot = ( q.x0, q.y0, 1) * omega0
                // translation: Omega0_tra = (-v0.y, v0.x, 0)
                //
                // resulting Omega0 = Omega0_rot + Omega0_tra

                // case with initial translation
                u[1] = vec2dp(-cm_w_spd0.y, cm_w_spd0.x, cm_w_omega0);
            }

            void calc_rhs()
            {

                // Create mdspan views for rhs calculation
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto rhs = mdspan<vec2dp, dextents<size_t, 1>>(rhs_mem.data(), 2);


                // get current state
                [[maybe_unused]] vec2dp B = u[0]; // position transformation B is in u[0]
                vec2dp Omega = u[1]; // velocity trafo d(B)/dt = Omega is in u[1]

                // forces and torques to change linear and angular acceleration

                // no force f acting the body frame at the origin (F_b = O_2dp ^ f)
                // auto F_b = bivec2dp{0.0, 0.0, 0.0};

                // gravitation acting through center of mass in in the body system
                // with center of mass being the origin O_2dp of the body system
                // gravitation f = (x=0, y=-m_tot*g, z=0) acting at O_2dp
                // where the y-axis is assumed pointing upwards
                auto F_b = wdg(O_2dp, vec2dp{0.0, -m_tot * 9.81, 0.0});

                // Set right-hand side for ODE system:
                // u[0]' = velocity trafo Omega = d(B)/dt (linear and angular)
                // u[1]' = acceleration trafo d(Omega)/dt (linear and angular)
                rhs[0] = Omega;
                rhs[1] = compute_omega_dot(I_inv, F_b, Omega, I);
            }

            void calc_rkstep(double dt)
            {
                // Create mdspan views for RK4 integration (for all npts points)
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto uh = mdspan<vec2dp, dextents<size_t, 2>>(uh_mem.data(), 2, 2);
                auto rhs = mdspan<vec2dp const, dextents<size_t, 1>>(rhs_mem.data(), 2);

                // Perform RK4 integration (4 sub-steps)
                for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
                    calc_rhs();
                    rk4_step(u, uh, rhs, dt, rk_step);
                }
            }

            void print_sim(double t)
            {
                // Create mdspan view for printing
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);


                fmt::println("t = {:>-7.3f}:", t);

                for (size_t n = 0; n < npts; ++n) {

                    // get current state (= current positional transformation bivector)
                    vec2dp B = u[0];     // B = Omega * t + B0 (from integration)
                    vec2dp Omega = u[1]; // dB/dt = Omega = dB^2/dt^2 * t + Omega0

                    // calculate current position from B via M = 0.5 * B) ⟇ M0
                    // and via pts(t) = M ⟇ pts(t0) ⟇ rrev(M)
                    //
                    auto M = rexp(0.5 * B);
                    auto pt = move2dp(pts[n], M);
                    // fmt::println(
                    //     "    n = {}, B = {:>-7.3f}, Omega = {:>-7.3f}, M = {:>-7.3f}",
                    //     n, B, Omega, M);
                    fmt::println("    n = {}, pos. = {:>-7.3f}, vel. = {:>-7.3f}", n, pt,
                                 rcmt(Omega, pt));
                }
            }

          private:

            // number of independent points to solve system for
            size_t npts;

            // positions of points in the body frame (fixed for now to get started)
            std::vector<vec2dp> pts; // point positions (fixed for now)
            std::vector<value_t> m;  // mass of points

            // initial position and speed of center of mass relative to world system
            vec2dp cm_w_pos0; // initial position of center of mass (unitized)
            vec2dp cm_w_spd0; // initial speed of center of mass (.z == 0)

            // initial angular position of cm (body vs. world frame, e12 defines pos phi)
            [[maybe_unused]] value_t cm_w_phi0;

            // initial angular velocity of cm (body vs. world frame)
            value_t cm_w_omega0;

            // RK4 integration state for point n with system order = 2
            // => [n+0: position, n+1: velocity])
            std::vector<vec2dp> u_mem;   // [n+0: position, n+1: velocity]
            std::vector<vec2dp> uh_mem;  // helper for integration
            std::vector<vec2dp> rhs_mem; // right-hand side values

            // inertia map and its inverse (calculated from descrete input values)
            // calculation is done in body frame
            inertia2dp I;
            inertia2dp I_inv;

            value_t m_tot; // total mass (needed for force calculation, e.g. gravity)
            vec2dp cg;     // center of gravity (in body frame == O_b)
        };

        auto constexpr num_points = 2;
        // (initial) positions and masses of each mass point in the body frame
        std::vector<vec2dp> pos_b = {vec2dp{0, 1, 1}, vec2dp{0, -1, 1}}; // position
        std::vector<value_t> m = {1.0, 1.0};                             // mass

        // initial position and movement of body frame vs. world frame (in world frame)
        auto cm_pos = O_2dp;           // world and body frame coincide at t=0
        auto cm_spd = vec2dp{1, 0, 0}; // body frame is moving with v0 vs. world frame

        // initial rotation and anguar velocity of body frame vs. world frame
        auto cm_phi = 0.0;   // no initial rotation
        auto cm_omega = 0.0; // no initial angular velocity

        sim_ode_rigid_pga2dp sim(num_points, pos_b, m, cm_pos, cm_spd, cm_phi, cm_omega);

        // time range (from, to, number of steps)
        auto t_rng = discrete_range(0.0, 2.0, 20);

        sim.set_initial_values();
        sim.print_sim(t_rng.min());

        for (size_t n = 1; n <= t_rng.steps(); n++) {

            // integration from t to t + dt
            double t = t_rng.min() + n * t_rng.delta();
            sim.calc_rkstep(t_rng.delta());

            // print sim status at t+dt
            sim.print_sim(t);
        }

        fmt::println("");
    }


    TEST_CASE("ega3d: circular motion (multi-body system; force+moment in vec3d)")
    {
        fmt::println("ega3d: circular motion (multi-body system; force+moment in vec3d)");

        // hint: rigid body is simulated as n points (integrated separately)
        //       where all see adequate forces to assure movement as rigid body
        //
        // TODO: Supports independent multi-body problem. Thus could be simplified to just
        //       integrate movement of center of gravity and its orientation to solve with
        //       minimal effort.
        //       However, the idea here is just to get a reference solution for
        //       comparision with pga2dp, not to get a solution with minimal effort

        class sim_ode_rigid_ega3d { // model 2nd order ode by a 1st order system
                                    // integrate a -> v -> s,
                                    // i.e. from acceleration via velocity to position
          public:

            sim_ode_rigid_ega3d(size_t npts_in) :
                npts(npts_in), u_mem(npts * 2), uh_mem(2 * npts * 2), rhs_mem(npts * 2)
            {
            }

            void set_initial_values()
            {

                // Create mdspan view for setting initial values
                auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);

                // use vec3d to encode 2d position/velocity in x- and y-component
                // of 3d vector and angle/angular velocity in z-component of 3d vector

                // Set initial state for all n points n = [0,n) and their components:
                // u[n,0] = position (x  ,   y,   phi)
                // u[n,1] = velocity (v_x, v_y, omega)

                size_t n = 0;
                u[n, 0] = vec3d{0, 1, 0}; // initial position (x,   y,   phi)
                u[n, 1] = vec3d{1, 0, 0}; // initial velocity (v_x, v_y, omega)

                n = 1;
                u[n, 0] = vec3d{0, -1, 0}; // initial position
                u[n, 1] = vec3d{-1, 0, 0}; // initial velocity
            }

            void calc_rhs()
            {

                // Create mdspan views for rhs calculation
                auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);
                auto rhs = mdspan<vec3d, dextents<size_t, 2>>(rhs_mem.data(), npts, 2);

                for (size_t n = 0; n < npts; ++n) {

                    // get current state
                    vec3d position = u[n, 0]; // position is stored as point in u[n,0]
                    vec3d velocity = u[n, 1]; // velocity is stored as vector in u[n,1]

                    // calculate forces and torques to change linear and angular
                    // acceleration

                    // HINT: this implements a multi-body solution with internal forces,
                    //       would not be required to this extend for a rigid body system

                    // equilibrium force for circular motion at given v and r:
                    // F = m * a_r = m * r * omega^2 = m * r * (v / r)^2 = m * v^2 / r
                    //
                    // T = 2*pi/omega = 2*pi/(v/r) = 2*pi*r/v

                    // keep difference to origin for comparison with 2dp case later on
                    auto r = nrm(position - O_3d);
                    auto force = -normalize(position - O_3d) * m * nrm_sq(velocity) / r;
                    auto torque = 0.0; // no torque (z)

                    // rhs of equation d2/dt^2(x) = a = F/m (linear)
                    // rhs of equation d2/dt^2(phi) = d/dt(omega) = M/J (angular)
                    vec3d acceleration;
                    acceleration.x = velocity.z * velocity.y;  // coupling (v_y, omega)
                    acceleration.x += force.x / m;             // impact of force_x
                    acceleration.y = -velocity.z * velocity.x; // coupling (v_x, omega)
                    acceleration.y += force.y / m;             // impact of force_y

                    // J = m * r^2 with r being the distance from rot-axis
                    auto J = m * r * r;
                    acceleration.z = torque / J; // angular part in (z)
                                                 // (only valid for const inertia J!)

                    // Set right-hand side for ODE system
                    // u[n,0]' = velocity (linear and angular)
                    // u[n,1]' = acceleration (linear and angular)
                    rhs[n, 0] = velocity;
                    rhs[n, 1] = acceleration;
                }
            }

            void calc_rkstep(double dt)
            {
                // Create mdspan views for RK4 integration (for all npts points)
                auto u = mdspan<vec3d, dextents<size_t, 1>>(u_mem.data(), npts * 2);
                auto uh = mdspan<vec3d, dextents<size_t, 2>>(uh_mem.data(), 2, npts * 2);
                auto rhs =
                    mdspan<vec3d const, dextents<size_t, 1>>(rhs_mem.data(), npts * 2);

                // Perform RK4 integration (4 sub-steps)
                for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
                    calc_rhs();
                    rk4_step(u, uh, rhs, dt, rk_step);
                }
            }

            void print_sim(double t)
            {
                // Create mdspan view for printing
                auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);

                if (npts == 1) {

                    fmt::print("t = {:>-7.3f}", t);

                    // get current state
                    vec3d position = u[0, 0]; // position is stored as point in u[n,0]
                    vec3d velocity = u[0, 1]; // velocity is stored as vector in u[n,1]

                    fmt::println(", pos. = {:>-7.3f}, vel. = {:>-7.3f}", position,
                                 velocity);
                }
                else {

                    fmt::println("t = {:>-7.3f}:", t);

                    for (size_t n = 0; n < npts; ++n) {

                        // get current state
                        vec3d position = u[n, 0];
                        vec3d velocity = u[n, 1];

                        fmt::println("    n = {}, pos. = {:>-7.3f}, vel. = {:>-7.3f}", n,
                                     position, velocity);
                    }
                }
            }

          private:

            // number of independent points to solve system for
            size_t npts;

            // RK4 integration state for point n with system order = 2
            // => [n+0: position, n+1: velocity])
            std::vector<vec3d> u_mem;   // [n+0: position, n+1: velocity]
            std::vector<vec3d> uh_mem;  // helper for integration
            std::vector<vec3d> rhs_mem; // right-hand side values

            // mass (= linear inertia)
            double m = 1.0;
        };

        auto constexpr num_points = 2;
        sim_ode_rigid_ega3d sim(num_points);

        // time range (1 full revolution)
        auto t_rng = discrete_range(0.0, 6.2832, 60);

        sim.set_initial_values();
        sim.print_sim(t_rng.min());

        for (size_t n = 1; n <= t_rng.steps(); n++) {

            // integration from t to t + dt
            double t = t_rng.min() + n * t_rng.delta();
            sim.calc_rkstep(t_rng.delta());

            // print sim status at t+dt
            sim.print_sim(t);
        }

        fmt::println("");
    }


    TEST_CASE("pga2dp: circular motion (rigid body system)")
    {
        fmt::println("pga2dp: circular motion (rigid body system)");

        class sim_ode_rigid_pga2dp { // model 2nd order ode by a 1st order system

            // actually only for the integration of the center of mass and its orientation
            // all points constituting the rigid body will be transformed accordingly

          public:

            sim_ode_rigid_pga2dp(size_t npts_in, std::vector<vec2dp> pos_in,
                                 std::vector<value_t> m_in, vec2dp const& cm_pos_in,
                                 vec2dp const& cm_spd_in, value_t cm_phi_in,
                                 value_t cm_omega_in) :
                npts(npts_in), pts(std::move(pos_in)), m(std::move(m_in)),
                cm_w_pos0(cm_pos_in), cm_w_spd0(cm_spd_in), cm_w_phi0(cm_phi_in),
                cm_w_omega0(cm_omega_in), u_mem(2), uh_mem(2 * 2), rhs_mem(2)
            {
                if (npts < 2) {
                    throw std::invalid_argument("sim_ode: rigid body simulation requires"
                                                " at least two discrete points.");
                }
                if ((pts.size() != npts) || (m.size() != npts)) {
                    throw std::invalid_argument(
                        "sim_ode: position and mass vectors require same size.");
                }

                fmt::println("sim_ode_rigid_pga2dp: circlular motion.");
                fmt::println("points:");
                for (size_t n = 0; n < npts; ++n) {
                    fmt::println("n = {}, m = {}, pt = {:>-7.3f}", n, m[n], pts[n]);
                }
                fmt::println("point inertia:");
                // build inertia map from discrete points
                // maps rate of change (a vector in 2D) to the momentum bivector
                for (size_t n = 0; n < npts; ++n) {
                    auto Ip = get_point_inertia(m[n], pts[n]);
                    fmt::println("n = {}, Ip = {:>-7.3f}", n, Ip);
                    I += Ip;
                }

                // get inverse inertia map
                // maps a momentum bivector to the rate of change (a vector in 2D)
                I_inv = get_inertia_inverse(I);

                fmt::println("system inertia:");
                fmt::println("I     = {:>-7.3f}", I);
                fmt::println("I_inv = {:>-7.3f}", I_inv);
            }

            void set_initial_values()
            {
                // HINT: currently input variables for initial position NOT used;
                //       for now assume "no initial transformation for position"

                // Create mdspan view for setting initial values
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);


                // initial transformation of "position" encoded in B-vector
                // encoding:
                //  B_rot=(x0_fix, y0_phi, 1) * phi0
                //  B_tra=(-y0_trans, x0_trans, 0)
                //  B    = B_rot + B_tra
                //
                //  CAUTION (generators add only at the velocity level):
                //  B = B_rot + B_tra is the additive split of ONE twist. rexp(0.5*B) of
                //  that single B is ONE finite screw (3D) / rotation about a shifted
                //  pivot (2D) -- NOT "translate, then rotate". To compose two SEPARATE
                //  finite motions, MULTIPLY their motors: M = rgpr(M2, M1), which is
                //  order-dependent (non-commutative). Adding generators == multiplying
                //  motors only when they commute: for translation+rotation that is
                //  solely the 3D screw case, never in 2D.
                //
                //  typical starting value is "no initial transformation":
                //  B0=(0,0,0) => M0 = rexp(0.5 * B0) = pscalar(1)
                //                (=identity transformation at t=0)
                //
                u[0] = vec2dp(0.0, 0.0, 0.0);

                // initial rate of change transformation of "velocity" dB/dt = Omega
                // encoding:
                // rotation: Omega_rot = (q.x0, q.y0, 1) * omega0 (Q is the
                // fixed-point) translation: Omega_tra = (-v0.y, v0.x, 0)
                //
                // resulting Omega = Omega_rot + Omega_tra

                // case with initial rotation
                //
                u[1] = wdg(O_2dp, scalar2dp(cm_w_omega0));
                // fmt::println("n = {}, u = {}", n, u[n, 1]);
            }

            void calc_rhs()
            {

                // Create mdspan views for rhs calculation
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto rhs = mdspan<vec2dp, dextents<size_t, 1>>(rhs_mem.data(), 2);

                // get current state
                [[maybe_unused]] vec2dp B =
                    u[0];            // position transformation B is in u[n,0]
                vec2dp Omega = u[1]; // velocity trafo d(B)/dt = Omega is in u[n,1]

                // forces and torques to change linear and angular acceleration

                // no force f acting the body frame at the origin (F_b = O_2dp ^ f)
                // internal forces not relevant for rigid body
                auto F_b = bivec2dp{0.0, 0.0, 0.0};

                // gravitation acting through center of mass in in the body system
                // with center of mass being the origin O_2dp of the body system
                // gravitation f = (x=0, y=-m*g, z=0) acting at O_2dp
                // where the y-axis is assumed pointing upwards
                // auto F_b = wdg(O_2dp, vec2dp{0.0, -m[n] * 9.81, 0.0});

                // Set right-hand side for ODE system:
                // u[0]' = velocity trafo Omega = d(B)/dt (linear and angular)
                // u[1]' = acceleration trafo d(Omega)/dt (linear and angular)
                rhs[0] = Omega;
                rhs[1] = compute_omega_dot(I_inv, F_b, Omega, I);
            }

            void calc_rkstep(double dt)
            {
                // Create mdspan views for RK4 integration (for all npts points)
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto uh = mdspan<vec2dp, dextents<size_t, 2>>(uh_mem.data(), 2, 2);
                auto rhs = mdspan<vec2dp const, dextents<size_t, 1>>(rhs_mem.data(), 2);

                // Perform RK4 integration (4 sub-steps)
                for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
                    calc_rhs();
                    rk4_step(u, uh, rhs, dt, rk_step);
                }
            }

            void print_sim(double t)
            {
                // Create mdspan view for printing
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);


                fmt::println("t = {:>-7.3f}:", t);

                for (size_t n = 0; n < npts; ++n) {

                    // get current state (= current positional transformation bivector)
                    vec2dp B = u[0];     // B = Omega * t + B0 (from integration)
                    vec2dp Omega = u[1]; // dB/dt = Omega = dB^2/dt^2 * t + Omega0

                    // calculate current position from B via M = rexp(0.5 * B) ⟇ M0
                    // and via pts(t) = M ⟇ pts(t0) ⟇ rrev(M) = move2dp(pts(t0),M)
                    // and d(pts(t))/dt = rcmt(Omega, pts(t))

                    auto M = rexp(0.5 * B);
                    auto pt = move2dp(pts[n], M);
                    // fmt::println(
                    //     "    n = {}, B = {:>-7.3f}, Omega = {:>-7.3f}, M = {:>-7.3f}",
                    //     n, B, Omega, M);
                    fmt::println("    n = {}, pos. = {:>-7.3f}, vel. = {:>-7.3f}", n, pt,
                                 rcmt(Omega, pt));
                }
            }

          private:

            // number of independent points to solve system for
            size_t npts;

            // positions of points in the body frame (fixed for now to get started)
            std::vector<vec2dp> pts; // point positions (fixed for now)
            std::vector<value_t> m;  // mass of points

            // initial position and speed of center of mass relative to world system
            vec2dp cm_w_pos0; // initial position of center of mass (unitized)
            vec2dp cm_w_spd0; // initial speed of center of mass (.z == 0)

            // initial angular position of cm (body vs. world frame, e12 defines pos phi)
            [[maybe_unused]] value_t cm_w_phi0;

            // initial angular velocity of cm (body vs. world frame)
            value_t cm_w_omega0;

            // RK4 integration state for point n with system order = 2
            // => [n+0: position, n+1: velocity])
            std::vector<vec2dp> u_mem;   // [n+0: position, n+1: velocity]
            std::vector<vec2dp> uh_mem;  // helper for integration
            std::vector<vec2dp> rhs_mem; // right-hand side values

            // inertia map and its inverse (calculated from descrete input values)
            // calculation is done in body frame
            inertia2dp I;
            inertia2dp I_inv;
        };

        auto constexpr num_points = 2;
        // (initial) positions and masses of each mass point in the body frame
        std::vector<vec2dp> pos_b = {vec2dp{0, 1, 1}, vec2dp{0, -1, 1}}; // position
        std::vector<value_t> m = {1.0, 1.0};                             // mass

        // initial position and movement of body frame vs. world frame (in world frame)
        auto cm_pos = O_2dp;           // world and body frame coincide at t=0
        auto cm_spd = vec2dp{0, 0, 0}; // body frame is moving with v0 vs. world frame

        // initial rotation and anguar velocity of body frame vs. world frame
        auto cm_phi = 0.0;    // no initial rotation
        auto cm_omega = -1.0; // initial angular velocity (T = 2*pi/omega = 6.2832s)
                              // turns clockwise, i.e. againt direction of e12

        sim_ode_rigid_pga2dp sim(num_points, pos_b, m, cm_pos, cm_spd, cm_phi, cm_omega);

        // time range (from, to, number of steps)
        // time range (1 full revolution)
        auto t_rng = discrete_range(0.0, 6.2832, 60);

        sim.set_initial_values();
        sim.print_sim(t_rng.min());

        for (size_t n = 1; n <= t_rng.steps(); n++) {

            // integration from t to t + dt
            double t = t_rng.min() + n * t_rng.delta();
            sim.calc_rkstep(t_rng.delta());

            // print sim status at t+dt
            sim.print_sim(t);
        }

        fmt::println("");
    }


    ///////////////////////////////////////////////////////////////////////
    // Kinematic frame-transformation test cases
    //
    // Purpose: demonstrate and verify the body-frame / world-frame
    // transformation framework using purely kinematic setups (constant
    // velocity fields, no forces, no ODE integration).
    //
    // Notation (PGA2DP):
    //   B_b   - accumulated position generator in body frame (integrated)
    //   B_w   - accumulated position generator in world frame (integrated)
    //   Omega_b - velocity generator in body frame (constant here)
    //   Omega_w - velocity generator in world frame (= move2dp(Omega_b, M))
    //   M0    - initial motor: maps body frame to world frame at t=0
    //   M(t)  - current motor: M0 ⟇ rexp(½ B_b(t))  [body-frame formulation]
    //   P_b   - a fixed point in the body frame
    //   P_w   - P_b expressed in world frame: move2dp(P_b, M)
    //
    // Key identity that reveals the bug in the original code:
    //   move2dp(B, rexp(½B)) == B   always (B commutes with its own exponential)
    //   => B_w differs from B_b only when M0 != ps (non-trivial initial state)
    //
    // Extensibility: body-frame velocity vectors v_b = vec2dp{vx, vy, 0} (z=0)
    //   transform to world via the same sandwich:  v_w = move2dp(v_b, M)
    ///////////////////////////////////////////////////////////////////////

    // Helper: print one time-step for kinematic cases
    auto print_kine_step = [](double t, mvec2dp_u const& /*M0*/, vec2dp const& B_b,
                              vec2dp const& B_w, mvec2dp_u const& /*M*/,
                              vec2dp const& P_b, vec2dp const& P_w, vec2dp const& /*Q_b*/,
                              vec2dp const& Q_w_check) {
        fmt::println("  t={:>6.3f}: B_b=({:>7.3f},{:>7.3f},{:>7.3f})"
                     " B_w=({:>7.3f},{:>7.3f},{:>7.3f})"
                     " P_b->P_w=({:>7.3f},{:>7.3f})->({:>7.3f},{:>7.3f})"
                     " Q_w_check=({:>7.3f},{:>7.3f})",
                     t, B_b.x, B_b.y, B_b.z, B_w.x, B_w.y, B_w.z, P_b.x, P_b.y, P_w.x,
                     P_w.y, Q_w_check.x, Q_w_check.y);
    };

    TEST_CASE("pga2dp: kinematics - pure translation")
    {
        fmt::println("pga2dp: kinematics - pure translation");
        fmt::println("  Body moves with velocity v=(+1,0) in world frame.");
        fmt::println("  Omega_b = (-vy, vx, 0) = (0, 1, 0)  [translation encoding]");

        // body-frame reference point (e.g. body origin at body-frame (0,0,1))
        vec2dp const P_b{0.0, 0.0, 1.0}; // body origin
        // body-frame "pivot" anchor (z=0 means ideal point / pure direction)
        vec2dp const Q_b{0.0, 0.0, 0.0}; // unused for translation; use O_b

        // constant velocity generator in body frame:
        // translation in +x: encoding = (-vy, vx, 0) with vx=1, vy=0
        vec2dp const Omega_b{0.0, 1.0, 0.0};

        double const T = 3.0;
        double const dt = T / 4.0;
        int const N = 4;

        fmt::println("  --- Case A: M0 = identity (body starts at world origin) ---");
        {
            mvec2dp_u const M0 = rexp(0.5 * vec2dp{0.0, 0.0, 0.0}); // ps (identity)
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, rexp(0.5 * B_b));
                vec2dp const Omega_w = move2dp(Omega_b, M);
                vec2dp const P_w = move2dp(P_b, M);
                // Q_w_check: body origin in world (should equal P_w here)
                vec2dp const Q_w = move2dp(P_b, M);
                print_kine_step(t, M0, B_b, B_w, M, P_b, P_w, Q_b, Q_w);

                if (i < N) {
                    B_b = B_b + Omega_b * dt;
                    B_w = B_w + Omega_w * dt;
                }
            }
        }

        fmt::println("  --- Case B: M0 = translation to (2,1) (body starts offset) ---");
        {
            // M0: pure translation motor placing body origin at world (2, 1)
            // Translation encoding: P_tra = (-ty, tx, 0) = (-1, 2, 0)
            vec2dp const cm0{2.0, 1.0, 1.0};
            mvec2dp_u const M0 = rexp(0.5 * vec2dp{-cm0.y, cm0.x, 0.0});
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, rexp(0.5 * B_b));
                vec2dp const Omega_w = move2dp(Omega_b, M);
                vec2dp const P_w = move2dp(P_b, M);
                vec2dp const Q_w = P_w; // body origin in world
                print_kine_step(t, M0, B_b, B_w, M, P_b, P_w, Q_b, Q_w);

                if (i < N) {
                    B_b = B_b + Omega_b * dt;
                    B_w = B_w + Omega_w * dt;
                }
            }
        }
        fmt::println(
            "  Note: B_w == B_b always for pure translation (regardless of M0).");
        fmt::println(
            "        Translation velocity (z=0) is a direction vector: invariant under");
        fmt::println(
            "        any motor, so Omega_w == Omega_b and B_w accumulates identically.");
        fmt::println("");
    }

    TEST_CASE("pga2dp: kinematics - pure rotation")
    {
        fmt::println("pga2dp: kinematics - pure rotation");
        fmt::println(
            "  Body rotates at omega=1 rad/s about body-frame pivot Q_b=(1,0,1).");

        // body-frame pivot (e.g. the TR corner of a plate at (1,0) from body origin)
        vec2dp const Q_b{1.0, 0.0, 1.0};
        // a different body-frame reference point to track
        vec2dp const P_b{-1.0, 0.5, 1.0};

        double const omega = 1.0; // rad/s
        // constant velocity generator in body frame: rotation about Q_b
        vec2dp const Omega_b = omega * Q_b;

        double const T = 2.0 * pi; // one full revolution
        double const dt = T / 8.0;
        int const N = 8;

        fmt::println("  --- Case A: M0 = identity (Q_b is at world origin) ---");
        {
            mvec2dp_u const M0 = rexp(0.5 * vec2dp{0.0, 0.0, 0.0});
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};
            // pivot world position (should stay fixed)
            [[maybe_unused]] vec2dp const Q_world_expected = move2dp(Q_b, M0);

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, rexp(0.5 * B_b));
                vec2dp const Omega_w = move2dp(Omega_b, M);
                vec2dp const P_w = move2dp(P_b, M);
                vec2dp const Q_w_chk = move2dp(Q_b, M); // pivot must stay fixed
                print_kine_step(t, M0, B_b, B_w, M, P_b, P_w, Q_b, Q_w_chk);

                if (i < N) {
                    B_b = B_b + Omega_b * dt;
                    B_w = B_w + Omega_w * dt;
                }
            }
        }

        fmt::println("  --- Case B: M0 = translation placing Q_b at world (2,1) ---");
        {
            // We want Q_b to be at world position (2, 1).
            // M0 is a pure translation by (2-Q_b.x, 1-Q_b.y) = (1, 1).
            // Translation encoding: P_tra = (-ty, tx, 0) = (-1, 1, 0)
            vec2dp const Q_world{2.0, 1.0, 1.0};
            vec2dp const trans{Q_world.x - Q_b.x, Q_world.y - Q_b.y, 1.0};
            mvec2dp_u const M0 = rexp(0.5 * vec2dp{-trans.y, trans.x, 0.0});
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, rexp(0.5 * B_b));
                vec2dp const Omega_w = move2dp(Omega_b, M);
                vec2dp const P_w = move2dp(P_b, M);
                vec2dp const Q_w_chk = move2dp(Q_b, M); // pivot must stay at (2,1)
                print_kine_step(t, M0, B_b, B_w, M, P_b, P_w, Q_b, Q_w_chk);

                if (i < N) {
                    B_b = B_b + Omega_b * dt;
                    B_w = B_w + Omega_w * dt;
                }
            }
        }
        fmt::println("  Note: Q_w_check should remain constant (pivot is fixed).");
        fmt::println(
            "        B_w != B_b when M0 != identity: world-frame encoding differs.");
        fmt::println(
            "        B_b.z == B_w.z always: rotation angle is frame-independent in 2D.");
        fmt::println("");
    }

    TEST_CASE("pga2dp: kinematics - combined translation and rotation")
    {
        fmt::println("pga2dp: kinematics - combined translation and rotation");
        fmt::println(
            "  Body rotates about Q_b=(1,0.5,1) and simultaneously translates (+x).");

        vec2dp const Q_b{1.0, 0.5, 1.0};  // body-frame pivot
        vec2dp const P_b{-1.0, 0.0, 1.0}; // body-frame reference point

        double const omega = 1.0; // rad/s
        // combined: rotation about Q_b + translation in +x
        // translation encoding: (-vy, vx, 0) = (0, 0.5, 0) for vx=0.5, vy=0
        vec2dp const Omega_b = omega * Q_b + vec2dp{0.0, 0.5, 0.0};

        double const T = 4.0;
        double const dt = T / 8.0;
        int const N = 8;

        fmt::println("  --- Case A: M0 = identity ---");
        {
            mvec2dp_u const M0 = rexp(0.5 * vec2dp{0.0, 0.0, 0.0});
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, rexp(0.5 * B_b));
                vec2dp const Omega_w = move2dp(Omega_b, M);
                vec2dp const P_w = move2dp(P_b, M);
                vec2dp const Q_w_chk = move2dp(Q_b, M);
                print_kine_step(t, M0, B_b, B_w, M, P_b, P_w, Q_b, Q_w_chk);

                if (i < N) {
                    B_b = B_b + Omega_b * dt;
                    B_w = B_w + Omega_w * dt;
                }
            }
        }

        fmt::println("  --- Case B: M0 = translation to (3,2), initial rotation "
                     "phi0=pi/4 in B_b ---");
        {
            // M0 encodes the REST position (phi=0) of the body in world frame.
            // The initial rotation phi0 is encoded in B_b(0) = phi0 * Q_b.
            // By the conjugation identity: M0_tra ⟇ rexp(½ phi0 Q_b)
            //   == (rotation by phi0 about Q_world) ⟇ M0_tra
            // so this correctly places the plate at the phi0-rotated rest position.
            vec2dp const cm0{3.0, 2.0, 1.0}; // rest cm world position (phi=0)
            mvec2dp_u const M0 =
                rexp(0.5 * vec2dp{-cm0.y, cm0.x, 0.0}); // pure translation

            double const phi0 = pi / 4.0;
            vec2dp B_b = phi0 * Q_b; // initial B_b: phi0 rotation about body-frame Q_b
            vec2dp B_w = phi0 * Q_b; // B_w starts equal to B_b (same reference pose)

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, rexp(0.5 * B_b));
                vec2dp const Omega_w = move2dp(Omega_b, M);
                vec2dp const P_w = move2dp(P_b, M);
                vec2dp const Q_w_chk = move2dp(Q_b, M);
                print_kine_step(t, M0, B_b, B_w, M, P_b, P_w, Q_b, Q_w_chk);

                if (i < N) {
                    B_b = B_b + Omega_b * dt;
                    B_w = B_w + Omega_w * dt;
                }
            }
        }
        fmt::println(
            "  Note: Omega_b is constant; Omega_w varies (body frame rotates in world).");
        fmt::println(
            "        B_w and B_b accumulate differently because Omega_w(t) != Omega_b.");
        fmt::println("");
    }

    TEST_CASE("pga2dp: combined motion (rigid 1-body system)")
    {
        fmt::println("pga2dp: combined motion (rigid 1-body system)");

        class sim_ode_plate_pga2dp { // model 2nd order ode by a 1st order system


          public:

            sim_ode_plate_pga2dp(value_t m_in, value_t w_in, value_t h_in,
                                 vec2dp const& cm_pos_in, vec2dp const& cm_spd_in,
                                 value_t cm_phi_in, value_t cm_omega_in) :
                m(m_in), width(w_in), height(h_in), cm_w_pos0(cm_pos_in),
                cm_w_spd0(cm_spd_in), cm_w_phi0(cm_phi_in), cm_w_omega0(cm_omega_in),
                u_mem(2), uh_mem(2 * 2), rhs_mem(2)
            {

                fmt::println("sim_ode_plate_pga2dp: combined motion.");

                // Body frame origin = CENTER OF MASS (plate is symmetric about cm).
                // Pivot = TR corner Q_b = (hw, hh, 1) in body frame.
                // get_plate_inertia(m, w, h, Q_b) applies the parallel-axis correction
                // (Steiner) so I is directly about Q_b — no manual scaling needed in
                // calc_rhs. And the total kinetic energy is K_rot
                value_t const hw = width / 2.0;
                value_t const hh = height / 2.0;
                vec2dp const Q_b{hw, hh, 1.0};
                fmt::println(
                    "inertia of plate in body system (about pivot Q_b = TR corner):");
                I = get_plate_inertia(m, width, height, Q_b);
                // I = get_plate_inertia(m, width, height);

                // get inverse inertia map (maps momentum bivector to velocity vector)
                I_inv = get_inertia_inverse(I);

                fmt::println("system inertia (about pivot Q_b):");
                fmt::println("I     = {:>-7.3f}", I);
                fmt::println("I_inv = {:>-7.3f}", I_inv);
            }

            void set_initial_values()
            {
                // Create mdspan view for setting initial values
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);

                value_t const hw = width / 2.0;
                value_t const hh = height / 2.0;

                // Body origin = CM. Pivot = TR corner = Q_b = (hw, hh, 1) in body frame.
                // M0: pure translation placing body origin (= cm) at cm_w_pos0.
                // Encoding: translation by (tx, ty) -> motor rexp(0.5*(-ty, tx, 0))
                M0 = rexp(0.5 * vec2dp{-cm_w_pos0.y, cm_w_pos0.x, 0.0});

                // pivot_w: world position of body-frame pivot Q_b via M0 (stays fixed)
                vec2dp const Q_b{hw, hh, 1.0};
                pivot_w = move2dp(Q_b, M0);
                fmt::println("pivot_w = {:>-7.3f}  (world position of pivot Q_b via M0)",
                             pivot_w);

                // B_b(0): initial rotation by phi0 about body-frame pivot Q_b.
                // move2dp(Q_b, rexp(½*phi*Q_b)) = Q_b for all phi (pivot invariant).
                u[0] = cm_w_phi0 * Q_b;

                // Omega_b(0): initial angular velocity about pivot Q_b.
                // For a pivot-constrained body, translational velocity of cm is
                // determined by omega (not independent): Omega_b = omega * Q_b.
                u[1] = cm_w_omega0 * Q_b;
            }

            void calc_rhs()
            {

                // Create mdspan views for rhs calculation
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto rhs = mdspan<vec2dp, dextents<size_t, 1>>(rhs_mem.data(), 2);

                // get current state
                vec2dp const B = u[0];     // position transformation B is in u[0]
                vec2dp const Omega = u[1]; // velocity trafo d(B)/dt = Omega is in u[1]

                // current motor: M(t) = M0 ⟇ rexp(½ B_b(t))  [body-frame formulation]
                auto const M = rgpr(M0, rexp(0.5 * B));

                // cm world position: body origin = cm = (0,0,1) in body frame
                vec2dp const O_b{0.0, 0.0, 1.0};
                auto const cm_w = move2dp(O_b, M);

                // fmt::println("    cm_w = {:>-7.3f}, pivot_w = {:>-7.3f}", cm_w,
                // pivot_w);

                auto const g = vec2dp{0.0, -9.81, 0.0}; // 9.81 m/s^2 in neg. y-direction

                // x-direction: Force couple in world frame with
                //              x-component resulting from splitting y-acceleration
                //              into centrifugal acceleration a_cf = r * omega^2 and
                auto const omega = Omega.z; // angular speed in z-component
                auto const a_cf = (cm_w - pivot_w) * omega * omega; // centrifugal accel.
                auto const f_cm = m * (a_cf + g);                   // total force at cm

                // fmt::println("    a_cf = {:>-7.3f}, f_cm = {:>-7.3f}", a_cf, f_cm);

                auto const F_cm = wdg(cm_w, f_cm);
                auto const F_piv = wdg(pivot_w, -f_cm);
                auto const F_b = move2dp(F_cm + F_piv, rrev(M)); // res. in body frame

                // fmt::println("    F_b = {:>-7.3f}", F_b);

                // Enforce pivot constraint explicitly: extract scalar alpha from
                // dOmega.z and project back onto Q_b so that Omega_b = omega*Q_b
                // is preserved exactly, preventing numerical drift of off-constraint
                // components (e31/e32) during RK4 integration.
                // In principle correct without projection: rhs[1] =
                // compute_omega_dot(I_inv, F_b, Omega, I);
                value_t const hw_r = width / 2.0;
                value_t const hh_r = height / 2.0;
                vec2dp const Q_b_r{hw_r, hh_r, 1.0};
                auto const dOmega = compute_omega_dot(I_inv, F_b, Omega, I);
                value_t const alpha = dOmega.z; // = tau / I_pivot[2,2]

                rhs[0] = Omega;         // dB_b/dt = Omega_b  (= omega * Q_b)
                rhs[1] = alpha * Q_b_r; // dOmega_b/dt = alpha * Q_b  (pivot constraint)
            }

            void calc_rkstep(double dt)
            {
                // Create mdspan views for RK4 integration (for all npts points)
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto uh = mdspan<vec2dp, dextents<size_t, 2>>(uh_mem.data(), 2, 2);
                auto rhs = mdspan<vec2dp const, dextents<size_t, 1>>(rhs_mem.data(), 2);

                // Perform RK4 integration (4 sub-steps)
                for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
                    calc_rhs();
                    rk4_step(u, uh, rhs, dt, rk_step);
                }
            }

            void print_sim(double t)
            {
                // Create mdspan view for printing
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);


                fmt::println("t = {:>-7.3f}:", t);

                // get current state (= current positional transformation bivector)

                vec2dp const B = u[0]; // B_b = phi * Q_b (accumulated rotation)

                // current motor: M(t) = M0 ⟇ rexp(½ B_b(t))  [body-frame formulation]
                auto const M = rgpr(M0, rexp(0.5 * B));

                // B_b = phi * Q_b, so B.z = phi (rotation angle about Q_b)
                value_t const phi_b = B.z;

                // cm world position: body origin = cm = (0,0,1) in body frame
                vec2dp const O_b{0.0, 0.0, 1.0};
                auto const cm_w = move2dp(O_b, M);

                // B_w: B_b expressed in world frame via M0 (not full M).
                // For pure couple torque, move2dp(B,M) = B (commutator identity).
                // move2dp(B, M0) gives the world-frame encoding when M0 != ps.
                auto const B_w = move2dp(B, M0);
                value_t const phi_w = B_w.z; // rotation angle is frame-invariant

                fmt::println("    cm_w = {:>-6.3f}, phi_b = {:>-6.3f}, phi_w = {:>-6.3f},"
                             " B_w = {:>-6.3f}, B = {:>-6.3f}",
                             cm_w, phi_b, phi_w, B_w, B);
            }

            // Accessors for validation checks
            vec2dp get_cm_world() const
            {
                auto u = mdspan<vec2dp const, dextents<size_t, 1>>(u_mem.data(), 2);
                auto const M = rgpr(M0, rexp(0.5 * u[0]));
                return move2dp(vec2dp{0.0, 0.0, 1.0}, M);
            }

            value_t get_omega() const
            {
                // Omega = omega * Q_b,  Q_b.z = 1  =>  omega = Omega.z
                auto u = mdspan<vec2dp const, dextents<size_t, 1>>(u_mem.data(), 2);
                return u[1].z;
            }

            value_t get_I_zz_pivot() const { return I.view()[2, 2]; }

          private:


            value_t m;      // mass of plate
            value_t width;  // width of plate (e1-direction in body frame), sym. to O_b
            value_t height; // height of plate (e2-direction in body frame), sym. to O_b

            // initial position and speed of center of mass relative to world system
            vec2dp cm_w_pos0; // initial position of center of mass (unitized)
            vec2dp cm_w_spd0; // initial speed of center of mass (.z == 0)

            // initial angular position of cm (body vs. world frame, e12 defines pos phi)
            value_t cm_w_phi0;

            // initial angular velocity of cm (body vs. world frame)
            value_t cm_w_omega0;

            mvec2dp_u M0;   // initial motor: maps body frame to world at t=0
            vec2dp pivot_w; // world position of body-frame TR corner (fixed pivot)

            // RK4 integration state for point n with system order = 2
            // => [n+0: position, n+1: velocity])
            std::vector<vec2dp> u_mem;   // [n+0: position, n+1: velocity]
            std::vector<vec2dp> uh_mem;  // helper for integration
            std::vector<vec2dp> rhs_mem; // right-hand side values

            // inertia map about pivot Q_b (parallel-axis corrected via get_plate_inertia)
            inertia2dp I;
            inertia2dp I_inv;
        };

        value_t m = 1.0; // mass
        value_t w = 2.0; // width
        value_t h = 2.0; // height

        // initial position and movement of body frame vs. world frame (in world frame)
        auto cm_pos = vec2dp{-1, -1, 1}; // world and body frame coincide at t=0
        auto cm_spd = vec2dp{0, 0, 0};   // body frame is moving with v0 vs. world frame

        // initial rotation and anguar velocity of body frame vs. world frame
        auto cm_phi = 0.0;   // no initial rotation
        auto cm_omega = 0.0; // initial angular velocity

        sim_ode_plate_pga2dp sim(m, w, h, cm_pos, cm_spd, cm_phi, cm_omega);

        // time range (from, to, number of steps)
        // time range (1 full revolution)
        auto t_rng = discrete_range(0.0, 2.0, 60);

        sim.set_initial_values();
        sim.print_sim(t_rng.min());

        // --- Static validation: inertia analytical check -----------------------
        // I[2,2] = m*(w^2+h^2)/12 + m*(hw^2+hh^2)
        //        = 1*(4+4)/12 + 1*(1+1) = 2/3 + 2 = 8/3
        value_t const hw = w / 2.0;
        value_t const hh = h / 2.0;
        CHECK(sim.get_I_zz_pivot() == doctest::Approx(8.0 / 3.0));

        // --- Dynamic validation setup ------------------------------------------
        // Initial energy: KE=0, PE = m*g*cm_w.y = 9.81*(-1) = -9.81 J
        value_t const g = 9.81;
        value_t const E_0 = m * g * cm_pos.y; // = -9.81 J

        // Pivot at world origin (0,0) => pivot distance^2 = hw^2 + hh^2 = 2.0
        value_t const pivot_dist_sq = hw * hw + hh * hh; // = 2.0

        for (size_t n = 1; n <= t_rng.steps(); n++) {

            // integration from t to t + dt
            double t = t_rng.min() + n * t_rng.delta();
            sim.calc_rkstep(t_rng.delta());

            // print sim status at t+dt
            sim.print_sim(t);

            vec2dp const cm_w = sim.get_cm_world();
            value_t const omega = sim.get_omega();

            // Energy conservation: T + V = E_0.
            // T = 0.5 * I_pivot * omega^2 (rotation about fixed pivot, parallel-axis I)
            // V = m * g * cm_w.y  (gravitational PE, y=0 datum)
            value_t const T_kin = 0.5 * sim.get_I_zz_pivot() * omega * omega;
            value_t const V_pot = m * g * cm_w.y;
            CHECK((T_kin + V_pot) == doctest::Approx(E_0).epsilon(1e-3));

            // Pivot constraint: distance from pivot (0,0) to cm = rsqrt(hw^2+hh^2)
            value_t const r_sq = cm_w.x * cm_w.x + cm_w.y * cm_w.y;
            CHECK(r_sq == doctest::Approx(pivot_dist_sq).epsilon(1e-3));
        }

        fmt::println("");
    }


    TEST_CASE("pga2dp: multibody system transformation - simple translation")
    {
        fmt::println("pga2dp: multibody system transformation - simple translation)");

        static_frame2dp rf;
        fmt::println("parent_frame = {:.3f}", rf); // nested spec forwarding

        auto phi = deg2rad(0);         // rotation vs. parent
        auto origin = vec2dp{2, 1, 1}; // origin in parent coordinates
        static_frame2dp cf("C"s, origin, phi);
        fmt::println("child_frame  = {:.3f}", cf);

        fmt::println("");
        fmt::println("transformation from parent to child (rf -> cf)");

        static_system2dp ds;
        ds.add_frame(rf);
        ds.add_frame(cf);
        fmt::println("ds = {:.3f}", ds);

        auto p0_rf = vec2dp(3, 2, 1);
        auto M = ds.get_pos_trafo(0, 1); // trafo from (idx: 0) to cf (idx: 1)

        auto p0_cf = move2dp(p0_rf, M);
        auto e1_rf = e1_2dp;
        auto e1_cf = move2dp(e1_rf, M);

        fmt::println("M = {:.3f}", M);
        fmt::println("p0_rf = {:.3f}", p0_rf);
        fmt::println("p0_cf = {:.3f}", p0_cf);
        fmt::println("e1_rf = {:.3f}", e1_rf);
        fmt::println("e1_cf = {:.3f}", e1_cf);

        CHECK(p0_cf == vec2dp(1, 1, 1));
        CHECK(e1_cf == e1_rf); // for simple cs translation, a vector is not transformed
                               // (a point at infinity remains at infinity under transl.)

        fmt::println("");
        fmt::println("transformation from child to parent (cf -> rf)");

        auto Mr = ds.get_pos_trafo(1, 0); // trafo from (idx: 1) to cf (idx: 0)

        auto p0_rf_t = move2dp(p0_cf, Mr);
        auto e1_rf_t = move2dp(e1_cf, Mr);

        fmt::println("Mr = {:.3f}", Mr);
        fmt::println("p0_rf_t = {:.3f}", p0_rf_t);
        fmt::println("e1_rf_t = {:.3f}", e1_rf_t);

        CHECK(p0_rf_t == p0_rf);
        CHECK(e1_rf_t == e1_rf);

        fmt::println("");
    }

    TEST_CASE("pga2dp: multibody system transformation - simple rotation")
    {
        fmt::println("pga2dp: multibody system transformation - simple rotation)");

        static_frame2dp rf;
        fmt::println("parent_frame = {:.3f}", rf); // nested spec forwarding

        auto phi = deg2rad(30);        // rotation vs. parent
        auto origin = vec2dp{0, 0, 1}; // origin in parent coordinates
        static_frame2dp cf("C"s, origin, phi);
        fmt::println("child_frame  = {:.3f}", cf);

        fmt::println("");
        fmt::println("transformation from parent to child (rf -> cf)");

        static_system2dp ds;
        ds.add_frame(rf);
        ds.add_frame(cf);
        fmt::println("ds = {:.3f}", ds);

        auto p0_rf = vec2dp(3, 2, 1);
        auto M = ds.get_pos_trafo(0, 1); // trafo from (idx: 0) to cf (idx: 1)

        auto p0_cf = move2dp(p0_rf, M);
        auto e1_rf = e1_2dp;
        auto e1_cf = move2dp(e1_rf, M);

        fmt::println("M = {:.3f}", M);
        fmt::println("p0_rf = {:.3f}", p0_rf);
        fmt::println("p0_cf = {:.3f}", p0_cf);
        fmt::println("e1_rf = {:.3f}", e1_rf);
        fmt::println("e1_cf = {:.3f}", e1_cf);

        // for a simple rotation point and vector coordinates are transformed
        auto l = to_val(bulk_nrm(p0_rf));
        auto rphi_cf = angle(p0_rf, e1_rf) - phi; // relative angle towards p0 in cf

        auto constexpr cmp_eps = 1e-10;
        CHECK(rphi_cf == doctest::Approx(angle(p0_cf, e1_2dp))
                             .epsilon(cmp_eps)); // rel. angle == angle vs. e1 in cf

        fmt::println("relative angle in cf = {:.3f}", rad2deg(rphi_cf));

        // component-wise comparison with selectable accuracy (the default eps of
        // vec2dp::operator== is too strict after the wedge/division in move2dp)

        CHECK(p0_cf.x == doctest::Approx(l * std::cos(rphi_cf)).epsilon(cmp_eps));
        CHECK(p0_cf.y == doctest::Approx(l * std::sin(rphi_cf)).epsilon(cmp_eps));
        CHECK(p0_cf.z == doctest::Approx(1.0).epsilon(cmp_eps));
        CHECK(e1_cf == vec2dp{std::cos(phi), -std::sin(phi), 0.0});


        fmt::println("");
        fmt::println("transformation from child to parent (cf -> rf)");

        auto Mr = ds.get_pos_trafo(1, 0); // trafo from (idx: 1) to cf (idx: 0)

        auto p0_rf_t = move2dp(p0_cf, Mr);
        auto e1_rf_t = move2dp(e1_cf, Mr);

        fmt::println("Mr = {:.3f}", Mr);
        fmt::println("p0_rf_t = {:.3f}", p0_rf_t);
        fmt::println("e1_rf_t = {:.3f}", e1_rf_t);

        CHECK(p0_rf_t == p0_rf);
        CHECK(e1_rf_t == e1_rf);

        fmt::println("");
    }

    TEST_CASE("pga2dp: multibody system transformation - combined translation/rotation")
    {
        fmt::println(
            "pga2dp: multibody system transformation - combined translation/rotation)");

        static_frame2dp rf;
        fmt::println("parent_frame = {:.3f}", rf); // nested spec forwarding

        auto phi = deg2rad(8);           // rotation vs. parent
        auto origin = vec2dp{2, 1.5, 1}; // origin in parent coordinates
        static_frame2dp cf("C"s, origin, phi);
        fmt::println("child_frame  = {:.3f}", cf);

        fmt::println("");
        fmt::println("transformation from parent to child (rf -> cf)");

        static_system2dp ds;
        ds.add_frame(rf);
        ds.add_frame(cf);
        fmt::println("ds = {:.3f}", ds);

        auto p0_rf = vec2dp(4, 5, 1);
        auto M = ds.get_pos_trafo(0, 1); // trafo from (idx: 0) to cf (idx: 1)

        auto p0_cf = move2dp(p0_rf, M);
        auto e1_rf = e1_2dp;
        auto e1_cf = move2dp(e1_rf, M);

        fmt::println("M = {:.3f}", M);
        fmt::println("p0_rf = {:.3f}", p0_rf);
        fmt::println("p0_cf = {:.3f}", p0_cf);
        fmt::println("e1_rf = {:.3f}", e1_rf);
        fmt::println("e1_cf = {:.3f}", e1_cf);

        auto rl = to_val(bulk_nrm(p0_rf - origin));
        auto rphi_cf =
            angle(p0_rf - origin, e1_rf) - phi; // relative angle towards p0 in cf

        fmt::println("relative angle in cf = {:.3f}", rad2deg(rphi_cf));

        // component-wise comparison with selectable accuracy (the default eps of
        // vec2dp::operator== is too strict after the wedge/division in move2dp)
        auto constexpr cmp_eps = 1e-10;
        CHECK(p0_cf.x == doctest::Approx(rl * std::cos(rphi_cf)).epsilon(cmp_eps));
        CHECK(p0_cf.y == doctest::Approx(rl * std::sin(rphi_cf)).epsilon(cmp_eps));
        CHECK(p0_cf.z == doctest::Approx(1.0).epsilon(cmp_eps));
        CHECK(e1_cf == vec2dp{std::cos(phi), -std::sin(phi), 0.0});


        fmt::println("");
        fmt::println("transformation from child to parent (cf -> rf)");

        auto Mr = ds.get_pos_trafo(1, 0); // trafo from (idx: 1) to cf (idx: 0)

        auto p0_rf_t = move2dp(p0_cf, Mr);
        auto e1_rf_t = move2dp(e1_cf, Mr);

        fmt::println("Mr = {:.3f}", Mr);
        fmt::println("p0_rf_t = {:.3f}", p0_rf_t);
        fmt::println("e1_rf_t = {:.3f}", e1_rf_t);

        // round-trip: component-wise with selectable accuracy (strict operator== eps is
        // too tight after the forward+backward move2dp accumulates rounding error)
        CHECK(p0_rf_t.x == doctest::Approx(p0_rf.x).epsilon(cmp_eps));
        CHECK(p0_rf_t.y == doctest::Approx(p0_rf.y).epsilon(cmp_eps));
        CHECK(p0_rf_t.z == doctest::Approx(p0_rf.z).epsilon(cmp_eps));
        CHECK(e1_rf_t.x == doctest::Approx(e1_rf.x).epsilon(cmp_eps));
        CHECK(e1_rf_t.y == doctest::Approx(e1_rf.y).epsilon(cmp_eps));
        CHECK(e1_rf_t.z == doctest::Approx(e1_rf.z).epsilon(cmp_eps));

        fmt::println("");
    }

    TEST_CASE("pga2dp: multibody system transformation - chain of three systems")
    {
        fmt::println("pga2dp: multibody system transformation - chain of three systems)");

        // Chain of three frames; each child is placed RELATIVE to its parent (the
        // previously added frame), forming the dependency chain rf -> cf1 -> cf2.
        static_frame2dp rf; // inertial root frame (identity)
        fmt::println("parent_frame = {:.3f}", rf);

        auto phi1 = deg2rad(8);           // cf1 rotation vs. rf
        auto origin1 = vec2dp{2, 1.5, 1}; // cf1 origin in rf coordinates
        static_frame2dp cf1("C1"s, origin1, phi1);
        fmt::println("child_frame1 = {:.3f}", cf1);

        auto phi2 = deg2rad(12);           // cf2 rotation vs. cf1
        auto origin2 = vec2dp{1, -0.5, 1}; // cf2 origin in cf1 coordinates
        static_frame2dp cf2("C2"s, origin2, phi2);
        fmt::println("child_frame2 = {:.3f}", cf2);

        static_system2dp ds;
        ds.add_frame(rf);
        ds.add_frame(cf1);
        ds.add_frame(cf2);
        fmt::println("ds = {:.3f}", ds);

        auto p0_rf = vec2dp(4, 5, 1);

        // reference: passive coordinate change one frame deeper into the chain,
        // p_child = R(-phi) * (p_parent - origin_rel)
        auto to_child = [](vec2dp const& p, vec2dp const& o, value_t phi) {
            auto const c = std::cos(phi);
            auto const s = std::sin(phi);
            auto const dx = p.x - o.x;
            auto const dy = p.y - o.y;
            return vec2dp(c * dx + s * dy, -s * dx + c * dy, 1.0);
        };

        // forward transforms via the chain
        auto M01 = ds.get_pos_trafo(0, 1); // rf  -> cf1
        auto M02 = ds.get_pos_trafo(0, 2); // rf  -> cf2 (composed over cf1)
        auto M12 = ds.get_pos_trafo(1, 2); // cf1 -> cf2

        auto p0_cf1 = move2dp(p0_rf, M01);
        auto p0_cf2 = move2dp(p0_rf, M02);

        auto e1_rf = e1_2dp;
        auto e1_cf2 = move2dp(e1_rf, M02);

        fmt::println("M01 = {:.3f}", M01);
        fmt::println("M02 = {:.3f}", M02);
        fmt::println("p0_rf  = {:.3f}", p0_rf);
        fmt::println("p0_cf1 = {:.3f}", p0_cf1);
        fmt::println("p0_cf2 = {:.3f}", p0_cf2);
        fmt::println("e1_cf2 = {:.3f}", e1_cf2);

        auto constexpr cmp_eps = 1e-10;

        // link rf -> cf1 reproduces the standalone two-frame result
        auto ref_cf1 = to_child(p0_rf, origin1, phi1);
        CHECK(p0_cf1.x == doctest::Approx(ref_cf1.x).epsilon(cmp_eps));
        CHECK(p0_cf1.y == doctest::Approx(ref_cf1.y).epsilon(cmp_eps));
        CHECK(p0_cf1.z == doctest::Approx(ref_cf1.z).epsilon(cmp_eps));

        // composed rf -> cf2 equals chaining the two passive coordinate changes
        auto ref_cf2 = to_child(to_child(p0_rf, origin1, phi1), origin2, phi2);
        CHECK(p0_cf2.x == doctest::Approx(ref_cf2.x).epsilon(cmp_eps));
        CHECK(p0_cf2.y == doctest::Approx(ref_cf2.y).epsilon(cmp_eps));
        CHECK(p0_cf2.z == doctest::Approx(ref_cf2.z).epsilon(cmp_eps));

        // a direction is only rotated -- by the accumulated angle (phi1 + phi2);
        // translation along the chain does not affect points at infinity
        CHECK(e1_cf2.x == doctest::Approx(std::cos(phi1 + phi2)).epsilon(cmp_eps));
        CHECK(e1_cf2.y == doctest::Approx(-std::sin(phi1 + phi2)).epsilon(cmp_eps));
        CHECK(e1_cf2.z == doctest::Approx(0.0).epsilon(cmp_eps));

        // composition consistency: M02 acts as the two consecutive single steps,
        // move2dp(p, M02) == move2dp(move2dp(p, M01), M12)
        auto p0_cf2_via_cf1 = move2dp(move2dp(p0_rf, M01), M12);
        CHECK(p0_cf2_via_cf1.x == doctest::Approx(p0_cf2.x).epsilon(cmp_eps));
        CHECK(p0_cf2_via_cf1.y == doctest::Approx(p0_cf2.y).epsilon(cmp_eps));
        CHECK(p0_cf2_via_cf1.z == doctest::Approx(p0_cf2.z).epsilon(cmp_eps));

        fmt::println("");
        fmt::println("transformation back along the chain (cf2 -> rf)");

        auto M20 = ds.get_pos_trafo(2, 0); // cf2 -> rf (inverse of the chain)
        auto p0_rf_t = move2dp(p0_cf2, M20);
        auto e1_rf_t = move2dp(e1_cf2, M20);

        fmt::println("M20 = {:.3f}", M20);
        fmt::println("p0_rf_t = {:.3f}", p0_rf_t);
        fmt::println("e1_rf_t = {:.3f}", e1_rf_t);

        // round-trip cf2 -> rf recovers the original point and direction
        CHECK(p0_rf_t.x == doctest::Approx(p0_rf.x).epsilon(cmp_eps));
        CHECK(p0_rf_t.y == doctest::Approx(p0_rf.y).epsilon(cmp_eps));
        CHECK(p0_rf_t.z == doctest::Approx(p0_rf.z).epsilon(cmp_eps));
        CHECK(e1_rf_t.x == doctest::Approx(e1_rf.x).epsilon(cmp_eps));
        CHECK(e1_rf_t.y == doctest::Approx(e1_rf.y).epsilon(cmp_eps));
        CHECK(e1_rf_t.z == doctest::Approx(e1_rf.z).epsilon(cmp_eps));

        fmt::println("");
    }

    TEST_CASE("pga2dp: kinematic_system2dp - opening a box lid (static setup + dynamics)")
    {
        fmt::println(
            "pga2dp: kinematic_system2dp - opening a box lid (static setup + dynamics)");

        // Square of side 2 with frame "S" at its lower-left corner (3,2). S is
        // axis-aligned with the world, so the square spans x in [3,5], y in [2,4] and the
        // top-left corner is at (3,4). The lid "L" is hinged at that top-left corner, has
        // length 2, and is opened by an angle vs. the top side (which runs along +x).
        // Point P sits at the end of the lid. Frame tree: W (world) -> S (box) -> L
        // (lid).
        auto side = value_t(2.0);
        auto lid_len = value_t(2.0);
        auto open0 = deg2rad(15); // initial opening angle
        auto constexpr cmp_eps = 1e-10;

        // build W -> S -> L for a given lid opening angle (all frames initially at rest)
        auto build = [&](value_t open_angle) {
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp("W"s)); // world (root)
            ks.add_frame(
                static_frame2dp("S"s, vec2dp{3, 2, 1}, 0.0));      // lower-left corner
            ks.add_frame(static_frame2dp("L"s, vec2dp{0, side, 1}, // hinge (0,2) in S
                                         open_angle));             // parent defaults to S
            return ks;
        };

        auto P_L = vec2dp{lid_len, 0, 1}; // P at the end of the lid, in L coordinates

        // ===== static setup: lid open at 15 deg -> vector to P in W, S, L =====
        auto ks = build(open0);
        auto P_in_L = P_L; // by construction
        auto P_in_S = move2dp(P_L, ks.get_pos_trafo("L", "S"));
        auto P_in_W = move2dp(P_L, ks.get_pos_trafo("L", "W"));

        fmt::println("lid at 15 deg:  P_L={:.3f}  P_S={:.3f}  P_W={:.3f}", P_in_L, P_in_S,
                     P_in_W);

        auto c0 = std::cos(open0), s0 = std::sin(open0);
        CHECK(P_in_L.x == doctest::Approx(lid_len).epsilon(cmp_eps));
        CHECK(P_in_L.y == doctest::Approx(0.0).epsilon(cmp_eps));
        CHECK(P_in_S.x == doctest::Approx(lid_len * c0).epsilon(cmp_eps));
        CHECK(P_in_S.y == doctest::Approx(lid_len * s0 + side).epsilon(cmp_eps));
        CHECK(P_in_W.x == doctest::Approx(lid_len * c0 + 3.0).epsilon(cmp_eps));
        CHECK(P_in_W.y == doctest::Approx(lid_len * s0 + 4.0).epsilon(cmp_eps));

        // ===== add the dynamics: open 15 -> 30 deg in 1 s at constant angular speed
        // =====
        auto open1 = deg2rad(30);
        auto dt = value_t(1.0);
        auto lid_omega = (open1 - open0) / dt; // relative angular speed = 15 deg/s

        ks.set_state("L", kin_state2dp{.omega = lid_omega}); // drive the lid joint

        fmt::println("lid relative angular speed = {:.3f} rad/s = {:.3f} deg/s",
                     lid_omega, rad2deg(lid_omega));

        // relative angular speed of the lid (a scalar; in L, S, or W it is the same
        // value).
        CHECK(ks.relative_twist("L").z == doctest::Approx(lid_omega).epsilon(cmp_eps));
        // S does not rotate, so the lid's WORLD angular speed equals its relative speed
        CHECK(ks.twist_world("L").z == doctest::Approx(lid_omega).epsilon(cmp_eps));
        // instantaneous tip speed = omega * lid_len (tangential to the lid)
        auto v_P = ks.point_velocity(P_in_W, "L");
        CHECK(to_val(bulk_nrm(v_P)) ==
              doctest::Approx(lid_omega * lid_len).epsilon(cmp_eps));

        // ===== new position of P after the lid reached 30 deg =====
        // the angular displacement over dt lands exactly at 30 deg (constant speed)
        CHECK(open0 + lid_omega * dt == doctest::Approx(open1).epsilon(cmp_eps));

        auto ks2 = build(open1);
        auto P_in_L_new = P_L;
        auto P_in_S_new = move2dp(P_L, ks2.get_pos_trafo("L", "S"));
        auto P_in_W_new = move2dp(P_L, ks2.get_pos_trafo("L", "W"));

        fmt::println("lid at 30 deg:  P_L={:.3f}  P_S={:.3f}  P_W={:.3f}", P_in_L_new,
                     P_in_S_new, P_in_W_new);

        auto c1 = std::cos(open1), s1 = std::sin(open1);
        CHECK(P_in_L_new.x ==
              doctest::Approx(lid_len).epsilon(cmp_eps)); // unchanged in L
        CHECK(P_in_L_new.y == doctest::Approx(0.0).epsilon(cmp_eps));
        CHECK(P_in_S_new.x == doctest::Approx(lid_len * c1).epsilon(cmp_eps));
        CHECK(P_in_S_new.y == doctest::Approx(lid_len * s1 + side).epsilon(cmp_eps));
        CHECK(P_in_W_new.x == doctest::Approx(lid_len * c1 + 3.0).epsilon(cmp_eps));
        CHECK(P_in_W_new.y == doctest::Approx(lid_len * s1 + 4.0).epsilon(cmp_eps));

        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // Cross-checks vs. the ga_view body-frame / world-frame scenes
    // (get_frame_trafo_scenes() in ga_view/src/w_mainwindow.cpp, driven by
    // active_frame_trafo). Each scene's body is modelled as a kinematic_system2dp frame
    // "B" at the scene's M0 pose with the scene's body twist Omega_b. We verify our API
    // reproduces the two quantities ga_view tracks:
    //   - the body->world motor  M = M0          (ga_view: M(t) = M0 (x) rexp(0.5 B_b))
    //   - the world generator    B_w = move2dp(Omega_b, M0)   (ga_view's Omega_w)
    /////////////////////////////////////////////////////////////////////////////////////

    // ga_view's M0 builders (mirror make_M0_trans / make_M0_rot in w_mainwindow.cpp)
    inline mvec2dp_u gv_M0_trans(double tx, double ty)
    {
        return rexp(0.5 * vec2dp{-ty, tx, 0.0});
    }
    inline mvec2dp_u gv_M0_rot(double phi) { return rexp(0.5 * vec2dp{0.0, 0.0, phi}); }

    inline void cross_check_scene(vec2dp const& origin, value_t phi,
                                  twist2dp const& Omega_b, mvec2dp_u const& M0_ref,
                                  double cmp_eps = 1e-10)
    {
        kinematic_system2dp ks;
        ks.add_frame(static_frame2dp("W"s));              // world (root)
        ks.add_frame(static_frame2dp("B"s, origin, phi)); // body at M0 pose (B_b = 0)
        ks.set_twist("B", Omega_b);                       // ga_view's body twist Omega_b

        // our body->world motor reproduces ga_view's M0 (compare action on probe points)
        auto const M = ks.get_pos_trafo("B", "W");
        for (auto const& p :
             {vec2dp{0, 0, 1}, vec2dp{1, 0, 1}, vec2dp{0, 1, 1}, vec2dp{2, -3, 1}}) {
            auto const a = move2dp(p, M);
            auto const b = move2dp(p, M0_ref);
            CHECK(a.x == doctest::Approx(b.x).epsilon(cmp_eps));
            CHECK(a.y == doctest::Approx(b.y).epsilon(cmp_eps));
            CHECK(a.z == doctest::Approx(b.z).epsilon(cmp_eps));
        }

        // our world generator B_w reproduces ga_view's Omega_w = move2dp(Omega_b, M0)
        auto const Bw = ks.twist_world("B");
        auto const Bw_ref = move2dp(Omega_b, M0_ref);
        CHECK(Bw.x == doctest::Approx(Bw_ref.x).epsilon(cmp_eps));
        CHECK(Bw.y == doctest::Approx(Bw_ref.y).epsilon(cmp_eps));
        CHECK(Bw.z == doctest::Approx(Bw_ref.z).epsilon(cmp_eps));
    }

    TEST_CASE("pga2dp: ga_view cross-check - Scene 1 (pure translation)")
    {
        // ga_view Scene 1: body slides in +x at vx = 1, so Omega_b = (0,1,0).
        // Upper sub-item M0 = identity; lower M0 = translation to (0,-2).
        fmt::println("pga2dp: ga_view cross-check - Scene 1 (pure translation)");

        auto const Omega_b = twist2dp{0.0, 1.0, 0.0}; // vx = 1, vy = 0

        cross_check_scene(vec2dp{0, 0, 1}, 0.0, Omega_b,
                          gv_M0_rot(0.0)); // upper: M0 = id
        cross_check_scene(vec2dp{0, -2, 1}, 0.0, Omega_b,
                          gv_M0_trans(0.0, -2.0)); // lower
        fmt::println("");
    }

    TEST_CASE("pga2dp: ga_view cross-check - Scene 2 (pure rotation)")
    {
        // ga_view Scene 2: rotation about pivot Q_b = (0.8,0) at omega = pi rad/s, so
        // Omega_b = omega * Q_b = (0.8*pi, 0, pi). Upper M0 = id; lower M0 = T(0,-2).
        fmt::println("pga2dp: ga_view cross-check - Scene 2 (pure rotation)");

        auto const omega = deg2rad(180);                    // pi rad/s
        auto const Omega_b = omega * vec2dp{0.8, 0.0, 1.0}; // = (0.8*pi, 0, pi)

        cross_check_scene(vec2dp{0, 0, 1}, 0.0, Omega_b,
                          gv_M0_rot(0.0)); // upper: M0 = id
        cross_check_scene(vec2dp{0, -2, 1}, 0.0, Omega_b,
                          gv_M0_trans(0.0, -2.0)); // lower
        fmt::println("");
    }

    TEST_CASE("pga2dp: ga_view cross-check - Scene 3 (combined M0 = trans + rot)")
    {
        // ga_view Scene 3: same Omega_b as Scene 2; lower M0 = T(0,-2) (x) R(pi/4), i.e.
        // body home at (0,-2) oriented at 45 deg. Upper M0 = id (reference).
        fmt::println("pga2dp: ga_view cross-check - Scene 3 (combined M0)");

        auto const omega = deg2rad(180);
        auto const Omega_b = omega * vec2dp{0.8, 0.0, 1.0};

        cross_check_scene(vec2dp{0, 0, 1}, 0.0, Omega_b,
                          gv_M0_rot(0.0)); // upper: M0 = id
        cross_check_scene(vec2dp{0, -2, 1}, deg2rad(45), Omega_b,
                          rgpr(gv_M0_trans(0.0, -2.0), gv_M0_rot(deg2rad(45)))); // lower
        fmt::println("");
    }

    TEST_CASE("pga2dp: ga_view cross-check - Scene 4 (world-frame driven)")
    {
        // ga_view Scene 4: upper is a body-driven pure translation (M0 = T(-3,1),
        // Omega_b = (0,1,0)). Lower is WORLD-frame driven: M(t) = T(O(t)) (x) R(omega*t)
        // with a constant world generator Omega_w = (o_z*ox - vy, o_z*oy + vx, o_z). We
        // check the t = 0 configuration: the world-driven motion maps to a body twist
        // Omega_b = move2dp(Omega_w, rrev(M0)) that our system transports back to
        // Omega_w.
        fmt::println("pga2dp: ga_view cross-check - Scene 4 (world-frame driven)");

        // --- upper: body-driven translation ---
        cross_check_scene(vec2dp{-3, 1, 1}, 0.0, twist2dp{0.0, 1.0, 0.0},
                          gv_M0_trans(-3.0, 1.0));

        // --- lower: world-driven, evaluated at t = 0 ---
        auto const cm_ox = -3.0, cm_oy = -1.5; // initial CM position (world)
        auto const cm_vx = 1.0, cm_vy = 0.0;   // CM velocity (world)
        auto const o_z = deg2rad(120);         // spin rate cm_omega = 4*pi/6

        auto const M0_lower = gv_M0_trans(cm_ox, cm_oy); // M(0) = T(O(0)) (x) R(0)
        // ga_view's world generator at t = 0 (tick(), world_frame_drive branch)
        auto const Omega_w0 = vec2dp{o_z * cm_ox - cm_vy, o_z * cm_oy + cm_vx, o_z};
        // the equivalent body twist ga_view derives: Omega_b = move2dp(Omega_w, rrev(M))
        auto const Omega_b0 = move2dp(Omega_w0, rrev(M0_lower));

        cross_check_scene(vec2dp{cm_ox, cm_oy, 1}, 0.0, Omega_b0, M0_lower);

        // and explicitly: our world generator equals ga_view's Omega_w formula
        kinematic_system2dp ks;
        ks.add_frame(static_frame2dp("W"s));
        ks.add_frame(static_frame2dp("B"s, vec2dp{cm_ox, cm_oy, 1}, 0.0));
        ks.set_twist("B", Omega_b0);
        auto const Bw = ks.twist_world("B");
        CHECK(Bw.x == doctest::Approx(Omega_w0.x).epsilon(1e-10));
        CHECK(Bw.y == doctest::Approx(Omega_w0.y).epsilon(1e-10));
        CHECK(Bw.z == doctest::Approx(Omega_w0.z).epsilon(1e-10));
        fmt::println("");
    }

    TEST_CASE("pga2dp: multibody system transformation - branching tree of frames")
    {
        fmt::println(
            "pga2dp: multibody system transformation - branching tree of frames)");

        // Tree topology with explicit parents:

        // Tree topology with explicit parents:
        //
        //              rf[0] (root)
        //              /          \
        //         cfA[1]           cfB[2]   (cfB branches off the root, NOT off cfA)
        //            |
        //         cfC[3]
        //
        static_frame2dp rf; // inertial root frame (identity)

        auto phiA = deg2rad(8);
        auto originA = vec2dp{2, 1.5, 1}; // cfA pose vs. root
        static_frame2dp cfA("A"s, originA, phiA);

        auto phiB = deg2rad(-15);
        auto originB = vec2dp{-1, 3, 1}; // cfB pose vs. root
        static_frame2dp cfB("B"s, originB, phiB);

        auto phiC = deg2rad(20);
        auto originC = vec2dp{0.5, -1, 1}; // cfC pose vs. cfA
        static_frame2dp cfC("C"s, originC, phiC);

        static_system2dp ds;
        ds.add_frame(rf);     // idx 0: root (parent = self)
        ds.add_frame(cfA);    // idx 1: parent defaults to previous -> root (0)
        ds.add_frame(cfB, 0); // idx 2: EXPLICIT parent 0 -> branch off the root
        ds.add_frame(cfC, 1); // idx 3: EXPLICIT parent 1 -> hangs off cfA
        fmt::println("ds = {:.3f}", ds);

        CHECK(!ds.is_linear_chain()); // genuine tree, not a plain chain

        // passive coordinate change one frame deeper: p_child = R(-phi)*(p_parent - o)
        auto to_child = [](vec2dp const& p, vec2dp const& o, value_t phi) {
            auto const c = std::cos(phi), s = std::sin(phi);
            auto const dx = p.x - o.x, dy = p.y - o.y;
            return vec2dp(c * dx + s * dy, -s * dx + c * dy, 1.0);
        };

        auto p_root = vec2dp(4, 5, 1);
        auto constexpr cmp_eps = 1e-10;

        // cfB hangs off the ROOT: a single step, independent of cfA's pose
        auto p_cfB = move2dp(p_root, ds.get_pos_trafo(0, 2));
        auto ref_cfB = to_child(p_root, originB, phiB);
        CHECK(p_cfB.x == doctest::Approx(ref_cfB.x).epsilon(cmp_eps));
        CHECK(p_cfB.y == doctest::Approx(ref_cfB.y).epsilon(cmp_eps));
        CHECK(p_cfB.z == doctest::Approx(ref_cfB.z).epsilon(cmp_eps));

        // cfC hangs off cfA: root -> cfA -> cfC
        auto p_cfC = move2dp(p_root, ds.get_pos_trafo(0, 3));
        auto ref_cfC = to_child(to_child(p_root, originA, phiA), originC, phiC);
        CHECK(p_cfC.x == doctest::Approx(ref_cfC.x).epsilon(cmp_eps));
        CHECK(p_cfC.y == doctest::Approx(ref_cfC.y).epsilon(cmp_eps));
        CHECK(p_cfC.z == doctest::Approx(ref_cfC.z).epsilon(cmp_eps));

        // cross-branch transform cfB -> cfC: walk up to the LCA (root), then down to cfC
        auto M_BC = ds.get_pos_trafo(2, 3);
        auto p_cfC_from_B = move2dp(p_cfB, M_BC);
        CHECK(p_cfC_from_B.x == doctest::Approx(p_cfC.x).epsilon(cmp_eps));
        CHECK(p_cfC_from_B.y == doctest::Approx(p_cfC.y).epsilon(cmp_eps));
        CHECK(p_cfC_from_B.z == doctest::Approx(p_cfC.z).epsilon(cmp_eps));

        // a direction is path-independent: cfB -> cfC composed equals the direct root ->
        // cfC rotation (translation along either branch does not affect points at inf.)
        auto e1_cfB = move2dp(e1_2dp, ds.get_pos_trafo(0, 2));
        auto e1_cfC_from_B = move2dp(e1_cfB, M_BC);
        auto e1_cfC_direct = move2dp(e1_2dp, ds.get_pos_trafo(0, 3));
        CHECK(e1_cfC_from_B.x == doctest::Approx(e1_cfC_direct.x).epsilon(cmp_eps));
        CHECK(e1_cfC_from_B.y == doctest::Approx(e1_cfC_direct.y).epsilon(cmp_eps));
        CHECK(e1_cfC_from_B.z == doctest::Approx(e1_cfC_direct.z).epsilon(cmp_eps));

        // round-trip across branches: cfB -> cfC -> cfB recovers the original point
        auto p_cfB_rt = move2dp(p_cfC_from_B, ds.get_pos_trafo(3, 2));
        CHECK(p_cfB_rt.x == doctest::Approx(p_cfB.x).epsilon(cmp_eps));
        CHECK(p_cfB_rt.y == doctest::Approx(p_cfB.y).epsilon(cmp_eps));
        CHECK(p_cfB_rt.z == doctest::Approx(p_cfB.z).epsilon(cmp_eps));

        fmt::println("");
    }

    TEST_CASE("pga2dp: multibody system transformation - merry-go-round (platform + 3 "
              "turntables)")
    {
        fmt::println("pga2dp: multibody system transformation - merry-go-round)");

        // Tree topology: a rotating main platform carrying three rotating turntables.
        //
        //                 ground[0] (root)
        //                     |
        //                 platform[1]              (rotates about its centre by theta)
        //               /     |      \
        //          tt0[2]   tt1[3]   tt2[4]        (mounted at r, 120 deg apart, own
        //          spin)
        //
        static_frame2dp ground; // inertial root

        auto theta = deg2rad(30);      // platform rotation vs. ground
        auto centre = vec2dp{1, 2, 1}; // platform centre in ground coordinates
        static_frame2dp platform("P"s, centre, theta);

        auto r = value_t(2.0);        // mounting radius on the platform
        auto mount = [r](value_t a) { // turntable origin in PLATFORM coordinates
            return vec2dp(r * std::cos(a), r * std::sin(a), 1.0);
        };
        auto psi0 = deg2rad(10); // individual turntable spins vs. the platform
        auto psi1 = deg2rad(25);
        auto psi2 = deg2rad(-5);
        static_frame2dp tt0("T0"s, mount(deg2rad(0)), psi0);
        static_frame2dp tt1("T1"s, mount(deg2rad(120)), psi1);
        static_frame2dp tt2("T2"s, mount(deg2rad(240)), psi2);

        static_system2dp ds;
        ds.add_frame(ground);   // idx 0: root
        ds.add_frame(platform); // idx 1: parent defaults to ground
        ds.add_frame(tt0, 1);   // idx 2: mounted on the platform
        ds.add_frame(tt1, 1);   // idx 3: mounted on the platform
        ds.add_frame(tt2, 1);   // idx 4: mounted on the platform
        fmt::println("ds = {:.3f}", ds);

        CHECK(!ds.is_linear_chain()); // three children off one parent -> a tree

        auto constexpr cmp_eps = 1e-10;

        // passive (parent -> child) and active (child -> parent) coordinate changes
        auto to_child = [](vec2dp const& p, vec2dp const& o, value_t phi) {
            auto const c = std::cos(phi), s = std::sin(phi);
            auto const dx = p.x - o.x, dy = p.y - o.y;
            return vec2dp(c * dx + s * dy, -s * dx + c * dy, 1.0);
        };
        auto to_parent = [](vec2dp const& q, vec2dp const& o, value_t phi) {
            auto const c = std::cos(phi), s = std::sin(phi);
            return vec2dp(c * q.x - s * q.y + o.x, s * q.x + c * q.y + o.y, 1.0);
        };

        // each turntable origin (its (0,0,1) point) expressed in ground coordinates,
        // checked against R(theta) * mount_point + centre
        auto tt_origin_ground = [&](size_t idx, value_t mount_angle) {
            auto o_ground = move2dp(O_2dp, ds.get_pos_trafo(idx, 0));
            auto ref = to_parent(mount(mount_angle), centre, theta);
            CHECK(o_ground.x == doctest::Approx(ref.x).epsilon(cmp_eps));
            CHECK(o_ground.y == doctest::Approx(ref.y).epsilon(cmp_eps));
            CHECK(o_ground.z == doctest::Approx(ref.z).epsilon(cmp_eps));
            return o_ground;
        };
        auto g0 = tt_origin_ground(2, deg2rad(0));
        auto g1 = tt_origin_ground(3, deg2rad(120));
        auto g2 = tt_origin_ground(4, deg2rad(240));

        fmt::println("tt0 origin in ground = {:.3f}", g0);
        fmt::println("tt1 origin in ground = {:.3f}", g1);
        fmt::println("tt2 origin in ground = {:.3f}", g2);

        // merry-go-round invariants: every turntable sits at radius r from the centre
        // (the platform rotation preserves it), and 120 deg spacing makes adjacent
        // turntables equidistant with chord length r * rsqrt(3)
        CHECK(to_val(bulk_nrm(g0 - centre)) == doctest::Approx(r).epsilon(cmp_eps));
        CHECK(to_val(bulk_nrm(g1 - centre)) == doctest::Approx(r).epsilon(cmp_eps));
        CHECK(to_val(bulk_nrm(g2 - centre)) == doctest::Approx(r).epsilon(cmp_eps));

        auto chord = r * std::sqrt(3.0);
        CHECK(to_val(bulk_nrm(g1 - g0)) == doctest::Approx(chord).epsilon(cmp_eps));
        CHECK(to_val(bulk_nrm(g2 - g1)) == doctest::Approx(chord).epsilon(cmp_eps));
        CHECK(to_val(bulk_nrm(g0 - g2)) == doctest::Approx(chord).epsilon(cmp_eps));

        // turntable-to-turntable transform routes up to the platform (LCA) and back
        // down: a point given in tt0 coordinates, expressed in tt2 coordinates
        auto p_tt0 = vec2dp(0.3, -0.4, 1);
        auto p_tt2 = move2dp(p_tt0, ds.get_pos_trafo(2, 4));
        auto ref_tt2 = to_child(to_parent(p_tt0, mount(deg2rad(0)), psi0),
                                mount(deg2rad(240)), psi2);
        CHECK(p_tt2.x == doctest::Approx(ref_tt2.x).epsilon(cmp_eps));
        CHECK(p_tt2.y == doctest::Approx(ref_tt2.y).epsilon(cmp_eps));
        CHECK(p_tt2.z == doctest::Approx(ref_tt2.z).epsilon(cmp_eps));

        // round-trip tt0 -> tt2 -> tt0
        auto p_tt0_rt = move2dp(p_tt2, ds.get_pos_trafo(4, 2));
        CHECK(p_tt0_rt.x == doctest::Approx(p_tt0.x).epsilon(cmp_eps));
        CHECK(p_tt0_rt.y == doctest::Approx(p_tt0.y).epsilon(cmp_eps));
        CHECK(p_tt0_rt.z == doctest::Approx(p_tt0.z).epsilon(cmp_eps));

        fmt::println("");
    }

    TEST_CASE("pga2dp: static_system2dp - name-addressed frames")
    {
        fmt::println("pga2dp: static_system2dp - name-addressed frames");

        static_frame2dp ground;                                       // name "I"
        static_frame2dp platform("P"s, vec2dp{1, 2, 1}, deg2rad(30)); // name "P"
        static_frame2dp turntable("T0"s, vec2dp{2, 0, 1}, deg2rad(10));

        static_system2dp ds;
        ds.add_frame(ground);
        ds.add_frame(platform);
        ds.add_frame(turntable, 1);

        // index_of resolves names to the add order
        CHECK(ds.index_of("I") == 0);
        CHECK(ds.index_of("P") == 1);
        CHECK(ds.index_of("T0") == 2);

        // name-addressed get_pos_trafo matches the index-addressed result
        auto p = vec2dp(0.4, -0.7, 1);
        auto by_name = move2dp(p, ds.get_pos_trafo("T0", "I"));
        auto by_idx = move2dp(p, ds.get_pos_trafo(2, 0));
        auto constexpr cmp_eps = 1e-10;
        CHECK(by_name.x == doctest::Approx(by_idx.x).epsilon(cmp_eps));
        CHECK(by_name.y == doctest::Approx(by_idx.y).epsilon(cmp_eps));
        CHECK(by_name.z == doctest::Approx(by_idx.z).epsilon(cmp_eps));

        // unknown name throws
        CHECK_THROWS_AS(ds.index_of("nope"), std::runtime_error);

        // duplicate frame name throws on add
        static_system2dp ds2;
        ds2.add_frame(static_frame2dp("A"s));
        CHECK_THROWS_AS(ds2.add_frame(static_frame2dp("A"s)), std::runtime_error);

        fmt::println("");
    }

    TEST_CASE("pga2dp: kinematic_system2dp - velocity of a spinning wheel (intuitive)")
    {
        // A single wheel spinning counter-clockwise at omega = 2 rad/s about the world
        // origin. Physics 101: a point at radius r moves TANGENTIALLY with speed omega*r.
        fmt::println("pga2dp: kinematic_system2dp - velocity of a spinning wheel");

        auto omega = value_t(2.0); // [rad/s], CCW positive

        kinematic_system2dp ks;
        ks.add_frame(static_frame2dp{});            // ground (root)
        ks.add_frame(static_frame2dp("wheel"s),     // wheel centred at the origin ...
                     kin_state2dp{.omega = omega}); // ... spinning at 2 rad/s

        auto constexpr cmp_eps = 1e-12;

        // a point 3 to the RIGHT of the centre moves straight UP at speed omega*3 = 6
        auto v_right = ks.point_velocity(vec2dp{3, 0, 1}, "wheel");
        fmt::println("  point (3,0): v = {:.3f}   expect (0, 6)", v_right);
        CHECK(v_right.x == doctest::Approx(0.0).epsilon(cmp_eps));
        CHECK(v_right.y == doctest::Approx(6.0).epsilon(cmp_eps));

        // a point 3 ABOVE the centre moves straight LEFT at speed 6
        auto v_top = ks.point_velocity(vec2dp{0, 3, 1}, "wheel");
        fmt::println("  point (0,3): v = {:.3f}   expect (-6, 0)", v_top);
        CHECK(v_top.x == doctest::Approx(-6.0).epsilon(cmp_eps));
        CHECK(v_top.y == doctest::Approx(0.0).epsilon(cmp_eps));

        // the centre itself does not move
        auto v_centre = ks.point_velocity(vec2dp{0, 0, 1}, "wheel");
        fmt::println("  point (0,0): v = {:.3f}   expect (0, 0)", v_centre);
        CHECK(v_centre.x == doctest::Approx(0.0).epsilon(cmp_eps));
        CHECK(v_centre.y == doctest::Approx(0.0).epsilon(cmp_eps));

        // speed grows linearly with radius: at r = 5, speed = omega*5 = 10
        auto v_far = ks.point_velocity(vec2dp{5, 0, 1}, "wheel");
        CHECK(to_val(bulk_nrm(v_far)) == doctest::Approx(10.0).epsilon(cmp_eps));

        fmt::println("");
    }

    TEST_CASE(
        "pga2dp: kinematic_system2dp - acceleration of a spinning wheel (intuitive)")
    {
        // Same wheel, now the acceleration. Uniform spin (omega const) => pure
        // CENTRIPETAL acceleration a = -omega^2 * r, pointing toward the centre. Spinning
        // up (alpha) adds a TANGENTIAL term alpha * r.
        fmt::println("pga2dp: kinematic_system2dp - acceleration of a spinning wheel");

        auto omega = value_t(2.0); // [rad/s]
        auto constexpr cmp_eps = 1e-12;

        // --- constant spin: pure centripetal, magnitude omega^2 * r = 4 * r ---
        {
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp{});
            ks.add_frame(static_frame2dp("wheel"s), kin_state2dp{.omega = omega});

            // point (3,0): centripetal points to centre (-x), magnitude 4*3 = 12
            auto a_right = ks.point_acceleration(vec2dp{3, 0, 1}, "wheel");
            fmt::println("  point (3,0): a = {:.3f}   expect (-12, 0)", a_right);
            CHECK(a_right.x == doctest::Approx(-12.0).epsilon(cmp_eps));
            CHECK(a_right.y == doctest::Approx(0.0).epsilon(cmp_eps));

            // point (0,3): centripetal points to centre (-y): (0, -12)
            auto a_top = ks.point_acceleration(vec2dp{0, 3, 1}, "wheel");
            fmt::println("  point (0,3): a = {:.3f}   expect (0, -12)", a_top);
            CHECK(a_top.x == doctest::Approx(0.0).epsilon(cmp_eps));
            CHECK(a_top.y == doctest::Approx(-12.0).epsilon(cmp_eps));

            // the centre has no acceleration
            auto a_centre = ks.point_acceleration(vec2dp{0, 0, 1}, "wheel");
            CHECK(a_centre.x == doctest::Approx(0.0).epsilon(cmp_eps));
            CHECK(a_centre.y == doctest::Approx(0.0).epsilon(cmp_eps));
        }

        // --- spinning up: alpha = 5 rad/s^2 adds a tangential term alpha * r ---
        {
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp{});
            ks.add_frame(static_frame2dp("wheel"s),
                         kin_state2dp{.omega = omega, .alpha = 5.0});

            // point (3,0): centripetal (-12,0) + tangential (alpha*r in +y = 15) =
            // (-12,15)
            auto a = ks.point_acceleration(vec2dp{3, 0, 1}, "wheel");
            fmt::println("  point (3,0), alpha=5: a = {:.3f}   expect (-12, 15)", a);
            CHECK(a.x == doctest::Approx(-12.0).epsilon(cmp_eps));
            CHECK(a.y == doctest::Approx(15.0).epsilon(cmp_eps));
        }

        fmt::println("");
    }

    TEST_CASE("pga2dp: kinematic_system2dp - acceleration on a merry-go-round (Coriolis)")
    {
        // A turntable spinning on a rotating platform: the point's acceleration now mixes
        // two centripetal terms AND a Coriolis coupling between the two spins. Validated
        // against the second central difference of the true motor evolution -- no
        // analytic shortcut, so this exercises the full Newton-Euler propagation (the
        // bracket term).
        fmt::println("pga2dp: kinematic_system2dp - acceleration on a merry-go-round");

        auto Omega = deg2rad(40);        // platform angular velocity  [rad/s]
        auto omega = deg2rad(90);        // turntable spin vs. platform [rad/s]
        auto centre = vec2dp{1, 2, 1};   // platform centre in ground coords
        auto mount = vec2dp{2, 0, 1};    // turntable origin in platform coords
        auto X_tt = vec2dp{0.3, 0.1, 1}; // a point fixed on the turntable

        // world position of X_tt at time t (poses advanced by the constant spin rates)
        auto world_pos = [&](double t) {
            static_system2dp ds;
            ds.add_frame(static_frame2dp{});
            ds.add_frame(static_frame2dp("P"s, centre, Omega * t));
            ds.add_frame(static_frame2dp("T"s, mount, omega * t), 1);
            return move2dp(X_tt, ds.get_pos_trafo(2, 0));
        };

        // momentary kinematic system at t = 0
        kinematic_system2dp ks;
        ks.add_frame(static_frame2dp{});
        ks.add_frame(static_frame2dp("P"s, centre, 0.0), kin_state2dp{.omega = Omega});
        ks.add_frame(static_frame2dp("T"s, mount, 0.0), kin_state2dp{.omega = omega}, 1);

        auto X_world0 = move2dp(X_tt, ks.get_pos_trafo(2, 0));
        auto a_formula = ks.point_acceleration(X_world0, "T");

        // second central difference: a ~ (p(+dt) - 2 p(0) + p(-dt)) / dt^2
        double const dt = 1e-4;
        auto p0 = world_pos(0.0);
        auto pp = world_pos(dt);
        auto pm = world_pos(-dt);
        auto a_fd = vec2dp{(pp.x - 2 * p0.x + pm.x) / (dt * dt),
                           (pp.y - 2 * p0.y + pm.y) / (dt * dt), 0.0};

        fmt::println("  a_formula = {:.3f}", a_formula);
        fmt::println("  a_fd      = {:.3f}", a_fd);

        // match the finite-difference reference (loose eps: 2nd difference is O(dt^2))
        CHECK(a_formula.x == doctest::Approx(a_fd.x).epsilon(1e-4));
        CHECK(a_formula.y == doctest::Approx(a_fd.y).epsilon(1e-4));

        fmt::println("");
    }

    TEST_CASE("pga2dp: kinematic_system2dp - turntable-to-turntable kinematics")
    {
        // Merry-go-round: platform spins at Omega; turntables tt0, tt1, tt2 (120 deg
        // apart, all children of the platform) each spin at their own rate. For a point P
        // on tt0:
        //   (a) its world velocity and its position expressed in tt2 (cross-branch, LCA =
        //   P); (b) the velocity an observer RIDING tt2 measures for P (transport
        //   theorem),
        //       validated against the finite difference of P's position in tt2 over time.
        fmt::println("pga2dp: kinematic_system2dp - turntable-to-turntable kinematics");

        auto Omega = deg2rad(40); // platform angular velocity
        auto w0 = deg2rad(90);    // tt0 spin vs. platform
        auto w1 = deg2rad(60);    // tt1 spin vs. platform
        auto w2 = deg2rad(-30);   // tt2 spin vs. platform
        auto pP = deg2rad(25); // initial angles (so observer-frame rotation is exercised)
        auto p0 = deg2rad(10), p1 = deg2rad(5), p2 = deg2rad(40);
        auto centre = vec2dp{1, 2, 1};
        auto r = value_t(2.0);
        auto mount = [r](value_t a) {
            return vec2dp(r * std::cos(a), r * std::sin(a), 1.0);
        };
        auto P_tt0 = vec2dp{0.3, 0.1, 1}; // point fixed on tt0

        auto constexpr cmp_eps = 1e-10;

        // momentary kinematic system at t = 0
        kinematic_system2dp ks;
        ks.add_frame(static_frame2dp("W"s));
        ks.add_frame(static_frame2dp("P"s, centre, pP), kin_state2dp{.omega = Omega});
        ks.add_frame(static_frame2dp("T0"s, mount(deg2rad(0)), p0),
                     kin_state2dp{.omega = w0}, ks.index_of("P"));
        ks.add_frame(static_frame2dp("T1"s, mount(deg2rad(120)), p1),
                     kin_state2dp{.omega = w1}, ks.index_of("P"));
        ks.add_frame(static_frame2dp("T2"s, mount(deg2rad(240)), p2),
                     kin_state2dp{.omega = w2}, ks.index_of("P"));

        // (a) world velocity of P and its position expressed in tt2 (cross-branch)
        // --------
        auto P_world = move2dp(P_tt0, ks.get_pos_trafo("T0", "W"));
        auto v_world = ks.point_velocity(P_world, "T0");
        auto P_in_tt2 = move2dp(P_tt0, ks.get_pos_trafo("T0", "T2"));
        fmt::println("  P world pos = {:.3f}, world vel = {:.3f}", P_world, v_world);
        fmt::println("  P in tt2    = {:.3f}", P_in_tt2);

        // relative velocity w.r.t. the inertial world equals the absolute world velocity
        auto v_rel_world = ks.relative_point_velocity(P_tt0, "T0", "W");
        CHECK(v_rel_world.x == doctest::Approx(v_world.x).epsilon(cmp_eps));
        CHECK(v_rel_world.y == doctest::Approx(v_world.y).epsilon(cmp_eps));
        // an observer reads zero for a point fixed in their own frame
        auto v_self = ks.relative_point_velocity(P_tt0, "T0", "T0");
        CHECK(to_val(bulk_nrm(v_self)) == doctest::Approx(0.0).epsilon(cmp_eps));

        // (b) velocity of P as measured by an observer riding tt2, vs finite difference
        // ----
        auto v_rel_tt2 = ks.relative_point_velocity(P_tt0, "T0", "T2");

        // P's position in tt2 coordinates, with all poses advanced by time t
        auto P_in_tt2_at = [&](double t) {
            static_system2dp ds;
            ds.add_frame(static_frame2dp("W"s));
            ds.add_frame(static_frame2dp("P"s, centre, pP + Omega * t));
            ds.add_frame(static_frame2dp("T0"s, mount(deg2rad(0)), p0 + w0 * t),
                         ds.index_of("P"));
            ds.add_frame(static_frame2dp("T1"s, mount(deg2rad(120)), p1 + w1 * t),
                         ds.index_of("P"));
            ds.add_frame(static_frame2dp("T2"s, mount(deg2rad(240)), p2 + w2 * t),
                         ds.index_of("P"));
            return move2dp(P_tt0, ds.get_pos_trafo("T0", "T2"));
        };
        double const dt = 1e-5;
        auto pp = P_in_tt2_at(dt);
        auto pm = P_in_tt2_at(-dt);
        auto v_fd = vec2dp{(pp.x - pm.x) / (2 * dt), (pp.y - pm.y) / (2 * dt), 0.0};

        fmt::println("  v_rel (tt2 observer) = {:.3f}, finite-diff = {:.3f}", v_rel_tt2,
                     v_fd);
        CHECK(v_rel_tt2.x == doctest::Approx(v_fd.x).epsilon(1e-6));
        CHECK(v_rel_tt2.y == doctest::Approx(v_fd.y).epsilon(1e-6));

        // (c) acceleration of P as measured by the tt2 observer, vs 2nd central
        // difference --
        auto a_rel_tt2 = ks.relative_point_acceleration(P_tt0, "T0", "T2");
        double const dt_a =
            1e-4; // a 2nd difference needs a larger step than a 1st difference
        auto qp = P_in_tt2_at(dt_a);
        auto q0 = P_in_tt2_at(0.0);
        auto qm = P_in_tt2_at(-dt_a);
        auto a_fd = vec2dp{(qp.x - 2 * q0.x + qm.x) / (dt_a * dt_a),
                           (qp.y - 2 * q0.y + qm.y) / (dt_a * dt_a), 0.0};
        fmt::println("  a_rel (tt2 observer) = {:.3f}, finite-diff = {:.3f}", a_rel_tt2,
                     a_fd);
        CHECK(a_rel_tt2.x == doctest::Approx(a_fd.x).epsilon(1e-4));
        CHECK(a_rel_tt2.y == doctest::Approx(a_fd.y).epsilon(1e-4));

        // an observer reads zero acceleration for a point fixed in their own frame
        auto a_self = ks.relative_point_acceleration(P_tt0, "T0", "T0");
        CHECK(to_val(bulk_nrm(a_self)) == doctest::Approx(0.0).epsilon(cmp_eps));

        fmt::println("");
    }

    TEST_CASE("pga2dp: kinematic_system2dp - step(dt) time evolution")
    {
        // step(dt) advances each frame's relative pose on the motor manifold by
        // rexp(0.5 B dt) and ramps its velocity by the relative acceleration. For
        // CONSTANT twist the pose evolution is EXACT: N steps of dt match a single
        // rebuild at t = N*dt to ~machine precision.
        fmt::println("pga2dp: kinematic_system2dp - step(dt) time evolution");
        auto constexpr cmp_eps = 1e-9;

        // (a) constant spin: stepping a wheel == rebuilding at the accumulated angle
        // -------
        {
            auto omega = deg2rad(50);
            auto P_B = vec2dp{0.5, -0.3, 1};
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp("W"s));
            ks.add_frame(static_frame2dp("B"s, vec2dp{2, 1, 1}, deg2rad(10)),
                         kin_state2dp{.omega = omega});

            double const T = 1.0;
            int const N = 100;
            for (int k = 0; k < N; ++k)
                ks.step(T / N);

            // reference: rebuild at the final angle 10 deg + omega*T (origin stays fixed)
            static_system2dp ref;
            ref.add_frame(static_frame2dp("W"s));
            ref.add_frame(
                static_frame2dp("B"s, vec2dp{2, 1, 1}, deg2rad(10) + omega * T));

            auto P_step = move2dp(P_B, ks.get_pos_trafo("B", "W"));
            auto P_ref = move2dp(P_B, ref.get_pos_trafo("B", "W"));
            fmt::println("  stepped P = {:.3f}, rebuilt P = {:.3f}", P_step, P_ref);
            CHECK(P_step.x == doctest::Approx(P_ref.x).epsilon(cmp_eps));
            CHECK(P_step.y == doctest::Approx(P_ref.y).epsilon(cmp_eps));
            CHECK(ks.relative_twist("B").z == doctest::Approx(omega).epsilon(cmp_eps));
        }

        // (b) opening a lid: step 15 -> 30 deg over 1 s and check the tip
        // ------------------
        {
            auto side = value_t(2.0), lid_len = value_t(2.0);
            auto lid_omega = deg2rad(15); // 15 deg/s -> 30 deg after 1 s
            auto P_L = vec2dp{lid_len, 0, 1};

            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp("W"s));
            ks.add_frame(static_frame2dp("S"s, vec2dp{3, 2, 1}, 0.0));
            ks.add_frame(static_frame2dp("L"s, vec2dp{0, side, 1}, deg2rad(15)),
                         kin_state2dp{.omega = lid_omega});

            double const T = 1.0;
            int const N = 200;
            for (int k = 0; k < N; ++k)
                ks.step(T / N);

            auto P_W = move2dp(P_L, ks.get_pos_trafo("L", "W"));
            auto c = std::cos(deg2rad(30)), s = std::sin(deg2rad(30));
            fmt::println("  lid tip after opening: P_W = {:.3f}", P_W);
            CHECK(P_W.x == doctest::Approx(lid_len * c + 3.0).epsilon(cmp_eps));
            CHECK(P_W.y == doctest::Approx(lid_len * s + 4.0).epsilon(cmp_eps));
        }

        // (c) constant angular acceleration ramps the velocity exactly
        // --------------------
        {
            auto alpha = deg2rad(20);
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp("W"s));
            ks.add_frame(static_frame2dp("B"s),
                         kin_state2dp{.omega = 0.0, .alpha = alpha});

            double const T = 1.0;
            int const N = 500;
            for (int k = 0; k < N; ++k)
                ks.step(T / N);

            // omega(T) = alpha*T is an exact sum of equal Euler increments
            CHECK(ks.relative_twist("B").z ==
                  doctest::Approx(alpha * T).epsilon(cmp_eps));
            // the angle ~ 0.5*alpha*T^2 (explicit Euler, first order in dt) -> loose
            // check
            CHECK(rad2deg(ks.frame(1).get_pose().phi) ==
                  doctest::Approx(rad2deg(0.5 * alpha * T * T)).epsilon(2e-2));
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: kinematic_system2dp - inertial vs non-inertial frames (intuitive)")
    {
        // An inertial frame neither rotates nor accelerates. Constant linear velocity is
        // fine; any rotation or any acceleration makes a frame non-inertial.
        fmt::println("pga2dp: kinematic_system2dp - inertial vs non-inertial frames");

        // the ground / root is the inertial reference
        {
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp{}); // ground
            CHECK(ks.is_inertial_frame("I"));
            CHECK(ks.is_inertial_system());
        }

        // a frame gliding at CONSTANT velocity is still inertial (no spin, no
        // acceleration)
        {
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp{});
            ks.add_frame(static_frame2dp("glider"s, vec2dp{1, 0, 1}),
                         kin_state2dp{.vel = vec2dp{2, 0, 0}}); // constant 2 in +x
            CHECK(ks.is_inertial_frame("glider"));
            CHECK(ks.is_inertial_system()); // ground AND glider are inertial
        }

        // a frame that ACCELERATES linearly is non-inertial
        {
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp{});
            ks.add_frame(static_frame2dp("rocket"s, vec2dp{1, 0, 1}),
                         kin_state2dp{.acc = vec2dp{0, 3, 0}});
            CHECK_FALSE(ks.is_inertial_frame("rocket"));
            CHECK_FALSE(ks.is_inertial_system());
        }

        // a spinning frame is non-inertial (rotating frames always are)
        {
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp{});
            ks.add_frame(static_frame2dp("wheel"s), kin_state2dp{.omega = 2.0});
            CHECK_FALSE(ks.is_inertial_frame("wheel"));
        }

        // a turntable that does NOT spin relative to a rotating platform is STILL
        // non-inertial: it inherits the platform's rotation (and its origin orbits)
        {
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp{});
            ks.add_frame(static_frame2dp("P"s, vec2dp{0, 0, 1}),
                         kin_state2dp{.omega = 1.0});
            ks.add_frame(static_frame2dp("T"s, vec2dp{2, 0, 1}),
                         1); // mounted, no own spin
            CHECK_FALSE(ks.is_inertial_frame("P"));
            CHECK_FALSE(ks.is_inertial_frame("T")); // inherits the platform's rotation
        }

        fmt::println("");
    }

    TEST_CASE("pga2dp: kinematic_system2dp - world twist & point velocity")
    {
        fmt::println("pga2dp: kinematic_system2dp - world twist & point velocity");

        auto constexpr cmp_eps = 1e-10;

        // --- pure translation: a body frame sliding at v=(1,0) relative to ground ---
        {
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp{});                     // ground (root, at rest)
            ks.add_frame(static_frame2dp("B"s, vec2dp{2, 0, 1}), // body at (2,0)
                         kin_state2dp{.vel = vec2dp{1.0, 0.0, 0.0}}); // sliding v = (1,0)

            // a translating frame has zero angular rate ...
            CHECK(ks.twist_world(1).z == doctest::Approx(0.0).epsilon(cmp_eps));

            // ... and every attached point moves with the same velocity (1,0)
            auto vA = ks.point_velocity(vec2dp{5, 7, 1}, 1);
            auto vB = ks.point_velocity(vec2dp{-3, 2, 1}, "B");
            CHECK(vA.x == doctest::Approx(1.0).epsilon(cmp_eps));
            CHECK(vA.y == doctest::Approx(0.0).epsilon(cmp_eps));
            CHECK(vB.x == doctest::Approx(1.0).epsilon(cmp_eps));
            CHECK(vB.y == doctest::Approx(0.0).epsilon(cmp_eps));
        }

        // --- merry-go-round: a point on a turntable mounted on a rotating platform,
        //     validated by finite-differencing the actual motor evolution over time ---
        {
            auto Omega = deg2rad(40);        // platform angular velocity  [rad/s]
            auto omega = deg2rad(90);        // turntable spin vs. platform [rad/s]
            auto centre = vec2dp{1, 2, 1};   // platform centre in ground coords
            auto mount = vec2dp{2, 0, 1};    // turntable origin in platform coords
            auto X_tt = vec2dp{0.3, 0.1, 1}; // a point fixed on the turntable

            // world position of X_tt at time t: poses advanced by the constant rates
            // (platform spins about its centre, turntable spins about its own origin)
            auto world_pos = [&](double t) {
                static_system2dp ds;
                ds.add_frame(static_frame2dp{});
                ds.add_frame(static_frame2dp("P"s, centre, Omega * t));
                ds.add_frame(static_frame2dp("T"s, mount, omega * t), 1);
                return move2dp(X_tt, ds.get_pos_trafo(2, 0));
            };

            // momentary kinematic system at t = 0
            kinematic_system2dp ks;
            ks.add_frame(static_frame2dp{});
            ks.add_frame(static_frame2dp("P"s, centre, 0.0),
                         kin_state2dp{.omega = Omega}); // platform spins about its centre
            ks.add_frame(static_frame2dp("T"s, mount, 0.0), kin_state2dp{.omega = omega},
                         1); // turntable spins about origin

            // formatters: print a single momentary state and the whole configured system
            // (each line shows the frame pose plus its decoded relative kinematic state)
            fmt::println("a single state: {:.3f}",
                         kin_state2dp{.vel = vec2dp{1.0, -0.5, 0.0}, .omega = 2.0});
            fmt::println("ks = {:.4f}", ks);

            auto X_world0 = move2dp(X_tt, ks.get_pos_trafo(2, 0));
            auto v_formula = ks.point_velocity(X_world0, "T");

            // central finite difference of the true motor evolution
            double const dt = 1e-5;
            auto p_plus = world_pos(dt);
            auto p_minus = world_pos(-dt);
            auto v_fd = vec2dp{(p_plus.x - p_minus.x) / (2 * dt),
                               (p_plus.y - p_minus.y) / (2 * dt), 0.0};

            fmt::println("X_world0  = {:.3f}", X_world0);
            fmt::println("v_formula = {:.3f}", v_formula);
            fmt::println("v_fd      = {:.3f}", v_fd);

            // world angular rate is the sum of the two spins (omega + Omega)
            CHECK(ks.twist_world("T").z ==
                  doctest::Approx(Omega + omega).epsilon(cmp_eps));

            // point_velocity matches the finite-difference reference (loose eps: O(dt^2))
            CHECK(v_formula.x == doctest::Approx(v_fd.x).epsilon(1e-6));
            CHECK(v_formula.y == doctest::Approx(v_fd.y).epsilon(1e-6));
        }

        fmt::println("");
    }

    // TODO: show more complex setups
    //
    //       - combined rotation and translation
    //       - movement of rigid bodies with continuously distributed mass
    //       - two-paticles systems with internal coupling forces
    //       - coupled systems in general with more complex inertia

    //////////////////////////////////////
    // save a version for integrating multiple points (mainly to show u memory acces)
    //////////////////////////////////////

    // TEST_CASE("pga2dp: multi_body")
    // {
    //     fmt::println("pga2dp: multi_body");


    //     class sim_ode_multi_body_pga2dp { // model 2nd order ode by a 1st order system

    //       public:

    //         sim_ode_multi_body_pga2dp(size_t npts_in, std::vector<vec2dp> pos_in,
    //         std::vector<value_t> m_in,
    //                 vec2dp const& cm_pos_in, vec2dp const& cm_spd_in, value_t
    //                 cm_phi_in, value_t cm_omega_in) :
    //             npts(npts_in), pts(std::move(pos_in)), m(std::move(m_in)),
    //             cm_w_pos0(cm_pos_in), cm_w_spd0(cm_spd_in),
    //             cm_w_phi0(cm_phi_in), cm_w_omega0(cm_omega_in), u_mem(npts *
    //             2), uh_mem(2 * npts * 2), rhs_mem(npts * 2)
    //         {
    //             if (npts < 2) {
    //                 throw std::invalid_argument("sim_ode: rigid body simulation
    //                 requires"
    //                                             " at least two discrete points.");
    //             }
    //             if ((pts.size() != npts) || (m.size() != npts)) {
    //                 throw std::invalid_argument(
    //                     "sim_ode: position and mass vectors require same size.");
    //             }

    //             // build inertia map from discrete points
    //             // maps rate of change (a vector in 2D) to the momentum bivector
    //             for (size_t n = 0; n < npts; ++n) {
    //                 I += get_point_inertia(m[n], pts[n]);
    //             }

    //             // get inverse inertia map
    //             // maps a momentum bivector to the rate of change (a vector in 2D)
    //             I_inv = get_inertia_inverse(I);
    //         }

    //         void set_initial_values()
    //         {

    //             // Create mdspan view for setting initial values
    //             auto u = mdspan<vec2dp, dextents<size_t, 2>>(u_mem.data(), npts, 2);

    //             for (size_t n = 0; n < npts; ++n) {

    //                 // initial transformation of "position" encoded in B-vector
    //                 // encoding:
    //                 // B_rot=(x0_fix, y0_phi, 1) * phi0
    //                 // B_tra=(- y0_trans, y0_phi*phi0 + x0_trans, phi0)
    //                 // typical starting value is "no initial transformation":
    //                 // B0=(0,0,0) => M0 = pscalar2dp(1)
    //                 //               (=identiy transformation at t=0)
    //                 u[n, 0] = vec2dp(0.0, 0.0, 0.0);

    //                 // initial rate of change transformation of "velocity" dB/dt =
    //                 Omega
    //                 // encoding:
    //                 // rotation: Omega_rot = (qx0, qy0, 1) * omega0
    //                 // translation: Omega_tra = (-v0.y, v0.x, 0)
    //                 //
    //                 // resulting Omega = Omega_rot + Omega_tra

    //                 // case with initial translation
    //                 u[n, 1] = vec2dp(-cm_w_spd0.y, cm_w_spd0.x,
    //                 cm_w_omega0);
    //             }
    //         }

    //         void calc_rhs()
    //         {

    //             // Create mdspan views for rhs calculation
    //             auto u = mdspan<vec2dp, dextents<size_t, 2>>(u_mem.data(), npts, 2);
    //             auto rhs = mdspan<vec2dp, dextents<size_t, 2>>(rhs_mem.data(), npts,
    //             2);

    //             for (size_t n = 0; n < npts; ++n) {

    //                 // get current state
    //                 vec2dp B = u[n, 0];     // position transformation B is in u[n,0]
    //                 vec2dp Omega = u[n, 1]; // velocity trafo d(B)/dt = Omega is in
    //                 u[n,1]

    //                 // forces and torques to change linear and angular acceleration

    //                 // no force f acting the body frame at the origin (F_b = O_2dp ^ f)
    //                 // auto F_b = bivec2dp{0.0, 0.0, 0.0};

    //                 // gravitation acting through center of mass in in the body system
    //                 // with center of mass being the origin O_2dp of the body system
    //                 // gravitation f = (x=0, y=-m*g, z=0) acting at O_2dp
    //                 // where the y-axis is assumed pointing upwards
    //                 auto F_b = wdg(O_2dp, vec2dp{0.0, -m[n] * 9.81, 0.0});

    //                 // Set right-hand side for ODE system:
    //                 // u[n,0]' = velocity trafo Omega = d(B)/dt (linear and angular)
    //                 // u[n,1]' = acceleration trafo d(Omega)/dt (linear and angular)
    //                 rhs[n, 0] = Omega;
    //                 rhs[n, 1] = compute_omega_dot(I_inv, F_b, Omega, I);
    //             }
    //         }

    //         void calc_rkstep(double dt)
    //         {
    //             // Create mdspan views for RK4 integration (for all npts points)
    //             auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), npts * 2);
    //             auto uh = mdspan<vec2dp, dextents<size_t, 2>>(uh_mem.data(), 2, npts *
    //             2); auto rhs =
    //                 mdspan<vec2dp const, dextents<size_t, 1>>(rhs_mem.data(), npts *
    //                 2);

    //             // Perform RK4 integration (4 sub-steps)
    //             for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
    //                 calc_rhs();
    //                 rk4_step(u, uh, rhs, dt, rk_step);
    //             }
    //         }

    //         void print_sim(double t)
    //         {
    //             // Create mdspan view for printing
    //             auto u = mdspan<vec2dp, dextents<size_t, 2>>(u_mem.data(), npts, 2);


    //             fmt::println("t = {:>-7.3f}:", t);

    //             for (size_t n = 0; n < npts; ++n) {

    //                 // get current state (= current positional transformation bivector)
    //                 vec2dp B = u[n, 0];     // B = Omega * t + B0 (from integration)
    //                 vec2dp Omega = u[n, 1]; // dB/dt = Omega = dB^2/dt^2 * t + Omega0

    //                 // calculate current position from B via M = rexp(0.5 * B) ⟇ M0
    //                 // and via pts(t) = M ⟇ pts(t0) ⟇ rrev(M)
    //                 //
    //                 auto M = rexp(0.5 * B);
    //                 fmt::println(
    //                     "    n = {}, B = {:>-7.3f}, Omega = {:>-7.3f}, M = {:>-7.3f}",
    //                     n, B, Omega, M);
    //                 fmt::println("    n = {}, pts(t=0) = {:>-7.3f}, pts(t) =
    //                 {:>-7.3f}",
    //                              n, pts[n], gr1(rgpr(rgpr(M, pts[n]), rrev(M))));
    //                 fmt::println("");
    //             }
    //         }

    //       private:

    //         // number of independent points to solve system for
    //         size_t npts;

    //         // positions of points in the body frame (fixed for now to get started)
    //         std::vector<vec2dp> pts; // point positions (fixed for now)
    //         std::vector<value_t> m;  // mass of points

    //         // initial position and speed of center of mass relative to world system
    //         vec2dp cm_w_pos0; // initial position of center of mass (unitized)
    //         vec2dp cm_w_spd0; // initial speed of center of mass (.z == 0)

    //         // initial angular position of cm (body vs. world frame, e12 defines pos
    //         phi) value_t cm_w_phi0;

    //         // initial angular velocity of cm (body vs. world frame)
    //         value_t cm_w_omega0;

    //         // RK4 integration state for point n with system order = 2
    //         // => [n+0: position, n+1: velocity])
    //         std::vector<vec2dp> u_mem;   // [n+0: position, n+1: velocity]
    //         std::vector<vec2dp> uh_mem;  // helper for integration
    //         std::vector<vec2dp> rhs_mem; // right-hand side values

    //         // inertia map and its inverse (calculated from descrete input values)
    //         // calculation is done in body frame
    //         inertia2dp I;
    //         inertia2dp I_inv;
    //     };

    //     auto constexpr num_points = 2;
    //     // (initial) positions and masses of each mass point in the body frame
    //     std::vector<vec2dp> pos_b = {vec2dp{0, 1, 1}, vec2dp{0, -1, 1}}; // position
    //     std::vector<value_t> m = {1.0, 1.0};                             // mass

    //     // initial position and movement of body frame vs. world frame (in world frame)
    //     auto cm_pos = O_2dp;           // world and body frame coincide at t=0
    //     auto cm_spd = vec2dp{1, 0, 0}; // body frame is moving with v0 vs. world frame

    //     // initial rotation and anguar velocity of body frame vs. world frame
    //     auto cm_phi = 0.0;   // no initial rotation
    //     auto cm_omega = 0.0; // no initial angular velocity

    //     sim_ode_multi_body_pga2dp sim(num_points, pos_b, m, cm_pos, cm_spd, cm_phi,
    //     cm_omega);

    //     // time range (from, to, number of steps)
    //     auto t_rng = discrete_range(0.0, 2.0, 20);

    //     sim.set_initial_values();
    //     sim.print_sim(t_rng.min());

    //     for (size_t n = 1; n <= t_rng.steps(); n++) {

    //         // integration from t to t + dt
    //         double t = t_rng.min() + n * t_rng.delta();
    //         sim.calc_rkstep(t_rng.delta());

    //         // print sim status at t+dt
    //         sim.print_sim(t);
    //     }

    //     fmt::println("");
    // }

    /////////////////////////////////////////////////////////////////////////////////////
    // dynamic_system2dp -- Phase B, Milestone 1: free rigid body under gravity.
    // Validation criterion (per the plan): energy conservation. The cm of a 2D rigid body
    // in a uniform field moves as a point mass (parabola, decoupled from spin), and with
    // no torque about the cm the angular velocity is conserved.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga2dp: dynamic_system2dp - free rigid body energy conservation (M1)")
    {
        fmt::println(
            "pga2dp: dynamic_system2dp - free rigid body energy conservation (M1)");

        // a uniform 1x1 plate (mass 2) launched from (0,5) with velocity (1,0) and spin 2
        dynamic_system2dp sys;
        sys.add_frame(static_frame2dp("W")); // inertial root
        sys.add_body(static_frame2dp("B", vec2dp{0.0, 5.0, 1.0}, 0.0),
                     make_plate_body(2.0, 1.0, 1.0),
                     kin_state2dp{.vel = vec2dp{1.0, 0.0, 0.0}, .omega = 2.0});

        value_t const dt = 0.0005;
        size_t const N = 2000; // -> t = 1.0 s

        value_t const E0 = sys.total_energy();
        value_t Emin = E0, Emax = E0;
        for (size_t n = 0; n < N; ++n) {
            sys.step(dt);
            value_t const E = sys.total_energy();
            if (E < Emin) Emin = E;
            if (E > Emax) Emax = E;
        }
        value_t const drift = (Emax - Emin) / std::abs(E0);

        // 1. energy conservation (the criterion): RK4 drift must be tiny
        CHECK(drift < 1.0e-6);

        // 2. cm follows the projectile parabola (decoupled from spin):
        //    x(1) = 0 + 1*1 = 1 ;  y(1) = 5 - 1/2 * 9.81 * 1^2 = 0.095
        vec2dp const cm = move2dp(O_2dp, sys.get_pos_trafo("B", "W"));
        CHECK(cm.x == doctest::Approx(1.0).epsilon(1e-4));
        CHECK(cm.y == doctest::Approx(5.0 - 0.5 * 9.81).epsilon(1e-4));

        // 3. no torque about the cm -> angular velocity is conserved
        CHECK(sys.twist_world("B").z == doctest::Approx(2.0).epsilon(1e-6));

        fmt::println("  E0 = {:.6f}, dE/E = {:.2e}, cm = {:>-7.3f}, omega = {:.4f}", E0,
                     drift, cm, sys.twist_world("B").z);
        fmt::println("");
    }

    TEST_CASE("pga2dp: free body integration converges at FOURTH order on the group")
    {
        fmt::println("pga2dp: dynamic_system2dp - free body convergence order (dexp^-1)");

        // Same gate as the 3D case: the generator B and the body twist Omega do not
        // commute for a spinning body whose body-frame velocity rotates, so without the
        // dexp^-1 brackets in step_free_body the scheme is second order on the pose
        // (measured ratio 4 per halving of dt; with both brackets 16). Gravity on, so the
        // body-frame velocity changes as well. Reference at 32x the finest resolution.
        auto make = [] {
            dynamic_system2dp sys;
            sys.add_frame(static_frame2dp("W"));
            sys.add_body(static_frame2dp("B", vec2dp{0.0, 5.0, 1.0}, 0.0),
                         make_plate_body(2.0, 1.0, 0.3),
                         kin_state2dp{.vel = vec2dp{3.0, 1.0, 0.0}, .omega = 4.0});
            return sys;
        };
        value_t const T = 1.0;
        auto const nrm3 = [](vec2dp const& D) {
            return std::sqrt(D.x * D.x + D.y * D.y + D.z * D.z);
        };

        auto ref = make();
        for (size_t i = 0; i < 19200; ++i)
            ref.step(T / 19200.0);
        auto const Mref = ref.get_pos_trafo(1, 0);

        value_t err_prev = 0.0;
        for (size_t N : {150, 300, 600}) {
            auto sys = make();
            for (size_t i = 0; i < N; ++i)
                sys.step(T / value_t(N));
            value_t const err = nrm3(rlog(rgpr(rrev(Mref), sys.get_pos_trafo(1, 0))));
            fmt::println("  N = {:4}: pose err = {:.3e} (ratio {:5.2f})", N, err,
                         err_prev / err);
            if (err_prev > 0.0) CHECK(err_prev / err > 14.0); // fourth order: 16
            err_prev = err;
        }
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // dynamic_system2dp -- Phase B, Milestone 2: single pivoted plate (revolute joint).
    // The plate hangs from a body-fixed hinge Q_b and swings under gravity (a compound
    // pendulum). Validation: (1) the initial angular acceleration equals the analytic
    // tau/I_hinge, (2) energy is conserved across several swings.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga2dp: dynamic_system2dp - pivoted plate pendulum (M2)")
    {
        fmt::println("pga2dp: dynamic_system2dp - pivoted plate pendulum (M2)");

        // square plate w=h=2, mass 1; cm at rest at the world origin; revolute hinge at
        // the body corner Q_b = (1,1) -> world hinge (1,1). Released from rest (phi0 =
        // 0).
        value_t const m = 1.0, w = 2.0, h = 2.0;
        vec2dp const Q_b{1.0, 1.0, 1.0};

        dynamic_system2dp sys;
        sys.add_frame(static_frame2dp("W")); // inertial root
        sys.add_revolute_body(static_frame2dp("B", vec2dp{0.0, 0.0, 1.0}, 0.0),
                              make_plate_body(m, w, h), Q_b, /*phi0*/ 0.0,
                              /*omega0*/ 0.0);
        size_t const B = sys.index_of("B");

        // 1. initial angular acceleration vs the analytic compound-pendulum value
        //    I_hinge = m(w^2+h^2)/12 + m|Q_b|^2 = 8/12 + 2 = 2.66667
        //    tau about hinge (1,1) from gravity at cm (0,0): r=(-1,-1) -> tau_z = 9.81
        value_t const I_hinge = m * (w * w + h * h) / 12.0 + m * (1.0 + 1.0);
        value_t const alpha0 = 9.81 / I_hinge;
        value_t const alpha_init = sys.joint_accel(B);
        CHECK(alpha_init == doctest::Approx(alpha0).epsilon(1e-12));

        // 2. energy conservation over several swing periods
        //    T = 2*pi*rsqrt(I_hinge / (m g d)), d = |Q_b| = rsqrt(2)  ->  T ~ 2.76 s
        value_t const dt = 0.0005;
        size_t const N = 8000; // 4 s (~1.5 periods)
        value_t const E0 = sys.total_energy();
        value_t Emin = E0, Emax = E0, KEmax = 0.0;
        for (size_t n = 0; n < N; ++n) {
            sys.step(dt);
            value_t const E = sys.total_energy();
            if (E < Emin) Emin = E;
            if (E > Emax) Emax = E;
            value_t const ke = sys.kinetic_energy();
            if (ke > KEmax) KEmax = ke;
        }
        value_t const drift = (Emax - Emin) / KEmax; // E0 ~ 0 here, so scale by KEmax
        CHECK(drift < 1e-5);

        // 3. the lowest point (cm directly below the hinge) reaches the energy-predicted
        //    kinetic energy: KE_low = m g dh, dh = |Q_b| - 1 = rsqrt(2) - 1
        CHECK(KEmax == doctest::Approx(m * 9.81 * (std::sqrt(2.0) - 1.0)).epsilon(1e-3));

        // 4. the relative twist is exactly the hinge rotation (encoding consistency)
        CHECK(sys.twist_world("B").z ==
              doctest::Approx(sys.joint_omega(B)).epsilon(1e-12));

        fmt::println(
            "  alpha_init = {:.5f} (analytic {:.5f}), dE/KEmax = {:.2e}, KEmax = {:.4f}",
            alpha_init, alpha0, drift, KEmax);
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // dynamic_system2dp -- Phase B, Milestone 3: the DOUBLE PENDULUM (two coupled
    // revolute plates). Plate 2 is hinged to plate 1's opposite corner. This is the
    // coupled articulated case, solved by the joint-space forward dynamics. The motion is
    // chaotic (no analytic trajectory), so correctness is checked by: (1) the joint-space
    // mass matrix being consistent with the kinetic energy, and (2) energy conservation
    // across the swing.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga2dp: dynamic_system2dp - double pendulum energy conservation (M3)")
    {
        fmt::println(
            "pga2dp: dynamic_system2dp - double pendulum energy conservation (M3)");

        // two square plates (w=h=2, m=1). Plate 1 hinged to the world at its corner
        // Q=(1,1) -> world pivot (1,1). Plate 2 hinged to plate 1's opposite corner; its
        // own connecting corner is also Q=(1,1). Rest poses (relative to parent): cm1 at
        // (0,0), cm2 at (-2,-2). Released from an energetic, asymmetric state.
        value_t const m = 1.0, w = 2.0, h = 2.0;
        vec2dp const Q{1.0, 1.0, 1.0};

        dynamic_system2dp sys;
        sys.add_frame(static_frame2dp("W")); // inertial root
        sys.add_revolute_body(static_frame2dp("B1", vec2dp{0.0, 0.0, 1.0}, 0.0),
                              make_plate_body(m, w, h), Q, /*phi1*/ 2.0, /*omega1*/ 1.0);
        sys.add_revolute_body(static_frame2dp("B2", vec2dp{-2.0, -2.0, 1.0}, 0.0),
                              make_plate_body(m, w, h), Q, /*phi2*/ -1.5, /*omega2*/ -2.0,
                              sys.index_of("B1"));
        size_t const B1 = sys.index_of("B1"), B2 = sys.index_of("B2");

        // 1. joint-space mass matrix is consistent with the kinetic energy:
        //    1/2 * qdot^T M(q) qdot == kinetic_energy()   (validates the
        //    inertia/Jacobian)
        std::vector<value_t> const Mm = sys.mass_matrix();
        value_t const w1 = sys.joint_omega(B1), w2 = sys.joint_omega(B2);
        value_t const ke_M =
            0.5 * (w1 * (Mm[0] * w1 + Mm[1] * w2) + w2 * (Mm[2] * w1 + Mm[3] * w2));
        CHECK(ke_M == doctest::Approx(sys.kinetic_energy()).epsilon(1e-12));

        // 2. energy conservation across the coupled (chaotic) swing
        value_t const dt = 0.0002;
        size_t const N = 25000; // 5 s
        value_t const E0 = sys.total_energy();
        value_t Emin = E0, Emax = E0, KEmax = sys.kinetic_energy();
        for (size_t n = 0; n < N; ++n) {
            sys.step(dt);
            value_t const E = sys.total_energy();
            if (E < Emin) Emin = E;
            if (E > Emax) Emax = E;
            value_t const ke = sys.kinetic_energy();
            if (ke > KEmax) KEmax = ke;
        }
        value_t scale = KEmax;
        if (std::abs(E0) > scale) scale = std::abs(E0);
        value_t const drift = (Emax - Emin) / scale;
        CHECK(drift < 1e-4);

        fmt::println(
            "  E0 = {:.6f}, KEmax = {:.4f}, dE/scale = {:.2e}, phi=({:.3f},{:.3f})", E0,
            KEmax, drift, sys.joint_phi(B1), sys.joint_phi(B2));
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // dynamic_system2dp -- PRISMATIC joint: the translational DoF, to show the PGA
    // unification. A prismatic slider runs through the SAME machinery as the revolute
    // hinge: only the joint screw generator differs (an ideal point / direction with z=0
    // instead of a finite point with z=1). A body on a frictionless 45-degree rail under
    // gravity should accelerate at -g*sin(45) along the rail, with energy conserved.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga2dp: dynamic_system2dp - prismatic slider (unification of DoF)")
    {
        fmt::println("pga2dp: dynamic_system2dp - prismatic slider (unification of DoF)");

        value_t const m = 2.0;
        value_t const theta = pi / 4.0;
        vec2dp const dir{std::cos(theta), std::sin(theta), 0.0}; // unit rail direction

        dynamic_system2dp sys;
        sys.add_frame(static_frame2dp("W"));
        sys.add_prismatic_body(static_frame2dp("B", vec2dp{0.0, 0.0, 1.0}, 0.0),
                               make_plate_body(m, 1.0, 1.0), dir, /*s0*/ 0.0, /*v0*/ 0.0);
        size_t const B = sys.index_of("B");

        // 1. slide acceleration = component of gravity along the rail = -g sin(theta)
        CHECK(sys.joint_accel(B) ==
              doctest::Approx(-9.81 * std::sin(theta)).epsilon(1e-12));

        // 2. the joint-space mass for a PURE translation is just the mass (the angular
        //    I_cm term drops out because the prismatic joint screw has S.z = 0)
        CHECK(sys.mass_matrix()[0] == doctest::Approx(m).epsilon(1e-12));

        // 3. energy conservation as it slides (constant accel -> RK4 exact)
        value_t const dt = 0.0005;
        size_t const N = 4000; // 2 s
        value_t const E0 = sys.total_energy();
        value_t Emin = E0, Emax = E0, KEmax = 0.0;
        for (size_t n = 0; n < N; ++n) {
            sys.step(dt);
            value_t const E = sys.total_energy();
            if (E < Emin) Emin = E;
            if (E > Emax) Emax = E;
            value_t const ke = sys.kinetic_energy();
            if (ke > KEmax) KEmax = ke;
        }
        CHECK((Emax - Emin) / KEmax < 1e-12);

        fmt::println(
            "  a_slide = {:.5f} (analytic {:.5f}), M[0]={:.4f}, dE/KEmax = {:.2e}",
            sys.joint_accel(B), -9.81 * std::sin(theta), sys.mass_matrix()[0],
            (Emax - Emin) / KEmax);
        fmt::println("");
    }

    TEST_CASE("pga2dp: dynamic_system2dp - spring/damper damped oscillator (Phase A)")
    {
        fmt::println("pga2dp: dynamic_system2dp - joint spring/damper (Phase A)");

        // 2D twin of the 3D check: a plate on a horizontal prismatic joint carrying a
        // linear spring (k) + damper (c) on its coordinate, no gravity. The reduced
        // inertia is the pure mass (M[0] = m), so the joint obeys m q'' + c q' + k q = 0
        // -- a damped harmonic oscillator with a closed-form solution. Validates the
        // additive spring/damper force elements (set_joint_spring_damper) folded into
        // tau.
        value_t const m = 2.0, k = 8.0, c = 0.8;
        value_t const wn = std::sqrt(k / m);                  // natural frequency 2 rad/s
        value_t const zeta = c / (2.0 * std::sqrt(k * m));    // damping ratio 0.1
        value_t const wd = wn * std::sqrt(1.0 - zeta * zeta); // damped frequency
        value_t const A0 = 0.5;                               // initial q, zero velocity

        dynamic_system2dp sys;
        sys.set_gravity(vec2dp{0.0, 0.0, 0.0}); // isolate the spring oscillator
        sys.add_frame(static_frame2dp("W"));
        sys.add_prismatic_body(static_frame2dp("B", vec2dp{0.0, 0.0, 1.0}, 0.0),
                               make_plate_body(m, 1.0, 1.0), vec2dp{1.0, 0.0, 0.0}, A0,
                               0.0); // slide along e1; q0 = A0, v0 = 0
        size_t const B = sys.index_of("B");
        sys.set_joint_spring_damper(B, k, c); // q_rest = 0

        // reduced inertia is the pure mass; initial accel = -k*A0/m (spring only, v = 0)
        CHECK(sys.mass_matrix()[0] == doctest::Approx(m));
        CHECK(sys.joint_accel(B) == doctest::Approx(-k * A0 / m));

        // closed-form underdamped free response (q(0) = A0, q'(0) = 0):
        //   q(t) = e^{-zeta wn t} [ A0 cos(wd t) + (zeta wn A0 / wd) sin(wd t) ]
        auto q_exact = [&](value_t t) {
            return std::exp(-zeta * wn * t) *
                   (A0 * std::cos(wd * t) + (zeta * wn * A0 / wd) * std::sin(wd * t));
        };

        value_t const dt = 1.0e-4;
        value_t t = 0.0, max_err = 0.0;
        for (size_t nstep = 1; nstep <= 50000; ++nstep) { // 5 s, ~1.6 damped periods
            sys.step(dt);
            t += dt;
            max_err = std::max(max_err, std::abs(sys.joint_phi(B) - q_exact(t)));
        }
        fmt::println(
            "  wn = {:.4f}, zeta = {:.3f}, wd = {:.4f}; max|q - q_exact| = {:.2e}", wn,
            zeta, wd, max_err);
        CHECK(max_err < 1e-9); // RK4 matches the closed form to roundoff

        // undamped variant: KE + spring potential 1/2 k q^2 conserved (spring PE term in
        // potential_energy)
        dynamic_system2dp usys;
        usys.set_gravity(vec2dp{0.0, 0.0, 0.0});
        usys.add_frame(static_frame2dp("W"));
        usys.add_prismatic_body(static_frame2dp("B", vec2dp{0.0, 0.0, 1.0}, 0.0),
                                make_plate_body(m, 1.0, 1.0), vec2dp{1.0, 0.0, 0.0}, A0,
                                0.0);
        usys.set_joint_spring_damper(usys.index_of("B"), k, 0.0);
        value_t const E0 = usys.total_energy(); // = 1/2 k A0^2
        CHECK(E0 == doctest::Approx(0.5 * k * A0 * A0));
        value_t max_dE = 0.0, kemax = 0.0;
        for (size_t nstep = 1; nstep <= 40000; ++nstep) { // 4 s
            usys.step(dt);
            max_dE = std::max(max_dE, std::abs(usys.total_energy() - E0));
            kemax = std::max(kemax, usys.kinetic_energy());
        }
        fmt::println("  undamped: E0 = {:.5f}, KEmax = {:.5f}, dE/KEmax = {:.2e}", E0,
                     kemax, max_dE / kemax);
        CHECK(max_dE / kemax < 1e-11);

        fmt::println("");
    }

    TEST_CASE("pga2dp: applied wrench + time threading - forced oscillator (Phase A.2)")
    {
        fmt::println("pga2dp: dynamic_system2dp - applied wrench / forced (Phase A.2)");

        // 2D twin of the Phase-A.2 check: a time-varying applied wrench
        // (set_applied_wrench)
        // + RK4 sub-step time threading drive a single prismatic joint (mass mt) with a
        // spring/damper, giving the forced damped oscillator mt q'' + R q' + K q =
        // F0 cos(Omega t). The steady-state amplitude matches the closed form
        // A = F0 / rsqrt((Omega R)^2 + (K - mt Omega^2)^2) -- B.1's Eq. (4) with F0 = m e
        // Omega^2 (the 2D radial plane is e1-e2).
        value_t const mt = 7.35, K = 2000.0, R = 5.0;
        value_t const me = 4.9 * 0.008336; // unbalance m*e [kg m]
        vec2dp const e1{1.0, 0.0, 0.0};

        // sign/magnitude check: a CONSTANT force F0 along e1 gives joint accel F0/mt
        {
            dynamic_system2dp s;
            s.set_gravity(vec2dp{0.0, 0.0, 0.0});
            s.add_frame(static_frame2dp("W"));
            s.add_prismatic_body(static_frame2dp("B", vec2dp{0.0, 0.0, 1.0}, 0.0),
                                 make_plate_body(mt, 1.0, 1.0), e1);
            value_t const F0 = 3.0;
            s.set_applied_wrench(
                1, [F0](value_t) { return wdg(O_2dp, vec2dp{F0, 0.0, 0.0}); });
            CHECK(s.joint_accel(1) == doctest::Approx(F0 / mt));
        }

        // forced sweep at the resonance and just off it; compare to the closed form
        auto amp_cf = [&](value_t Om) {
            value_t const d = K - mt * Om * Om;
            return me * Om * Om / std::sqrt((Om * R) * (Om * R) + d * d);
        };
        auto run = [&](value_t Om) {
            dynamic_system2dp s;
            s.set_gravity(vec2dp{0.0, 0.0, 0.0});
            s.add_frame(static_frame2dp("W"));
            s.add_prismatic_body(static_frame2dp("B", vec2dp{0.0, 0.0, 1.0}, 0.0),
                                 make_plate_body(mt, 1.0, 1.0), e1);
            s.set_joint_spring_damper(1, K, R);
            s.set_applied_wrench(1, [me, Om](value_t t) {
                return wdg(O_2dp, vec2dp{me * Om * Om * std::cos(Om * t), 0.0, 0.0});
            });
            value_t const dt = 1.0e-3, t_end = 28.0, t_meas = 22.0;
            value_t qmin = std::numeric_limits<value_t>::max();
            value_t qmax = std::numeric_limits<value_t>::lowest();
            for (value_t t = 0.0; t < t_end; t += dt) {
                s.step(dt);
                if (s.time() >= t_meas) {
                    qmin = std::min(qmin, s.joint_phi(1));
                    qmax = std::max(qmax, s.joint_phi(1));
                }
            }
            return 0.5 * (qmax - qmin);
        };

        value_t const om1 = std::sqrt(K / mt);
        for (value_t Om : {om1, 24.0}) {
            value_t const a = run(Om), cf = amp_cf(Om);
            fmt::println("  Omega = {:.3f}: A_lib = {:.4e}, A_cf = {:.4e}, err% = {:.3f}",
                         Om, a, cf, 100.0 * std::abs(a - cf) / cf);
            CHECK(a == doctest::Approx(cf).epsilon(0.01));
        }

        fmt::println("");
    }

    TEST_CASE(
        "pga2dp: driven joint - spinning radial slider (centrifugal eq., Phase A.3)")
    {
        fmt::println(
            "pga2dp: dynamic_system2dp - driven spin + radial slider (Phase A.3)");

        // 2D twin: a driven revolute base spins at constant Omega about the origin; a
        // radial prismatic slider (along the rotating e1) carries a mass m on a spring
        // (k) + damper (c) with rest at radius L. The mass settles where the spring
        // balances the centrifugal force m Omega^2 (L+q): q_eq = m Omega^2 L / (k - m
        // Omega^2). The centrifugal force EMERGES from the driven base's bias -- nothing
        // adds it by hand.
        value_t const m = 1.0, L = 1.0, Om = 2.0, k = 20.0, c = 2.0;
        auto const slider = make_plate_body(m, 0.05, 0.05); // near point mass
        dynamic_system2dp sys;
        sys.set_gravity(vec2dp{0.0, 0.0, 0.0});
        sys.add_frame(static_frame2dp("W"));
        sys.add_revolute_body(static_frame2dp("S"), make_plate_body(1.0, 0.1, 0.1),
                              O_2dp); // spin hinge about the origin
        sys.add_prismatic_body(static_frame2dp("B", vec2dp{L, 0.0, 1.0}, 0.0), slider,
                               vec2dp{1.0, 0.0, 0.0}); // radial slider along e1
        size_t const B = sys.index_of("B");
        sys.set_joint_spring_damper(B, k, c);
        sys.set_driven_rate(sys.index_of("S"), Om);

        value_t const dt = 5.0e-4;
        for (int n = 0; n < 20000; ++n)
            sys.step(dt); // 10 s -> steady state
        value_t const q_eq = m * Om * Om * L / (k - m * Om * Om);
        fmt::println("  q_settled = {:.5f} (analytic q_eq = {:.5f})", sys.joint_phi(B),
                     q_eq);
        CHECK(sys.joint_phi(B) == doctest::Approx(q_eq).epsilon(0.01));
        fmt::println("");
    }

    TEST_CASE("pga2dp: grounded spring/damper - stiffness emergence (Phase C.1)")
    {
        fmt::println(
            "pga2dp: dynamic_system2dp - grounded spring stiffness emergence (C.1)");

        // 2D twin of the 3D C.1 gate: the grounded spatial spring (add_grounded_spring)
        // yields a translational stiffness, and -- via its lever arm -- a tilt stiffness.

        // Gate 1: translational stiffness == k, free response == x0 cos(wn t).
        {
            value_t const m = 2.0, k = 800.0, x0 = 0.01;
            value_t const wn = std::sqrt(k / m);
            dynamic_system2dp sys;
            sys.set_gravity(vec2dp{0.0, 0.0, 0.0});
            sys.add_frame(static_frame2dp("W"));
            sys.add_prismatic_body(static_frame2dp("B"), make_plate_body(m, 0.1, 0.1),
                                   vec2dp{1.0, 0.0, 0.0}, x0, 0.0);
            sys.add_grounded_spring(1, O_2dp, O_2dp, vec2dp{k, 0.0, 0.0}, 0.0);

            value_t const k_emergent = -sys.mass_matrix()[0] * sys.joint_accel(1) / x0;
            fmt::println("  Gate1 translational: k_emergent = {:.4f} (k = {:.1f})",
                         k_emergent, k);
            CHECK(k_emergent == doctest::Approx(k));

            value_t const dt = 1.0e-4;
            value_t t = 0.0, max_err = 0.0;
            for (size_t nstep = 1; nstep <= 30000; ++nstep) {
                sys.step(dt);
                t += dt;
                max_err =
                    std::max(max_err, std::abs(sys.joint_phi(1) - x0 * std::cos(wn * t)));
            }
            fmt::println("  Gate1 free response: max|q - x0 cos(wn t)| = {:.2e}",
                         max_err);
            CHECK(max_err < 1e-9);
        }

        // Gate 2: tilt stiffness from the lever arm == k * l^2 -- a perpendicular spring
        // at radial offset l on a body that rotates about the origin gives k_theta = k
        // l^2.
        {
            value_t const k = 1000.0, l = 0.1, th0 = 1.0e-4;
            dynamic_system2dp sys;
            sys.set_gravity(vec2dp{0.0, 0.0, 0.0});
            sys.add_frame(static_frame2dp("W"));
            sys.add_revolute_body(static_frame2dp("B"), make_plate_body(1.0, 0.4, 0.05),
                                  O_2dp, th0, 0.0);
            sys.add_grounded_spring(1, vec2dp{l, 0.0, 1.0}, vec2dp{l, 0.0, 1.0},
                                    vec2dp{0.0, k, 0.0}, 0.0);

            value_t const ktheta_emergent =
                -sys.mass_matrix()[0] * sys.joint_accel(1) / th0;
            fmt::println("  Gate2 tilt: k_theta_emergent = {:.5f} (k l^2 = {:.5f})",
                         ktheta_emergent, k * l * l);
            CHECK(ktheta_emergent == doctest::Approx(k * l * l).epsilon(1e-3));
        }

        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // closed_loop_system2dp -- Phase 1 (position-level assembly): the planar FOUR-BAR
    // linkage, the canonical 1-DOF closed loop. As a spanning tree it is two branches off
    // the ground: crank -> coupler (open 3-link chain O2 -> A -> B), and the rocker (O4
    // -> B') as a separate branch. The loop is closed by a point-coincidence constraint
    // requiring the coupler tip B to coincide with the rocker tip B'. assemble() drives
    // the crank and solves the coupler + rocker joints so g(q) = 0 (Newton on the
    // constraint Jacobian, whose columns are the reused spatial-Jacobian velocity_field
    // partials). Validated against an analytically-known closed configuration.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga2dp: closed_loop_system2dp - four-bar assembly (Phase 1)")
    {
        fmt::println("pga2dp: closed_loop_system2dp - four-bar assembly (Phase 1)");

        // four-bar lengths chosen so a clean closed configuration exists at crank = pi/2:
        //   a (crank) = 2, b (coupler) = 3, c (rocker) = rsqrt(5), d (ground) = 4
        //   O2 = (0,0), O4 = (4,0); at theta2 = pi/2: A = (0,2), B = B' = (3,2)
        //   -> dependent angles theta3 = -pi/2 (coupler vs crank), theta4 = atan2(2,-1)
        value_t const a = 2.0, b = 3.0, c = std::sqrt(5.0), d = 4.0;
        value_t const theta2 = pi / 2.0;        // crank angle (driver, held fixed)
        value_t const theta3_exact = -pi / 2.0; // coupler vs crank (ground truth)
        value_t const theta4_exact = std::atan2(2.0, -1.0); // rocker (ground truth)

        // nominal body (its inertia/cm is irrelevant to the position-level assembly)
        auto const link = make_plate_body(1.0, 1.0, 1.0);

        closed_loop_system2dp cl;
        cl.add_frame(static_frame2dp("W")); // inertial root / ground

        // branch 1: crank (revolute about O2) then coupler (revolute about A)
        cl.add_revolute_body(static_frame2dp("CR", vec2dp{0.0, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, /*phi0*/ theta2, /*omega0*/ 0.0,
                             cl.index_of("W"));
        cl.add_revolute_body(static_frame2dp("CO", vec2dp{a, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, /*phi0 (perturbed)*/ -1.4,
                             /*omega0*/ 0.0, cl.index_of("CR"));
        // branch 2: rocker (revolute about O4)
        cl.add_revolute_body(static_frame2dp("RO", vec2dp{d, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, /*phi0 (perturbed)*/ 2.2,
                             /*omega0*/ 0.0, cl.index_of("W"));

        size_t const CR = cl.index_of("CR"), CO = cl.index_of("CO"),
                     RO = cl.index_of("RO");

        // close the loop: coupler tip B (= (b,0) in CO) coincides with rocker tip B'
        // (= (c,0) in RO)
        cl.add_loop_constraint(loop_constraint2dp{CO, vec2dp{b, 0.0, 1.0}, RO,
                                                  vec2dp{c, 0.0, 1.0},
                                                  constraint2dp::coincidence});

        // 1. the loop is initially OPEN (coupler/rocker tips do not coincide)
        value_t const g0 = cl.residual_norm();
        CHECK(g0 > 0.1);

        // 2. assemble: drive the crank, solve coupler + rocker (square Newton: 2
        // dependent
        //    joints, 2 coincidence equations) -> residual driven to ~0
        value_t const gfin = cl.assemble(/*driven*/ {CR});
        CHECK(gfin < 1e-12);
        CHECK(cl.residual_norm() < 1e-12);

        // 3. the solved dependent joint angles match the analytic closed configuration
        value_t const theta3_solved = cl.system().joint_phi(CO);
        value_t const theta4_solved = cl.system().joint_phi(RO);
        CHECK(theta3_solved == doctest::Approx(theta3_exact).epsilon(1e-9));
        CHECK(theta4_solved == doctest::Approx(theta4_exact).epsilon(1e-9));
        // the driven crank angle was held fixed
        CHECK(cl.system().joint_phi(CR) == doctest::Approx(theta2).epsilon(1e-15));

        // 4. the closure point B reached from both branches lands at (3,2)
        auto const Bc = move2dp(vec2dp{b, 0.0, 1.0}, cl.system().get_pos_trafo(CO, 0));
        auto const Br = move2dp(vec2dp{c, 0.0, 1.0}, cl.system().get_pos_trafo(RO, 0));
        CHECK(Bc.x / Bc.z == doctest::Approx(3.0).epsilon(1e-9));
        CHECK(Bc.y / Bc.z == doctest::Approx(2.0).epsilon(1e-9));
        CHECK(Br.x / Br.z == doctest::Approx(3.0).epsilon(1e-9));
        CHECK(Br.y / Br.z == doctest::Approx(2.0).epsilon(1e-9));

        // 5. sweep the crank over a range and re-assemble at each step: the loop stays
        //    closed (continuous assembly within one branch of the four-bar)
        value_t gmax = 0.0;
        for (int k = 1; k <= 20; ++k) {
            value_t const th = theta2 + 0.3 * std::sin(0.3 * k); // wiggle the crank
            cl.set_joint(CR, th);                                // drive the crank to th
            value_t const gk = cl.assemble(/*driven*/ {CR});
            gmax = std::max(gmax, gk);
        }
        CHECK(gmax < 1e-10);

        fmt::println(
            "  g_open = {:.4f} -> g_assembled = {:.2e}; theta3 = {:.6f} "
            "(exact {:.6f}), theta4 = {:.6f} (exact {:.6f}); g_sweep_max = {:.2e}",
            g0, gfin, theta3_solved, theta3_exact, theta4_solved, theta4_exact, gmax);
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // closed_loop_system2dp -- Phase 2 (kinematic closed loop): velocity and acceleration
    // distribution on the same planar four-bar. Driving the crank at a prescribed rate,
    // the dependent (coupler, rocker) joint rates follow from the velocity-level loop
    // closure G q-dot = 0, and their accelerations from the acceleration-level closure
    // G q-ddot = -G-dot q-dot. Validated against analytic four-bar ratios at the closed
    // configuration theta2 = pi/2 (A = (0,2), B = (3,2)):
    //
    //   velocity closure  v_B(crank+coupler) = v_B(rocker):
    //     omega4 = omega2,            omega3(rel) = -4/3 omega2
    //   acceleration closure (alpha2 = 0):
    //     alpha4 = 2/3 omega2^2,      alpha3(rel) = -2/9 omega2^2
    //
    // (derived by equating the two branch expressions for the velocity / acceleration of
    // the shared coupler tip B; omega3/alpha3 are the joint RELATIVE coupler-vs-crank
    // rates, which is what the solver returns.)
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga2dp: closed_loop_system2dp - four-bar velocity/acceleration (Phase 2)")
    {
        fmt::println(
            "pga2dp: closed_loop_system2dp - four-bar velocity/acceleration (Phase 2)");

        value_t const a = 2.0, b = 3.0, c = std::sqrt(5.0), d = 4.0;
        value_t const theta2 = pi / 2.0;
        auto const link = make_plate_body(1.0, 1.0, 1.0);

        closed_loop_system2dp cl;
        cl.add_frame(static_frame2dp("W"));
        cl.add_revolute_body(static_frame2dp("CR", vec2dp{0.0, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, theta2, 0.0, cl.index_of("W"));
        cl.add_revolute_body(static_frame2dp("CO", vec2dp{a, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, -1.4, 0.0, cl.index_of("CR"));
        cl.add_revolute_body(static_frame2dp("RO", vec2dp{d, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, 2.2, 0.0, cl.index_of("W"));
        size_t const CR = cl.index_of("CR"), CO = cl.index_of("CO"),
                     RO = cl.index_of("RO");
        cl.add_loop_constraint(loop_constraint2dp{CO, vec2dp{b, 0.0, 1.0}, RO,
                                                  vec2dp{c, 0.0, 1.0},
                                                  constraint2dp::coincidence});

        // assemble the closed configuration first (g ~ 0)
        cl.assemble(/*driven*/ {CR});
        REQUIRE(cl.residual_norm() < 1e-12);

        // world-coordinate (unitized) closure point reached from each branch
        auto unit = [](vec2dp const& p) { return vec2dp{p.x / p.z, p.y / p.z, 1.0}; };
        auto const Pa =
            unit(move2dp(vec2dp{b, 0.0, 1.0}, cl.system().get_pos_trafo(CO, 0)));
        auto const Pb =
            unit(move2dp(vec2dp{c, 0.0, 1.0}, cl.system().get_pos_trafo(RO, 0)));

        // --- velocity level: drive the crank at omega2 = 1, distribute to coupler/rocker
        value_t const omega2 = 1.0;
        cl.set_joint_rate(CR, omega2);
        cl.solve_velocities(/*driven*/ {CR});

        CHECK(cl.joint_rate(CO) == doctest::Approx(-4.0 / 3.0 * omega2).epsilon(1e-9));
        CHECK(cl.joint_rate(RO) == doctest::Approx(omega2).epsilon(1e-9));

        // velocity closure: the shared tip B has the SAME velocity from both branches
        auto const vA = cl.system().point_velocity(Pa, CO);
        auto const vB = cl.system().point_velocity(Pb, RO);
        CHECK(std::abs(vA.x - vB.x) < 1e-9);
        CHECK(std::abs(vA.y - vB.y) < 1e-9);

        // --- acceleration level: constant crank rate (alpha2 = 0) -> coupler/rocker
        // accel
        auto const adep = cl.solve_accelerations(/*driven*/ {CR}, /*alpha_drv*/ {0.0});
        // adep is in dependent-joint order = [CO, RO]
        CHECK(adep[0] == doctest::Approx(-2.0 / 9.0 * omega2 * omega2).epsilon(1e-9));
        CHECK(adep[1] == doctest::Approx(2.0 / 3.0 * omega2 * omega2).epsilon(1e-9));

        // acceleration closure: the shared tip B has the SAME acceleration from both
        // branches (the defining property G q-ddot + G-dot q-dot = 0)
        auto const aA = cl.system().point_acceleration(Pa, CO);
        auto const aB = cl.system().point_acceleration(Pb, RO);
        CHECK(std::abs(aA.x - aB.x) < 1e-9);
        CHECK(std::abs(aA.y - aB.y) < 1e-9);

        fmt::println("  omega: coupler = {:.6f} (exact {:.6f}), rocker = {:.6f} (exact "
                     "{:.6f})",
                     cl.joint_rate(CO), -4.0 / 3.0, cl.joint_rate(RO), 1.0);
        fmt::println("  alpha: coupler = {:.6f} (exact {:.6f}), rocker = {:.6f} (exact "
                     "{:.6f})",
                     adep[0], -2.0 / 9.0, adep[1], 2.0 / 3.0);
        fmt::println("  closure: |dv| = ({:.1e},{:.1e}), |da| = ({:.1e},{:.1e})",
                     std::abs(vA.x - vB.x), std::abs(vA.y - vB.y), std::abs(aA.x - aB.x),
                     std::abs(aA.y - aB.y));
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // closed_loop_system2dp -- Phase 3 (dynamic closed loop): constrained forward
    // dynamics of the four-bar under gravity. The three joints (crank, coupler, rocker)
    // are all free and coupled by the 2 coincidence equations (a 1-DOF mechanism),
    // integrated by RK4 with the acceleration-level KKT solve at each sub-step and
    // post-step projection stabilisation. The headline correctness metric is ENERGY
    // CONSERVATION: the constraint forces (Lagrange multipliers) do no work, so total
    // energy is conserved to integrator tolerance, and the closure error ‖g‖ stays
    // bounded (no drift). Each link is a plate whose centre of mass sits at its frame
    // origin (the joint), so the coupler mass at the moving crank tip A is the
    // gravitating element that drives the loop.
    /////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("pga2dp: closed_loop_system2dp - four-bar energy conservation (Phase 3)")
    {
        fmt::println(
            "pga2dp: closed_loop_system2dp - four-bar energy conservation (Phase 3)");

        value_t const a = 2.0, b = 3.0, c = std::sqrt(5.0), d = 4.0;
        // start the crank away from theta2 = pi/2 (where the coupler mass at A sits
        // directly above O2 and gravity is radial -> a zero-torque unstable equilibrium
        // that would never leave rest); 1.2 rad gives a real driving torque.
        value_t const theta2 = 1.2;
        auto const link = make_plate_body(1.0, 1.0, 1.0); // unit-mass plate, cm at origin

        closed_loop_system2dp cl;
        cl.add_frame(static_frame2dp("W"));
        cl.add_revolute_body(static_frame2dp("CR", vec2dp{0.0, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, theta2, 0.0, cl.index_of("W"));
        cl.add_revolute_body(static_frame2dp("CO", vec2dp{a, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, -1.4, 0.0, cl.index_of("CR"));
        cl.add_revolute_body(static_frame2dp("RO", vec2dp{d, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, 2.2, 0.0, cl.index_of("W"));
        size_t const CR = cl.index_of("CR");
        cl.add_loop_constraint(loop_constraint2dp{cl.index_of("CO"), vec2dp{b, 0.0, 1.0},
                                                  cl.index_of("RO"), vec2dp{c, 0.0, 1.0},
                                                  constraint2dp::coincidence});

        // start from the assembled configuration, released from rest (all rates 0)
        cl.assemble(/*driven*/ {CR});
        REQUIRE(cl.residual_norm() < 1e-12);

        // 1. constraint forces do no work -> energy conserved over the swing
        value_t const dt = 0.0005;
        size_t const N = 6000; // 3 s
        value_t const E0 = cl.system().total_energy();
        value_t Emin = E0, Emax = E0, KEmax = 0.0, gmax = 0.0;
        for (size_t n = 0; n < N; ++n) {
            cl.step(dt);
            value_t const E = cl.system().total_energy();
            Emin = std::min(Emin, E);
            Emax = std::max(Emax, E);
            KEmax = std::max(KEmax, cl.system().kinetic_energy());
            gmax = std::max(gmax, cl.residual_norm());
        }
        value_t scale = KEmax;
        if (std::abs(E0) > scale) scale = std::abs(E0);
        value_t const drift = (Emax - Emin) / scale;
        CHECK(drift <
              1e-4); // observed ~1e-6; same tolerance class as the open-loop tests

        // 2. the loop stays closed across the whole run (projection controls drift)
        CHECK(gmax < 1e-9);

        // 3. the mechanism actually moved (non-trivial dynamics, not a frozen config)
        CHECK(KEmax > 0.05);

        fmt::println(
            "  E0 = {:.6f}, KEmax = {:.4f}, dE/scale = {:.2e}, max||g|| = {:.2e}", E0,
            KEmax, drift, gmax);
        fmt::println("");
    }

    /////////////////////////////////////////////////////////////////////////////////////
    // closed_loop_system2dp -- planar 5-bar (2-RRR) parallel manipulator: the 2D analogue
    // of a delta / Stewart-Gough robot, and the mechanism behind the active_planar_delta
    // ga_view scene. Two 2-link arms rise from two fixed, actuated shoulders (SA, SB);
    // the forearm tips are pinned to a shared end-effector by ONE point-coincidence. 4
    // revolute joints (2 shoulders + 2 elbows) - 2 coincidence equations -> 2 DOF. Unlike
    // the single-input four-bar, BOTH shoulders are driven: assemble({SA, SB}) solves the
    // two dependent elbows (forward kinematics of the closed loop). The parameters mirror
    // the ga_view scene defaults, so this also guards that the demo stays in the
    // well-conditioned interior of the workspace (away from elbow full-extension
    // singularities) across the whole driven sweep.
    /////////////////////////////////////////////////////////////////////////////////////
    TEST_CASE("pga2dp: closed_loop_system2dp - planar 5-bar parallel manipulator")
    {
        fmt::println("pga2dp: closed_loop_system2dp - planar 5-bar parallel manipulator");

        // mirrors planar_delta_params (ga_view active_planar_delta scene defaults)
        vec2dp const base_a{-1.0, -0.6, 1.0}, base_b{1.0, -0.6, 1.0};
        value_t const upper = 1.5, fore = 1.8;
        value_t const sa0 = 1.9198622, sb0 = 1.2217305; // 110 / 70 deg
        value_t const amp = 0.25, phase = pi / 2.0, period = 5.0;

        // upper circle-circle intersection -> consistent start configuration
        auto circ2 = [](vec2dp const& c0, value_t r0, vec2dp const& c1,
                        value_t r1) -> vec2dp {
            value_t const dx = c1.x - c0.x, dy = c1.y - c0.y;
            value_t const dd = std::hypot(dx, dy);
            value_t const aa = (r0 * r0 - r1 * r1 + dd * dd) / (2.0 * dd);
            value_t const h2 = r0 * r0 - aa * aa;
            value_t const h = h2 > 0.0 ? std::sqrt(h2) : 0.0;
            value_t const px = c0.x + aa * dx / dd, py = c0.y + aa * dy / dd;
            return vec2dp{px - h * dy / dd, py + h * dx / dd, 1.0};
        };

        auto const link = make_plate_body(1.0, 1.0, 1.0);

        closed_loop_system2dp cl;
        cl.add_frame(static_frame2dp("W")); // inertial root / ground

        vec2dp const EA0{base_a.x + upper * std::cos(sa0),
                         base_a.y + upper * std::sin(sa0), 1.0};
        vec2dp const EB0{base_b.x + upper * std::cos(sb0),
                         base_b.y + upper * std::sin(sb0), 1.0};
        vec2dp const tip = circ2(EA0, fore, EB0, fore);
        value_t const tA0 = std::atan2(tip.y - EA0.y, tip.x - EA0.x) - sa0;
        value_t const tB0 = std::atan2(tip.y - EB0.y, tip.x - EB0.x) - sb0;

        // arm A: shoulder SA (revolute about base_a) then elbow EA
        cl.add_revolute_body(static_frame2dp("SA", base_a, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, /*phi0*/ sa0, /*omega0*/ 0.0,
                             cl.index_of("W"));
        cl.add_revolute_body(static_frame2dp("EA", vec2dp{upper, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, /*phi0*/ tA0, /*omega0*/ 0.0,
                             cl.index_of("SA"));
        // arm B: shoulder SB (revolute about base_b) then elbow EB
        cl.add_revolute_body(static_frame2dp("SB", base_b, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, /*phi0*/ sb0, /*omega0*/ 0.0,
                             cl.index_of("W"));
        cl.add_revolute_body(static_frame2dp("EB", vec2dp{upper, 0.0, 1.0}, 0.0), link,
                             vec2dp{0.0, 0.0, 1.0}, /*phi0*/ tB0, /*omega0*/ 0.0,
                             cl.index_of("SB"));

        size_t const SA = cl.index_of("SA"), EA = cl.index_of("EA"),
                     SB = cl.index_of("SB"), EB = cl.index_of("EB");

        // close the loop: forearm tip of arm A coincides with forearm tip of arm B
        cl.add_loop_constraint(loop_constraint2dp{EA, vec2dp{fore, 0.0, 1.0}, EB,
                                                  vec2dp{fore, 0.0, 1.0},
                                                  constraint2dp::coincidence});

        // 1. assemble the start configuration (the closed-form guess is already
        // consistent,
        //    so the residual starts small and is polished to ~0)
        value_t const gfin = cl.assemble(/*driven*/ {SA, SB});
        CHECK(gfin < 1e-12);
        CHECK(cl.residual_norm() < 1e-12);

        // 2. drive BOTH shoulders with phase-shifted sinusoids and re-assemble each step:
        //    the loop stays closed and the two forearm tips coincide throughout
        value_t gmax = 0.0, tip_mismatch_max = 0.0;
        value_t Px_min = tip.x, Px_max = tip.x, Py_min = tip.y, Py_max = tip.y;
        value_t const dt = 0.016;
        for (value_t t = dt; t <= 12.0; t += dt) {
            cl.set_joint(SA, sa0 + amp * std::sin(2.0 * pi * t / period));
            cl.set_joint(SB, sb0 + amp * std::sin(2.0 * pi * t / period + phase));
            value_t const gk = cl.assemble(/*driven*/ {SA, SB});
            gmax = std::max(gmax, gk);

            // effector reached from each arm independently (should coincide)
            auto Pa = move2dp(vec2dp{fore, 0.0, 1.0}, cl.system().get_pos_trafo(EA, 0));
            auto Pb = move2dp(vec2dp{fore, 0.0, 1.0}, cl.system().get_pos_trafo(EB, 0));
            value_t const ax = Pa.x / Pa.z, ay = Pa.y / Pa.z;
            value_t const bx = Pb.x / Pb.z, by = Pb.y / Pb.z;
            tip_mismatch_max = std::max(tip_mismatch_max, std::hypot(ax - bx, ay - by));
            Px_min = std::min(Px_min, ax);
            Px_max = std::max(Px_max, ax);
            Py_min = std::min(Py_min, ay);
            Py_max = std::max(Py_max, ay);
        }
        CHECK(gmax < 1e-10);            // well-conditioned across the whole driven sweep
        CHECK(tip_mismatch_max < 1e-9); // the shared effector is single-valued

        // 3. the effector actually traces a 2D curve (both shoulders contribute -> the
        //    motion is genuinely 2-DOF, not a degenerate 1-DOF line)
        value_t const span_x = Px_max - Px_min, span_y = Py_max - Py_min;
        CHECK(span_x > 0.2);
        CHECK(span_y > 0.2);

        fmt::println("  start g = {:.2e}; driven-sweep max||g|| = {:.2e}, max tip "
                     "mismatch = {:.2e}; effector span = ({:.3f} x {:.3f})",
                     gfin, gmax, tip_mismatch_max, span_x, span_y);
        fmt::println("");
    }

    TEST_CASE("pga2dp: closed_loop_system2dp - an inactive constraint is absent (D1)")
    {
        fmt::println("pga2dp: closed_loop_system2dp - inactive constraint == absent");

        // a 3R arm whose hand is pinned to the world; system A carries a SECOND pin at
        // the elbow that is switched off, system B never had it. Every solve must agree
        // to the bit -- an inactive constraint contributes nothing -- and switching the
        // extra pin on must change the answer (so the flag is not a no-op).
        auto const seg = make_plate_body(1.5, 0.5, 0.05);
        auto build = [&](bool extra) {
            closed_loop_system2dp cl;
            cl.add_frame(static_frame2dp("W"));
            cl.add_revolute_body(static_frame2dp("L1", vec2dp{0.25, 0.0, 1.0}, 0.0), seg,
                                 vec2dp{-0.25, 0.0, 1.0}, 0.30, 0.0, cl.index_of("W"));
            cl.add_revolute_body(static_frame2dp("L2", vec2dp{0.5, 0.0, 1.0}, 0.0), seg,
                                 vec2dp{-0.25, 0.0, 1.0}, -0.40, 0.0, cl.index_of("L1"));
            cl.add_revolute_body(static_frame2dp("L3", vec2dp{0.5, 0.0, 1.0}, 0.0), seg,
                                 vec2dp{-0.25, 0.0, 1.0}, 0.50, 0.0, cl.index_of("L2"));
            size_t const L2 = cl.index_of("L2"), L3 = cl.index_of("L3"),
                         W = cl.index_of("W");
            auto unit = [](vec2dp const& p) { return vec2dp{p.x / p.z, p.y / p.z, 1.0}; };
            vec2dp const hand =
                unit(move2dp(vec2dp{0.25, 0.0, 1.0}, cl.system().get_pos_trafo(L3, 0)));
            cl.add_loop_constraint(loop_constraint2dp{L3, vec2dp{0.25, 0.0, 1.0}, W, hand,
                                                      constraint2dp::coincidence});
            if (extra) {
                vec2dp const elbow = unit(
                    move2dp(vec2dp{0.25, 0.0, 1.0}, cl.system().get_pos_trafo(L2, 0)));
                size_t const c = cl.add_loop_constraint(loop_constraint2dp{
                    L2, vec2dp{0.25, 0.0, 1.0}, W, elbow, constraint2dp::coincidence});
                cl.set_loop_active(c, false);
            }
            return cl;
        };

        auto A = build(true);
        auto B = build(false);
        CHECK(A.loop_count() == 2);
        CHECK(A.active_loop_count() == 1);
        CHECK(A.constraint_rows() == 2);

        std::vector<value_t> lA, lB;
        auto const qA = A.joint_accelerations(&lA);
        auto const qB = B.joint_accelerations(&lB);
        CHECK(lA.size() == 2);
        for (size_t k = 0; k < qA.size(); ++k)
            CHECK(qA[k] == qB[k]); // bit-identical, not approximately
        for (size_t k = 0; k < lA.size(); ++k)
            CHECK(lA[k] == lB[k]);

        for (int i = 0; i < 100; ++i) {
            A.step(2.0e-3);
            B.step(2.0e-3);
        }
        for (size_t j = 1; j < 4; ++j) {
            CHECK(A.joint_phi(j) == B.joint_phi(j));
            CHECK(A.joint_rate(j) == B.joint_rate(j));
        }

        // switch the extra pin on: the elbow is now held too, and the answer changes
        A.set_loop_active(1, true);
        CHECK(A.constraint_rows() == 4);
        A.assemble();
        std::vector<value_t> lA2;
        auto const qA2 = A.joint_accelerations(&lA2);
        CHECK(lA2.size() == 4);
        bool differs = false;
        for (size_t k = 0; k < qA2.size(); ++k)
            differs = differs || (qA2[k] != qA[k]);
        CHECK(differs);

        // and with EVERY constraint off the closed-loop system is the open chain: its
        // step must agree with dynamic_system2dp's to integrator round-off
        auto C = build(false);
        C.set_loop_active(0, false);
        CHECK(C.constraint_rows() == 0);
        dynamic_system2dp open;
        open.add_frame(static_frame2dp("W"));
        open.add_revolute_body(static_frame2dp("L1", vec2dp{0.25, 0.0, 1.0}, 0.0), seg,
                               vec2dp{-0.25, 0.0, 1.0}, 0.30, 0.0, open.index_of("W"));
        open.add_revolute_body(static_frame2dp("L2", vec2dp{0.5, 0.0, 1.0}, 0.0), seg,
                               vec2dp{-0.25, 0.0, 1.0}, -0.40, 0.0, open.index_of("L1"));
        open.add_revolute_body(static_frame2dp("L3", vec2dp{0.5, 0.0, 1.0}, 0.0), seg,
                               vec2dp{-0.25, 0.0, 1.0}, 0.50, 0.0, open.index_of("L2"));
        for (int i = 0; i < 100; ++i) {
            C.step(2.0e-3);
            open.step(2.0e-3);
        }
        for (size_t j = 1; j < 4; ++j)
            CHECK(C.joint_phi(j) == doctest::Approx(open.joint_phi(j)).epsilon(1e-12));

        fmt::println("  inactive: bit-identical to the system without the pin; all off: "
                     "the open chain; active: {} rows, answer changes",
                     A.constraint_rows());
        fmt::println("");
    }

    TEST_CASE("pga2dp: closed_loop_system2dp - distance and frame kinds (D1)")
    {
        fmt::println("pga2dp: closed_loop_system2dp - distance (rod) and frame (weld)");

        value_t const g = 9.81;
        auto const seg = make_plate_body(2.0, 0.6, 0.06); // 0.6 long along x
        auto unit = [](vec2dp const& p) { return vec2dp{p.x / p.z, p.y / p.z, 1.0}; };
        auto torque = [](vec2dp const& Q, bivec2dp const& L) {
            return moment_about(Q, L).z;
        };

        // -- distance: a hinged bar whose tip is tied by a rod to a fixed point. One dof,
        //    one row: static. The multiplier is the rod tension, pinned against the GA
        //    moment balance about the hinge, moment_about(hinge, weight line + rod line)
        //    = 0, with the KKT layer's sign convention: the constraint force on anchor a
        //    is -lambda * n, n the unit direction from b to a, so lambda > 0 is tension.
        {
            closed_loop_system2dp cl;
            cl.add_frame(static_frame2dp("W"));
            cl.add_revolute_body(static_frame2dp("bar", vec2dp{0.3, 0.0, 1.0}, 0.0), seg,
                                 vec2dp{-0.3, 0.0, 1.0}, -0.40, 0.0, cl.index_of("W"));
            size_t const bar = cl.index_of("bar"), W = cl.index_of("W");
            vec2dp const tip_b{0.3, 0.0, 1.0};
            vec2dp const tip = unit(move2dp(tip_b, cl.system().get_pos_trafo(bar, 0)));
            vec2dp const anchor{tip.x + 0.2, tip.y + 0.7, 1.0}; // the rod's other end
            vec2dp const d{tip.x - anchor.x, tip.y - anchor.y, 0.0};
            value_t const L = std::sqrt(d.x * d.x + d.y * d.y);
            cl.add_loop_constraint(
                loop_constraint2dp{bar, tip_b, W, anchor, constraint2dp::distance, L});
            CHECK(cl.residual_norm() < 1e-12);
            CHECK(cl.constraint_rows() == 1);

            std::vector<value_t> lambda;
            auto const qdd = cl.joint_accelerations(&lambda);
            CHECK(std::abs(qdd[0]) < 1e-9); // static
            CHECK(lambda.size() == 1);

            vec2dp const hinge =
                unit(move2dp(vec2dp{-0.3, 0.0, 1.0}, cl.system().get_pos_trafo(bar, 0)));
            vec2dp const cm = unit(move2dp(O_2dp, cl.system().get_pos_trafo(bar, 0)));
            vec2dp const n{d.x / L, d.y / L, 0.0};
            value_t const Mw = torque(hinge, wdg(cm, vec2dp{0.0, -2.0 * g, 0.0}));
            value_t const Mr = torque(hinge, wdg(tip, n)); // per unit rod force
            value_t const lambda_pred = Mw / Mr;           // Mw - lambda * Mr = 0
            fmt::println("  rod: lambda = {: .5f}  statics = {: .5f}", lambda[0],
                         lambda_pred);
            CHECK(lambda[0] == doctest::Approx(lambda_pred).epsilon(1e-9));

            // dynamic: free the bar's hinge by adding a second bar, so the rod-tied
            // chain has one dof, and check the energy over a swing
            closed_loop_system2dp dyn;
            dyn.add_frame(static_frame2dp("W"));
            dyn.add_revolute_body(static_frame2dp("b1", vec2dp{0.3, 0.0, 1.0}, 0.0), seg,
                                  vec2dp{-0.3, 0.0, 1.0}, -0.40, 0.0, dyn.index_of("W"));
            dyn.add_revolute_body(static_frame2dp("b2", vec2dp{0.6, 0.0, 1.0}, 0.0), seg,
                                  vec2dp{-0.3, 0.0, 1.0}, 0.80, 0.0, dyn.index_of("b1"));
            size_t const b2 = dyn.index_of("b2");
            vec2dp const t2 = unit(move2dp(tip_b, dyn.system().get_pos_trafo(b2, 0)));
            vec2dp const a2{t2.x + 0.5, t2.y + 0.3, 1.0};
            value_t const L2 = std::hypot(t2.x - a2.x, t2.y - a2.y);
            dyn.add_loop_constraint(loop_constraint2dp{b2, tip_b, dyn.index_of("W"), a2,
                                                       constraint2dp::distance, L2});
            value_t const E0 = dyn.system().total_energy();
            value_t max_dE = 0.0, max_g = 0.0;
            for (int i = 0; i < 1500; ++i) {
                dyn.step(1.0e-3);
                max_dE = std::max(max_dE, std::abs(dyn.system().total_energy() - E0));
                max_g = std::max(max_g, dyn.residual_norm());
            }
            fmt::println("  rod-tied double bar: dE/E = {:.2e}, max |g| = {:.2e}",
                         max_dE / std::abs(E0), max_g);
            CHECK(max_dE / std::abs(E0) < 1e-6);
            CHECK(max_g < 1e-9);
        }

        // -- frame: a 3R arm whose hand is WELDED to the world (3 rows, 3 dof: static).
        //    lambda is minus the wrench the wall exerts on the hand, i.e. the total
        //    gravity wrench of the arm about the hand: force = sum of the weights,
        //    moment balance about every joint pivot, in GA (moment_about).
        {
            closed_loop_system2dp cl;
            cl.add_frame(static_frame2dp("W"));
            value_t const q0[3] = {0.3, -0.5, 0.7};
            size_t parent = cl.index_of("W");
            for (int i = 0; i < 3; ++i) {
                cl.add_revolute_body(
                    static_frame2dp("L" + std::to_string(i), vec2dp{0.6, 0.0, 1.0}, 0.0),
                    seg, vec2dp{-0.3, 0.0, 1.0}, q0[i], 0.0, parent);
                parent = cl.index_of("L" + std::to_string(i));
            }
            size_t const hand_f = parent, W = cl.index_of("W");
            vec2dp const hand_b{0.3, 0.0, 1.0};
            vec2dp const hand =
                unit(move2dp(hand_b, cl.system().get_pos_trafo(hand_f, 0)));
            // a world-fixed frame at the hand's current pose, to weld against
            cl.add_frame(static_frame2dp("wall"), W);
            cl.system().set_pose(
                cl.index_of("wall"),
                pose2dp_from_motor(cl.system().get_pos_trafo(hand_f, 0)));
            cl.add_loop_constraint(loop_constraint2dp{hand_f, hand_b, cl.index_of("wall"),
                                                      hand_b, constraint2dp::frame});
            CHECK(cl.residual_norm() < 1e-10);
            CHECK(cl.constraint_rows() == 3);

            std::vector<value_t> lambda;
            auto const qdd = cl.joint_accelerations(&lambda);
            CHECK(lambda.size() == 3);
            for (auto v : qdd)
                CHECK(std::abs(v) < 1e-8); // welded: static

            // Passive joints carry no torque, so about EVERY pivot j the moment of the
            // distal weight lines and of the wall wrench vanishes. The wall acts on the
            // hand with the force -lambda_F at the hand point and the couple -lambda_M
            // (the KKT convention: the constraint force is -G^T lambda). This is NOT the
            // gravity resultant reduced at the hand -- the base pin carries part of the
            // load, so the weld reaction is fixed by the three joint balances, not by the
            // force sum (measured: lambda_F = (38.3, -37.2) N against a weight of 58.9
            // N).
            fmt::println("  weld: lambda = ({: .4f}, {: .4f} | {: .4f})", lambda[0],
                         lambda[1], lambda[2]);
            vec2dp const lam_F{lambda[0], lambda[1], 0.0};
            for (int j = 0; j < 3; ++j) {
                auto const Mj =
                    cl.system().get_pos_trafo(cl.index_of("L" + std::to_string(j)), 0);
                vec2dp const pivot = unit(move2dp(vec2dp{-0.3, 0.0, 1.0}, Mj));
                bivec2dp Wg{}; // the weight lines distal of joint j
                for (int i = j; i < 3; ++i) {
                    auto const Mi = cl.system().get_pos_trafo(
                        cl.index_of("L" + std::to_string(i)), 0);
                    Wg = Wg + wdg(unit(move2dp(O_2dp, Mi)), vec2dp{0.0, -2.0 * g, 0.0});
                }
                value_t const r =
                    torque(pivot, Wg) - torque(pivot, wdg(hand, lam_F)) - lambda[2];
                fmt::println("    joint {}: moment balance residual = {: .2e}", j, r);
                CHECK(std::abs(r) < 1e-7);
            }

            // dynamic: a 4R arm welded at the hand keeps one dof (a four-bar with a
            // welded coupler); energy conserved and the weld holds
            closed_loop_system2dp dyn;
            dyn.add_frame(static_frame2dp("W"));
            value_t const q1[4] = {0.3, -0.5, 0.7, 0.2};
            parent = dyn.index_of("W");
            for (int i = 0; i < 4; ++i) {
                dyn.add_revolute_body(
                    static_frame2dp("L" + std::to_string(i), vec2dp{0.6, 0.0, 1.0}, 0.0),
                    seg, vec2dp{-0.3, 0.0, 1.0}, q1[i], 0.0, parent);
                parent = dyn.index_of("L" + std::to_string(i));
            }
            size_t const h4 = parent;
            dyn.add_frame(static_frame2dp("wall"), dyn.index_of("W"));
            dyn.system().set_pose(dyn.index_of("wall"),
                                  pose2dp_from_motor(dyn.system().get_pos_trafo(h4, 0)));
            dyn.add_loop_constraint(loop_constraint2dp{h4, hand_b, dyn.index_of("wall"),
                                                       hand_b, constraint2dp::frame});
            value_t const E0 = dyn.system().total_energy();
            value_t max_dE = 0.0, max_g = 0.0;
            for (int i = 0; i < 1500; ++i) {
                dyn.step(1.0e-3);
                max_dE = std::max(max_dE, std::abs(dyn.system().total_energy() - E0));
                max_g = std::max(max_g, dyn.residual_norm());
            }
            fmt::println("  welded 4R: dE/E = {:.2e}, max |g| = {:.2e}",
                         max_dE / std::abs(E0), max_g);
            CHECK(max_dE / std::abs(E0) < 1e-6);
            CHECK(max_g < 1e-9);
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: closed_loop_system2dp - two arms welded hand to hand (D1)")
    {
        fmt::println("pga2dp: closed_loop_system2dp - weld between two MOVING frames");

        // The weld tests above close against a world-fixed frame, which has no joint
        // ancestors -- the b side of the frame kind's rows never contributes there. Here
        // both sides move: arm A (2R) and arm B (1R, then 2R) are welded hand to hand.
        // Exact closure by construction: arm B is built under its own base frame W2,
        // whose pose is then set to M_handA (x) rrev(M_handB_rel), so B's hand lands on
        // A's hand in position AND orientation.
        value_t const g = 9.81;
        auto const seg = make_plate_body(2.0, 0.6, 0.06);
        // arm B's links are SHORTER: a 1-link arm B of A's length, welded with A's hand
        // orientation, would sit exactly on top of A's last link (its base pivot on A's
        // elbow) -- a degenerate structure with a singular KKT matrix. Measured before
        // this was understood: lambda ~ 1e16.
        auto const segB = make_plate_body(1.5, 0.4, 0.06);
        auto unit = [](vec2dp const& p) { return vec2dp{p.x / p.z, p.y / p.z, 1.0}; };
        auto torque = [](vec2dp const& Q, bivec2dp const& L) {
            return moment_about(Q, L).z;
        };
        vec2dp const hand_b{0.3, 0.0, 1.0}, piv_b{-0.3, 0.0, 1.0};
        vec2dp const hand_bB{0.2, 0.0, 1.0}, piv_bB{-0.2, 0.0, 1.0};
        vec2dp const w{0.0, -2.0 * g, 0.0}, wB{0.0, -1.5 * g, 0.0};

        // build: nA links on arm A from W, nB links on arm B from W2; returns the system
        // with the weld registered and closed
        auto build = [&](int nA, std::vector<value_t> const& qA, int nB,
                         std::vector<value_t> const& qB, std::vector<size_t>& A,
                         std::vector<size_t>& B) {
            closed_loop_system2dp cl;
            cl.add_frame(static_frame2dp("W"));
            size_t parent = cl.index_of("W");
            for (int i = 0; i < nA; ++i) {
                cl.add_revolute_body(
                    static_frame2dp("A" + std::to_string(i), vec2dp{0.6, 0.0, 1.0}, 0.0),
                    seg, piv_b, qA[i], 0.0, parent);
                parent = cl.index_of("A" + std::to_string(i));
                A.push_back(parent);
            }
            cl.add_frame(static_frame2dp("W2"), cl.index_of("W"));
            parent = cl.index_of("W2");
            for (int i = 0; i < nB; ++i) {
                cl.add_revolute_body(
                    static_frame2dp("B" + std::to_string(i), vec2dp{0.4, 0.0, 1.0}, 0.0),
                    segB, piv_bB, qB[i], 0.0, parent);
                parent = cl.index_of("B" + std::to_string(i));
                B.push_back(parent);
            }
            // place W2 so that B's hand frame coincides with A's hand frame (the two
            // hand anchors differ in body coordinates but land on the same world point)
            auto const M_hA = cl.system().get_pos_trafo(A.back(), 0);
            auto const M_hB = cl.system().get_pos_trafo(B.back(), cl.index_of("W2"));
            auto const T_ab = rgpr(M_hA, rrev(M_hB)); // B's hand frame onto A's
            // both anchors are their link's tip, (0.3, 0) on A and (0.2, 0) on B, so with
            // the frames coincident the tips differ by 0.1 along the shared x: shift W2
            // by that difference (in world) so the two tips coincide exactly
            pose2dp p2 = pose2dp_from_motor(T_ab);
            vec2dp const tipA = unit(move2dp(hand_b, M_hA));
            vec2dp const tipB = unit(move2dp(hand_bB, rgpr(T_ab, M_hB)));
            p2.origin =
                vec2dp{p2.origin.x + tipA.x - tipB.x, p2.origin.y + tipA.y - tipB.y, 1.0};
            cl.system().set_pose(cl.index_of("W2"), p2);
            cl.add_loop_constraint(loop_constraint2dp{A.back(), hand_b, B.back(), hand_bB,
                                                      constraint2dp::frame});
            return cl;
        };

        // -- static: 2R + 1R = 3 dof, 3 rows. The interaction wrench is pinned by the
        //    moment balance about every pivot of BOTH arms: on A's hand the force is
        //    -lambda_F, the couple -lambda_M; on B's hand +lambda_F, +lambda_M.
        {
            std::vector<size_t> A, B;
            auto cl = build(2, {0.4, -0.7}, 1, {0.9}, A, B);
            CHECK(cl.residual_norm() < 1e-12);
            std::vector<value_t> lambda;
            auto const qdd = cl.joint_accelerations(&lambda);
            for (auto v : qdd)
                CHECK(std::abs(v) < 1e-8);
            vec2dp const hand =
                unit(move2dp(hand_b, cl.system().get_pos_trafo(A.back(), 0)));
            vec2dp const lam_F{lambda[0], lambda[1], 0.0};
            auto balance = [&](std::vector<size_t> const& arm, vec2dp const& piv,
                               vec2dp const& wt, value_t sgn) {
                for (size_t j = 0; j < arm.size(); ++j) {
                    vec2dp const pivot =
                        unit(move2dp(piv, cl.system().get_pos_trafo(arm[j], 0)));
                    bivec2dp Wg{};
                    for (size_t i = j; i < arm.size(); ++i)
                        Wg = Wg + wdg(unit(move2dp(O_2dp,
                                                   cl.system().get_pos_trafo(arm[i], 0))),
                                      wt);
                    value_t const r = torque(pivot, Wg) +
                                      sgn * (torque(pivot, wdg(hand, lam_F)) + lambda[2]);
                    CHECK(std::abs(r) < 1e-7);
                }
            };
            balance(A, piv_b, w, -1.0);
            balance(B, piv_bB, wB, +1.0);
            fmt::println("  2R + 1R welded: static, both arms' joint balances hold, "
                         "lambda = ({: .3f}, {: .3f} | {: .3f})",
                         lambda[0], lambda[1], lambda[2]);
        }

        // -- dynamic: 2R + 2R = 4 dof, 3 rows: one dof. Energy conserved and the weld
        //    holds through a swing (both sides of the rows exercised every sub-step).
        //    Joint A1 passes through -pi at t ~ 0.65 s: its motor crosses to the far
        //    sheet of the double cover, and without the sign normalisation in
        //    relative_rotation the orientation residual reads a full turn there and the
        //    projection spins joints through 2 pi -- measured as dE/E jumping to 0.17
        //    and then 0.53. This case is the regression gate for that.
        {
            std::vector<size_t> A, B;
            auto cl = build(2, {0.4, -0.7}, 2, {0.9, 0.5}, A, B);
            CHECK(cl.residual_norm() < 1e-12);
            value_t const E0 = cl.system().total_energy();
            value_t max_dE = 0.0, max_g = 0.0;
            for (int i = 0; i < 1500; ++i) {
                cl.step(1.0e-3);
                max_dE = std::max(max_dE, std::abs(cl.system().total_energy() - E0));
                max_g = std::max(max_g, cl.residual_norm());
            }
            fmt::println("  2R + 2R welded: dE/E = {:.2e}, max |g| = {:.2e}",
                         max_dE / std::abs(E0), max_g);
            CHECK(max_dE / std::abs(E0) < 1e-6);
            CHECK(max_g < 1e-9);
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: dynamic_system2dp - a massless link may carry a joint (D2)")
    {
        fmt::println("pga2dp: dynamic_system2dp - massless carrier link in a chain");

        // massive -> MASSLESS -> massive: the carrier holds a joint but no inertia. It
        // used to be refused at construction (make_plate_body inverted a singular
        // inertia); now it builds, and the chain integrates with energy conserved. The
        // carrier's own joint moves the massive link below it, so the mass matrix stays
        // regular.
        auto const heavy = make_plate_body(2.0, 0.6, 0.06);
        auto const carrier = make_plate_body(0.0, 0.3, 0.03); // massless
        auto const light = make_plate_body(1.5, 0.5, 0.05);
        CHECK(carrier.mass == 0.0);

        dynamic_system2dp sys;
        sys.add_frame(static_frame2dp("W"));
        sys.add_revolute_body(static_frame2dp("A", vec2dp{0.3, 0.0, 1.0}, 0.0), heavy,
                              vec2dp{-0.3, 0.0, 1.0}, 0.4, 0.0, sys.index_of("W"));
        sys.add_revolute_body(static_frame2dp("C", vec2dp{0.45, 0.0, 1.0}, 0.0), carrier,
                              vec2dp{-0.15, 0.0, 1.0}, -0.7, 0.0, sys.index_of("A"));
        sys.add_revolute_body(static_frame2dp("B", vec2dp{0.4, 0.0, 1.0}, 0.0), light,
                              vec2dp{-0.25, 0.0, 1.0}, 0.5, 0.0, sys.index_of("C"));
        auto const M = sys.mass_matrix();
        CHECK(M.size() == 9);
        for (size_t j = 0; j < 3; ++j)
            CHECK(M[j * 3 + j] > 0.0); // every joint moves some inertia

        value_t const E0 = sys.total_energy();
        value_t max_dE = 0.0;
        for (int i = 0; i < 1500; ++i) {
            sys.step(1.0e-3);
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));
        }
        fmt::println("  heavy -> massless -> light: dE/E = {:.2e}",
                     max_dE / std::abs(E0));
        CHECK(max_dE / std::abs(E0) < 1e-6);

        // a massless LEAF with a dof joint moves no inertia: refused with the cause --
        // unless the joint is driven, which takes it out of the coupled solve
        dynamic_system2dp leaf;
        leaf.add_frame(static_frame2dp("W"));
        leaf.add_revolute_body(static_frame2dp("A", vec2dp{0.3, 0.0, 1.0}, 0.0), heavy,
                               vec2dp{-0.3, 0.0, 1.0}, 0.4, 0.0, leaf.index_of("W"));
        leaf.add_revolute_body(static_frame2dp("C", vec2dp{0.45, 0.0, 1.0}, 0.0), carrier,
                               vec2dp{-0.15, 0.0, 1.0}, -0.7, 0.0, leaf.index_of("A"));
        CHECK_THROWS_WITH_AS(leaf.step(1.0e-3), doctest::Contains("moves no inertia"),
                             std::runtime_error);
        leaf.set_driven_rate(leaf.index_of("C"), 2.0);
        CHECK_NOTHROW(leaf.step(1.0e-3));

        // a massless FREE body is inert: skipped by the free-body integrator, no throw
        dynamic_system2dp free;
        free.add_frame(static_frame2dp("W"));
        free.add_body(static_frame2dp("ghost", vec2dp{1.0, 1.0, 1.0}, 0.0), carrier);
        CHECK_NOTHROW(free.step(1.0e-3));
        fmt::println("  massless leaf joint refused with the cause; driven it runs; a "
                     "massless free body is inert");
        fmt::println("");
    }

    TEST_CASE("pga2dp: make_body_from_inertia and make_point_body2dp")
    {
        fmt::println(
            "pga2dp: body builders from a moment of inertia and for a point mass");

        // the moment builder reproduces the plate builder byte for byte
        value_t const m = 2.0, w = 0.6, h = 0.2;
        auto const pl = make_plate_body(m, w, h);
        auto const gen = make_body_from_inertia(m, m * (w * w + h * h) / 12.0);
        for (size_t i = 0; i < 9; ++i) {
            CHECK(gen.I.data[i] == pl.I.data[i]);
            CHECK(gen.I_inv.data[i] == pl.I_inv.data[i]);
        }

        // point mass on a massless rod = the simple pendulum: alpha0 = -g sin(theta)/L
        // exactly, energy conserved over the swing
        value_t const g = 9.81, Lr = 0.8, th0 = 0.6;
        dynamic_system2dp sys;
        sys.add_frame(static_frame2dp("W"));
        // the pendulum is ONE revolute body -- the point mass -- with its hinge L away in
        // body coordinates: the massless rod is the pivot offset, not a body of its own
        // (a second, unjointed body on the path would be a FREE body, refused until F)
        sys.add_revolute_body(static_frame2dp("bob", vec2dp{0.0, -Lr, 1.0}, 0.0),
                              make_point_body2dp(1.7), vec2dp{0.0, Lr, 1.0}, th0, 0.0,
                              sys.index_of("W"));
        CHECK(sys.joint_accel(sys.index_of("bob")) ==
              doctest::Approx(-g * std::sin(th0) / Lr).epsilon(1e-12));
        value_t const E0 = sys.total_energy();
        value_t max_dE = 0.0;
        for (int i = 0; i < 2000; ++i) {
            sys.step(1.0e-3);
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));
        }
        fmt::println("  point bob on a massless rod: alpha0 = -g sin(th)/L exact, dE/E = "
                     "{:.2e}",
                     max_dE / std::abs(E0));
        CHECK(max_dE / std::abs(E0) < 1e-8);
        fmt::println("");
    }

    TEST_CASE("pga2dp: screw_axis and jacobian (D3)")
    {
        fmt::println("pga2dp: screw_axis (centre, angle) and the Jacobian vs finite "
                     "differences");

        // a rotation about a point: the generator IS the centre scaled by the angle, so
        // screw_axis recovers the centre by unitizing and the full angle from 2 rlog
        vec2dp const C{0.7, -0.4, 1.0};
        value_t const theta = 0.9;
        auto const M = get_motor(C, theta);
        auto const s = screw_axis(M);
        CHECK(s.angle == doctest::Approx(theta).epsilon(1e-12));
        CHECK(s.centre.x == doctest::Approx(C.x).epsilon(1e-12));
        CHECK(s.centre.y == doctest::Approx(C.y).epsilon(1e-12));
        CHECK(is_same_motion(get_motor(s.centre, s.angle), M));
        // a translation: angle 0 and the displacement vector
        auto const T = get_motor(vec2dp{0.3, -0.4, 0.0});
        auto const st = screw_axis(T);
        CHECK(st.angle == 0.0);
        CHECK(st.translation.x == doctest::Approx(0.3).epsilon(1e-12));
        CHECK(st.translation.y == doctest::Approx(-0.4).epsilon(1e-12));

        // the Jacobian of a 3R arm's hand: body columns vs central finite differences
        // of the hand motor's log, and the space form as the body form transported
        auto const seg = make_plate_body(2.0, 0.6, 0.06);
        value_t const q0[3] = {0.3, -0.5, 0.7};
        dynamic_system2dp sys;
        sys.add_frame(static_frame2dp("W"));
        size_t parent = sys.index_of("W");
        for (int i = 0; i < 3; ++i) {
            sys.add_revolute_body(
                static_frame2dp("L" + std::to_string(i), vec2dp{0.6, 0.0, 1.0}, 0.0), seg,
                vec2dp{-0.3, 0.0, 1.0}, q0[i], 0.0, parent);
            parent = sys.index_of("L" + std::to_string(i));
        }
        size_t const hand = parent;
        auto const rj = sys.dof_joints();
        auto const Jb = sys.jacobian_columns(hand, true);
        auto const Js = sys.jacobian_columns(hand, false);
        mvec2dp_u const Mf = sys.get_pos_trafo(hand, 0);
        value_t const h = 1.0e-5;
        value_t max_err = 0.0;
        for (size_t j = 0; j < 3; ++j) {
            value_t const q = sys.joint_phi(rj[j]);
            sys.set_joint(rj[j], q + h);
            mvec2dp_u const Mp = sys.get_pos_trafo(hand, 0);
            sys.set_joint(rj[j], q - h);
            mvec2dp_u const Mm = sys.get_pos_trafo(hand, 0);
            sys.set_joint(rj[j], q);
            twist2dp const fd = (1.0 / h) * rlog(rgpr(rrev(Mm), Mp));
            twist2dp const d = fd - Jb[j];
            max_err = std::max(max_err, std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z));
            CHECK(is_close(move2dp(Jb[j], Mf), Js[j]));
        }
        fmt::println("  max |J_b column - FD of rlog| = {:.2e}", max_err);
        CHECK(max_err < 1e-8);

        // V = J q-dot reproduces twist_world; flat layout and subset
        value_t const qd[3] = {0.5, -0.3, 0.8};
        for (size_t j = 0; j < 3; ++j) {
            sys.set_joint_rate(rj[j], qd[j]);
        }
        twist2dp V{};
        for (size_t j = 0; j < 3; ++j)
            V = V + qd[j] * Js[j];
        CHECK(is_close(V, sys.twist_world(hand)));
        auto const J = sys.jacobian(hand);
        CHECK(J.size() == 9);
        CHECK(J[2 * 3 + 1] == Js[1].z);
        auto const Jsub = sys.jacobian(hand, std::vector<size_t>{rj[2]});
        CHECK(Jsub.size() == 3);
        CHECK(Jsub[1] == Js[2].y);
        // every revolute column is a unit-rate rotation about its world pivot: its screw
        // axis is that pivot with angle 1
        for (size_t j = 0; j < 3; ++j) {
            auto const sj = screw_axis(Js[j]);
            vec2dp const piv =
                unitize(move2dp(vec2dp{-0.3, 0.0, 1.0}, sys.get_pos_trafo(rj[j], 0)));
            CHECK(sj.angle == doctest::Approx(1.0));
            CHECK(sj.centre.x == doctest::Approx(piv.x).epsilon(1e-12));
            CHECK(sj.centre.y == doctest::Approx(piv.y).epsilon(1e-12));
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: screw_system - span and rcmt-closure in the plane (W1)")
    {
        fmt::println("pga2dp: screw_system - span vs. Lie closure of planar screws");

        // planar screws: a rotation about a point is the point scaled by the rate; two
        // distinct centres bracket to a translation and close se(2) (3); a single
        // rotation is its own subalgebra; two translations commute (closure 2)
        vec2dp const R1{0.0, 0.0, 1.0}, R2{1.0, 0.0, 1.0};
        vec2dp const T1{1.0, 0.0, 0.0}, T2{0.0, 1.0, 0.0};
        auto const one = screw_system({R1});
        CHECK(one.span == 1);
        CHECK(one.closure == 1);
        auto const two = screw_system({R1, R2});
        CHECK(two.span == 2);
        CHECK(two.closure == 3);
        auto const trans = screw_system({T1, T2});
        CHECK(trans.span == 2);
        CHECK(trans.closure == 2); // translations commute: a subalgebra
        auto const dep = screw_system({R1, 3.0 * R1});
        CHECK(dep.span == 1);
        fmt::println("  one rotation: ({},{})  two centres: ({},{})  two translations: "
                     "({},{})",
                     one.span, one.closure, two.span, two.closure, trans.span,
                     trans.closure);
        fmt::println("");
    }

    TEST_CASE("pga2dp: joint torque - an actuator as a generalised force")
    {
        fmt::println("pga2dp: dynamic_system2dp - joint torque (actuator)");

        // a gravity pendulum: a plate of mass m whose centre hangs l below the pivot
        value_t const m = 2.0, l = 0.5, g = 9.81, w = 0.05;
        value_t const I_piv = m * (w * w + 4.0 * l * l) / 12.0 + m * l * l;
        auto make = [&](value_t gg) {
            dynamic_system2dp s;
            s.set_gravity(vec2dp{0.0, -gg, 0.0});
            s.add_frame(static_frame2dp("W"));
            s.add_revolute_body(static_frame2dp("P", vec2dp{0.0, -l, 1.0}, 0.0),
                                make_plate_body(m, w, 2.0 * l), vec2dp{0.0, l, 1.0});
            return s;
        };

        // a torque of zero changes nothing, bit for bit
        {
            auto a = make(g), b = make(g);
            b.set_joint_torque(1, [](value_t) { return 0.0; });
            a.set_joint(1, 0.3);
            b.set_joint(1, 0.3);
            for (int i = 0; i < 200; ++i) {
                a.step(1e-3);
                b.step(1e-3);
            }
            CHECK(a.joint_phi(1) == b.joint_phi(1));
            CHECK(a.joint_omega(1) == b.joint_omega(1));
        }
        // sign and magnitude: at rest a torque tau gives q-ddot = tau / I_pivot
        {
            auto s = make(0.0);
            s.set_joint_torque(1, [](value_t) { return 0.7; });
            CHECK(s.joint_accel(1) == doctest::Approx(0.7 / I_piv).epsilon(1e-12));
        }
        // a constant torque against gravity settles (damped) at m g l sin(q) = tau
        {
            auto s = make(g);
            value_t const tau = 0.3 * m * g * l;
            s.set_joint_torque(1, [tau](value_t) { return tau; });
            s.set_joint_spring_damper(1, 0.0, 2.0);
            for (int i = 0; i < 20000; ++i)
                s.step(1e-3);
            CHECK(s.joint_phi(1) == doctest::Approx(std::asin(0.3)).epsilon(1e-7));
        }
        // an actuator is not conservative: the change of total energy is its work,
        // integral(tau q-dot dt). The trapezoid work integral is second order, so the
        // mismatch must fall 4x per halving of dt (a gate by convergence order)
        {
            auto work_err = [&](value_t dt) {
                auto s = make(g);
                s.set_joint_torque(1, [](value_t t) { return 1.5 * std::sin(3.0 * t); });
                value_t const E0 = s.total_energy();
                value_t W = 0.0, p_prev = 0.0;
                for (int i = 0; i * dt < 2.0; ++i) {
                    value_t const p = 1.5 * std::sin(3.0 * s.time()) * s.joint_omega(1);
                    s.step(dt);
                    value_t const p_new =
                        1.5 * std::sin(3.0 * s.time()) * s.joint_omega(1);
                    W += 0.5 * (p + p_new) * dt;
                    p_prev = p_new;
                }
                (void)p_prev;
                return std::abs(s.total_energy() - E0 - W);
            };
            value_t const e1 = work_err(2e-3), e2 = work_err(1e-3);
            CHECK(e1 > 1e-9); // the mismatch is a real discretisation error ...
            CHECK(e1 / e2 == doctest::Approx(4.0).epsilon(0.1)); // ... of second order
        }
        // a driven joint ignores a torque: its coordinate is prescribed
        {
            auto s = make(g);
            s.set_driven_rate(1, 0.5);
            s.set_joint_torque(1, [](value_t) { return 100.0; });
            for (int i = 0; i < 1000; ++i)
                s.step(1e-3);
            CHECK(s.joint_phi(1) == doctest::Approx(0.5).epsilon(1e-12));
        }
        // the motor-joint form: a free body (3 screws) with one generalised force per
        // screw; the solve returns M q-ddot = tau (the fold lands on the right index)
        {
            dynamic_system2dp s;
            s.set_gravity(vec2dp{0.0, 0.0, 0.0});
            s.add_frame(static_frame2dp("W"));
            s.add_body(static_frame2dp("B", vec2dp{0.3, 0.2, 1.0}, 0.4),
                       make_plate_body(3.0, 0.2, 0.5));
            std::vector<value_t> const tau{0.0, 2.0, -0.5};
            s.set_joint_torque(1, [tau](value_t) { return tau; });
            auto const M = s.mass_matrix();
            auto const qdd = s.joint_accelerations();
            REQUIRE(qdd.size() == 3);
            for (size_t r = 0; r < 3; ++r) {
                value_t acc = 0.0;
                for (size_t c = 0; c < 3; ++c)
                    acc += M[3 * r + c] * qdd[c];
                CHECK(acc == doctest::Approx(tau[r]).epsilon(1e-9).scale(1.0));
            }
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: mass_bias / constraint_jacobian - the EoM read off (L4)")
    {
        fmt::println("pga2dp: mass_bias() and constraint_jacobian() accessors");

        // the public read-off of  M(q) q-ddot = RHS + tau - G^T lambda  for an external
        // control law (computed torque). Falsified during development by forcing
        // with_joint_torques = true in mass_bias() (the exclusion gate below fails) and
        // by flipping the sign of one FD column (the Jacobian gate fails).

        value_t const m = 2.0, l = 0.5, g = 9.81, w = 0.05;

        // -- open chain: a pendulum carrying every PASSIVE force element (gravity,
        //    joint spring/damper, applied wrench). mass_bias() must close the equation
        //    of motion: lu_solve(M, RHS) = joint_accelerations().
        {
            dynamic_system2dp s;
            s.set_gravity(vec2dp{0.0, -g, 0.0});
            s.add_frame(static_frame2dp("W"));
            s.add_revolute_body(static_frame2dp("P", vec2dp{0.0, -l, 1.0}, 0.0),
                                make_plate_body(m, w, 2.0 * l), vec2dp{0.0, l, 1.0}, 0.4,
                                0.7);
            s.set_joint_spring_damper(1, 3.0, 0.2, 0.1);
            s.set_applied_wrench(1, [](value_t) { return bivec2dp{0.4, 0.0, -0.6}; });

            auto const [M, RHS] = s.mass_bias();
            REQUIRE(M.size() == 1);
            REQUIRE(RHS.size() == 1);

            // M is the mass matrix (no driven joints here: the two paths agree)
            auto const Mm = s.mass_matrix();
            CHECK(M[0] == doctest::Approx(Mm[0]).epsilon(1e-13));

            // the equation of motion closes on the public surface
            auto const qdd = lu_solve(M, RHS, 1);
            auto const qdd_sys = s.joint_accelerations();
            CHECK(qdd[0] == doctest::Approx(qdd_sys[0]).epsilon(1e-12));

            // the EXCLUSION: registering an actuator torque leaves RHS bit-identical
            // (a feedback law reading the accessor must not see its own output), while
            // the system's accelerations shift by exactly M^-1 tau
            value_t const tau = 0.8;
            s.set_joint_torque(1, [tau](value_t) { return tau; });
            auto const [M2, RHS2] = s.mass_bias();
            CHECK(M2[0] == M[0]);
            CHECK(RHS2[0] == RHS[0]); // bit-identical: tau is NOT in the bias
            std::vector<value_t> rhs_tau{RHS2[0] + tau};
            auto const qdd_tau = lu_solve(M2, rhs_tau, 1);
            CHECK(qdd_tau[0] ==
                  doctest::Approx(s.joint_accelerations()[0]).epsilon(1e-12));
            fmt::println("  pendulum: qdd = {:.6f}, with tau = {:.6f}", qdd_sys[0],
                         qdd_tau[0]);
        }

        // -- a kinematically DRIVEN child carrying inertia: mass_bias()'s M is the
        //    matrix forward dynamics solves (moving-base inertia retained), which
        //    mass_matrix() does not include -- the documented difference
        {
            dynamic_system2dp s;
            s.set_gravity(vec2dp{0.0, -g, 0.0});
            s.add_frame(static_frame2dp("W"));
            s.add_revolute_body(static_frame2dp("P", vec2dp{0.0, -l, 1.0}, 0.0),
                                make_plate_body(m, w, 2.0 * l), vec2dp{0.0, l, 1.0}, 0.3,
                                0.0);
            s.add_revolute_body(static_frame2dp("D", vec2dp{0.0, -l, 1.0}, 0.0),
                                make_plate_body(1.5, 0.3, 0.3), vec2dp{0.0, 0.0, 1.0},
                                0.0, 0.0, 1);
            s.set_driven_rate(2, 5.0);
            auto const [M, RHS] = s.mass_bias();
            REQUIRE(M.size() == 1); // one dof: the driven joint is not a coordinate
            CHECK(M[0] > s.mass_matrix()[0]); // the driven body's inertia is in M
            auto const qdd = lu_solve(M, RHS, 1);
            CHECK(qdd[0] == doctest::Approx(s.joint_accelerations()[0]).epsilon(1e-12));
        }

        // -- closed loop: constraint_jacobian() against central finite differences of
        //    the residual, coordinate by coordinate, with all three constraint kinds
        //    active at once (coincidence 2 rows, distance 1, frame 3)
        {
            auto const seg = make_plate_body(2.0, 0.6, 0.06);
            auto unit = [](vec2dp const& p) { return vec2dp{p.x / p.z, p.y / p.z, 1.0}; };
            closed_loop_system2dp cl;
            cl.add_frame(static_frame2dp("W"));
            cl.add_revolute_body(static_frame2dp("b1", vec2dp{0.3, 0.0, 1.0}, 0.0), seg,
                                 vec2dp{-0.3, 0.0, 1.0}, 0.35, 0.0, cl.index_of("W"));
            cl.add_revolute_body(static_frame2dp("b2", vec2dp{0.6, 0.0, 1.0}, 0.0), seg,
                                 vec2dp{-0.3, 0.0, 1.0}, -0.55, 0.0, cl.index_of("b1"));
            cl.add_revolute_body(static_frame2dp("b3", vec2dp{0.6, 0.0, 1.0}, 0.0), seg,
                                 vec2dp{-0.3, 0.0, 1.0}, 0.25, 0.0, cl.index_of("b2"));
            size_t const W = cl.index_of("W");
            size_t const b2 = cl.index_of("b2"), b3 = cl.index_of("b3");
            vec2dp const tip_b{0.3, 0.0, 1.0};
            // anchors at the CURRENT geometry so the residual starts near zero (not
            // required for the derivative test, but keeps the rod axis well-defined)
            vec2dp const t3 = unit(move2dp(tip_b, cl.system().get_pos_trafo(b3, 0)));
            vec2dp const t2 = unit(move2dp(tip_b, cl.system().get_pos_trafo(b2, 0)));
            vec2dp const a2{t2.x + 0.4, t2.y + 0.5, 1.0};
            value_t const L2 = std::hypot(t2.x - a2.x, t2.y - a2.y);
            cl.add_loop_constraint(
                loop_constraint2dp{b3, tip_b, W, t3, constraint2dp::coincidence});
            cl.add_loop_constraint(
                loop_constraint2dp{b2, tip_b, W, a2, constraint2dp::distance, L2});
            cl.add_loop_constraint(
                loop_constraint2dp{b3, tip_b, W, t3, constraint2dp::frame});

            size_t const mrows = cl.constraint_rows();
            auto const rj = cl.system().dof_joints();
            size_t const n = rj.size();
            REQUIRE(mrows == 6); // 2 + 1 + 3
            REQUIRE(n == 3);
            auto const G = cl.constraint_jacobian();
            REQUIRE(G.size() == mrows * n);

            value_t const h = 1e-6;
            value_t max_err = 0.0;
            for (size_t c = 0; c < n; ++c) {
                value_t const q0 = cl.system().joint_phi(rj[c]);
                cl.system().set_joint(rj[c], q0 + h);
                auto const gp = cl.residual();
                cl.system().set_joint(rj[c], q0 - h);
                auto const gm = cl.residual();
                cl.system().set_joint(rj[c], q0); // restore
                for (size_t r = 0; r < mrows; ++r) {
                    value_t const fd = (gp[r] - gm[r]) / (2.0 * h);
                    max_err = std::max(max_err, std::abs(G[r * n + c] - fd));
                }
            }
            fmt::println("  constraint_jacobian vs FD of the residual: max err = {:.2e}",
                         max_err);
            CHECK(max_err < 1e-8);
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: dynamic_system2dp - total mass, centre of mass, gravity wrench")
    {
        fmt::println("pga2dp: mass distribution -- total_mass / centre_of_mass / "
                     "gravity_wrench");
        // a hub (free body, 4 kg at (1, 2)) carrying a revolute link (2 kg, its origin
        // 0.5 further along x at q = 0): the centre of mass by hand, gravity as one
        // wrench whose moment about the centre of mass vanishes and whose moment about
        // the origin is the classic r x f
        dynamic_system2dp sys;
        sys.add_frame(static_frame2dp("W"));
        sys.add_body(static_frame2dp("hub", vec2dp{1.0, 2.0, 1.0}, 0.0),
                     make_plate_body(4.0, 0.2, 0.2));
        sys.add_revolute_body(static_frame2dp("link", vec2dp{0.5, 0.0, 1.0}, 0.0),
                              make_plate_body(2.0, 0.5, 0.05), vec2dp{-0.25, 0.0, 1.0},
                              0.0, 0.0, sys.index_of("hub"));
        CHECK(sys.total_mass() == doctest::Approx(6.0).epsilon(1e-12));
        vec2dp const c = sys.centre_of_mass();
        CHECK(c.x == doctest::Approx((4.0 * 1.0 + 2.0 * 1.5) / 6.0).epsilon(1e-12));
        CHECK(c.y == doctest::Approx(2.0).epsilon(1e-12));
        bivec2dp const W = sys.gravity_wrench();
        CHECK(std::abs(moment_about(c, W).z) < 1e-12); // no moment about the CoM
        value_t const g = -sys.gravity().y;
        CHECK(moment_about(O_2dp, W).z == doctest::Approx(-c.x * 6.0 * g).epsilon(1e-12));
        // the identity behind it: the summed force lines have the same moment about ANY
        // point as ONE force line of the total weight through the centre of mass
        bivec2dp const W1 = wdg(c, 6.0 * sys.gravity());
        vec2dp const R{-2.0, 0.5, 1.0};
        CHECK(moment_about(R, W).z ==
              doctest::Approx(moment_about(R, W1).z).epsilon(1e-12));
        fmt::println("  m = {:.1f} kg, CoM = ({:.4f}, {:.4f}), M_O = {:.3f} N m", 6.0,
                     c.x, c.y, moment_about(O_2dp, W).z);
        fmt::println("");
    }

    TEST_CASE("pga2dp: closed_loop_system2dp - the impact map (W2)")
    {
        fmt::println("pga2dp: closed_loop_system2dp - the impact map (W2)");

        // A rod (free body, cm at the origin, length 1) falls at v = 1 with no rotation
        // and no gravity; its left end reaches a fixed pin. The plastic impact pins the
        // end: angular momentum about the pin is conserved through the impact, so
        //
        //     m v L/2 = (I_cm + m L^2/4) w   ->   w = 3 v / (2 L) for a thin rod,
        //
        // T+ / T- = 3/4, and the impulse on the rod is its momentum change m (w L/2 - v)
        // = -m v/4. The plate has a thickness t, so I_cm = m (L^2 + t^2) / 12 is used
        // as it is (w = 1.49985, not 1.5 -- the gate is at 1e-12, not at "thin").
        value_t const m = 2.0, L = 1.0, t = 0.02, v = 1.0;
        auto const rod = make_plate_body(m, L, t);
        value_t const I_pin = m * (L * L + t * t) / 12.0 + m * L * L / 4.0;
        value_t const w_ref = -m * v * (L / 2.0) / I_pin; // clockwise: the pin is at -x
        auto build = [&](value_t vy) {
            closed_loop_system2dp cl;
            cl.system().set_gravity(vec2dp{0.0, 0.0, 0.0});
            cl.add_frame(static_frame2dp("W"));
            cl.add_body(static_frame2dp("rod", vec2dp{0.0, 0.0, 1.0}, 0.0), rod,
                        kin_state2dp{.vel = vec2dp{0.0, vy, 0.0}, .omega = 0.0});
            size_t const c = cl.add_loop_constraint(loop_constraint2dp{
                cl.index_of("rod"), vec2dp{-0.5, 0.0, 1.0}, cl.index_of("W"),
                vec2dp{-0.5, 0.0, 1.0}, constraint2dp::coincidence});
            cl.set_loop_active(c, false); // in flight: no rows
            return std::pair{std::move(cl), c};
        };
        auto momentum = [](closed_loop_system2dp& cl) {
            bivec2dp P{};
            for (size_t i = 1; i < cl.system().size(); ++i)
                P = P + cl.system().momentum_world(i);
            return P;
        };
        auto end_velocity = [](closed_loop_system2dp& cl) {
            size_t const r = cl.index_of("rod");
            vec2dp const P =
                move2dp(vec2dp{-0.5, 0.0, 1.0}, cl.system().get_pos_trafo(r, 0));
            return cl.system().point_velocity(vec2dp{P.x / P.z, P.y / P.z, 1.0}, r);
        };

        { // plastic touchdown: the analytic post-impact state
            auto [cl, c] = build(-v);
            value_t const T0 = cl.system().total_energy();
            bivec2dp const P0 = momentum(cl);
            std::vector<value_t> const Lam = cl.activate_loop_with_impact(c, 0.0);
            REQUIRE(Lam.size() == 2);
            vec2dp const ve = end_velocity(cl);
            CHECK(std::abs(ve.x) < 1e-12); // G q-dot+ = 0: the end is at rest
            CHECK(std::abs(ve.y) < 1e-12);
            // the cm now moves as a point of the rod rotating about the pin: v_cm = w L/2
            vec2dp const vc =
                cl.system().point_velocity(vec2dp{0.0, 0.0, 1.0}, cl.index_of("rod"));
            value_t const w = vc.y / (L / 2.0);
            CHECK(w == doctest::Approx(w_ref).epsilon(1e-12));
            value_t const T1 = cl.system().total_energy();
            CHECK(T1 == doctest::Approx(0.5 * I_pin * w * w).epsilon(1e-12));
            CHECK(T1 / T0 == doctest::Approx(0.75).epsilon(1e-3)); // thin-rod value
            bivec2dp const dP = momentum(cl) - P0; // impulse = momentum change = -Λ
            CHECK(-Lam[0] == doctest::Approx(dP.x).epsilon(1e-12));
            CHECK(-Lam[1] == doctest::Approx(dP.y).epsilon(1e-12));
            // Λ_y < 0 for a foot pressed on the ground, like its λ_y (the impulse on
            // the rod is upward, -Λ): -m (w L/2 + v) with v downward
            CHECK(Lam[1] == doctest::Approx(-m * (w * L / 2.0 + v)).epsilon(1e-12));
            // afterwards: a pinned rod, energy conserved between impacts
            for (int i = 0; i < 500; ++i)
                cl.step(1.0e-3);
            CHECK(std::abs(cl.system().total_energy() - T1) / T1 < 1e-9);
            fmt::println("  plastic: w = {:.4f} (3v/2L), T+/T- = {:.4f}, Λ_y = {:.4f}", w,
                         T1 / T0, Lam[1]);
        }
        { // restitution e = 1/2: the end rebounds at half its closing speed
            auto [cl, c] = build(-v);
            cl.activate_loop_with_impact(c, 0.5);
            vec2dp const ve = end_velocity(cl);
            CHECK(ve.y == doctest::Approx(0.5 * v).epsilon(1e-12));
            CHECK(std::abs(ve.x) < 1e-12);
            CHECK(cl.system().total_energy() < 0.5 * m * v * v); // still dissipative
        }
        { // zero closing velocity: no impulse, the state untouched to the bit
            auto [cl, c] = build(0.0);
            auto const before = cl.system().joint_rates(cl.index_of("rod"));
            std::vector<value_t> const Lam = cl.activate_loop_with_impact(c, 0.0);
            for (value_t const l : Lam)
                CHECK(l == 0.0);
            auto const after = cl.system().joint_rates(cl.index_of("rod"));
            for (size_t k = 0; k < before.size(); ++k)
                CHECK(after[k] == before[k]);
        }
        { // mass weighting: a heavy hub carrying a light link whose tip hits a pin. The
          // link swings, the hub barely notices -- the Euclidean projection of the drift
          // stabiliser would have corrected the hub's coordinates as much as the link's
            auto const hub = make_plate_body(40.0, 0.4, 0.4);
            auto const link = make_plate_body(0.1, 0.5, 0.02);
            closed_loop_system2dp cl;
            cl.system().set_gravity(vec2dp{0.0, 0.0, 0.0});
            cl.add_frame(static_frame2dp("W"));
            cl.add_body(static_frame2dp("hub", vec2dp{0.0, 0.0, 1.0}, 0.0), hub,
                        kin_state2dp{.vel = vec2dp{0.0, -v, 0.0}, .omega = 0.0});
            cl.add_revolute_body(static_frame2dp("link", vec2dp{0.45, 0.0, 1.0}, 0.0),
                                 link, vec2dp{-0.25, 0.0, 1.0}, 0.0, 0.0,
                                 cl.index_of("hub"));
            size_t const c = cl.add_loop_constraint(loop_constraint2dp{
                cl.index_of("link"), vec2dp{0.25, 0.0, 1.0}, cl.index_of("W"),
                vec2dp{0.7, 0.0, 1.0}, constraint2dp::coincidence});
            cl.set_loop_active(c, false);
            value_t const T0 = cl.system().total_energy();
            cl.activate_loop_with_impact(c, 0.0);
            vec2dp const vh =
                cl.system().point_velocity(vec2dp{0.0, 0.0, 1.0}, cl.index_of("hub"));
            value_t const dv_hub = std::sqrt(vh.x * vh.x + (vh.y + v) * (vh.y + v));
            CHECK(dv_hub < 0.02 * v);
            CHECK(std::abs(cl.joint_rate(cl.index_of("link"))) > 1.0); // the link swings
            CHECK(cl.system().total_energy() < T0);
            fmt::println(
                "  mass weighting: hub |dv| = {:.2e} (v = 1), link rate = {:.3f}", dv_hub,
                cl.joint_rate(cl.index_of("link")));
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: ground_contact2dp - touchdown, lift-off, slope, flat foot, guard")
    {
        fmt::println("pga2dp: ground_contact2dp - unilateral contact by events");

        // A rod (free body, cm at height h, tilted by th0 so its left end is the lower
        // one, at rest) falls onto the floor y = 0 with both ends as contact points.
        // The left end's touchdown must be located at its free-fall time
        // t* = sqrt(2 h_end / g) to the bisection tolerance, then the plastic impact of
        // W2 pins the end: w = m v (L/2 cos th0) / I_pin with v = g t*; the rod swings
        // down and the far end lands too -- the rod lies on the ground.
        value_t const m = 2.0, L = 1.0, t = 0.02, g = 9.81, h = 0.5;
        auto const rod = make_plate_body(m, L, t);
        value_t const I_pin = m * (L * L + t * t) / 12.0 + m * L * L / 4.0;
        auto const floor =
            ground_contact2dp::ground_line(vec2dp{0.0, 0.0, 1.0}, vec2dp{1.0, 0.0, 1.0});
        value_t const th0 = 0.1; // ccw: the left end is lower by 0.5 sin(th0)
        auto build = [&](value_t y0, value_t th = 0.0) {
            closed_loop_system2dp cl;
            cl.system().set_gravity(vec2dp{0.0, -g, 0.0});
            cl.add_frame(static_frame2dp("W"));
            cl.add_body(static_frame2dp("rod", vec2dp{0.0, y0, 1.0}, th), rod);
            return cl;
        };
        {
            closed_loop_system2dp cl = build(h, th0);
            ground_contact2dp gc(cl, floor);
            CHECK(gc.height(vec2dp{3.0, 0.7, 1.0}) ==
                  doctest::Approx(0.7).epsilon(1e-12));
            CHECK(gc.normal().y == doctest::Approx(1.0).epsilon(1e-12));
            size_t const c = gc.add({cl.index_of("rod"), vec2dp{-0.5, 0.0, 1.0}});
            size_t const c2 = gc.add({cl.index_of("rod"), vec2dp{0.5, 0.0, 1.0}});
            CHECK(gc.count() == 2);
            CHECK(!gc.active(c));
            value_t const h_end = gc.contact_height(c);
            CHECK(h_end == doctest::Approx(h - 0.5 * std::sin(th0)).epsilon(1e-12));
            value_t const dt = 1.0e-3;
            for (int i = 0; i < 400 && gc.events().empty(); ++i)
                gc.step(dt);
            REQUIRE(gc.events().size() == 1);
            auto const& ev = gc.events()[0];
            value_t const t_star = std::sqrt(2.0 * h_end / g);
            CHECK(ev.touchdown);
            CHECK(ev.t == doctest::Approx(t_star).epsilon(1e-9)); // exact in time
            CHECK(std::abs(ev.point.y) < 1e-12);
            CHECK(ev.point.x == doctest::Approx(-0.5 * std::cos(th0)).epsilon(1e-9));
            CHECK(gc.active(c));
            CHECK(std::abs(gc.contact_height(c)) < 1e-9);
            // the impact of W2 at the event: the impulse -Lambda is the rod's momentum
            // change, m (v+ - v-) with v- = -g t*, so the post-impact cm velocity and
            // hence the rotation rate about the pinned end follow from the recorded
            // impulse alone (the state has moved on by the rest of the step since)
            value_t const v = g * t_star;
            value_t const arm = 0.5 * L * std::cos(th0); // the cm's lever about the pin
            value_t const w_ref = -m * v * arm / I_pin;
            REQUIRE(ev.impulse.size() == 2);
            value_t const v_plus = -v - ev.impulse[1] / m; // -Lambda_y / m
            value_t const w_ev = v_plus / arm;
            CHECK(w_ev == doctest::Approx(w_ref).epsilon(1e-9));
            value_t const w = cl.system().joint_rates(cl.index_of("rod"))[2];
            CHECK(w == doctest::Approx(w_ref).epsilon(5e-3)); // + < 1 ms of swing
            fmt::println("  touchdown at t = {:.9f} (free fall {:.9f}), w from the "
                         "impulse = {:.6f} (rod-onto-pin {:.6f})",
                         ev.t, t_star, w_ev, w_ref);
            // the pinned rod swings down: its far end (the higher one) reaches the
            // floor and LANDS too -- the rod then lies on the ground on both ends, at
            // rest, each end carrying half the weight
            for (int i = 0; i < 900; ++i)
                gc.step(dt);
            CHECK(gc.active(c));
            CHECK(gc.active(c2));
            CHECK(std::abs(gc.contact_height(c2)) < 1e-9);
            CHECK(gc.normal_force(c) == doctest::Approx(0.5 * m * g).epsilon(1e-6));
            CHECK(gc.normal_force(c2) == doctest::Approx(0.5 * m * g).epsilon(1e-6));
            fmt::println("  the far end landed too: the rod lies flat, {:.3f} + {:.3f} N",
                         gc.normal_force(c), gc.normal_force(c2));
        }

        // lift-off: the pinned rod is lifted by an upward force larger than its weight
        // -- the reaction turns to a pull and the contact is released; without the
        // force it comes down again and lands (a second touchdown)
        {
            closed_loop_system2dp cl = build(0.0);
            ground_contact2dp gc(cl, floor);
            size_t const c = gc.add({cl.index_of("rod"), vec2dp{-0.5, 0.0, 1.0}});
            gc.step(1.0e-3);     // the rod placed ON the floor starts to fall: that IS a
            CHECK(gc.active(c)); // touchdown (height 0, moving down), at t ~ 0
            // drop it from 1 mm so it lands
            closed_loop_system2dp cl2 = build(1.0e-3);
            ground_contact2dp gc2(cl2, floor);
            size_t const c2 = gc2.add({cl2.index_of("rod"), vec2dp{-0.5, 0.0, 1.0}});
            for (int i = 0; i < 50; ++i)
                gc2.step(1.0e-3);
            REQUIRE(gc2.active(c2));
            value_t const t_up = cl2.system().time();
            cl2.system().set_applied_wrench(cl2.index_of("rod"), [&](value_t tt) {
                bool const on = tt >= t_up && tt < t_up + 0.25;
                vec2dp const P = unitize(
                    move2dp(vec2dp{-0.5, 0.0, 1.0}, cl2.system().get_pos_trafo(1, 0)));
                return on ? wdg(P, vec2dp{0.0, 3.0 * m * g, 0.0}) : bivec2dp{};
            });
            size_t n_before = gc2.events().size();
            for (int i = 0; i < 1500; ++i)
                gc2.step(1.0e-3);
            auto const& evs = gc2.events();
            REQUIRE(evs.size() >= n_before + 2);
            CHECK(!evs[n_before].touchdown);      // released ...
            CHECK(evs[n_before].t < t_up + 0.01); // ... as soon as the pull appears
            CHECK(evs[n_before + 1].touchdown);   // ... and landed again
            CHECK(gc2.active(c2));
            CHECK(gc2.normal_force(c2) > 0.0);
            fmt::println("  lifted: released at t = {:.4f}, landed again at t = {:.4f}",
                         evs[n_before].t - t_up, evs[n_before + 1].t - t_up);
        }

        // a slope: the ground line through (0, 0) and (1, -0.2); the rod lands ON the
        // line (height 0 to round-off) and the reaction is along the slope's normal
        {
            closed_loop_system2dp cl = build(h);
            auto const slope = ground_contact2dp::ground_line(vec2dp{0.0, 0.0, 1.0},
                                                              vec2dp{1.0, -0.2, 1.0});
            ground_contact2dp gc(cl, slope);
            CHECK(gc.normal().x == doctest::Approx(0.2 / std::sqrt(1.04)).epsilon(1e-12));
            size_t const c = gc.add({cl.index_of("rod"), vec2dp{-0.5, 0.0, 1.0}});
            for (int i = 0; i < 500 && gc.events().empty(); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            CHECK(std::abs(gc.contact_height(c)) < 1e-9);
            vec2dp const P = gc.contact_point(c);
            CHECK(P.y == doctest::Approx(-0.2 * P.x).epsilon(1e-9)); // on the line
            fmt::println("  slope: landed at ({:.4f}, {:.4f}), height {:.1e}", P.x, P.y,
                         gc.contact_height(c));
        }

        // a flat foot: the rod's end welded to the floor on landing (the impact absorbs
        // the rotation too) -- a horizontal cantilever whose ankle moment is m g L/2,
        // so the centre of pressure sits at L/2 = 0.5 > half_length: tipping
        {
            closed_loop_system2dp cl = build(h, th0);
            ground_contact2dp gc(cl, floor);
            size_t const c = gc.add(
                {cl.index_of("rod"), vec2dp{-0.5, 0.0, 1.0}, contact_kind2dp::flat, 0.1});
            for (int i = 0; i < 500 && gc.events().empty(); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            CHECK(cl.constraint_rows() == 3);
            auto const rates = cl.system().joint_rates(cl.index_of("rod"));
            for (auto r : rates)
                CHECK(std::abs(r) < 1e-9); // welded: at rest entirely
            CHECK(gc.normal_force(c) == doctest::Approx(m * g).epsilon(1e-9));
            CHECK(std::abs(gc.moment(c)) ==
                  doctest::Approx(m * g * 0.5 * std::cos(th0)).epsilon(1e-9));
            CHECK(std::abs(gc.cop(c)) ==
                  doctest::Approx(0.5 * std::cos(th0)).epsilon(1e-9));
            CHECK(gc.tipping(c));
            fmt::println("  flat foot: N = {:.3f} N, M = {:.3f} N m, cop = {:.3f} -> "
                         "tipping; as a pin it swings",
                         gc.normal_force(c), gc.moment(c), gc.cop(c));
            // back to a point contact: the angular row goes, the rod swings
            gc.set_kind(c, contact_kind2dp::point);
            CHECK(cl.constraint_rows() == 2);
            for (int i = 0; i < 100; ++i)
                gc.step(1.0e-3);
            CHECK(std::abs(cl.system().joint_rates(cl.index_of("rod"))[2]) > 0.1);
        }
        fmt::println("");
    }

    TEST_CASE("pga2dp: free-floating chain conserves momentum and energy (F)")
    {
        fmt::println("pga2dp: floating base - a planar space station: momentum + energy");

        // The 2D twin of the decisive Phase F gate: a heavy hub (free root joint)
        // carrying two two-link arms swinging with opposite rates, no gravity, no
        // external wrench -- total linear and angular momentum conserved in the world
        // frame, energy alongside, over a long integration. Same description as in 3D,
        // only the types differ.
        auto const hub = make_plate_body(40.0, 1.2, 0.8);
        auto const link = make_plate_body(3.0, 0.8, 0.08);

        auto build = [&] {
            dynamic_system2dp sys;
            sys.set_gravity(vec2dp{0.0, 0.0, 0.0});
            sys.add_frame(static_frame2dp("W"));
            sys.add_body(static_frame2dp("hub", vec2dp{0.0, 0.0, 1.0}, 0.0), hub,
                         kin_state2dp{.vel = vec2dp{0.1, -0.05, 0.0}, .omega = 0.3});
            sys.add_revolute_body(static_frame2dp("A1", vec2dp{1.0, 0.0, 1.0}, 0.0), link,
                                  vec2dp{-0.4, 0.0, 1.0}, 0.3, 1.5, sys.index_of("hub"));
            sys.add_revolute_body(static_frame2dp("A2", vec2dp{0.8, 0.0, 1.0}, 0.0), link,
                                  vec2dp{-0.4, 0.0, 1.0}, -0.5, -2.0, sys.index_of("A1"));
            sys.add_revolute_body(static_frame2dp("B1", vec2dp{-1.0, 0.0, 1.0}, 0.0),
                                  link, vec2dp{0.4, 0.0, 1.0}, -0.2, -1.0,
                                  sys.index_of("hub"));
            sys.add_revolute_body(static_frame2dp("B2", vec2dp{-0.8, 0.0, 1.0}, 0.0),
                                  link, vec2dp{0.4, 0.0, 1.0}, 0.4, 1.2,
                                  sys.index_of("B1"));
            return sys;
        };
        dynamic_system2dp sys = build();
        CHECK(sys.dof_coords().size() == 7); // 3 + 4

        auto total_momentum = [&] {
            bivec2dp P{};
            for (size_t i = 1; i < sys.size(); ++i)
                P = P + sys.momentum_world(i);
            return P;
        };
        bivec2dp const P0 = total_momentum();
        value_t const E0 = sys.total_energy();
        value_t const P0n = std::sqrt(P0.x * P0.x + P0.y * P0.y + P0.z * P0.z);
        value_t max_dP = 0.0, max_dE = 0.0;
        for (int i = 0; i < 3000; ++i) {
            sys.step(1.0e-3);
            bivec2dp const d = total_momentum() - P0;
            max_dP = std::max(max_dP, std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z));
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));
        }
        fmt::println("  |P| = {:.4f}: max |dP|/|P| = {:.2e}, dE/E = {:.2e}", P0n,
                     max_dP / P0n, max_dE / std::abs(E0));
        CHECK(max_dP / P0n < 1e-9);
        CHECK(max_dE / std::abs(E0) < 1e-9);
        bivec2dp const dh = sys.momentum_world(sys.index_of("hub")) - P0;
        CHECK(std::sqrt(dh.x * dh.x + dh.y * dh.y) > 1e-3); // the hub reacted

        // Conservation alone does not see the root's pose retraction (in the plane a
        // pose error leaves momentum and energy intact), so the hub's pose at t = 1 s
        // must also converge at FOURTH order under dt-halving -- the dexp^-1 terms
        // in the motor-joint integrator are what make it do so (2nd order without)
        auto hub_at = [&](value_t dt) {
            dynamic_system2dp s = build();
            for (int i = 0; i < int(std::round(1.0 / dt)); ++i)
                s.step(dt);
            return unitize(move2dp(O_2dp, s.get_pos_trafo(s.index_of("hub"), 0)));
        };
        vec2dp const h1 = hub_at(4.0e-3), h2 = hub_at(2.0e-3), h3 = hub_at(1.0e-3);
        value_t const e1 = std::abs(h1.x - h2.x) + std::abs(h1.y - h2.y);
        value_t const e2 = std::abs(h2.x - h3.x) + std::abs(h2.y - h3.y);
        fmt::println("  hub pose under dt-halving: {:.2e} -> {:.2e}, ratio {:.1f}", e1,
                     e2, e1 / e2);
        CHECK(e1 / e2 > 14.0);
        fmt::println("");
    }

} // TEST_SUITE("PGA2DP: physics tests implementation")

/////////////////////////////////////////////////////////////////////////////////////////
// Geodetic coordinates in the meridian section (the 2D case of ga_usr_geodesics.hpp)
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA2DP: geodesics in the meridian section")
{

    TEST_CASE("pga2dp: meridian section vs. the 3D case at longitude zero")
    {
        fmt::println("pga2dp: meridian section vs. the 3D case at longitude zero");

        auto const near2 = [](vec2dp const& a, vec2dp const& b, value_t tol) {
            return std::abs(a.x - b.x) < tol && std::abs(a.y - b.y) < tol &&
                   std::abs(a.z - b.z) < tol;
        };

        /////////////////////////////////////////////////////////////////////////////////
        // the claim the shared core rests on: 3D geodesy IS 2D geodesy plus a rotation
        /////////////////////////////////////////////////////////////////////////////////
        //
        // Set the longitude to zero in three dimensions and the meridian plane is the
        // e1-e3 plane. Everything the 2D case computes must then reproduce the 3D result
        // component for component, with (e1, e2) standing for (X, Z).

        for (int ilat = -90; ilat <= 90; ilat += 5) {
            for (auto const h : {0.0, 250.0, 8848.0}) {

                auto const p2 = geo_pos2dp{deg2rad(double(ilat)), h};
                auto const p3 = geo_pos{deg2rad(double(ilat)), 0.0, h};

                auto const P2 = geo_to_ecef(p2);
                auto const P3 = geo_to_ecef(p3);

                // the position: (r, z) is (X, Z), and Y vanishes on the prime meridian
                CHECK(P2.x == doctest::Approx(P3.x).epsilon(1e-13));
                CHECK(P2.y == doctest::Approx(P3.z).epsilon(1e-13));
                CHECK(std::abs(P3.y) < 1e-9);

                // ... and so does the distance from the geocenter
                CHECK(distance_from_geocenter(p2) ==
                      doctest::Approx(distance_from_geocenter(p3)).epsilon(1e-13));

                // the local frame: 2D up/north are the 3D up/north with Y dropped
                auto const F2 = un_at(p2);
                auto const F3 = enu_at(p3);

                CHECK(F2.up.x == doctest::Approx(F3.up.x).scale(1.0).epsilon(1e-12));
                CHECK(F2.up.y == doctest::Approx(F3.up.z).scale(1.0).epsilon(1e-12));
                CHECK(F2.north.x ==
                      doctest::Approx(F3.north.x).scale(1.0).epsilon(1e-12));
                CHECK(F2.north.y ==
                      doctest::Approx(F3.north.z).scale(1.0).epsilon(1e-12));
                CHECK(std::abs(F3.up.y) < 1e-12);    // the 3D frame's out-of-plane
                CHECK(std::abs(F3.north.y) < 1e-12); // components vanish here
            }
        }

        /////////////////////////////////////////////////////////////////////////////////
        // the position, pinned against an independently computed reference
        /////////////////////////////////////////////////////////////////////////////////

        auto const Berlin = to_geo_pos(geo_pos_dms2dp{"52°31'12\"N", 35});

        CHECK(Berlin.lat == doctest::Approx(deg2rad(52.52)).epsilon(1e-14));
        CHECK(Berlin.height == 35.0);

        auto const P_B = geo_to_ecef(Berlin);
        CHECK(P_B.x == doctest::Approx(3889225.181433).epsilon(1e-11));
        CHECK(P_B.y == doctest::Approx(5038246.874625).epsilon(1e-11));
        CHECK(P_B.z == 1.0); // a point, not a direction
        CHECK(distance_from_geocenter(Berlin) ==
              doctest::Approx(6364746.977026).epsilon(1e-13));

        // the two analytic pins: the equator is the equatorial radius along e1, the pole
        // the polar radius along e2
        CHECK(near2(geo_to_ecef(geo_pos2dp{0.0, 0.0}), vec2dp{wgs84.r_equator, 0.0, 1.0},
                    1e-9));
        CHECK(near2(geo_to_ecef(geo_pos2dp{0.5 * pi, 0.0}),
                    vec2dp{0.0, wgs84.r_pole, 1.0}, 1e-6));

        // the inverse inverts, poles included
        for (int ilat = -90; ilat <= 90; ilat += 5) {
            auto const p = geo_pos2dp{deg2rad(double(ilat)), 250.0};
            auto const q = ecef_to_geo(geo_to_ecef(p));
            CHECK(q.lat == doctest::Approx(p.lat).scale(1.0).epsilon(1e-12));
            CHECK(q.height == doctest::Approx(p.height).scale(1.0).epsilon(1e-9));
        }

        // the far half of the ellipse is the 180-degree meridian and is refused rather
        // than folded onto a wrong latitude
        CHECK_THROWS_AS(ecef_to_geo(vec2dp{-6.0e6, 1.0e6, 1.0}), std::invalid_argument);

        /////////////////////////////////////////////////////////////////////////////////
        // the local frame -- and why its axis order is (up, north)
        /////////////////////////////////////////////////////////////////////////////////

        auto const F = un_at(Berlin);

        CHECK(near2(F.up, vec2dp{0.608484459368082, 0.793565789778978, 0.0}, 1e-14));
        CHECK(near2(F.north, vec2dp{-0.793565789778978, 0.608484459368082, 0.0}, 1e-14));

        // the geometric and the trigonometric construction agree everywhere. Note the
        // sweep runs to the poles: unlike 3D there is NO degenerate configuration here,
        // because the degeneracy in 3D comes from the longitude, which does not exist in
        // a meridian section.
        for (int ilat = -90; ilat <= 90; ilat += 5) {
            auto const p = geo_pos2dp{deg2rad(double(ilat)), 0.0};
            CHECK(near2(un_at(p).up, un_basis_at(p).up, 1e-13));
            CHECK(near2(un_at(p).north, un_basis_at(p).north, 1e-13));
        }
        CHECK_THROWS_AS(un_at(vec2dp{0.0, 0.0, 1.0}), std::invalid_argument);

        // orthonormal, and ORIENTED: (up, north) is the positive pair, (north, up) the
        // negative one. That is what forces the axis order -- a motor is a rotation, so
        // it could never carry (e1, e2) onto (north, up), which is a reflection.
        auto const u2 = vec2d{F.up.x, F.up.y};
        auto const n2 = vec2d{F.north.x, F.north.y};
        CHECK(nrm(u2) == doctest::Approx(1.0).epsilon(1e-14));
        CHECK(nrm(n2) == doctest::Approx(1.0).epsilon(1e-14));
        CHECK(value_t(dot(u2, n2)) == doctest::Approx(0.0).scale(1.0).epsilon(1e-14));
        CHECK(value_t(wdg(u2, n2)) == doctest::Approx(1.0).epsilon(1e-14));  // +e12
        CHECK(value_t(wdg(n2, u2)) == doctest::Approx(-1.0).epsilon(1e-14)); // -e12

        /////////////////////////////////////////////////////////////////////////////////
        // the motor: one rotation, and it is the latitude itself
        /////////////////////////////////////////////////////////////////////////////////

        auto const M = un_to_ecef_motor(Berlin);

        CHECK(near2(unitize(move2dp(O_2dp, M)), P_B, 1e-6));
        CHECK(near2(move2dp(x_dir_2dp, M), F.up, 1e-13));    // local e1 -> up
        CHECK(near2(move2dp(y_dir_2dp, M), F.north, 1e-13)); // local e2 -> north

        // the rotation part IS the rotation by the latitude about the geocenter
        CHECK(is_same_motion(rgpr(get_motor(vec2dp{-P_B.x, -P_B.y, 0.0}), M),
                             get_motor(O_2dp, Berlin.lat), 1e-9));

        auto const M_inv = rrev(M);
        CHECK(near2(unitize(move2dp(P_B, M_inv)), O_2dp, 1e-6));

        // a point 100 m straight up sits at (0, 100) locally
        auto const P_up = vec2dp{P_B.x + 100.0 * F.up.x, P_B.y + 100.0 * F.up.y, 1.0};
        CHECK(near2(unitize(move2dp(P_up, M_inv)), vec2dp{100.0, 0.0, 1.0}, 1e-6));

        // ... which is the LOCAL x-axis, since up is the first axis here
        CHECK(near2(unitize(move2dp(geo_to_ecef(geo_pos2dp{Berlin.lat, 135.0}), M_inv)),
                    vec2dp{100.0, 0.0, 1.0}, 1e-6));

        // the point-based motor agrees with the position-based one over the sweep
        for (int ilat = -90; ilat <= 90; ilat += 10) {
            auto const p = geo_pos2dp{deg2rad(double(ilat)), 250.0};
            CHECK(is_same_motion(un_motor_at(geo_to_ecef(p)), un_to_ecef_motor(p), 1e-9));
        }

        /////////////////////////////////////////////////////////////////////////////////
        // the same station read in 2D and in 3D
        /////////////////////////////////////////////////////////////////////////////////

        // a second station in the same meridian, seen from Berlin: the local coordinates
        // must match the 3D answer at longitude zero, component for component
        {
            auto const Cairo2 = to_geo_pos(geo_pos_dms2dp{"30°02'N", 23});
            auto const Cairo3 = geo_pos{Cairo2.lat, 0.0, Cairo2.height};
            auto const Berlin3 = geo_pos{Berlin.lat, 0.0, Berlin.height};

            auto const q2 = unitize(move2dp(geo_to_ecef(Cairo2), M_inv));
            auto const q3 =
                unitize(move3dp(geo_to_ecef(Cairo3), rrev(enu_to_ecef_motor(Berlin3))));

            // 2D (up, north) against 3D (east, north, up): up <-> z, north <-> y
            CHECK(q2.x == doctest::Approx(q3.z).epsilon(1e-9));
            CHECK(q2.y == doctest::Approx(q3.y).epsilon(1e-9));
            CHECK(std::abs(q3.x) < 1e-6); // nothing to the east, same meridian

            fmt::println("   Cairo seen from Berlin, meridian section (up, north) ="
                         " ({:.1f}, {:.1f}) km",
                         q2.x / 1e3, q2.y / 1e3);
            fmt::println("   the 3D answer at lon 0 (east, north, up)          ="
                         " ({:.1f}, {:.1f}, {:.1f}) km",
                         q3.x / 1e3, q3.y / 1e3, q3.z / 1e3);
        }

        /////////////////////////////////////////////////////////////////////////////////
        // the types print, and dms output closes the parser's round trip
        /////////////////////////////////////////////////////////////////////////////////

        geo_pos_dms2dp const Berlin_dms{"52°31'12\"N", 35};

        CHECK(fmt::format("{}", Berlin_dms) ==
              "geo_pos_dms2dp(lat = 52°31'12\"N, height = 35 m)");
        CHECK(fmt::format("{:.4f}", Berlin) ==
              "geo_pos2dp(lat = 0.9166 rad, height = 35.0000 m)");
        CHECK(fmt::format("{:.3f}", F).starts_with("un_frame(up = Vec2dp("));

        // deg2dms is the inverse of dms2deg: what was read can be written back
        CHECK(deg2dms(52.52, geo_angle::latitude) == "52°31'12.0\"N");
        CHECK(deg2dms(-(33.0 + 52.0 / 60.0), geo_angle::latitude, 0) == "33°52'00\"S");
        CHECK(deg2dms(-(3.0 + 43.0 / 60.0), geo_angle::longitude, 0) == "3°43'00\"W");
        CHECK(rad2dms(deg2rad(13.41), geo_angle::longitude) == "13°24'36.0\"E");

        // and the round trip holds to the requested resolution, both ways
        for (auto const& str :
             {"52°31'12.0\"N", "0°00'00.0\"N", "89°59'59.9\"N", "33°52'00.0\"S"}) {
            CHECK(deg2dms(dms2deg(str, geo_angle::latitude), geo_angle::latitude) == str);
        }
        for (value_t d = -89.0; d <= 89.0; d += 7.3) {
            auto const round =
                dms2deg(deg2dms(d, geo_angle::latitude, 4), geo_angle::latitude);
            CHECK(round == doctest::Approx(d).scale(1.0).epsilon(1e-11));
        }

        // rounding carries instead of printing an impossible 60"
        CHECK(deg2dms(52.0 + 59.0 / 60.0 + 59.99 / 3600.0, geo_angle::latitude, 1) ==
              "53°00'00.0\"N");

        fmt::println("   {}", Berlin_dms);
        fmt::println("   {:.6f}", Berlin);
        fmt::println("   written back out: {} / {}",
                     rad2dms(Berlin.lat, geo_angle::latitude, 3),
                     rad2dms(deg2rad(13.41), geo_angle::longitude, 3));

        fmt::println("");
    }

} // TEST_SUITE("PGA2DP: geodesics in the meridian section")
