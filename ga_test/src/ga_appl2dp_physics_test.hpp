// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <cmath>
#include <iostream>

#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

#include "hd/hd_determinant.hpp"

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

        auto det_I_X1 = hd::det(I_X1.view());
        auto det_I_X2 = hd::det(I_X2.view());
        auto det_I_tot = hd::det(I_tot.view());

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

        // auto v = vec2dp{1, 1, 0};   // v_x = 1 m/s, v_y = 1 m/s v' = |v|*sqrt(2))
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
        fmt::println("M(1s) = exp(0.5 * B_1s) = {:>-7.3f}", M_1s);
        fmt::println("M(1s)_alt               = {:>-7.3f}", M_1s_alt);
        fmt::println("M(2s) = exp(0.5 * B_2s) = {:>-7.3f}", M_2s);
        fmt::println("M(2s)_alt               = {:>-7.3f}", M_2s_alt);
        fmt::println("");
        fmt::println("gr1(rgpr(rgpr(M, v), rrev(M)))) = {:>-7.3f}", move2dp(X1, M_1s));
        fmt::println("gr1(rgpr(rgpr(M, v), rrev(M)))) = {:>-7.3f}", move2dp(X1, M_2s));
        fmt::println("");

        CHECK(X11 == move2dp(X1, M_1s));
        CHECK(X12 == move2dp(X1, M_2s));
        CHECK(M_1s == exp(0.5 * B_1s));
        CHECK(M_2s == exp(0.5 * B_2s));

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

        auto const det_I_grid = hd::det(I_grid.view());
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
                //  typical starting value is "no initial transformation":
                //  B0=(0,0,0) => M0 = exp(0.5 * B0) = pscalar2dp(1)
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
                    auto M = exp(0.5 * B);
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
                //  typical starting value is "no initial transformation":
                //  B0=(0,0,0) => M0=0.5 * B0) = pscalar(1)
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
                [[maybe_unused]] vec2dp B = u[0]; // position transformation B is in u[n,0]
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

                    // calculate current position from B via M = exp(0.5 * B) ⟇ M0
                    // and via pts(t) = M ⟇ pts(t0) ⟇ rrev(M) = move2dp(pts(t0),M)
                    // and d(pts(t))/dt = rcmt(Omega, pts(t))

                    auto M = exp(0.5 * B);
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
    //   M(t)  - current motor: M0 ⟇ exp(½ B_b(t))  [body-frame formulation]
    //   P_b   - a fixed point in the body frame
    //   P_w   - P_b expressed in world frame: move2dp(P_b, M)
    //
    // Key identity that reveals the bug in the original code:
    //   move2dp(B, exp(½B)) == B   always (B commutes with its own exponential)
    //   => B_w differs from B_b only when M0 != ps (non-trivial initial state)
    //
    // Extensibility: body-frame velocity vectors v_b = vec2dp{vx, vy, 0} (z=0)
    //   transform to world via the same sandwich:  v_w = move2dp(v_b, M)
    ///////////////////////////////////////////////////////////////////////

    // Helper: print one time-step for kinematic cases
    auto print_kine_step = [](double t, mvec2dp_u const& /*M0*/, vec2dp const& B_b,
                              vec2dp const& B_w, mvec2dp_u const& /*M*/, vec2dp const& P_b,
                              vec2dp const& P_w, vec2dp const& /*Q_b*/,
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
            mvec2dp_u const M0 = exp(0.5 * vec2dp{0.0, 0.0, 0.0}); // ps (identity)
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, exp(0.5 * B_b));
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
            mvec2dp_u const M0 = exp(0.5 * vec2dp{-cm0.y, cm0.x, 0.0});
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, exp(0.5 * B_b));
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

        double const T = 2.0 * M_PI; // one full revolution
        double const dt = T / 8.0;
        int const N = 8;

        fmt::println("  --- Case A: M0 = identity (Q_b is at world origin) ---");
        {
            mvec2dp_u const M0 = exp(0.5 * vec2dp{0.0, 0.0, 0.0});
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};
            // pivot world position (should stay fixed)
            [[maybe_unused]] vec2dp const Q_world_expected = move2dp(Q_b, M0);

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, exp(0.5 * B_b));
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
            mvec2dp_u const M0 = exp(0.5 * vec2dp{-trans.y, trans.x, 0.0});
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, exp(0.5 * B_b));
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
            mvec2dp_u const M0 = exp(0.5 * vec2dp{0.0, 0.0, 0.0});
            vec2dp B_b{0.0, 0.0, 0.0};
            vec2dp B_w{0.0, 0.0, 0.0};

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, exp(0.5 * B_b));
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
            // By the conjugation identity: M0_tra ⟇ exp(½ phi0 Q_b)
            //   == (rotation by phi0 about Q_world) ⟇ M0_tra
            // so this correctly places the plate at the phi0-rotated rest position.
            vec2dp const cm0{3.0, 2.0, 1.0}; // rest cm world position (phi=0)
            mvec2dp_u const M0 =
                exp(0.5 * vec2dp{-cm0.y, cm0.x, 0.0}); // pure translation

            double const phi0 = M_PI / 4.0;
            vec2dp B_b = phi0 * Q_b; // initial B_b: phi0 rotation about body-frame Q_b
            vec2dp B_w = phi0 * Q_b; // B_w starts equal to B_b (same reference pose)

            for (int i = 0; i <= N; ++i) {
                double const t = i * dt;
                mvec2dp_u const M = rgpr(M0, exp(0.5 * B_b));
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
                // Encoding: translation by (tx, ty) -> motor exp(0.5*(-ty, tx, 0))
                M0 = exp(0.5 * vec2dp{-cm_w_pos0.y, cm_w_pos0.x, 0.0});

                // pivot_w: world position of body-frame pivot Q_b via M0 (stays fixed)
                vec2dp const Q_b{hw, hh, 1.0};
                pivot_w = move2dp(Q_b, M0);
                fmt::println("pivot_w = {:>-7.3f}  (world position of pivot Q_b via M0)",
                             pivot_w);

                // B_b(0): initial rotation by phi0 about body-frame pivot Q_b.
                // move2dp(Q_b, exp(½*phi*Q_b)) = Q_b for all phi (pivot invariant).
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

                // current motor: M(t) = M0 ⟇ exp(½ B_b(t))  [body-frame formulation]
                auto const M = rgpr(M0, exp(0.5 * B));

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
                // In principle correct without projection: rhs[1] = compute_omega_dot(I_inv, F_b, Omega, I);
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

                // current motor: M(t) = M0 ⟇ exp(½ B_b(t))  [body-frame formulation]
                auto const M = rgpr(M0, exp(0.5 * B));

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
                auto const M = rgpr(M0, exp(0.5 * u[0]));
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

            // Pivot constraint: distance from pivot (0,0) to cm = sqrt(hw^2+hh^2)
            value_t const r_sq = cm_w.x * cm_w.x + cm_w.y * cm_w.y;
            CHECK(r_sq == doctest::Approx(pivot_dist_sq).epsilon(1e-3));
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
    // save a version for multiple points, but still symmetric,
    // i.e. reducable to one eqn. via symmetry of the internal forces
    // typically done via reduces masses and forces by weighting with position
    //////////////////////////////////////

    // TEST_CASE("ega3d: circular motion II (symmetric multi-body system)")
    // {
    //     fmt::println("ega3d: circular motion II (symmetric multi-body system)");

    //     // TODO: translate to system of symmetric central forces,
    //     //       e.g. NFCM, p. 230 (two particle system)
    //     //
    //     //       and then reduce again to only use center of mass and orientaion

    //     class sim_ode_multibody_ega3d { // model 2nd order ode by a 1st order system
    //                                     // integrate a -> v -> s,
    //                                     // i.e. from acceleration via velocity to
    //                                     position

    //       public:

    //         sim_ode_multibody_ega3d(size_t npts_in) :
    //             npts(npts_in), u_mem(npts * 2), uh_mem(2 * npts * 2), rhs_mem(npts * 2)
    //         {
    //         }

    //         void set_initial_values()
    //         {

    //             // Create mdspan view for setting initial values
    //             auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);

    //             // use vec3d to encode 2d position/velocity in x- and y-component
    //             // of 3d vector and angle/angular velocity in z-component of 3d vector

    //             // Set initial state for all n points n = [0,n) and their components:
    //             // u[n,0] = position (x  ,   y,   phi)
    //             // u[n,1] = velocity (v_x, v_y, omega)

    //             size_t n = 0;
    //             u[n, 0] = vec3d{0, 1, 0};     // initial position
    //             u[n, 1] = vec3d{1, 0, 0.0}    // circumferentential initial velocity
    //                       + vec3d{0, 0.1, 0}; // velocity away from cg
    //                                           //    + vec3d{1, 0, 0}; // velocity of cg

    //             n = 1;
    //             u[n, 0] = vec3d{0, -1, 0};     // initial position
    //             u[n, 1] = vec3d{-1, 0, 0.0}    // circumferentential initial velocity
    //                       + vec3d{0, -0.1, 0}; // velocity away from cg
    //                                            //   + vec3d{1, 0, 0}; // velocity of cg

    //             // we ignore phi and omega, since we only have internal, central forces
    //         }

    //         void calc_rhs()
    //         {

    //             // Create mdspan views for rhs calculation
    //             auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);
    //             auto rhs = mdspan<vec3d, dextents<size_t, 2>>(rhs_mem.data(), npts, 2);


    //             double c_sp = 46.3;     // spring constant in N/m  => F_sp = -c*(l-l0)
    //             double l0 = 2.0 * 0.95; // initial spring length in m


    //             // get current state
    //             vec3d position0 = u[0, 0]; // position is stored as point in u[n,0]
    //             vec3d velocity0 = u[0, 1]; // velocity is stored as vector in u[n,1]

    //             vec3d position1 = u[1, 0];
    //             vec3d velocity1 = u[1, 1];

    //             // calculate forces and torques to change linear and angular
    //             // acceleration

    //             // HINT: this implements a multi-body solution with internal forces,
    //             //       would not be required to this extend for a rigid body system

    //             // equilibrium force for circular motion at given v and r:
    //             // F_eq = m * a_r = m * r * omega^2 = m * r * (v / r)^2 = m * v^2 / r
    //             // F_sp = -c_sp * (l - l0) = -c_sp *(2*r - l0)
    //             //
    //             // T = 2*pi/omega = 2*pi/(v/r) = 2*pi*r/v

    //             // keep difference to origin for comparison with 2dp case later on
    //             auto r_hat = position1 - position0;
    //             auto r = nrm(r_hat);
    //             r_hat /= r; // normalize r_hat

    //             auto acceleration0 = c_sp / m * (r - l0) * r_hat;
    //             auto acceleration1 = -c_sp / m * (r - l0) * r_hat;


    //             // J = m * r^2 with r being the distance from rot-axis
    //             // from dL/dt = J*omega = M = 0 (no torque) => alpha = - dJ/dt*omega/J
    //             // with dJ/dt = dJ/dr*dr/dt = dJ/dr*v_r = 2*m*r*v_r
    //             // alpha = -2*omega*v_r/r  (=angular acceleration)
    //             acceleration0.z = 0.0; // angular acceleration
    //             acceleration1.z = 0.0; // angular acceleration
    //                                    // (implicitly covered in coupled equations)

    //             // Set right-hand side for ODE system
    //             // u[n,0]' = velocity (linear and angular)
    //             // u[n,1]' = acceleration (linear and angular)
    //             rhs[0, 0] = velocity0;
    //             rhs[0, 1] = acceleration0;

    //             rhs[1, 0] = velocity1;
    //             rhs[1, 1] = acceleration1;
    //         }

    //         void calc_rkstep(double dt)
    //         {
    //             // Create mdspan views for RK4 integration (for all npts points)
    //             auto u = mdspan<vec3d, dextents<size_t, 1>>(u_mem.data(), npts * 2);
    //             auto uh = mdspan<vec3d, dextents<size_t, 2>>(uh_mem.data(), 2, npts *
    //             2); auto rhs =
    //                 mdspan<vec3d const, dextents<size_t, 1>>(rhs_mem.data(), npts * 2);

    //             // Perform RK4 integration (4 sub-steps)
    //             for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
    //                 calc_rhs();
    //                 rk4_step(u, uh, rhs, dt, rk_step);
    //             }
    //         }

    //         void print_sim(double t)
    //         {
    //             // Create mdspan view for printing
    //             auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);


    //             fmt::println("t = {:>-7.3f}, dist = nrm(r1 - r0) = {}:", t,
    //                          nrm(u[1, 0] - u[2, 0]));

    //             for (size_t n = 0; n < npts; ++n) {

    //                 // get current state
    //                 vec3d position = u[n, 0];
    //                 vec3d velocity = u[n, 1];

    //                 fmt::println("    n = {}, pos. = {:>-7.3f}, vel. = {:>-7.3f}", n,
    //                              position, velocity);
    //             }
    //         }

    //       private:

    //         // number of independent points to solve system for
    //         size_t npts;

    //         // RK4 integration state for point n with system order = 2
    //         // => [n+0: position, n+1: velocity])
    //         std::vector<vec3d> u_mem;   // [n+0: position, n+1: velocity]
    //         std::vector<vec3d> uh_mem;  // helper for integration
    //         std::vector<vec3d> rhs_mem; // right-hand side values

    //         // mass (= linear inertia)
    //         double m = 1.0;
    //     };

    //     auto constexpr num_points = 2;
    //     sim_ode_multibody_ega3d sim(num_points);

    //     // time range (1 full revolution)
    //     auto t_rng = discrete_range(0.0, 6.2832, 60);

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

    //                 // calculate current position from B via M = exp(0.5 * B) ⟇ M0
    //                 // and via pts(t) = M ⟇ pts(t0) ⟇ rrev(M)
    //                 //
    //                 auto M = exp(0.5 * B);
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

} // TEST_SUITE("PGA2DP: physics tests implementation")