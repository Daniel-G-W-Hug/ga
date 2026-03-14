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

        fmt::println("X1                     = {:>-7.3f}", X1);
        fmt::println("v                      = {:>-7.3f}", v);
        fmt::println("omega_tra              = {:>-7.3f}", omega_tra);
        fmt::println("B(1s) = omega_tra * 1s = {:>-7.3f}", B_1s);
        fmt::println("B(2s) = omega_tra * 2s = {:>-7.3f}", B_2s);
        fmt::println("M(1s) = exp(B(1s))     = {:>-7.3f}", M_1s);
        fmt::println("M(1s)_alt              = {:>-7.3f}", M_1s_alt);
        fmt::println("M(2s) = exp(B(2s))     = {:>-7.3f}", M_2s);
        fmt::println("M(2s)_alt              = {:>-7.3f}", M_2s_alt);
        fmt::println("");
        fmt::println("gr1(rgpr(rgpr(M, v), rrev(M)))) = {:>-7.3f}", move2dp(X1, M_1s));
        fmt::println("gr1(rgpr(rgpr(M, v), rrev(M)))) = {:>-7.3f}", move2dp(X1, M_2s));
        fmt::println("");

        CHECK(X11 == move2dp(X1, M_1s));
        CHECK(X12 == move2dp(X1, M_2s));

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

} // TEST_SUITE("PGA2DP: physics tests prep")


/////////////////////////////////////////////////////////////////////////////////////////
// PGA3DP physics tests preparation - Inertia matrix for rigid body dynamics
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: physics tests prep")
{

    TEST_CASE("pga3dp: get_point_inertia - single point mass")
    {
        fmt::println("pga3dp: get_point_inertia - single point mass");

        // Point mass at (1, 0, 0) with mass 2.0
        double m = 2.0;
        Vec3dp<double> X{1.0, 0.0, 0.0, 1.0}; // unitized point at (1,0,0)

        auto I = get_point_inertia(m, X);
        auto v = I.view();

        // For X = (1, 0, 0, 1), the matrix simplifies significantly
        // Many terms become zero

        // Row 0: [0, 0, 0, 1, 0, 0] * m = [0, 0, 0, 2, 0, 0]
        CHECK(v[0, 0] == doctest::Approx(0.0));
        CHECK(v[0, 1] == doctest::Approx(0.0));
        CHECK(v[0, 2] == doctest::Approx(0.0));
        CHECK(v[0, 3] == doctest::Approx(2.0)); // Xw^2 * m
        CHECK(v[0, 4] == doctest::Approx(0.0));
        CHECK(v[0, 5] == doctest::Approx(0.0));

        // Row 1: [0, 0, 1, 0, 1, 0] * m = [0, 0, 2, 0, 2, 0]
        CHECK(v[1, 0] == doctest::Approx(0.0));
        CHECK(v[1, 1] == doctest::Approx(0.0));
        CHECK(v[1, 2] == doctest::Approx(2.0)); // Xx*Xw * m
        CHECK(v[1, 3] == doctest::Approx(0.0));
        CHECK(v[1, 4] == doctest::Approx(2.0)); // Xw^2 * m
        CHECK(v[1, 5] == doctest::Approx(0.0));

        // Row 2: [0, -1, 0, 0, 0, 1] * m = [0, -2, 0, 0, 0, 2]
        CHECK(v[2, 0] == doctest::Approx(0.0));
        CHECK(v[2, 1] == doctest::Approx(-2.0)); // -Xx*Xw * m
        CHECK(v[2, 2] == doctest::Approx(0.0));
        CHECK(v[2, 3] == doctest::Approx(0.0));
        CHECK(v[2, 4] == doctest::Approx(0.0));
        CHECK(v[2, 5] == doctest::Approx(2.0)); // Xw^2 * m
    }

    TEST_CASE("pga3dp: Inertia3dp operator+= - accumulate multiple point masses")
    {
        fmt::println("pga3dp: Inertia3dp operator+= - accumulate multiple point masses");

        // Three point masses
        double m1 = 1.0;
        Vec3dp<double> X1{1.0, 0.0, 0.0, 1.0};

        double m2 = 2.0;
        Vec3dp<double> X2{0.0, 1.0, 0.0, 1.0};

        double m3 = 1.5;
        Vec3dp<double> X3{0.0, 0.0, 1.0, 1.0};

        auto I1 = get_point_inertia(m1, X1);
        auto I2 = get_point_inertia(m2, X2);
        auto I3 = get_point_inertia(m3, X3);

        // Accumulate
        Inertia3dp<double> I_total{};
        I_total += I1;
        I_total += I2;
        I_total += I3;

        // Verify it equals sum of individual matrices
        for (size_t i = 0; i < 36; ++i) {
            CHECK(I_total.data[i] ==
                  doctest::Approx(I1.data[i] + I2.data[i] + I3.data[i]));
        }
    }

    TEST_CASE("pga3dp: Inertia3dp operator() - matrix-vector multiplication")
    {
        fmt::println("pga3dp: Inertia3dp operator() - matrix-vector multiplication");

        // Point at (1, 0, 0)
        double m = 1.0;
        Vec3dp<double> X{1.0, 0.0, 0.0, 1.0};

        auto I = get_point_inertia(m, X);

        // Apply to a simple rate of change
        // BiVec3dp: (vx=e41, vy=e42, vz=e43, mx=e23, my=e31, mz=e12)
        BiVec3dp<double> Omega{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}; // rotation about e23

        BiVec3dp<double> result = I(Omega);

        // Manual verification based on matrix structure
        // The result should be non-trivial
        // For X=(1,0,0,1), row 3 is [0, 0, 0, 0, 0, 0] (Xy^2+Xz^2=0)
        // So result.mx should be 0
        CHECK(result.mx == doctest::Approx(0.0));
    }

    TEST_CASE("pga3dp: get_inertia_inverse - verify I * I_inv = Identity")
    {
        fmt::println("pga3dp: get_inertia_inverse - verify I * I_inv = Identity");

        // Create a non-trivial inertia matrix from multiple points
        // Need at least 3 non-collinear points for full rank
        Inertia3dp<double> I{};
        I += get_point_inertia(1.0, Vec3dp<double>{1.0, 0.0, 0.0, 1.0});
        I += get_point_inertia(2.0, Vec3dp<double>{0.0, 1.0, 0.0, 1.0});
        I += get_point_inertia(1.5, Vec3dp<double>{0.0, 0.0, 1.0, 1.0});
        I += get_point_inertia(0.5, Vec3dp<double>{1.0, 1.0, 1.0, 1.0});

        auto I_inv = get_inertia_inverse(I);

        // Verify I * I_inv ≈ Identity by checking I(I_inv(e_i)) ≈ e_i
        BiVec3dp<double> e1{1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        BiVec3dp<double> e2{0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
        BiVec3dp<double> e3{0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
        BiVec3dp<double> e4{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        BiVec3dp<double> e5{0.0, 0.0, 0.0, 0.0, 1.0, 0.0};
        BiVec3dp<double> e6{0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

        auto r1 = I(I_inv(e1));
        auto r2 = I(I_inv(e2));
        auto r3 = I(I_inv(e3));
        auto r4 = I(I_inv(e4));
        auto r5 = I(I_inv(e5));
        auto r6 = I(I_inv(e6));

        double tol = 1e-10;

        CHECK(r1.vx == doctest::Approx(1.0).epsilon(tol));
        CHECK(r1.vy == doctest::Approx(0.0).epsilon(tol));
        CHECK(r1.vz == doctest::Approx(0.0).epsilon(tol));
        CHECK(r1.mx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r1.my == doctest::Approx(0.0).epsilon(tol));
        CHECK(r1.mz == doctest::Approx(0.0).epsilon(tol));

        CHECK(r2.vx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r2.vy == doctest::Approx(1.0).epsilon(tol));
        CHECK(r2.vz == doctest::Approx(0.0).epsilon(tol));
        CHECK(r2.mx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r2.my == doctest::Approx(0.0).epsilon(tol));
        CHECK(r2.mz == doctest::Approx(0.0).epsilon(tol));

        CHECK(r3.vx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r3.vy == doctest::Approx(0.0).epsilon(tol));
        CHECK(r3.vz == doctest::Approx(1.0).epsilon(tol));
        CHECK(r3.mx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r3.my == doctest::Approx(0.0).epsilon(tol));
        CHECK(r3.mz == doctest::Approx(0.0).epsilon(tol));

        CHECK(r4.vx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r4.vy == doctest::Approx(0.0).epsilon(tol));
        CHECK(r4.vz == doctest::Approx(0.0).epsilon(tol));
        CHECK(r4.mx == doctest::Approx(1.0).epsilon(tol));
        CHECK(r4.my == doctest::Approx(0.0).epsilon(tol));
        CHECK(r4.mz == doctest::Approx(0.0).epsilon(tol));

        CHECK(r5.vx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r5.vy == doctest::Approx(0.0).epsilon(tol));
        CHECK(r5.vz == doctest::Approx(0.0).epsilon(tol));
        CHECK(r5.mx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r5.my == doctest::Approx(1.0).epsilon(tol));
        CHECK(r5.mz == doctest::Approx(0.0).epsilon(tol));

        CHECK(r6.vx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r6.vy == doctest::Approx(0.0).epsilon(tol));
        CHECK(r6.vz == doctest::Approx(0.0).epsilon(tol));
        CHECK(r6.mx == doctest::Approx(0.0).epsilon(tol));
        CHECK(r6.my == doctest::Approx(0.0).epsilon(tol));
        CHECK(r6.mz == doctest::Approx(1.0).epsilon(tol));
    }

    TEST_CASE("pga3dp: compute_omega_dot - ODE right-hand side")
    {
        fmt::println("pga3dp: compute_omega_dot - ODE right-hand side");

        // Create inertia from multiple point masses
        Inertia3dp<double> I{};
        I += get_point_inertia(1.0, Vec3dp<double>{1.0, 0.0, 0.0, 1.0});
        I += get_point_inertia(1.0, Vec3dp<double>{-1.0, 0.0, 0.0, 1.0});
        I += get_point_inertia(1.0, Vec3dp<double>{0.0, 1.0, 0.0, 1.0});
        I += get_point_inertia(1.0, Vec3dp<double>{0.0, -1.0, 0.0, 1.0});

        auto I_inv = get_inertia_inverse(I);

        // Zero force, zero omega -> zero omega_dot
        BiVec3dp<double> F_zero{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        BiVec3dp<double> Omega_zero{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        auto Omega_dot = compute_omega_dot(I_inv, F_zero, Omega_zero, I);

        double tol = 1e-10;
        CHECK(Omega_dot.vx == doctest::Approx(0.0).epsilon(tol));
        CHECK(Omega_dot.vy == doctest::Approx(0.0).epsilon(tol));
        CHECK(Omega_dot.vz == doctest::Approx(0.0).epsilon(tol));
        CHECK(Omega_dot.mx == doctest::Approx(0.0).epsilon(tol));
        CHECK(Omega_dot.my == doctest::Approx(0.0).epsilon(tol));
        CHECK(Omega_dot.mz == doctest::Approx(0.0).epsilon(tol));

        // Non-zero force, zero omega -> non-zero omega_dot (linear response)
        BiVec3dp<double> F{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}; // torque about e12
        Omega_dot = compute_omega_dot(I_inv, F, Omega_zero, I);

        // Result should be I_inv(F) since rcmt(0, ...) = 0
        auto expected = I_inv(F);
        CHECK(Omega_dot.vx == doctest::Approx(expected.vx).epsilon(tol));
        CHECK(Omega_dot.vy == doctest::Approx(expected.vy).epsilon(tol));
        CHECK(Omega_dot.vz == doctest::Approx(expected.vz).epsilon(tol));
        CHECK(Omega_dot.mx == doctest::Approx(expected.mx).epsilon(tol));
        CHECK(Omega_dot.my == doctest::Approx(expected.my).epsilon(tol));
        CHECK(Omega_dot.mz == doctest::Approx(expected.mz).epsilon(tol));
    }

} // TEST_SUITE("PGA3DP: physics tests")


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
                cm_pos_init_w(cm_pos_in), cm_spd_init_w(cm_spd_in),
                cm_phi_init_w(cm_phi_in), cm_omega_init_w(cm_omega_in), u_mem(2),
                uh_mem(2 * 2), rhs_mem(2)
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
                //  B0=(0,0,0) => M0=exp(B0)=pscalar(1) (=identity transformation at t=0)
                //
                u[0] = vec2dp(0.0, 0.0, 0.0);

                // initial rate of change transformation of "velocity" dB(t0)/dt = Omega0
                // encoding:
                // rotation:    Omega0_rot = ( q.x0, q.y0, 1) * omega0
                // translation: Omega0_tra = (-v0.y, v0.x, 0)
                //
                // resulting Omega0 = Omega0_rot + Omega0_tra

                // case with initial translation
                u[1] = vec2dp(-cm_spd_init_w.y, cm_spd_init_w.x, cm_omega_init_w);
            }

            void calc_rhs()
            {

                // Create mdspan views for rhs calculation
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto rhs = mdspan<vec2dp, dextents<size_t, 1>>(rhs_mem.data(), 2);


                // get current state
                vec2dp B = u[0];     // position transformation B is in u[0]
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

                    // calculate current position from B via M = exp(B) ⟇ M0
                    // and via pts(t) = M ⟇ pts(t0) ⟇ rrev(M)
                    //
                    auto M = exp(B);
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
            vec2dp cm_pos_init_w; // initial position of center of mass (unitized)
            vec2dp cm_spd_init_w; // initial speed of center of mass (.z == 0)

            // initial angular position of cm (body vs. world frame, e12 defines pos phi)
            value_t cm_phi_init_w;

            // initial angular velocity of cm (body vs. world frame)
            value_t cm_omega_init_w;

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
                cm_pos_init_w(cm_pos_in), cm_spd_init_w(cm_spd_in),
                cm_phi_init_w(cm_phi_in), cm_omega_init_w(cm_omega_in), u_mem(2),
                uh_mem(2 * 2), rhs_mem(2)
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
                //  B0=(0,0,0) => M0=exp(B0)=pscalar(1) (=identity transformation at t=0)
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
                u[1] = wdg(O_2dp, scalar2dp(cm_omega_init_w));
                // fmt::println("n = {}, u = {}", n, u[n, 1]);
            }

            void calc_rhs()
            {

                // Create mdspan views for rhs calculation
                auto u = mdspan<vec2dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto rhs = mdspan<vec2dp, dextents<size_t, 1>>(rhs_mem.data(), 2);

                // get current state
                vec2dp B = u[0];     // position transformation B is in u[n,0]
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

                    // calculate current position from B via M = exp(B) ⟇ M0
                    // and via pts(t) = M ⟇ pts(t0) ⟇ rrev(M) = move2dp(pts(t0),M)
                    // and d(pts(t))/dt = rcmt(Omega, pts(t))

                    auto M = exp(B);
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
            vec2dp cm_pos_init_w; // initial position of center of mass (unitized)
            vec2dp cm_spd_init_w; // initial speed of center of mass (.z == 0)

            // initial angular position of cm (body vs. world frame, e12 defines pos phi)
            value_t cm_phi_init_w;

            // initial angular velocity of cm (body vs. world frame)
            value_t cm_omega_init_w;

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
    //             cm_pos_init_w(cm_pos_in), cm_spd_init_w(cm_spd_in),
    //             cm_phi_init_w(cm_phi_in), cm_omega_init_w(cm_omega_in), u_mem(npts *
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
    //                 // B0=(0,0,0) => M0=pscalar(1) (=identiy transformation at t=0)
    //                 u[n, 0] = vec2dp(0.0, 0.0, 0.0);

    //                 // initial rate of change transformation of "velocity" dB/dt =
    //                 Omega
    //                 // encoding:
    //                 // rotation: Omega_rot = (qx0, qy0, 1) * omega0
    //                 // translation: Omega_tra = (-v0.y, v0.x, 0)
    //                 //
    //                 // resulting Omega = Omega_rot + Omega_tra

    //                 // case with initial translation
    //                 u[n, 1] = vec2dp(-cm_spd_init_w.y, cm_spd_init_w.x,
    //                 cm_omega_init_w);
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

    //                 // calculate current position from B via M = exp(B) ⟇ M0
    //                 // and via pts(t) = M ⟇ pts(t0) ⟇ rrev(M)
    //                 //
    //                 auto M = exp(B);
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
    //         vec2dp cm_pos_init_w; // initial position of center of mass (unitized)
    //         vec2dp cm_spd_init_w; // initial speed of center of mass (.z == 0)

    //         // initial angular position of cm (body vs. world frame, e12 defines pos
    //         phi) value_t cm_phi_init_w;

    //         // initial angular velocity of cm (body vs. world frame)
    //         value_t cm_omega_init_w;

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

} // TEST_SUITE("PGA2DP: physics tests implementation")


/////////////////////////////////////////////////////////////////////////////////////////
// PGA3DP: discrete vs. continuous inertia
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: physics tests implementation")
{

    TEST_CASE("pga3dp: discrete inertia of rectangular cuboid (25x25x25 grid)")
    {
        fmt::println("pga3dp: discrete inertia of rectangular cuboid (25x25x25 grid)");
        fmt::println("");
        fmt::println("Cuboid: e1 in [-2, 2], e2 in [-1, 1], e3 in [-0.5, 0.5]");
        fmt::println("Total mass M = 1.0, uniform distribution");
        fmt::println("");

        // Grid parameters
        int const nx = 25;
        int const ny = 25;
        int const nz = 25;
        double const x_min = -2.0, x_max = 2.0; // w = 4
        double const y_min = -1.0, y_max = 1.0; // h = 2
        double const z_min = -0.5, z_max = 0.5; // d = 1
        double const dx = (x_max - x_min) / (nx - 1);
        double const dy = (y_max - y_min) / (ny - 1);
        double const dz = (z_max - z_min) / (nz - 1);

        // Total mass = 1.0 distributed uniformly over all grid points
        double const M_total = 1.0;
        double const m_pt = M_total / (nx * ny * nz); // mass per point

        // Accumulate inertia over all grid points
        Inertia3dp<double> I_grid{};
        for (int ix = 0; ix < nx; ++ix) {
            double const x = x_min + ix * dx;
            for (int iy = 0; iy < ny; ++iy) {
                double const y = y_min + iy * dy;
                for (int iz = 0; iz < nz; ++iz) {
                    double const z = z_min + iz * dz;
                    Vec3dp<double> const X{x, y, z, 1.0}; // unitized point (Xw = 1)
                    I_grid += get_point_inertia(m_pt, X);
                }
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

        // Analytical continuous limit for get_cuboid_inertia(M, w, h, d):
        //   w = 4 (e1), h = 2 (e2), d = 1 (e3)
        //   Upper-right block: I[0,3]=I[1,4]=I[2,5] = M  (mass, Newton p=mv)
        //   Lower-left diagonal (classical rectangle moments):
        //     I[3,0] = M*(h^2+d^2)/12 = (4+1)/12 = 5/12  ~ 0.41667
        //     I[4,1] = M*(w^2+d^2)/12 = (16+1)/12 = 17/12 ~ 1.41667
        //     I[5,2] = M*(w^2+h^2)/12 = (16+4)/12 = 5/3  ~ 1.66667
        double const W = x_max - x_min;                           // 4.0
        double const H = y_max - y_min;                           // 2.0
        double const D = z_max - z_min;                           // 1.0
        double const Ixx_cont = M_total * (H * H + D * D) / 12.0; // 5/12
        double const Iyy_cont = M_total * (W * W + D * D) / 12.0; // 17/12
        double const Izz_cont = M_total * (W * W + H * H) / 12.0; // 5/3

        // Exact continuous result from get_cuboid_inertia() - no discretization error
        auto const I_cuboid = get_cuboid_inertia(M_total, W, H, D);
        auto const I_cuboid_inv = get_inertia_inverse(I_cuboid);

        fmt::println("Continuous limit (uniform rectangular cuboid):");
        fmt::println("  I[0,3]=I[1,4]=I[2,5] = M      = {:.6f}  (discrete: {:.6f})",
                     M_total, I_grid.view()[0, 3]);
        fmt::println("  I[3,0] = M*(h^2+d^2)/12 = {:.6f}  (discrete: {:.6f})", Ixx_cont,
                     I_grid.view()[3, 0]);
        fmt::println("  I[4,1] = M*(w^2+d^2)/12 = {:.6f}  (discrete: {:.6f})", Iyy_cont,
                     I_grid.view()[4, 1]);
        fmt::println("  I[5,2] = M*(w^2+h^2)/12 = {:.6f}  (discrete: {:.6f})", Izz_cont,
                     I_grid.view()[5, 2]);
        fmt::println(
            "  (discrete overestimates moments by ~{:.1f}% due to endpoint sampling)",
            100.0 * (I_grid.view()[5, 2] - Izz_cont) / Izz_cont);
        fmt::println("");
        fmt::println("get_cuboid_inertia() (exact continuous limit):");
        fmt::println("I_cuboid     = {:>-10.6f}", I_cuboid);
        fmt::println("");
        fmt::println("I_cuboid_inv = {:>-10.6f}", I_cuboid_inv);
        fmt::println("");

        // get_cuboid_inertia() matches the continuous limit exactly
        // Upper-right block: mass terms
        CHECK(I_cuboid.view()[0, 3] == doctest::Approx(M_total).epsilon(1e-15));
        CHECK(I_cuboid.view()[1, 4] == doctest::Approx(M_total).epsilon(1e-15));
        CHECK(I_cuboid.view()[2, 5] == doctest::Approx(M_total).epsilon(1e-15));
        // Lower-left block: moments of inertia
        CHECK(I_cuboid.view()[3, 0] == doctest::Approx(Ixx_cont).epsilon(1e-15));
        CHECK(I_cuboid.view()[4, 1] == doctest::Approx(Iyy_cont).epsilon(1e-15));
        CHECK(I_cuboid.view()[5, 2] == doctest::Approx(Izz_cont).epsilon(1e-15));
        // All other entries must be zero
        for (size_t r = 0; r < 6; ++r) {
            for (size_t c = 0; c < 6; ++c) {
                bool const is_mass_entry = (r < 3 && c == r + 3);
                bool const is_moment_entry = (r >= 3 && c == r - 3);
                if (!is_mass_entry && !is_moment_entry) {
                    CHECK(I_cuboid.view()[r, c] == doctest::Approx(0.0).epsilon(1e-15));
                }
            }
        }

        // Discrete grid overestimates all moments vs. the exact continuous cuboid
        CHECK(I_grid.view()[3, 0] > I_cuboid.view()[3, 0]);
        CHECK(I_grid.view()[4, 1] > I_cuboid.view()[4, 1]);
        CHECK(I_grid.view()[5, 2] > I_cuboid.view()[5, 2]);

        // Exact properties of the discrete grid: mass terms and symmetry
        CHECK(I_grid.view()[0, 3] == doctest::Approx(M_total).epsilon(1e-10));
        CHECK(I_grid.view()[1, 4] == doctest::Approx(M_total).epsilon(1e-10));
        CHECK(I_grid.view()[2, 5] == doctest::Approx(M_total).epsilon(1e-10));

        // Exact discrete moments for 25-pt endpoint-inclusive grid:
        //   mean(x^2) for [-2,2]:   13/9
        //   mean(y^2) for [-1,1]:   13/36
        //   mean(z^2) for [-0.5,0.5]: 13/144
        //   I[3,0] = M*(mean(y^2)+mean(z^2)) = 13/36 + 13/144 = 65/144
        //   I[4,1] = M*(mean(x^2)+mean(z^2)) = 13/9  + 13/144 = 221/144
        //   I[5,2] = M*(mean(x^2)+mean(y^2)) = 13/9  + 13/36  = 65/36
        CHECK(I_grid.view()[3, 0] == doctest::Approx(65.0 / 144.0).epsilon(1e-6));
        CHECK(I_grid.view()[4, 1] == doctest::Approx(221.0 / 144.0).epsilon(1e-6));
        CHECK(I_grid.view()[5, 2] == doctest::Approx(65.0 / 36.0).epsilon(1e-6));
    }

} // TEST_SUITE("PGA3DP: physics tests implementation")
