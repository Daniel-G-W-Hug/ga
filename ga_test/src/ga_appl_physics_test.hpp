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
// PGA2DP physics tests - Inertia matrix for rigid body dynamics
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA2DP: physics tests")
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

        Vec2dp<double> result = I(Omega);

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

        // Zero force, zero omega -> zero omega_dot
        Vec2dp<double> F_zero{0.0, 0.0, 0.0};
        Vec2dp<double> Omega_zero{0.0, 0.0, 0.0};

        auto Omega_dot = compute_omega_dot(I_inv, F_zero, Omega_zero, I);

        CHECK(Omega_dot.x == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(Omega_dot.y == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(Omega_dot.z == doctest::Approx(0.0).epsilon(1e-10));

        // Non-zero force, zero omega -> non-zero omega_dot (linear response)
        Vec2dp<double> F{0.0, 0.0, 1.0}; // torque about origin
        Omega_dot = compute_omega_dot(I_inv, F, Omega_zero, I);

        // Result should be I_inv(F) since rcmt(0, ...) = 0
        auto expected = I_inv(F);
        CHECK(Omega_dot.x == doctest::Approx(expected.x).epsilon(1e-10));
        CHECK(Omega_dot.y == doctest::Approx(expected.y).epsilon(1e-10));
        CHECK(Omega_dot.z == doctest::Approx(expected.z).epsilon(1e-10));
    }

    TEST_CASE("pga2dp: compute_motor_dot - motor derivative")
    {
        fmt::println("pga2dp: compute_motor_dot - motor derivative");

        // Identity motor (no rotation, no translation)
        MVec2dp_U<double> M{Vec2dp<double>{0.0, 0.0, 0.0}, PScalar2dp<double>{1.0}};

        // Zero rate of change -> zero motor derivative
        Vec2dp<double> Omega_zero{0.0, 0.0, 0.0};
        auto M_dot = compute_motor_dot(M, Omega_zero);

        // All components should be zero
        CHECK(M_dot.c0 == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(M_dot.c1 == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(M_dot.c2 == doctest::Approx(0.0).epsilon(1e-10));
        CHECK(M_dot.c3 == doctest::Approx(0.0).epsilon(1e-10));
    }

} // TEST_SUITE("PGA2DP: physics tests")


/////////////////////////////////////////////////////////////////////////////////////////
// PGA3DP physics tests - Inertia matrix for rigid body dynamics
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: physics tests")
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

    TEST_CASE("pga3dp: compute_motor_dot - motor derivative")
    {
        fmt::println("pga3dp: compute_motor_dot - motor derivative");

        // Identity motor (no rotation, no translation)
        // MVec3dp_E has components: c0=scalar, c1-c6=bivector, c7=pseudoscalar
        MVec3dp_E<double> M{Scalar3dp<double>{0.0},
                            BiVec3dp<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                            PScalar3dp<double>{1.0}};

        // Zero rate of change -> zero motor derivative
        BiVec3dp<double> Omega_zero{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        auto M_dot = compute_motor_dot(M, Omega_zero);

        double tol = 1e-10;
        // All components should be zero
        CHECK(M_dot.c0 == doctest::Approx(0.0).epsilon(tol));
        CHECK(M_dot.c1 == doctest::Approx(0.0).epsilon(tol));
        CHECK(M_dot.c2 == doctest::Approx(0.0).epsilon(tol));
        CHECK(M_dot.c3 == doctest::Approx(0.0).epsilon(tol));
        CHECK(M_dot.c4 == doctest::Approx(0.0).epsilon(tol));
        CHECK(M_dot.c5 == doctest::Approx(0.0).epsilon(tol));
        CHECK(M_dot.c6 == doctest::Approx(0.0).epsilon(tol));
        CHECK(M_dot.c7 == doctest::Approx(0.0).epsilon(tol));
    }

    TEST_CASE("ega3d for 2D case (force + moment in vec3d): linear motion")
    {
        fmt::println("ega3d for 2D case (force + moment in vec3d): linear motion");


        class sim_ode { // model 2nd order ode by a 1st order system

          public:

            sim_ode(size_t npts_in) :
                npts(npts_in), u_mem(npts * 2), uh_mem(2 * npts * 2), rhs_mem(npts * 2)
            {
            }

            void set_initial_values()
            {

                // Create mdspan view for setting initial values
                auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);

                // use vec3d to encode 2d position/velocity in x- and y-component
                // of 3d vector and angle/angular velocity in z-component of 3d vector

                // Set initial state for all n points n = [0,n):
                // u[n,0] = position
                // u[n,1] = velocity

                size_t n = 0;
                u[n, 0] = vec3d{0, 1, 0}; // initial position;
                u[n, 1] = vec3d{1, 0, 0}; // initial velocity;

                n = 1;
                u[n, 0] = vec3d{0, -1, 0}; // initial position;
                u[n, 1] = vec3d{1, 0, 0};  // initial velocity;
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
                    auto force = vec2d{0.0, 0.0}; // no force (x,y)
                    // auto force = vec2d{0.0, m * 9.81}; // gravitation (x=0, y=m*g)
                    auto torque = 0.0; // no torque (z)

                    vec3d acceleration;
                    acceleration.x = force.x / m; // linear acceleratin in (x,y)
                    acceleration.y = force.y / m;

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
        sim_ode sim(num_points);


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


    TEST_CASE("ega3d for 2D case (force + moment in vec3d): circular motion")
    {
        fmt::println("ega3d for 2D case (force + moment in vec3d): circular motion");


        class sim_ode { // model 2nd order ode by a 1st order system

          public:

            sim_ode(size_t npts_in) :
                npts(npts_in), u_mem(npts * 2), uh_mem(2 * npts * 2), rhs_mem(npts * 2)
            {
            }

            void set_initial_values()
            {

                // Create mdspan view for setting initial values
                auto u = mdspan<vec3d, dextents<size_t, 2>>(u_mem.data(), npts, 2);

                // use vec3d to encode 2d position/velocity in x- and y-component
                // of 3d vector and angle/angular velocity in z-component of 3d vector

                // Set initial state for all n points n = [0,n):
                // u[n,0] = position
                // u[n,1] = velocity

                size_t n = 0;
                u[n, 0] = vec3d{0, 1, 0}; // initial position;
                u[n, 1] = vec3d{1, 0, 0}; // initial velocity;

                n = 1;
                u[n, 0] = vec3d{0, -1, 0}; // initial position;
                u[n, 1] = vec3d{-1, 0, 0}; // initial velocity;
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
                    acceleration.x = force.x / m; // linear part in (x,y)
                    acceleration.y = force.y / m;

                    // J = m * r^2 with r being the distance from rot-axis
                    auto J = m * nrm_sq(position - O_3d);
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
        sim_ode sim(num_points);


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

} // TEST_SUITE("PGA3DP: physics tests")
