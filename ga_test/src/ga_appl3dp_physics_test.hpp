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

    //////////////////////////////////////////////////////////////////////////////////////
    // PGA3DP: discrete vs. continuous inertia
    //////////////////////////////////////////////////////////////////////////////////////

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

} // TEST_SUITE("PGA3DP: physics tests prep")


/////////////////////////////////////////////////////////////////////////////////////////
// PGA3DP physics tests application - 3D cuboid pendulum
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: physics tests application")
{

    TEST_CASE("pga3dp: ODE cuboid pendulum (pivot edge at TR corner)")
    {
        fmt::println("pga3dp: ODE cuboid pendulum (pivot edge at TR corner)");

        class sim_ode_cuboid_pga3dp { // model 2nd order ode by a 1st order system

            // Models a rigid cuboid acting as a 3D pendulum fixed at its
            // top-right z-edge (the pivot line L_b).
            //
            // PGA3DP body-frame formulation: M(t) = M0 ⟇ exp(½ B_b(t))
            // - M0: pure translation motor placing body origin (= cm) at initial world pos
            //       M0 = exp(0.5 * bivec3dp{0,0,0,cx,cy,cz})
            //       (translation uses the bulk/e23,e31,e12 components, not e41,e42,e43)
            // - B_b = phi * L_b: rotation bivector, phi = rotation angle about L_b
            // - Pivot line L_b = wdg(P_z0, P_z1) = bivec3dp{0, 0, 1, hh, -hw, 0}
            //   where P_z0 = (hw,hh,0,1), P_z1 = (hw,hh,1,1) in body frame
            //   => direction vz = +1 (along e3/z axis), moment = (hh, -hw, 0)
            // - Constraint projection: alpha = dOmega.vz  (since L_b.vz = 1)
            // - Pivot invariance: move3dp(L_b, exp(½ phi * L_b)) = L_b for all phi
            //
            // get_cuboid_inertia(m,w,h,d,L_pivot) applies the scalar parallel-axis
            // (Steiner) correction: only I[3,0], I[4,1], I[5,2] are updated by
            // m*(Py²+Pz²), m*(Px²+Pz²), m*(Px²+Py²) respectively, where P_foot is
            // the foot of perpendicular from the body origin to L_pivot.

          public:

            sim_ode_cuboid_pga3dp(value_t m_in, value_t w_in, value_t h_in, value_t d_in,
                                  vec3dp const& cm_pos_in,
                                  value_t cm_phi_in, value_t cm_omega_in) :
                m(m_in), width(w_in), height(h_in), depth(d_in), cm_w_pos0(cm_pos_in),
                cm_w_phi0(cm_phi_in), cm_w_omega0(cm_omega_in),
                u_mem(2), uh_mem(2 * 2), rhs_mem(2)
            {
                fmt::println("sim_ode_cuboid_pga3dp: 3D pendulum (pivot = z-edge at TR).");

                // Body frame origin = CENTER OF MASS (cuboid symmetric about cm).
                // Pivot = z-edge through TR corner at (hw, hh, *) in body frame.
                // L_b = wdg(P_z0, P_z1) with P_z0=(hw,hh,0,1) and P_z1=(hw,hh,1,1)
                //      = bivec3dp{0, 0, 1, hh, -hw, 0}  (L_b.vz = 1, unit direction)
                value_t const hw = width / 2.0;
                value_t const hh = height / 2.0;
                bivec3dp const L_b{0.0, 0.0, 1.0, hh, -hw, 0.0};
                fmt::println("pivot line L_b (body frame, z-edge at TR corner):");
                fmt::println("  L_b = {:>-7.3f}", L_b);

                // get_cuboid_inertia with parallel-axis correction about L_b.
                fmt::println("inertia of cuboid in body frame (about pivot L_b):");
                I = get_cuboid_inertia(m, width, height, depth, L_b);
                I_inv = get_inertia_inverse(I);
                fmt::println("I     = {:>-7.3f}", I);
                fmt::println("I_inv = {:>-7.3f}", I_inv);
            }

            void set_initial_values()
            {
                auto u = mdspan<bivec3dp, dextents<size_t, 1>>(u_mem.data(), 2);

                value_t const hw = width / 2.0;
                value_t const hh = height / 2.0;

                // Body origin = CM. Pivot = z-edge at (hw, hh, *) in body frame.
                // M0: pure translation placing cm at cm_w_pos0.
                // In PGA3DP, translate by (cx,cy,cz): motor = exp(0.5 * bivec3dp{0,0,0,cx,cy,cz})
                // Translation is encoded in the bulk (e23,e31,e12) components,
                // not the weight (e41,e42,e43) components.
                M0 = exp(0.5 * bivec3dp{0.0, 0.0, 0.0,
                                        cm_w_pos0.x, cm_w_pos0.y, cm_w_pos0.z});

                // pivot_pt_w: world position of representative pivot point (hw,hh,0,1)
                // in body frame, mapped via M0.  Stays fixed during rotation since
                // move3dp(P_on_L_b, M0 ⟇ exp(½ phi * L_b)) = M0-image of P_on_L_b.
                pivot_pt_w = move3dp(vec3dp{hw, hh, 0.0, 1.0}, M0);
                fmt::println("pivot_pt_w = {:>-7.3f}  (world position of (hw,hh,0) via M0)",
                             pivot_pt_w);

                // B_b(0): initial rotation by phi0 about body-frame pivot line L_b.
                bivec3dp const L_b{0.0, 0.0, 1.0, hh, -hw, 0.0};
                u[0] = cm_w_phi0 * L_b;

                // Omega_b(0): initial angular velocity about pivot line L_b.
                u[1] = cm_w_omega0 * L_b;
            }

            void calc_rhs()
            {
                auto u   = mdspan<bivec3dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto rhs = mdspan<bivec3dp, dextents<size_t, 1>>(rhs_mem.data(), 2);

                bivec3dp const B     = u[0]; // position bivector B_b = phi * L_b
                bivec3dp const Omega = u[1]; // velocity bivector Omega_b = omega * L_b

                // Current motor: M(t) = M0 ⟇ exp(½ B_b(t))  [body-frame formulation]
                auto const M = rgpr(M0, exp(0.5 * B));

                // CM world position: body origin = (0,0,0,1) in body frame
                vec3dp const O_b{0.0, 0.0, 0.0, 1.0};
                auto const cm_w = move3dp(O_b, M);

                // Force couple in world frame (net force = 0, only torque):
                // - Gravity at cm (downward -y):       F_dn = wdg(cm_w, (0,-mg,0,0))
                // - Reaction at pivot (upward +y):     F_up = wdg(pivot_pt_w, (0,+mg,0,0))
                // Using a single representative pivot POINT (z=0 on pivot edge) here.
                // The torque from the reaction force about L_b vanishes by definition.
                auto const F_dn_w = wdg(cm_w,       vec3dp{0.0, -m * 9.81, 0.0, 0.0});
                auto const F_up_w = wdg(pivot_pt_w, vec3dp{0.0,  m * 9.81, 0.0, 0.0});
                auto const F_b = move3dp(F_dn_w + F_up_w, rrev(M)); // to body frame

                // Constrained rotation about pivot line L_b = {0,0,1, hh,-hw,0}.
                // I is about pivot L_b (parallel-axis corrected in ctor).
                // compute_omega_dot gives the unconstrained angular acceleration dOmega.
                // Maintain constraint B_b ∝ L_b: project dOmega onto L_b direction.
                // Since L_b.vz = 1, the scalar rotation rate is alpha = dOmega.vz.
                value_t const hw = width / 2.0;
                value_t const hh = height / 2.0;
                bivec3dp const L_b{0.0, 0.0, 1.0, hh, -hw, 0.0};
                auto const dOmega = compute_omega_dot(I_inv, F_b, Omega, I);
                value_t const alpha = dOmega.vz; // L_b.vz = 1

                rhs[0] = Omega;       // dB_b/dt = Omega_b  (= omega * L_b)
                rhs[1] = alpha * L_b; // dOmega_b/dt = alpha * L_b (pivot constraint)
            }

            void calc_rkstep(double dt)
            {
                auto u   = mdspan<bivec3dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto uh  = mdspan<bivec3dp, dextents<size_t, 2>>(uh_mem.data(), 2, 2);
                auto rhs = mdspan<bivec3dp const, dextents<size_t, 1>>(rhs_mem.data(), 2);

                for (size_t rk_step = 1; rk_step <= 4; ++rk_step) {
                    calc_rhs();
                    rk4_step(u, uh, rhs, dt, rk_step);
                }
            }

            void print_sim(double t)
            {
                auto u = mdspan<bivec3dp, dextents<size_t, 1>>(u_mem.data(), 2);

                fmt::println("t = {:>-7.3f}:", t);

                bivec3dp const B = u[0]; // B_b = phi * L_b

                // Current motor: M(t) = M0 ⟇ exp(½ B_b(t))
                auto const M = rgpr(M0, exp(0.5 * B));

                // phi = B.vz  (since B_b = phi * L_b and L_b.vz = 1)
                value_t const phi_b = B.vz;

                // CM world position
                vec3dp const O_b{0.0, 0.0, 0.0, 1.0};
                auto const cm_w = move3dp(O_b, M);

                fmt::println("    cm_w = {:>-6.3f}, phi_b = {:>-6.3f}, B = {:>-6.3f}",
                             cm_w, phi_b, B);
            }

            // Accessors for validation checks
            vec3dp get_cm_world() const
            {
                auto u = mdspan<bivec3dp const, dextents<size_t, 1>>(u_mem.data(), 2);
                auto const M = rgpr(M0, exp(0.5 * u[0]));
                return move3dp(vec3dp{0.0, 0.0, 0.0, 1.0}, M);
            }

            value_t get_omega() const
            {
                // Omega_b = omega * L_b, L_b.vz = 1  =>  omega = Omega_b.vz
                auto u = mdspan<bivec3dp const, dextents<size_t, 1>>(u_mem.data(), 2);
                return u[1].vz;
            }

            value_t get_I_zz_pivot() const { return I.view()[5, 2]; }

          private:

            value_t m;      // mass of cuboid [kg]
            value_t width;  // extent in e1/x direction in body frame [m]
            value_t height; // extent in e2/y direction in body frame [m]
            value_t depth;  // extent in e3/z direction in body frame [m]

            vec3dp cm_w_pos0;    // initial world position of cm
            value_t cm_w_phi0;   // initial rotation angle about pivot line
            value_t cm_w_omega0; // initial angular velocity about pivot line

            mvec3dp_e M0;        // translation motor: cm to initial world position
            vec3dp pivot_pt_w;   // world pos of body-frame pivot representative point

            // RK4 state: u[0] = B_b (position bivector), u[1] = Omega_b (velocity)
            std::vector<bivec3dp> u_mem;   // [B, Omega]
            std::vector<bivec3dp> uh_mem;  // RK4 helper (2 * SYS_SIZE entries)
            std::vector<bivec3dp> rhs_mem; // right-hand side [dB/dt, dOmega/dt]

            // Inertia map about pivot L_b (parallel-axis corrected via get_cuboid_inertia)
            inertia3dp I;
            inertia3dp I_inv;
        };

        value_t m = 1.0; // mass [kg]
        value_t w = 2.0; // width  [m] (e1/x in body frame)
        value_t h = 2.0; // height [m] (e2/y in body frame)
        value_t d = 1.0; // depth  [m] (e3/z in body frame)

        // Initial position: cm at (-1,-1,0) so pivot edge (at body (1,1,*)) is at
        // world origin (0,0,*). No initial rotation or angular velocity.
        auto cm_pos   = vec3dp{-1.0, -1.0, 0.0, 1.0};
        auto cm_phi   = 0.0; // no initial rotation
        auto cm_omega = 0.0; // no initial angular velocity

        sim_ode_cuboid_pga3dp sim(m, w, h, d, cm_pos, cm_phi, cm_omega);

        // time range: 0..2 s, 60 steps (analogous to 2D plate pendulum test)
        auto t_rng = discrete_range(0.0, 2.0, 60);

        sim.set_initial_values();
        sim.print_sim(t_rng.min());

        // --- Static validation: inertia analytical check -----------------------
        // I[5,2] = moment of inertia about z-axis (pivot line direction) at pivot.
        // Analytically (Steiner): I_cm_zz + m*(hw^2+hh^2)
        //   = m*(w^2+h^2)/12 + m*(hw^2+hh^2)
        //   = 1*(4+4)/12 + 1*(1+1) = 2/3 + 2 = 8/3
        value_t const hw = w / 2.0;
        value_t const hh = h / 2.0;
        CHECK(sim.get_I_zz_pivot() == doctest::Approx(8.0 / 3.0));

        // Cross-check: 3D cuboid I[5,2] matches 2D plate I[2,2] with same w,h,pivot.
        // (The z-extent d=1 does not enter the in-plane moment of inertia.)
        auto const I_plate = get_plate_inertia(m, w, h, vec2dp{hw, hh, 1.0});
        CHECK(I_plate.view()[2, 2] == doctest::Approx(sim.get_I_zz_pivot()));

        // --- Dynamic validation setup ------------------------------------------
        // Initial mechanical energy (t=0): kinetic=0, potential = m*g*cm_w.y
        value_t const g = 9.81;
        value_t const E_0 = m * g * cm_pos.y; // = -9.81 J (cm starts at y=-1)

        // Pivot is at world origin (0,0,0) => pivot distance squared = hw^2+hh^2 = 2.0
        value_t const pivot_dist_sq = hw * hw + hh * hh; // = 2.0

        for (size_t n = 1; n <= t_rng.steps(); n++) {
            double t = t_rng.min() + n * t_rng.delta();
            sim.calc_rkstep(t_rng.delta());
            sim.print_sim(t);

            vec3dp const cm_w = sim.get_cm_world();
            value_t const omega = sim.get_omega();

            // Energy conservation: T + V = E_0  (RK4 drift tolerance)
            value_t const T_kin = 0.5 * sim.get_I_zz_pivot() * omega * omega;
            value_t const V_pot = m * g * cm_w.y;
            CHECK((T_kin + V_pot) == doctest::Approx(E_0).epsilon(1e-3));

            // Pivot constraint: |cm_w|_xy^2 = hw^2 + hh^2 = 2.0
            value_t const r_sq = cm_w.x * cm_w.x + cm_w.y * cm_w.y;
            CHECK(r_sq == doctest::Approx(pivot_dist_sq).epsilon(1e-3));

            // Planar motion: z-component of cm stays zero
            CHECK(cm_w.z == doctest::Approx(0.0).epsilon(1e-10));
        }

        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: physics tests application")
