// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "doctest/doctest.h"

#include <cmath>
#include <iostream>
#include <memory>

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
        // BiVec3dp: (vx=e41, vy=e42, vz=e43, mx=e23, my=e31, mz=e12); the bulk slot mx
        // (zero weight) is a LINEAR velocity along e23 (a translation generator)
        BiVec3dp<double> Omega{0.0, 0.0, 0.0,
                               1.0, 0.0, 0.0}; // unit linear velocity (e23)

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

        auto const det_I_grid = hd::ga::det(I_grid.view());
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

    TEST_CASE("rk4_step: mdspan and std::vector overloads produce identical trajectories")
    {
        // Cross-check the two rk4_step overloads by integrating the same
        // problem with both and asserting bit-identical results at every
        // sub-step. Problem: 1-D harmonic oscillator du/dt = -u (closed
        // form: u(t) = u0 * exp(-t)), state encoded in a single vec3d.

        constexpr size_t n = 1;
        constexpr value_t dt = 0.01;
        constexpr size_t n_steps = 100;

        // mdspan side: backing storage + spans.
        std::vector<vec3d> u_md_mem{vec3d{1.0, 2.0, 3.0}};
        std::vector<vec3d> uh_md_mem(2 * n);
        std::vector<vec3d> rhs_md_mem(n);
        auto u_md = mdspan<vec3d, dextents<size_t, 1>>(u_md_mem.data(), n);
        auto uh_md = mdspan<vec3d, dextents<size_t, 2>>(uh_md_mem.data(), 2, n);
        auto rhs_md_mut = mdspan<vec3d, dextents<size_t, 1>>(rhs_md_mem.data(), n);
        auto rhs_md = mdspan<vec3d const, dextents<size_t, 1>>(rhs_md_mem.data(), n);

        // vector side: same initial state as the mdspan side.
        std::vector<vec3d> u_vec{vec3d{1.0, 2.0, 3.0}};
        std::array<std::vector<vec3d>, 2> uh_vec{std::vector<vec3d>(n),
                                                 std::vector<vec3d>(n)};
        std::vector<vec3d> rhs_vec(n);

        for (size_t step = 0; step < n_steps; ++step) {
            for (size_t rk = 1; rk <= 4; ++rk) {
                // rhs(u) = -u for both sides.
                rhs_md_mut[0] = vec3d{-u_md[0].x, -u_md[0].y, -u_md[0].z};
                rhs_vec[0] = vec3d{-u_vec[0].x, -u_vec[0].y, -u_vec[0].z};

                // Sanity: rhs values match before the rk4_step call.
                CHECK(rhs_md[0].x == doctest::Approx(rhs_vec[0].x));
                CHECK(rhs_md[0].y == doctest::Approx(rhs_vec[0].y));
                CHECK(rhs_md[0].z == doctest::Approx(rhs_vec[0].z));

                rk4_step(u_md, uh_md, rhs_md, dt, rk);
                std::tie(u_vec, uh_vec) =
                    rk4_step(std::move(u_vec), std::move(uh_vec), rhs_vec, dt, rk);

                // Both overloads must agree to the bit on every sub-step.
                CHECK(u_md[0].x == u_vec[0].x);
                CHECK(u_md[0].y == u_vec[0].y);
                CHECK(u_md[0].z == u_vec[0].z);
                CHECK(uh_md[0, 0].x == uh_vec[0][0].x);
                CHECK(uh_md[1, 0].x == uh_vec[1][0].x);
            }
        }

        // Sanity: the trajectory is approximately u0 * exp(-T) at the end.
        value_t const T = n_steps * dt;
        value_t const decay = std::exp(-T);
        CHECK(u_vec[0].x == doctest::Approx(1.0 * decay).epsilon(1e-6));
        CHECK(u_vec[0].y == doctest::Approx(2.0 * decay).epsilon(1e-6));
        CHECK(u_vec[0].z == doctest::Approx(3.0 * decay).epsilon(1e-6));
    }

    TEST_CASE("rk4_step (vector overload): input-validation errors")
    {
        std::vector<vec3d> u{vec3d{1, 2, 3}};
        std::array<std::vector<vec3d>, 2> uh{std::vector<vec3d>(1),
                                             std::vector<vec3d>(1)};
        std::vector<vec3d> const rhs{vec3d{0, 0, 0}};

        // rk_step out of [1,4] must throw.
        CHECK_THROWS_AS(rk4_step(u, uh, rhs, 0.01, 0), std::invalid_argument);
        CHECK_THROWS_AS(rk4_step(u, uh, rhs, 0.01, 5), std::invalid_argument);

        // Length mismatch must throw.
        std::vector<vec3d> const rhs_bad{vec3d{0, 0, 0}, vec3d{0, 0, 0}};
        CHECK_THROWS_AS(rk4_step(u, uh, rhs_bad, 0.01, 1), std::invalid_argument);
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
            // - M0: pure translation motor placing body origin (= cm) at initial world
            // pos
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
                                  vec3dp const& cm_pos_in, value_t cm_phi_in,
                                  value_t cm_omega_in) :
                m(m_in), width(w_in), height(h_in), depth(d_in), cm_w_pos0(cm_pos_in),
                cm_w_phi0(cm_phi_in), cm_w_omega0(cm_omega_in), u_mem(2), uh_mem(2 * 2),
                rhs_mem(2)
            {
                fmt::println(
                    "sim_ode_cuboid_pga3dp: 3D pendulum (pivot = z-edge at TR).");

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
                // In PGA3DP, translate by (cx,cy,cz): motor = exp(0.5 *
                // bivec3dp{0,0,0,cx,cy,cz}) Translation is encoded in the bulk
                // (e23,e31,e12) components, not the weight (e41,e42,e43) components.
                M0 = exp(0.5 *
                         bivec3dp{0.0, 0.0, 0.0, cm_w_pos0.x, cm_w_pos0.y, cm_w_pos0.z});

                // pivot_pt_w: world position of representative pivot point (hw,hh,0,1)
                // in body frame, mapped via M0.  Stays fixed during rotation since
                // move3dp(P_on_L_b, M0 ⟇ exp(½ phi * L_b)) = M0-image of P_on_L_b.
                pivot_pt_w = move3dp(vec3dp{hw, hh, 0.0, 1.0}, M0);
                fmt::println(
                    "pivot_pt_w = {:>-7.3f}  (world position of (hw,hh,0) via M0)",
                    pivot_pt_w);

                // B_b(0): initial rotation by phi0 about body-frame pivot line L_b.
                bivec3dp const L_b{0.0, 0.0, 1.0, hh, -hw, 0.0};
                u[0] = cm_w_phi0 * L_b;

                // Omega_b(0): initial angular velocity about pivot line L_b.
                u[1] = cm_w_omega0 * L_b;
            }

            void calc_rhs()
            {
                auto u = mdspan<bivec3dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto rhs = mdspan<bivec3dp, dextents<size_t, 1>>(rhs_mem.data(), 2);

                bivec3dp const B = u[0];     // position bivector B_b = phi * L_b
                bivec3dp const Omega = u[1]; // velocity bivector Omega_b = omega * L_b

                // Current motor: M(t) = M0 ⟇ exp(½ B_b(t))  [body-frame formulation]
                auto const M = rgpr(M0, exp(0.5 * B));

                // CM world position: body origin = (0,0,0,1) in body frame
                vec3dp const O_b{0.0, 0.0, 0.0, 1.0};
                auto const cm_w = move3dp(O_b, M);

                // Force couple in world frame (net force = 0, only torque):
                // - Gravity at cm (downward -y):       F_dn = wdg(cm_w, (0,-mg,0,0))
                // - Reaction at pivot (upward +y):     F_up = wdg(pivot_pt_w,
                // (0,+mg,0,0)) Using a single representative pivot POINT (z=0 on pivot
                // edge) here. The torque from the reaction force about L_b vanishes by
                // definition.
                auto const F_dn_w = wdg(cm_w, vec3dp{0.0, -m * 9.81, 0.0, 0.0});
                auto const F_up_w = wdg(pivot_pt_w, vec3dp{0.0, m * 9.81, 0.0, 0.0});
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
                auto u = mdspan<bivec3dp, dextents<size_t, 1>>(u_mem.data(), 2);
                auto uh = mdspan<bivec3dp, dextents<size_t, 2>>(uh_mem.data(), 2, 2);
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

            mvec3dp_e M0;      // translation motor: cm to initial world position
            vec3dp pivot_pt_w; // world pos of body-frame pivot representative point

            // RK4 state: u[0] = B_b (position bivector), u[1] = Omega_b (velocity)
            std::vector<bivec3dp> u_mem;   // [B, Omega]
            std::vector<bivec3dp> uh_mem;  // RK4 helper (2 * SYS_SIZE entries)
            std::vector<bivec3dp> rhs_mem; // right-hand side [dB/dt, dOmega/dt]

            // Inertia map about pivot L_b (parallel-axis corrected via
            // get_cuboid_inertia)
            inertia3dp I;
            inertia3dp I_inv;
        };

        value_t m = 1.0; // mass [kg]
        value_t w = 2.0; // width  [m] (e1/x in body frame)
        value_t h = 2.0; // height [m] (e2/y in body frame)
        value_t d = 1.0; // depth  [m] (e3/z in body frame)

        // Initial position: cm at (-1,-1,0) so pivot edge (at body (1,1,*)) is at
        // world origin (0,0,*). No initial rotation or angular velocity.
        auto cm_pos = vec3dp{-1.0, -1.0, 0.0, 1.0};
        auto cm_phi = 0.0;   // no initial rotation
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


/////////////////////////////////////////////////////////////////////////////////////////
// PGA3DP dynamic_system3dp - Milestone 1: free rigid bodies
//   (the first 3D deliverable of the 2D->3D lift; validates the ported system classes,
//    the pose3dp constrained motor log, and the se(3) Euler equation against energy +
//    angular-momentum conservation and the classic intermediate-axis instability)
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: dynamic_system3dp (M1)")
{

    TEST_CASE("pga3dp: free rigid body energy + cm parabola (M1)")
    {
        fmt::println(
            "pga3dp: dynamic_system3dp - free rigid body energy + parabola (M1)");

        // a uniform cube (symmetric inertia) launched under gravity with an initial
        // linear velocity and a spin about e3. Symmetric => omega is constant; gravity
        // exerts no torque about the cm => the cm follows an exact parabola and energy is
        // conserved.
        value_t const m = 2.0;
        auto const cube =
            make_cuboid_body(m, 1.0, 1.0, 1.0); // I_cm = m/6 on the diagonal

        dynamic_system3dp sys;
        sys.add_frame(static_frame3dp("W")); // inertial root at the world origin

        vec3dp const v0{1.0, 2.0, 0.0, 0.0}; // initial cm velocity
        vec3dp const w0{0.0, 0.0, 3.0, 0.0}; // initial spin about e3
        sys.add_body(static_frame3dp("B", vec3dp{0.0, 0.0, 0.0, 1.0}), cube,
                     kin_state3dp{.vel = v0, .omega = w0});

        value_t const g = 9.81; // default gravity is (0,-g,0)
        value_t const E0 = sys.total_energy();
        // analytic: KE = 1/2 m|v|^2 + 1/2 I_zz w^2; PE = 0 (cm at y = 0)
        value_t const Izz = m / 6.0; // cube: m*(1+1)/12
        value_t const KE_analytic = 0.5 * m * (1.0 + 4.0) + 0.5 * Izz * 9.0;
        CHECK(E0 == doctest::Approx(KE_analytic));

        value_t const dt = 1.0e-3;
        size_t const N = 1500; // 1.5 s
        value_t max_dE = 0.0, max_cm_err = 0.0, max_w_err = 0.0;
        for (size_t n = 1; n <= N; ++n) {
            sys.step(dt);
            value_t const t = n * dt;

            // energy conservation (RK4 drift)
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));

            // cm follows the exact parabola  cm(t) = v0*t + 1/2 g t^2 (absolute error:
            // the body-frame twist representation couples the spin into the RK4
            // truncation)
            vec3dp const cm = move3dp(O_3dp, sys.get_pos_trafo("B", "W"));
            max_cm_err =
                std::max({max_cm_err, std::abs(cm.x - 1.0 * t),
                          std::abs(cm.y - (2.0 * t - 0.5 * g * t * t)), std::abs(cm.z)});

            // symmetric body, gravity has no torque about the cm => the body spin is
            // EXACTLY constant (omega_dot = 0, so RK4 leaves it untouched)
            twist3dp const Om = sys.relative_twist("B"); // omega = (vx,vy,vz)
            max_w_err = std::max(
                {max_w_err, std::abs(Om.vx), std::abs(Om.vy), std::abs(Om.vz - 3.0)});
        }
        fmt::println("  E0 = {:.6f}, dE/E = {:.2e}, max cm err = {:.2e}, max omega err = "
                     "{:.2e}",
                     E0, max_dE / std::abs(E0), max_cm_err, max_w_err);
        CHECK(max_dE / std::abs(E0) < 1e-4); // RK4 drift (spin-coupled)
        CHECK(max_cm_err < 1e-3);
        CHECK(max_w_err < 1e-12); // exact: omega_dot == 0 for this case
        fmt::println("");
    }

    TEST_CASE("pga3dp: torque-free intermediate-axis instability (Dzhanibekov, M1)")
    {
        fmt::println("pga3dp: dynamic_system3dp - torque-free Dzhanibekov flip (M1)");

        // an asymmetric cuboid (w<h<d) spun about its INTERMEDIATE principal axis (e2):
        //   I_xx = m(h^2+d^2)/12 = 13/12,  I_yy = m(w^2+d^2)/12 = 10/12 (intermediate),
        //   I_zz = m(w^2+h^2)/12 =  5/12.
        // The tennis-racket / Dzhanibekov theorem: spin about the intermediate axis is
        // unstable, so a small perturbation makes the spin axis periodically FLIP. With
        // no gravity, energy AND the world angular-momentum vector are exactly conserved
        // -- the unambiguous correctness test of the se(3) Euler equation.
        value_t const m = 1.0;
        auto const body = make_cuboid_body(m, 1.0, 2.0, 3.0);

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0}); // torque-free
        sys.add_frame(static_frame3dp("W"));

        // spin mostly about e2 (intermediate) with a small e1 perturbation; no
        // translation
        vec3dp const w0{0.3, 5.0, 0.0, 0.0};
        sys.add_body(static_frame3dp("B"), body, kin_state3dp{.omega = w0});

        value_t const E0 = sys.total_energy(); // kinetic only
        bivec3dp const L0 = sys.momentum_world(1);
        value_t const L0n = std::sqrt(L0.mx * L0.mx + L0.my * L0.my +
                                      L0.mz * L0.mz); // |angular momentum|

        value_t const dt = 1.0e-3;
        size_t const N = 12000; // 12 s -- long enough for at least one flip
        value_t max_dE = 0.0, max_dL = 0.0, min_wy = w0.y, max_abs_wx = 0.0, max_cm = 0.0;
        for (size_t n = 1; n <= N; ++n) {
            sys.step(dt);
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));

            bivec3dp const L = sys.momentum_world(1);
            max_dL = std::max({max_dL, std::abs(L.mx - L0.mx), std::abs(L.my - L0.my),
                               std::abs(L.mz - L0.mz)});

            twist3dp const Om = sys.relative_twist("B"); // body-frame omega = (vx,vy,vz)
            min_wy = std::min(min_wy, Om.vy);
            max_abs_wx = std::max(max_abs_wx, std::abs(Om.vx));

            // cm stays put (no linear momentum)
            vec3dp const cm = move3dp(O_3dp, sys.get_pos_trafo(1, 0));
            max_cm = std::max({max_cm, std::abs(cm.x), std::abs(cm.y), std::abs(cm.z)});
        }
        fmt::println("  E0 = {:.6f}, dE/E = {:.2e}, |L|0 = {:.4f}, max|dL|/|L| = {:.2e}",
                     E0, max_dE / std::abs(E0), L0n, max_dL / L0n);
        fmt::println(
            "  spin about e2: w_y in [{:.3f}, {:.3f}], max|w_x| = {:.3f}, max cm = "
            "{:.2e}",
            min_wy, w0.y, max_abs_wx, max_cm);

        // conservation: energy and the world angular-momentum vector
        CHECK(max_dE / std::abs(E0) < 1e-6);
        CHECK(max_dL / L0n < 1e-5);
        CHECK(max_cm < 1e-8); // cm fixed (pure rotation, no linear momentum)
        // the flip: w_y reverses sign (drops from +5 to clearly negative), and the
        // perturbation about e1 grows large during the flip
        CHECK(min_wy < -3.0);
        CHECK(max_abs_wx > 2.0);
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: dynamic_system3dp (M1)")


/////////////////////////////////////////////////////////////////////////////////////////
// PGA3DP dynamic_system3dp - Milestone 2: single revolute joint (spatial pendulum)
//   validates the ported joint tier + inertia-map forward_dynamics against the analytic
//   compound-pendulum reduced inertia, restoring law, small-oscillation frequency, and
//   energy conservation.
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: dynamic_system3dp (M2)")
{

    TEST_CASE("pga3dp: revolute joint - compound spatial pendulum (M2)")
    {
        fmt::println("pga3dp: dynamic_system3dp - spatial pendulum (revolute joint, M2)");

        // A cube hinged at a fixed horizontal axis (e3), its cm a distance L below the
        // axis, swinging under gravity -> a compound pendulum. The hinge moment of
        // inertia emerges (Steiner) from the spatial Jacobian: I_hinge = I_cm_zz + m L^2
        // = m s^2/6 + m L^2. Small-oscillation frequency w0 = sqrt(m g L / I_hinge).
        value_t const m = 1.0, s = 1.0, L = 1.0, g = 9.81;
        auto const cube = make_cuboid_body(m, s, s, s);

        value_t const I_hinge = m * s * s / 6.0 + m * L * L;         // = 7/6
        value_t const w0 = std::sqrt(g * L / (s * s / 6.0 + L * L)); // rad/s
        value_t const T0 = 2.0 * std::acos(-1.0) / w0;               // period

        // build a pendulum at initial angle phi0: rest cm at (0,-L,0), hinge line through
        // the body point (0,L,0) along e3 (so the world hinge sits at the origin)
        auto build = [&](value_t phi0) {
            auto sys = std::make_unique<dynamic_system3dp>();
            sys->add_frame(static_frame3dp("W"));
            sys->add_revolute_body(static_frame3dp("B", vec3dp{0.0, -L, 0.0, 1.0}), cube,
                                   vec3dp{0.0, L, 0.0, 1.0}, vec3dp{0.0, 0.0, 1.0, 0.0},
                                   phi0);
            return sys;
        };

        // --- a) reduced inertia + the pendulum restoring law at a small angle ----------
        {
            value_t const phi0 = 0.05;
            auto sys = build(phi0);
            value_t const alpha0 = -(g * L / (s * s / 6.0 + L * L)) * std::sin(phi0);
            fmt::println(
                "  I_hinge = {:.5f} (analytic {:.5f}); alpha0 = {:.5f} (analytic "
                "{:.5f})",
                sys->mass_matrix()[0], I_hinge, sys->joint_accel(1), alpha0);
            CHECK(sys->mass_matrix()[0] == doctest::Approx(I_hinge)); // reduced inertia
            CHECK(sys->joint_accel(1) == doctest::Approx(alpha0));    // restoring law
        }

        // --- b) small-oscillation period via the first phi zero-crossing (T0/4) --------
        {
            auto sys = build(0.05);
            value_t const dt = 1.0e-4;
            value_t t_cross = 0.0, phi_prev = sys->joint_phi(1);
            for (size_t n = 1; n <= 20000; ++n) {
                sys->step(dt);
                value_t const phi = sys->joint_phi(1);
                if (phi_prev > 0.0 && phi <= 0.0) { // first downward zero crossing
                    t_cross = ((n - 1) + phi_prev / (phi_prev - phi)) * dt; // interp
                    break;
                }
                phi_prev = phi;
            }
            value_t const T_meas = 4.0 * t_cross;
            fmt::println("  T0 = {:.5f} s (w0 = {:.5f}); T_measured = {:.5f} s", T0, w0,
                         T_meas);
            CHECK(T_meas == doctest::Approx(T0).epsilon(2e-3)); // anharmonic + RK4
        }

        // --- c) energy conservation + M-vs-KE consistency over a moderate swing --------
        {
            auto sys = build(0.6); // moderate amplitude (non-chaotic single pendulum)
            value_t const E0 = sys->total_energy();
            value_t const dt = 5.0e-4;
            value_t max_dE = 0.0, kemax = 0.0, max_mke = 0.0;
            for (size_t n = 1; n <= 6000; ++n) { // 3 s
                sys->step(dt);
                max_dE = std::max(max_dE, std::abs(sys->total_energy() - E0));
                value_t const ke = sys->kinetic_energy();
                kemax = std::max(kemax, ke);
                // 1/2 qdot^T M(q) qdot == kinetic_energy()
                value_t const w = sys->joint_omega(1);
                max_mke =
                    std::max(max_mke, std::abs(0.5 * sys->mass_matrix()[0] * w * w - ke));
            }
            fmt::println("  E0 = {:.5f}, KEmax = {:.5f}, dE/KEmax = {:.2e}, M-vs-KE = "
                         "{:.2e}",
                         E0, kemax, max_dE / kemax, max_mke);
            CHECK(max_dE / kemax < 1e-6); // RK4 energy drift
            CHECK(max_mke < 1e-9);        // mass-matrix / kinetic-energy identity
        }

        fmt::println("");
    }

    TEST_CASE("pga3dp: prismatic joint - slider on an incline (M2)")
    {
        fmt::println("pga3dp: dynamic_system3dp - prismatic slider (M2)");

        // A cube on a prismatic joint sliding along a fixed 45-deg direction in the x-y
        // plane under gravity. A slider has NO angular DOF, so the reduced inertia is
        // just the mass (M[0] = m) and the acceleration is the gravity component along
        // the slide direction: q-ddot = dir . g. This exercises the rotation/translation
        // unification: the SAME inertia-map forward_dynamics handles it with no angular
        // term -- the 3D analogue of the 2D prismatic check that caught the dropped
        // S.z*S.z term (a trap that no longer exists in the inertia-map form).
        value_t const m = 2.0, g = 9.81;
        auto const cube = make_cuboid_body(m, 1.0, 1.0, 1.0);
        value_t const c = 1.0 / std::sqrt(2.0);
        vec3dp const dir{c, -c, 0.0, 0.0}; // unit 45-deg direction (down-right)

        dynamic_system3dp sys;
        sys.add_frame(static_frame3dp("W"));
        sys.add_prismatic_body(static_frame3dp("B"), cube, dir);

        value_t const a_analytic = dir.y * (-g); // dir . g  (g = (0,-g,0)) = c*g
        fmt::println("  M[0] = {:.5f} (analytic {:.5f}); a_slide = {:.5f} (analytic "
                     "{:.5f})",
                     sys.mass_matrix()[0], m, sys.joint_accel(1), a_analytic);
        CHECK(sys.mass_matrix()[0] == doctest::Approx(m)); // pure mass, NO angular term
        CHECK(sys.joint_accel(1) == doctest::Approx(a_analytic)); // q-ddot = dir . g

        // energy conservation over the slide (constant accel => RK4 near-exact)
        value_t const E0 = sys.total_energy();
        value_t const dt = 5.0e-4;
        value_t max_dE = 0.0, kemax = 0.0;
        for (size_t n = 1; n <= 4000; ++n) { // 2 s
            sys.step(dt);
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));
            kemax = std::max(kemax, sys.kinetic_energy());
            // the slide stays straight: the cm only translates along dir (no rotation)
            CHECK(sys.joint_accel(1) == doctest::Approx(a_analytic)); // a stays constant
        }
        fmt::println("  E0 = {:.5f}, KEmax = {:.5f}, dE/KEmax = {:.2e}", E0, kemax,
                     max_dE / kemax);
        CHECK(max_dE / kemax < 1e-10);

        fmt::println("");
    }

    TEST_CASE("pga3dp: prismatic joint - spring/damper damped oscillator (Phase A)")
    {
        fmt::println("pga3dp: dynamic_system3dp - joint spring/damper (Phase A)");

        // A cube on a horizontal prismatic joint carrying a linear spring (k) + damper
        // (c) on its coordinate, no gravity: the reduced inertia is the pure mass (M[0] =
        // m), so the joint obeys  m q'' + c q' + k q = 0  -- a damped harmonic oscillator
        // with a closed-form solution. Validates the additive spring/damper force
        // elements (set_joint_spring_damper) folded into tau.
        value_t const m = 2.0, k = 8.0, c = 0.8;
        value_t const wn = std::sqrt(k / m);                  // natural frequency 2 rad/s
        value_t const zeta = c / (2.0 * std::sqrt(k * m));    // damping ratio 0.1
        value_t const wd = wn * std::sqrt(1.0 - zeta * zeta); // damped frequency
        value_t const A0 = 0.5;                               // initial q, zero velocity

        auto const cube = make_cuboid_body(m, 1.0, 1.0, 1.0);

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0}); // isolate the spring oscillator
        sys.add_frame(static_frame3dp("W"));
        sys.add_prismatic_body(static_frame3dp("B"), cube, vec3dp{1.0, 0.0, 0.0, 0.0}, A0,
                               0.0);          // slide along e1; q0 = A0, v0 = 0
        sys.set_joint_spring_damper(1, k, c); // q_rest = 0

        // reduced inertia is the pure mass; initial accel = -k*A0/m (spring only, v = 0)
        CHECK(sys.mass_matrix()[0] == doctest::Approx(m));
        CHECK(sys.joint_accel(1) == doctest::Approx(-k * A0 / m));

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
            max_err = std::max(max_err, std::abs(sys.joint_phi(1) - q_exact(t)));
        }
        fmt::println(
            "  wn = {:.4f}, zeta = {:.3f}, wd = {:.4f}; max|q - q_exact| = {:.2e}", wn,
            zeta, wd, max_err);
        CHECK(max_err < 1e-9); // RK4 matches the closed form to roundoff (~1e-12)

        fmt::println("");
    }

    TEST_CASE("pga3dp: prismatic joint - undamped spring energy conservation (Phase A)")
    {
        fmt::println("pga3dp: dynamic_system3dp - undamped spring energy (Phase A)");

        // Same rig with c = 0: KE + spring potential 1/2 k q^2 must be conserved, which
        // exercises the spring term added to potential_energy().
        value_t const m = 2.0, k = 8.0, A0 = 0.5;
        auto const cube = make_cuboid_body(m, 1.0, 1.0, 1.0);

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
        sys.add_frame(static_frame3dp("W"));
        sys.add_prismatic_body(static_frame3dp("B"), cube, vec3dp{1.0, 0.0, 0.0, 0.0}, A0,
                               0.0);
        sys.set_joint_spring_damper(1, k, 0.0); // spring only, no damping

        value_t const E0 = sys.total_energy(); // = 1/2 k A0^2
        CHECK(E0 == doctest::Approx(0.5 * k * A0 * A0));

        value_t const dt = 1.0e-4;
        value_t max_dE = 0.0, kemax = 0.0;
        for (size_t nstep = 1; nstep <= 40000; ++nstep) { // 4 s
            sys.step(dt);
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));
            kemax = std::max(kemax, sys.kinetic_energy());
        }
        fmt::println("  E0 = {:.5f}, KEmax = {:.5f}, dE/KEmax = {:.2e}", E0, kemax,
                     max_dE / kemax);
        CHECK(max_dE / kemax < 1e-11); // KE + spring PE conserved to roundoff

        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: dynamic_system3dp (M2)")


/////////////////////////////////////////////////////////////////////////////////////////
// PGA3DP dynamic_system3dp - Milestone 3: coupled spatial double pendulum
//   two revolute joints with NON-PARALLEL axes -> genuine 3D coupling. This is where the
//   spatial-form forward_dynamics earns its keep: the angular bias-acceleration AND the
//   se(3) gyroscopic generalised force (both zero for a single joint) are now non-zero
//   and required. Validated by energy conservation over the chaotic swing + the
//   mass-matrix / kinetic-energy identity.
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: dynamic_system3dp (M3)")
{

    TEST_CASE("pga3dp: spatial double pendulum - non-parallel axes (M3)")
    {
        fmt::println("pga3dp: dynamic_system3dp - spatial double pendulum (M3)");

        // two cubes, two revolute joints with NON-PARALLEL axes: B1 hinges about e3, B2
        // (child of B1) hinges about e1. The axes are non-parallel and time-varying as B1
        // turns, so B2 swings out of B1's plane -> genuine 3D coupling (the gyroscopic +
        // angular-bias terms are exercised). No analytic trajectory (chaotic); validate
        // by energy conservation and the identity 1/2 qdot^T M(q) qdot ==
        // kinetic_energy().
        value_t const m = 1.0, sgeo = 0.6;
        auto const cube = make_cuboid_body(m, sgeo, sgeo, sgeo);

        dynamic_system3dp sys;
        sys.add_frame(static_frame3dp("W"));
        // B1: hinge about e3 at the world origin, cm 1 below; given a brisk initial spin
        // so it drags B2 around -> the gyroscopic / out-of-plane coupling is strongly
        // excited
        sys.add_revolute_body(static_frame3dp("B1", vec3dp{0.0, -1.0, 0.0, 1.0}), cube,
                              vec3dp{0.0, 1.0, 0.0, 1.0}, vec3dp{0.0, 0.0, 1.0, 0.0}, 1.2,
                              3.0);
        // B2: child of B1, hinge about e1 (NON-parallel to e3), cm 1 below B1's cm
        sys.add_revolute_body(static_frame3dp("B2", vec3dp{0.0, -1.0, 0.0, 1.0}), cube,
                              vec3dp{0.0, 1.0, 0.0, 1.0}, vec3dp{1.0, 0.0, 0.0, 0.0}, 0.6,
                              0.0, 1 /* parent = B1 */);

        value_t const E0 = sys.total_energy();
        value_t const dt = 2.0e-4;
        size_t const N = 25000; // 5 s
        value_t max_dE = 0.0, kemax = 0.0, max_mke = 0.0, min_phi2 = 0.4, max_phi2 = 0.4;
        for (size_t n = 1; n <= N; ++n) {
            sys.step(dt);
            value_t const ke = sys.kinetic_energy();
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));
            kemax = std::max(kemax, ke);
            // 1/2 qdot^T M(q) qdot == kinetic_energy()  (M is symmetric: M[1] == M[2])
            auto const M = sys.mass_matrix();
            value_t const w1 = sys.joint_omega(1), w2 = sys.joint_omega(2);
            value_t const qMq =
                0.5 * (M[0] * w1 * w1 + (M[1] + M[2]) * w1 * w2 + M[3] * w2 * w2);
            max_mke = std::max(max_mke, std::abs(qMq - ke));
            min_phi2 = std::min(min_phi2, sys.joint_phi(2));
            max_phi2 = std::max(max_phi2, sys.joint_phi(2));
        }
        fmt::println("  E0 = {:.5f}, KEmax = {:.5f}, dE/KEmax = {:.2e}, M-vs-KE = {:.2e}",
                     E0, kemax, max_dE / kemax, max_mke);
        fmt::println("  phi2 in [{:.3f}, {:.3f}] (swings out of plane)", min_phi2,
                     max_phi2);

        CHECK(max_dE / kemax < 1e-6); // energy conserved over the chaotic 3D swing
        CHECK(max_mke < 1e-9);        // mass-matrix / kinetic-energy identity
        // genuine 3D motion: the second joint actually swings appreciably
        CHECK(max_phi2 - min_phi2 > 0.3);

        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: dynamic_system3dp (M3)")


/////////////////////////////////////////////////////////////////////////////////////////
// closed_loop_system3dp -- the 3D lift of the closed-loop layer, exercised on a
// genuinely-spatial closed mechanism: a TWO-ARM REACHER whose two hands are pinned
// together.
//
// Geometry (top view, looking down the vertical e3 axis; gravity acts along -e2):
//
//                              tip (shared hand)
//                                    o                 e2 (up)
//                          forearm  / \  forearm        ^
//                                  /   \                 |
//                          elbow  E1    E2  elbow        +---> e1
//                                /       \
//                       upperarm/         \upperarm
//                              /           \
//                        S1  o             o  S2     <- shoulders (fixed to the world)
//                       (-1.2,0,0)      (1.2,0,0)
//
// Two identical 2R arms hang from fixed shoulders on the e1 axis (1.2 to each side of the
// origin). Each arm has two revolute joints and two unit-length links:
//
//   shoulder (S1/S2): YAW   about the vertical e3 axis  -> swings the arm in the e1-e2
//   plane elbow    (E1/E2): PITCH about the e2 axis           -> lifts the forearm OUT of
//   that
//                                                          plane (into +/- e3)
//
// The loop is closed by a single point-coincidence constraint: the two forearm tips (a
// point L along each elbow frame's local e1) must occupy the same world point -- a
// spatial spherical ("ball") joint joining the two hands. That is 3 scalar equations (x,
// y, z).
//
// DOF: 4 revolute joints - 3 coincidence equations = 1 DOF. The redundancy that IS the
// one remaining freedom lives in the two elbows: the two shoulders (yaw) set where the
// shared hand sits in the horizontal e1-e2 plane (2 dof for the 2 in-plane equations),
// and ONE elbow pitch is enough to satisfy the vertical (e3) equation -- so the OTHER
// elbow pitch is free, and folding one elbow forces the other to follow. Driving an elbow
// (E1) and solving {S1, S2, E2} therefore gives a well-posed, full-rank assembly; driving
// a shoulder instead would leave two elbows fighting over the single vertical equation
// while one shoulder alone cannot span the horizontal plane -- a rank-deficient
// (singular) Jacobian.
//
// The out-of-plane elbow pitch is also what keeps the *constraint* full-rank in 3D: a
// flat (planar) linkage would make the vertical equation permanently inactive, leaving
// the KKT saddle-point matrix singular. Here all three coincidence components are live.
//
// Validated across all three tiers on the pga3dp base (sharing the dimension-agnostic
// lstsq_solve / kkt_solve kernels with the 2D layer): Phase 1 assembly (drive g -> 0),
// Phase 2 consistent velocity distribution (the shared hand has one well-defined
// velocity), and Phase 3 energy-conserving constrained dynamics.
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA3DP: closed_loop_system3dp")
{

    TEST_CASE("pga3dp: closed_loop_system3dp - spatial two-arm loop (assembly + energy)")
    {
        fmt::println("pga3dp: closed_loop_system3dp - spatial two-arm loop");

        // one unit-mass cube body reused for every link; L = link length (upper arm =
        // forearm = 1). A frame's rest origin is its centre of mass, so the energy
        // bookkeeping (potential = -m g . cm) reads the frame origins directly.
        value_t const m = 1.0, sgeo = 0.4, L = 1.0;
        auto const link = make_cuboid_body(m, sgeo, sgeo, sgeo);

        closed_loop_system3dp cl;
        cl.add_frame(static_frame3dp("W")); // inertial root (world)

        // --- arm 1 (left), hanging from the fixed shoulder S1 at (-1.2, 0, 0)
        // ----------- shoulder S1: yaw about e3 (pivot = its own origin, axis = e3);
        // start yawed inward (0.8 rad) so the hand reaches toward the centre
        cl.add_revolute_body(static_frame3dp("S1", vec3dp{-1.2, 0.0, 0.0, 1.0}), link,
                             vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{0.0, 0.0, 1.0, 0.0},
                             /*yaw0*/ 0.8, /*rate0*/ 0.0, cl.index_of("W"));
        // elbow E1: child of S1, its rest origin one link-length (L) out along S1's e1;
        // pitch about e2 (lifts the forearm out of the horizontal plane)
        cl.add_revolute_body(static_frame3dp("E1", vec3dp{L, 0.0, 0.0, 1.0}), link,
                             vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{0.0, 1.0, 0.0, 0.0},
                             /*pitch0*/ 0.0, /*rate0*/ 0.0, cl.index_of("S1"));

        // --- arm 2 (right), mirror image, shoulder S2 at (1.2, 0, 0)
        // --------------------- shoulder S2 yawed inward (2.34 rad ~ pi - 0.8) so its
        // hand also reaches the centre
        cl.add_revolute_body(static_frame3dp("S2", vec3dp{1.2, 0.0, 0.0, 1.0}), link,
                             vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{0.0, 0.0, 1.0, 0.0},
                             /*yaw0*/ 2.34, /*rate0*/ 0.0, cl.index_of("W"));
        cl.add_revolute_body(static_frame3dp("E2", vec3dp{L, 0.0, 0.0, 1.0}), link,
                             vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{0.0, 1.0, 0.0, 0.0},
                             /*pitch0*/ 0.0, /*rate0*/ 0.0, cl.index_of("S2"));

        size_t const E1 = cl.index_of("E1"), E2 = cl.index_of("E2");

        // close the loop: arm 1's hand (a point L along E1's local e1) must coincide with
        // arm 2's hand (L along E2's local e1) -- a spatial spherical joint joining them
        cl.add_loop_constraint(loop_constraint3dp{E1, vec3dp{L, 0.0, 0.0, 1.0}, E2,
                                                  vec3dp{L, 0.0, 0.0, 1.0},
                                                  constraint3dp::coincidence});

        // 1. (Phase 1) the two hands start apart (loop open); assemble() folds the joints
        //    until they meet. Drive an ELBOW (the 1-DOF freedom) and solve {S1, S2, E2}:
        //    that dependent set spans the 3D hand-gap, so the Newton Jacobian is full
        //    rank.
        value_t const g_open = cl.residual_norm();
        CHECK(g_open > 0.05);
        cl.assemble(/*driven*/ {E1});
        CHECK(cl.residual_norm() < 1e-10);

        // 2. (Phase 2) give the mechanism a consistent initial motion: spin the driver
        //    elbow at 1.5 rad/s and distribute that to the other joints so G q-dot = 0.
        cl.set_joint_rate(E1, 1.5);
        cl.solve_velocities(/*driven*/ {E1});
        // velocity closure: with the loop respected, the single shared hand has ONE
        // velocity -- arm 1 and arm 2 must compute the same value for it.
        auto unit = [](vec3dp const& p) { return unitize(p); };
        auto const PtA =
            unit(move3dp(vec3dp{L, 0.0, 0.0, 1.0}, cl.system().get_pos_trafo(E1, 0)));
        auto const PtB =
            unit(move3dp(vec3dp{L, 0.0, 0.0, 1.0}, cl.system().get_pos_trafo(E2, 0)));
        auto const vA = cl.system().point_velocity(PtA, E1);
        auto const vB = cl.system().point_velocity(PtB, E2);
        CHECK(std::abs(vA.x - vB.x) < 1e-9);
        CHECK(std::abs(vA.y - vB.y) < 1e-9);
        CHECK(std::abs(vA.z - vB.z) < 1e-9);

        // 3. (Phase 3) let it swing: constrained forward dynamics under gravity (the KKT
        //    solve supplies the joint accelerations + constraint/leg forces; RK4
        //    advances; projection removes drift). The initial spin plus gravity exchange
        //    kinetic and potential energy, but the constraint forces do no work, so TOTAL
        //    energy is conserved to integrator tolerance while the loop stays closed.
        value_t const dt = 2.0e-4;
        size_t const N = 10000; // 2 s
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

        CHECK(drift < 1e-4); // energy conserved over the spatial swing
        CHECK(gmax < 1e-9);  // loop stays closed (projection controls drift)
        CHECK(KEmax > 0.1);  // genuinely moving (driven by the initial spin + gravity)

        fmt::println(
            "  g_open = {:.3f} -> assembled = {:.2e}; E0 = {:.5f}, KEmax = {:.4f}, "
            "dE/scale = {:.2e}, max||g|| = {:.2e}",
            g_open, cl.residual_norm(), E0, KEmax, drift, gmax);
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: closed_loop_system3dp")
