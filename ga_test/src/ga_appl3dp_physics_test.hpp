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

    TEST_CASE("pga3dp: get_disc_inertia / make_disc_body")
    {
        fmt::println("pga3dp: get_disc_inertia (cylinder, symmetry axis e3)");

        value_t const m = 3.0, r = 2.0, t = 0.5;
        auto const I = get_disc_inertia(m, r, t);
        value_t const I_axial = 0.5 * m * r * r;                  // I_zz = 6.0
        value_t const I_trans = m * (r * r / 4.0 + t * t / 12.0); // I_xx = I_yy = 3.0625

        // linear mass block (upper-right)
        CHECK(I.view()[0, 3] == doctest::Approx(m));
        CHECK(I.view()[1, 4] == doctest::Approx(m));
        CHECK(I.view()[2, 5] == doctest::Approx(m));
        // moments (lower-left): transverse about e1,e2; axial about e3
        CHECK(I.view()[3, 0] == doctest::Approx(I_trans));
        CHECK(I.view()[4, 1] == doctest::Approx(I_trans));
        CHECK(I.view()[5, 2] == doctest::Approx(I_axial));

        // I * I_inv == identity
        auto const I_inv = get_inertia_inverse(I);
        auto const P = I.view(), Q = I_inv.view();
        for (size_t i = 0; i < 6; ++i)
            for (size_t j = 0; j < 6; ++j) {
                value_t s = 0.0;
                for (size_t k = 0; k < 6; ++k)
                    s += P[i, k] * Q[k, j];
                CHECK(s == doctest::Approx(i == j ? 1.0 : 0.0).epsilon(1e-12));
            }

        // thin-disc limit (t -> 0): the perpendicular-axis theorem I_zz = I_xx + I_yy
        // holds
        auto const I_thin = get_disc_inertia(m, r, 0.0);
        CHECK(I_thin.view()[5, 2] ==
              doctest::Approx(I_thin.view()[3, 0] + I_thin.view()[4, 1]));

        // make_disc_body wires the inertia + its inverse + mass
        auto const body = make_disc_body(m, r, t);
        CHECK(body.mass == doctest::Approx(m));
        CHECK(body.I.view()[5, 2] == doctest::Approx(I_axial));
        fmt::println("");
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

    TEST_CASE(
        "pga3dp: applied wrench + time threading - Sommerfeld in-library (Phase A.2)")
    {
        fmt::println(
            "pga3dp: dynamic_system3dp - applied wrench / Sommerfeld (Phase A.2)");

        // Phase A.2: validate the time-varying applied-wrench force element folded into
        // tau (set_applied_wrench) + the RK4 sub-step time threading, by reproducing
        // B.1's Sommerfeld forced response INSIDE the library. B.1's two radial axes are
        // independent 1-DOF damped oscillators (Bisoi Eq. 2), so each is a single
        // prismatic joint (mass m+M) with a spring/damper (Phase A.1) and a rotating
        // unbalance force m*e*Omega^2 applied as a world-frame wrench wdg(O, F(t)). The
        // steady-state amplitude must match the closed-form Eq. (4) -- the same gate as
        // standalone B.1.
        value_t const mt = 4.9 + 2.45;        // m + M [kg]
        value_t const me = 4.9 * 0.008336;    // unbalance m*e [kg m] (Bisoi Table 1)
        value_t const Kx = 2000.0, Rx = 5.0;  // e2 (x) foundation
        value_t const Ky = 6000.0, Ry = 10.0; // e3 (y) foundation

        // sign/magnitude check: a CONSTANT force F0 along e2 gives joint accel F0/mt
        {
            dynamic_system3dp s;
            s.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            s.add_frame(static_frame3dp("W"));
            s.add_prismatic_body(static_frame3dp("B"),
                                 make_cuboid_body(mt, 1.0, 1.0, 1.0),
                                 vec3dp{0.0, 1.0, 0.0, 0.0}); // slide along e2
            value_t const F0 = 3.0;
            s.set_applied_wrench(
                1, [F0](value_t) { return wdg(O_3dp, vec3dp{0.0, F0, 0.0, 0.0}); });
            // at q = 0, v = 0 the only generalised force is the applied F0 along e2
            CHECK(s.joint_accel(1) == doctest::Approx(F0 / mt));
        }

        auto amp_cf = [&](value_t Om, value_t K, value_t R) {
            value_t const d = K - mt * Om * Om;
            return me * Om * Om / std::sqrt((Om * R) * (Om * R) + d * d);
        };

        // run a single forced axis to steady state and return its vibration amplitude
        auto run_axis = [&](value_t Om, vec3dp const& dir, value_t K, value_t R,
                            bool use_sin) {
            dynamic_system3dp s;
            s.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            s.add_frame(static_frame3dp("W"));
            s.add_prismatic_body(static_frame3dp("B"),
                                 make_cuboid_body(mt, 1.0, 1.0, 1.0), dir);
            s.set_joint_spring_damper(1, K, R);
            s.set_applied_wrench(1, [me, Om, dir, use_sin](value_t t) {
                value_t const ph = use_sin ? std::sin(Om * t) : std::cos(Om * t);
                vec3dp const F = (me * Om * Om * ph) * dir;
                return wdg(O_3dp, F);
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

        value_t const om1 = std::sqrt(Kx / mt), om2 = std::sqrt(Ky / mt);
        vec3dp const e2{0.0, 1.0, 0.0, 0.0}, e3{0.0, 0.0, 1.0, 0.0};

        // x-axis at its resonance, x off-resonance, y-axis at its resonance
        struct axis_case {
            char const* nm;
            value_t Om, K, R;
            vec3dp dir;
            bool use_sin;
        };
        axis_case const cases[] = {
            {"x@om1", om1, Kx, Rx, e2, false},
            {"x@24.0", 24.0, Kx, Rx, e2, false},
            {"y@om2", om2, Ky, Ry, e3, true},
        };
        for (auto const& c : cases) {
            value_t const a = run_axis(c.Om, c.dir, c.K, c.R, c.use_sin);
            value_t const cf = amp_cf(c.Om, c.K, c.R);
            fmt::println("  {}: A_lib = {:.4e}, A_cf = {:.4e}, err% = {:.3f}", c.nm, a,
                         cf, 100.0 * std::abs(a - cf) / cf);
            CHECK(a == doctest::Approx(cf).epsilon(0.01)); // matches B.1's Eq. (4)
        }

        fmt::println("");
    }

    TEST_CASE("pga3dp: driven joint - constant-rate spin kinematics (Phase A.3)")
    {
        fmt::println("pga3dp: dynamic_system3dp - driven joint kinematics (Phase A.3)");

        // A kinematically driven revolute joint spins at constant Omega about e1; after
        // time t its coordinate is exactly Omega*t and a point at radius L on a
        // co-rotating child frame moves at speed Omega*L. Validates set_driven_rate + the
        // prescribed- motion kinematics (the joint is excluded from the dynamic DOFs).
        value_t const Om = 1.7, L = 2.0;
        dynamic_system3dp sys;
        sys.add_frame(static_frame3dp("W"));
        sys.add_revolute_body(static_frame3dp("S"), make_cuboid_body(1.0, 0.1, 0.1, 0.1),
                              O_3dp, vec3dp{1.0, 0.0, 0.0, 0.0}); // hinge about e1
        sys.add_frame(static_frame3dp("P", vec3dp{0.0, L, 0.0, 1.0}),
                      sys.index_of("S")); // co-rotating point at radius L
        size_t const S = sys.index_of("S");
        sys.set_driven_rate(S, Om);

        value_t const dt = 1.0e-3;
        for (int n = 0; n < 500; ++n)
            sys.step(dt); // 0.5 s
        value_t const t = sys.time();
        CHECK(sys.joint_phi(S) == doctest::Approx(Om * t));
        CHECK(sys.joint_omega(S) == doctest::Approx(Om));

        auto const P = unitize(move3dp(O_3dp, sys.get_pos_trafo("P", "W")));
        auto const v = sys.point_velocity(P, "P");
        CHECK(to_val(bulk_nrm(v)) == doctest::Approx(Om * L));
        fmt::println("  phi = {:.5f} (Om*t = {:.5f}); |v_P| = {:.5f} (Om*L = {:.5f})",
                     sys.joint_phi(S), Om * t, to_val(bulk_nrm(v)), Om * L);
        fmt::println("");
    }

    TEST_CASE(
        "pga3dp: driven joint - spinning radial slider (centrifugal eq., Phase A.3)")
    {
        fmt::println(
            "pga3dp: dynamic_system3dp - driven spin + radial slider (Phase A.3)");

        // A driven revolute base spins at constant Omega about e1; a radial prismatic
        // slider (along the rotating e2) carries a mass m on a spring (k) + damper (c)
        // with rest at radius L. The mass settles where the spring balances the
        // centrifugal force m Omega^2 (L+q):  q_eq = m Omega^2 L / (k - m Omega^2). The
        // centrifugal force EMERGES from the driven base's Newton-Euler bias -- nothing
        // adds it by hand.
        value_t const m = 1.0, L = 1.0, Om = 2.0, k = 20.0, c = 2.0;
        auto const cube = make_cuboid_body(m, 0.05, 0.05, 0.05); // near point mass
        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
        sys.add_frame(static_frame3dp("W"));
        sys.add_revolute_body(static_frame3dp("S"), make_cuboid_body(1.0, 0.1, 0.1, 0.1),
                              O_3dp, vec3dp{1.0, 0.0, 0.0, 0.0}); // spin hinge about e1
        sys.add_prismatic_body(static_frame3dp("B", vec3dp{0.0, L, 0.0, 1.0}), cube,
                               vec3dp{0.0, 1.0, 0.0, 0.0}); // radial slider along e2
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

    TEST_CASE("pga3dp: grounded spring/damper - stiffness emergence (Phase C.1)")
    {
        fmt::println(
            "pga3dp: dynamic_system3dp - grounded spring stiffness emergence (C.1)");

        // Phase C.1: the grounded spatial spring/damper element (add_grounded_spring).
        // Two analytic gates establish the mechanism the Tao spindle relies on -- that a
        // single linear spring at a body-fixed point yields BOTH a translational and (via
        // its lever arm) a tilt stiffness, with no separate torsional spring.

        // --- Gate 1: translational stiffness == k, free response == closed form --------
        // A rotor on a prismatic e1-slider; a grounded spring at the body origin (cm)
        // with stiffness k along e1, rest at the origin. Displaced by x0 the spring is
        // the only force, so the reduced system is  m q'' + k q = 0.
        {
            value_t const m = 2.0, k = 800.0, x0 = 0.01;
            value_t const wn = std::sqrt(k / m); // 20 rad/s
            auto const cube = make_cuboid_body(m, 0.1, 0.1, 0.1);

            dynamic_system3dp sys;
            sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            sys.add_frame(static_frame3dp("W"));
            sys.add_prismatic_body(static_frame3dp("B"), cube, vec3dp{1.0, 0.0, 0.0, 0.0},
                                   x0, 0.0); // q0 = x0
            // grounded radial spring at the cm, rest = world origin (k along e1 only)
            sys.add_grounded_spring(1, O_3dp, O_3dp, vec3dp{k, 0.0, 0.0, 0.0}, 0.0);

            // emergent translational stiffness read from the lib's own reduced model:
            //   k_emergent = -M_eff * q'' / q   (M_eff = m here)
            value_t const k_emergent = -sys.mass_matrix()[0] * sys.joint_accel(1) / x0;
            fmt::println("  Gate1 translational: k_emergent = {:.4f} (k = {:.1f})",
                         k_emergent, k);
            CHECK(k_emergent == doctest::Approx(k));
            CHECK(sys.total_energy() == doctest::Approx(0.5 * k * x0 * x0));

            // free undamped response q(t) = x0 cos(wn t) through the full RK4 path
            value_t const dt = 1.0e-4;
            value_t t = 0.0, max_err = 0.0;
            for (size_t nstep = 1; nstep <= 30000; ++nstep) { // 3 s, ~9.5 periods
                sys.step(dt);
                t += dt;
                max_err =
                    std::max(max_err, std::abs(sys.joint_phi(1) - x0 * std::cos(wn * t)));
            }
            fmt::println("  Gate1 free response: max|q - x0 cos(wn t)| = {:.2e}",
                         max_err);
            CHECK(max_err < 1e-9);
        }

        // --- Gate 2: tilt stiffness from the lever arm == k * l^2 (THE headline) -------
        // The SAME kind of grounded radial spring, now placed at an axial offset l on a
        // body that tilts about e1 through the origin: tilting by theta moves the anchor
        // laterally by ~l*theta, and the lateral spring k produces a restoring moment
        // k l^2 theta. So a torsional stiffness k_theta = k l^2 EMERGES geometrically --
        // exactly how Tao's tilt stiffness (k_x+k_y)(l1^2+l2^2) arises from the radial
        // bearings at the spindle ends.
        {
            value_t const k = 1000.0, l = 0.1, th0 = 1.0e-4;
            auto const cube = make_cuboid_body(1.0, 0.05, 0.05, 0.4);

            dynamic_system3dp sys;
            sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            sys.add_frame(static_frame3dp("W"));
            // revolute hinge about e1 through the origin (the tilt DOF), start tilted th0
            sys.add_revolute_body(static_frame3dp("B"), cube, O_3dp,
                                  vec3dp{1.0, 0.0, 0.0, 0.0}, th0, 0.0);
            // grounded spring at body point (0,0,l), lateral (e2) stiffness k, rest there
            sys.add_grounded_spring(1, vec3dp{0.0, 0.0, l, 1.0}, vec3dp{0.0, 0.0, l, 1.0},
                                    vec3dp{0.0, k, 0.0, 0.0}, 0.0);

            value_t const ktheta_emergent =
                -sys.mass_matrix()[0] * sys.joint_accel(1) / th0;
            fmt::println("  Gate2 tilt: k_theta_emergent = {:.5f} (k l^2 = {:.5f})",
                         ktheta_emergent, k * l * l);
            CHECK(ktheta_emergent == doctest::Approx(k * l * l).epsilon(1e-3));
        }

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


/////////////////////////////////////////////////////////////////////////////////////////
// PGA3DP dynamic_system3dp - Phase C: the Tao wheel-spindle vibration model
//   (Tao et al., Int. J. Mech. Sci. 232 (2022) 107620). A 6-joint serial stack on the
//   inertial root reproduces the paper's 5-DOF spindle: 3 prismatic (x,y,z) + 2 revolute
//   (theta,phi) vibration DOFs (massless frames) carrying a driven revolute SPIN joint
//   with the rotor body (mass m, principal inertias Jx=Jy, Jz, cm offset by e from the
//   spin axis). Aerostatic bearings = GROUNDED spatial springs on the NON-spinning
//   housing (the phi frame): two radial springs at axial +-Lb give the full radial
//   stiffness AND
//   -- via the lever arm -- the emergent tilt stiffness; one axial spring at the cm gives
//   the axial stiffness. All parameters are Table-1 values in strict SI. Validated
//   against the closed-form characteristic frequencies Eqs (15)-(19) and the centrifugal
//   / gyroscopic structure of Eq. (13).
//
// Lb is a CALIBRATION, not a physical bearing offset -- and the reason is now understood
// (derived in TODO/tao_eq13_derivation.md): Tao's lumped radial tilt stiffness
// (kx+ky)(l1^2+l2^2) is a factor of 4 too large vs a consistent two-bearing model
// (kx/2)(l1^2+l2^2). The 4x = 2 [Tao adds ky, but a tilt about y is resisted only by the
// x-springs] x 2 [kx is the TOTAL radial stiffness of the two parallel bearings (pinned
// by the f_x cross-check), so each bearing carries kx/2, not kx]. At the table's
// l1=l2=100 mm a consistent model gives f_theta ~34% low (~4430 vs ~6700 Hz). We absorb
// the over-count by inflating the lever to Lb = sqrt(K_tilt/kx) ~ 231 mm (~2.3x physical)
// with two kx/2 springs at +-Lb, so the emergent tilt stiffness kx Lb^2 = K_tilt
// reproduces Fig. 4. So Lb encodes Tao's 4x; the physical tilt frequency is really ~34%
// lower. Confirmed NOT caused by units (strict SI throughout) or by e (the +-e in l1,l2
// shifts l1^2+l2^2 by <1e-7). The GA model itself is correct -- it is Tao's coefficient
// that does not reduce from rigid-body mechanics.
//
// FALSIFIABLE PREDICTIONS (where the corrected model parts from the paper -- both errors
// are tilt-only, so x/y/z and all of section 3/4 are untouched; full note in
// TODO/tao_eq13_derivation.md): (a) Fig. 4(c,d) and Fig. 6(b) -- a SINGLE tilt peak, not
// the twin 6680/6727: the gyroscopic split is (Jz/Jt)*f_s ~ 1.76 Hz (unresolvable), and
// the center sits ~30-40% lower (~3.7-4.4 kHz) for a consistent tilt stiffness. (b) Tao's
// +-f_s/2 makes the peak SEPARATION grow with spin; his own measurement (sec 4.2.1: peaks
// "the same under 2250 and 2650 r/min") says it does NOT -- which supports the corrected
// (negligible-gyroscopic) picture and refutes Eq.13's gyroscopic origin. (c) Topography
// results / experimental wavelengths (Table 3) are UNCHANGED: that half feeds the surface
// model with the MEASURED f_b = 6253 Hz, never Eq.13's predicted frequency.
/////////////////////////////////////////////////////////////////////////////////////////

namespace {

// Tao Table-1 spindle performance parameters, in strict SI (m, kg, s, N). The table's
// stiffness rows read "N/um" (180, 950); a compliance of 180 um/N would be absurd and the
// frequency cross-checks only close as 180 N/um = 1.8e8 N/m. e is 2 um (not 2 mm).
struct tao_spindle_params {
    value_t m = 0.8;     // rotor mass [kg]
    value_t kx = 1.8e8;  // radial stiffness kx=ky [N/m]  (180 N/um)
    value_t kz = 9.5e8;  // axial stiffness   [N/m]       (950 N/um)
    value_t Jx = 5.4e-3; // transverse inertia Jx=Jy [kg m^2]  (5.4 g m^2)
    value_t Jz = 1.9e-4; // spin-axis inertia [kg m^2]         (0.19 g m^2)
    value_t e = 2.0e-6;  // mass eccentricity [m]              (2 um)
    value_t l1 = 0.1;    // bearing axial stations [m] (the +-e correction is negligible:
    value_t l2 = 0.1;    //   l1^2+l2^2 = 0.02 to 8 sig figs, so dropped here)
    value_t Rm = 0.05;   // rotor radius (axial-bearing lever) [m]
    value_t Rw = 0.15;   // grinding-wheel radius [m]
    value_t l3 = 0.18;   // wheel-surface-to-rotor-center axial distance [m]

    // Tao's lumped tilt stiffness (Eq.12/18/19) and the CALIBRATED radial-bearing offset
    // that makes the emergent tilt stiffness (2 * (kx/2) * Lb^2 = kx Lb^2) equal it. Lb
    // is a fit to Fig. 4, NOT a physical offset: at the table's l1=l2 a consistent
    // two-bearing model gives 1/4 of Tao's radial K_tilt -- the factor-4 over-count is
    // derived in
    // TODO/tao_eq13_derivation.md (Tao adds ky and uses the total kx per bearing).
    value_t K_tilt() const { return (kx + kx) * (l1 * l1 + l2 * l2) + kz * Rm * Rm; }
    value_t Lb() const { return std::sqrt(K_tilt() / kx); }

    // closed-form characteristic frequencies [Hz], Eqs (16)-(19) (no-spin tilt = f_th0)
    value_t f_x() const { return std::sqrt(kx / m) / (2.0 * pi); }
    value_t f_z() const { return std::sqrt(kz / m) / (2.0 * pi); }
    value_t f_th0() const { return std::sqrt(K_tilt() / Jx) / (2.0 * pi); }
};

// rotor body3dp with principal moments (Jx=Jy, Jz) about the cm, via a cuboid whose
// extents reproduce them: w=h with Izz=m(w^2+h^2)/12=Jz, Ixx=m(h^2+d^2)/12=Jx.
inline body3dp make_rotor_body(value_t m, value_t Jx, value_t Jz)
{
    value_t const w = std::sqrt(6.0 * Jz / m);
    value_t const d = std::sqrt(12.0 * Jx / m - 6.0 * Jz / m);
    return make_cuboid_body(m, w, w, d);
}

struct spindle_idx {
    size_t jx, jy, jz, jth, jph, spin, housing;
};

// Build the 6-joint Tao spindle into `sys`. Initial vibration offsets q0 = {x,y,z,th,ph};
// if spin_rate != 0 the spin joint is kinematically driven at that constant rate.
inline spindle_idx build_tao_spindle(dynamic_system3dp& sys, tao_spindle_params const& p,
                                     std::array<value_t, 5> const& q0 = {},
                                     value_t spin_rate = 0.0,
                                     std::array<value_t, 5> const& qdot0 = {},
                                     std::array<value_t, 2> const& tilt_rest = {},
                                     value_t c_bearing = 0.0)
{
    auto const massless = body3dp{};
    auto const rotor = make_rotor_body(p.m, p.Jx, p.Jz);
    sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
    sys.add_frame(static_frame3dp("W"));
    sys.add_prismatic_body(static_frame3dp("Jx"), massless, vec3dp{1, 0, 0, 0}, q0[0],
                           qdot0[0]);
    sys.add_prismatic_body(static_frame3dp("Jy"), massless, vec3dp{0, 1, 0, 0}, q0[1],
                           qdot0[1]);
    sys.add_prismatic_body(static_frame3dp("Jz"), massless, vec3dp{0, 0, 1, 0}, q0[2],
                           qdot0[2]);
    sys.add_revolute_body(static_frame3dp("Jth"), massless, O_3dp, vec3dp{1, 0, 0, 0},
                          q0[3], qdot0[3]);
    sys.add_revolute_body(static_frame3dp("Jph"), massless, O_3dp, vec3dp{0, 1, 0, 0},
                          q0[4], qdot0[4]);
    // Spin joint: axis e3 through (e,0,0) so the rotor cm (body origin) orbits at radius
    // e. If spin_rate != 0 the spin is kinematically DRIVEN (motor-clamped, Omega-ddot =
    // 0, matching Tao) -- see set_driven_rate below; assemble_mass_bias now retains a
    // driven joint's inertia as a moving base, so the clamped gyroscopic dynamics are
    // produced.
    sys.add_revolute_body(static_frame3dp("rotor"), rotor, vec3dp{p.e, 0, 0, 1},
                          vec3dp{0, 0, 1, 0}, 0.0, spin_rate);

    spindle_idx ix;
    ix.jx = sys.index_of("Jx");
    ix.jy = sys.index_of("Jy");
    ix.jz = sys.index_of("Jz");
    ix.jth = sys.index_of("Jth");
    ix.jph = sys.index_of("Jph");
    ix.spin = sys.index_of("rotor");
    ix.housing = ix.jph; // bearings on the non-spinning phi frame

    value_t const Lb = p.Lb(), kr = p.kx / 2.0;
    // Optional STATIC bearing-anchor misalignment (Zhou-flatness bridge, Phase Dyn.1):
    // offset the two radial ground anchors so the unstressed (zero-force) configuration
    // is tilted by (alpha about e1 -> theta, beta about e2 -> phi). Under that small tilt
    // the
    // +-Lb body points tip (linearized) to (+-Lb*beta, -+Lb*alpha, +-Lb); placing each
    // ground anchor there makes the tilted pose the spring rest, so the spindle SETTLES
    // to (theta,phi) = (alpha,beta) -- the static alignment deviation, now produced by
    // the dynamics rather than hand-set. Default {0,0} keeps the anchors at +-Lb (Phase C
    // byte-unchanged). c_bearing adds isotropic point damping so the settling is
    // non-oscillatory (default 0 -> Phase C byte-unchanged).
    value_t const ta = tilt_rest[0], tb = tilt_rest[1];
    sys.add_grounded_spring(ix.housing, vec3dp{0, 0, Lb, 1},
                            vec3dp{Lb * tb, -Lb * ta, Lb, 1}, vec3dp{kr, kr, 0, 0},
                            c_bearing);
    sys.add_grounded_spring(ix.housing, vec3dp{0, 0, -Lb, 1},
                            vec3dp{-Lb * tb, Lb * ta, -Lb, 1}, vec3dp{kr, kr, 0, 0},
                            c_bearing);
    sys.add_grounded_spring(ix.housing, O_3dp, O_3dp, vec3dp{0, 0, p.kz, 0}, c_bearing);
    if (spin_rate != 0.0) sys.set_driven_rate(ix.spin, spin_rate); // motor-clamped spin
    return ix;
}

// Measure the oscillation frequency [Hz] of joint `jidx` by downward zero-crossings (with
// linear sub-step interpolation) of its coordinate over `nsteps` of size `dt`.
inline value_t measure_freq(dynamic_system3dp& sys, size_t jidx, value_t dt, int nsteps)
{
    std::vector<value_t> cross;
    value_t qprev = sys.joint_phi(jidx), t = 0.0;
    for (int n = 0; n < nsteps; ++n) {
        sys.step(dt);
        t += dt;
        value_t const q = sys.joint_phi(jidx);
        if (qprev > 0.0 && q <= 0.0) // downward crossing, interpolated
            cross.push_back(t - dt * q / (q - qprev));
        qprev = q;
    }
    value_t T = 0.0;
    for (size_t i = 1; i < cross.size(); ++i)
        T += cross[i] - cross[i - 1];
    return (cross.size() - 1.0) / T; // mean period -> frequency
}

} // namespace

TEST_SUITE("PGA3DP: Tao wheel-spindle (Phase C)")
{

    TEST_CASE("pga3dp: Tao spindle - radial/axial/tilt free frequencies (Phase C.2)")
    {
        fmt::println("pga3dp: Tao spindle - characteristic frequencies (Phase C.2)");

        // No-spin free-vibration frequencies of the assembled spindle. Each vibration DOF
        // is released from a small offset (all others at rest, no spin, no damping) and
        // its frequency measured. f_x, f_z are FIRST-PRINCIPLES (no calibration): the
        // radial/axial bearing stiffnesses are the table values, so matching Eqs
        // (16),(17) genuinely validates the assembled model. f_th0 instead checks that
        // the CALIBRATED Lb reproduces Tao's tilt stiffness -- matched by construction,
        // not derived (the ~factor-4 gap to a consistent l1=l2 two-spring model is
        // unresolved; see the suite header). With spin off the DOFs decouple (gyroscopic
        // coupling vanishes).
        tao_spindle_params const p;
        fmt::println("  Lb (calibrated radial-bearing offset) = {:.4f} m (l1 = {:.3f})",
                     p.Lb(), p.l1);

        // radial f_x: release x
        {
            dynamic_system3dp sys;
            build_tao_spindle(sys, p, {1.0e-7, 0, 0, 0, 0});
            value_t const f = measure_freq(sys, sys.index_of("Jx"), 2.0e-8, 60000);
            fmt::println("  f_x = {:.1f} Hz  (Eq.16: {:.1f}, paper 2385)", f, p.f_x());
            CHECK(f == doctest::Approx(p.f_x()).epsilon(2e-3));
        }
        // axial f_z: release z
        {
            dynamic_system3dp sys;
            build_tao_spindle(sys, p, {0, 0, 1.0e-7, 0, 0});
            value_t const f = measure_freq(sys, sys.index_of("Jz"), 1.0e-8, 60000);
            fmt::println("  f_z = {:.1f} Hz  (Eq.17: {:.1f}, paper 5490)", f, p.f_z());
            CHECK(f == doctest::Approx(p.f_z()).epsilon(2e-3));
        }
        // no-spin tilt f_th0: release theta (validates the emergent calibrated tilt
        // stiff.)
        {
            dynamic_system3dp sys;
            build_tao_spindle(sys, p, {0, 0, 0, 1.0e-7, 0});
            value_t const f = measure_freq(sys, sys.index_of("Jth"), 1.0e-8, 50000);
            fmt::println("  f_th0 = {:.1f} Hz  (Eq.18/19 at f_s=0: {:.1f})", f,
                         p.f_th0());
            CHECK(f == doctest::Approx(p.f_th0()).epsilon(3e-3));
        }

        fmt::println("");
    }

    TEST_CASE("pga3dp: Tao spindle - gyroscopic + centrifugal emergence (Phase C.3)")
    {
        fmt::println("pga3dp: Tao spindle - gyroscopic + centrifugal (Phase C.3)");

        // The two "for free" payoffs of building on rigid-body GA dynamics: with the spin
        // joint DRIVEN at omega, BOTH the gyroscopic tilt coupling (J_y-J_z)omega phi-dot
        // (Eq.13) and the centrifugal unbalance forcing m e omega^2 EMERGE from the
        // Newton-Euler velocity-product bias of the spinning offset-cm rotor -- nothing
        // is added by hand. Validated instantaneously (forward dynamics at a chosen
        // state) against the closed-form coefficients, the same altitude as the A.3
        // centrifugal slider.
        tao_spindle_params const p;
        value_t const f_s = 50.0;             // spin frequency [Hz] (N_w = 3000 r/min)
        value_t const omega = 2.0 * pi * f_s; // spin rate [rad/s]

        // --- Gyroscopic coupling == textbook polar-inertia value (clamped spin) -------
        // With the spin motor-CLAMPED (driven, Omega-ddot = 0, as Tao assumes), a tilt
        // rate theta-dot induces phi-dd through the gyroscopic velocity-product bias. The
        // standard rotordynamics result is that the coupling moment is the POLAR angular
        // momentum J_z*omega, so  phi-dd = (J_z / J_x) * omega * theta-dot  (forward /
        // backward whirl split = (J_z/J_x)*omega). The GA rigid-body model reproduces
        // this textbook value exactly. (assemble_mass_bias now retains the driven rotor's
        // inertia, so the clamped case runs; the value equals the free-spin case here
        // because the spin-DOF coupling to phi vanishes at this config.)
        //
        // DEVIATION (now resolved -- derived in TODO/tao_eq13_derivation.md): Tao's Eq.13
        // uses (J_y - J_z)*omega ~ J_x*omega (the TRANSVERSE inertia), giving his +-f_s/2
        // ~ 25 Hz split. That disagrees with both standard rotordynamics and this model,
        // which give the polar-inertia whirl split (J_z/J_x)*omega ~ 1.76 Hz. It is a
        // PAPER ERROR, not an artefact of his Eq.1/2 angles: projecting the body-frame
        // Euler equations onto theta,phi, the time derivative d/dt(J*omega_c) contributes
        // a kinematic frame-rotation term +J_t*omega*phi-dot that combines with the
        // cross-product term (J_z - J_t)*omega*phi-dot to leave the correct coupling
        // J_z*omega. Tao dropped the +J_t*omega*phi-dot term, leaving (J_t - J_z) = (J_y
        // - J_z). The GA model keeps the full derivative, so it is correct.
        {
            value_t const thd0 = 1.0; // theta-dot [rad/s]
            dynamic_system3dp sys;
            auto const ix = build_tao_spindle(sys, p, {}, omega, {0, 0, 0, thd0, 0});
            value_t const phidd = sys.joint_accel(ix.jph);
            value_t const coupling =
                p.Jz / p.Jx * omega;                     // textbook (J_z/J_t)*omega [1/s]
            value_t const split = coupling / (2.0 * pi); // whirl split [Hz]
            fmt::println(
                "  gyroscopic: phi-dd/theta-dot = {:.4f} (textbook (Jz/Jx)*omega "
                "= {:.4f}); whirl split = {:.2f} Hz  [Tao Eq.13 -> ~{:.0f} Hz; "
                "deviation flagged]",
                phidd / thd0, coupling, split, f_s);
            CHECK(std::abs(phidd / thd0) == doctest::Approx(coupling).epsilon(1e-3));
        }

        // --- Centrifugal unbalance forcing: x-dd at rest, spin on ---------------------
        // Eq.13 x-equation: m x-dd + ... = m e omega^2 sin(omega t) + F_x. At t=0 (cm
        // offset along the spin-axis radial, x=0, x-dot=0, no grinding force) the radial
        // acceleration the lib returns must be the centrifugal forcing / mass = e
        // omega^2.
        {
            dynamic_system3dp sys;
            auto const ix = build_tao_spindle(sys, p, {}, omega); // only spin, all q=0
            value_t const xdd = sys.joint_accel(ix.jx);
            value_t const ydd = sys.joint_accel(ix.jy);
            value_t const pred = p.e * omega * omega; // m e omega^2 / m
            fmt::println("  centrifugal: x-dd = {:.5f} m/s^2 (predicted e*omega^2 = "
                         "{:.5f}); y-dd = {:.2e}",
                         xdd, pred, ydd);
            CHECK(std::abs(xdd) == doctest::Approx(pred).epsilon(1e-3));
            CHECK(std::abs(ydd) < 1e-3 * pred); // forcing is along x at this spin phase
        }

        fmt::println("");
    }

    TEST_CASE("pga3dp: Tao spindle - axial runout z_b = z - R_w*phi (Phase C.4)")
    {
        fmt::println("pga3dp: Tao spindle - axial runout z_b (Eq.14, Phase C.4)");

        // Eq.14: the total wheel-wafer axial displacement at the grinding rim is
        //   z_b = z - R_w * phi
        // -- the axial DOF z plus the rim's axial swing under tilt phi (a point at radius
        // R_w tips out of the grinding plane by R_w*phi). z_b is THE surface-error driver
        // in Tao's model (it sets the grinding-mark depth, Eq.26). We put the spindle in
        // a known (z, phi) configuration and confirm the world axial position of a rim
        // point at radius R_w, computed by the GA model, equals z - R_w*phi.
        tao_spindle_params const p;
        value_t const z0 = 1.0e-6, phi0 = 5.0e-6; // small axial offset + tilt [m, rad]
        dynamic_system3dp sys;
        build_tao_spindle(sys, p, {0, 0, z0, 0, phi0}); // z on Jz, phi on Jph
        size_t const jz = sys.index_of("Jz"), jph = sys.index_of("Jph");

        // world axial position of a rim point at radius R_w on the (non-spinning) housing
        vec3dp const rim_w =
            unitize(move3dp(vec3dp{p.Rw, 0.0, 0.0, 1.0}, sys.get_pos_trafo("Jph", "W")));
        value_t const zb_model = rim_w.z; // from the GA pose
        value_t const zb_eq14 = sys.joint_phi(jz) - p.Rw * std::sin(sys.joint_phi(jph));
        value_t const zb_lin = z0 - p.Rw * phi0; // Eq.14 (small angle)
        fmt::println("  z = {:.3e}, phi = {:.3e} -> z_b(model) = {:.4e}, z_b(Eq.14) = "
                     "{:.4e} (linear {:.4e})",
                     z0, phi0, zb_model, zb_eq14, zb_lin);
        CHECK(zb_model ==
              doctest::Approx(zb_eq14).epsilon(1e-9)); // GA pose == z - R_w sin(phi)
        CHECK(zb_model ==
              doctest::Approx(zb_lin).epsilon(1e-3)); // == Eq.14 to small angle

        fmt::println("");
    }

    // -------------------------------------------------------------------------------------
    // Phase Dyn.1 -- the static-deviation run: the Zhou-geometry <-> Tao-dynamics BRIDGE.
    //
    // Until now the two error sources lived in separate worlds: the global FLATNESS came
    // from a HAND-SET kinematic wafer tilt (Zhou Figs 5-7, ga_appl3dp_appl_test.hpp), and
    // the spindle DYNAMICS (Phase C) only produced free vibrations about an ideal pose.
    // This case fuses them. The alignment misalignment (alpha,beta) is injected as a
    // STATIC bearing-anchor offset (build_tao_spindle's tilt_rest), so under its own
    // springs/dampers the spindle SETTLES to a tilted equilibrium
    // (theta,phi)=(alpha,beta). That settled tilt -- now an OUTPUT of the dynamics, not
    // an input -- is then fed into the Zhou single-grain cutting-path sampler and shown
    // to reproduce the same cone/dome/bowl global-flatness shapes.
    //
    // Axis relabeling: the spindle's two independent tilt DOFs are theta (about e1) and
    // phi (about e2); the Zhou sampler tilts the wafer by (a about e2, b about e3). Each
    // spindle tilt DOF maps onto one of Zhou's two perpendicular tilt slots -- theta -> a
    // (the SYMMETRIC slot, along the grain's radial line -> CONE), phi -> b (the
    // perpendicular slot -> sign-dependent DOME/BOWL). The bridge is the
    // magnitude-and-shape correspondence: a dynamically-settled tilt of magnitude alpha
    // carves Zhou's cone, one of magnitude beta carves the dome.
    //
    // A constant infeed (prescribed axial translation on Jz) is a UNIFORM z offset of the
    // whole path -- it sets the engagement depth but not the flatness SHAPE -- so it is
    // omitted from this shape gate and enters in the Step-2 dynamic-vibration run (where
    // the time-varying z_b adds the waviness on top of this mean tilt).
    TEST_CASE(
        "pga3dp: Tao spindle - static tilt deviation -> Zhou flatness (Phase Dyn.1)")
    {
        fmt::println("pga3dp: Tao spindle - static tilt deviation -> Zhou flatness "
                     "(Phase Dyn.1)");

        tao_spindle_params const p;
        value_t const deg = pi / 180.0;
        value_t const alpha = 0.1 * deg, beta = 0.1 * deg; // target misalignment [rad]
        value_t const c_bearing = 4000.0; // ~critical tilt damping -> monotone settle
        value_t const dt = 2.0e-8;
        int const nsteps = 60000; // 1.2 ms ~ 8 tilt periods: fully settled

        // Release the spindle from the IDEAL (theta=phi=0) pose with the anchors offset
        // to the misaligned rest; integrate until it settles, and report the equilibrium
        // tilt.
        auto settle = [&](value_t a, value_t b) {
            dynamic_system3dp sys;
            auto const ix = build_tao_spindle(sys, p, {}, 0.0, {}, {a, b}, c_bearing);
            for (int n = 0; n < nsteps; ++n)
                sys.step(dt);
            struct {
                value_t th, ph, rad;
            } r;
            r.th = sys.joint_phi(ix.jth);
            r.ph = sys.joint_phi(ix.jph);
            // the radial/axial DOFs are NOT driven by a pure tilt offset (the two anchor
            // forces cancel) -- confirm they stay at the ideal
            r.rad =
                std::max({std::abs(sys.joint_phi(ix.jx)), std::abs(sys.joint_phi(ix.jy)),
                          std::abs(sys.joint_phi(ix.jz))});
            return r;
        };

        // (theta,phi) = (alpha,0): settles to the alpha-tilt equilibrium
        auto const ce = settle(alpha, 0.0);
        fmt::println("  settle(alpha,0): theta_eq = {:.5e} (target {:.5e}), phi_eq = "
                     "{:.2e}, |x,y,z| <= {:.1e} m",
                     ce.th, alpha, ce.ph, ce.rad);
        CHECK(ce.th == doctest::Approx(alpha).epsilon(1e-3)); // dynamics reproduce alpha
        CHECK(std::abs(ce.ph) < 1e-3 * alpha);                // no phi cross-coupling
        CHECK(ce.rad < 1e-7);                                 // radial/axial undisturbed

        // (theta,phi) = (0,beta): settles to the beta-tilt equilibrium
        auto const de = settle(0.0, beta);
        fmt::println("  settle(0,beta): phi_eq   = {:.5e} (target {:.5e}), theta_eq = "
                     "{:.2e}",
                     de.ph, beta, de.th);
        CHECK(de.ph == doctest::Approx(beta).epsilon(1e-3)); // dynamics reproduce beta
        CHECK(std::abs(de.th) < 1e-3 * beta);

        // Confirm (alpha,beta) is a genuine EQUILIBRIUM: with the spindle placed exactly
        // at the misaligned rest (no spin, no velocity), the tilt accelerations vanish --
        // so the settled state above is the static fixed point, not a slow transient.
        // Compare to the restoring acceleration from the ideal (theta=phi=0) pose.
        {
            dynamic_system3dp eq, off;
            auto const ie =
                build_tao_spindle(eq, p, {0, 0, 0, alpha, beta}, 0.0, {}, {alpha, beta});
            auto const io = build_tao_spindle(off, p, {}, 0.0, {}, {alpha, beta});
            value_t const a_eq = std::abs(eq.joint_accel(ie.jth));
            value_t const a_off =
                std::abs(off.joint_accel(io.jth)); // restoring from ideal
            fmt::println(
                "  equilibrium: |theta-ddot| at rest = {:.3e} rad/s^2 (restoring "
                "from ideal = {:.3e}); ratio {:.1e}",
                a_eq, a_off, a_eq / a_off);
            CHECK(a_eq <
                  1e-5 * a_off); // negligible vs the restoring scale -> equilibrium
        }

        // ---- Flatness consequence: feed the settled tilt into the Zhou sampler
        // ---------- One wheel grain (radius r1=R, axis offset L=R, half-overlap) sampled
        // in the tilted, spinning wafer frame over one wafer revolution; z = wafer-axis
        // component = profile height. Identical to the committed Zhou Figs 5-7 case, but
        // the tilt now comes from the dynamics (theta_eq -> a, phi_eq -> b). Geometry in
        // mm; only the dimensionless settled ANGLES cross over from the SI dynamics.
        double const R = 150.0, L = R, r1 = R;
        double const w1 = rpm2radps(1500.0), w2 = rpm2radps(50.0); // n2/n1 = 1/30
        auto profile = [&](double a, double b) {
            kinematic_system3dp sys;
            sys.add_frame(static_frame3dp("W"));
            sys.add_frame(static_frame3dp("wheel_rot", vec3dp{0.0, L, 0.0, 1.0}),
                          kin_state3dp{.omega = vec3dp{w1, 0.0, 0.0, 0.0}},
                          sys.index_of("W"));
            sys.add_frame(static_frame3dp("grain", vec3dp{0.0, -r1, 0.0, 1.0}),
                          kin_state3dp{}, sys.index_of("wheel_rot"));
            sys.add_frame(static_frame3dp("wafer_tilt", O_3dp, vec3dp{0.0, a, b, 0.0}),
                          kin_state3dp{}, sys.index_of("W"));
            sys.add_frame(static_frame3dp("wafer_rot"),
                          kin_state3dp{.omega = vec3dp{w2, 0.0, 0.0, 0.0}},
                          sys.index_of("wafer_tilt"));
            size_t const wheel = sys.index_of("wheel_rot"),
                         wafer = sys.index_of("wafer_rot");
            double const T = 2.0 * pi / w2;
            int const N = 6000;
            struct {
                double zmin = 1e9, zmax = -1e9, z_ctr = 0.0, rmax = 0.0;
            } pr;
            for (int k = 0; k <= N; ++k) {
                double const t = T * double(k) / double(N);
                sys.set_pose(wheel, pose3dp{vec3dp{0.0, L, 0.0, 1.0},
                                            vec3dp{w1 * t, 0.0, 0.0, 0.0}});
                sys.set_pose(wafer, pose3dp{O_3dp, vec3dp{w2 * t, 0.0, 0.0, 0.0}});
                auto const g =
                    unitize(move3dp(O_3dp, sys.get_pos_trafo("grain", "wafer_rot")));
                double const z = g.x, r = std::sqrt(g.y * g.y + g.z * g.z);
                if (r > R + 1e-6) continue;
                pr.zmin = std::min(pr.zmin, z);
                pr.zmax = std::max(pr.zmax, z);
                if (r < 3.0) pr.z_ctr = z;
                pr.rmax = std::max(pr.rmax, r);
            }
            return pr;
        };

        double const amp = R * std::tan(0.1 * deg); // 0.2618 mm unclipped rim amplitude

        auto const flat = profile(0.0, 0.0); // ideal alignment -> FLAT
        CHECK(std::abs(flat.zmin) < 1e-9);
        CHECK(std::abs(flat.zmax) < 1e-9);

        auto const cone = profile(ce.th, 0.0); // settled alpha -> CONE (symmetric)
        CHECK(std::abs(cone.zmax + cone.zmin) < 1e-3 * amp);
        CHECK(std::abs(cone.z_ctr) < 1e-3);
        CHECK(cone.zmax > 0.5 * amp);
        CHECK(cone.zmax < 1.01 * amp);
        CHECK(cone.rmax > 0.99 * R);

        auto const dome = profile(0.0, de.ph); // settled +beta -> CONVEX dome
        CHECK(dome.zmin == doctest::Approx(0.0).epsilon(1e-4));
        CHECK(dome.zmax > 0.3 * amp);
        CHECK(dome.zmax < 0.75 * cone.zmax); // ~half the cone (Zhou: beta ~ half alpha)

        auto const bowl = profile(0.0, -de.ph); // settled -beta -> CONCAVE bowl (mirror)
        CHECK(bowl.zmax == doctest::Approx(0.0).epsilon(1e-4));
        CHECK(bowl.zmin == doctest::Approx(-dome.zmax).epsilon(1e-4));

        fmt::println("  flat z in [{:+.4f},{:+.4f}] | cone [{:+.4f},{:+.4f}] | dome "
                     "[{:+.4f},{:+.4f}] | bowl [{:+.4f},{:+.4f}] mm",
                     flat.zmin, flat.zmax, cone.zmin, cone.zmax, dome.zmin, dome.zmax,
                     bowl.zmin, bowl.zmax);
        fmt::println("  -> dynamically-settled tilt reproduces Zhou's cone/dome/bowl "
                     "flatness (dome/cone = {:.2f})",
                     dome.zmax / cone.zmax);
        fmt::println("");
    }

    // -------------------------------------------------------------------------------------
    // Phase Dyn.2 -- the dynamic-vibration run: ONE live spindle run yields BOTH error
    // sources at once. On top of the Phase-Dyn.1 static misalignment (the MEAN tilt ->
    // global flatness), the spindle now spins (driven Omega, mass unbalance e) and its
    // lightly-damped tilt DOFs RING at their natural frequency, so the axial runout
    //
    //     z_b = z - R_w * sin(phi)      (Eq.14)
    //
    // becomes time-varying. The live z_b then carries both signatures of Tao's model: its
    // MEAN is the static tilt (flatness, Zhou cone/dome of Dyn.1), and its OSCILLATION at
    // the tilt natural frequency f_b (Tao Eq.22 = f_th0) carves the WAVINESS whose
    // spatial wavelength is (surface speed)/f_b (Phase D.1). The headline: f_b is now an
    // OUTPUT of the live dynamics, not the measured 6253 Hz fed into the Phase-D.1
    // surface model.
    //
    // The constant grinding infeed (feed ~ 30 um/min) is a UNIFORM slow z-ramp -- over
    // the ~1.6 ms vibration window it advances < 1 nm, i.e. a constant offset of the mean
    // cutting depth, invisible to the ripple. It sets the engagement depth (and combines
    // with z_b for the carve), not the waviness frequency, so the axial DOF is left free
    // here (z ~ 0) and the feed enters the depth in the Step-3 visualization.
    TEST_CASE(
        "pga3dp: Tao spindle - dynamic vibration: flatness + waviness (Phase Dyn.2)")
    {
        fmt::println("pga3dp: Tao spindle - dynamic vibration: flatness + waviness "
                     "(Phase Dyn.2)");

        tao_spindle_params const p;
        value_t const deg = pi / 180.0;
        value_t const beta = 0.1 * deg; // static mean tilt -> global flatness
        value_t const delta = 1.0e-5;   // tilt perturbation -> z_b ripple [rad]
        value_t const f_s = 50.0;       // 3000 r/min driven spin
        value_t const Omega = 2.0 * pi * f_s;
        value_t const dt = 2.0e-8;
        int const nsteps = 80000; // 1.6 ms ~ 11 tilt periods -> ~10 z_b crossings

        // Static misalignment via the anchor offset (mean flatness); UNDAMPED so the tilt
        // rings cleanly at its natural freq over the window (the table damping ratio
        // ~5e-5 is negligible here); driven spin + unbalance on; phi released `delta`
        // above its tilted equilibrium so z_b ripples about the mean tilt.
        dynamic_system3dp sys;
        auto const ix = build_tao_spindle(sys, p, {0, 0, 0, 0.0, beta + delta}, Omega, {},
                                          {0.0, beta}, 0.0);

        std::vector<value_t> zb, ts;
        zb.reserve(nsteps + 1);
        ts.reserve(nsteps + 1);
        value_t phi_sum = 0.0, t = 0.0;
        for (int n = 0; n <= nsteps; ++n) {
            value_t const z = sys.joint_phi(ix.jz);
            value_t const phi = sys.joint_phi(ix.jph);
            zb.push_back(z - p.Rw * std::sin(phi));
            ts.push_back(t);
            phi_sum += phi;
            if (n < nsteps) {
                sys.step(dt);
                t += dt;
            }
        }
        value_t const phi_mean = phi_sum / value_t(nsteps + 1);

        // (1) MEAN tilt == the static beta -> the GLOBAL FLATNESS survives the vibration:
        // averaging the live phi recovers the prescribed alignment tilt (Dyn.1
        // cone/dome).
        fmt::println("  mean tilt phi = {:.5e} rad (static beta = {:.5e}) -> flatness",
                     phi_mean, beta);
        CHECK(phi_mean == doctest::Approx(beta).epsilon(2e-2));

        // (2) z_b OSCILLATION frequency == the tilt natural freq f_b (Tao Eq.22 = f_th0),
        // measured from the live run by downward zero-crossings of (z_b - mean), linearly
        // interpolated.
        value_t zb_mean = 0.0;
        for (auto v : zb)
            zb_mean += v;
        zb_mean /= value_t(zb.size());
        std::vector<value_t> cross;
        for (size_t i = 1; i < zb.size(); ++i) {
            value_t const a = zb[i - 1] - zb_mean, b = zb[i] - zb_mean;
            if (a > 0.0 && b <= 0.0) cross.push_back(ts[i - 1] + dt * a / (a - b));
        }
        value_t Tsum = 0.0;
        for (size_t i = 1; i < cross.size(); ++i)
            Tsum += cross[i] - cross[i - 1];
        value_t const f_b = (cross.size() - 1.0) / Tsum;
        value_t const A_b = 0.5 * (*std::max_element(zb.begin(), zb.end()) -
                                   *std::min_element(zb.begin(), zb.end())); // [m]
        fmt::println(
            "  z_b ripple: f_b = {:.1f} Hz (Eq.22 f_th0 = {:.1f}), A_b = {:.3f} um", f_b,
            p.f_th0(), A_b * 1e6);
        CHECK(f_b == doctest::Approx(p.f_th0()).epsilon(5e-3));

        // (3) WAVINESS wavelength lambda = (surface speed)/f_b (Phase D.1), now with f_b
        // from the LIVE dynamics. WMD lambda_m uses the wheel surface speed v_w =
        // n_w*R_w; WCD lambda_c uses the wafer surface speed v_s = n_s*r at a
        // representative radius.
        value_t const v_w = rpm2radps(3000.0) * p.Rw; // wheel surface speed [m/s]
        value_t const v_s = rpm2radps(265.0) * 0.030; // wafer @ r=30 mm [m/s]
        value_t const lambda_m = v_w / f_b * 1e3;     // WMD wavelength [mm]
        value_t const lambda_c = v_s / f_b * 1e3;     // WCD wavelength [mm]
        fmt::println("  waviness from live f_b: lambda_m = {:.3f} mm (WMD), lambda_c = "
                     "{:.1f} um (WCD)",
                     lambda_m, lambda_c * 1e3);
        CHECK(lambda_m > 5.0); // WMD ~ mm scale (Tao Table 3: 5-7 mm)
        CHECK(lambda_m < 10.0);
        CHECK(lambda_c * 1e3 > 50.0); // WCD ~ 100 um scale (Tao Table 3: 127-163 um)
        CHECK(lambda_c * 1e3 < 300.0);

        fmt::println(
            "  -> one live run: MEAN tilt = flatness (Zhou cone/dome), z_b ripple "
            "at f_b = waviness (lambda = v/f_b)");
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: Tao wheel-spindle (Phase C)")

// =================================================================================
// PGA3DP dynamic_system3dp - Phase D.2b: closing the grinding force loop
//
// D.2a (ga_appl3dp_appl_test.hpp) computed the grinding force FEED-FORWARD from a fixed
// engagement. D.2b promotes the contact element into dynamic_system3dp (set_contact_force
// + the swappable grinding_law3dp + the assemble_mass_bias fold-in) and CLOSES the loop:
// the spindle's own axial compliance sets the engagement, the force law turns it into a
// wrench, the wrench feeds back. The closed loop has an analytic equilibrium, which is
// the gate.
//
// SIMPLIFIED RIG (the full Tao 5-DOF spindle of Phase C is the next refinement): a single
// axial DOF -- a prismatic joint along the wafer normal e1 carrying the rotor mass m,
// with an axial bearing spring k_spring + damper c to ground -- plus a DRIVEN wheel spin
// about e1 (massless, purely kinematic) so the contact has a real relative-sliding speed
// n_w*R_w. A static wafer frame supplies the wafer side of the relative velocity. The
// in-plane (tangential) part of the force is perpendicular to the axial DOF, so only the
// axial F_z drives the loop -- the equation of motion is
//
//     m z'' + c z' + k_spring z = F_z = k_grind (delta0 - z)        [z_ref = 0 at start]
//
// i.e. a damped oscillator with effective stiffness (k_spring + k_grind) and an offset
// drive k_grind*delta0, settling to
//
//     z_eq = k_grind delta0 / (k_spring + k_grind),   F_eq = k_grind (delta0 - z_eq).
//
// The HEADLINE: closing the loop pulls the force BELOW the feed-forward nominal
// (k_grind*delta0 = 100 N) because the spindle deflects and backs out of the cut; in the
// weak-coupling / stiff-spindle limit (k_grind << k_spring) F_eq -> the feed-forward
// value, recovering D.2a. Reaction on the wafer (chuck side) stays OFF here -- topology 1
// (D.2c).
// =================================================================================
TEST_SUITE("PGA3DP: grinding force loop (Phase D.2b)")
{

    TEST_CASE("pga3dp: closed grinding-force loop, axial spindle (Phase D.2b)")
    {
        fmt::println("");
        fmt::println("pga3dp: closed grinding-force loop, axial spindle (Phase D.2b)");
        fmt::println("");

        // SI units. Calibration from D.2a: delta0 = 1 um, k_grind = 1e8 N/m -> 100 N
        // nominal; mu = 0.3536 (Tao in-plane). Axial bearing k_spring = Tao k_z ~ 950
        // N/um.
        value_t const m = 0.8;                 // rotor mass [kg] (Tao Table 1)
        value_t const R_w = 0.15;              // grinding-wheel radius [m]
        value_t const n_w = rpm2radps(3000.0); // wheel spin [rad/s]
        value_t const delta0 = 1.0e-6;         // nominal engagement [m]
        value_t const k_grind = 1.0e8; // depth-law stiffness [N/m] (-> 100 N at delta0)
        value_t const mu = std::sqrt(2.0) * 25.0 / 100.0; // in-plane / axial = 0.3536
        value_t const k_spring = 9.5e8; // axial bearing stiffness [N/m] (Tao k_z)
        value_t const k_tot = k_spring + k_grind;
        value_t const c = 2.0 * 0.7 * std::sqrt(k_tot * m); // damping, zeta_eff = 0.7

        // analytic closed-loop equilibrium
        value_t const z_eq = k_grind * delta0 / k_tot;
        value_t const F_eq = k_grind * (delta0 - z_eq);
        value_t const F_ff = k_grind * delta0; // feed-forward nominal (D.2a) = 100 N

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0}); // isolate the axial loop
        sys.add_frame(static_frame3dp("W"));         // inertial root
        // static wafer surface (wafer side of the relative-sliding velocity)
        sys.add_frame(static_frame3dp("wafer"), sys.index_of("W"));
        // axial DOF: prismatic along e1, carries the rotor mass, bearing spring + damper
        sys.add_prismatic_body(static_frame3dp("spindle_z"),
                               make_cuboid_body(m, 0.05, 0.05, 0.05),
                               vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys.index_of("W"));
        size_t const zj = sys.index_of("spindle_z");
        sys.set_joint_spring_damper(zj, k_spring, c);
        // driven wheel spin about e1: a negligible-mass kinematic flywheel (1e-6 kg, just
        // enough for a non-singular inertia) -- it supplies the relative-sliding speed
        // n_w*R_w without materially loading the axial DOF (M[0] stays ~ m).
        sys.add_revolute_body(static_frame3dp("wheel"),
                              make_cuboid_body(1.0e-6, R_w, R_w, 0.01), O_3dp,
                              vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, zj);
        size_t const wj = sys.index_of("wheel");
        sys.set_driven_rate(wj, n_w);

        // register the grinding contact: rim point on the wheel vs the wafer, depth law
        vec3dp const rim_b{0.0, R_w, 0.0, 1.0}; // rim point in the wheel body frame
        sys.set_contact_force(
            wj, sys.index_of("wafer"),
            [=](contact_state3dp const& s) {
                return grinding_force_depth(s, k_grind, mu);
            },
            vec3dp{1.0, 0.0, 0.0, 0.0}, delta0, 30.0e-6 / 60.0, 0.8e-3, rim_b);

        // GATE 1: at the nominal pose (z = 0) the only force is the contact normal force
        // k_grind*delta0 = 100 N on the mass m -> initial axial accel 125 m/s^2
        CHECK(sys.mass_matrix()[0] == doctest::Approx(m));
        CHECK(sys.joint_accel(zj) == doctest::Approx(F_ff / m));
        // GATE 2: the contact has a real relative-sliding speed (wheel-dominated) ==
        // n_w*R_w
        vec3dp const P0 = unitize(move3dp(rim_b, sys.get_pos_trafo(wj, 0)));
        CHECK(to_val(bulk_nrm(sys.point_velocity(P0, wj))) == doctest::Approx(n_w * R_w));

        // integrate to steady state (zeta_eff = 0.7 -> settles in a few periods)
        value_t const dt = 1.0e-6;
        for (size_t nstep = 1; nstep <= 20000; ++nstep) // 20 ms
            sys.step(dt);

        value_t const z_ss = sys.joint_phi(zj);
        value_t const F_ss = sys.contact_force(wj).x; // axial component (n_hat = e1)
        value_t const d_ss = sys.contact_engagement(wj);
        fmt::println("  feed-forward nominal force : {:>7.3f} N  (k_grind*delta0)", F_ff);
        fmt::println("  closed-loop equilibrium    : z_eq = {:.4f} um, F_eq = {:>7.3f} N "
                     "(analytic)",
                     z_eq * 1.0e6, F_eq);
        fmt::println(
            "  simulated steady state     : z_ss = {:.4f} um, F_ss = {:>7.3f} N, "
            "delta_ss = {:.4f} um",
            z_ss * 1.0e6, F_ss, d_ss * 1.0e6);
        fmt::println("  -> spindle compliance backs the wheel out of the cut: "
                     "{:.1f} N -> {:.1f} N ({:+.1f}%)",
                     F_ff, F_ss, 100.0 * (F_ss - F_ff) / F_ff);

        // GATE 3: the closed loop settles to the analytic equilibrium
        CHECK(z_ss == doctest::Approx(z_eq).epsilon(1e-4));
        CHECK(F_ss == doctest::Approx(F_eq).epsilon(1e-4));
        CHECK(d_ss == doctest::Approx(delta0 - z_eq).epsilon(1e-4));
        // GATE 4: closing the loop reduced the force below the feed-forward nominal
        CHECK(F_ss < F_ff);
        CHECK(F_eq == doctest::Approx(F_ff * k_spring / k_tot)); // == feed-forward as
                                                                 // k_spring -> infinity

        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: grinding force loop (Phase D.2b)")

// =================================================================================
// PGA3DP dynamic_system3dp - Phase D.2c: wafer thinning + optional chuck-side reaction
//
// Two pieces of infrastructure on top of the D.2b loop:
//  (1) WAFER THINNING -- a quasi-static wafer thickness decremented between feed
//  macro-steps
//      by a SIMPLE removal law MRR = k_mrr * delta * v_rel (one tunable constant k_mrr,
//      lumping contact width / removal efficiency). The thickness is exposed in
//      contact_state so a thickness-dependent force law or the D.2e feed control can read
//      it.
//  (2) CHUCK-SIDE REACTION -- the optional, user-configurable equal/opposite force line
//      -wdg(P, F) on the wafer frame (set_contact_reaction). DEFAULT off = rigid wafer;
//      the realistic reaction physics is a later refinement, so this just builds + tests
//      the infrastructure (Newton's third law for a simple two-DOF case).
// =================================================================================
TEST_SUITE("PGA3DP: grinding wafer thinning + reaction (Phase D.2c)")
{

    TEST_CASE("pga3dp: quasi-static wafer thinning from MRR (Phase D.2c)")
    {
        fmt::println("");
        fmt::println("pga3dp: quasi-static wafer thinning from MRR (Phase D.2c)");
        fmt::println("");

        // the D.2b axial loop + a removal constant k_mrr. The wafer thins
        // quasi-statically (held constant within the kHz dynamics, decremented between
        // feed macro-steps).
        value_t const m = 0.8, R_w = 0.15, n_w = rpm2radps(3000.0);
        value_t const delta0 = 1.0e-6, k_grind = 1.0e8;
        value_t const mu = std::sqrt(2.0) * 25.0 / 100.0;
        value_t const k_spring = 9.5e8, k_tot = k_spring + k_grind;
        value_t const c = 2.0 * 0.7 * std::sqrt(k_tot * m);
        value_t const tw0 = 0.8e-3; // initial wafer thickness [m]
        // tunable removal constant [1/m]: MRR = k_mrr*delta*v_rel is the
        // thickness-reduction RATE [m/s]; chosen so the thinning is ~ Tao's feed scale
        // (30 um/min ~ 0.4 um/s).
        value_t const k_mrr = 0.01;

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
        sys.add_frame(static_frame3dp("W"));
        sys.add_frame(static_frame3dp("wafer"), sys.index_of("W"));
        sys.add_prismatic_body(static_frame3dp("spindle_z"),
                               make_cuboid_body(m, 0.05, 0.05, 0.05),
                               vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys.index_of("W"));
        size_t const zj = sys.index_of("spindle_z");
        sys.set_joint_spring_damper(zj, k_spring, c);
        sys.add_revolute_body(static_frame3dp("wheel"),
                              make_cuboid_body(1.0e-6, R_w, R_w, 0.01), O_3dp,
                              vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, zj);
        size_t const wj = sys.index_of("wheel");
        sys.set_driven_rate(wj, n_w);
        vec3dp const rim_b{0.0, R_w, 0.0, 1.0};
        sys.set_contact_force(
            wj, sys.index_of("wafer"),
            [=](contact_state3dp const& s) {
                return grinding_force_depth(s, k_grind, mu);
            },
            vec3dp{1.0, 0.0, 0.0, 0.0}, delta0, 30.0e-6 / 60.0, tw0, rim_b, k_mrr);

        // settle to the axial equilibrium, then thin quasi-statically
        value_t const dt = 1.0e-6;
        for (size_t nstep = 1; nstep <= 20000; ++nstep)
            sys.step(dt);

        value_t const d_eq = sys.contact_engagement(wj);
        value_t const rate = sys.removal_rate(wj); // MRR = k_mrr * delta_eq * v_rel
        CHECK(rate == doctest::Approx(k_mrr * d_eq * n_w *
                                      R_w)); // v_rel = n_w*R_w (static wafer)
        fmt::println("  delta_eq = {:.4f} um, v_rel = {:.2f} m/s, MRR = {:.4e} m/s",
                     d_eq * 1e6, n_w * R_w, rate);

        value_t const dt_macro = 1.0; // feed macro-step [s]
        value_t prev = sys.wafer_thickness(wj);
        CHECK(prev == doctest::Approx(tw0));
        fmt::println("  wafer thickness over 5 feed macro-steps (dt = {:.2f} s):",
                     dt_macro);
        for (int i = 0; i < 5; ++i) {
            sys.update_wafer_thinning(dt_macro);
            value_t const tw = sys.wafer_thickness(wj);
            CHECK(tw < prev);                                     // monotone decrease
            CHECK(prev - tw == doctest::Approx(rate * dt_macro)); // exact decrement
            fmt::println("    step {}: tw = {:.3f} um  (removed {:.4f} um)", i + 1,
                         tw * 1e6, (prev - tw) * 1e6);
            prev = tw;
        }
        fmt::println("");
    }

    TEST_CASE("pga3dp: optional chuck-side reaction wrench (Phase D.2c)")
    {
        fmt::println("");
        fmt::println("pga3dp: optional chuck-side reaction wrench (Phase D.2c)");
        fmt::println("");

        // Newton's third law for the contact: with the reaction ON, the equal/opposite
        // force line -wdg(P, F) drives a COMPLIANT wafer DOF opposite the spindle; with
        // it OFF (default rigid wafer) the wafer stays put. Two prismatic DOFs along e1
        // (spindle mass m1, wafer mass m2) facing across the contact; at t = 0 (both at
        // rest, mu = 0, delta = delta0) the contact normal force is F = k_grind*delta0 =
        // 100 N.
        value_t const m1 = 0.8, m2 = 2.0;
        value_t const delta0 = 1.0e-6, k_grind = 1.0e8,
                      mu = 0.0;             // pure normal force here
        value_t const F = k_grind * delta0; // = 100 N

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
        sys.add_frame(static_frame3dp("W"));
        sys.add_prismatic_body(static_frame3dp("spindle"),
                               make_cuboid_body(m1, 0.05, 0.05, 0.05),
                               vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys.index_of("W"));
        sys.add_prismatic_body(static_frame3dp("wafer"),
                               make_cuboid_body(m2, 0.05, 0.05, 0.05),
                               vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys.index_of("W"));
        size_t const si = sys.index_of("spindle"), wi = sys.index_of("wafer");
        sys.set_contact_force(
            si, wi,
            [=](contact_state3dp const& s) {
                return grinding_force_depth(s, k_grind, mu);
            },
            vec3dp{1.0, 0.0, 0.0, 0.0}, delta0);

        // reaction OFF (default rigid wafer): only the spindle accelerates
        CHECK(sys.joint_accel(si) == doctest::Approx(F / m1));
        CHECK(sys.joint_accel(wi) == doctest::Approx(0.0));

        // reaction ON: the wafer accelerates opposite, ratio -m1/m2 (Newton's third law)
        sys.set_contact_reaction(si, true);
        value_t const a_s = sys.joint_accel(si), a_w = sys.joint_accel(wi);
        CHECK(a_s == doctest::Approx(F / m1));
        CHECK(a_w == doctest::Approx(-F / m2));
        CHECK(a_w / a_s == doctest::Approx(-m1 / m2));
        fmt::println(
            "  reaction ON : spindle accel = {:>8.2f}, wafer accel = {:>8.2f} m/s^2"
            "  (ratio {:.3f} = -m1/m2)",
            a_s, a_w, a_w / a_s);

        // toggling back OFF restores the rigid-wafer behaviour
        sys.set_contact_reaction(si, false);
        CHECK(sys.joint_accel(wi) == doctest::Approx(0.0));
        fmt::println("  reaction OFF: wafer accel = 0 (rigid wafer)");
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: grinding wafer thinning + reaction (Phase D.2c)")

// =================================================================================
// PGA3DP dynamic_system3dp - Phase D.2d-2: integrator selector on the grinding loop
//
// dynamic_system3dp::set_integrator(rk4 | abm2) picks the time integrator for
// coupled_step (RK4 default, byte-identical to the old loop; ABM2 =
// Adams-Bashforth-Moulton 2nd order, a drop-in sharing the same forward-dynamics rhs).
// Run the D.2b axial grinding loop with both and (1) confirm both reach the analytic
// equilibrium and agree, and (2) MEASURE whether the loop is actually stiff -- the reason
// this step comes before any adaptive/implicit work (TODO/grinding.md Phase D.2d).
// Verdict: the working dt sits far inside both integrators' stability regions, so dt is
// ACCURACY-bound, not stability-bound -> NOT classically stiff at the current contact
// stiffness; an implicit solver (D.2d-4) is only needed if a near-rigid contact penalty
// injects a mode far above the ~6 kHz physics band.
// =================================================================================
TEST_SUITE("PGA3DP: grinding loop integrator (Phase D.2d-2)")
{

    TEST_CASE(
        "pga3dp: RK4 vs ABM2 on the grinding loop + stiffness verdict (Phase D.2d-2)")
    {
        fmt::println("");
        fmt::println("pga3dp: RK4 vs ABM2 on the grinding loop (Phase D.2d-2)");
        fmt::println("");

        value_t const m = 0.8, R_w = 0.15, n_w = rpm2radps(3000.0);
        value_t const delta0 = 1.0e-6, k_grind = 1.0e8;
        value_t const mu = std::sqrt(2.0) * 25.0 / 100.0;
        value_t const k_spring = 9.5e8, k_tot = k_spring + k_grind;
        value_t const c = 2.0 * 0.7 * std::sqrt(k_tot * m);
        value_t const z_eq = k_grind * delta0 / k_tot;
        value_t const F_eq = k_grind * (delta0 - z_eq); // analytic equilibrium force

        auto build = [&](integrator_kind ik) {
            auto sys = std::make_unique<dynamic_system3dp>();
            sys->set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            sys->add_frame(static_frame3dp("W"));
            sys->add_frame(static_frame3dp("wafer"), sys->index_of("W"));
            sys->add_prismatic_body(
                static_frame3dp("spindle_z"), make_cuboid_body(m, 0.05, 0.05, 0.05),
                vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, sys->index_of("W"));
            size_t const zj = sys->index_of("spindle_z");
            sys->set_joint_spring_damper(zj, k_spring, c);
            sys->add_revolute_body(static_frame3dp("wheel"),
                                   make_cuboid_body(1.0e-6, R_w, R_w, 0.01), O_3dp,
                                   vec3dp{1.0, 0.0, 0.0, 0.0}, 0.0, 0.0, zj);
            size_t const wj = sys->index_of("wheel");
            sys->set_driven_rate(wj, n_w);
            sys->set_contact_force(
                wj, sys->index_of("wafer"),
                [=](contact_state3dp const& s) {
                    return grinding_force_depth(s, k_grind, mu);
                },
                vec3dp{1.0, 0.0, 0.0, 0.0}, delta0, 30.0e-6 / 60.0, 0.8e-3,
                vec3dp{0.0, R_w, 0.0, 1.0});
            sys->set_integrator(ik);
            return sys;
        };

        value_t const dt = 1.0e-6;
        size_t const N = 20000;
        auto run = [&](integrator_kind ik) {
            auto sys = build(ik);
            for (size_t i = 0; i < N; ++i)
                sys->step(dt);
            return sys->contact_force(sys->index_of("wheel")).x; // axial force
        };
        value_t const F_rk = run(integrator_kind::rk4);
        value_t const F_ab = run(integrator_kind::abm2);

        fmt::println("  analytic equilibrium F_eq = {:.4f} N", F_eq);
        fmt::println("  RK4  : F = {:.4f} N", F_rk);
        fmt::println("  ABM2 : F = {:.4f} N   (|F_rk - F_ab| / F_eq = {:.2e})", F_ab,
                     std::abs(F_rk - F_ab) / F_eq);

        // both integrators reach the analytic equilibrium and agree to <= 1% (D.2d gate)
        CHECK(F_rk == doctest::Approx(F_eq).epsilon(1e-4));
        CHECK(F_ab == doctest::Approx(F_eq).epsilon(1e-3)); // ABM2 lower order, looser
        CHECK(std::abs(F_rk - F_ab) / F_eq < 0.01);

        // stiffness verdict: compare the working dt to the explicit-stability dt limits.
        // If dt is far below them, stability is NOT the binding constraint (dt is
        // accuracy-set to resolve the ~6 kHz vibration) -> the loop is not classically
        // stiff.
        value_t const w_eff = std::sqrt(k_tot / m); // ~36229 rad/s
        value_t const dt_rk = 2.78 / w_eff;         // RK4 stability limit (order of mag.)
        value_t const dt_ab = 1.0 / w_eff;          // ABM2 PECE stability limit (smaller)
        fmt::println("  w_eff = {:.0f} rad/s; stability dt: RK4 ~{:.2e}, ABM2 ~{:.2e}; "
                     "working dt = {:.0e}",
                     w_eff, dt_rk, dt_ab, dt);
        fmt::println(
            "  -> working dt is {:.0f}x below the ABM2 stability limit: dt is "
            "ACCURACY-bound, the loop is NOT classically stiff (k_grind = {:.0e}).",
            dt_ab / dt, k_grind);
        CHECK(dt < 0.1 * dt_ab); // stability is not binding -> not stiff at this k_grind
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: grinding loop integrator (Phase D.2d-2)")

// =================================================================================
// PGA3DP - Phase D.2e: grinding feed control (force limiting) + the throughput cost
//
// grinding_controller3dp is the optional control tier on top of the force loop: it holds
// the grinding force within a target envelope by commanding the infeed rate. This runs
// the PROCESS-level (macro-timescale) loop -- the quasi-static reduction of the fast
// force loop, in which the steady force is ~proportional to the feed (removal balance
// v_feed = MRR = k_mrr*delta*v_rel -> delta proportional to feed, F = k_grind*delta).
//
// The point of the test (and the answer to "is reducing the force to 70% just feed ->
// 70%?"): cutting the force target to a fraction cuts the feed to the same fraction AND
// the removal rate to the same fraction, so the wafer thins slower and the process takes
// ~1/fraction longer -- the throughput cost of force-limiting. Two schedules: (a)
// constant 70% cap; (b) thickness-scheduled taper 100% -> 70% between 150 and 80 um.
// =================================================================================
TEST_SUITE("PGA3DP: grinding feed control (Phase D.2e)")
{

    TEST_CASE("pga3dp: feed control -- force limiting + throughput cost (Phase D.2e)")
    {
        fmt::println("");
        fmt::println("pga3dp: grinding feed control -- force limiting + throughput cost "
                     "(Phase D.2e)");
        fmt::println("");

        value_t const F_ref = 100.0;          // reference grinding force [N] (Tao F_z)
        value_t const v_nom = 30.0e-6 / 60.0; // nominal feed 30 um/min -> [m/s]
        value_t const gain =
            F_ref / v_nom; // quasi-static F <-> feed gain (removal balance)
        value_t const tw0 = 200.0e-6, tw_end = 80.0e-6; // grind 200 um -> 80 um

        // step the slow process: feed = ctrl.feed_command(tw); steady F = gain*feed; the
        // removal rate equals the feed (steady balance) so tw -= feed*dt. Returns
        // {process time, peak force}.
        auto grind = [&](grinding_controller3dp const& ctrl) {
            value_t tw = tw0, t = 0.0, F_max = 0.0;
            value_t const dt = 1.0; // macro-step [s]
            while (tw > tw_end) {
                value_t const feed = ctrl.feed_command(tw);
                F_max = std::max(F_max, gain * feed);
                tw -= feed * dt; // removal rate = feed (steady-state)
                t += dt;
            }
            return std::pair{t, F_max};
        };

        grinding_controller3dp base{v_nom, F_ref, 1.0, 0.0, 0.0};      // no cap (100%)
        grinding_controller3dp capA{v_nom, F_ref, 0.7, 0.0, 0.0};      // (a) constant 70%
        grinding_controller3dp capB{v_nom, F_ref, 0.7, 150e-6, 80e-6}; // (b) taper
        auto const [t_base, F_base] = grind(base);
        auto const [t_a, F_a] = grind(capA);
        auto const [t_b, F_b] = grind(capB);

        fmt::println("  baseline (100%):    F = {:>5.1f} N, process = {:>4.0f} s", F_base,
                     t_base);
        fmt::println("  option a (70% cap): F = {:>5.1f} N, process = {:>4.0f} s  "
                     "({:.2f}x longer -- the throughput cost)",
                     F_a, t_a, t_a / t_base);
        fmt::println("  option b (taper):   F_max = {:>5.1f} N (full above 150um), "
                     "process = {:>4.0f} s ({:.2f}x)",
                     F_b, t_b, t_b / t_base);

        // option (a): force held at 70%, and the process is ~1/0.7 = 1.43x longer
        CHECK(F_a == doctest::Approx(0.7 * F_ref));
        CHECK(t_a / t_base == doctest::Approx(1.0 / 0.7).epsilon(0.01));
        // option (b): peak force is still 100% (above 150 um); the taper slows the late
        // phase -> longer than baseline but cheaper than the constant cap
        CHECK(F_b == doctest::Approx(F_ref));
        CHECK(t_b > t_base);
        CHECK(t_b < t_a);

        // the taper schedule itself: full above tw_hi, capped at/below tw_lo, linear
        // between
        CHECK(capB.target_fraction(200e-6) == doctest::Approx(1.0));  // above 150
        CHECK(capB.target_fraction(150e-6) == doctest::Approx(1.0));  // at 150
        CHECK(capB.target_fraction(115e-6) == doctest::Approx(0.85)); // midpoint -> 85%
        CHECK(capB.target_fraction(80e-6) == doctest::Approx(0.7));   // at 80 -> capped
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: grinding feed control (Phase D.2e)")
