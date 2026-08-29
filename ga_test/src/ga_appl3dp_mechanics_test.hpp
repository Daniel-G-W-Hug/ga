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
        // form: u(t) = u0 * rexp(-t)), state encoded in a single vec3d.

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

        // Sanity: the trajectory is approximately u0 * rexp(-T) at the end.
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

    TEST_CASE("pga3dp: moment_about - torque of a force line about a point")
    {
        fmt::println("pga3dp: moment_about - torque of a force line about a point");

        // a force f applied at P is the line F = wdg(P, f); its torque about R is the
        // documented bulk[(Q - R) ^ f] with Q any point on the line. moment_about uses
        // bulk(F) - r ^ att(F) instead, which needs no point on the line -- pin the two
        // against each other, and against the classical r x f.
        vec3dp const P{1.5, -0.7, 2.2, 1.0};
        vec3dp const f{0.3, 1.1, -0.4, 0.0};
        vec3dp const R{-2.0, 0.8, 0.5, 1.0};
        bivec3dp const F = wdg(P, f);

        bivec3dp const M = moment_about(R, F);
        bivec3dp const M_doc = bulk(wdg(vec3dp{P.x - R.x, P.y - R.y, P.z - R.z, 0.0}, f));
        CHECK(is_close(M, M_doc));

        // classical: (P - R) x f, read from the bulk slots (mx, my, mz)
        vec3dp const d{P.x - R.x, P.y - R.y, P.z - R.z, 0.0};
        CHECK(M.mx == doctest::Approx(d.y * f.z - d.z * f.y));
        CHECK(M.my == doctest::Approx(d.z * f.x - d.x * f.z));
        CHECK(M.mz == doctest::Approx(d.x * f.y - d.y * f.x));
        // a free bivector: no weight (force) part
        CHECK(std::abs(M.vx) + std::abs(M.vy) + std::abs(M.vz) < 1e-14);

        // independent of which point on the line built F, and of R's scale
        bivec3dp const F2 =
            wdg(vec3dp{P.x + 3.0 * f.x, P.y + 3.0 * f.y, P.z + 3.0 * f.z, 1.0}, f);
        CHECK(is_close(moment_about(R, F2), M));
        CHECK(is_close(moment_about(vec3dp{-4.0, 1.6, 1.0, 2.0}, F), M));

        // zero on the line of action itself
        bivec3dp const M0 = moment_about(
            vec3dp{P.x + 0.5 * f.x, P.y + 0.5 * f.y, P.z + 0.5 * f.z, 1.0}, F);
        CHECK(std::abs(M0.mx) + std::abs(M0.my) + std::abs(M0.mz) < 1e-14);

        // force lines add: the torque of a sum is the sum of the torques
        vec3dp const P2{0.2, 0.4, -1.0, 1.0};
        vec3dp const f2{-0.9, 0.2, 0.7, 0.0};
        CHECK(is_close(moment_about(R, F + wdg(P2, f2)),
                       moment_about(R, F) + moment_about(R, wdg(P2, f2))));
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
            // PGA3DP body-frame formulation: M(t) = M0 ⟇ rexp(½ B_b(t))
            // - M0: pure translation motor placing body origin (= cm) at initial world
            // pos
            //       M0 = rexp(0.5 * bivec3dp{0,0,0,cx,cy,cz})
            //       (translation uses the bulk/e23,e31,e12 components, not e41,e42,e43)
            // - B_b = phi * L_b: rotation bivector, phi = rotation angle about L_b
            // - Pivot line L_b = wdg(P_z0, P_z1) = bivec3dp{0, 0, 1, hh, -hw, 0}
            //   where P_z0 = (hw,hh,0,1), P_z1 = (hw,hh,1,1) in body frame
            //   => direction vz = +1 (along e3/z axis), moment = (hh, -hw, 0)
            // - Constraint projection: alpha = dOmega.vz  (since L_b.vz = 1)
            // - Pivot invariance: move3dp(L_b, rexp(½ phi * L_b)) = L_b for all phi
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
                // In PGA3DP, translate by (cx,cy,cz): motor = rexp(0.5 *
                // bivec3dp{0,0,0,cx,cy,cz}) Translation is encoded in the bulk
                // (e23,e31,e12) components, not the weight (e41,e42,e43) components.
                M0 = rexp(0.5 *
                          bivec3dp{0.0, 0.0, 0.0, cm_w_pos0.x, cm_w_pos0.y, cm_w_pos0.z});

                // pivot_pt_w: world position of representative pivot point (hw,hh,0,1)
                // in body frame, mapped via M0.  Stays fixed during rotation since
                // move3dp(P_on_L_b, M0 ⟇ rexp(½ phi * L_b)) = M0-image of P_on_L_b.
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

                // Current motor: M(t) = M0 ⟇ rexp(½ B_b(t))  [body-frame formulation]
                auto const M = rgpr(M0, rexp(0.5 * B));

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

                // Current motor: M(t) = M0 ⟇ rexp(½ B_b(t))
                auto const M = rgpr(M0, rexp(0.5 * B));

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
                auto const M = rgpr(M0, rexp(0.5 * u[0]));
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

            // cm follows the exact parabola  cm(t) = v0*t + 1/2 g t^2. The body-frame
            // twist couples the spin into the generator rate; with the dexp^-1 brackets
            // in step_free_body that coupling is integrated at fourth order (the
            // tolerances below were 1e-4 / 1e-3 while the scheme was second order)
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
        CHECK(max_dE / std::abs(E0) < 1e-9); // RK4 drift (measured 8.6e-12)
        CHECK(max_cm_err < 1e-8);            // measured 2.9e-11
        CHECK(max_w_err < 1e-12);            // exact: omega_dot == 0 for this case
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
        CHECK(max_dE / std::abs(E0) < 1e-9); // measured 1.9e-13
        CHECK(max_dL / L0n < 1e-9);          // measured 8.8e-12 (was 1e-5 at 2nd order)
        CHECK(max_cm < 1e-8); // cm fixed (pure rotation, no linear momentum)
        // the flip: w_y reverses sign (drops from +5 to clearly negative), and the
        // perturbation about e1 grows large during the flip
        CHECK(min_wy < -3.0);
        CHECK(max_abs_wx > 2.0);
        fmt::println("");
    }

    TEST_CASE("pga3dp: free body integration converges at FOURTH order on the group")
    {
        fmt::println("pga3dp: dynamic_system3dp - free body convergence order (dexp^-1)");

        // The order gate that the conservation tolerances above cannot provide: a
        // tolerance passes a second-order scheme at a small enough dt, a convergence
        // ratio does not. The tumbling body has B and Omega non-commuting throughout, so
        // the dexp^-1 bracket terms in step_free_body are exercised. Reference: the
        // same scheme at 32x the finest resolution, so its own error is ~1e-6 of the
        // coarsest cell's. Falsified: dropping the 1/12 double bracket makes the ratio
        // 8, dropping both brackets makes it 4 (measured).
        auto const body = make_cuboid_body(2.0, 0.1, 0.3, 0.6); // three distinct moments
        auto make = [&] {
            dynamic_system3dp sys;
            sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            sys.add_frame(static_frame3dp("W"));
            sys.add_body(static_frame3dp("B"), body,
                         kin_state3dp{.omega = vec3dp{0.05, 3.0, 0.02, 0.0}});
            return sys;
        };
        value_t const T = 3.0;
        auto const nrm6 = [](bivec3dp const& D) {
            return std::sqrt(D.vx * D.vx + D.vy * D.vy + D.vz * D.vz + D.mx * D.mx +
                             D.my * D.my + D.mz * D.mz);
        };

        auto ref = make();
        for (size_t i = 0; i < 19200; ++i)
            ref.step(T / 19200.0);
        auto const Mref = ref.get_pos_trafo(1, 0);

        value_t err_prev = 0.0, drift_prev = 0.0;
        for (size_t N : {150, 300, 600}) {
            auto sys = make();
            bivec3dp const L0 = sys.momentum_world(1);
            value_t drift = 0.0;
            for (size_t i = 0; i < N; ++i) {
                sys.step(T / value_t(N));
                drift = std::max(drift, nrm6(sys.momentum_world(1) - L0));
            }
            // pose error as the generator between the two motors
            value_t const err = nrm6(rlog(rgpr(rrev(Mref), sys.get_pos_trafo(1, 0))));
            fmt::println("  N = {:4}: pose err = {:.3e} (ratio {:5.2f}), max |dL| = "
                         "{:.3e} (ratio {:5.2f})",
                         N, err, err_prev / err, drift, drift_prev / drift);
            if (err_prev > 0.0) {
                CHECK(err_prev / err > 14.0);     // fourth order: 16 per halving
                CHECK(drift_prev / drift > 14.0); // the momentum drift shrinks alike
            }
            err_prev = err;
            drift_prev = drift;
        }
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
        // = m s^2/6 + m L^2. Small-oscillation frequency w0 = rsqrt(m g L / I_hinge).
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

    TEST_CASE("pga3dp: applied wrench + time threading - forced response in-library "
              "(Phase A.2)")
    {
        fmt::println(
            "pga3dp: dynamic_system3dp - applied wrench / forced response (Phase A.2)");

        // Phase A.2: validate the time-varying applied-wrench force element folded into
        // tau (set_applied_wrench) + the RK4 sub-step time threading, by reproducing
        // the standalone rotating-unbalance forced response INSIDE the library. Its two
        // radial axes are independent 1-DOF damped oscillators (closed-form reference),
        // so each is a single prismatic joint (mass m+M) with a spring/damper (Phase A.1)
        // and a rotating unbalance force m*e*Omega^2 applied as a world-frame wrench
        // wdg(O, F(t)). The steady-state amplitude must match the closed-form reference
        // -- the same gate as the standalone reference case.
        value_t const mt = 4.9 + 2.45;     // m + M [kg]
        value_t const me = 4.9 * 0.008336; // unbalance m*e [kg m] (reference parameters)
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
        // Two analytic gates establish the mechanism the application relies on -- that a
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
        // exactly how the application's tilt stiffness (k_x+k_y)(l1^2+l2^2) arises from
        // the radial springs at the rotor ends.
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

    TEST_CASE("pga3dp: a free body and a jointed chain must not share a path")
    {
        fmt::println("pga3dp: a free body carrying a jointed chain is a floating base");

        // History: a free body was integrated by step_free_body() and the 1-dof joints
        // by coupled_step(), and neither saw the other's inertia -- a 100 kg free body
        // carrying one revolute link produced a mass matrix byte-identical to the same
        // chain with the body absent, and fell at exactly -1/2 g t^2 while its child hung
        // off it. A guard then refused the configuration. Phase F: the free body IS a
        // 6-dof joint in the coupled solve (Featherstone 9.3), so the root's mass is
        // counted -- the case the user asked for ("not neglecting the mass of the free
        // floating root body").

        auto const heavy = make_cuboid_body(100.0, 1.0, 1.0, 1.0);
        auto const link = make_cuboid_body(1.0, 0.1, 0.1, 0.4);

        // a) jointed body BELOW a free body: 6 + 1 coordinates, and the mass matrix
        //    carries the base's 100 kg (the base's own diagonal block is its inertia)
        {
            dynamic_system3dp sys;
            sys.add_frame(static_frame3dp("I"));
            sys.add_body(static_frame3dp("base", vec3dp{0, 0, 0, 1}), heavy);
            CHECK_NOTHROW(
                sys.add_revolute_body(static_frame3dp("link", vec3dp{0, 0, 0.5, 1}), link,
                                      O_3dp, vec3dp{0, 1, 0, 0}));
            CHECK(sys.dof_coords().size() == 7);
            auto const M = sys.mass_matrix();
            REQUIRE(M.size() == 49);
            // translational block of the base (coordinates 3..5 = the bulk screws):
            // total mass 101 kg on the diagonal
            CHECK(M[3 * 7 + 3] == doctest::Approx(101.0));
            CHECK(M[4 * 7 + 4] == doctest::Approx(101.0));
            CHECK(M[5 * 7 + 5] == doctest::Approx(101.0));
            // and the whole thing falls at g: with no joint torque the base's vertical
            // acceleration is -g, and the momentum-form cross terms are consistent
            auto const qdd = sys.joint_accelerations();
            REQUIRE(qdd.size() == 7);
            CHECK(qdd[4] == doctest::Approx(-9.81).epsilon(1e-10)); // vertical (e42 slot)
        }

        // b) the same mix built in the other order: also a floating assembly now
        {
            dynamic_system3dp sys;
            sys.add_frame(static_frame3dp("I"));
            CHECK_NOTHROW(
                sys.add_revolute_body(static_frame3dp("link", vec3dp{0, 0, 0.5, 1}), link,
                                      O_3dp, vec3dp{0, 1, 0, 0}));
            CHECK_NOTHROW(sys.add_body(static_frame3dp("base", vec3dp{0, 0, 1.0, 1}),
                                       heavy, kin_state3dp{}, sys.index_of("link")));
            CHECK(sys.dof_coords().size() == 7);
        }

        // c) a MASSLESS free frame on the path is inert (step_free_body skips it), so it
        //    is allowed -- the guard tests mass, not type alone
        {
            dynamic_system3dp sys;
            sys.add_frame(static_frame3dp("I"));
            sys.add_body(static_frame3dp("carrier", vec3dp{0, 0, 0, 1}), body3dp{});
            CHECK_NOTHROW(
                sys.add_revolute_body(static_frame3dp("link", vec3dp{0, 0, 0.5, 1}), link,
                                      O_3dp, vec3dp{0, 1, 0, 0}));
        }

        // d) two INDEPENDENT subsystems hanging off the root are fine -- the check is per
        //    path, not per system. A free body on one branch, a chain on another.
        {
            dynamic_system3dp sys;
            sys.add_frame(static_frame3dp("I"));
            sys.add_body(static_frame3dp("floater", vec3dp{2, 0, 0, 1}), heavy, {}, 0);
            CHECK_NOTHROW(
                sys.add_revolute_body(static_frame3dp("arm", vec3dp{0, 0, 0.5, 1}), link,
                                      O_3dp, vec3dp{0, 1, 0, 0}, 0.0, 0.0, 0));
        }

        // e) the guard does not disturb the ordinary cases: a standalone free body, and
        //    a plain jointed chain, both still build and step
        {
            dynamic_system3dp sys;
            sys.add_frame(static_frame3dp("I"));
            CHECK_NOTHROW(sys.add_body(static_frame3dp("B", vec3dp{0, 0, 0, 1}), heavy));
            CHECK_NOTHROW(sys.step(0.001));
        }
        {
            dynamic_system3dp sys;
            sys.add_frame(static_frame3dp("I"));
            CHECK_NOTHROW(sys.add_revolute_body(static_frame3dp("l1"), link, O_3dp,
                                                vec3dp{0, 1, 0, 0}, 0.3));
            CHECK_NOTHROW(
                sys.add_revolute_body(static_frame3dp("l2", vec3dp{0, 0, 0.4, 1}), link,
                                      O_3dp, vec3dp{1, 0, 0, 0}, 0.2));
            CHECK_NOTHROW(sys.step(0.001));
        }
        fmt::println("");
    }

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

    TEST_CASE("pga3dp: dynamic_system3dp - a massless link may carry a joint (D2)")
    {
        fmt::println("pga3dp: dynamic_system3dp - massless carrier link in a chain");

        // massive -> MASSLESS -> massive with non-parallel axes: the carrier holds a
        // joint but no inertia (the composition of a universal joint from two revolutes
        // with a massless link between them). Builds, and integrates with energy
        // conserved; every joint moves some inertia so the mass matrix is regular.
        auto const heavy = make_cuboid_body(2.0, 0.6, 0.06, 0.06);
        auto const carrier = make_cuboid_body(0.0, 0.1, 0.03, 0.03); // massless
        auto const light = make_cuboid_body(1.5, 0.5, 0.05, 0.05);
        CHECK(carrier.mass == 0.0);
        vec3dp const ez{0.0, 0.0, 1.0, 0.0}, ey{0.0, 1.0, 0.0, 0.0};

        dynamic_system3dp sys;
        sys.add_frame(static_frame3dp("W"));
        sys.add_revolute_body(static_frame3dp("A", vec3dp{0.3, 0.0, 0.0, 1.0}), heavy,
                              vec3dp{-0.3, 0.0, 0.0, 1.0}, ez, 0.4, 0.0,
                              sys.index_of("W"));
        sys.add_revolute_body(static_frame3dp("C", vec3dp{0.35, 0.0, 0.0, 1.0}), carrier,
                              vec3dp{-0.05, 0.0, 0.0, 1.0}, ey, -0.7, 0.0,
                              sys.index_of("A"));
        sys.add_revolute_body(static_frame3dp("B", vec3dp{0.3, 0.0, 0.0, 1.0}), light,
                              vec3dp{-0.25, 0.0, 0.0, 1.0}, ez, 0.5, 0.0,
                              sys.index_of("C"));
        auto const M = sys.mass_matrix();
        for (size_t j = 0; j < 3; ++j)
            CHECK(M[j * 3 + j] > 0.0);

        value_t const E0 = sys.total_energy();
        value_t max_dE = 0.0;
        for (int i = 0; i < 1500; ++i) {
            sys.step(1.0e-3);
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));
        }
        fmt::println("  heavy -> massless -> light: dE/E = {:.2e}",
                     max_dE / std::abs(E0));
        CHECK(max_dE / std::abs(E0) < 1e-6);

        // a massless LEAF with a dof joint is refused with the cause; driven it runs
        dynamic_system3dp leaf;
        leaf.add_frame(static_frame3dp("W"));
        leaf.add_revolute_body(static_frame3dp("A", vec3dp{0.3, 0.0, 0.0, 1.0}), heavy,
                               vec3dp{-0.3, 0.0, 0.0, 1.0}, ez, 0.4, 0.0,
                               leaf.index_of("W"));
        leaf.add_revolute_body(static_frame3dp("C", vec3dp{0.35, 0.0, 0.0, 1.0}), carrier,
                               vec3dp{-0.05, 0.0, 0.0, 1.0}, ey, -0.7, 0.0,
                               leaf.index_of("A"));
        CHECK_THROWS_WITH_AS(leaf.step(1.0e-3), doctest::Contains("moves no inertia"),
                             std::runtime_error);
        leaf.set_driven_rate(leaf.index_of("C"), 2.0);
        CHECK_NOTHROW(leaf.step(1.0e-3));

        // a massless FREE body is inert
        dynamic_system3dp free;
        free.add_frame(static_frame3dp("W"));
        free.add_body(static_frame3dp("ghost", vec3dp{1.0, 1.0, 0.0, 1.0}), carrier);
        CHECK_NOTHROW(free.step(1.0e-3));
        fmt::println("  massless leaf joint refused with the cause; driven it runs; a "
                     "massless free body is inert");
        fmt::println("");
    }

    TEST_CASE("pga3dp: make_body_from_inertia and make_point_body3dp")
    {
        fmt::println("pga3dp: body builders from a tensor and for a point mass");

        // the tensor builder reproduces the shape builder from its principal moments,
        // byte for byte
        value_t const m = 2.0, w = 0.6, h = 0.2, d = 0.1;
        auto const cub = make_cuboid_body(m, w, h, d);
        auto const gen = make_body_from_inertia(m, m * (h * h + d * d) / 12.0,
                                                m * (w * w + d * d) / 12.0,
                                                m * (w * w + h * h) / 12.0);
        for (size_t i = 0; i < 36; ++i) {
            CHECK(gen.I.data[i] == cub.I.data[i]);
            CHECK(gen.I_inv.data[i] == cub.I_inv.data[i]);
        }

        // with products of inertia the angular block is the classical L = J omega
        value_t const Ixx = 0.5, Iyy = 0.7, Izz = 0.9, Ixy = 0.1, Ixz = -0.05, Iyz = 0.2;
        auto const b = make_body_from_inertia(3.0, Ixx, Iyy, Izz, Ixy, Ixz, Iyz);
        twist3dp const Om{0.3, -1.1, 0.7, 0.0, 0.0, 0.0}; // pure angular velocity
        bivec3dp const L = b.I(Om);
        CHECK(L.mx == doctest::Approx(Ixx * Om.vx - Ixy * Om.vy - Ixz * Om.vz));
        CHECK(L.my == doctest::Approx(-Ixy * Om.vx + Iyy * Om.vy - Iyz * Om.vz));
        CHECK(L.mz == doctest::Approx(-Ixz * Om.vx - Iyz * Om.vy + Izz * Om.vz));
        CHECK(std::abs(L.vx) + std::abs(L.vy) + std::abs(L.vz) < 1e-15); // no linear part
        // and the cached inverse is the inverse
        twist3dp const back = b.I_inv(L);
        CHECK(back.vx == doctest::Approx(Om.vx));
        CHECK(back.vy == doctest::Approx(Om.vy));
        CHECK(back.vz == doctest::Approx(Om.vz));

        // point mass on a massless rod = the simple pendulum: the initial angular
        // acceleration is exactly -g sin(theta) / L, and the swing conserves energy.
        // The point body carries the mass, the rod (a massless carrier) the hinge.
        value_t const g = 9.81, Lr = 0.8, th0 = 0.6;
        dynamic_system3dp sys;
        sys.add_frame(static_frame3dp("W"));
        // the pendulum is ONE revolute body -- the point mass -- with its hinge L away in
        // body coordinates: the massless rod is the pivot offset, not a body of its own.
        // (A rigidly attached mass is the link's body; adding it as a second, unjointed
        // body on the path would be a FREE body, which the tier refuses until Phase F.)
        sys.add_revolute_body(static_frame3dp("bob", vec3dp{0.0, -Lr, 0.0, 1.0}),
                              make_point_body3dp(1.7), vec3dp{0.0, Lr, 0.0, 1.0},
                              vec3dp{0.0, 0.0, 1.0, 0.0}, th0, 0.0, sys.index_of("W"));
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

    TEST_CASE("pga3dp: screw_axis - the screw of a motor and of a twist (D3)")
    {
        fmt::println("pga3dp: screw_axis - Chasles' axis, angle and distance");

        // round trip through the builder: get_motor(l, theta, dist) -> screw_axis gives
        // (l, theta, dist) back -- the FULL angle and distance (rlog carries the halves)
        bivec3dp const l =
            unitize(wdg(vec3dp{0.3, -0.2, 0.5, 1.0}, vec3dp{0.6, 0.8, 0.0, 0.0}));
        value_t const theta = 0.9, dist = 0.35;
        auto const M = get_motor(l, theta, dist);
        auto const s = screw_axis(M);
        CHECK(s.angle == doctest::Approx(theta).epsilon(1e-12));
        CHECK(s.dist == doctest::Approx(dist).epsilon(1e-12));
        CHECK(is_close(s.axis, l));
        CHECK(s.pitch() == doctest::Approx(dist / theta).epsilon(1e-12));
        // and forward again
        CHECK(is_same_motion(get_motor(s.axis, s.angle, s.dist), M));

        // a joint screw is a twist: a revolute's is its pivot line at unit rate, pitch 0
        dynamic_system3dp sys;
        sys.add_frame(static_frame3dp("W"));
        vec3dp const piv{-0.3, 0.1, 0.0, 1.0}, ax{0.0, 0.6, 0.8, 0.0};
        sys.add_revolute_body(static_frame3dp("L", vec3dp{0.5, 0.2, -0.1, 1.0}),
                              make_cuboid_body(1.0, 0.6, 0.1, 0.1), piv, ax, 0.7, 0.0,
                              sys.index_of("W"));
        auto const js = screw_axis(sys.joint_props(1).screw_b);
        CHECK(js.angle == doctest::Approx(1.0));
        CHECK(std::abs(js.dist) < 1e-14);
        CHECK(is_close(js.axis, unitize(wdg(piv, ax))));

        // a pure translation has no finite axis: angle 0, dist = the length, axis = the
        // direction at infinity
        auto const T = get_motor(vec3dp{0.3, -0.4, 1.2, 0.0});
        auto const st = screw_axis(T);
        CHECK(st.angle == 0.0);
        CHECK(st.dist == doctest::Approx(std::sqrt(0.09 + 0.16 + 1.44)).epsilon(1e-12));
        fmt::println("  motor round trip exact; revolute screw = pivot line, pitch 0; "
                     "translation: angle 0, dist = |delta|");
        fmt::println("");
    }

    TEST_CASE("pga3dp: jacobian - analytic columns vs finite differences of rlog (D3)")
    {
        fmt::println("pga3dp: jacobian - analytic vs finite differences, body and space");

        // the 6R arm of the weld test; the hand's body Jacobian column j must equal the
        // finite difference of the hand motor's log w.r.t. q_j:
        //
        //     J_b[:, j] = d/dq_j [ 2 rlog( rrev(M_f(q)) (x) M_f(q + h e_j) ) ]      (h ->
        //     0)
        //
        // central differences, so the truncation error is O(h^2); and the space form is
        // the body form transported by the hand motor, J_s = Ad(M_f) J_b.
        auto const seg = make_cuboid_body(2.0, 0.6, 0.06, 0.06);
        vec3dp const ez{0.0, 0.0, 1.0, 0.0}, ey{0.0, 1.0, 0.0, 0.0};
        value_t const q0[6] = {0.3, -0.5, 0.7, 0.4, -0.6, 0.2};
        dynamic_system3dp sys;
        sys.add_frame(static_frame3dp("W"));
        size_t parent = sys.index_of("W");
        for (int i = 0; i < 6; ++i) {
            sys.add_revolute_body(
                static_frame3dp("L" + std::to_string(i), vec3dp{0.6, 0.0, 0.0, 1.0}), seg,
                vec3dp{-0.3, 0.0, 0.0, 1.0}, (i % 2 == 0) ? ez : ey, q0[i], 0.0, parent);
            parent = sys.index_of("L" + std::to_string(i));
        }
        size_t const hand = parent;
        auto const rj = sys.dof_joints();
        CHECK(rj.size() == 6);

        auto const Jb = sys.jacobian_columns(hand, true);
        auto const Js = sys.jacobian_columns(hand, false);
        mvec3dp_e const Mf = sys.get_pos_trafo(hand, 0);
        value_t const h = 1.0e-5;
        value_t max_err = 0.0;
        for (size_t j = 0; j < 6; ++j) {
            value_t const q = sys.joint_phi(rj[j]);
            sys.set_joint(rj[j], q + h);
            mvec3dp_e const Mp = sys.get_pos_trafo(hand, 0);
            sys.set_joint(rj[j], q - h);
            mvec3dp_e const Mm = sys.get_pos_trafo(hand, 0);
            sys.set_joint(rj[j], q);
            // body twist over the step: the log of the relative motor, central
            twist3dp const fd = (1.0 / h) * rlog(rgpr(rrev(Mm), Mp)); // 2 rlog / (2 h)
            twist3dp const d = fd - Jb[j];
            max_err =
                std::max(max_err, std::sqrt(d.vx * d.vx + d.vy * d.vy + d.vz * d.vz +
                                            d.mx * d.mx + d.my * d.my + d.mz * d.mz));
            // space column = body column transported by the hand motor
            CHECK(is_close(move3dp(Jb[j], Mf), Js[j]));
        }
        fmt::println("  max |J_b column - FD of rlog| = {:.2e} (h = {:.0e})", max_err, h);
        CHECK(max_err < 1e-8);

        // V = J q-dot: with rates set, the space Jacobian reproduces twist_world
        value_t const qd[6] = {0.5, -0.3, 0.8, 0.1, -0.7, 0.4};
        for (size_t j = 0; j < 6; ++j) {
            sys.set_joint_rate(rj[j], qd[j]);
        }
        auto const J = sys.jacobian(hand);
        twist3dp V{};
        for (size_t j = 0; j < 6; ++j)
            V = V + qd[j] * Js[j];
        CHECK(is_close(V, sys.twist_world(hand)));
        CHECK(J.size() == 36);
        CHECK(J[0 * 6 + 2] == Js[2].vx); // row-major: row = component, column = joint
        // a column subset: joints 1 and 4 only, 6 x 2
        auto const Jsub = sys.jacobian(hand, std::vector<size_t>{rj[1], rj[4]});
        CHECK(Jsub.size() == 12);
        CHECK(Jsub[3 * 2 + 1] == Js[4].mx);
        // a joint that is no ancestor of the queried frame has a zero column
        auto const J2 = sys.jacobian_columns(rj[2]);
        CHECK(std::abs(J2[5].vx) + std::abs(J2[5].mx) == 0.0);
        fmt::println("");
    }

    TEST_CASE("pga3dp: screw_system - span and rcmt-closure, mobility (W1)")
    {
        fmt::println("pga3dp: screw_system - span vs. Lie closure of joint screws");

        // Selig ch. 8 / Featherstone 8.10: the closure of a screw system under the
        // bracket tells momentary from full-cycle mobility. The bracket is rcmt.
        vec3dp const O{0.0, 0.0, 0.0, 1.0}, P{1.0, 0.0, 0.0, 1.0}, Q{0.0, 1.0, 0.5, 1.0};
        vec3dp const ez{0.0, 0.0, 1.0, 0.0}, ex{1.0, 0.0, 0.0, 0.0};
        bivec3dp const L1 = wdg(O, ez), L2 = wdg(P, ez), L3 = wdg(Q, ez), L4 = wdg(Q, ex);

        auto const one = screw_system({L1});
        CHECK(one.span == 1);
        CHECK(one.closure == 1); // a single screw is its own subalgebra

        // two parallel axes: span 2, bracket = the translation perpendicular to the
        // plane through them -> the planar subalgebra, closure 3, and it stops there
        auto const par = screw_system({L1, L2});
        CHECK(par.span == 2);
        CHECK(par.closure == 3);
        // three parallel axes (a planar linkage) already span the planar subalgebra
        auto const planar = screw_system({L1, L2, L3});
        CHECK(planar.span == 3);
        CHECK(planar.closure == 3);
        // two skew screws generate all of se(3)
        auto const skew = screw_system({L1, L4});
        CHECK(skew.span == 2);
        CHECK(skew.closure == 6);
        // a dependent set: the same axis twice counts once
        auto const dep = screw_system({L1, 2.0 * L1});
        CHECK(dep.span == 1);
        fmt::println("  one: ({},{})  parallel pair: ({},{})  planar triple: ({},{})  "
                     "skew pair: ({},{})",
                     one.span, one.closure, par.span, par.closure, planar.span,
                     planar.closure, skew.span, skew.closure);
        fmt::println("");
    }

    TEST_CASE("pga3dp: joint torque - an actuator as a generalised force")
    {
        fmt::println("pga3dp: dynamic_system3dp - joint torque (actuator)");

        // a gravity pendulum about the z axis: a bar of mass m, centre l below the pivot
        value_t const m = 2.0, l = 0.5, g = 9.81, w = 0.05;
        value_t const I_piv = m * (w * w + 4.0 * l * l) / 12.0 + m * l * l;
        auto make = [&](value_t gg) {
            dynamic_system3dp s;
            s.set_gravity(vec3dp{0.0, -gg, 0.0, 0.0});
            s.add_frame(static_frame3dp("W"));
            s.add_revolute_body(static_frame3dp("P", vec3dp{0.0, -l, 0.0, 1.0},
                                                vec3dp{0.0, 0.0, 0.0, 0.0}),
                                make_cuboid_body(m, w, 2.0 * l, w),
                                vec3dp{0.0, l, 0.0, 1.0}, vec3dp{0.0, 0.0, 1.0, 0.0});
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
        // the change of total energy is the actuator's work integral(tau q-dot dt); the
        // trapezoid integral is second order: the mismatch falls 4x per halving of dt
        {
            auto work_err = [&](value_t dt) {
                auto s = make(g);
                s.set_joint_torque(1, [](value_t t) { return 1.5 * std::sin(3.0 * t); });
                value_t const E0 = s.total_energy();
                value_t W = 0.0;
                for (int i = 0; i * dt < 2.0; ++i) {
                    value_t const p = 1.5 * std::sin(3.0 * s.time()) * s.joint_omega(1);
                    s.step(dt);
                    value_t const p_new =
                        1.5 * std::sin(3.0 * s.time()) * s.joint_omega(1);
                    W += 0.5 * (p + p_new) * dt;
                }
                return std::abs(s.total_energy() - E0 - W);
            };
            value_t const e1 = work_err(2e-3), e2 = work_err(1e-3);
            CHECK(e1 > 1e-9);
            CHECK(e1 / e2 == doctest::Approx(4.0).epsilon(0.1));
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
        // the motor-joint form: a spherical joint (3 screws) with one generalised force
        // per screw; the solve returns M q-ddot = tau (the fold lands on the right index)
        {
            dynamic_system3dp s;
            s.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            s.add_frame(static_frame3dp("W"));
            s.add_spherical_body(static_frame3dp("B", vec3dp{0.0, -0.4, 0.0, 1.0},
                                                 vec3dp{0.0, 0.0, 0.0, 0.0}),
                                 make_cuboid_body(3.0, 0.2, 0.8, 0.1),
                                 vec3dp{0.0, 0.4, 0.0, 1.0});
            std::vector<value_t> const tau{0.3, -2.0, 0.5};
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

    TEST_CASE("pga3dp: mass_bias / constraint_jacobian - the EoM read off (L4)")
    {
        fmt::println("pga3dp: mass_bias() and constraint_jacobian() accessors");

        // the public read-off of  M(q) q-ddot = RHS + tau - G^T lambda  for an external
        // control law (computed torque) -- the 3D mirror of the 2D case. Falsified the
        // same way (forcing with_joint_torques = true; flipping an FD column sign).

        value_t const m = 2.0, l = 0.5, g = 9.81, w = 0.05;

        // -- open chain: a pendulum carrying every PASSIVE force element
        {
            dynamic_system3dp s;
            s.set_gravity(vec3dp{0.0, -g, 0.0, 0.0});
            s.add_frame(static_frame3dp("W"));
            s.add_revolute_body(static_frame3dp("P", vec3dp{0.0, -l, 0.0, 1.0},
                                                vec3dp{0.0, 0.0, 0.0, 0.0}),
                                make_cuboid_body(m, w, 2.0 * l, w),
                                vec3dp{0.0, l, 0.0, 1.0}, vec3dp{0.0, 0.0, 1.0, 0.0}, 0.4,
                                0.7);
            s.set_joint_spring_damper(1, 3.0, 0.2, 0.1);
            s.set_applied_wrench(
                1, [](value_t) { return bivec3dp{0.4, 0.0, -0.6, 0.1, 0.2, 0.0}; });

            auto const [M, RHS] = s.mass_bias();
            REQUIRE(M.size() == 1);
            REQUIRE(RHS.size() == 1);

            auto const Mm = s.mass_matrix();
            CHECK(M[0] == doctest::Approx(Mm[0]).epsilon(1e-13));

            auto const qdd = lu_solve(M, RHS, 1);
            auto const qdd_sys = s.joint_accelerations();
            CHECK(qdd[0] == doctest::Approx(qdd_sys[0]).epsilon(1e-12));

            // the EXCLUSION: an actuator torque leaves RHS bit-identical, the
            // accelerations shift by exactly M^-1 tau
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

        // -- a kinematically DRIVEN child carrying inertia: mass_bias()'s M retains the
        //    moving-base inertia, which mass_matrix() does not include
        {
            dynamic_system3dp s;
            s.set_gravity(vec3dp{0.0, -g, 0.0, 0.0});
            s.add_frame(static_frame3dp("W"));
            s.add_revolute_body(static_frame3dp("P", vec3dp{0.0, -l, 0.0, 1.0},
                                                vec3dp{0.0, 0.0, 0.0, 0.0}),
                                make_cuboid_body(m, w, 2.0 * l, w),
                                vec3dp{0.0, l, 0.0, 1.0}, vec3dp{0.0, 0.0, 1.0, 0.0}, 0.3,
                                0.0);
            s.add_revolute_body(static_frame3dp("D", vec3dp{0.0, -l, 0.0, 1.0},
                                                vec3dp{0.0, 0.0, 0.0, 0.0}),
                                make_cuboid_body(1.5, 0.3, 0.3, 0.05),
                                vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{0.0, 1.0, 0.0, 0.0},
                                0.0, 0.0, 1);
            s.set_driven_rate(2, 5.0);
            auto const [M, RHS] = s.mass_bias();
            REQUIRE(M.size() == 1);
            CHECK(M[0] > s.mass_matrix()[0]);
            auto const qdd = lu_solve(M, RHS, 1);
            CHECK(qdd[0] == doctest::Approx(s.joint_accelerations()[0]).epsilon(1e-12));
        }

        // -- closed loop: constraint_jacobian() against central finite differences of
        //    the residual, coordinate by coordinate, with all three constraint kinds
        //    active at once (coincidence 3 rows, distance 1, frame 6) on a spatial
        //    3R chain whose axes do not align
        {
            auto const seg = make_cuboid_body(2.0, 0.6, 0.06, 0.06);
            auto unit = [](vec3dp const& p) {
                return vec3dp{p.x / p.w, p.y / p.w, p.z / p.w, 1.0};
            };
            closed_loop_system3dp cl;
            cl.add_frame(static_frame3dp("W"));
            cl.add_revolute_body(static_frame3dp("b1", vec3dp{0.3, 0.0, 0.0, 1.0},
                                                 vec3dp{0.0, 0.0, 0.0, 0.0}),
                                 seg, vec3dp{-0.3, 0.0, 0.0, 1.0},
                                 vec3dp{0.0, 0.0, 1.0, 0.0}, 0.35, 0.0, cl.index_of("W"));
            cl.add_revolute_body(static_frame3dp("b2", vec3dp{0.6, 0.0, 0.0, 1.0},
                                                 vec3dp{0.0, 0.0, 0.0, 0.0}),
                                 seg, vec3dp{-0.3, 0.0, 0.0, 1.0},
                                 vec3dp{0.0, 1.0, 0.0, 0.0}, -0.55, 0.0,
                                 cl.index_of("b1"));
            cl.add_revolute_body(static_frame3dp("b3", vec3dp{0.6, 0.0, 0.0, 1.0},
                                                 vec3dp{0.0, 0.0, 0.0, 0.0}),
                                 seg, vec3dp{-0.3, 0.0, 0.0, 1.0},
                                 vec3dp{1.0, 0.0, 0.0, 0.0}, 0.25, 0.0,
                                 cl.index_of("b2"));
            size_t const W = cl.index_of("W");
            size_t const b2 = cl.index_of("b2"), b3 = cl.index_of("b3");
            vec3dp const tip_b{0.3, 0.0, 0.0, 1.0};
            vec3dp const t3 = unit(move3dp(tip_b, cl.system().get_pos_trafo(b3, 0)));
            vec3dp const t2 = unit(move3dp(tip_b, cl.system().get_pos_trafo(b2, 0)));
            vec3dp const a2{t2.x + 0.4, t2.y + 0.5, t2.z - 0.3, 1.0};
            value_t const L2 =
                std::sqrt((t2.x - a2.x) * (t2.x - a2.x) + (t2.y - a2.y) * (t2.y - a2.y) +
                          (t2.z - a2.z) * (t2.z - a2.z));
            cl.add_loop_constraint(
                loop_constraint3dp{b3, tip_b, W, t3, constraint3dp::coincidence});
            cl.add_loop_constraint(
                loop_constraint3dp{b2, tip_b, W, a2, constraint3dp::distance, L2});
            // the weld's ground side is a static frame AT b3's current pose, so the
            // orientation residual (a rotation VECTOR, the log) is zero at the test
            // configuration -- G is the velocity-level Jacobian, and d(log)/dq equals
            // the relative angular rate exactly at closure (dexp^-1(0) = id); away from
            // closure the two differ by the dexp factor, which is G's operating
            // condition (GGL keeps g ~ 0), not a defect
            auto const pw = pose3dp_from_motor(cl.system().get_pos_trafo(b3, 0));
            cl.add_frame(static_frame3dp("Wg", pw.origin, pw.rot), W);
            cl.add_loop_constraint(loop_constraint3dp{b3, tip_b, cl.index_of("Wg"), tip_b,
                                                      constraint3dp::frame});

            size_t const mrows = cl.constraint_rows();
            auto const rj = cl.system().dof_joints();
            size_t const n = rj.size();
            REQUIRE(mrows == 10); // 3 + 1 + 6
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

    TEST_CASE("pga3dp: dynamic_system3dp - total mass, centre of mass, gravity wrench")
    {
        fmt::println("pga3dp: mass distribution -- total_mass / centre_of_mass / "
                     "gravity_wrench");
        dynamic_system3dp sys;
        sys.add_frame(static_frame3dp("W"));
        sys.add_body(static_frame3dp("hub", vec3dp{1.0, 2.0, 3.0, 1.0}),
                     make_cuboid_body(4.0, 0.2, 0.2, 0.2));
        sys.add_revolute_body(static_frame3dp("link", vec3dp{0.5, 0.0, 0.0, 1.0}),
                              make_cuboid_body(2.0, 0.5, 0.05, 0.05),
                              vec3dp{-0.25, 0.0, 0.0, 1.0}, vec3dp{0.0, 0.0, 1.0, 0.0},
                              0.0, 0.0, sys.index_of("hub"));
        CHECK(sys.total_mass() == doctest::Approx(6.0).epsilon(1e-12));
        vec3dp const c = sys.centre_of_mass();
        CHECK(c.x == doctest::Approx((4.0 * 1.0 + 2.0 * 1.5) / 6.0).epsilon(1e-12));
        CHECK(c.y == doctest::Approx(2.0).epsilon(1e-12));
        CHECK(c.z == doctest::Approx(3.0).epsilon(1e-12));
        bivec3dp const W = sys.gravity_wrench();
        bivec3dp const Mc = moment_about(c, W);
        CHECK(std::abs(Mc.mx) + std::abs(Mc.my) + std::abs(Mc.mz) < 1e-12);
        // the identity behind the centre of mass: the summed force lines have the same
        // moment about ANY point as ONE force line of the total weight through the CoM
        bivec3dp const W1 = wdg(c, 6.0 * sys.gravity());
        for (vec3dp const& R : {O_3dp, vec3dp{-2.0, 0.5, 4.0, 1.0}}) {
            bivec3dp const Mo = moment_about(R, W), M1 = moment_about(R, W1);
            CHECK(Mo.mx == doctest::Approx(M1.mx).epsilon(1e-12));
            CHECK(Mo.my == doctest::Approx(M1.my).epsilon(1e-12));
            CHECK(Mo.mz == doctest::Approx(M1.mz).epsilon(1e-12));
        }
        fmt::println("  m = 6 kg, CoM = ({:.4f}, {:.4f}, {:.4f})", c.x, c.y, c.z);
        fmt::println("");
    }

    TEST_CASE("pga3dp: closed_loop_system3dp - the impact map (W2)")
    {
        fmt::println("pga3dp: closed_loop_system3dp - the impact map (W2)");

        // The 2D rod-onto-pin impact in space: a slender cuboid (free body, cm at the
        // origin, length L along x) falls along -y with no rotation and no gravity; its
        // -x end reaches a fixed ball joint (a coincidence constraint: three rows). The
        // plastic impact pins the end and the rod rotates about z through the pin --
        // angular momentum about the pin is conserved through the impact:
        //
        //     m v L/2 = (I_zz + m L^2/4) w,     v_cm+ = w L/2,     T+ = 1/2 I_pin w^2,
        //
        // the impulse is the momentum change m (v_cm+ - v_cm-) = -Λ, and the other two
        // rotations stay zero (no angular momentum about x or y to conserve).
        value_t const m = 2.0, L = 1.0, h = 0.02, v = 1.0;
        auto const rod = make_cuboid_body(m, L, h, h);
        value_t const I_pin = m * (L * L + h * h) / 12.0 + m * L * L / 4.0;
        value_t const w_ref = -m * v * (L / 2.0) / I_pin;
        vec3dp const end_b{-0.5, 0.0, 0.0, 1.0};
        auto build = [&](value_t vy) {
            closed_loop_system3dp cl;
            cl.system().set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            cl.add_frame(static_frame3dp("W"));
            cl.add_body(static_frame3dp("rod", vec3dp{0.0, 0.0, 0.0, 1.0}), rod,
                        kin_state3dp{.vel = vec3dp{0.0, vy, 0.0, 0.0}});
            size_t const c = cl.add_loop_constraint(
                loop_constraint3dp{cl.index_of("rod"), end_b, cl.index_of("W"), end_b,
                                   constraint3dp::coincidence});
            cl.set_loop_active(c, false);
            return std::pair{std::move(cl), c};
        };
        auto vel_at = [](closed_loop_system3dp& cl, vec3dp const& Pb) {
            size_t const r = cl.index_of("rod");
            vec3dp const P = unitize(move3dp(Pb, cl.system().get_pos_trafo(r, 0)));
            return cl.system().point_velocity(P, r);
        };

        { // plastic touchdown: the analytic post-impact state
            auto [cl, c] = build(-v);
            value_t const T0 = cl.system().total_energy();
            std::vector<value_t> const Lam = cl.activate_loop_with_impact(c, 0.0);
            REQUIRE(Lam.size() == 3);
            vec3dp const ve = vel_at(cl, end_b);
            CHECK(std::abs(ve.x) + std::abs(ve.y) + std::abs(ve.z) < 1e-12); // pinned
            vec3dp const vc = vel_at(cl, O_3dp);
            value_t const w = vc.y / (L / 2.0);
            CHECK(w == doctest::Approx(w_ref).epsilon(1e-12));
            CHECK(std::abs(vc.x) + std::abs(vc.z) < 1e-12);
            value_t const T1 = cl.system().total_energy();
            CHECK(T1 == doctest::Approx(0.5 * I_pin * w * w).epsilon(1e-12));
            CHECK(T1 / T0 == doctest::Approx(0.75).epsilon(1e-3));
            // impulse = momentum change = -Λ, componentwise (v_cm- = (0, -v, 0))
            CHECK(-Lam[0] == doctest::Approx(m * vc.x).epsilon(1e-12));
            CHECK(-Lam[1] == doctest::Approx(m * (vc.y + v)).epsilon(1e-12));
            CHECK(-Lam[2] == doctest::Approx(m * vc.z).epsilon(1e-12));
            CHECK(Lam[1] < 0.0); // pressed on the pin, like a loaded foot's λ_y
            for (int i = 0; i < 500; ++i)
                cl.step(1.0e-3);
            CHECK(std::abs(cl.system().total_energy() - T1) / T1 < 1e-9);
            fmt::println("  plastic: w = {:.4f}, T+/T- = {:.4f}, Λ = ({:.3f}, {:.3f}, "
                         "{:.3f})",
                         w, T1 / T0, Lam[0], Lam[1], Lam[2]);
        }
        { // restitution e = 1/2: the end rebounds at half its closing speed
            auto [cl, c] = build(-v);
            cl.activate_loop_with_impact(c, 0.5);
            vec3dp const ve = vel_at(cl, end_b);
            CHECK(ve.y == doctest::Approx(0.5 * v).epsilon(1e-12));
            CHECK(std::abs(ve.x) + std::abs(ve.z) < 1e-12);
            CHECK(cl.system().total_energy() < 0.5 * m * v * v);
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
        fmt::println("");
    }

    TEST_CASE("pga3dp: ground_contact3dp - touchdown, slope plane, flat foot, guard")
    {
        fmt::println("pga3dp: ground_contact3dp - unilateral contact by events (3D)");

        // The 2D chapter's rod in space: a slender cuboid (free body, cm at height h
        // above the floor z = 0, along x, at rest) falls with its -x end as the contact
        // point. Touchdown at the free-fall time to the bisection tolerance, the impact
        // of W2 pins the end (a ball joint, 3 rows) and the rod rotates about y.
        value_t const m = 2.0, L = 1.0, h_ = 0.02, g = 9.81, h = 0.5;
        auto const rod = make_cuboid_body(m, L, h_, h_);
        value_t const I_pin = m * (L * L + h_ * h_) / 12.0 + m * L * L / 4.0;
        vec3dp const end_b{-0.5, 0.0, 0.0, 1.0};
        auto const floor = ground_contact3dp::ground_plane(vec3dp{0.0, 0.0, 0.0, 1.0},
                                                           vec3dp{1.0, 0.0, 0.0, 1.0},
                                                           vec3dp{0.0, 1.0, 0.0, 1.0});
        // tilted about y by th0 < 0 so the -x end is the lower one (a rotation about +y
        // by th maps (x, 0, 0) to (x cos th, 0, -x sin th))
        value_t const th0 = -0.1;
        auto build = [&](value_t z0, value_t th = 0.0) {
            closed_loop_system3dp cl;
            cl.system().set_gravity(vec3dp{0.0, 0.0, -g, 0.0});
            cl.add_frame(static_frame3dp("W"));
            cl.add_body(static_frame3dp("rod", vec3dp{0.0, 0.0, z0, 1.0},
                                        vec3dp{0.0, th, 0.0, 0.0}),
                        rod);
            return cl;
        };
        {
            closed_loop_system3dp cl = build(h, th0);
            ground_contact3dp gc(cl, floor);
            CHECK(gc.height(vec3dp{3.0, -1.0, 0.7, 1.0}) ==
                  doctest::Approx(0.7).epsilon(1e-12));
            CHECK(gc.normal().z == doctest::Approx(1.0).epsilon(1e-12));
            size_t const c = gc.add({cl.index_of("rod"), end_b});
            size_t const c2 = gc.add({cl.index_of("rod"), vec3dp{0.5, 0.0, 0.0, 1.0}});
            value_t const h_end = gc.contact_height(c);
            CHECK(h_end == doctest::Approx(h - 0.5 * std::sin(-th0)).epsilon(1e-12));
            for (int i = 0; i < 400 && gc.events().empty(); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.events().size() == 1);
            auto const& ev = gc.events()[0];
            value_t const t_star = std::sqrt(2.0 * h_end / g);
            CHECK(ev.touchdown);
            CHECK(ev.t == doctest::Approx(t_star).epsilon(1e-9));
            CHECK(std::abs(ev.point.z) < 1e-12);
            CHECK(gc.active(c));
            CHECK(cl.constraint_rows() == 3);
            CHECK(std::abs(gc.contact_height(c)) < 1e-9);
            // the impact: -Lambda is the momentum change, v- = -g t* along z
            REQUIRE(ev.impulse.size() == 3);
            value_t const v = g * t_star;
            value_t const arm = 0.5 * L * std::cos(th0); // the cm's lever about the pin
            value_t const w_ref = -m * v * arm / I_pin;  // about y, sign by r x v
            value_t const v_plus = -v - ev.impulse[2] / m;
            value_t const w_ev = v_plus / arm;
            CHECK(std::abs(w_ev) == doctest::Approx(std::abs(w_ref)).epsilon(1e-9));
            fmt::println("  touchdown at t = {:.9f} (free fall {:.9f}), |w| from the "
                         "impulse = {:.6f} ({:.6f})",
                         ev.t, t_star, std::abs(w_ev), std::abs(w_ref));
            for (int i = 0; i < 900; ++i)
                gc.step(1.0e-3);
            CHECK(gc.active(c)); // the far end landed too: the rod lies on the floor
            CHECK(gc.active(c2));
            CHECK(std::abs(gc.contact_height(c2)) < 1e-9);
            CHECK(gc.normal_force(c) + gc.normal_force(c2) ==
                  doctest::Approx(m * g).epsilon(1e-6));
        }
        // a slope plane z = -0.2 x (normal up): the rod lands ON it
        {
            closed_loop_system3dp cl = build(h);
            auto const slope = ground_contact3dp::ground_plane(
                vec3dp{0.0, 0.0, 0.0, 1.0}, vec3dp{1.0, 0.0, -0.2, 1.0},
                vec3dp{0.0, 1.0, 0.0, 1.0});
            ground_contact3dp gc(cl, slope);
            CHECK(gc.normal().z > 0.0);
            size_t const c = gc.add({cl.index_of("rod"), end_b});
            for (int i = 0; i < 500 && gc.events().empty(); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            CHECK(std::abs(gc.contact_height(c)) < 1e-9);
            vec3dp const P = gc.contact_point(c);
            CHECK(P.z == doctest::Approx(-0.2 * P.x).epsilon(1e-9));
            fmt::println("  slope: landed at ({:.4f}, {:.4f}, {:.4f}), height {:.1e}",
                         P.x, P.y, P.z, gc.contact_height(c));
        }
        // a flat foot: welded on landing (6 rows), a cantilever along x -- the moment
        // has magnitude m g L/2 about y, the centre of pressure sits at x = L/2 in the
        // foot's axes: outside a 0.1 x 0.1 foot, tipping; as a ball joint it swings
        {
            closed_loop_system3dp cl = build(h, th0);
            ground_contact3dp gc(cl, floor);
            size_t const c =
                gc.add({cl.index_of("rod"), end_b, contact_kind3dp::flat, 0.1, 0.1});
            for (int i = 0; i < 500 && gc.events().empty(); ++i)
                gc.step(1.0e-3);
            REQUIRE(gc.active(c));
            CHECK(cl.constraint_rows() == 6);
            for (auto r : cl.system().joint_rates(cl.index_of("rod")))
                CHECK(std::abs(r) < 1e-9);
            CHECK(gc.normal_force(c) == doctest::Approx(m * g).epsilon(1e-9));
            vec3dp const M = gc.moment(c);
            CHECK(std::sqrt(M.x * M.x + M.y * M.y + M.z * M.z) ==
                  doctest::Approx(m * g * 0.5 * std::cos(th0)).epsilon(1e-9));
            CHECK(std::abs(gc.cop(c).x) ==
                  doctest::Approx(0.5 * std::cos(th0)).epsilon(1e-9));
            CHECK(std::abs(gc.cop(c).y) < 1e-9);
            CHECK(gc.tipping(c));
            fmt::println("  flat foot: N = {:.3f} N, |M| = {:.3f} N m, cop = ({:.3f}, "
                         "{:.3f}) -> tipping",
                         gc.normal_force(c), std::sqrt(M.x * M.x + M.y * M.y + M.z * M.z),
                         gc.cop(c).x, gc.cop(c).y);
            gc.set_kind(c, contact_kind3dp::point);
            CHECK(cl.constraint_rows() == 3);
            for (int i = 0; i < 100; ++i)
                gc.step(1.0e-3);
            auto const rates = cl.system().joint_rates(cl.index_of("rod"));
            value_t rot = 0.0;
            for (auto r : rates)
                rot = std::max(rot, std::abs(r));
            CHECK(rot > 0.1);
        }
        fmt::println("");
    }

    TEST_CASE("pga3dp: free-floating chain conserves momentum and energy (F)")
    {
        fmt::println("pga3dp: floating base - a space station: momentum and energy");

        // The decisive gate of Phase F. A free-floating assembly -- a heavy hub (the
        // free root joint) carrying two arms of two revolute links each, the arms set
        // swinging with opposite initial rates -- with no gravity and no external wrench
        // must conserve its total linear AND angular momentum, in the world frame, over a
        // long integration. An incorrect floating base leaks momentum through the root;
        // nothing tests it more sharply. Energy is conserved alongside.
        auto const hub = make_cuboid_body(40.0, 1.2, 0.8, 0.6);
        auto const link = make_cuboid_body(3.0, 0.8, 0.08, 0.08);
        vec3dp const ez{0.0, 0.0, 1.0, 0.0}, ey{0.0, 1.0, 0.0, 0.0};

        dynamic_system3dp sys;
        sys.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
        sys.add_frame(static_frame3dp("W"));
        sys.add_body(static_frame3dp("hub", vec3dp{0.0, 0.0, 0.0, 1.0}), hub,
                     kin_state3dp{.vel = vec3dp{0.1, 0.0, -0.05, 0.0},
                                  .omega = vec3dp{0.02, 0.3, -0.1, 0.0}});
        // arm A off +x: two links, hinges about z then y
        sys.add_revolute_body(static_frame3dp("A1", vec3dp{1.0, 0.0, 0.0, 1.0}), link,
                              vec3dp{-0.4, 0.0, 0.0, 1.0}, ez, 0.3, 1.5,
                              sys.index_of("hub"));
        sys.add_revolute_body(static_frame3dp("A2", vec3dp{0.8, 0.0, 0.0, 1.0}), link,
                              vec3dp{-0.4, 0.0, 0.0, 1.0}, ey, -0.5, -2.0,
                              sys.index_of("A1"));
        // arm B off -x
        sys.add_revolute_body(static_frame3dp("B1", vec3dp{-1.0, 0.0, 0.0, 1.0}), link,
                              vec3dp{0.4, 0.0, 0.0, 1.0}, ez, -0.2, -1.0,
                              sys.index_of("hub"));
        sys.add_revolute_body(static_frame3dp("B2", vec3dp{-0.8, 0.0, 0.0, 1.0}), link,
                              vec3dp{0.4, 0.0, 0.0, 1.0}, ey, 0.4, 1.2,
                              sys.index_of("B1"));
        CHECK(sys.dof_coords().size() == 10); // 6 + 4

        auto total_momentum = [&] {
            bivec3dp P{};
            for (size_t i = 1; i < sys.size(); ++i)
                P = P + sys.momentum_world(i);
            return P;
        };
        bivec3dp const P0 = total_momentum();
        value_t const E0 = sys.total_energy();
        value_t const P0n = std::sqrt(P0.vx * P0.vx + P0.vy * P0.vy + P0.vz * P0.vz +
                                      P0.mx * P0.mx + P0.my * P0.my + P0.mz * P0.mz);
        value_t max_dP = 0.0, max_dE = 0.0;
        for (int i = 0; i < 3000; ++i) { // 3 s of tumbling and swinging
            sys.step(1.0e-3);
            bivec3dp const d = total_momentum() - P0;
            max_dP = std::max(max_dP, std::sqrt(d.vx * d.vx + d.vy * d.vy + d.vz * d.vz +
                                                d.mx * d.mx + d.my * d.my + d.mz * d.mz));
            max_dE = std::max(max_dE, std::abs(sys.total_energy() - E0));
        }
        fmt::println("  |P| = {:.4f}: max |dP|/|P| = {:.2e}, dE/E = {:.2e}", P0n,
                     max_dP / P0n, max_dE / std::abs(E0));
        CHECK(max_dP / P0n < 1e-9);
        CHECK(max_dE / std::abs(E0) < 1e-9);

        // the hub actually moved and turned (the arms' motion reacts on it): no leak
        // means the hub's momentum changed while the total did not
        bivec3dp const Ph = sys.momentum_world(sys.index_of("hub"));
        bivec3dp const dh = Ph - P0;
        CHECK(std::sqrt(dh.vx * dh.vx + dh.vy * dh.vy + dh.vz * dh.vz) > 1e-3);
        fmt::println("");
    }

    TEST_CASE("pga3dp: motor joints agree with their composed 1-dof twins (F1)")
    {
        fmt::println(
            "pga3dp: spherical / cylindrical / planar vs. composed 1-dof chains");

        // Each multi-dof joint must move exactly like the composition of 1-dof joints
        // with massless carriers between them (the oracle the compose route provides):
        // same pendulum, same energy, same trajectory of the body's cm. The composed
        // chains carry the same initial rates on the corresponding screws.
        auto const bob = make_cuboid_body(2.0, 0.3, 0.2, 0.5);
        auto const ghost = make_cuboid_body(0.0, 0.05, 0.05, 0.05); // massless carrier
        vec3dp const ex{1.0, 0.0, 0.0, 0.0}, ey{0.0, 1.0, 0.0, 0.0},
            ez{0.0, 0.0, 1.0, 0.0};
        vec3dp const piv{0.0, 0.6, 0.0, 1.0}; // the pivot 0.6 above the body origin

        auto trajectory_error = [&](dynamic_system3dp& a, dynamic_system3dp& b, size_t fa,
                                    size_t fb, value_t dt = 1.0e-3) {
            value_t err = 0.0;
            for (int i = 0; i < int(std::round(0.8 / dt)); ++i) {
                a.step(dt);
                b.step(dt);
                vec3dp const pa = unitize(move3dp(O_3dp, a.get_pos_trafo(fa, 0)));
                vec3dp const pb = unitize(move3dp(O_3dp, b.get_pos_trafo(fb, 0)));
                err = std::max(err, std::abs(pa.x - pb.x) + std::abs(pa.y - pb.y) +
                                        std::abs(pa.z - pb.z));
            }
            return err;
        };

        // -- spherical: a ball-jointed pendulum vs. three revolutes (x, y, z through the
        //    pivot) with massless carriers. Rates: spin about the pivot's y and a swing.
        //    The two do NOT agree to round-off at one dt: the 3R chain integrates three
        //    angle coordinates whose RK4 truncation error is ~1e5 times that of the
        //    motor joint (measured: at dt = 1e-3 the chain moves 3.9e-7 under dt-halving,
        //    the motor joint 2.4e-12). So the gate is convergence: the difference must
        //    fall by 2^4 per halving, which says both approach the same trajectory
        {
            auto build = [&]() {
                dynamic_system3dp sph;
                sph.add_frame(static_frame3dp("W"));
                sph.add_spherical_body(
                    static_frame3dp("bob", vec3dp{0.0, -0.6, 0.0, 1.0}), bob, piv,
                    sph.index_of("W"));
                sph.set_joint_rates(sph.index_of("bob"), {0.7, 2.0, -0.4});

                dynamic_system3dp cmp;
                cmp.add_frame(static_frame3dp("W"));
                cmp.add_revolute_body(static_frame3dp("gx"), ghost, O_3dp, ex, 0.0, 0.7,
                                      cmp.index_of("W"));
                cmp.add_revolute_body(static_frame3dp("gy"), ghost, O_3dp, ey, 0.0, 2.0,
                                      cmp.index_of("gx"));
                cmp.add_revolute_body(static_frame3dp("bob", vec3dp{0.0, -0.6, 0.0, 1.0}),
                                      bob, piv, ez, 0.0, -0.4, cmp.index_of("gy"));
                return std::pair{sph, cmp};
            };
            auto [sph, cmp] = build();
            value_t const E1 = sph.total_energy(), E2 = cmp.total_energy();
            CHECK(E1 == doctest::Approx(E2).epsilon(1e-10));
            value_t const err1 =
                trajectory_error(sph, cmp, sph.index_of("bob"), cmp.index_of("bob"));
            CHECK(std::abs(sph.total_energy() - E1) / std::abs(E1) < 1e-8);
            auto [sph2, cmp2] = build();
            value_t const err2 = trajectory_error(sph2, cmp2, sph2.index_of("bob"),
                                                  cmp2.index_of("bob"), 5.0e-4);
            fmt::println("  spherical vs 3R: cm trajectory error {:.2e} (dt 1e-3) -> "
                         "{:.2e} (dt 5e-4), ratio {:.1f} (4th order: 16)",
                         err1, err2, err1 / err2);
            CHECK(err1 < 1e-6);
            CHECK(err1 / err2 == doctest::Approx(16.0).epsilon(0.1));
        }

        // -- cylindrical: rotation about + translation along y through the pivot, vs. a
        //    revolute (y) carrying a prismatic (y) with a massless carrier
        {
            dynamic_system3dp cyl;
            cyl.add_frame(static_frame3dp("W"));
            cyl.add_cylindrical_body(static_frame3dp("bob", vec3dp{0.0, -0.6, 0.0, 1.0}),
                                     bob, piv, ey, cyl.index_of("W"));
            cyl.set_joint_rates(cyl.index_of("bob"), {1.5, 0.3});

            dynamic_system3dp cmp;
            cmp.add_frame(static_frame3dp("W"));
            cmp.add_revolute_body(static_frame3dp("gy"), ghost, O_3dp, ey, 0.0, 1.5,
                                  cmp.index_of("W"));
            cmp.add_prismatic_body(static_frame3dp("bob", vec3dp{0.0, -0.6, 0.0, 1.0}),
                                   bob, ey, 0.0, 0.3, cmp.index_of("gy"));
            value_t const err =
                trajectory_error(cyl, cmp, cyl.index_of("bob"), cmp.index_of("bob"));
            fmt::println("  cylindrical vs R+P: cm trajectory error {:.2e}", err);
            CHECK(err < 1e-8);
        }

        // -- planar: the x-z plane through the pivot (normal y): vs. prismatic x,
        //    prismatic z and a revolute about y with massless carriers
        {
            dynamic_system3dp pl;
            pl.add_frame(static_frame3dp("W"));
            pl.add_planar_body(static_frame3dp("bob", vec3dp{0.0, -0.6, 0.0, 1.0}), bob,
                               piv, ey, pl.index_of("W"));
            // the planar screws are (t1, t2, rotation) with t1, t2 an orthonormal pair
            // perpendicular to y; read them back to set matching composed rates
            auto const& S = pl.joint_props(pl.index_of("bob")).screws;
            REQUIRE(S.size() == 3);
            pl.set_joint_rates(pl.index_of("bob"), {0.4, -0.2, 1.1});

            dynamic_system3dp cmp;
            cmp.add_frame(static_frame3dp("W"));
            cmp.add_prismatic_body(static_frame3dp("g1"), ghost,
                                   vec3dp{S[0].mx, S[0].my, S[0].mz, 0.0}, 0.0, 0.4,
                                   cmp.index_of("W"));
            cmp.add_prismatic_body(static_frame3dp("g2"), ghost,
                                   vec3dp{S[1].mx, S[1].my, S[1].mz, 0.0}, 0.0, -0.2,
                                   cmp.index_of("g1"));
            cmp.add_revolute_body(static_frame3dp("bob", vec3dp{0.0, -0.6, 0.0, 1.0}),
                                  bob, piv, ey, 0.0, 1.1, cmp.index_of("g2"));
            value_t const err =
                trajectory_error(pl, cmp, pl.index_of("bob"), cmp.index_of("bob"));
            fmt::println("  planar vs P+P+R: cm trajectory error {:.2e}", err);
            CHECK(err < 1e-8);
        }

        // -- helical: one coordinate turns and advances; pitch h means the body climbs
        //    h * q along the axis -- a screw joint driven at constant rate is a helix
        {
            dynamic_system3dp hel;
            hel.set_gravity(vec3dp{0.0, 0.0, 0.0, 0.0});
            hel.add_frame(static_frame3dp("W"));
            value_t const pitch = 0.05;
            // the screw axis is the WORLD z axis, given in the nut's own frame: the nut
            // sits at x = 0.3, so the axis passes through (-0.3, 0, 0) of the nut
            hel.add_helical_body(static_frame3dp("nut", vec3dp{0.3, 0.0, 0.0, 1.0}), bob,
                                 vec3dp{-0.3, 0.0, 0.0, 1.0}, ez, pitch, 0.0, 0.0,
                                 hel.index_of("W"));
            hel.set_driven_rate(hel.index_of("nut"), 2.0);
            for (int i = 0; i < 500; ++i)
                hel.step(1.0e-3); // q = 1.0 rad
            vec3dp const p =
                unitize(move3dp(O_3dp, hel.get_pos_trafo(hel.index_of("nut"), 0)));
            CHECK(p.z == doctest::Approx(pitch * 1.0).epsilon(1e-10)); // advanced h q
            CHECK(p.x == doctest::Approx(0.3 * std::cos(1.0)).epsilon(1e-10)); // turned q
            CHECK(p.y == doctest::Approx(0.3 * std::sin(1.0)).epsilon(1e-10));
            fmt::println("  helical: q = 1 rad -> z = {:.4f} (= pitch), turned by 1 rad",
                         p.z);
        }
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

    TEST_CASE("pga3dp: closed_loop_system3dp - an inactive constraint is absent (D1)")
    {
        fmt::println("pga3dp: closed_loop_system3dp - inactive constraint == absent");

        // a 3R arm whose hand is pinned to the world; system A carries a SECOND pin at
        // the elbow that is switched off, system B never had it. Every solve must agree
        // to the bit -- an inactive constraint contributes nothing -- and switching the
        // extra pin on must change the answer (so the flag is not a no-op).
        auto const seg = make_cuboid_body(1.5, 0.5, 0.05, 0.05);
        auto build = [&](bool extra) {
            closed_loop_system3dp cl;
            cl.add_frame(static_frame3dp("W"));
            vec3dp const ez{0.0, 0.0, 1.0, 0.0}, ey{0.0, 1.0, 0.0, 0.0};
            cl.add_revolute_body(static_frame3dp("L1", vec3dp{0.25, 0.0, 0.0, 1.0}), seg,
                                 vec3dp{-0.25, 0.0, 0.0, 1.0}, ez, 0.30, 0.0,
                                 cl.index_of("W"));
            cl.add_revolute_body(static_frame3dp("L2", vec3dp{0.5, 0.0, 0.0, 1.0}), seg,
                                 vec3dp{-0.25, 0.0, 0.0, 1.0}, ey, -0.40, 0.0,
                                 cl.index_of("L1"));
            cl.add_revolute_body(static_frame3dp("L3", vec3dp{0.5, 0.0, 0.0, 1.0}), seg,
                                 vec3dp{-0.25, 0.0, 0.0, 1.0}, ez, 0.50, 0.0,
                                 cl.index_of("L2"));
            size_t const L2 = cl.index_of("L2"), L3 = cl.index_of("L3"),
                         W = cl.index_of("W");
            vec3dp const hand = unitize(
                move3dp(vec3dp{0.25, 0.0, 0.0, 1.0}, cl.system().get_pos_trafo(L3, 0)));
            cl.add_loop_constraint(loop_constraint3dp{L3, vec3dp{0.25, 0.0, 0.0, 1.0}, W,
                                                      hand, constraint3dp::coincidence});
            if (extra) {
                vec3dp const elbow = unitize(move3dp(vec3dp{0.25, 0.0, 0.0, 1.0},
                                                     cl.system().get_pos_trafo(L2, 0)));
                size_t const c = cl.add_loop_constraint(
                    loop_constraint3dp{L2, vec3dp{0.25, 0.0, 0.0, 1.0}, W, elbow,
                                       constraint3dp::coincidence});
                cl.set_loop_active(c, false);
            }
            return cl;
        };

        auto A = build(true);
        auto B = build(false);
        CHECK(A.loop_count() == 2);
        CHECK(A.active_loop_count() == 1);
        CHECK(A.constraint_rows() == 3);
        CHECK(B.constraint_rows() == 3);

        std::vector<value_t> lA, lB;
        auto const qA = A.joint_accelerations(&lA);
        auto const qB = B.joint_accelerations(&lB);
        CHECK(lA.size() == 3);
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

        // switch the extra pin on: the elbow is now held too. In space the 3R arm with
        // its hand pinned is already immobile (3 rows on 3 joints), so q-ddot stays
        // zero -- what changes is the LOAD PATH: the elbow pin now carries part of it,
        // and the hand's multipliers change (the multipliers of a redundant closure
        // are the minimum-norm split)
        A.set_loop_active(1, true);
        CHECK(A.active_loop_count() == 2);
        CHECK(A.constraint_rows() == 6);
        A.assemble(); // the elbow has moved since t = 0; close the loop where it is now
        std::vector<value_t> lA2;
        auto const qA2 = A.joint_accelerations(&lA2);
        CHECK(lA2.size() == 6);
        for (size_t k = 0; k < qA2.size(); ++k)
            CHECK(std::abs(qA2[k]) < 1e-9); // still immobile
        bool differs = false;
        for (size_t k = 0; k < 3; ++k)
            differs = differs || (std::abs(lA2[k] - lA[k]) > 1e-9);
        CHECK(differs);

        fmt::println("  inactive: q-ddot / lambda / 100 steps bit-identical to the "
                     "system without the pin; active: {} rows, answer changes",
                     A.constraint_rows());
        fmt::println("");
    }

    TEST_CASE("pga3dp: closed_loop_system3dp - distance and frame kinds (D1)")
    {
        fmt::println("pga3dp: closed_loop_system3dp - distance (rod) and frame (weld)");

        value_t const g = 9.81;
        auto const seg = make_cuboid_body(2.0, 0.6, 0.06, 0.06); // 0.6 long along x

        // -- distance: a hinged bar whose tip is tied by a rod to a fixed point. One dof,
        //    one row: static. The multiplier is the rod tension, pinned against the
        //    moment balance about the hinge in GA -- moment_about(hinge, weight line +
        //    rod line) = 0 -- with the sign convention of the KKT layer: the constraint
        //    force on anchor a is -lambda * n, n the unit direction from b to a, so
        //    lambda > 0 is tension.
        {
            closed_loop_system3dp cl;
            cl.add_frame(static_frame3dp("W"));
            vec3dp const ez{0.0, 0.0, 1.0, 0.0};
            cl.add_revolute_body(static_frame3dp("bar", vec3dp{0.3, 0.0, 0.0, 1.0}), seg,
                                 vec3dp{-0.3, 0.0, 0.0, 1.0}, ez, -0.40, 0.0,
                                 cl.index_of("W"));
            size_t const bar = cl.index_of("bar"), W = cl.index_of("W");
            vec3dp const tip_b{0.3, 0.0, 0.0, 1.0};
            vec3dp const tip = unitize(move3dp(tip_b, cl.system().get_pos_trafo(bar, 0)));
            vec3dp const anchor{tip.x + 0.2, tip.y + 0.7, tip.z + 0.1, 1.0}; // rod end
            vec3dp const d{tip.x - anchor.x, tip.y - anchor.y, tip.z - anchor.z, 0.0};
            value_t const L = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
            cl.add_loop_constraint(
                loop_constraint3dp{bar, tip_b, W, anchor, constraint3dp::distance, L});
            CHECK(cl.residual_norm() < 1e-12);
            CHECK(cl.constraint_rows() == 1);

            std::vector<value_t> lambda;
            auto const qdd = cl.joint_accelerations(&lambda);
            CHECK(std::abs(qdd[0]) < 1e-9); // static
            CHECK(lambda.size() == 1);

            // GA statics about the hinge: the weight line at the cm and the rod line at
            // the tip must have zero total moment about the hinge axis (e_z here)
            vec3dp const hinge = unitize(
                move3dp(vec3dp{-0.3, 0.0, 0.0, 1.0}, cl.system().get_pos_trafo(bar, 0)));
            vec3dp const cm = unitize(move3dp(O_3dp, cl.system().get_pos_trafo(bar, 0)));
            vec3dp const n{d.x / L, d.y / L, d.z / L, 0.0};
            bivec3dp const Mw =
                moment_about(hinge, wdg(cm, vec3dp{0.0, -2.0 * g, 0.0, 0.0}));
            bivec3dp const Mr = moment_about(hinge, wdg(tip, n)); // per unit rod force
            // force on the bar's tip = -lambda n:  Mw.mz - lambda * Mr.mz = 0
            value_t const lambda_pred = Mw.mz / Mr.mz;
            fmt::println("  rod: lambda = {: .5f}  statics = {: .5f}", lambda[0],
                         lambda_pred);
            CHECK(lambda[0] == doctest::Approx(lambda_pred).epsilon(1e-9));
        }

        // -- frame: a 6R arm whose hand is WELDED to the world (6 rows, 6 dof: static).
        //    lambda is minus the wrench the wall exerts on the hand; it is pinned by the
        //    moment balance about every joint axis, in GA (moment_about).
        {
            closed_loop_system3dp cl;
            cl.add_frame(static_frame3dp("W"));
            vec3dp const ez{0.0, 0.0, 1.0, 0.0}, ey{0.0, 1.0, 0.0, 0.0};
            value_t const q0[6] = {0.3, -0.5, 0.7, 0.4, -0.6, 0.2};
            size_t parent = cl.index_of("W");
            for (int i = 0; i < 6; ++i) {
                cl.add_revolute_body(
                    static_frame3dp("L" + std::to_string(i), vec3dp{0.6, 0.0, 0.0, 1.0}),
                    seg, vec3dp{-0.3, 0.0, 0.0, 1.0}, (i % 2 == 0) ? ez : ey, q0[i], 0.0,
                    parent);
                parent = cl.index_of("L" + std::to_string(i));
            }
            size_t const hand_f = parent, W = cl.index_of("W");
            vec3dp const hand_b{0.3, 0.0, 0.0, 1.0};
            vec3dp const hand =
                unitize(move3dp(hand_b, cl.system().get_pos_trafo(hand_f, 0)));
            // the world anchor frame must carry the hand's CURRENT orientation for the
            // weld to close exactly: use the hand frame itself against a world-fixed
            // copy -- i.e. constrain the hand frame to a frame added at its pose
            cl.add_frame(static_frame3dp("wall"), W);
            cl.system().set_pose(
                cl.index_of("wall"),
                pose3dp_from_motor(cl.system().get_pos_trafo(hand_f, 0)));
            cl.add_loop_constraint(loop_constraint3dp{hand_f, hand_b, cl.index_of("wall"),
                                                      hand_b, constraint3dp::frame});
            CHECK(cl.residual_norm() < 1e-10);
            CHECK(cl.constraint_rows() == 6);

            std::vector<value_t> lambda;
            auto const qdd = cl.joint_accelerations(&lambda);
            CHECK(lambda.size() == 6);
            for (auto v : qdd)
                CHECK(std::abs(v) < 1e-8); // welded: static

            // Passive joints carry no torque, so about EVERY joint axis the moment of the
            // distal weight lines and of the wall wrench vanishes: axis . [moment about
            // the pivot] = 0. The wall acts on the hand with the force -lambda_F at the
            // hand point and the couple -lambda_M (KKT convention: the constraint force
            // is -G^T lambda). NOT the gravity resultant reduced at the hand -- the base
            // pin carries part of the load, so the weld reaction is fixed by the six
            // axis balances, not by the force sum.
            fmt::println(
                "  weld: lambda = ({: .4f}, {: .4f}, {: .4f} | {: .4f}, {: .4f}, "
                "{: .4f})",
                lambda[0], lambda[1], lambda[2], lambda[3], lambda[4], lambda[5]);
            vec3dp const lam_F{lambda[0], lambda[1], lambda[2], 0.0};
            for (int j = 0; j < 6; ++j) {
                auto const Mj =
                    cl.system().get_pos_trafo(cl.index_of("L" + std::to_string(j)), 0);
                vec3dp const pivot = unitize(move3dp(vec3dp{-0.3, 0.0, 0.0, 1.0}, Mj));
                vec3dp const axis = move3dp((j % 2 == 0) ? ez : ey, Mj); // world axis
                bivec3dp Wg{}; // the weight lines distal of joint j
                for (int i = j; i < 6; ++i) {
                    auto const Mi = cl.system().get_pos_trafo(
                        cl.index_of("L" + std::to_string(i)), 0);
                    Wg = Wg + wdg(unitize(move3dp(O_3dp, Mi)),
                                  vec3dp{0.0, -2.0 * g, 0.0, 0.0});
                }
                bivec3dp const Mtot =
                    moment_about(pivot, Wg) - moment_about(pivot, wdg(hand, lam_F));
                value_t const r = axis.x * (Mtot.mx - lambda[3]) +
                                  axis.y * (Mtot.my - lambda[4]) +
                                  axis.z * (Mtot.mz - lambda[5]);
                fmt::println("    joint {}: axis moment balance residual = {: .2e}", j,
                             r);
                CHECK(std::abs(r) < 1e-7);
            }
        }
        fmt::println("");
    }

    TEST_CASE("pga3dp: closed_loop_system3dp - two arms welded hand to hand (D1)")
    {
        fmt::println("pga3dp: closed_loop_system3dp - weld between two MOVING frames");

        // Both sides of the weld move (the wall-weld above has no joints on its b side):
        // arm A (3R) and arm B (3R, then 4R) welded hand to hand. Exact closure by
        // construction: B is built under its own base frame W2, then W2 is placed at
        // M_handA (x) rrev(M_handB_rel) so B's hand lands on A's hand, pose and all.
        value_t const g = 9.81;
        auto const seg = make_cuboid_body(2.0, 0.6, 0.06, 0.06);
        vec3dp const ez{0.0, 0.0, 1.0, 0.0}, ey{0.0, 1.0, 0.0, 0.0};
        vec3dp const hand_b{0.3, 0.0, 0.0, 1.0}, piv_b{-0.3, 0.0, 0.0, 1.0};
        vec3dp const w{0.0, -2.0 * g, 0.0, 0.0};

        auto build = [&](int nA, std::vector<value_t> const& qA, int nB,
                         std::vector<value_t> const& qB, std::vector<size_t>& A,
                         std::vector<size_t>& B) {
            closed_loop_system3dp cl;
            cl.add_frame(static_frame3dp("W"));
            size_t parent = cl.index_of("W");
            for (int i = 0; i < nA; ++i) {
                cl.add_revolute_body(
                    static_frame3dp("A" + std::to_string(i), vec3dp{0.6, 0.0, 0.0, 1.0}),
                    seg, piv_b, (i % 2 == 0) ? ez : ey, qA[i], 0.0, parent);
                parent = cl.index_of("A" + std::to_string(i));
                A.push_back(parent);
            }
            cl.add_frame(static_frame3dp("W2"), cl.index_of("W"));
            parent = cl.index_of("W2");
            for (int i = 0; i < nB; ++i) {
                cl.add_revolute_body(
                    static_frame3dp("B" + std::to_string(i), vec3dp{0.6, 0.0, 0.0, 1.0}),
                    seg, piv_b, (i % 2 == 0) ? ey : ez, qB[i], 0.0, parent);
                parent = cl.index_of("B" + std::to_string(i));
                B.push_back(parent);
            }
            auto const M_hA = cl.system().get_pos_trafo(A.back(), 0);
            auto const M_hB = cl.system().get_pos_trafo(B.back(), cl.index_of("W2"));
            cl.system().set_pose(cl.index_of("W2"),
                                 pose3dp_from_motor(rgpr(M_hA, rrev(M_hB))));
            cl.add_loop_constraint(loop_constraint3dp{A.back(), hand_b, B.back(), hand_b,
                                                      constraint3dp::frame});
            return cl;
        };

        // -- static: 3R + 3R = 6 dof, 6 rows. The interaction wrench is pinned by the
        //    moment balance about every joint axis of BOTH arms: on A's hand the force
        //    is -lambda_F and the couple -lambda_M, on B's hand +lambda_F, +lambda_M.
        {
            std::vector<size_t> A, B;
            auto cl = build(3, {0.4, -0.7, 0.5}, 3, {0.9, -0.3, 0.6}, A, B);
            CHECK(cl.residual_norm() < 1e-12);
            std::vector<value_t> lambda;
            auto const qdd = cl.joint_accelerations(&lambda);
            for (auto v : qdd)
                CHECK(std::abs(v) < 1e-8);
            vec3dp const hand =
                unitize(move3dp(hand_b, cl.system().get_pos_trafo(A.back(), 0)));
            vec3dp const lam_F{lambda[0], lambda[1], lambda[2], 0.0};
            auto balance = [&](std::vector<size_t> const& arm, bool a_side, value_t sgn) {
                for (size_t j = 0; j < arm.size(); ++j) {
                    auto const Mj = cl.system().get_pos_trafo(arm[j], 0);
                    vec3dp const pivot = unitize(move3dp(piv_b, Mj));
                    bool const even = (j % 2 == 0);
                    vec3dp const axis =
                        move3dp(a_side ? (even ? ez : ey) : (even ? ey : ez), Mj);
                    bivec3dp Wg{};
                    for (size_t i = j; i < arm.size(); ++i)
                        Wg = Wg + wdg(unitize(move3dp(
                                          O_3dp, cl.system().get_pos_trafo(arm[i], 0))),
                                      w);
                    bivec3dp const Mt = moment_about(pivot, Wg) +
                                        sgn * moment_about(pivot, wdg(hand, lam_F));
                    value_t const r = axis.x * (Mt.mx + sgn * lambda[3]) +
                                      axis.y * (Mt.my + sgn * lambda[4]) +
                                      axis.z * (Mt.mz + sgn * lambda[5]);
                    CHECK(std::abs(r) < 1e-7);
                }
            };
            balance(A, true, -1.0);
            balance(B, false, +1.0);
            fmt::println("  3R + 3R welded: static, both arms' axis balances hold");
        }

        // -- dynamic: 3R + 4R = 7 dof, 6 rows: one dof. Energy conserved and the weld
        //    holds through a swing.
        {
            std::vector<size_t> A, B;
            auto cl = build(3, {0.4, -0.7, 0.5}, 4, {0.9, -0.3, 0.6, 0.2}, A, B);
            CHECK(cl.residual_norm() < 1e-12);
            value_t const E0 = cl.system().total_energy();
            value_t max_dE = 0.0, max_g = 0.0;
            for (int i = 0; i < 1500; ++i) {
                cl.step(1.0e-3);
                max_dE = std::max(max_dE, std::abs(cl.system().total_energy() - E0));
                max_g = std::max(max_g, cl.residual_norm());
            }
            fmt::println("  3R + 4R welded: dE/E = {:.2e}, max |g| = {:.2e}",
                         max_dE / std::abs(E0), max_g);
            CHECK(max_dE / std::abs(E0) < 1e-6);
            CHECK(max_g < 1e-9);
        }
        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: closed_loop_system3dp")

TEST_SUITE("PGA3DP: coordinate transformation")
{

    TEST_CASE("pga3dp: coordinate transformation ECEF to ENU")
    {
        fmt::println("pga3dp: coordinate transformation ECEF to ENU");

        // ECEF: Earth-Centered, Earth-Fixed
        // -> cs rotates with earth, point fixed on earth surface has fixed coordinates
        // -> origin of cs in center of mass of earth (geocenter) at (0,0,0,1)
        // -> e1 axis (X) in equatorial plane through prime meridian (0° longitude)
        // -> e2 axis (Y) in equatorial plane 90° east of e1 (90° East longitude)
        // -> e3 axis (Z) aligns with the earth's mean rotational axis, pointing north
        //
        //
        // ENU: East-North-Up
        // -> cs local Cartesian system representing a each point on earth
        // -> e1 axis (X) points locally east
        // -> e2 axis (Y) points locally north
        // -> e3 axis (Z) points locally upwards

        // HINT 1: both systems are right handed systems
        // HINT 2: this also works directly in EGA3D if only directions are relevant

        /////////////////////////////////////////////////////////////////////////////////
        // the input form: what a user copies out of an encyclopedia entry
        /////////////////////////////////////////////////////////////////////////////////

        // the spellings a user might realistically paste in; each is only an EXAMPLE of a
        // notation, they are not all meant to name the same spot
        geo_pos_dms Berlin0 = {"52°31'12.0\"N", "13°24'36.0\"E", 35};
        geo_pos_dms Berlin1 = {"52°31'12\"N", "13°24'36\"E", 35};
        geo_pos_dms Berlin3 = {"52°31'N", "13°24'E", 35};
        geo_pos_dms Berlin4 = {"52°31' N", "13°24' E", 35};
        geo_pos_dms Berlin5 = {"52.51667N", "13.4E", 35};
        geo_pos_dms Berlin6 = {"52.15167N", "13.4E", 35};

        // nothing is parsed on construction -- the strings are kept as given, and
        // to_geo_pos() below is the one place the conversion happens
        CHECK(Berlin0.lat == "52°31'12.0\"N");
        CHECK(Berlin0.height == 35.0);

        // what the parser has to reproduce:
        //
        //     Berlin0,1  : 52°31'12" N, 13°24'36" E  (with seconds)  -> 52.52   , 13.41
        //     Berlin3,4,5: 52°31'    N, 13°24'    E  (no seconds)    -> 52.51666, 13.4
        //     Berlin6    : decimal degrees, approximate              -> 52.15167, 13.4
        //
        // dropping the 12" moves the position by ~370 m, so the first two groups differ;
        // Berlin6 is a deliberately approximate value and lands elsewhere again.

        CHECK(dms2deg(Berlin0.lat, geo_angle::latitude) ==
              doctest::Approx(52.52).epsilon(1e-12));
        CHECK(dms2deg(Berlin0.lon, geo_angle::longitude) ==
              doctest::Approx(13.41).epsilon(1e-12));

        CHECK(dms2deg(Berlin1.lat, geo_angle::latitude) == // seconds w/o decimal point
              dms2deg(Berlin0.lat, geo_angle::latitude));
        CHECK(dms2deg(Berlin1.lon, geo_angle::longitude) ==
              dms2deg(Berlin0.lon, geo_angle::longitude));

        CHECK(dms2deg(Berlin3.lat, geo_angle::latitude) ==
              doctest::Approx(52.0 + 31.0 / 60.0).epsilon(1e-12));
        CHECK(dms2deg(Berlin4.lat, geo_angle::latitude) == // blank before the letter
              dms2deg(Berlin3.lat, geo_angle::latitude));
        CHECK(dms2deg(Berlin4.lon, geo_angle::longitude) ==
              dms2deg(Berlin3.lon, geo_angle::longitude));

        // decimal degrees, no delimiter at all -- the same position as Berlin3/4 to the
        // five decimals given (~1 m)
        CHECK(dms2deg(Berlin5.lat, geo_angle::latitude) ==
              doctest::Approx(dms2deg(Berlin3.lat, geo_angle::latitude)).epsilon(1e-7));
        CHECK(dms2deg(Berlin6.lat, geo_angle::latitude) ==
              doctest::Approx(52.15167).epsilon(1e-12));
        CHECK(dms2deg(Berlin6.lon, geo_angle::longitude) ==
              doctest::Approx(13.4).epsilon(1e-12));

        // The German 'O' ("Ost") is deliberately NOT accepted -- only N/S/E/W. 'O' reads
        // as "Ouest" (= West) in French, so accepting it would silently flip the sign for
        // some users; the parser fails loudly instead of guessing.
        CHECK_THROWS_AS(dms2deg("13°24'36\"O", geo_angle::longitude),
                        std::invalid_argument);
        CHECK_THROWS_AS(dms2deg("13°24'O", geo_angle::longitude), std::invalid_argument);
        CHECK_THROWS_AS(dms2deg("13°24' O", geo_angle::longitude), std::invalid_argument);

        // signed decimal degrees are accepted as well, but not together with a
        // hemisphere letter; the hemisphere letter must match the angle asked for
        CHECK(dms2deg("-33°51'54\"", geo_angle::latitude) ==
              doctest::Approx(-(33.0 + 51.0 / 60.0 + 54.0 / 3600.0)).epsilon(1e-12));
        CHECK(dms2deg("151°12'36\"E", geo_angle::longitude) ==
              doctest::Approx(151.0 + 12.0 / 60.0 + 36.0 / 3600.0).epsilon(1e-12));
        CHECK(dms2rad("90°N", geo_angle::latitude) ==
              doctest::Approx(0.5 * pi).epsilon(1e-14));
        CHECK_THROWS_AS(dms2deg("-52°31'N", geo_angle::latitude), std::invalid_argument);
        CHECK_THROWS_AS(dms2deg("52°31'E", geo_angle::latitude), std::invalid_argument);
        CHECK_THROWS_AS(dms2deg("52°31'N", geo_angle::longitude), std::invalid_argument);
        CHECK_THROWS_AS(dms2deg("52°70'N", geo_angle::latitude), std::invalid_argument);
        CHECK_THROWS_AS(dms2deg("95°N", geo_angle::latitude), std::invalid_argument);
        CHECK_THROWS_AS(dms2deg("52'31°N", geo_angle::latitude), std::invalid_argument);
        CHECK_THROWS_AS(dms2deg("north", geo_angle::latitude), std::invalid_argument);

        /////////////////////////////////////////////////////////////////////////////////
        // the reference ellipsoid
        /////////////////////////////////////////////////////////////////////////////////

        // the radius runs from the equatorial to the polar radius, and is monotone in
        // between -- the two radii ARE the bases the shape is built on
        CHECK(wgs84.radius(0.0) == doctest::Approx(wgs84.r_equator).epsilon(1e-15));
        CHECK(wgs84.radius(0.5 * pi) == doctest::Approx(wgs84.r_pole).epsilon(1e-15));
        CHECK(wgs84.radius(deg2rad(45.0)) < wgs84.radius(deg2rad(30.0)));
        CHECK(wgs84.radius(deg2rad(60.0)) < wgs84.radius(deg2rad(45.0)));
        CHECK(wgs84.radius(deg2rad(-52.52)) == // symmetric about the equator
              doctest::Approx(wgs84.radius(deg2rad(52.52))).epsilon(1e-15));

        CHECK(wgs84.flattening() == doctest::Approx(1.0 / 298.257223563).epsilon(1e-12));
        CHECK(wgs84.e_sq() == doctest::Approx(0.00669437999014132).epsilon(1e-12));

        // geodetic vs geocentric latitude: equal at the equator and the pole, and
        // furthest apart in between (~11.5' at 45°)
        CHECK(wgs84.geocentric_lat(0.0) ==
              doctest::Approx(0.0).scale(1.0).epsilon(1e-14));
        CHECK(wgs84.geocentric_lat(0.5 * pi) == doctest::Approx(0.5 * pi).epsilon(1e-14));
        CHECK(rad2deg(deg2rad(45.0) - wgs84.geocentric_lat(deg2rad(45.0))) * 60.0 ==
              doctest::Approx(11.5).epsilon(2e-2)); // in arc minutes

        /////////////////////////////////////////////////////////////////////////////////
        // the calculation form
        /////////////////////////////////////////////////////////////////////////////////

        auto const B0 = to_geo_pos(Berlin0);

        CHECK(B0.lat == doctest::Approx(deg2rad(52.52)).epsilon(1e-14));
        CHECK(B0.lon == doctest::Approx(deg2rad(13.41)).epsilon(1e-14));
        CHECK(B0.height == 35.0); // passed through unchanged, along the normal

        // the "total radius" is derived, not stored: adding the height to the ellipsoid
        // radius is only exact ON the surface, because off it the height runs along the
        // normal rather than along the geocentric ray (0.18 mm apart already at 35 m)
        CHECK(distance_from_geocenter(B0) ==
              doctest::Approx(6364746.977026).epsilon(1e-13));
        CHECK(distance_from_geocenter(geo_pos{B0.lat, B0.lon, 0.0}) ==
              doctest::Approx(wgs84.radius(B0.lat)).epsilon(1e-13));
        CHECK(std::abs(distance_from_geocenter(B0) - (wgs84.radius(B0.lat) + B0.height)) <
              1e-3); // the naive sum, off by well under a millimetre at this height

        /////////////////////////////////////////////////////////////////////////////////
        // geodetic -> ECEF
        /////////////////////////////////////////////////////////////////////////////////

        // component-wise comparison with an EXPLICIT tolerance -- operator== compares
        // against an absolute eps, which is meaningless once coordinates are of the
        // order of the earth's radius (3.8e6 m: one ulp is already ~5e-10 m)
        //
        // NOT named `near`: on MSVC <windows.h> arrives via doctest's implementation
        // section and defines the legacy `near`/`far` memory-model keywords as EMPTY
        // macros, so `auto const near = ...` would expand to `auto const = ...`.
        auto const is_near = [](vec3dp const& a, vec3dp const& b, value_t tol) {
            return std::abs(a.x - b.x) < tol && std::abs(a.y - b.y) < tol &&
                   std::abs(a.z - b.z) < tol && std::abs(a.w - b.w) < tol;
        };

        // two positions where the ellipsoid fixes the answer analytically: on the
        // equator at the prime meridian the point is the equatorial radius, at the north
        // pole it is the polar radius
        auto const P_equ = geo_to_ecef(geo_pos{0.0, 0.0, 0.0});
        CHECK(is_near(P_equ, vec3dp{wgs84.r_equator, 0.0, 0.0, 1.0}, 1e-9));

        auto const P_pole = geo_to_ecef(geo_pos{0.5 * pi, 0.0, 0.0});
        CHECK(is_near(P_pole, vec3dp{0.0, 0.0, wgs84.r_pole, 1.0}, 1e-6));

        // Berlin0 against an independently computed reference [m] (pinned to the tenth
        // of a millimetre, i.e. ~1.5e-11 relative)
        auto const P_B0 = geo_to_ecef(B0);
        CHECK(P_B0.x == doctest::Approx(3783187.073822).epsilon(1e-11));
        CHECK(P_B0.y == doctest::Approx(901980.086455).epsilon(1e-11));
        CHECK(P_B0.z == doctest::Approx(5038246.874625).epsilon(1e-11));
        CHECK(P_B0.w == 1.0); // a point, not a direction

        CHECK(std::sqrt(P_B0.x * P_B0.x + P_B0.y * P_B0.y + P_B0.z * P_B0.z) ==
              doctest::Approx(distance_from_geocenter(B0)).epsilon(1e-13));

        // ecef_to_geo (Bowring) inverts it -- a genuine cross-check, since the two are
        // computed by different closed-form expressions
        for (auto const& p : {B0, to_geo_pos(Berlin6), geo_pos{0.0, 0.0, 0.0},
                              to_geo_pos(geo_pos_dms{"33°51'54\"S", "151°12'36\"E", 58}),
                              to_geo_pos(geo_pos_dms{"45°N", "179°54'W", 4000}),
                              geo_pos{deg2rad(89.999), deg2rad(-120.0), 4000.0}}) {
            auto const q = ecef_to_geo(geo_to_ecef(p));
            CHECK(q.lat == doctest::Approx(p.lat).scale(1.0).epsilon(1e-12));
            CHECK(q.lon == doctest::Approx(p.lon).scale(1.0).epsilon(1e-12));
            CHECK(q.height == doctest::Approx(p.height).scale(1.0).epsilon(1e-8));
        }

        /////////////////////////////////////////////////////////////////////////////////
        // the local ENU frame -- built from the geometry, and cross-checked by trig
        /////////////////////////////////////////////////////////////////////////////////

        // enu_at() derives the frame from the meridian plane N = wdg(e3, g) alone, with
        // g the ellipsoid gradient at the position: up = g, north = (g << N),
        // east = dual(N). Up is the ELLIPSOID NORMAL, i.e. the local vertical.
        auto const F = enu_at(B0);

        // the same frame written out from the angles -- an independent implementation,
        // so agreement between the two is a real cross-check, not a tautology
        auto const G = enu_basis_at(B0);

        CHECK(is_near(F.east, G.east, 1e-14));
        CHECK(is_near(F.north, G.north, 1e-14));
        CHECK(is_near(F.up, G.up, 1e-14));

        // and both against an independently computed reference
        CHECK(is_near(F.east, vec3dp{-0.231917681383160, 0.972735415753872, 0.0, 0.0},
                      1e-14));
        CHECK(is_near(
            F.north,
            vec3dp{-0.771929548448704, -0.184041937990537, 0.608484459368082, 0.0},
            1e-14));
        CHECK(is_near(
            F.up, vec3dp{0.591894383563182, 0.141118304974331, 0.793565789778978, 0.0},
            1e-14));

        // up is the ellipsoid NORMAL, not the ray from the geocenter: it is parallel to
        // the surface gradient 2*(x/a^2, y/a^2, z/b^2) at the foot point ...
        {
            auto const S = geo_to_ecef(geo_pos{B0.lat, B0.lon, 0.0});
            auto const aa = wgs84.r_equator * wgs84.r_equator;
            auto const bb = wgs84.r_pole * wgs84.r_pole;
            CHECK(is_congruent(vec3dp{S.x / aa, S.y / aa, S.z / bb, 0.0}, F.up));
        }

        // ... and it is measurably NOT the station ray -- 11.2' apart at this latitude,
        // which is the whole reason the normal is the one used in practice
        {
            auto const ray = vec3dp{P_B0.x, P_B0.y, P_B0.z, 0.0};
            auto const r = std::sqrt(ray.x * ray.x + ray.y * ray.y + ray.z * ray.z);
            auto const cos_a = (ray.x * F.up.x + ray.y * F.up.y + ray.z * F.up.z) / r;
            auto const arcmin = rad2deg(std::acos(std::min(1.0, cos_a))) * 60.0;
            CHECK(arcmin == doctest::Approx(11.2).epsilon(2e-2));
        }

        // taking the gradient AT the position instead of at its foot point is exact on
        // the surface, and drifts only slowly with height (see enu_at)
        {
            auto const p0 = geo_pos{B0.lat, B0.lon, 0.0};
            CHECK(is_near(enu_at(geo_to_ecef(p0)).up, enu_basis_at(p0).up, 1e-15));

            auto const p1 = geo_pos{B0.lat, B0.lon, 1000.0};
            CHECK(is_near(enu_at(geo_to_ecef(p1)).up, enu_basis_at(p1).up, 1e-6));
            CHECK(!is_near(enu_at(geo_to_ecef(p1)).up, enu_basis_at(p1).up, 1e-9));
        }

        // orthonormal and right-handed: east x north = up (the EGA3D cross product,
        // reachable here because only directions are involved -- HINT 2 above)
        auto const e3d = vec3d{F.east.x, F.east.y, F.east.z};
        auto const n3d = vec3d{F.north.x, F.north.y, F.north.z};
        auto const u3d = vec3d{F.up.x, F.up.y, F.up.z};
        CHECK(nrm(e3d) == doctest::Approx(1.0).epsilon(1e-14));
        CHECK(nrm(n3d) == doctest::Approx(1.0).epsilon(1e-14));
        CHECK(nrm(u3d) == doctest::Approx(1.0).epsilon(1e-14));
        CHECK(value_t(dot(e3d, n3d)) == doctest::Approx(0.0).scale(1.0).epsilon(1e-14));
        CHECK(value_t(dot(n3d, u3d)) == doctest::Approx(0.0).scale(1.0).epsilon(1e-14));
        CHECK(value_t(dot(u3d, e3d)) == doctest::Approx(0.0).scale(1.0).epsilon(1e-14));
        CHECK(cross(e3d, n3d) == u3d);

        // The frame must hold up over the whole latitude range, not just at the one
        // station pinned above. The degeneracy guards are the fragile part: the gradient
        // carries a 1/a^2 factor, so a guard written against an absolute eps starts
        // rejecting perfectly ordinary high-latitude positions long before the pole.
        // Swept on the surface, where the gradient IS the normal, so the two forms may be
        // held to full precision (off the surface they part by the documented drift --
        // checked separately above).
        for (int i = -89; i <= 89; ++i) {
            auto const p = geo_pos{deg2rad(double(i)), deg2rad(0.37 * double(i)), 0.0};
            auto const Fi = enu_at(geo_to_ecef(p));
            auto const Gi = enu_basis_at(p);
            CHECK(is_near(Fi.east, Gi.east, 1e-13));
            CHECK(is_near(Fi.north, Gi.north, 1e-13));
            CHECK(is_near(Fi.up, Gi.up, 1e-13));
        }

        // The MOTOR is what every later transformation runs on, and so far it has only
        // been exercised at one northern station. The frame sweep above cannot cover it:
        // a wrong sign in Rz(lon + 90°) or Rx(90° - lat) lives in enu_to_ecef_motor
        // alone. So sweep the motor against the frame over both hemispheres and all four
        // longitude quadrants -- if the two agree everywhere, no sign is left to guess.
        // the poles are included: the geo_pos overload has no degenerate case, since it
        // builds the meridian plane from the longitude rather than from the position
        for (int ilat = -90; ilat <= 90; ilat += 5) {
            for (int ilon = -180; ilon < 180; ilon += 15) {

                auto const p =
                    geo_pos{deg2rad(double(ilat)), deg2rad(double(ilon)), 250.0};
                auto const Mi = enu_to_ecef_motor(p);
                auto const Fi = enu_basis_at(p);

                CHECK(is_near(move3dp(x_dir_3dp, Mi), Fi.east, 1e-13));
                CHECK(is_near(move3dp(y_dir_3dp, Mi), Fi.north, 1e-13));
                CHECK(is_near(move3dp(z_dir_3dp, Mi), Fi.up, 1e-13));

                // the GA construction agrees with both, poles included
                auto const Gi = enu_at(p);
                CHECK(is_near(Gi.east, Fi.east, 1e-13));
                CHECK(is_near(Gi.north, Fi.north, 1e-13));
                CHECK(is_near(Gi.up, Fi.up, 1e-13));

                // the motor puts the ENU origin on the station ...
                CHECK(is_near(unitize(move3dp(O_3dp, Mi)), geo_to_ecef(p), 1e-6));
                // ... and its inverse brings the station back to the ENU origin
                CHECK(is_near(unitize(move3dp(geo_to_ecef(p), rrev(Mi))), O_3dp, 1e-6));
            }
        }

        // AT THE POLE only one axis is determined. Up is the polar axis, but every
        // horizontal direction is south (north pole), so the frame is free to rotate
        // about up -- and the right-hand rule cannot close that gap, since it fixes the
        // third axis only once a horizontal one is given. Which frame comes out is
        // therefore a matter of WHICH MERIDIAN the pole is approached along, and the two
        // overloads differ in what they can know about it:
        {
            // a geo_pos names a meridian, so the frame follows it continuously ...
            for (auto const lon_deg : {0.0, 45.0, 90.0, -90.0, 180.0}) {
                auto const at_pole = geo_pos{0.5 * pi, deg2rad(lon_deg), 0.0};
                auto const nearly = geo_pos{deg2rad(89.999999), deg2rad(lon_deg), 0.0};

                CHECK(is_near(enu_at(at_pole).east, enu_at(nearly).east, 1e-7));
                CHECK(is_near(enu_at(at_pole).north, enu_at(nearly).north, 1e-7));
                CHECK(is_near(enu_at(at_pole).up, vec3dp{0.0, 0.0, 1.0, 0.0}, 1e-13));

                // still an orthonormal right-handed triad, exactly at the pole
                auto const Fp = enu_at(at_pole);
                CHECK(cross(vec3d{Fp.east.x, Fp.east.y, Fp.east.z},
                            vec3d{Fp.north.x, Fp.north.y, Fp.north.z}) ==
                      vec3d{Fp.up.x, Fp.up.y, Fp.up.z});
            }

            // ... whereas a bare ECEF point on the axis carries no longitude, so the
            // prime meridian is taken by convention -- the same choice ecef_to_geo makes
            // when it reports lon == 0 there
            auto const F_np = enu_at(P_pole);
            CHECK(is_near(F_np.east, vec3dp{0.0, 1.0, 0.0, 0.0}, 1e-13));   // +e2
            CHECK(is_near(F_np.north, vec3dp{-1.0, 0.0, 0.0, 0.0}, 1e-13)); // -e1
            CHECK(is_near(F_np.up, vec3dp{0.0, 0.0, 1.0, 0.0}, 1e-13));     // +e3
            CHECK(is_near(F_np.east, enu_at(geo_pos{0.5 * pi, 0.0, 0.0}).east, 1e-13));

            auto const F_sp = enu_at(geo_to_ecef(geo_pos{-0.5 * pi, 0.0, 0.0}));
            CHECK(is_near(F_sp.east, vec3dp{0.0, 1.0, 0.0, 0.0}, 1e-13));  // +e2
            CHECK(is_near(F_sp.north, vec3dp{1.0, 0.0, 0.0, 0.0}, 1e-13)); // +e1
            CHECK(is_near(F_sp.up, vec3dp{0.0, 0.0, -1.0, 0.0}, 1e-13));   // -e3
        }

        // enu_motor_at(): the motor for a station given as an ECEF point. It must agree
        // with the geo_pos form everywhere -- including lon = 90 deg, where the tempting
        // "build the motor by aligning the reference frame onto the local one"
        // construction collapses (that rotation is by pi there, so its axis is
        // undetermined). See the note at enu_motor_at.
        for (int ilat = -80; ilat <= 80; ilat += 10) {
            for (int ilon = -180; ilon < 180; ilon += 15) {
                auto const p =
                    geo_pos{deg2rad(double(ilat)), deg2rad(double(ilon)), 250.0};
                auto const Mp = enu_motor_at(geo_to_ecef(p));
                CHECK(is_same_motion(Mp, enu_to_ecef_motor(p), 1e-9));
            }
        }

        // The poles are excluded above on purpose, and this is the reason: a point ON
        // the axis carries no longitude, so the point-based form can only produce the
        // prime-meridian frame -- the same convention enu_at(vec3dp) takes, and the same
        // one ecef_to_geo takes when it reports lon == 0 there. It therefore agrees with
        // the geo_pos form only for a station whose longitude IS zero.
        for (auto const lat : {0.5 * pi, -0.5 * pi}) {
            auto const at_zero = geo_pos{lat, 0.0, 250.0};
            CHECK(is_same_motion(enu_motor_at(geo_to_ecef(at_zero)),
                                 enu_to_ecef_motor(at_zero), 1e-9));

            auto const at_90 = geo_pos{lat, 0.5 * pi, 250.0};
            CHECK(!is_same_motion(enu_motor_at(geo_to_ecef(at_90)),
                                  enu_to_ecef_motor(at_90), 1e-9));
            // ... and what it does produce is the lon == 0 motor
            CHECK(is_same_motion(enu_motor_at(geo_to_ecef(at_90)),
                                 enu_to_ecef_motor(at_zero), 1e-9));
        }

        /////////////////////////////////////////////////////////////////////////////////
        // the height datum: a quoted elevation is orthometric, geo_pos is ellipsoidal
        /////////////////////////////////////////////////////////////////////////////////

        {
            geo_pos_dms const town = {"52°31'12\"N", "13°24'36\"E", 35};

            // an illustrative undulation for central Europe [m]; the real field is
            // measured and depends on longitude too, which is why it is supplied
            value_t const N_geoid = 45.0;
            auto const plain = to_geo_pos(town);
            auto const corrected = to_geo_pos(town, N_geoid);

            // h = H + N: only the height moves, the angles are untouched, so the local
            // frame's orientation is unaffected
            CHECK(corrected.height == doctest::Approx(35.0 + N_geoid).epsilon(1e-14));
            CHECK(corrected.lat == plain.lat);
            CHECK(corrected.lon == plain.lon);
            CHECK(is_near(enu_at(corrected).up, enu_at(plain).up, 1e-15));
            CHECK(is_near(enu_at(corrected).east, enu_at(plain).east, 1e-15));

            // and the position moves by exactly N along the local vertical -- which is
            // what "measured along the ellipsoid normal" means
            auto const P0 = geo_to_ecef(plain);
            auto const P1 = geo_to_ecef(corrected);
            auto const up = enu_at(plain).up;
            CHECK(is_near(P1,
                          vec3dp{P0.x + N_geoid * up.x, P0.y + N_geoid * up.y,
                                 P0.z + N_geoid * up.z, 1.0},
                          1e-6));

            // read back in the uncorrected station's own ENU frame it is straight up
            auto const q = unitize(move3dp(P1, rrev(enu_to_ecef_motor(plain))));
            CHECK(is_near(q, vec3dp{0.0, 0.0, N_geoid, 1.0}, 1e-6));

            // the default really is zero, i.e. the quoted number passes through
            CHECK(plain.height == 35.0);
        }

        // The geocenter is the one case that CANNOT be resolved by a convention: there
        // is no local vertical there either, so no axis of the frame is determined -- a
        // different situation from the pole, where one axis survives. It stays a throw.
        CHECK_THROWS_AS(enu_at(vec3dp{0.0, 0.0, 0.0, 1.0}), std::invalid_argument);

        /////////////////////////////////////////////////////////////////////////////////
        // the ECEF <-> ENU motor
        /////////////////////////////////////////////////////////////////////////////////

        // the motor must reproduce that frame: it takes the ENU origin to the station
        // point and the ENU axis directions to east/north/up
        auto const M = enu_to_ecef_motor(B0);

        // 1e-6 m = one micrometre, i.e. ~1.5e-13 relative on the earth's radius and
        // roughly 3 ulp of the coordinates involved
        CHECK(is_near(unitize(move3dp(O_3dp, M)), P_B0, 1e-6));
        CHECK(is_near(move3dp(x_dir_3dp, M), F.east, 1e-14));
        CHECK(is_near(move3dp(y_dir_3dp, M), F.north, 1e-14));
        CHECK(is_near(move3dp(z_dir_3dp, M), F.up, 1e-14));

        // and the inverse direction ECEF -> ENU is rrev(M)
        auto const M_inv = rrev(M);

        // the station is the ENU origin
        CHECK(is_near(unitize(move3dp(P_B0, M_inv)), O_3dp, 1e-6));

        // a point 100 m straight up over the station has ENU coordinates (0, 0, 100)
        auto const P_up = vec3dp{P_B0.x + 100.0 * F.up.x, P_B0.y + 100.0 * F.up.y,
                                 P_B0.z + 100.0 * F.up.z, 1.0};
        CHECK(is_near(unitize(move3dp(P_up, M_inv)), vec3dp{0.0, 0.0, 100.0, 1.0}, 1e-6));

        // the same station quoted 100 m higher sits at the same ENU spot -- this is the
        // height being measured along the ellipsoid normal, i.e. along ENU's up
        CHECK(
            is_near(unitize(move3dp(
                        geo_to_ecef(geo_pos{B0.lat, B0.lon, B0.height + 100.0}), M_inv)),
                    vec3dp{0.0, 0.0, 100.0, 1.0}, 1e-6));

        /////////////////////////////////////////////////////////////////////////////////
        // a second station seen from the first
        /////////////////////////////////////////////////////////////////////////////////

        // a rigid motion preserves distances: the ENU coordinates of a second station
        // must have the same length as the ECEF chord between the two stations
        {
            geo_pos_dms Madrid_dms = {"40°25'N", "3°43'W", 650};
            auto const Madrid = to_geo_pos(Madrid_dms);
            auto const P_M = geo_to_ecef(Madrid);
            auto const Q_M = unitize(move3dp(P_M, M_inv)); // Madrid in Berlin's ENU

            auto const dx = P_M.x - P_B0.x;
            auto const dy = P_M.y - P_B0.y;
            auto const dz = P_M.z - P_B0.z;
            auto const chord = std::sqrt(dx * dx + dy * dy + dz * dz);

            CHECK(std::sqrt(Q_M.x * Q_M.x + Q_M.y * Q_M.y + Q_M.z * Q_M.z) ==
                  doctest::Approx(chord).epsilon(1e-12));

            // Madrid is south-west of Berlin and, being far around the curve, well below
            // the local horizon plane
            CHECK(Q_M.x < 0.0); // west
            CHECK(Q_M.y < 0.0); // south
            CHECK(Q_M.z < 0.0); // below the local horizon plane

            fmt::println("   Madrid seen from Berlin in ENU: ({:.1f}, {:.1f}, {:.1f}) m,"
                         " chord = {:.1f} m",
                         Q_M.x, Q_M.y, Q_M.z, chord);
        }

        /////////////////////////////////////////////////////////////////////////////////
        // the types print
        /////////////////////////////////////////////////////////////////////////////////

        CHECK(fmt::format("{}", Berlin0) ==
              "geo_pos_dms(lat = 52°31'12.0\"N, lon = 13°24'36.0\"E, height = 35 m)");
        CHECK(fmt::format("{:.3f}", wgs84) ==
              "ellipsoid(r_equator = 6378137.000, r_pole = 6356752.314)");
        CHECK(fmt::format("{:.4f}", B0) ==
              "geo_pos(lat = 0.9166 rad, lon = 0.2340 rad, height = 35.0000 m)");
        CHECK(fmt::format("{:.3f}", F).starts_with("enu_frame(east = Vec3dp("));

        fmt::println("   {}", Berlin0);
        fmt::println("   {:.6f}", B0);
        fmt::println("   {:.6f}", F);

        fmt::println("");
    }


    /////////////////////////////////////////////////////////////////////////////////////
    // Satellites seen from the ground: ECEF and the two local ENU frames
    /////////////////////////////////////////////////////////////////////////////////////
    //
    // ORBIT MODEL (stated up front, because it fixes every number below).
    //
    // Each orbit is treated as a CIRCLE FIXED IN ECEF, sampled geometrically at 60 points
    // around it. The earth is not advanced between samples: the sweep walks the orbit's
    // geometry, it is not a time history. That is what makes the sweep useful for
    // catching sign errors -- it visits all 60 orientations of every orbit -- and it is
    // why the polar satellite stays in the lon = 0° plane throughout, as specified.
    //
    // The VELOCITY at each sampled point is nevertheless the physical one. A body on a
    // circular orbit of radius r runs at the mean motion n = rsqrt(mu/r^3) in the
    // INERTIAL frame, so its ECEF velocity carries the frame-rotation term:
    //
    //     v_inertial = n * (h x p)                h = orbit normal, p = position
    //     v_ecef     = v_inertial - omega x p     omega = earth rotation rate about e3
    //
    // Equivalently the whole configuration rotates rigidly in ECEF at
    //
    //     Omega_ecef = n * h - omega * e3
    //
    // which is the form the frame-tree block below hands to kinematic_system3dp.
    //
    // Two consequences worth watching, both asserted below:
    //
    // - the geostationary ring has n == omega, so v_ecef vanishes at EVERY sampled point
    // --
    //   the satellite stands still over the ground, which is the definition;
    // - the polar satellite's POSITION stays in the lon = 0° plane while its ECEF
    // VELOCITY
    //   does not: the earth turning underneath adds a westward component. A real polar
    //   orbit would precess out of that plane over time; here the earth is frozen, so the
    //   plane is the orbit's orientation at the sampled instant.
    //
    // The satellite body frame is as specified: x along the flight direction, z towards
    // the geocenter (nadir), y completing a right-handed triad, i.e. y = z x x. "Flight
    // direction" is the INERTIAL velocity -- for the geostationary satellite the ECEF
    // velocity is zero and would not define a direction at all.
    /////////////////////////////////////////////////////////////////////////////////////

    // WGS84 gravitational parameter and earth rotation rate
    inline value_t const mu_earth = 3.986004418e14;  // [m^3/s^2]
    inline value_t const omega_earth = 7.2921150e-5; // [rad/s]

    // a circular orbit as an earth-fixed circle: the plane is spanned by the orthonormal
    // pair (u0, u1) and the body runs from u0 towards u1
    struct circular_orbit {
        value_t r;    // radius from the geocenter [m]
        vec3d u0, u1; // orthonormal spanning pair, motion runs u0 -> u1
        value_t n;    // mean motion [rad/s]
    };

    inline circular_orbit make_orbit(value_t r, vec3d const& u0, vec3d const& u1)
    {
        return circular_orbit{r, u0, u1, std::sqrt(mu_earth / (r * r * r))};
    }

    // full state at orbit angle t, plus the body frame (x = flight, z = nadir, y = z x x)
    struct sat_state {
        vec3dp pos;        // ECEF position (w = 1)
        vec3dp v_ecef;     // ECEF velocity (w = 0)
        vec3dp v_in;       // inertial velocity (w = 0)
        vec3dp bx, by, bz; // body axes in ECEF (w = 0)
        vec3dp Omega;      // ECEF angular velocity of the whole configuration (w = 0)
    };

    inline vec3dp as_dir(vec3d const& d) { return vec3dp{d.x, d.y, d.z, 0.0}; }
    inline vec3d as_vec3d(vec3dp const& d) { return vec3d{d.x, d.y, d.z}; }

    inline sat_state state_at(circular_orbit const& o, value_t t)
    {
        auto const c = std::cos(t);
        auto const s = std::sin(t);

        auto const p = o.r * (c * o.u0 + s * o.u1);
        auto const vi = (o.n * o.r) * (-s * o.u0 + c * o.u1);

        // the frame-rotation term: omega x p with omega = omega_earth * e3
        auto const wxp = cross(omega_earth * e3_3d, p);

        auto const bz = normalize(-p); // nadir
        auto const bx = normalize(vi); // flight direction (inertial)
        auto const by = cross(bz, bx); // right-handed completion

        auto const h = cross(o.u0, o.u1);              // orbit normal
        auto const Om = o.n * h - omega_earth * e3_3d; // ECEF angular velocity

        return sat_state{vec3dp{p.x, p.y, p.z, 1.0},
                         as_dir(vi - wxp),
                         as_dir(vi),
                         as_dir(bx),
                         as_dir(by),
                         as_dir(bz),
                         as_dir(Om)};
    }

    // look angles of a satellite read straight off its ENU coordinates
    inline value_t elevation_deg(vec3dp const& q)
    {
        return rad2deg(std::asin(q.z / std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z)));
    }

    inline value_t azimuth_deg(vec3dp const& q)
    {
        auto const a = rad2deg(std::atan2(q.x, q.y)); // from north, towards east
        return (a < 0.0) ? a + 360.0 : a;
    }

    // The body->parent motor of a frame given by its origin and two of its axes: rotate
    // e1 onto ex, then turn about ex until e2 lands on ey, then translate out to the
    // origin. Both steps are ordinary get_motor() rotations about a line through the
    // origin.
    inline mvec3dp_e motor_from_frame(vec3dp const& origin, vec3d const& ex,
                                      vec3d const& ey)
    {
        auto const line_through_origin = [](vec3d const& d) {
            return wdg(O_3dp, vec3dp{d.x, d.y, d.z, 0.0});
        };

        // step 1: e1 -> ex
        mvec3dp_e R1 = I_3dp_mv_e;
        auto const c1 = cross(e1_3d, ex);
        if (nrm_sq(c1) > eps) {
            R1 = get_motor(line_through_origin(normalize(c1)), angle(e1_3d, ex));
        }
        else if (dot(e1_3d, ex) < 0.0) { // antiparallel: any perpendicular axis will do
            R1 = get_motor(line_through_origin(e3_3d), pi);
        }

        // step 2: turn about ex until the carried e2 lands on ey (signed angle about ex)
        auto const y1 = as_vec3d(move3dp(y_dir_3dp, R1));
        auto const R2 =
            get_motor(line_through_origin(ex),
                      std::atan2(value_t(dot(cross(y1, ey), ex)), value_t(dot(y1, ey))));

        auto const T = get_motor(vec3dp{origin.x, origin.y, origin.z, 0.0});
        return rgpr(T, rgpr(R2, R1)); // rotate first, then translate
    }

    // the three orbits, shared by both blocks below
    inline std::vector<std::pair<char const*, circular_orbit>> the_orbits()
    {
        value_t const r_leo = wgs84.r_equator + 360.0e3;

        // the geostationary radius follows from the definition n == omega_earth, giving a
        // height of ~35786 km ("about 36000 km")
        value_t const r_geo = std::cbrt(mu_earth / (omega_earth * omega_earth));

        return {
            // 1) equatorial, eastward, currently at lon 0: the circle runs e1 -> e2
            {"LEO equatorial 360 km", make_orbit(r_leo, e1_3d, e2_3d)},
            // 2) geostationary, same plane and sense, at the geostationary radius
            {"geostationary 35786 km", make_orbit(r_geo, e1_3d, e2_3d)},
            // 3) polar, in the plane containing lon 0° (the ECEF e1-e3 plane), running
            //    north over lon 0: the circle runs e1 -> e3
            {"LEO polar 360 km", make_orbit(r_leo, e1_3d, e3_3d)}};
    }

    inline size_t const n_samples = 60;


    TEST_CASE("pga3dp: satellites in ECEF and local ENU -- direct transformation")
    {
        fmt::println("pga3dp: satellites in ECEF and local ENU -- direct transformation");

        auto const Berlin = to_geo_pos(geo_pos_dms{"52°31'12\"N", "13°24'36\"E", 35});
        auto const Madrid = to_geo_pos(geo_pos_dms{"40°25'N", "3°43'W", 650});

        auto const M_B = enu_to_ecef_motor(Berlin); // Berlin ENU -> ECEF
        auto const M_M = enu_to_ecef_motor(Madrid); // Madrid ENU -> ECEF
        auto const B_inv = rrev(M_B);               // ECEF -> Berlin ENU
        auto const M_inv = rrev(M_M);               // ECEF -> Madrid ENU

        auto const P_B = geo_to_ecef(Berlin);
        auto const P_M = geo_to_ecef(Madrid);

        auto const orbits = the_orbits();

        value_t const r_leo = orbits[0].second.r;
        value_t const r_geo = orbits[1].second.r;

        CHECK(r_leo == doctest::Approx(6738137.0).epsilon(1e-15));
        CHECK(r_geo == doctest::Approx(42164172.931157).epsilon(1e-12));
        CHECK((r_geo - wgs84.r_equator) / 1000.0 ==
              doctest::Approx(35786.036).epsilon(1e-6)); // ~36000 km, as specified
        CHECK(orbits[1].second.n == doctest::Approx(omega_earth).epsilon(1e-12));
        CHECK(orbits[0].second.n * r_leo == doctest::Approx(7691.286152).epsilon(1e-10));

        // ---- the sweep: 60 points around each orbit -----------------------------------

        for (auto const& [name, orbit] : orbits) {

            for (size_t k = 0; k < n_samples; ++k) {

                value_t const t = 2.0 * pi * double(k) / double(n_samples);
                auto const S = state_at(orbit, t);

                // ---- the body frame is a right-handed orthonormal triad ---------------
                auto const bx = as_vec3d(S.bx), by = as_vec3d(S.by), bz = as_vec3d(S.bz);
                CHECK(nrm(bx) == doctest::Approx(1.0).epsilon(1e-14));
                CHECK(nrm(by) == doctest::Approx(1.0).epsilon(1e-14));
                CHECK(nrm(bz) == doctest::Approx(1.0).epsilon(1e-14));
                CHECK(value_t(dot(bx, by)) ==
                      doctest::Approx(0.0).scale(1.0).epsilon(1e-14));
                CHECK(value_t(dot(by, bz)) ==
                      doctest::Approx(0.0).scale(1.0).epsilon(1e-14));
                CHECK(value_t(dot(bz, bx)) ==
                      doctest::Approx(0.0).scale(1.0).epsilon(1e-14));
                CHECK(cross(bx, by) == bz); // right-handed

                // z points at the geocenter, x is perpendicular to the radius (circular)
                CHECK(is_congruent(S.bz, vec3dp{-S.pos.x, -S.pos.y, -S.pos.z, 0.0}));
                CHECK(value_t(dot(bx, vec3d{S.pos.x, S.pos.y, S.pos.z})) ==
                      doctest::Approx(0.0).scale(orbit.r).epsilon(1e-14));

                // the rigid ECEF rotation reproduces the velocity: v = Omega x p
                CHECK(value_t(bulk_nrm(S.v_ecef - as_dir(cross(as_vec3d(S.Omega),
                                                               vec3d{S.pos.x, S.pos.y,
                                                                     S.pos.z})))) < 1e-8);

                // ---- transformation into the two local frames -------------------------
                auto const qB = unitize(move3dp(S.pos, B_inv)); // position in Berlin ENU
                auto const qM = unitize(move3dp(S.pos, M_inv)); // position in Madrid ENU

                auto const vB = move3dp(S.v_ecef, B_inv); // velocity in Berlin ENU
                auto const vM = move3dp(S.v_ecef, M_inv); // velocity in Madrid ENU

                // ECEF and both ENU frames are earth-fixed and differ by a rigid motion,
                // so the speed is the SAME number in all three systems
                CHECK(value_t(bulk_nrm(vB)) ==
                      doctest::Approx(value_t(bulk_nrm(S.v_ecef)))
                          .scale(1.0)
                          .epsilon(1e-11));
                CHECK(value_t(bulk_nrm(vM)) ==
                      doctest::Approx(value_t(bulk_nrm(S.v_ecef)))
                          .scale(1.0)
                          .epsilon(1e-11));

                // ... and so is the distance from each station to the satellite
                auto const d_ecef = [&](vec3dp const& S0) {
                    return std::sqrt((S.pos.x - S0.x) * (S.pos.x - S0.x) +
                                     (S.pos.y - S0.y) * (S.pos.y - S0.y) +
                                     (S.pos.z - S0.z) * (S.pos.z - S0.z));
                };
                CHECK(value_t(bulk_nrm(qB)) ==
                      doctest::Approx(d_ecef(P_B)).epsilon(1e-11));
                CHECK(value_t(bulk_nrm(qM)) ==
                      doctest::Approx(d_ecef(P_M)).epsilon(1e-11));

                // the body axes carry over as directions (rotation only, length kept)
                CHECK(value_t(bulk_nrm(move3dp(S.bx, B_inv))) ==
                      doctest::Approx(1.0).epsilon(1e-13));
                CHECK(cross(as_vec3d(move3dp(S.bx, B_inv)),
                            as_vec3d(move3dp(S.by, B_inv))) ==
                      as_vec3d(move3dp(S.bz, B_inv))); // still right-handed locally

                // round trip back to ECEF
                CHECK(std::abs(unitize(move3dp(qB, M_B)).x - S.pos.x) < 1e-6);
                CHECK(std::abs(unitize(move3dp(qB, M_B)).z - S.pos.z) < 1e-6);
            }
        }

        // ---- per-orbit invariants that make the sweep worth running
        // --------------------

        for (size_t k = 0; k < n_samples; ++k) {

            value_t const t = 2.0 * pi * double(k) / double(n_samples);

            // 1) equatorial: stays in the equatorial plane, purely eastward ground track
            // at
            //    the ground-relative speed (n - omega) * r
            auto const S1 = state_at(orbits[0].second, t);
            CHECK(std::abs(S1.pos.z) < 1e-6);
            CHECK(std::abs(S1.v_ecef.z) < 1e-9);
            CHECK(value_t(bulk_nrm(S1.v_ecef)) ==
                  doctest::Approx((orbits[0].second.n - omega_earth) * r_leo)
                      .epsilon(1e-11));
            CHECK(value_t(bulk_nrm(S1.v_ecef)) ==
                  doctest::Approx(7199.933453).epsilon(1e-10));

            // 2) geostationary: motionless over the ground at EVERY sampled point -- the
            //    inertial orbital motion is exactly cancelled by the earth's rotation
            auto const S2 = state_at(orbits[1].second, t);
            CHECK(value_t(bulk_nrm(S2.v_ecef)) < 1e-8);
            CHECK(value_t(bulk_nrm(S2.v_in)) ==
                  doctest::Approx(3074.659979).epsilon(1e-10));
            CHECK(std::abs(S2.pos.z) < 1e-6); // equatorial

            // 3) polar: the POSITION stays in the lon 0° plane (y == 0) ...
            auto const S3 = state_at(orbits[2].second, t);
            CHECK(std::abs(S3.pos.y) < 1e-6);

            //    ... while the ECEF VELOCITY leaves it: the earth turning eastward under
            //    the orbit shows up as a westward component -omega * x
            CHECK(S3.v_ecef.y == doctest::Approx(-omega_earth * S3.pos.x).epsilon(1e-10));
            CHECK(value_t(bulk_nrm(S3.v_in)) ==
                  doctest::Approx(7691.286152).epsilon(1e-10));

            // over the poles the radius is parallel to the earth's axis, so the rotation
            // term vanishes and the ECEF speed equals the inertial one
            if (std::abs(S3.pos.x) < 1e-6) {
                CHECK(value_t(bulk_nrm(S3.v_ecef)) ==
                      doctest::Approx(value_t(bulk_nrm(S3.v_in))).epsilon(1e-11));
            }
        }

        // ---- the geostationary look angles, against an independent reference
        // ------------

        {
            auto const S = state_at(orbits[1].second, 0.0); // at lon 0°
            auto const qB = unitize(move3dp(S.pos, B_inv));
            auto const qM = unitize(move3dp(S.pos, M_inv));

            // classic spherical-earth formula for a geostationary look angle:
            //     el = atan((cos(gamma) - Re/r) / sin(gamma)),
            //     cos(gamma) = cos(lat) cos(lon_sat - lon_station)
            // It ignores the flattening, so it agrees only to a few hundredths of a
            // degree
            // -- which is precisely why it is an independent check on the ellipsoidal
            // result.
            auto const el_spherical = [&](geo_pos const& g) {
                auto const cg = std::cos(g.lat) * std::cos(0.0 - g.lon);
                auto const gamma = std::acos(cg);
                return rad2deg(
                    std::atan((cg - wgs84.r_equator / r_geo) / std::sin(gamma)));
            };

            CHECK(elevation_deg(qB) == doctest::Approx(28.6956).epsilon(1e-5));
            CHECK(azimuth_deg(qB) == doctest::Approx(196.7323).epsilon(1e-5));
            CHECK(elevation_deg(qB) ==
                  doctest::Approx(el_spherical(Berlin)).epsilon(2e-3));

            CHECK(elevation_deg(qM) == doctest::Approx(43.1320).epsilon(1e-5));
            CHECK(azimuth_deg(qM) == doctest::Approx(174.2741).epsilon(1e-5));
            CHECK(elevation_deg(qM) ==
                  doctest::Approx(el_spherical(Madrid)).epsilon(2e-3));

            // Both stations so far are NORTHERN, where a geostationary satellite always
            // sits to the south -- so a flipped north axis would show up merely as an
            // azimuth mirrored about due south, and both stations would agree on the lie.
            // Southern and equatorial stations break that symmetry: the same satellite is
            // then to the NORTH, and near the equator the north component is what almost
            // the whole azimuth is made of.
            {
                struct station {
                    char const* name;
                    geo_pos_dms pos;
                    value_t el, az; // independently computed reference [deg]
                    bool visible;   // is it above this station's horizon?
                };

                std::vector<station> const southern{
                    // southern, east of Greenwich: the satellite sits to the north-WEST
                    {"Cape Town",
                     {"33°55'S", "18°25'E", 25},
                     45.920953204,
                     329.150962372,
                     true},
                    // 0.2° south of the equator and 78° west, at 2850 m: almost due EAST
                    // and barely above the horizon
                    {"Quito",
                     {"0°13'S", "78°31'W", 2850},
                     2.789244143,
                     89.956208300,
                     true},
                    // 151° east: the satellite is on the far side of the earth
                    {"Sydney",
                     {"33°52'S", "151°13'E", 58},
                     -52.018583403,
                     224.564204018,
                     false}};

                for (auto const& st : southern) {
                    auto const g = to_geo_pos(st.pos);
                    auto const q = unitize(move3dp(S.pos, rrev(enu_to_ecef_motor(g))));

                    CHECK(elevation_deg(q) == doctest::Approx(st.el).epsilon(1e-8));
                    CHECK(azimuth_deg(q) == doctest::Approx(st.az).epsilon(1e-8));
                    CHECK((q.z > 0.0) == st.visible);

                    // The visibility test IS the scalar product of the local up with the
                    // line of sight -- expressing that line in ENU coordinates makes its
                    // third component exactly dot(up, los), because up IS the third ENU
                    // axis. So q.z already carries the test; computing the scalar
                    // product separately is the same number, as pinned here.
                    auto const P_st = geo_to_ecef(g);
                    auto const los =
                        vec3d{S.pos.x - P_st.x, S.pos.y - P_st.y, S.pos.z - P_st.z};
                    CHECK(value_t(dot(as_vec3d(enu_at(g).up), los)) ==
                          doctest::Approx(q.z).epsilon(1e-12));

                    // (this is the geometric horizon of the local tangent plane; a
                    // station at altitude actually sees slightly further, since its true
                    // horizon dips below horizontal by ~1.7° at Quito's 2850 m)

                    // For a southern station that can SEE it, the satellite lies to the
                    // north -- the sign simply unavailable from a northern station, where
                    // it is always to the south. Sydney is 151° of longitude away and
                    // would be looking through the earth, so its line of sight points
                    // south instead; that asymmetry is worth pinning too.
                    if (st.visible) {
                        CHECK(q.y > 0.0);
                    }
                    else {
                        CHECK(q.y < 0.0);
                    }

                    fmt::println("   {:<10} sees the geostationary satellite at"
                                 " el = {:7.3f}°, az = {:8.3f}°  ({})",
                                 st.name, elevation_deg(q), azimuth_deg(q),
                                 (q.z > 0.0) ? "above the horizon" : "BELOW the horizon");
                }

                // Cape Town looks north-west, Quito almost due east, Sydney cannot see it
                auto const q_cpt = unitize(
                    move3dp(S.pos, rrev(enu_to_ecef_motor(to_geo_pos(southern[0].pos)))));
                auto const q_uio = unitize(
                    move3dp(S.pos, rrev(enu_to_ecef_motor(to_geo_pos(southern[1].pos)))));
                auto const q_syd = unitize(
                    move3dp(S.pos, rrev(enu_to_ecef_motor(to_geo_pos(southern[2].pos)))));

                CHECK(q_cpt.x < 0.0); // west
                CHECK(q_cpt.z > 0.0); // above the horizon
                CHECK(q_uio.x > 0.0); // east
                CHECK(q_uio.z > 0.0);
                CHECK(std::abs(q_uio.y) < 0.01 * q_uio.x); // nearly on the equator
                CHECK(q_syd.z < 0.0); // below the horizon, as it must be
            }

            // Berlin (13.4°E) sees it south-west, Madrid (3.7°W) south-east and higher up
            CHECK(qB.x < 0.0);                            // west
            CHECK(qB.y < 0.0);                            // south
            CHECK(qM.x > 0.0);                            // east
            CHECK(qM.y < 0.0);                            // south
            CHECK(elevation_deg(qM) > elevation_deg(qB)); // Madrid is nearer the equator
        }

        // ---- report
        // --------------------------------------------------------------------

        fmt::println("");
        for (auto const& [name, orbit] : orbits) {
            fmt::println(
                "   {} : r = {:.3f} km, period = {:.1f} min, v_inertial = {:.3f} m/s",
                name, orbit.r / 1000.0, 2.0 * pi / orbit.n / 60.0, orbit.n * orbit.r);

            for (auto const t : {0.0, 0.5 * pi, pi, 1.5 * pi}) {
                auto const S = state_at(orbit, t);
                auto const qB = unitize(move3dp(S.pos, B_inv));
                auto const qM = unitize(move3dp(S.pos, M_inv));
                auto const vB = move3dp(S.v_ecef, B_inv);
                auto const vM = move3dp(S.v_ecef, M_inv);

                fmt::println("     u = {:5.1f}°", rad2deg(t));
                fmt::println("       ECEF   p = ({:11.1f},{:11.1f},{:11.1f}) km"
                             "  v = ({:9.3f},{:9.3f},{:9.3f}) |v| = {:8.3f} m/s",
                             S.pos.x / 1e3, S.pos.y / 1e3, S.pos.z / 1e3, S.v_ecef.x,
                             S.v_ecef.y, S.v_ecef.z, value_t(bulk_nrm(S.v_ecef)));
                fmt::println("       Berlin p = ({:11.1f},{:11.1f},{:11.1f}) km"
                             "  v = ({:9.3f},{:9.3f},{:9.3f}) |v| = {:8.3f} m/s",
                             qB.x / 1e3, qB.y / 1e3, qB.z / 1e3, vB.x, vB.y, vB.z,
                             value_t(bulk_nrm(vB)));
                fmt::println("       Madrid p = ({:11.1f},{:11.1f},{:11.1f}) km"
                             "  v = ({:9.3f},{:9.3f},{:9.3f}) |v| = {:8.3f} m/s",
                             qM.x / 1e3, qM.y / 1e3, qM.z / 1e3, vM.x, vM.y, vM.z,
                             value_t(bulk_nrm(vM)));
                if (qB.z > 0.0 || qM.z > 0.0) {
                    fmt::println(
                        "       above the horizon: Berlin el = {:6.2f}° az = {:6.2f}°"
                        "  |  Madrid el = {:6.2f}° az = {:6.2f}°",
                        elevation_deg(qB), azimuth_deg(qB), elevation_deg(qM),
                        azimuth_deg(qM));
                }
            }
            fmt::println("");
        }
    }


    TEST_CASE("pga3dp: satellites via static_system3dp / kinematic_system3dp")
    {
        fmt::println("pga3dp: satellites via static_system3dp / kinematic_system3dp");

        // The same three orbits, but nothing is transformed by hand. The scene is
        // declared as a frame tree hanging off ECEF,
        //
        //     ECEF (root)
        //      +-- Berlin   pose = the ENU motor of the station
        //      +-- Madrid   pose = the ENU motor of the station
        //      +-- Sat      pose = position + body attitude,  state = its ECEF motion
        //
        // and every quantity the direct block computed is then READ OFF the tree:
        // get_pos_trafo(a, b) is the transformation between any two of the frames, and
        // point_velocity() is the velocity of the satellite as a point of the moving
        // frame. The two blocks must agree to rounding -- that is the point of running
        // both.

        auto const Berlin = to_geo_pos(geo_pos_dms{"52°31'12\"N", "13°24'36\"E", 35});
        auto const Madrid = to_geo_pos(geo_pos_dms{"40°25'N", "3°43'W", 650});

        auto const M_B = enu_to_ecef_motor(Berlin);
        auto const M_M = enu_to_ecef_motor(Madrid);

        // the stations enter the tree as poses -- the ENU motor decoded into origin +
        // axis times angle, which is what a static_frame3dp carries
        auto const pose_B = pose3dp_from_motor(M_B);
        auto const pose_M = pose3dp_from_motor(M_M);

        kinematic_system3dp sys;
        sys.add_frame(static_frame3dp("ECEF")); // root, index 0
        sys.add_frame(static_frame3dp("Berlin", pose_B.origin, pose_B.rot), 0);
        sys.add_frame(static_frame3dp("Madrid", pose_M.origin, pose_M.rot), 0);
        sys.add_frame(static_frame3dp("Sat"), kin_state3dp{}, 0); // re-posed per sample

        size_t const i_sat = sys.index_of("Sat");

        // The tree reproduces the motors the direct block used. Two motors are compared
        // by what they DO, not with operator== : a motor's translation lives in
        // components of order 1e6 here, where an absolute-eps equality cannot resolve
        // anything. A rigid motion is fixed by where it sends the origin and the axes.
        // is_same_motion() now lives in the library (ga_pga3dp_ops.hpp) -- motors
        // double-cover the rigid motions, so operator== compares representations while
        // this compares motions.

        CHECK(is_same_motion(sys.get_pos_trafo("Berlin", "ECEF"), M_B));
        CHECK(is_same_motion(sys.get_pos_trafo("ECEF", "Berlin"), rrev(M_B)));
        CHECK(is_same_motion(sys.get_pos_trafo("Madrid", "ECEF"), M_M));

        // and the comparison is not vacuous -- it separates the two stations
        CHECK(!is_same_motion(sys.get_pos_trafo("Berlin", "ECEF"), M_M));

        // Berlin and Madrid seen from each other -- a frame-to-frame transformation that
        // never mentions ECEF at the call site
        {
            auto const M_BM = sys.get_pos_trafo("Berlin", "Madrid");
            auto const B_in_M =
                unitize(move3dp(O_3dp, M_BM)); // Berlin's origin, in Madrid
            auto const M_in_B =
                unitize(move3dp(O_3dp, sys.get_pos_trafo("Madrid", "Berlin")));

            CHECK(value_t(bulk_nrm(B_in_M)) ==
                  doctest::Approx(value_t(bulk_nrm(M_in_B))).epsilon(1e-12));
            CHECK(B_in_M.y > 0.0); // Berlin is north-east of Madrid
            CHECK(B_in_M.x > 0.0);
            CHECK(M_in_B.y < 0.0); // and Madrid south-west of Berlin
            CHECK(M_in_B.x < 0.0);

            fmt::println("");
            fmt::println("   Berlin seen from Madrid (ENU): ({:.1f}, {:.1f}, {:.1f}) km",
                         B_in_M.x / 1e3, B_in_M.y / 1e3, B_in_M.z / 1e3);
            fmt::println("   Madrid seen from Berlin (ENU): ({:.1f}, {:.1f}, {:.1f}) km",
                         M_in_B.x / 1e3, M_in_B.y / 1e3, M_in_B.z / 1e3);
        }

        auto const orbits = the_orbits();

        for (auto const& [name, orbit] : orbits) {

            for (size_t k = 0; k < n_samples; ++k) {

                value_t const t = 2.0 * pi * double(k) / double(n_samples);
                auto const S = state_at(orbit, t);

                // ---- put the satellite into the tree ----------------------------------
                //
                // pose: its ECEF position with the body axes as orientation
                auto const M_sat =
                    motor_from_frame(S.pos, as_vec3d(S.bx), as_vec3d(S.by));
                sys.set_pose(i_sat, pose3dp_from_motor(M_sat));

                // the motor really does carry the body frame
                CHECK(value_t(bulk_nrm(move3dp(x_dir_3dp, M_sat) - S.bx)) < 1e-12);
                CHECK(value_t(bulk_nrm(move3dp(y_dir_3dp, M_sat) - S.by)) < 1e-12);
                CHECK(value_t(bulk_nrm(move3dp(z_dir_3dp, M_sat) - S.bz)) < 1e-12);
                CHECK(value_t(bulk_nrm(unitize(move3dp(O_3dp, M_sat)) - S.pos)) < 1e-6);

                // state: the configuration rotates rigidly at Omega_ecef about the
                // geocenter, expressed in the satellite's own axes (a relative twist is
                // carried in body coordinates)
                auto const Om = as_vec3d(S.Omega);
                auto const bx = as_vec3d(S.bx), by = as_vec3d(S.by), bz = as_vec3d(S.bz);
                auto const v = as_vec3d(S.v_ecef);

                sys.set_state(
                    i_sat, kin_state3dp{vec3dp{value_t(dot(v, bx)), value_t(dot(v, by)),
                                               value_t(dot(v, bz)), 0.0},
                                        vec3dp{0.0, 0.0, 0.0, 0.0},
                                        vec3dp{value_t(dot(Om, bx)), value_t(dot(Om, by)),
                                               value_t(dot(Om, bz)), 0.0},
                                        vec3dp{0.0, 0.0, 0.0, 0.0}});

                // ---- read everything off the tree -------------------------------------

                auto const M_sat_to_B = sys.get_pos_trafo("Sat", "Berlin");
                auto const M_sat_to_M = sys.get_pos_trafo("Sat", "Madrid");

                // the satellite's position is its own origin, seen in the station's frame
                auto const qB_sys = unitize(move3dp(O_3dp, M_sat_to_B));
                auto const qM_sys = unitize(move3dp(O_3dp, M_sat_to_M));

                // the direct block's answer, for comparison
                auto const qB_dir = unitize(move3dp(S.pos, rrev(M_B)));
                auto const qM_dir = unitize(move3dp(S.pos, rrev(M_M)));

                CHECK(value_t(bulk_nrm(qB_sys - qB_dir)) < 1e-6);
                CHECK(value_t(bulk_nrm(qM_sys - qM_dir)) < 1e-6);

                // the body axes in the station frames agree too
                CHECK(value_t(bulk_nrm(move3dp(x_dir_3dp, M_sat_to_B) -
                                       move3dp(S.bx, rrev(M_B)))) < 1e-12);
                CHECK(value_t(bulk_nrm(move3dp(z_dir_3dp, M_sat_to_B) -
                                       move3dp(S.bz, rrev(M_B)))) < 1e-12);

                // ---- the velocity comes out of the kinematic layer
                // ---------------------
                //
                // point_velocity gives the world velocity of the satellite as a point of
                // the moving frame; it must reproduce v_ecef, earth-rotation term
                // included
                auto const v_sys = sys.point_velocity(S.pos, i_sat);
                CHECK(value_t(bulk_nrm(v_sys - S.v_ecef)) < 1e-8);

                // and in the station frames (a pure rotation, so the magnitude is kept)
                auto const vB_sys = move3dp(v_sys, rrev(M_B));
                auto const vM_sys = move3dp(v_sys, rrev(M_M));
                CHECK(value_t(bulk_nrm(vB_sys)) ==
                      doctest::Approx(value_t(bulk_nrm(S.v_ecef)))
                          .scale(1.0)
                          .epsilon(1e-9));
                CHECK(value_t(bulk_nrm(vM_sys)) ==
                      doctest::Approx(value_t(bulk_nrm(S.v_ecef)))
                          .scale(1.0)
                          .epsilon(1e-9));
                CHECK(value_t(bulk_nrm(vB_sys - move3dp(S.v_ecef, rrev(M_B)))) < 1e-8);

                // the geostationary satellite is at rest in the tree as well
                if (std::abs(orbit.n - omega_earth) < 1e-12) {
                    CHECK(value_t(bulk_nrm(v_sys)) < 1e-7);
                }
            }
        }

        // ---- report: the same table, produced entirely through the frame tree
        // -----------

        fmt::println("");
        for (auto const& [name, orbit] : orbits) {
            fmt::println("   {} (via the frame tree)", name);

            for (auto const t : {0.0, 0.5 * pi}) {
                auto const S = state_at(orbit, t);
                sys.set_pose(i_sat, pose3dp_from_motor(motor_from_frame(
                                        S.pos, as_vec3d(S.bx), as_vec3d(S.by))));

                auto const Om = as_vec3d(S.Omega);
                auto const bx = as_vec3d(S.bx), by = as_vec3d(S.by), bz = as_vec3d(S.bz);
                auto const v = as_vec3d(S.v_ecef);
                sys.set_state(
                    i_sat, kin_state3dp{vec3dp{value_t(dot(v, bx)), value_t(dot(v, by)),
                                               value_t(dot(v, bz)), 0.0},
                                        vec3dp{0.0, 0.0, 0.0, 0.0},
                                        vec3dp{value_t(dot(Om, bx)), value_t(dot(Om, by)),
                                               value_t(dot(Om, bz)), 0.0},
                                        vec3dp{0.0, 0.0, 0.0, 0.0}});

                auto const qB =
                    unitize(move3dp(O_3dp, sys.get_pos_trafo("Sat", "Berlin")));
                auto const qM =
                    unitize(move3dp(O_3dp, sys.get_pos_trafo("Sat", "Madrid")));
                auto const vw = sys.point_velocity(S.pos, i_sat);
                auto const vB = move3dp(vw, rrev(M_B));
                auto const vM = move3dp(vw, rrev(M_M));

                fmt::println("     u = {:5.1f}°", rad2deg(t));
                fmt::println("       Berlin p = ({:11.1f},{:11.1f},{:11.1f}) km"
                             "  v = ({:9.3f},{:9.3f},{:9.3f}) |v| = {:8.3f} m/s",
                             qB.x / 1e3, qB.y / 1e3, qB.z / 1e3, vB.x, vB.y, vB.z,
                             value_t(bulk_nrm(vB)));
                fmt::println("       Madrid p = ({:11.1f},{:11.1f},{:11.1f}) km"
                             "  v = ({:9.3f},{:9.3f},{:9.3f}) |v| = {:8.3f} m/s",
                             qM.x / 1e3, qM.y / 1e3, qM.z / 1e3, vM.x, vM.y, vM.z,
                             value_t(bulk_nrm(vM)));
            }
            fmt::println("");
        }
    }


    TEST_CASE("pga3dp: is_close and is_same_motion")
    {
        fmt::println("pga3dp: is_close and is_same_motion");

        /////////////////////////////////////////////////////////////////////////////////
        // is_close: relative equality, where operator== cannot resolve anything
        /////////////////////////////////////////////////////////////////////////////////

        // at coordinates of the order of an earth radius one ulp is already ~5e-10 m, so
        // operator== (absolute eps) rejects values that agree to every digit a double
        // carries, while is_close accepts them
        auto const big = vec3dp{3783187.073822, 901980.086455, 5038246.874625, 1.0};
        auto const big_1ulp =
            vec3dp{std::nextafter(big.x, 1e9), std::nextafter(big.y, 1e9),
                   std::nextafter(big.z, 1e9), 1.0};

        CHECK(big != big_1ulp);         // operator== says they differ ...
        CHECK(is_close(big, big_1ulp)); // ... is_close says they agree, correctly

        // and it is not merely permissive: a millimetre apart at that magnitude is a
        // relative 1.6e-10, far above the 1e-12 default, so it is still rejected
        auto const big_mm = vec3dp{big.x + 1.0e-3, big.y, big.z, 1.0};
        CHECK(!is_close(big, big_mm));
        CHECK(is_close(big, big_mm, 1e-9)); // ... unless the tolerance says otherwise

        // near zero the floor of 1 keeps the test absolute, so a component that should
        // vanish is still measured against rel_tol and not against nothing
        CHECK(is_close(vec3dp{0.0, 0.0, 0.0, 0.0}, vec3dp{1e-15, 0.0, 0.0, 0.0}));
        CHECK(!is_close(vec3dp{0.0, 0.0, 0.0, 0.0}, vec3dp{1e-6, 0.0, 0.0, 0.0}));

        // distinct from is_congruent, which allows an arbitrary scale factor
        auto const dir = vec3dp{1.0, 2.0, 3.0, 0.0};
        CHECK(is_congruent(dir, vec3dp{2.0, 4.0, 6.0, 0.0})); // same subspace ...
        CHECK(!is_close(dir, vec3dp{2.0, 4.0, 6.0, 0.0})); // ... but not the same value

        /////////////////////////////////////////////////////////////////////////////////
        // is_same_motion: motors double-cover the rigid motions
        /////////////////////////////////////////////////////////////////////////////////

        auto const M =
            enu_to_ecef_motor(to_geo_pos(geo_pos_dms{"52°31'12\"N", "13°24'36\"E", 35}));
        auto const M_neg = mvec3dp_e{-M};

        // THE point of the function: -M is the very same motion, and operator== does not
        // know it. Every object is moved identically by the two.
        CHECK(M != M_neg);               // representations differ ...
        CHECK(is_same_motion(M, M_neg)); // ... the motion does not
        CHECK(is_close(move3dp(x_dir_3dp, M), move3dp(x_dir_3dp, M_neg)));
        CHECK(is_close(unitize(move3dp(O_3dp, M)), unitize(move3dp(O_3dp, M_neg))));

        // a 2*pi turn is the identity motion but returns the NEGATED motor -- the
        // textbook case where comparing representations gives the wrong answer
        auto const R_2pi = get_motor(z_axis_3dp, 2.0 * pi);
        CHECK(is_same_motion(R_2pi, I_3dp_mv_e));
        CHECK(is_same_motion(get_motor(z_axis_3dp, 4.0 * pi), I_3dp_mv_e));

        // it still separates motions that really differ
        CHECK(!is_same_motion(M, get_motor(z_axis_3dp, deg2rad(1.0))));
        CHECK(!is_same_motion(get_motor(z_axis_3dp, deg2rad(30.0)),
                              get_motor(z_axis_3dp, deg2rad(-30.0))));
        CHECK(!is_same_motion(get_motor(vec3dp{1.0, 0.0, 0.0, 0.0}),
                              get_motor(vec3dp{0.0, 1.0, 0.0, 0.0})));

        // a SCALED motor is not the same motion: points land where they were (the
        // sandwich scales by |M|^2, which unitize divides out) but directions stretch
        auto const M_scaled = mvec3dp_e{2.0 * M};
        CHECK(is_close(unitize(move3dp(O_3dp, M)), unitize(move3dp(O_3dp, M_scaled))));
        CHECK(!is_same_motion(M, M_scaled));

        /////////////////////////////////////////////////////////////////////////////////
        // the same two properties in pga2dp
        /////////////////////////////////////////////////////////////////////////////////

        // translate by (3,4), then turn 35 degrees about the origin
        auto const R2 = get_motor(O_2dp, deg2rad(35.0));
        auto const M2 = rgpr(get_motor(vec2dp{3.0, 4.0, 0.0}), R2);
        auto const M2_neg = mvec2dp_u{-M2};

        CHECK(M2 != M2_neg);
        CHECK(is_same_motion(M2, M2_neg));
        CHECK(!is_same_motion(M2, R2));
        CHECK(is_same_motion(get_motor(O_2dp, 2.0 * pi), I_2dp_mv_u));

        CHECK(
            is_close(vec2dp{1e6, 2e6, 1.0}, vec2dp{std::nextafter(1e6, 1e9), 2e6, 1.0}));

        fmt::println("");
    }


    TEST_CASE("pga3dp: the line of sight, computed locally and compared across frames")
    {
        fmt::println("pga3dp: the line of sight, computed locally and compared "
                     "across frames");

        // The question this answers: "what is the line from my location to the
        // satellite?" -- and, more to the point, whether the line a user computes IN
        // THEIR OWN FRAME (the only thing they can measure directly) is the same object
        // as the line another frame computes. In PGA a line is a first-class element, so
        // the comparison is between two bivectors, not between six numbers that have to
        // be re-derived per frame.

        auto const Berlin = to_geo_pos(geo_pos_dms{"52°31'12\"N", "13°24'36\"E", 35});
        auto const Madrid = to_geo_pos(geo_pos_dms{"40°25'N", "3°43'W", 650});

        auto const M_B = enu_to_ecef_motor(Berlin);
        auto const M_M = enu_to_ecef_motor(Madrid);

        auto const P_B = geo_to_ecef(Berlin);
        auto const P_M = geo_to_ecef(Madrid);

        // the geostationary satellite of the previous case, at longitude 0
        value_t const r_geo = std::cbrt(mu_earth / (omega_earth * omega_earth));
        auto const P_S = vec3dp{r_geo, 0.0, 0.0, 1.0};

        /////////////////////////////////////////////////////////////////////////////////
        // the same line, built in two different frames
        /////////////////////////////////////////////////////////////////////////////////

        // in ECEF: join the two points
        auto const L_ecef = join(P_B, P_S);

        // in Berlin's OWN frame: the satellite as Berlin sees it, joined with Berlin's
        // own origin -- no knowledge of ECEF enters this line at all
        auto const S_in_B = unitize(move3dp(P_S, rrev(M_B)));
        auto const L_local = join(O_3dp, S_in_B);

        // carrying the local line out to ECEF must reproduce the ECEF line. join() does
        // not normalize, so the two agree up to weight -- which is what is_congruent
        // asks, and after unitize() they agree as values.
        CHECK(is_congruent(move3dp(L_local, M_B), L_ecef));
        CHECK(is_close(unitize(move3dp(L_local, M_B)), unitize(L_ecef), 1e-9));

        // both endpoints lie on the line: joining a point already on a line with that
        // line gives no plane. With the line unitized the residual is metres-scale
        // rounding on coordinates of order 4e7.
        // A line is compared to a LINE, not probed with an incidence test: join() was
        // what built L_ecef out of these very points, so join(join(A,B), A) is the
        // repeated-factor identity and vanishes for any A and B. It could only ever have
        // measured cancellation on the coordinate products, never the geometry.

        /////////////////////////////////////////////////////////////////////////////////
        // what the line tells the user, read in the local frame
        /////////////////////////////////////////////////////////////////////////////////

        // the line's ATTITUDE is the direction of sight; from the local origin it is
        // just the direction to the satellite
        auto const d_local = att(L_local);
        auto const F = enu_at(Berlin);

        CHECK(is_congruent(d_local, vec3dp{S_in_B.x, S_in_B.y, S_in_B.z, 0.0}));

        // the elevation is the angle the line makes with the local horizon, i.e. 90 deg
        // minus the angle against the local vertical -- and it must be the same number
        // the previous case pinned from the coordinates
        auto const dir = ega::normalize(vec3d{d_local.x, d_local.y, d_local.z});
        auto const up = vec3d{F.up.x, F.up.y, F.up.z}; // = e3 in Berlin's own frame
        value_t const elev =
            90.0 - rad2deg(std::acos(std::clamp(value_t(dot(dir, vec3d{0.0, 0.0, 1.0})),
                                                value_t(-1.0), value_t(1.0))));
        CHECK(elev == doctest::Approx(28.6956).epsilon(1e-5));
        CHECK(nrm(up) == doctest::Approx(1.0).epsilon(1e-14));

        /////////////////////////////////////////////////////////////////////////////////
        // the same line seen by somebody else
        /////////////////////////////////////////////////////////////////////////////////

        // Madrid can express BERLIN's line of sight in its own frame -- it needs only
        // the motor between the frames, not a re-derivation
        auto const L_in_M = move3dp(L_ecef, rrev(M_M));

        // and Madrid could equally build it from the two endpoints as IT sees them
        auto const B_in_M = unitize(move3dp(P_B, rrev(M_M)));
        auto const S_in_M = unitize(move3dp(P_S, rrev(M_M)));
        CHECK(is_close(unitize(L_in_M), unitize(join(B_in_M, S_in_M)), 1e-9));

        // Berlin's line is NOT Madrid's line of sight -- different lines to the same
        // satellite, meeting only at the satellite itself
        auto const L_M_own = join(P_M, P_S);
        CHECK(!is_congruent(L_ecef, L_M_own));

        // the comparison discriminates: Madrid does not lie on Berlin's line of sight,
        // and neither station's line is the other's
        CHECK(!is_congruent(L_ecef, join(P_M, P_S)));
        CHECK(!is_congruent(unitize(L_ecef), unitize(L_M_own)));

        /////////////////////////////////////////////////////////////////////////////////
        // and the same, obtained through the frame tree instead of by hand
        /////////////////////////////////////////////////////////////////////////////////

        {
            auto const pose_B = pose3dp_from_motor(M_B);
            auto const pose_M = pose3dp_from_motor(M_M);

            static_system3dp sys;
            sys.add_frame(static_frame3dp("ECEF"));
            sys.add_frame(static_frame3dp("Berlin", pose_B.origin, pose_B.rot), 0);
            sys.add_frame(static_frame3dp("Madrid", pose_M.origin, pose_M.rot), 0);

            // a line transforms exactly like a point does -- same motor, same call
            auto const L_tree = move3dp(L_local, sys.get_pos_trafo("Berlin", "ECEF"));
            CHECK(is_close(unitize(L_tree), unitize(L_ecef), 1e-9));

            auto const L_tree_M = move3dp(L_local, sys.get_pos_trafo("Berlin", "Madrid"));
            CHECK(is_close(unitize(L_tree_M), unitize(L_in_M), 1e-9));
        }

        fmt::println("   line of sight Berlin -> geostationary satellite");
        fmt::println("     in Berlin's own frame : att = ({:.6f}, {:.6f}, {:.6f}),"
                     " elevation = {:.3f}°",
                     dir.x, dir.y, dir.z, elev);
        fmt::println("     the same line in ECEF : {:.3e}", unitize(L_ecef));
        fmt::println("");
    }


    /////////////////////////////////////////////////////////////////////////////////////////
    // Paths on the ellipsoid: the planar curve, and the true geodesic
    /////////////////////////////////////////////////////////////////////////////////////////
    //
    // On a SPHERE the shortest path is a plane section -- the plane through the two
    // points and the centre cuts a great circle -- so the whole problem is incidence and
    // GA states it in one join. On an oblate ellipsoid that collapses: geodesics are not
    // plane curves at all. Both are computed here so the difference is visible rather
    // than asserted.
    //
    // THE PLANAR CURVE (great ellipse). Still pure GA: the plane through the geocenter
    // and the two points is join(join(A, B), O), and the curve in it is swept by a ROTOR
    // built from that plane's own bivector -- rexp(-B_hat t/2) carried the start
    // direction round to the end direction. Each swept direction is pushed out to the
    // surface along its ray. It is a genuine curve on the ellipsoid; it is simply not the
    // shortest one.
    //
    // THE GEODESIC. Writing the ellipsoid as the quadratic form G = diag(1/a^2, 1/a^2,
    // 1/b^2)
    // -- the surface is P.(G P) = 1 and its gradient is G P -- the geodesic condition is
    // that the tangent turns ONLY towards the surface normal, never within the surface:
    //
    //     dP/ds = T
    //     dT/ds = kappa * U,     U = normalize(G P),   kappa = -(T.(G T)) / |G P|
    //
    // No auxiliary sphere, no meridian/prime-vertical radii, no series: the ellipsoid
    // enters only through G, and U is the same gradient enu_at() already uses for "up".
    // On a sphere (a == b == R) the expression reduces to kappa = -1/R, the familiar
    // great-circle result. Integrated with the library's own rk4_step; the A -> B
    // boundary value problem is solved by shooting on the initial azimuth.
    //
    // Clairaut's relation r*sin(alpha) = const is the conserved quantity of the
    // ellipsoid's rotational symmetry -- the ONLY thing that survives of "the constant
    // plane" -- and is checked along the integrated path.
    /////////////////////////////////////////////////////////////////////////////////////////

    // the ellipsoid as a quadratic form: G v = (v.x/a^2, v.y/a^2, v.z/b^2)
    inline vec3d ell_G(vec3d const& v, ellipsoid const& el)
    {
        value_t const aa = el.r_equator * el.r_equator;
        value_t const bb = el.r_pole * el.r_pole;
        return vec3d{v.x / aa, v.y / aa, v.z / bb};
    }

    // push a direction from the geocenter out to the surface along its own ray
    inline vec3d ell_project(vec3d const& d, ellipsoid const& el)
    {
        return d * (1.0 / std::sqrt(value_t(dot(d, ell_G(d, el)))));
    }

    inline vec3d ell_normal(vec3d const& P, ellipsoid const& el)
    {
        return normalize(ell_G(P, el));
    }

    inline vec3d ecef3(vec3dp const& P) { return vec3d{P.x, P.y, P.z}; }

    // the azimuth of a surface direction T at surface point P [rad, from north towards
    // east]
    inline value_t azimuth_at(vec3d const& P, vec3d const& T)
    {
        auto const F = enu_at(vec3dp{P.x, P.y, P.z, 1.0});
        auto const e = vec3d{F.east.x, F.east.y, F.east.z};
        auto const n = vec3d{F.north.x, F.north.y, F.north.z};
        return std::atan2(value_t(dot(T, e)), value_t(dot(T, n)));
    }

    // wrap an angle in degrees into [0, 360)
    inline value_t deg_wrap360(value_t deg)
    {
        return std::fmod(std::fmod(deg, 360.0) + 360.0, 360.0);
    }

    // Clairaut's invariant r * sin(alpha) at a point of a surface path
    inline value_t clairaut(vec3d const& P, vec3d const& T)
    {
        return std::sqrt(P.x * P.x + P.y * P.y) * std::sin(azimuth_at(P, T));
    }

    // ---- the planar curve
    // -----------------------------------------------------------------

    // arc length of the great ellipse from A to B, swept by a rotor in the plane's
    // bivector
    inline value_t great_ellipse_length(vec3dp const& A, vec3dp const& B, size_t n_seg,
                                        ellipsoid const& el = wgs84)
    {
        auto const ua = normalize(ecef3(A));
        auto const ub = normalize(ecef3(B));

        auto const Bpl = wdg(ua, ub);                      // the plane, as a bivector
        auto const Bhat = Bpl * (1.0 / value_t(nrm(Bpl))); // its unit bivector
        value_t const theta = angle(ua, ub);

        auto prev = ell_project(ua, el);
        value_t L = 0.0;
        for (size_t i = 1; i <= n_seg; ++i) {
            value_t const t = theta * value_t(i) / value_t(n_seg);
            auto const d = rotate(ua, exp(-Bhat * (0.5 * t))); // the rotor sweep
            auto const P = ell_project(d, el);
            L += value_t(nrm(P - prev));
            prev = P;
        }
        return L;
    }

    // ---- the geodesic
    // ---------------------------------------------------------------------

    // one RK4 step of dP/ds = T, dT/ds = kappa U, followed by a projection back onto the
    // surface (the constraint is not enforced by the ODE, so it is restored each step)
    inline void geodesic_step(vec3d & P, vec3d & T, value_t ds, ellipsoid const& el)
    {
        std::vector<vec3d> u_mem{P, T};
        std::vector<vec3d> uh_mem(4);
        std::vector<vec3d> rhs_mem(2);

        auto u = mdspan<vec3d, dextents<size_t, 1>>(u_mem.data(), 2);
        auto uh = mdspan<vec3d, dextents<size_t, 2>>(uh_mem.data(), 2, 2);
        auto rhs = mdspan<vec3d const, dextents<size_t, 1>>(rhs_mem.data(), 2);

        for (size_t rk = 1; rk <= 4; ++rk) {
            auto const g = ell_G(u_mem[0], el);
            value_t const ng = value_t(nrm(g));
            rhs_mem[0] = u_mem[1];
            rhs_mem[1] = g * (-value_t(dot(u_mem[1], ell_G(u_mem[1], el))) / (ng * ng));
            rk4_step(u, uh, rhs, ds, rk);
        }

        P = ell_project(u_mem[0], el); // back onto the surface ...
        auto const U = ell_normal(P, el);
        T = normalize(u_mem[1] -
                      U * value_t(dot(u_mem[1], U))); // ... and tangent into it
    }

    struct geodesic_shot {
        value_t length; // arc length at closest approach, corrected along the tangent [m]
        value_t side;   // signed cross-track miss (left of the path positive) [m]
        value_t miss;   // closest approach distance [m]
        vec3d P_end;    // the point of closest approach
        vec3d T_end;    // the tangent there
        value_t clairaut_spread; // max - min of r*sin(alpha) along the path [m]
    };

    // integrate from A on the given initial azimuth and record the closest approach to B
    inline geodesic_shot geodesic_shoot(vec3d const& A, vec3d const& B, value_t az,
                                        value_t ds, ellipsoid const& el = wgs84)
    {
        auto const F = enu_at(vec3dp{A.x, A.y, A.z, 1.0});
        auto const e = vec3d{F.east.x, F.east.y, F.east.z};
        auto const n = vec3d{F.north.x, F.north.y, F.north.z};

        vec3d P = A;
        vec3d T = normalize(n * std::cos(az) + e * std::sin(az));

        value_t s = 0.0;
        value_t cl_min = clairaut(P, T), cl_max = cl_min;
        geodesic_shot best{0.0, 0.0, 1e30, P, T, 0.0};

        value_t const s_max = 4.0e7;
        while (s < s_max) {
            geodesic_step(P, T, ds, el);
            s += ds;

            value_t const cl = clairaut(P, T);
            cl_min = std::min(cl_min, cl);
            cl_max = std::max(cl_max, cl);

            value_t const d = value_t(nrm(P - B));
            if (d < best.miss) {
                best = geodesic_shot{s, 0.0, d, P, T, 0.0};
            }
            else if (d > best.miss + 5.0e3) {
                break; // past the closest approach
            }
        }

        auto const U = ell_normal(best.P_end, el);
        auto const left = cross(U, best.T_end);
        best.side = value_t(dot(B - best.P_end, left));
        // the last fraction of a step, along the tangent
        best.length += value_t(dot(B - best.P_end, best.T_end));
        best.clairaut_spread = cl_max - cl_min;
        return best;
    }

    // solve A -> B by bisecting the initial azimuth on the signed cross-track miss
    inline geodesic_shot geodesic_solve(vec3d const& A, vec3d const& B, value_t az_guess,
                                        value_t ds_final, ellipsoid const& el = wgs84)
    {
        value_t lo = az_guess - deg2rad(5.0);
        value_t hi = az_guess + deg2rad(5.0);
        value_t f_lo = geodesic_shoot(A, B, lo, 1000.0, el).side;

        for (int i = 0; i < 45; ++i) {
            value_t const mid = 0.5 * (lo + hi);
            value_t const f_mid = geodesic_shoot(A, B, mid, 1000.0, el).side;
            if ((f_mid < 0.0) == (f_lo < 0.0)) {
                lo = mid;
                f_lo = f_mid;
            }
            else {
                hi = mid;
            }
        }

        auto out = geodesic_shoot(A, B, 0.5 * (lo + hi), ds_final, el);
        out.side = 0.5 * (lo + hi); // report the converged azimuth here
        return out;
    }


    TEST_CASE("pga3dp: Berlin to Madrid -- the planar curve and the true geodesic")
    {
        fmt::println(
            "pga3dp: Berlin to Madrid -- the planar curve and the true geodesic");

        // both endpoints ON the ellipsoid (height 0): a path along the surface
        auto const A = ecef3(geo_to_ecef(
            geo_pos{to_geo_pos(geo_pos_dms{"52°31'12\"N", "13°24'36\"E", 0}).lat,
                    to_geo_pos(geo_pos_dms{"52°31'12\"N", "13°24'36\"E", 0}).lon, 0.0}));
        auto const Bm = ecef3(geo_to_ecef(
            geo_pos{to_geo_pos(geo_pos_dms{"40°25'N", "3°43'W", 0}).lat,
                    to_geo_pos(geo_pos_dms{"40°25'N", "3°43'W", 0}).lon, 0.0}));

        auto const Ap = vec3dp{A.x, A.y, A.z, 1.0};
        auto const Bp = vec3dp{Bm.x, Bm.y, Bm.z, 1.0};

        /////////////////////////////////////////////////////////////////////////////////
        // the rotor sweep really is a curve from A to B on the surface
        /////////////////////////////////////////////////////////////////////////////////

        {
            auto const ua = normalize(A);
            auto const ub = normalize(Bm);
            auto const Bpl = wdg(ua, ub);
            auto const Bhat = Bpl * (1.0 / value_t(nrm(Bpl)));
            value_t const theta = angle(ua, ub);

            // at t = theta the rotor has carried the start direction onto the end one
            CHECK(rotate(ua, exp(-Bhat * (0.5 * theta))) == ub);
            CHECK(rotate(ua, exp(-Bhat * 0.0)) == ua);

            // and the plane is the one PGA joins out of the two points and the geocenter
            auto const plane = join(join(Ap, Bp), O_3dp);
            auto const nrm_dir = vec3d{plane.x, plane.y, plane.z};
            CHECK(std::abs(value_t(dot(normalize(nrm_dir), ua))) < 1e-12);
            CHECK(std::abs(value_t(dot(normalize(nrm_dir), ub))) < 1e-12);

            // every swept point lies ON the ellipsoid: P.(G P) == 1
            for (int i = 0; i <= 10; ++i) {
                value_t const t = theta * value_t(i) / 10.0;
                auto const P = ell_project(rotate(ua, exp(-Bhat * (0.5 * t))), wgs84);
                CHECK(value_t(dot(P, ell_G(P, wgs84))) ==
                      doctest::Approx(1.0).epsilon(1e-14));
            }
        }

        // the length converges as the sweep is refined
        value_t const L_planar = great_ellipse_length(Ap, Bp, 200000);
        CHECK(L_planar ==
              doctest::Approx(great_ellipse_length(Ap, Bp, 100000)).epsilon(1e-10));
        CHECK(L_planar == doctest::Approx(1872384.247).epsilon(1e-8));

        /////////////////////////////////////////////////////////////////////////////////
        // the geodesic, integrated
        /////////////////////////////////////////////////////////////////////////////////

        auto const g = geodesic_solve(A, Bm, deg2rad(231.0), 25.0);
        value_t const az_start = deg_wrap360(rad2deg(g.side));

        // against Vincenty's inverse solution for WGS84
        CHECK(g.length == doctest::Approx(1872384.2).epsilon(1e-7));
        CHECK(az_start == doctest::Approx(230.8862).epsilon(1e-5));
        CHECK(g.miss < 30.0); // the sampled closest approach, before the tangent step

        value_t const az_end = deg_wrap360(rad2deg(azimuth_at(g.P_end, g.T_end)));
        CHECK(az_end == doctest::Approx(218.3559).epsilon(1e-5));

        // CLAIRAUT: the conserved quantity of the ellipsoid's rotational symmetry. It is
        // what survives of "the constant plane" once the sphere is left behind, and it
        // holds to a fraction of a micrometre over 1872 km of integration.
        CHECK(g.clairaut_spread < 1.0e-5);

        // the integrated path stays ON the surface (the projection each step is what
        // keeps it there -- the ODE alone does not enforce the constraint)
        CHECK(value_t(dot(g.P_end, ell_G(g.P_end, wgs84))) ==
              doctest::Approx(1.0).epsilon(1e-12));
        CHECK(std::abs(value_t(dot(g.T_end, ell_normal(g.P_end, wgs84)))) < 1e-12);

        /////////////////////////////////////////////////////////////////////////////////
        // the difference, which is the point of computing both
        /////////////////////////////////////////////////////////////////////////////////

        // the geodesic is SHORTER -- it must be, it is the shortest path
        CHECK(g.length < L_planar);
        CHECK(L_planar - g.length == doctest::Approx(0.036).epsilon(0.3));

        fmt::println("   Berlin -> Madrid");
        fmt::println("     great ellipse (planar, rotor-swept) : {:15.3f} m", L_planar);
        fmt::println("     geodesic      (integrated)          : {:15.3f} m", g.length);
        fmt::println("     planar longer by                    : {:15.3f} m",
                     L_planar - g.length);
        fmt::println("     azimuth {:.4f}° at Berlin, {:.4f}° arriving; Clairaut spread"
                     " {:.2e} m",
                     az_start, az_end, g.clairaut_spread);

        /////////////////////////////////////////////////////////////////////////////////
        // a route where the planar curve is visibly wrong
        /////////////////////////////////////////////////////////////////////////////////

        // Berlin -> Sydney: at 16000 km the two curves part by tens of metres, where at
        // Berlin -> Madrid they agree to a few centimetres. The planar approximation is
        // not "wrong at continental range" -- it is wrong at range.
        {
            auto const S = to_geo_pos(geo_pos_dms{"33°52'S", "151°13'E", 0});
            auto const Sy = ecef3(geo_to_ecef(geo_pos{S.lat, S.lon, 0.0}));
            auto const Sp = vec3dp{Sy.x, Sy.y, Sy.z, 1.0};

            value_t const L2 = great_ellipse_length(Ap, Sp, 200000);
            auto const g2 = geodesic_solve(A, Sy, deg2rad(78.0), 25.0);

            CHECK(g2.length == doctest::Approx(16088403.8).epsilon(1e-6));
            CHECK(g2.length < L2);
            CHECK(L2 - g2.length == doctest::Approx(24.28).epsilon(0.2));
            CHECK(g2.clairaut_spread < 1.0e-4);

            fmt::println("   Berlin -> Sydney");
            fmt::println("     great ellipse (planar, rotor-swept) : {:15.3f} m", L2);
            fmt::println("     geodesic      (integrated)          : {:15.3f} m",
                         g2.length);
            fmt::println("     planar longer by                    : {:15.3f} m",
                         L2 - g2.length);
        }

        fmt::println("");
    }

} // TEST_SUITE("PGA3DP: coordinate transformation")
