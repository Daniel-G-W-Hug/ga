#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Physics operations require hd_solver.hpp from the hd utility library.
// If not available, this header provides no functionality.
#if __has_include("hd/hd_solver.hpp")
#define HD_GA_PGA3DP_HAS_PHYSICS_OPS 1
#include "hd/hd_determinant.hpp" // includes hd_solver.hpp
#else
#define HD_GA_PGA3DP_HAS_PHYSICS_OPS 0
#endif

#include "ga_pga3dp_ops.hpp"
#include "ga_value_t.hpp" // for value_t used in convenience type alias

#include <array>
#include <mdspan>
#include <stdexcept>


namespace hd::ga::pga {

#if HD_GA_PGA3DP_HAS_PHYSICS_OPS

////////////////////////////////////////////////////////////////////////////////
// Inertia3dp: Inertia matrix for 3D projective GA (6x6 matrix)
//
// Used for rigid body dynamics in PGA3DP. The inertia map I[Omega] maps the
// rate of change Omega (a BiVec3dp) to momentum (also a BiVec3dp).
//
// From ga_docu/5_ga_physics_modelling.tex eq. 604-611:
// I_3D = m * [  0              Xz*Xw       -Xy*Xw        Xw^2     0        0       ]
//            [ -Xz*Xw          0            Xx*Xw        0        Xw^2     0       ]
//            [  Xy*Xw         -Xx*Xw        0            0        0        Xw^2    ]
//            [  Xy^2+Xz^2     -Xx*Xy       -Xx*Xz        0       -Xz*Xw    Xy*Xw   ]
//            [ -Xx*Xy          Xx^2+Xz^2   -Xy*Xz        Xz*Xw    0       -Xx*Xw   ]
//            [ -Xx*Xz         -Xy*Xz        Xx^2+Xy^2   -Xy*Xw    Xx*Xw    0       ]
//
// where X = (Xx, Xy, Xz, Xw) is the position vector with Xw as homogeneous coord.
//
// BiVec3dp components: (vx=e41, vy=e42, vz=e43, mx=e23, my=e31, mz=e12)
// The matrix maps: [Omega.vx, Omega.vy, Omega.vz, Omega.mx, Omega.my, Omega.mz]
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct Inertia3dp {
    std::array<T, 36> data{}; // row-major storage (6x6 matrix)

    // Default constructor (zero matrix)
    constexpr Inertia3dp() = default;

    // Accumulate inertia from another point mass
    constexpr Inertia3dp& operator+=(Inertia3dp const& other)
    {
        for (size_t i = 0; i < 36; ++i) {
            data[i] += other.data[i];
        }
        return *this;
    }

    // mdspan accessor for 2D indexing (mutable)
    auto view() { return std::mdspan<T, std::extents<size_t, 6, 6>>{data.data()}; }

    // mdspan accessor for 2D indexing (const)
    auto view() const
    {
        return std::mdspan<T const, std::extents<size_t, 6, 6>>{data.data()};
    }

    // Apply inertia map: I[Omega] = I * Omega
    // BiVec3dp components: vx=e41, vy=e42, vz=e43, mx=e23, my=e31, mz=e12
    BiVec3dp<T> operator()(BiVec3dp<T> const& Omega) const
    {
        auto I = view();
        return BiVec3dp<T>{
            I[0, 0] * Omega.vx + I[0, 1] * Omega.vy + I[0, 2] * Omega.vz +
                I[0, 3] * Omega.mx + I[0, 4] * Omega.my + I[0, 5] * Omega.mz,
            I[1, 0] * Omega.vx + I[1, 1] * Omega.vy + I[1, 2] * Omega.vz +
                I[1, 3] * Omega.mx + I[1, 4] * Omega.my + I[1, 5] * Omega.mz,
            I[2, 0] * Omega.vx + I[2, 1] * Omega.vy + I[2, 2] * Omega.vz +
                I[2, 3] * Omega.mx + I[2, 4] * Omega.my + I[2, 5] * Omega.mz,
            I[3, 0] * Omega.vx + I[3, 1] * Omega.vy + I[3, 2] * Omega.vz +
                I[3, 3] * Omega.mx + I[3, 4] * Omega.my + I[3, 5] * Omega.mz,
            I[4, 0] * Omega.vx + I[4, 1] * Omega.vy + I[4, 2] * Omega.vz +
                I[4, 3] * Omega.mx + I[4, 4] * Omega.my + I[4, 5] * Omega.mz,
            I[5, 0] * Omega.vx + I[5, 1] * Omega.vy + I[5, 2] * Omega.vz +
                I[5, 3] * Omega.mx + I[5, 4] * Omega.my + I[5, 5] * Omega.mz};
    }
};


// Create inertia matrix for a point mass at position X with mass m
// Pre: X should be unitized (X.w = 1 for finite points)
template <typename T>
    requires(std::floating_point<T>)
Inertia3dp<T> get_point_inertia(T m, Vec3dp<T> const& X)
{
    // Matrix from tex eq. 604-611
    // X components: x=Xx, y=Xy, z=Xz, w=Xw (homogeneous)
    T const Xx = X.x;
    T const Xy = X.y;
    T const Xz = X.z;
    T const Xw = X.w;

    Inertia3dp<T> I;
    auto v = I.view();

    // Row 0: [0, Xz*Xw, -Xy*Xw, Xw^2, 0, 0]
    v[0, 0] = T{0};
    v[0, 1] = m * Xz * Xw;
    v[0, 2] = -m * Xy * Xw;
    v[0, 3] = m * Xw * Xw;
    v[0, 4] = T{0};
    v[0, 5] = T{0};

    // Row 1: [-Xz*Xw, 0, Xx*Xw, 0, Xw^2, 0]
    v[1, 0] = -m * Xz * Xw;
    v[1, 1] = T{0};
    v[1, 2] = m * Xx * Xw;
    v[1, 3] = T{0};
    v[1, 4] = m * Xw * Xw;
    v[1, 5] = T{0};

    // Row 2: [Xy*Xw, -Xx*Xw, 0, 0, 0, Xw^2]
    v[2, 0] = m * Xy * Xw;
    v[2, 1] = -m * Xx * Xw;
    v[2, 2] = T{0};
    v[2, 3] = T{0};
    v[2, 4] = T{0};
    v[2, 5] = m * Xw * Xw;

    // Row 3: [Xy^2+Xz^2, -Xx*Xy, -Xx*Xz, 0, -Xz*Xw, Xy*Xw]
    v[3, 0] = m * (Xy * Xy + Xz * Xz);
    v[3, 1] = -m * Xx * Xy;
    v[3, 2] = -m * Xx * Xz;
    v[3, 3] = T{0};
    v[3, 4] = -m * Xz * Xw;
    v[3, 5] = m * Xy * Xw;

    // Row 4: [-Xx*Xy, Xx^2+Xz^2, -Xy*Xz, Xz*Xw, 0, -Xx*Xw]
    v[4, 0] = -m * Xx * Xy;
    v[4, 1] = m * (Xx * Xx + Xz * Xz);
    v[4, 2] = -m * Xy * Xz;
    v[4, 3] = m * Xz * Xw;
    v[4, 4] = T{0};
    v[4, 5] = -m * Xx * Xw;

    // Row 5: [-Xx*Xz, -Xy*Xz, Xx^2+Xy^2, -Xy*Xw, Xx*Xw, 0]
    v[5, 0] = -m * Xx * Xz;
    v[5, 1] = -m * Xy * Xz;
    v[5, 2] = m * (Xx * Xx + Xy * Xy);
    v[5, 3] = -m * Xy * Xw;
    v[5, 4] = m * Xx * Xw;
    v[5, 5] = T{0};

    return I;
}


// Get inverse of inertia matrix using LU decomposition
// Solves I * I_inv = Identity by back-substitution for each column
// Throws std::invalid_argument if inertia matrix is singular (det = 0)
template <typename T>
    requires(std::floating_point<T>)
Inertia3dp<T> get_inertia_inverse(Inertia3dp<T> const& I)
{
    // Check determinant before attempting inversion
    T const d = hd::det(I.view());
    if (d == T{0}) {
        throw std::invalid_argument(
            "get_inertia_inverse: singular inertia matrix (determinant is zero)");
    }

    // Copy data for LU decomposition (modifies in place)
    std::array<double, 36> A_data;
    for (size_t i = 0; i < 36; ++i) {
        A_data[i] = static_cast<double>(I.data[i]);
    }
    std::array<int, 6> perm;

    auto A = std::mdspan<double, std::dextents<size_t, 2>>{A_data.data(), 6, 6};
    auto perm_view = std::mdspan<int, std::dextents<size_t, 1>>{perm.data(), 6};

    hd::lu_decomp(A, perm_view);

    // Solve for each column of identity matrix to get inverse
    Inertia3dp<T> I_inv;
    auto A_const =
        std::mdspan<double const, std::dextents<size_t, 2>>{A_data.data(), 6, 6};
    auto perm_const = std::mdspan<int const, std::dextents<size_t, 1>>{perm.data(), 6};

    for (size_t col = 0; col < 6; ++col) {
        std::array<double, 6> e = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        e[col] = 1.0;
        auto e_view = std::mdspan<double, std::dextents<size_t, 1>>{e.data(), 6};
        hd::lu_backsubs(A_const, perm_const, e_view);

        // Store column in row-major format
        for (size_t row = 0; row < 6; ++row) {
            I_inv.data[row * 6 + col] = static_cast<T>(e[row]);
        }
    }
    return I_inv;
}


////////////////////////////////////////////////////////////////////////////////
// ODE right-hand side helpers for 3D rigid body dynamics
//
// From ga_docu/5_ga_physics_modelling.tex eq. 495-499:
//   Omega_dot = I_inv[ F - rcmt(Omega, I[Omega]) ]
//   M_dot     = 0.5 * M regressive_dot Omega
////////////////////////////////////////////////////////////////////////////////

// Compute Omega_dot = I_inv[ F - rcmt(Omega, I[Omega]) ]
// where F is the applied force/torque (as BiVec3dp), Omega is the current rate
// of change, I is the inertia matrix, and I_inv is its inverse
template <typename T>
    requires(std::floating_point<T>)
BiVec3dp<T> compute_omega_dot(Inertia3dp<T> const& I_inv, BiVec3dp<T> const& F,
                              BiVec3dp<T> const& Omega, Inertia3dp<T> const& I)
{
    // Omega_dot = I_inv[ F - rcmt(Omega, I[Omega]) ]
    BiVec3dp<T> I_Omega = I(Omega);
    BiVec3dp<T> rhs = F - rcmt(Omega, I_Omega);
    return I_inv(rhs);
}


// Compute M_dot = 0.5 * M ∨ Omega (motor derivative using regressive product)
// where M is the current motor and Omega is the rate of change
template <typename T>
    requires(std::floating_point<T>)
MVec3dp_E<T> compute_motor_dot(MVec3dp_E<T> const& M, BiVec3dp<T> const& Omega)
{
    // M_dot = 0.5 * M ∨ Omega (regressive geometric product)
    return T{0.5} * rgpr(M, Omega);
}

#endif // HD_GA_PGA3DP_HAS_PHYSICS_OPS

} // namespace hd::ga::pga

/////////////////////////////////////////////////////////////////////////////////////////
// user convenience type for pga3dp physics (based on value_t)
/////////////////////////////////////////////////////////////////////////////////////////
namespace hd::ga {

#if HD_GA_PGA3DP_HAS_PHYSICS_OPS
using inertia3dp = pga::Inertia3dp<value_t>; // 6x6 inertia matrix for 3D rigid body
#endif

} // namespace hd::ga
