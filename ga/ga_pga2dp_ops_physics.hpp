#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Physics operations require hd_solver.hpp from the hd utility library.
// If not available, this header provides no functionality.
#if __has_include("hd/hd_solver.hpp")
#define HD_GA_PGA2DP_HAS_PHYSICS_OPS 1
#include "hd/hd_determinant.hpp" // includes hd_solver.hpp
#else
#define HD_GA_PGA2DP_HAS_PHYSICS_OPS 0
#endif

#include "ga_pga2dp_ops.hpp"
#include "ga_value_t.hpp" // for value_t used in convenience type alias

#include <array>
#include <mdspan>
#include <stdexcept>


namespace hd::ga::pga {

#if HD_GA_PGA2DP_HAS_PHYSICS_OPS

////////////////////////////////////////////////////////////////////////////////
// Inertia2dp: Inertia matrix for 2D projective GA (3x3 matrix)
//
// Used for rigid body dynamics in PGA2DP. The inertia map I[Omega] maps the
// rate of change Omega (a Vec2dp) to the momentum (a BiVec2dp).
//
// From ga_docu/5_ga_physics_modelling.tex eq. 539-546:
// I_2D = m * [  0        Xz^2      -Xy*Xz    ]
//            [ -Xz^2      0         Xx*Xz    ]
//            [ -Xx*Xz    -Xy*Xz    Xx^2+Xy^2 ]
//
// where X = (Xx, Xy, Xz) is the position vector with Xz as homogeneous coord.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
struct Inertia2dp {
    std::array<T, 9> data{}; // row-major storage (3x3 matrix)

    // Default constructor (zero matrix)
    constexpr Inertia2dp() = default;

    // Accumulate inertia from another point mass
    constexpr Inertia2dp& operator+=(Inertia2dp const& other)
    {
        for (size_t i = 0; i < 9; ++i) {
            data[i] += other.data[i];
        }
        return *this;
    }

    // mdspan accessor for 2D indexing (mutable)
    auto view() { return std::mdspan<T, std::extents<size_t, 3, 3>>{data.data()}; }

    // mdspan accessor for 2D indexing (const)
    auto view() const
    {
        return std::mdspan<T const, std::extents<size_t, 3, 3>>{data.data()};
    }

    // Apply inertia map: I[Omega]
    // (map rate of change vector Omega to momentum bivector in 2D)
    BiVec2dp<T> operator()(Vec2dp<T> const& Omega) const
    {
        auto I = view();
        return BiVec2dp<T>{I[0, 0] * Omega.x + I[0, 1] * Omega.y + I[0, 2] * Omega.z,
                           I[1, 0] * Omega.x + I[1, 1] * Omega.y + I[1, 2] * Omega.z,
                           I[2, 0] * Omega.x + I[2, 1] * Omega.y + I[2, 2] * Omega.z};
    }

    // Apply inverse inertia map: I_inv[arg]
    // (map momentum bivector arg to vector Omega in 2D)
    Vec2dp<T> operator()(BiVec2dp<T> const& arg) const
    {
        auto I = view();
        return Vec2dp<T>{I[0, 0] * arg.x + I[0, 1] * arg.y + I[0, 2] * arg.z,
                         I[1, 0] * arg.x + I[1, 1] * arg.y + I[1, 2] * arg.z,
                         I[2, 0] * arg.x + I[2, 1] * arg.y + I[2, 2] * arg.z};
    }
};


// Create inertia matrix for a point mass at position X with mass m
// Pre: X should be unitized (X.z = 1 for finite points)
template <typename T>
    requires(numeric_type<T>)
Inertia2dp<T> get_point_inertia(T m, Vec2dp<T> const& X)
{
    // Matrix from tex eq. 539-546
    T const Xx = X.x;
    T const Xy = X.y;
    T const Xz = X.z;

    Inertia2dp<T> I;
    auto v = I.view();

    // Row 0: [0, Xz^2, -Xy*Xz]
    v[0, 0] = T{0};
    v[0, 1] = m * Xz * Xz;
    v[0, 2] = -m * Xy * Xz;

    // Row 1: [-Xz^2, 0, Xx*Xz]
    v[1, 0] = -m * Xz * Xz;
    v[1, 1] = T{0};
    v[1, 2] = m * Xx * Xz;

    // Row 2: [-Xx*Xz, -Xy*Xz, Xx^2+Xy^2]
    v[2, 0] = -m * Xx * Xz;
    v[2, 1] = -m * Xy * Xz;
    v[2, 2] = m * (Xx * Xx + Xy * Xy);

    return I;
}


// Get inverse of inertia matrix using LU decomposition
// Solves I * I_inv = Identity by back-substitution for each column
// Throws std::invalid_argument if inertia matrix is singular (det = 0)
template <typename T>
    requires(numeric_type<T>)
Inertia2dp<T> get_inertia_inverse(Inertia2dp<T> const& I)
{
    // Check determinant before attempting inversion
    T const d = hd::det(I.view());
    if (d == T{0}) {
        throw std::invalid_argument(
            "get_inertia_inverse: singular inertia matrix (determinant is zero)");
    }

    // Copy data for LU decomposition (modifies in place)
    std::array<double, 9> A_data;
    for (size_t i = 0; i < 9; ++i) {
        A_data[i] = static_cast<double>(I.data[i]);
    }
    std::array<int, 3> perm;

    auto A = std::mdspan<double, std::dextents<size_t, 2>>{A_data.data(), 3, 3};
    auto perm_view = std::mdspan<int, std::dextents<size_t, 1>>{perm.data(), 3};

    hd::lu_decomp(A, perm_view);

    // Solve for each column of identity matrix to get inverse
    Inertia2dp<T> I_inv;
    auto A_const =
        std::mdspan<double const, std::dextents<size_t, 2>>{A_data.data(), 3, 3};
    auto perm_const = std::mdspan<int const, std::dextents<size_t, 1>>{perm.data(), 3};

    for (size_t col = 0; col < 3; ++col) {
        std::array<double, 3> e = {0.0, 0.0, 0.0};
        e[col] = 1.0;
        auto e_view = std::mdspan<double, std::dextents<size_t, 1>>{e.data(), 3};
        hd::lu_backsubs(A_const, perm_const, e_view);

        // Store column in row-major format
        for (size_t row = 0; row < 3; ++row) {
            I_inv.data[row * 3 + col] = static_cast<T>(e[row]);
        }
    }
    return I_inv;
}


////////////////////////////////////////////////////////////////////////////////
// ODE right-hand side helpers for 2D rigid body dynamics
//
// From ga_docu/5_ga_physics_modelling.tex eq. 495-499:
//   Omega_dot = I_inv[ F - rcmt(Omega, I[Omega]) ]
//   M_dot     = 0.5 * M regressive_dot Omega
////////////////////////////////////////////////////////////////////////////////

// Compute Omega_dot = I_inv[ F - rcmt(Omega, I[Omega]) ]
// where F is the applied force/torque F = wdg(Q,f), Omega is the current rate of change,
// I is the inertia matrix, and I_inv is its inverse
template <typename T>
    requires(numeric_type<T>)
Vec2dp<T> compute_omega_dot(Inertia2dp<T> const& I_inv, BiVec2dp<T> const& F,
                            Vec2dp<T> const& Omega, Inertia2dp<T> const& I)
{
    // Omega_dot = I_inv[ F - rcmt(Omega, I[Omega]) ]
    BiVec2dp<T> I_Omega = I(Omega);
    BiVec2dp<T> rhs = F - rcmt(Omega, I_Omega);
    return I_inv(rhs); // returns the change rate Omega
}


// Compute M_dot = 0.5 * M ∨ Omega (motor derivative using regressive product)
// where M is the current motor and Omega is the rate of change
template <typename T>
    requires(numeric_type<T>)
MVec2dp_U<T> compute_motor_dot(MVec2dp_U<T> const& M, Vec2dp<T> const& Omega)
{
    // M_dot = 0.5 * M ∨ Omega (regressive geometric product)
    return T{0.5} * rgpr(M, Omega);
}

#endif // HD_GA_PGA2DP_HAS_PHYSICS_OPS

} // namespace hd::ga::pga

/////////////////////////////////////////////////////////////////////////////////////////
// user convenience type for pga2dp physics (based on value_t)
/////////////////////////////////////////////////////////////////////////////////////////
namespace hd::ga {

#if HD_GA_PGA2DP_HAS_PHYSICS_OPS
using inertia2dp = pga::Inertia2dp<value_t>; // 3x3 inertia matrix for 2D rigid body
#endif

} // namespace hd::ga
