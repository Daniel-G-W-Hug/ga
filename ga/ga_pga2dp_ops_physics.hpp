#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_solver.hpp" // hd::ga::lu_decomp / lu_backsubs / det
#include "ga_pga2dp_ops.hpp"
#include "ga_value_t.hpp" // for value_t used in convenience type alias

#include <array>
#include <mdspan>
#include <stdexcept>


namespace hd::ga::pga {

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

    // Subtract inertia contribution
    constexpr Inertia2dp& operator-=(Inertia2dp const& other)
    {
        for (size_t i = 0; i < 9; ++i) {
            data[i] -= other.data[i];
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


// Create inertia matrix for a uniform flat rectangular plate centered at the origin.
//
// The plate has width w along e1 (x-direction) and height h along e2 (y-direction).
// It is assumed to be centered at the origin in the body frame (Xz = 1).
//
// This is the exact continuous limit of accumulating point masses over the plate area.
// Derivation: integrate get_point_inertia over the area element dA = dx dy,
// using mean(x^2) = w^2/12, mean(y^2) = h^2/12, mean(x) = mean(y) = 0 at origin.
//
// Base result about cm / body origin O_b = (0,0,1):
//   I_cm = m * [  0    1         0        ]
//              [ -1    0         0        ]
//              [  0    0    (w^2+h^2)/12  ]
//
// I[0,1] =  m : total mass (Newton p = mv, Hodge-crossed e2 -> e23)
// I[1,0] = -m : total mass (Newton p = mv, Hodge-crossed e1 -> e31, opposite sign)
// I[2,2] = m*(w^2+h^2)/12 : moment of inertia about centroid (classical rectangle formula)
//
// Optional P_pivot parameter (default = body origin = cm):
// When P_pivot != O_b, the scalar parallel-axis (Steiner) correction is applied:
//   I_pivot[2,2] = J_cm + m*(Px² + Py²)
// All other entries remain equal to I_cm (the upper-left 2×2 mass block and the
// zero off-diagonal coupling terms are unchanged). This is the correct form for a
// body constrained to rotate about a fixed pivot Q_b in 2D: the coupled Euler ODE
// compute_omega_dot then yields α = τ_pivot / I_pivot[2,2] exactly.
template <typename T>
    requires(numeric_type<T>)
Inertia2dp<T> get_plate_inertia(T m, T w, T h,
                                Vec2dp<T> const& P_pivot = Vec2dp<T>{T{0}, T{0}, T{1}})
{
    Inertia2dp<T> I;
    auto v = I.view();

    // Row 0: [0, m, 0]
    v[0, 0] = T{0};
    v[0, 1] = m;
    v[0, 2] = T{0};

    // Row 1: [-m, 0, 0]
    v[1, 0] = -m;
    v[1, 1] = T{0};
    v[1, 2] = T{0};

    // Row 2: [0, 0, m*(w^2+h^2)/12]
    v[2, 0] = T{0};
    v[2, 1] = T{0};
    v[2, 2] = m * (w * w + h * h) / T{12};

    // Apply parallel-axis (Steiner) correction if pivot differs from body origin (cm).
    // For a body constrained to rotate about a fixed pivot Q_b (body frame),
    // only the scalar moment I[2,2] gains the Steiner term m*r²:
    //   I[2,2] = J_cm + m*(Px² + Py²)
    // The off-diagonal coupling terms (I[0,2], I[1,2] etc.) must NOT be added here.
    // Adding them via get_point_inertia(m, P_pivot) would make I[Omega] = (0,0,J_cm*ω)
    // instead of (0,0,I_pivot_zz*ω), causing I_inv[2,2] = 1/J_cm instead of
    // 1/I_pivot_zz — a factor of ~4 error in angular acceleration for a square plate.
    Vec2dp<T> const O_b{T{0}, T{0}, T{1}};
    if (P_pivot.x != O_b.x || P_pivot.y != O_b.y) {
        v[2, 2] += m * (P_pivot.x * P_pivot.x + P_pivot.y * P_pivot.y);
    }

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
    T const d = hd::ga::det(I.view());
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

    hd::ga::lu_decomp(A, perm_view);

    // Solve for each column of identity matrix to get inverse
    Inertia2dp<T> I_inv;
    auto A_const =
        std::mdspan<double const, std::dextents<size_t, 2>>{A_data.data(), 3, 3};
    auto perm_const = std::mdspan<int const, std::dextents<size_t, 1>>{perm.data(), 3};

    for (size_t col = 0; col < 3; ++col) {
        std::array<double, 3> e = {0.0, 0.0, 0.0};
        e[col] = 1.0;
        auto e_view = std::mdspan<double, std::dextents<size_t, 1>>{e.data(), 3};
        hd::ga::lu_backsubs(A_const, perm_const, e_view);

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
//   Omega_dot = I_inv[ F - rcmt(Omega, I[Omega]) ]
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

} // namespace hd::ga::pga

/////////////////////////////////////////////////////////////////////////////////////////
// user convenience type for pga2dp physics (based on value_t)
/////////////////////////////////////////////////////////////////////////////////////////
namespace hd::ga {

using inertia2dp = pga::Inertia2dp<value_t>; // 3x3 inertia matrix for 2D rigid body

} // namespace hd::ga
