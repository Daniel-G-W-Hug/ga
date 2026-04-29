#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_solver.hpp" // hd::ga::lu_decomp / lu_backsubs / det
#include "ga_pga3dp_ops.hpp"
#include "ga_value_t.hpp" // for value_t used in convenience type alias

#include <array>
#include <mdspan>
#include <stdexcept>


namespace hd::ga::pga {

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

template <typename T>
    requires(numeric_type<T>)
struct Inertia3dp {
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

    // Subtract inertia contribution
    constexpr Inertia3dp& operator-=(Inertia3dp const& other)
    {
        for (size_t i = 0; i < 36; ++i) {
            data[i] -= other.data[i];
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

    // Apply inertia map: I[Omega]
    // (map rate of change bivector Omega to momentum bivector in 3D)
    // or
    // Apply inverse inertia map: I_inv[arg]
    // (map momentum bivector arg to bivector Omega in 3D)
    //
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


// Create inertia matrix for a uniform rectangular cuboid centered at the origin.
//
// The cuboid has width w along e1, height h along e2, depth d along e3.
// It is assumed to be centered at the origin in the body frame (Xw = 1).
//
// This is the exact continuous limit of accumulating point masses over the volume.
// Derivation: integrate get_point_inertia over volume element dV = dx dy dz,
// using mean(x^2) = w^2/12, mean(y^2) = h^2/12, mean(z^2) = d^2/12,
// and mean(x) = mean(y) = mean(z) = 0, mean(x*y) = mean(x*z) = mean(y*z) = 0
// (centered, independent variables).
//
// Base result about cm / body origin O_b = (0,0,0,1):
//
//   I_cm = m * [  0              0              0             1    0    0  ]
//              [  0              0              0             0    1    0  ]
//              [  0              0              0             0    0    1  ]
//              [ (h^2+d^2)/12    0              0             0    0    0  ]
//              [  0             (w^2+d^2)/12    0             0    0    0  ]
//              [  0              0             (w^2+h^2)/12   0    0    0  ]
//
// BiVec3dp index layout: (vx=e41, vy=e42, vz=e43, mx=e23, my=e31, mz=e12)
//   Indices 0-2 (vx,vy,vz): translational velocity / linear momentum  (ideal lines, e4*)
//   Indices 3-5 (mx,my,mz): angular velocity / angular momentum        (real lines,
//   e23/e31/e12)
//
// Block structure of the base matrix:
//   Upper-left  [0:3, 0:3] = 0:          no translational-to-translational coupling
//   Upper-right [0:3, 3:6] = m*Identity: angular velocity -> linear momentum (mass,
//   Newton p=mv) Lower-left  [3:6, 0:3] = m*J_rot:   translational vel -> angular
//   momentum (classical moments) Lower-right [3:6, 3:6] = 0:          no
//   rotational-to-rotational coupling
//
// Where J_rot is diagonal with the classical rectangle-rule moments of inertia:
//   I[3,0] = m*(h^2+d^2)/12  (rotation about e1-axis, depends on e2 and e3 extents)
//   I[4,1] = m*(w^2+d^2)/12  (rotation about e2-axis, depends on e1 and e3 extents)
//   I[5,2] = m*(w^2+h^2)/12  (rotation about e3-axis, depends on e1 and e2 extents)
//
// Why the coupling is "crossed" (upper-right / lower-left rather than diagonal):
//   In PGA3DP the Hodge complement swaps ideal lines (e41,e42,e43) with real lines
//   (e23,e31,e12). Angular velocity (e23/e31/e12, indices 3-5) couples to linear momentum
//   (e41/e42/e43, indices 0-2) via mass m, landing in the upper-right block.
//   Translational velocity (e41/e42/e43, indices 0-2) couples to angular momentum
//   (e23/e31/e12, indices 3-5) via moment of inertia, landing in the lower-left block.
//   This is the exact 3D analogue of the 2D case where mass appears off-diagonal.
//
// Optional L_pivot parameter (default = zero bivector = no correction):
// When L_pivot represents a line offset from the body origin, the scalar
// parallel-axis (Steiner) corrections are applied to the lower-left block only:
//   I[3,0] += m*(Py²+Pz²)   (moment about x-axis through origin)
//   I[4,1] += m*(Px²+Pz²)   (moment about y-axis through origin)
//   I[5,2] += m*(Px²+Py²)   (moment about z-axis through origin)
// where P_foot is the foot of the perpendicular from the body origin O_b to L_pivot:
//   n   = (L_pivot.vx, L_pivot.vy, L_pivot.vz)   (line direction, ideal part)
//   mom = (L_pivot.mx, L_pivot.my, L_pivot.mz)   (Plücker moment = p × n)
//   P_foot = (n × mom) / |n|²                    (closest point on L_pivot to O_b)
// Off-diagonal coupling terms (products of inertia, cross-blocks) must NOT be
// added here. Adding them via get_point_inertia(m, P_foot) would cancel the
// Steiner correction in I_inv (same mechanism as in the 2D case), causing
// alpha = τ/J_cm instead of τ/I_pivot (wrong by factor ~4 for a square plate).
// When L_pivot passes through the body origin (ideal part = 0), no correction.
// Pre: if ideal part != 0 then |n|^2 > 0 (L_pivot must be non-degenerate).
template <typename T>
    requires(std::floating_point<T>)
Inertia3dp<T> get_cuboid_inertia(T m, T w, T h, T d,
                                 BiVec3dp<T> const& L_pivot = BiVec3dp<T>{})
{
    Inertia3dp<T> I;
    auto v = I.view();

    // Upper-right block [0:3, 3:6]: angular velocity -> linear momentum (mass m)
    v[0, 3] = m; // mx -> vx: angular vel about e1 -> linear momentum e41
    v[1, 4] = m; // my -> vy: angular vel about e2 -> linear momentum e42
    v[2, 5] = m; // mz -> vz: angular vel about e3 -> linear momentum e43

    // Lower-left block [3:6, 0:3]: translational velocity -> angular momentum (moments)
    v[3, 0] = m * (h * h + d * d) / T{12}; // vx -> mx: I_xx = m*(h^2+d^2)/12
    v[4, 1] = m * (w * w + d * d) / T{12}; // vy -> my: I_yy = m*(w^2+d^2)/12
    v[5, 2] = m * (w * w + h * h) / T{12}; // vz -> mz: I_zz = m*(w^2+h^2)/12

    // Apply scalar parallel-axis (Steiner) corrections if pivot line is offset
    // from the body origin. A line through the origin has zero ideal part.
    T const v_sq =
        L_pivot.vx * L_pivot.vx + L_pivot.vy * L_pivot.vy + L_pivot.vz * L_pivot.vz;
    if (v_sq != T{0}) {
        // Foot of perpendicular from O_b to L_pivot: P_foot = (n × m_moment) / |n|²
        // where n = (vx,vy,vz) is the line direction (ideal part)
        // and m_moment = (mx,my,mz) is the Plücker moment (real part)
        T const inv_v_sq = T{1} / v_sq;
        T const Px = (L_pivot.vy * L_pivot.mz - L_pivot.vz * L_pivot.my) * inv_v_sq;
        T const Py = (L_pivot.vz * L_pivot.mx - L_pivot.vx * L_pivot.mz) * inv_v_sq;
        T const Pz = (L_pivot.vx * L_pivot.my - L_pivot.vy * L_pivot.mx) * inv_v_sq;
        // Update only the diagonal lower-left block (scalar moments of inertia).
        // Do NOT add off-diagonal coupling terms — they would cancel the Steiner
        // correction in I_inv, giving alpha = τ/J_cm instead of τ/I_pivot.
        v[3, 0] += m * (Py * Py + Pz * Pz); // I_xx += m*(Py²+Pz²)
        v[4, 1] += m * (Px * Px + Pz * Pz); // I_yy += m*(Px²+Pz²)
        v[5, 2] += m * (Px * Px + Py * Py); // I_zz += m*(Px²+Py²) = m*r²
    }

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
    T const d = hd::ga::det(I.view());
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

    hd::ga::lu_decomp(A, perm_view);

    // Solve for each column of identity matrix to get inverse
    Inertia3dp<T> I_inv;
    auto A_const =
        std::mdspan<double const, std::dextents<size_t, 2>>{A_data.data(), 6, 6};
    auto perm_const = std::mdspan<int const, std::dextents<size_t, 1>>{perm.data(), 6};

    for (size_t col = 0; col < 6; ++col) {
        std::array<double, 6> e = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        e[col] = 1.0;
        auto e_view = std::mdspan<double, std::dextents<size_t, 1>>{e.data(), 6};
        hd::ga::lu_backsubs(A_const, perm_const, e_view);

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
//   Omega_dot = I_inv[ F - rcmt(Omega, I[Omega]) ]
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

} // namespace hd::ga::pga

/////////////////////////////////////////////////////////////////////////////////////////
// user convenience type for pga3dp physics (based on value_t)
/////////////////////////////////////////////////////////////////////////////////////////
namespace hd::ga {

using inertia3dp = pga::Inertia3dp<value_t>; // 6x6 inertia matrix for 3D rigid body

} // namespace hd::ga
