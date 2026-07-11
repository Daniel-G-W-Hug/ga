#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_solver.hpp" // hd::ga::lu_decomp / lu_backsubs / det
#include "ga_pga3dp_ops.hpp"
#include "ga_usr_utilities.hpp" // hd::ga::rk4_step (shared RK4 integrator)
#include "ga_value_t.hpp"       // for value_t used in convenience type alias

#include <algorithm> // std::find, std::reverse
#include <array>
#include <cmath>      // std::abs
#include <functional> // std::function (time-varying applied wrench)
#include <limits>     // std::numeric_limits
#include <mdspan>
#include <optional> // std::optional (multistep integrator state)
#include <stdexcept>
#include <string>
#include <unordered_map> // std::unordered_map (frame name -> index)
#include <utility>       // std::pair, std::move (assemble_mass_bias return)
#include <vector>


namespace hd::ga::pga {

// Forward declaration of the (optional) closed-loop layer in
// ga_pga3dp_ops_constraints.hpp. dynamic_system3dp grants it friendship so the
// constrained KKT / assembly solver can reuse the tree's private assembly seam
// (assemble_mass_bias, the joint screws, dof_joints/is_ancestor, apply_joint_state)
// WITHOUT widening the open-loop public API. Open-loop users never include the
// constraints header.
class closed_loop_system3dp;

////////////////////////////////////////////////////////////////////////////////
// Inertia3dp: Inertia matrix for 3D projective GA (6x6 matrix)
//
// Used for rigid body dynamics in PGA3DP. The inertia map I[Omega] maps the
// rate of change Omega (a BiVec3dp) to momentum (also a BiVec3dp).
//
// From ga_docu/5_ga_modelling_mechanics.tex eq. 604-611:
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
// BiVec3dp index layout: (vx=e41, vy=e42, vz=e43, mx=e23, my=e31, mz=e12). The weight
// (vx,vy,vz) and bulk (mx,my,mz) play OPPOSITE physical roles in the input twist vs. the
// output momentum:
//
//   input twist Omega:  weight (vx,vy,vz) = ANGULAR velocity omega  (rotation lives in
//   the
//                       weight: exp() uses weight_nrm_sq over (vx,vy,vz) as the rotation
//                       angle), bulk (mx,my,mz) = LINEAR velocity v (a translation has
//                       zero weight, only bulk).
//   output I(Omega):    (vx,vy,vz) = LINEAR momentum p = m v,
//                       (mx,my,mz) = ANGULAR momentum L = J_rot omega.
//
// Block structure of the base matrix:
//
//   Upper-left  [0:3, 0:3] = 0:          no linear-velocity -> linear-momentum coupling
//   Upper-right [0:3, 3:6] = m*Identity: linear velocity (bulk) -> linear momentum
//                                        (vx,vy,vz), i.e. Newton p = m v
//   Lower-left  [3:6, 0:3] = J_rot:      angular velocity (vx,vy,vz) -> angular momentum
//                                        (mx,my,mz), i.e. L = J_rot omega
//   Lower-right [3:6, 3:6] = 0:          no angular-velocity -> angular-momentum (cm)
//   cross
//
// Where J_rot is diagonal with the classical rectangle-rule moments of inertia (acting on
// the angular velocity (vx,vy,vz)):
//
//   I[3,0] = m*(h^2+d^2)/12  (moment about e1-axis, depends on e2 and e3 extents)
//   I[4,1] = m*(w^2+d^2)/12  (moment about e2-axis, depends on e1 and e3 extents)
//   I[5,2] = m*(w^2+h^2)/12  (moment about e3-axis, depends on e1 and e2 extents)
//
// Why the inertia map is "crossed" (weight-input lands in the bulk output and vice
// versa):
//
//   the energy pairing <Omega, I(Omega)> = 2*KE is the regressive product
//   rwdg(Omega, I(Omega)), which pairs the WEIGHT of one bivector with the BULK of the
//   other (see rwdg(BiVec3dp, BiVec3dp)). So the map necessarily sends the angular
//   velocity (weight, vx,vy,vz) into the angular-momentum (bulk, mx,my,mz) slots and the
//   linear velocity (bulk) into the linear-momentum (weight) slots. This is the exact 3D
//   analogue of the 2D case where mass sits off-diagonal.
//
// Optional L_pivot parameter (default = zero bivector = no correction):
// When L_pivot represents a line offset from the body origin, the scalar
// parallel-axis (Steiner) corrections are applied to the lower-left block only:
//
//   I[3,0] += m*(Py²+Pz²)   (moment about x-axis through origin)
//   I[4,1] += m*(Px²+Pz²)   (moment about y-axis through origin)
//   I[5,2] += m*(Px²+Py²)   (moment about z-axis through origin)
//
// where P_foot is the foot of the perpendicular from the body origin O_b to L_pivot:
//
//   n   = (L_pivot.vx, L_pivot.vy, L_pivot.vz)   (line direction, ideal part)
//   mom = (L_pivot.mx, L_pivot.my, L_pivot.mz)   (Plücker moment = p × n)
//   P_foot = (n × mom) / |n|²                    (closest point on L_pivot to O_b)
//
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

    // Upper-right block [0:3, 3:6]: linear velocity (bulk) -> linear momentum (mass m)
    v[0, 3] = m; // mx -> vx: linear vel (e23) -> linear momentum e41, Newton p = m v
    v[1, 4] = m; // my -> vy: linear vel (e31) -> linear momentum e42
    v[2, 5] = m; // mz -> vz: linear vel (e12) -> linear momentum e43

    // Lower-left block [3:6, 0:3]: angular velocity (weight) -> angular momentum
    // (moments)
    v[3, 0] = m * (h * h + d * d) / T{12}; // vx -> mx: I_xx = m*(h^2+d^2)/12, L = J omega
    v[4, 1] = m * (w * w + d * d) / T{12}; // vy -> my: I_yy = m*(w^2+d^2)/12
    v[5, 2] = m * (w * w + h * h) / T{12}; // vz -> mz: I_zz = m*(w^2+h^2)/12

    // Apply scalar parallel-axis (Steiner) corrections if pivot line is offset
    // from the body origin. A line through the origin has zero moment.
    T const v_sq =
        L_pivot.vx * L_pivot.vx + L_pivot.vy * L_pivot.vy + L_pivot.vz * L_pivot.vz;
    if (v_sq != T{0}) {
        // Foot of perpendicular from O_b to L_pivot: P_foot = (n × m_moment) / |n|²
        // where n = (vx,vy,vz) is the line direction
        // and m_moment = (mx,my,mz) is the Plücker moment
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


// Inertia map of a uniform solid DISC / CYLINDER of radius r and thickness / height t,
// total mass m, with the SYMMETRY (spin) axis along e3 -- the disc lies in the e1-e2
// plane, thickness t along e3 -- and the body origin at the centre. Principal moments
// about the centroid:
//
//   axial  I_zz = m r^2 / 2           (about e3, the symmetry axis) transverse
//   I_xx = I_yy = m (r^2/4 + t^2/12)  (about e1 and e2)
//
// Same crossed layout (linear mass in the upper-right block, moments in the lower-left)
// and optional parallel-axis (Steiner) convention as get_cuboid_inertia -- only the
// lower-left diagonal block differs.
template <typename T>
    requires(std::floating_point<T>)
Inertia3dp<T> get_disc_inertia(T m, T r, T t, BiVec3dp<T> const& L_pivot = BiVec3dp<T>{})
{
    Inertia3dp<T> I;
    auto v = I.view();

    // Upper-right block: linear velocity (bulk) -> linear momentum (mass m)
    v[0, 3] = m;
    v[1, 4] = m;
    v[2, 5] = m;

    // Lower-left block: angular velocity (weight) -> angular momentum (moments). Symmetry
    // axis e3: e1, e2 are transverse, e3 is axial.
    T const I_trans = m * (r * r / T{4} + t * t / T{12});
    v[3, 0] = I_trans;          // I_xx (about e1, transverse)
    v[4, 1] = I_trans;          // I_yy (about e2, transverse)
    v[5, 2] = m * r * r / T{2}; // I_zz (about e3, axial)

    // Scalar parallel-axis (Steiner) corrections on the diagonal moments only (identical
    // to get_cuboid_inertia; off-diagonal coupling terms must NOT be added -- see there).
    T const v_sq =
        L_pivot.vx * L_pivot.vx + L_pivot.vy * L_pivot.vy + L_pivot.vz * L_pivot.vz;
    if (v_sq != T{0}) {
        T const inv_v_sq = T{1} / v_sq;
        T const Px = (L_pivot.vy * L_pivot.mz - L_pivot.vz * L_pivot.my) * inv_v_sq;
        T const Py = (L_pivot.vz * L_pivot.mx - L_pivot.vx * L_pivot.mz) * inv_v_sq;
        T const Pz = (L_pivot.vx * L_pivot.my - L_pivot.vy * L_pivot.mx) * inv_v_sq;
        v[3, 0] += m * (Py * Py + Pz * Pz);
        v[4, 1] += m * (Px * Px + Pz * Pz);
        v[5, 2] += m * (Px * Px + Py * Py);
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


/////////////////////////////////////////////////////////////////////////////////////////
// pose3dp: a rigid pose RELATIVE to a parent frame, the 3D twin of pose2dp{origin, phi}.
// In 2D the orientation is a scalar angle phi; in 3D it is a rotation VECTOR rot = axis *
// angle (both origin and rot are vec3dp -- origin a finite point with w = 1, rot a
// direction with w = 0). The pose is encoded into / decoded from the body->parent motor
// by the two converters below (the 3D analog of 2D's exp(vec2dp(0,0,phi)) build + atan2
// decode).
/////////////////////////////////////////////////////////////////////////////////////////

struct pose3dp {
    vec3dp origin{0.0, 0.0, 0.0, 1.0}; // origin of frame in parent coordinates (w = 1)
    vec3dp rot{0.0, 0.0, 0.0, 0.0};    // orientation vs. parent: axis * angle (w = 0)
};


// Build the body->parent motor M = translate(origin) (x) rotate(rot) from a pose. The
// rotation about the parent origin is exp(0.5 * {rot.x,rot.y,rot.z, 0,0,0}) (the weight /
// (vx,vy,vz) slots carry the rotation axis*angle); the translation is get_motor(origin as
// a direction). Inverse of pose3dp_from_motor (below).
inline mvec3dp_e motor_from_pose3dp(pose3dp const& p)
{
    auto const M_rot = exp(0.5 * bivec3dp{p.rot.x, p.rot.y, p.rot.z, 0.0, 0.0, 0.0});
    auto const M_tra = get_motor(vec3dp{p.origin.x, p.origin.y, p.origin.z, 0.0});
    return rgpr(M_tra,
                M_rot); // rotate about origin first, then translate -> body->parent
}

// Decode a body->parent motor M into a pose (origin, rot). The origin is where the body
// origin lands: unitize(move3dp(O, M)). The orientation is recovered by stripping the
// translation (M_rot = T(-origin) (x) M, a pure rotation about the origin) and reading
// the axis*angle off the rotor: the bivector weight (vx,vy,vz) is sin(angle/2)*axis and
// the pseudoscalar is cos(angle/2). A constrained motor log (no general screw log needed
// since the translation is split off via the origin). Inverse of motor_from_pose3dp.
inline pose3dp pose3dp_from_motor(mvec3dp_e const& M)
{
    auto const o = unitize(move3dp(O_3dp, M)); // body origin in parent (w = 1)

    // strip the translation: M_rot = T(-origin) (x) M is a pure rotation about the origin
    auto const T_inv = get_motor(vec3dp{-o.x, -o.y, -o.z, 0.0});
    auto const M_rot = rgpr(T_inv, M);

    auto const w = gr2(M_rot); // rotor bivector part
    value_t const wn = std::sqrt(w.vx * w.vx + w.vy * w.vy + w.vz * w.vz); // sin(angle/2)
    value_t const c = value_t(gr4(M_rot)); // cos(angle/2) (pseudoscalar)

    vec3dp rot{0.0, 0.0, 0.0, 0.0};
    if (wn > eps) {
        value_t const angle = 2.0 * std::atan2(wn, c);
        value_t const s = angle / wn; // angle / sin(angle/2): scales weight -> axis*angle
        rot = vec3dp{s * w.vx, s * w.vy, s * w.vz, 0.0};
    }
    return pose3dp{vec3dp{o.x, o.y, o.z, 1.0}, rot};
}


/////////////////////////////////////////////////////////////////////////////////////////
// static_frame3dp / static_system3dp: a tree of right-handed coordinate frames, each
// posed RELATIVE to its parent. The 3D twin of static_frame2dp/static_system2dp --
// algorithm identical (motor composition + LCA tree walk). The per-frame pose is a
// pose3dp{origin, rot} (the full 2D mirror of pose2dp{origin, phi}); step_pos_trafo
// BUILDS the relative motor from it (just like 2D builds it from origin + phi), and the
// dynamics round-trips through pose3dp via the constrained motor log above.
/////////////////////////////////////////////////////////////////////////////////////////

class static_frame3dp {

    // Static basis frame (a coordinate system without change over time). The parent frame
    // is assumed right-handed (e1,e2,e3) with e4 the homogeneous/origin direction.

    std::string name; // display name (default: 'I' = inertial frame)
    pose3dp pose;     // pose vs. parent coordinate system (ALWAYS relative to parent!)

  public:

    static_frame3dp(std::string frame_name = "I", vec3dp const& origin_in_parent = O_3dp,
                    vec3dp const& rot_vs_parent = vec3dp{0.0, 0.0, 0.0, 0.0}) :
        name(frame_name), pose(origin_in_parent, rot_vs_parent)
    {
        if (pose.origin.w != 1.0) {
            throw std::runtime_error(
                std::string("static_frame3dp: Unitized projective point expected. "
                            "Provided origin.w == ") +
                std::to_string(pose.origin.w));
        }
    }

    std::string const& get_name() const { return name; }
    pose3dp get_pose() const { return pose; }

    // reposition the frame relative to its parent (origin expected unitized, w = 1)
    void set_pose(pose3dp const& new_pose) { pose = new_pose; }
};


class static_system3dp {

    std::vector<static_frame3dp> vfr; // reference frames (index 0 == root)
    std::vector<size_t> parent_of;    // parent frame index per frame (root: itself)
    std::unordered_map<std::string, size_t> name_to_idx; // frame name -> index

  public:

    static_system3dp() = default; // create an empty system

    // sentinel for add_frame: select the default parent (the previously added frame)
    static size_t constexpr prev_frame = std::numeric_limits<size_t>::max();

    // Add a frame to the system. By default its parent is the previously added frame, so
    // a plain sequence of add_frame() calls builds a linear chain rf[0] -> rf[1] -> ...;
    // pass an explicit parent_idx to branch the tree off an earlier frame. A frame's pose
    // is always interpreted RELATIVE to its parent. The first frame added is the root;
    // its parent is itself (a self-loop terminating the upward walk in get_pos_trafo).
    void add_frame(static_frame3dp const& rf, size_t parent_idx = prev_frame)
    {
        size_t const new_idx = vfr.size();
        if (parent_idx == prev_frame) {
            parent_idx =
                (new_idx == 0) ? 0 : new_idx - 1; // default: previous (root: self)
        }
        else if (parent_idx >= new_idx) {
            throw std::runtime_error(
                std::string("static_system3dp: parent_idx must refer to an already added "
                            "frame in [0,") +
                std::to_string(new_idx) + std::string("), but provided parent_idx == ") +
                std::to_string(parent_idx));
        }
        if (name_to_idx.contains(rf.get_name())) {
            throw std::runtime_error(
                std::string("static_system3dp: frame name must be unique, but '") +
                rf.get_name() + std::string("' was already added"));
        }
        name_to_idx.emplace(rf.get_name(), new_idx);
        vfr.push_back(rf);
        parent_of.push_back(parent_idx);
    }

    // Look up a frame index by its name (throws if no such frame exists).
    size_t index_of(std::string const& frame_name) const
    {
        auto const it = name_to_idx.find(frame_name);
        if (it == name_to_idx.end()) {
            throw std::runtime_error(std::string("static_system3dp: no frame named '") +
                                     frame_name + std::string("'"));
        }
        return it->second;
    }

    // Point transformation M for a point p_from (or vector/bivector) in system from_idx
    // to system to_idx, by walking the frame tree: up from `from` to the lowest common
    // ancestor (LCA), then down to `to`. Identical algorithm to static_system2dp (motor
    // composition via rgpr; child->parent step is rrev(step_pos_trafo(child))).
    //
    //   To be used as "p_to = move3dp(p_from, M);".
    mvec3dp_e get_pos_trafo(size_t from_idx, size_t to_idx)
    {
        if (vfr.empty() || from_idx >= vfr.size() || to_idx >= vfr.size()) {
            throw std::runtime_error(
                std::string("static_system3dp: Index range of static_frame3dp must be "
                            "within permissible limits [0,") +
                std::to_string(vfr.size()) + std::string("), but provided from_idx == ") +
                std::to_string(from_idx) + std::string(", to_idx == ") +
                std::to_string(to_idx));
        }

        // identity transformation (M is the pseudoscalar, the neutral element of rgpr())
        if (from_idx == to_idx) return I_3dp_mv_e;

        auto const to_chain = ancestor_chain(to_idx); // [to, parent(to), ..., root]

        // M_up: from -> LCA. Each child -> parent step is rrev(step_pos_trafo(child)); a
        // further-up step is multiplied on the LEFT.
        mvec3dp_e M_up = I_3dp_mv_e;
        size_t node = from_idx;
        while (std::find(to_chain.begin(), to_chain.end(), node) == to_chain.end()) {
            M_up = rgpr(rrev(step_pos_trafo(node)), M_up);
            node = parent_of[node];
        }
        size_t const lca = node;

        // M_down: LCA -> to. Each parent -> child step is step_pos_trafo(child); the
        // deepest child (to) ends up on the far LEFT.
        mvec3dp_e M_down = I_3dp_mv_e;
        for (size_t j = 0; to_chain[j] != lca; ++j) {
            M_down = rgpr(M_down, step_pos_trafo(to_chain[j]));
        }

        return rgpr(M_down, M_up); // apply the up-segment first, then the down-segment
    }

    mvec3dp_e get_pos_trafo(std::string const& from_name, std::string const& to_name)
    {
        return get_pos_trafo(index_of(from_name), index_of(to_name));
    }

    // Single-step parent(child_idx) -> child_idx point transform, BUILT from child_idx's
    // pose3dp (stored relative to its parent). motor_from_pose3dp gives the body->parent
    // motor; its regressive reverse is the parent->child transform. Mirrors 2D's
    // step_pos_trafo (which builds rgpr(M_rot, M_tra) from origin + phi).
    mvec3dp_e step_pos_trafo(size_t child_idx) const
    {
        return rrev(motor_from_pose3dp(vfr[child_idx].get_pose())); // parent -> child
    }

    bool is_linear_chain() const
    {
        for (size_t i = 1; i < parent_of.size(); ++i) {
            if (parent_of[i] != i - 1) return false;
        }
        return true;
    }

    bool is_valid_for_transformations() const { return vfr.size() >= 2; }

    size_t size() const { return vfr.size(); }
    bool empty() const { return vfr.empty(); }
    static_frame3dp const& frame(size_t idx) const { return vfr[idx]; }
    size_t parent(size_t idx) const { return parent_of[idx]; }

    // reposition frame idx relative to its parent (origin expected unitized, w = 1)
    void set_pose(size_t idx, pose3dp const& p) { vfr[idx].set_pose(p); }

  private:

    // chain of ancestors from idx up to the root: [idx, parent(idx), ..., root]
    std::vector<size_t> ancestor_chain(size_t idx) const
    {
        std::vector<size_t> chain{idx};
        while (parent_of[idx] != idx) { // root is its own parent
            idx = parent_of[idx];
            chain.push_back(idx);
        }
        return chain;
    }
};


/////////////////////////////////////////////////////////////////////////////////////////
// kinematic_system3dp: a static_system3dp augmented with momentary RELATIVE velocity and
// acceleration twists per frame. The 3D twin of kinematic_system2dp. The physical input
// is a (linear velocity, angular velocity) pair (kin_state3dp); the dimension-specific
// to_twist() packs it into the abstract se(3) motor generator (a genuine BiVec3dp, grade
// 2
// -- where 2D used a grade-1 vec2dp with a scalar omega). world_VA propagates the twists
// root -> frame by the Newton-Euler adjoint + commutator (bracket) relations.
/////////////////////////////////////////////////////////////////////////////////////////

// A 3dp twist (instantaneous screw), the se(3) motor generator. exp(0.5 * twist) is the
// motor and move3dp(twist, M) is its adjoint. In the BiVec3dp{vx,vy,vz,mx,my,mz} encoding
// the (vx,vy,vz)=e41,e42,e43 (weight) slots carry the ANGULAR velocity (the rotation
// axis*rate, as in exp()'s weight_nrm_sq angle) and (mx,my,mz)=e23,e31,e12 (bulk) slots
// the LINEAR velocity -- decode as omega = (vx,vy,vz), v = (mx,my,mz). The alias
// documents intent at every signature and adds no overloads (it IS bivec3dp); 2D's
// twist2dp was a vec2dp -- see the 2D->3D notes.
using twist3dp = bivec3dp;


// Momentary kinematic state of a frame RELATIVE to its parent (physical inputs). Both
// quantities are direction vectors (w = 0). Carries NO pose -- the pose is held by the
// frame's static_frame3dp. (2D used a scalar omega; in 3D the angular velocity is a
// vector
// -- the natural lift.)
struct kin_state3dp {
    vec3dp vel{0.0, 0.0, 0.0, 0.0}; // linear velocity of the frame origin vs. parent
    vec3dp acc{0.0, 0.0, 0.0, 0.0}; // linear acceleration of the frame origin vs. parent
    vec3dp omega{0.0, 0.0, 0.0, 0.0}; // angular velocity vs. parent (axis * rate)
    vec3dp alpha{0.0, 0.0, 0.0, 0.0}; // angular acceleration vs. parent (axis * rate)
};


class kinematic_system3dp : public static_system3dp {

    // Relative velocity / acceleration twists per frame (the body-frame twists), kept in
    // sync with the base frame list via the add_frame overrides below.
    std::vector<twist3dp> rel_vtwist;
    std::vector<twist3dp> rel_atwist;

  public:

    // pack a (linear velocity, angular velocity) pair into a relative twist -- the only
    // dimension-specific builder (readiness rule #3). Angular velocity -> weight slots
    // (e41,e42,e43); linear velocity -> bulk slots (e23,e31,e12).
    static twist3dp to_twist(vec3dp const& v, vec3dp const& omega)
    {
        return twist3dp{omega.x, omega.y, omega.z, v.x, v.y, v.z};
    }

    // add a frame together with its momentary kinematic state (relative to the parent);
    // the frame's static_frame3dp still provides the initial pose
    void add_frame(static_frame3dp const& f, kin_state3dp const& k,
                   size_t parent_idx = prev_frame)
    {
        static_system3dp::add_frame(f, parent_idx);
        rel_vtwist.push_back(to_twist(k.vel, k.omega));
        rel_atwist.push_back(to_twist(k.acc, k.alpha));
    }

    // add a frame at rest (zero kinematic state) -- hides the base add_frame so
    // rel_vtwist/rel_atwist always stay in sync
    void add_frame(static_frame3dp const& f, size_t parent_idx = prev_frame)
    {
        add_frame(f, kin_state3dp{}, parent_idx);
    }

    // set / replace the momentary kinematic state (relative to the parent) of an
    // already-added frame -- e.g. build a static tree first, then drive a frame
    void set_state(size_t idx, kin_state3dp const& k)
    {
        rel_vtwist[idx] = to_twist(k.vel, k.omega);
        rel_atwist[idx] = to_twist(k.acc, k.alpha);
    }

    void set_state(std::string const& frame_name, kin_state3dp const& k)
    {
        set_state(index_of(frame_name), k);
    }

    // set the relative velocity twist directly (the Lie-algebra-native setter; the
    // physical alternative is set_state(...))
    void set_twist(size_t idx, twist3dp const& B) { rel_vtwist[idx] = B; }

    void set_twist(std::string const& frame_name, twist3dp const& B)
    {
        rel_vtwist[index_of(frame_name)] = B;
    }

    // read-only access to the stored relative twist of a frame (the body-frame twist;
    // decode as omega = (vx,vy,vz), v = (mx,my,mz))
    twist3dp relative_twist(size_t idx) const { return rel_vtwist[idx]; }

    twist3dp relative_twist(std::string const& frame_name) const
    {
        return relative_twist(index_of(frame_name));
    }

    // World-frame velocity twist of frame idx: the sum along the root -> idx path of each
    // relative twist transported into the world frame by that frame's world motor (the
    // adjoint move3dp). Spatial twists add (Lie algebra): V_i = V_parent + Ad(xi_i).
    twist3dp twist_world(size_t idx)
    {
        twist3dp V{};
        for (size_t n = idx; parent(n) != n; n = parent(n)) {
            V = V + move3dp(rel_vtwist[n], get_pos_trafo(n, 0)); // Ad_{M_n}(xi_n)
        }
        return V;
    }

    twist3dp twist_world(std::string const& frame_name)
    {
        return twist_world(index_of(frame_name));
    }

    // set the relative ACCELERATION twist directly (Lie-algebra-native; parallel to
    // set_twist). Used by the dynamics tier to set q-ddot per frame and to zero it for
    // the velocity-product (bias) pass of forward_dynamics.
    void set_accel_twist(size_t idx, twist3dp const& B) { rel_atwist[idx] = B; }

    twist3dp relative_accel_twist(size_t idx) const { return rel_atwist[idx]; }

    twist3dp relative_accel_twist(std::string const& frame_name) const
    {
        return relative_accel_twist(index_of(frame_name));
    }

    // Velocity field of a twist V at point X -- the PGA rate of change of a point:
    //
    //   Xdot = rcmt(V, X)      (3_ga_modelling_motion.tex). In 3D the twist is a bivector
    //
    // (twist3dp) and rcmt(bivec, vec) -> vec; argument ORDER matters:
    // rcmt(V,X)==-rcmt(X,V).
    static vec3dp velocity_field(twist3dp const& V, vec3dp const& X)
    {
        return rcmt(V, X);
    }

    // Acceleration field at point X of a rigid body with velocity twist V and
    // acceleration twist A (5_ga_modelling_mechanics.tex, "Moving coordinate systems"):
    //
    //   a(X) = rcmt(A, X)            [frame/Euler (alpha x r) + origin acceleration]
    //        + rcmt(V, rcmt(V, X))   [centripetal]
    static vec3dp accel_field(twist3dp const& V, twist3dp const& A, vec3dp const& X)
    {
        return rcmt(A, X) + rcmt(V, rcmt(V, X));
    }

    // Velocity of a world-space point X rigidly attached to frame idx.
    vec3dp point_velocity(vec3dp const& X_world, size_t idx)
    {
        return velocity_field(twist_world(idx), X_world);
    }

    vec3dp point_velocity(vec3dp const& X_world, std::string const& frame_name)
    {
        return point_velocity(X_world, index_of(frame_name));
    }

    // World-frame acceleration twist A = dV/dt of frame idx.
    twist3dp accel_twist_world(size_t idx) { return world_VA(idx).A; }

    twist3dp accel_twist_world(std::string const& frame_name)
    {
        return accel_twist_world(index_of(frame_name));
    }

    // Acceleration of a world-space point X rigidly attached to frame idx (the rigid-body
    // acceleration field of its world velocity + acceleration twists).
    vec3dp point_acceleration(vec3dp const& X_world, size_t idx)
    {
        auto const va = world_VA(idx);
        return accel_field(va.V, va.A, X_world);
    }

    vec3dp point_acceleration(vec3dp const& X_world, std::string const& frame_name)
    {
        return point_acceleration(X_world, index_of(frame_name));
    }

  private:

    // World velocity & acceleration twists of frame idx, propagated root -> idx by the
    // recursive Newton-Euler relations (twists transported to world by the adjoint
    // move3dp):
    //
    //   V_i = V_parent + Ad(xi_i)
    //   A_i = A_parent + Ad(xidot_i) + [V_i, Ad(xi_i)]   (Coriolis / centrifugal
    //   coupling)
    // The se(3) twist Lie bracket [.,.] is the regressive commutator
    // rcmt(BiVec3dp,BiVec3dp) (the 3D twin of the se(2) bracket, which is exactly rcmt of
    // the vec2dp twists).
    struct world_va3dp {
        twist3dp V; // velocity twist
        twist3dp A; // acceleration twist
    };

    world_va3dp world_VA(size_t idx)
    {
        std::vector<size_t> path; // root -> idx
        for (size_t n = idx;; n = parent(n)) {
            path.push_back(n);
            if (parent(n) == n) break; // reached the root (self-parent)
        }
        std::reverse(path.begin(), path.end());

        twist3dp V{};
        twist3dp A{};
        for (size_t n : path) {
            if (parent(n) == n) continue; // root contributes nothing
            auto const M = get_pos_trafo(n, 0);
            auto const zeta = move3dp(rel_vtwist[n], M);    // Ad(xi_n)     world rel. vel
            auto const zetadot = move3dp(rel_atwist[n], M); // Ad(xidot_n)  world rel. acc
            V = V + zeta;
            A = A + zetadot +
                rcmt(V, zeta); // rcmt = the se(3) twist Lie bracket [V, zeta]
        }
        return {V, A};
    }
};


/////////////////////////////////////////////////////////////////////////////////////////
// dynamic_system3dp: the forces/inertia tier on top of kinematic_system3dp. The 3D twin
// of dynamic_system2dp. A free body's relative acceleration twist is COMPUTED from the
// applied wrench via the se(3) Euler equation
//
//   Omega_dot = I^-1[ W - rcmt(Omega, I(Omega)) ]            (= compute_omega_dot)
//
// -- LITERALLY the same line as 2D (vec) with BiVec3dp instead. The pose is evolved on
// the motor manifold M(t) = M0 (x) exp(1/2 B) with RK4 on the Lie-algebra pair (B,
// Omega). Bodies may also be connected by 1-DOF JOINTS into articulated chains (coupled
// forward dynamics) with force elements (spring/damper, applied wrench, grounded spring,
// and subclass-provided elements via extra_wrenches()) feeding the generalised force.
/////////////////////////////////////////////////////////////////////////////////////////

// Rigid-body inertial properties of a frame (body frame, about the body origin = cm).
struct body3dp {
    Inertia3dp<value_t> I;     // inertia map: body twist -> body momentum
    Inertia3dp<value_t> I_inv; // its inverse (cached)
    value_t mass{0.0};         // total mass (gravity + energy)
};

// Build a body3dp for a uniform cuboid (extents w,h,d along e1,e2,e3) of total mass m,
// with the body origin at the centre of mass.
inline body3dp make_cuboid_body(value_t m, value_t w, value_t h, value_t d)
{
    auto const I =
        get_cuboid_inertia(m, w, h, d); // about cm (default pivot = body origin)
    return body3dp{I, get_inertia_inverse(I), m};
}

// Build a body3dp for a uniform disc/cylinder of radius r and thickness/height t
// (symmetry/spin axis along e3) of total mass m, body origin at the centre.
inline body3dp make_disc_body(value_t m, value_t r, value_t t)
{
    auto const I = get_disc_inertia(m, r, t); // about cm (default pivot = body origin)
    return body3dp{I, get_inertia_inverse(I), m};
}

// Joint type connecting a body to its parent (the reduced-coordinate degrees of freedom).
// Both 1-DOF kinds run through the SAME code -- only the screw generator differs (the PGA
// unification of rotation and translation): M(q) = rest (x) exp(1/2 q * screw).
//
//   free      : unconstrained 6-DOF rigid body; state lives in the base layer.
//
//   revolute  : 1-DOF hinge; screw = a LINE (the rotation axis); q rotates about it.
//
//   prismatic : 1-DOF slider; screw = an IDEAL line (translation dir); q translates it.
enum class joint3dp { free, revolute, prismatic };

// Per-frame joint state (parallel to the body[] list). Meaningful for 1-DOF joints; the
// free body ignores it and uses the base relative twist + pose directly.
struct joint_state3dp {
    joint3dp type{joint3dp::free};
    // joint screw generator (a twist3dp = BiVec3dp), in the body frame: a unit
    // rotation-axis LINE for a revolute hinge, an ideal line (translation generator) for
    // a prismatic slider. The same field drives both via the same exponential -- see
    // joint3dp above.
    twist3dp screw_b{};
    mvec3dp_e rest;     // body->parent motor at q = 0 (the reference pose)
    value_t phi{0.0};   // generalised coordinate q (revolute angle / prismatic distance)
    value_t omega{0.0}; // generalised rate q-dot

    // Optional linear force elements acting on the generalised coordinate q, folded into
    // the joint-space generalised force tau:
    //
    //     tau += -stiffness * (q - q_rest)  -  damping * q-dot
    //
    // A torsional spring/damper for a revolute joint, a linear spring/damper for a
    // prismatic slider (the PGA unification carries through). All default to zero, so the
    // gravity/bias-only path is byte-unchanged. The spring also contributes its potential
    // 1/2 k (q - q_rest)^2 to potential_energy(); the damper is dissipative.
    value_t stiffness{0.0}; // generalised spring constant k
    value_t damping{0.0};   // generalised damping constant c
    value_t q_rest{0.0};    // spring rest coordinate q0
};

// A GROUNDED spatial spring + damper: a linear spring/damper connecting a body-fixed
// attachment point P (anchor_b, in the body frame) to a fixed anchor P0 (p0_world) in the
// inertial frame -- an aerostatic bearing / bushing. Unlike joint_state3dp's
// spring-on-a-coordinate, this is a CONFIGURATION-dependent force element: at every
// sub-step it reads the live world position/velocity of P and produces the restoring
// force (world axes, anisotropic stiffness)
//
//     F = -[k.x dx, k.y dy, k.z dz] - c v_P ,   d = P_world - P0 ,   v_P = point velocity
//
// applied as the force line (wrench) wdg(P_world, F). Because the wrench acts at the
// physical attachment point, BOTH the translational stiffness AND the tilt stiffness
// (lever arm of P about the body origin) emerge geometrically from one element -- e.g.
// two radial springs at axial stations +-l give a tilt stiffness ~ k l^2 without
// any separate torsional spring. The spring also contributes its potential
// 1/2 (k.x dx^2 + k.y dy^2 + k.z dz^2) to potential_energy(); the damper is dissipative.
struct grounded_spring3dp {
    vec3dp anchor_b{O_3dp};       // body-fixed attachment point (w = 1, body frame)
    vec3dp p0_world{O_3dp};       // inertial anchor / rest position (w = 1, world frame)
    vec3dp k{0.0, 0.0, 0.0, 0.0}; // anisotropic stiffness along world e1/e2/e3
    value_t c{0.0};               // linear (isotropic) damping on the point velocity
};

// Time-integration scheme selectable on dynamic_system3dp (see set_integrator):
//   rk4  -- canonical 4th-order Runge-Kutta (default), 4 rhs evals/step
//   abm2 -- Adams-Bashforth-Moulton 2nd-order predictor-corrector, 2 rhs evals/step
// (both from ga/ga_usr_utilities.hpp). RK4 is the robust default; ABM2 trades order for
// fewer evaluations. Both are EXPLICIT -- neither is an implicit/stiff solver, so on a
// genuinely stiff system the step is bounded by stability, not accuracy.
enum class integrator_kind { rk4, abm2 };

class dynamic_system3dp : public kinematic_system3dp {

  public:

    // Virtual so an application subclass (e.g. a contact/penalty force model)
    // can be owned and deleted polymorphically.
    virtual ~dynamic_system3dp() = default;

  protected:

    // Extension point for a subclass: inject configuration-dependent world
    // wrenches into the force assembly, evaluated at each RK4 sub-step. Each
    // returned (frame_idx, wrench) pair is folded onto that frame's supporting
    // joints exactly like an applied wrench (same spatial_dot(S_j, W) pairing).
    // The generic base contributes none; a subclass overrides this to add its
    // own force elements without the base knowing anything about them.
    virtual std::vector<std::pair<size_t, bivec3dp>> extra_wrenches() { return {}; }

  private:


    // The optional closed-loop layer (ga_pga3dp_ops_constraints.hpp) composes a
    // dynamic_system3dp as its spanning tree and reuses this class's private assembly
    // seam (assemble_mass_bias, the joint screws, dof_joints/is_ancestor,
    // apply_joint_state) to build the loop-closure constraint Jacobian and the
    // constrained dynamics on top -- without those internals becoming public API.
    friend class closed_loop_system3dp;

    std::vector<body3dp> body;         // per-frame inertial properties (index = frame)
    std::vector<joint_state3dp> joint; // per-frame joint state (index = frame)
    // uniform gravity field [world frame, a direction with w = 0]
    vec3dp grav{0.0, -9.81, 0.0, 0.0};

    // Optional time-varying applied wrench per frame (world frame), folded into tau as an
    // external generalised force. Sparse map (frame index -> wrench function of time); a
    // frame without an entry has no applied wrench. The simulation clock `time_` is the
    // argument passed at each RK4 sub-step (set transiently by coupled_step). Only the
    // jointed (assemble_mass_bias) path consumes it; free bodies see gravity only.
    using wrench_fn = std::function<bivec3dp(value_t)>;
    std::unordered_map<size_t, wrench_fn> wrench_;
    value_t time_{0.0}; // simulation clock [s], advanced by step()

    // Optional KINEMATICALLY DRIVEN joints: a 1-DOF joint whose coordinate is PRESCRIBED
    // q(t) = q0 + rate*t (constant rate) rather than integrated from the dynamics -- a
    // motor-driven spin or a steady feed. A joint is "driven" iff it appears in this map
    // (it keeps its revolute/prismatic screw machinery, but is excluded from dof_joints
    // and instead re-evaluated from q(t) at each sub-step). The driven motion is a moving
    // base for the dynamic sub-chain below it: its velocity (omega*screw) propagates into
    // the descendants' Newton-Euler bias, so centrifugal/Coriolis/gyroscopic effects
    // emerge. (Constant rate => zero joint acceleration; a general prescribed q(t) with
    // q-ddot != 0 is a future extension.)
    struct driven_spec {
        value_t rate{0.0}; // prescribed q-dot
        value_t q0{0.0};   // q at t = 0
    };
    std::unordered_map<size_t, driven_spec> driven_;

    // Optional GROUNDED spatial springs/dampers per frame (a frame may carry several --
    // e.g. two radial springs + an axial spring). Configuration-dependent: their
    // restoring wrench is recomputed from the live pose/velocity at each RK4 sub-step, so
    // it is folded into assemble_mass_bias (the q-ddot-dependent path) like an applied
    // wrench, NOT prescribed as a function of time. See grounded_spring3dp.
    std::unordered_map<size_t, std::vector<grounded_spring3dp>> springs_;

    // Selectable time integrator for the coupled joint chain (coupled_step). RK4
    // (default, the canonical mdspan rk4_step) or ABM2 (Adams-Bashforth-Moulton 2nd-order
    // multistep). ABM2 is a MULTISTEP method, so its history (the previous derivative)
    // persists across step() calls in `abm_`, lazily sized to the dof count and reset on
    // a switch / dof change. Both share one rhs evaluation (forward dynamics), so RK4
    // stays byte-identical to the hand-rolled loop and ABM2 is a drop-in alternative. See
    // set_integrator.
    integrator_kind integ_{integrator_kind::rk4};
    std::optional<abm2_integrator> abm_;

  public:

    dynamic_system3dp() = default;

    // add a frame WITHOUT inertia (e.g. the inertial root); keeps body[]/joint[] in sync
    // with the base frame list. Mirrors the two base add_frame overloads.
    void add_frame(static_frame3dp const& f, kin_state3dp const& k,
                   size_t parent_idx = prev_frame)
    {
        kinematic_system3dp::add_frame(f, k, parent_idx);
        body.push_back(body3dp{});
        joint.push_back(joint_state3dp{});
    }

    void add_frame(static_frame3dp const& f, size_t parent_idx = prev_frame)
    {
        kinematic_system3dp::add_frame(f, parent_idx);
        body.push_back(body3dp{});
        joint.push_back(joint_state3dp{});
    }

    // add a free (6-DOF) dynamic rigid body: rest pose + inertial properties + initial
    // kinematic state (linear & angular velocity)
    void add_body(static_frame3dp const& f, body3dp const& b,
                  kin_state3dp const& k = kin_state3dp{}, size_t parent_idx = prev_frame)
    {
        add_frame(f, k, parent_idx);
        body.back() = b;
    }

    // add a revolute-jointed body: a 1-DOF hinge about the body-fixed axis line through
    // `pivot_b` (a finite point, w = 1) along the unit direction `axis_b` (w = 0). The
    // joint screw is that line, screw = wdg(pivot_b, axis_b); exp(1/2 q screw) is a
    // rotation about it. The effective hinge inertia (Steiner) emerges from the spatial
    // Jacobian in the forward dynamics; nothing extra to precompute.
    void add_revolute_body(static_frame3dp const& f, body3dp const& b,
                           vec3dp const& pivot_b, vec3dp const& axis_b,
                           value_t phi0 = 0.0, value_t omega0 = 0.0,
                           size_t parent_idx = prev_frame)
    {
        add_screw_joint(f, b, joint3dp::revolute, wdg(pivot_b, axis_b), phi0, omega0,
                        parent_idx);
    }

    // add a prismatic-jointed body: a 1-DOF slider along the body-fixed unit direction
    // `dir`. The joint screw is the translation generator (an ideal line, zero weight):
    // exp(1/2 s screw) is a pure translation by s along dir. The machinery is IDENTICAL
    // to the revolute joint -- only the generator differs (the PGA unification).
    void add_prismatic_body(static_frame3dp const& f, body3dp const& b, vec3dp const& dir,
                            value_t s0 = 0.0, value_t v0 = 0.0,
                            size_t parent_idx = prev_frame)
    {
        // translation generator along dir: the ideal line (bulk = linear) {0,0,0, dir}
        add_screw_joint(f, b, joint3dp::prismatic,
                        twist3dp{0.0, 0.0, 0.0, dir.x, dir.y, dir.z}, s0, v0, parent_idx);
    }

    void set_gravity(vec3dp const& g) { grav = g; } // g is a direction (w = 0)
    vec3dp gravity() const { return grav; }

    // Attach a time-varying applied wrench (world frame) to frame `idx`, folded into tau
    // as an external generalised force projected onto each supporting joint screw. The
    // wrench is a bivec3dp (force/torque line); for a pure force F through a point P use
    // wdg(P, F). Evaluated at each RK4 sub-step time. Pass an empty function to clear.
    void set_applied_wrench(size_t idx, wrench_fn fn)
    {
        if (fn) wrench_[idx] = std::move(fn);
        else wrench_.erase(idx);
    }

    value_t time() const { return time_; }  // simulation clock [s]
    void set_time(value_t t) { time_ = t; } // reset / seed the clock

    // Make a 1-DOF joint KINEMATICALLY DRIVEN at constant rate: q(t) = q0 + rate*t. The
    // joint must already exist (added via add_revolute_body / add_prismatic_body); it is
    // then excluded from the dynamic DOFs and prescribed instead, acting as a moving base
    // for the dynamic sub-chain below it. Use for a motor-driven spin or a steady feed.
    void set_driven_rate(size_t idx, value_t rate, value_t q0 = 0.0)
    {
        driven_[idx] = driven_spec{rate, q0};
        apply_driven_joints(); // seed the joint state at the current clock
    }

    void clear_driven_joint(size_t idx) { driven_.erase(idx); }
    bool is_driven(size_t idx) const { return driven_.count(idx) != 0; }

    // read-only access to a body's inertial properties
    body3dp const& body_props(size_t idx) const { return body[idx]; }

    value_t joint_phi(size_t idx) const { return joint[idx].phi; }     // joint coordinate
    value_t joint_omega(size_t idx) const { return joint[idx].omega; } // joint rate

    // Attach a linear spring + damper to a 1-DOF joint's generalised coordinate q:
    // tau += -k*(q - q0) - c*q-dot. Additive to the gravity/bias generalised forces; the
    // coupled forward dynamics and step() pick it up automatically. Only revolute /
    // prismatic joints consume it (a free/rigid frame ignores these fields).
    void set_joint_spring_damper(size_t idx, value_t k, value_t c, value_t q0 = 0.0)
    {
        joint[idx].stiffness = k;
        joint[idx].damping = c;
        joint[idx].q_rest = q0;
    }

    // Attach a GROUNDED spatial spring + damper to frame `idx`: a body-fixed point
    // `anchor_b` (body frame, w = 1) tied to the inertial anchor `p0_world` by
    // anisotropic stiffness `k` (world axes) and isotropic damping `c`. If `p0_world` is
    // omitted the current world position of the attachment point (the present
    // configuration) is taken as the rest position. A frame may carry several (returns
    // nothing; call repeatedly). See grounded_spring3dp for the force law and why tilt
    // stiffness emerges.
    void add_grounded_spring(size_t idx, vec3dp const& anchor_b, vec3dp const& k,
                             value_t c = 0.0)
    {
        vec3dp const p0 = unitize(move3dp(anchor_b, get_pos_trafo(idx, 0)));
        springs_[idx].push_back(grounded_spring3dp{anchor_b, p0, k, c});
    }

    void add_grounded_spring(size_t idx, vec3dp const& anchor_b, vec3dp const& p0_world,
                             vec3dp const& k, value_t c)
    {
        springs_[idx].push_back(grounded_spring3dp{anchor_b, p0_world, k, c});
    }

    // Select the time integrator for the coupled joint chain. RK4 (default) or ABM2
    // (Adams-Bashforth-Moulton 2nd-order). Switching resets the ABM2 multistep history so
    // the next step() self-starts cleanly. RK4 is byte-identical to the previous
    // hand-rolled loop; ABM2 is a drop-in alternative sharing the same forward-dynamics
    // rhs.
    void set_integrator(integrator_kind k)
    {
        integ_ = k;
        abm_.reset(); // restart the multistep history on any switch
    }
    integrator_kind get_integrator() const { return integ_; }

    void clear_grounded_springs(size_t idx) { springs_.erase(idx); }

    // Current acceleration of joint `idx`, from the COUPLED joint-space forward dynamics
    // at the present state (no integration).
    value_t joint_accel(size_t idx)
    {
        auto const rj = dof_joints();
        auto const qdd = forward_dynamics(rj);
        for (size_t k = 0; k < rj.size(); ++k)
            if (rj[k] == idx) return qdd[k];
        return 0.0;
    }

    // Recompute the joint accelerations from the applied forces (forward dynamics) and
    // write them into the per-frame relative acceleration twists, so accel_twist_world /
    // point_acceleration then return the ACTUAL dynamic accelerations (not just the
    // velocity-product bias left after step()).
    void sync_accelerations()
    {
        auto const rj = dof_joints();
        auto const qdd = forward_dynamics(rj); // zeroes rel_atwist internally
        for (size_t k = 0; k < rj.size(); ++k)
            set_accel_twist(rj[k], qdd[k] * joint[rj[k]].screw_b);
    }

    // Advance the system by dt. The 1-DOF joints form a COUPLED chain integrated together
    // in their reduced (joint) coordinates via the joint-space forward dynamics; free
    // bodies are integrated independently. RK4 throughout.
    void step(value_t dt)
    {
        auto const rj = dof_joints();
        if (!rj.empty())
            coupled_step(rj, dt); // uses time_ for sub-step wrench/drive eval
        for (size_t i = 1; i < size(); ++i)
            if (joint[i].type == joint3dp::free && body[i].mass > 0.0)
                step_free_body(i, dt);
        time_ += dt;           // advance the clock (coupled_step restores it to t0)
        apply_driven_joints(); // prescribe the driven joints at t + dt (final state)
    }

    // --- energy / momentum diagnostics (inertial / world frame) ----------------------

    // total kinetic energy: sum over bodies of 1/2 <V_body, I(V_body)>, the spatial
    // inertia-map quadratic form (= 1/2 m|v_cm|^2 + 1/2 omega.I.omega). V_body is the
    // world twist pulled into the body frame where the inertia map lives.
    // Dimension-agnostic form
    // -- LITERALLY the 2D code with the BiVec3dp rwdg overload.
    value_t kinetic_energy()
    {
        value_t ke = 0.0;
        for (size_t i = 1; i < size(); ++i) {
            auto const Minv = rrev(get_pos_trafo(i, 0));       // world -> body i
            twist3dp const Vb = move3dp(twist_world(i), Minv); // body twist of body i
            ke += 0.5 * spatial_dot(Vb, body[i].I(Vb));
        }
        return ke;
    }

    // total gravitational potential energy: sum over bodies of -m (g . r_cm)
    value_t potential_energy()
    {
        value_t pe = 0.0;
        for (size_t i = 1; i < size(); ++i) {
            vec3dp const cm_w = move3dp(O_3dp, get_pos_trafo(i, 0));
            pe += -body[i].mass * (grav.x * cm_w.x + grav.y * cm_w.y + grav.z * cm_w.z);
            // joint-spring potential 1/2 k (q - q0)^2 (zero unless a spring is attached)
            value_t const dq = joint[i].phi - joint[i].q_rest;
            pe += 0.5 * joint[i].stiffness * dq * dq;
        }
        // grounded-spring potential 1/2 (k.x dx^2 + k.y dy^2 + k.z dz^2)
        for (auto const& [fi, sps] : springs_) {
            auto const M = get_pos_trafo(fi, 0);
            for (auto const& sp : sps) {
                vec3dp const P = unitize(move3dp(sp.anchor_b, M));
                vec3dp const d = P - sp.p0_world;
                pe +=
                    0.5 * (sp.k.x * d.x * d.x + sp.k.y * d.y * d.y + sp.k.z * d.z * d.z);
            }
        }
        return pe;
    }

    value_t total_energy() { return kinetic_energy() + potential_energy(); }

    // World-frame momentum bivector of body idx: the body momentum I(V_body) transported
    // to the world frame by the adjoint. For a torque-free free body this is CONSERVED --
    // the classic Poinsot / Euler check that complements energy conservation.
    bivec3dp momentum_world(size_t idx)
    {
        auto const M = get_pos_trafo(idx, 0);                   // body -> world
        twist3dp const Vb = move3dp(twist_world(idx), rrev(M)); // body twist
        return move3dp(body[idx].I(Vb), M);                     // momentum back to world
    }

    // Joint-space mass matrix M(q) (n*n row-major, n = number of 1-DOF joints) at the
    // current configuration. A diagnostic / showcase quantity: the reduced inertia of the
    // articulated system, with the identity  1/2 * qdot^T M(q) qdot == kinetic_energy().
    std::vector<value_t> mass_matrix()
    {
        auto const rj = dof_joints();
        size_t const n = rj.size();
        std::vector<twist3dp> S(n); // world joint screws (unit rate)
        std::vector<mvec3dp_e> Minv(
            n); // world -> body i motor (pulls a twist into body i)
        for (size_t i = 0; i < n; ++i) {
            auto const M = get_pos_trafo(rj[i], 0); // body i -> world
            S[i] = move3dp(joint[rj[i]].screw_b, M);
            Minv[i] = rrev(M);
        }
        std::vector<value_t> Mmat(n * n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            auto const& I = body[rj[i]].I; // inertia map about body i's cm (body frame)
            for (size_t j = 0; j < n; ++j) {
                if (!is_ancestor(rj[j], rj[i])) continue;
                twist3dp const xj = move3dp(S[j], Minv[i]); // joint-j screw in body i
                for (size_t k = 0; k < n; ++k) {
                    if (!is_ancestor(rj[k], rj[i])) continue;
                    twist3dp const xk = move3dp(S[k], Minv[i]);
                    Mmat[j * n + k] += spatial_dot(xj, I(xk));
                }
            }
        }
        return Mmat;
    }

  private:

    // Spatial (reciprocal / Klein) pairing of a velocity twist with a momentum / wrench
    // bivector -> scalar: the dimension-agnostic inertia-map quadratic form, with
    // <xi, I(xi)> == 2x the kinetic energy of a body moving with body twist xi. In PGA3DP
    // the pairing is the regressive product -rwdg(xi, mom) -- IDENTICAL in form to the 2D
    // dynamic_system2dp::spatial_dot (only the rwdg overload differs: BiVec3dp here,
    // Vec2dp/BiVec2dp there). This is the form that makes the mass matrix lift unchanged.
    static value_t spatial_dot(twist3dp const& xi, bivec3dp const& mom)
    {
        return -value_t(rwdg(xi, mom));
    }

    // RK4-integrate one free rigid body (frame idx) over dt under gravity. The
    // integration state is the Lie-algebra pair (B, Omega): B is the relative generator
    // accumulated from the current relative pose M0 (so M(t) = M0 (x) exp(1/2 B)), Omega
    // the body twist. dB/dt = Omega; dOmega/dt = I^-1[ W_body - rcmt(Omega, I(Omega)) ].
    // For a torque-free body (grav = 0) this reduces to the pure se(3) Euler equation
    // (Poinsot / Dzhanibekov).
    void step_free_body(size_t idx, value_t dt)
    {
        auto const M0 = rrev(step_pos_trafo(idx)); // current body -> parent motor
        auto const& I = body[idx].I;
        auto const& I_inv = body[idx].I_inv;
        value_t const m = body[idx].mass;

        // body-frame twist rate from the gravity wrench acting at the cm (= body origin)
        auto omega_dot = [&](twist3dp const& B, twist3dp const& Om) -> twist3dp {
            auto const M = rgpr(M0, exp(0.5 * B));             // pose at this stage
            auto const W_w = wdg(move3dp(O_3dp, M), m * grav); // gravity wrench (world)
            auto const W_b = move3dp(W_w, rrev(M)); // pulled into the body frame
            return compute_omega_dot(I_inv, W_b, Om, I);
        };

        // RK4 (shared rk4_step) on the Lie-algebra pair u = (B, Omega): dB/dt = Omega,
        // dOmega/dt = omega_dot(B, Omega). B starts at 0 (M(t) = M0 (x) exp(1/2 B)).
        std::array<twist3dp, 2> u_mem{twist3dp{}, relative_twist(idx)};
        std::array<twist3dp, 4> uh_mem{};
        std::array<twist3dp, 2> rhs_mem{};
        auto u = std::mdspan<twist3dp, std::dextents<size_t, 1>>(u_mem.data(), 2);
        auto uh = std::mdspan<twist3dp, std::dextents<size_t, 2>>(uh_mem.data(), 2, 2);
        auto const rhs =
            std::mdspan<twist3dp const, std::dextents<size_t, 1>>(rhs_mem.data(), 2);
        for (size_t s = 1; s <= 4; ++s) {
            rhs_mem[0] = u[1];                  // dB/dt = Omega
            rhs_mem[1] = omega_dot(u[0], u[1]); // dOmega/dt
            rk4_step(u, uh, rhs, dt, s);
        }

        // decode the evolved pose (motor -> pose3dp via the constrained motor log) and
        // write pose + twist back into the base layers
        set_pose(idx, pose3dp_from_motor(rgpr(M0, exp(0.5 * u[0]))));
        set_twist(idx, u[1]);
    }

    // shared constructor for the 1-DOF screw joints (revolute / prismatic): the ONLY
    // difference between them is the screw generator passed in.
    void add_screw_joint(static_frame3dp const& f, body3dp const& b, joint3dp type,
                         twist3dp const& screw_b, value_t q0, value_t qdot0,
                         size_t parent_idx)
    {
        add_frame(f, parent_idx);
        size_t const idx = size() - 1;
        body[idx] = b;
        joint_state3dp js;
        js.type = type;
        js.screw_b = screw_b;
        js.rest = rrev(step_pos_trafo(idx)); // body->parent motor at the rest pose
        js.phi = q0;
        js.omega = qdot0;
        joint[idx] = js;
        apply_joint_state(idx); // write q0/qdot0 into the base pose + twist
    }

    // body->parent motor at generalised coordinate q: M(q) = rest (x) exp(1/2 q * screw).
    // The SAME exponential builds a rotation (revolute, screw = a line) or a translation
    // (prismatic, screw = an ideal line); the operators do not change.
    mvec3dp_e joint_motor(size_t idx, value_t q) const
    {
        return rgpr(joint[idx].rest, exp(0.5 * q * joint[idx].screw_b));
    }

    // Write the joint state (q, q-dot) into the base pose + relative twist, so all
    // kinematic queries (get_pos_trafo, twist_world, point_velocity, energy) stay valid.
    void apply_joint_state(size_t idx)
    {
        auto const M = joint_motor(idx, joint[idx].phi); // body->parent motor at q
        set_pose(idx, pose3dp_from_motor(M));
        set_twist(idx,
                  joint[idx].omega * joint[idx].screw_b); // rel twist = q-dot * screw
    }

    // Re-evaluate every kinematically driven joint at the current clock time_: write its
    // prescribed q(t) = q0 + rate*t and q-dot = rate into the joint state (pose +
    // relative twist). The relative acceleration stays zero (constant rate), so the
    // driven joint's velocity propagates into the descendants' bias while adding no
    // prescribed q-ddot.
    void apply_driven_joints()
    {
        for (auto const& [idx, d] : driven_) {
            joint[idx].phi = d.q0 + d.rate * time_;
            joint[idx].omega = d.rate;
            apply_joint_state(idx);
        }
    }

    // Frame indices of the DYNAMIC 1-DOF joints (revolute or prismatic) -- the
    // generalised coords integrated by the forward dynamics. Kinematically driven joints
    // (in driven_) are excluded: they are prescribed, not solved for.
    std::vector<size_t> dof_joints() const
    {
        std::vector<size_t> rj;
        for (size_t i = 1; i < size(); ++i)
            if ((joint[i].type == joint3dp::revolute ||
                 joint[i].type == joint3dp::prismatic) &&
                driven_.count(i) == 0)
                rj.push_back(i);
        return rj;
    }

    // Is joint frame `jf` on the path from body `bf` up to the root (inclusive)?
    bool is_ancestor(size_t jf, size_t bf) const
    {
        for (size_t n = bf;; n = parent(n)) {
            if (n == jf) return true;
            if (parent(n) == n) return false; // reached the root
        }
    }

    // Assemble the joint-space mass matrix M(q) and the generalised-force RHS for the
    // 1-DOF joint chain `rj`, by virtual work over the bodies in the dimension-agnostic
    // SPATIAL (screw) form:
    //
    //   M[j][k] = sum_i  spatial_dot( S_j^body_i , I_i( S_k^body_i ) )      (mass matrix)
    //   RHS[j]  = sum_i [ m_i vcm_i(S_j).g  -  spatial_dot( S_j^body_i, F_bias_i ) ]
    //
    // over bodies i with joint j (and k) as ancestor. S_j = move3dp(screw_j,
    // M_{j->world}) is the world joint screw; S_j^body_i = move3dp(S_j, rrev(M_i))
    // transports it into body i's frame (where its inertia map I_i lives). The gravity
    // term is a pure cm force (m_i vcm_j.g, no torque about the cm). The body-frame BIAS
    // WRENCH at q-ddot = 0 is the full spatial Newton-Euler bias:
    //
    //   F_bias_i = I_i(A_bias_body_i) + rcmt(V_body_i, I_i(V_body_i))
    //
    // where V_body_i / A_bias_body_i are the world velocity / bias-acceleration twists
    // pulled into body i, and the second term is the se(3) gyroscopic (velocity-product)
    // wrench. This is the GENERAL form: for a single joint both the angular bias and the
    // gyroscopic projection onto its own screw vanish (ad-invariance: <X,[X,Y]> == 0), so
    // it reduces to the 2D translational result; for coupled NON-PARALLEL axes the
    // angular bias + gyroscopic terms are non-zero and required. Pre: the joint state
    // (phi, omega) is already applied.
    //
    // SIDE EFFECT: this runs the bias pass, zeroing the chain's relative accel twists
    // (rel_atwist) so the world accel queries return only the q-ddot-independent part.
    //
    // Split out of forward_dynamics() as the reuse seam for the closed-loop layer
    // (ga_pga3dp_ops_constraints.hpp): the same spatial-Jacobian columns (velocity_field)
    // that build M and RHS also build the loop-closure constraint Jacobian, so the
    // constrained KKT solver assembles on top of this without duplicating the inertia-map
    // assembly. Private (the open-loop public surface is unchanged); the closed-loop
    // layer reaches it through friendship (see the forward declaration above). Returns {
    // Mmat (n*n, row-major), RHS (n) }.
    std::pair<std::vector<value_t>, std::vector<value_t>>
    assemble_mass_bias(std::vector<size_t> const& rj)
    {
        size_t const n = rj.size();

        // velocity-product (bias) pass: zero the relative accel twists so the world accel
        // queries return only the q-ddot-independent (Coriolis/gyroscopic) bias part.
        for (size_t k = 0; k < n; ++k)
            set_accel_twist(rj[k], twist3dp{});

        // World joint screws of the dof coordinates (the n unknowns).
        std::vector<twist3dp> S(n);
        for (size_t j = 0; j < n; ++j)
            S[j] = move3dp(joint[rj[j]].screw_b, get_pos_trafo(rj[j], 0));

        // Inertia-bearing bodies of the articulated chain: the dof joints AND the
        // kinematically DRIVEN joints. A driven joint is a MOVING BASE -- its body
        // inertia still loads its ancestor dof joints (mass matrix) and its prescribed
        // velocity feeds their Newton-Euler bias, so a spinning driven rotor produces the
        // gyroscopic / clamped dynamics (Omega-ddot = 0). Omitting driven bodies here
        // would silently drop their inertia (and, if a driven joint carries the only
        // inertia below a dof joint, make the mass matrix singular). Its relative
        // acceleration is zero (constant rate, apply_driven_joints), so the bias is pure
        // velocity-product.
        std::vector<size_t> bl(rj);
        for (auto const& [idx, d] : driven_)
            bl.push_back(idx);
        size_t const nb = bl.size();

        std::vector<vec3dp> cm(nb);      // body cm in world
        std::vector<mvec3dp_e> Minv(nb); // world -> body i motor
        std::vector<twist3dp> Fbias(nb); // body-frame spatial bias wrench
        std::vector<value_t> mass(nb);
        for (size_t i = 0; i < nb; ++i) {
            size_t const fb = bl[i];
            auto const M = get_pos_trafo(fb, 0);
            cm[i] = move3dp(O_3dp, M);
            Minv[i] = rrev(M);
            mass[i] = body[fb].mass;
            // body-frame velocity / bias-acceleration twists + the spatial bias wrench
            twist3dp const Vb = move3dp(twist_world(fb), Minv[i]);
            twist3dp const Ab = move3dp(accel_twist_world(fb), Minv[i]); // rel_atwist = 0
            auto const& I = body[fb].I;
            Fbias[i] = I(Ab) + rcmt(Vb, I(Vb)); // I*A_bias + gyroscopic V x* (I V)
        }

        std::vector<value_t> Mmat(n * n, 0.0), RHS(n, 0.0);
        for (size_t j = 0; j < n; ++j) {
            for (size_t i = 0; i < nb; ++i) { // contribution of body i to coordinate j
                if (!is_ancestor(rj[j], bl[i])) continue;
                vec3dp const vj = velocity_field(S[j], cm[i]); // cm velocity, unit rate j
                twist3dp const xj = move3dp(S[j], Minv[i]);    // joint-j screw in body i
                RHS[j] += mass[i] * (vj.x * grav.x + vj.y * grav.y + vj.z * grav.z) -
                          spatial_dot(xj, Fbias[i]);
                auto const& I = body[bl[i]].I; // inertia about body i's cm (body frame)
                for (size_t k = 0; k < n; ++k) {
                    if (!is_ancestor(rj[k], bl[i])) continue;
                    twist3dp const xk = move3dp(S[k], Minv[i]);
                    Mmat[j * n + k] += spatial_dot(xj, I(xk));
                }
            }
        }

        // linear spring/damper generalised forces on each joint coordinate (additive,
        // diagonal in joint space): tau_j += -k_j (q_j - q0_j) - c_j q-dot_j. Zero unless
        // a spring/damper was attached via set_joint_spring_damper.
        for (size_t j = 0; j < n; ++j) {
            auto const& js = joint[rj[j]];
            RHS[j] += -js.stiffness * (js.phi - js.q_rest) - js.damping * js.omega;
        }

        // applied external wrenches (world frame, evaluated at the current clock time_):
        // an applied wrench W on frame fi contributes the generalised force
        // spatial_dot(S_j, W) to every joint j that supports fi (j ancestor of fi). The
        // reciprocal pairing is the rate of work of W under unit joint rate -- the same
        // pairing that yields the gravity term. Zero unless a wrench was attached.
        for (auto const& [fi, fn] : wrench_) {
            if (!fn) continue;
            bivec3dp const W = fn(time_);
            for (size_t j = 0; j < n; ++j)
                if (is_ancestor(rj[j], fi)) RHS[j] += spatial_dot(S[j], W);
        }

        // grounded spatial springs/dampers: for each spring on frame fi, the live world
        // attachment point P and its velocity v_P give the restoring force F (world-axis
        // anisotropic stiffness + isotropic damping); the force line wdg(P, F) is the
        // wrench, projected onto every supporting joint screw. Recomputed from state here
        // (not a function of time) -- this is the configuration-dependent path.
        for (auto const& [fi, sps] : springs_) {
            auto const M = get_pos_trafo(fi, 0);
            twist3dp const Vw = twist_world(fi); // world velocity twist of frame fi
            for (auto const& sp : sps) {
                vec3dp const P = unitize(move3dp(sp.anchor_b, M)); // world point (w = 1)
                vec3dp const v = velocity_field(Vw, P);            // world point velocity
                vec3dp const F{-sp.k.x * (P.x - sp.p0_world.x) - sp.c * v.x,
                               -sp.k.y * (P.y - sp.p0_world.y) - sp.c * v.y,
                               -sp.k.z * (P.z - sp.p0_world.z) - sp.c * v.z, 0.0};
                bivec3dp const W = wdg(P, F);
                for (size_t j = 0; j < n; ++j)
                    if (is_ancestor(rj[j], fi)) RHS[j] += spatial_dot(S[j], W);
            }
        }

        // application-specific external wrenches contributed by a subclass (e.g. a
        // contact/penalty force model): each (frame, world-wrench) pair is projected onto
        // its frame's supporting joints exactly like the grounded-spring path above.
        // Empty in the generic base (extra_wrenches() returns nothing).
        for (auto const& [fi, W] : extra_wrenches()) {
            for (size_t j = 0; j < n; ++j)
                if (is_ancestor(rj[j], fi)) RHS[j] += spatial_dot(S[j], W);
        }
        return {std::move(Mmat), std::move(RHS)};
    }

    // Joint-space forward dynamics for the chain `rj`: returns the joint accelerations
    // q-ddot solving  M(q) q-ddot = RHS(q, q-dot). Thin wrapper over assemble_mass_bias()
    // (see there for the assembly and its bias-pass side effect) plus the shared LU
    // solve.
    std::vector<value_t> forward_dynamics(std::vector<size_t> const& rj)
    {
        auto const [Mmat, RHS] = assemble_mass_bias(rj);
        return hd::ga::lu_solve(Mmat, RHS, rj.size()); // shared LU (detail/ga_solver.hpp)
    }

    // RK4-integrate the coupled 1-DOF joint chain `rj` over dt in its joint coordinates
    // via the shared rk4_step. The state is u = [phi_0..phi_{n-1}, omega_0..omega_{n-1}];
    // the derivative (dphi, domega) = (omega, q-ddot) is recomputed each sub-step by
    // writing u into the joint state, refreshing the kinematic poses + twists, and
    // solving the coupled forward dynamics.
    void coupled_step(std::vector<size_t> const& rj, value_t dt)
    {
        size_t const n = rj.size();
        size_t const dim = 2 * n; // state = [phi_0..phi_{n-1}, omega_0..omega_{n-1}]
        std::vector<value_t> u_mem(dim);
        for (size_t k = 0; k < n; ++k) {
            u_mem[k] = joint[rj[k]].phi;
            u_mem[n + k] = joint[rj[k]].omega;
        }

        // write u into the joint state (+ refresh the kinematic poses/twists)
        auto write_u = [&](std::vector<value_t> const& u) {
            for (size_t k = 0; k < n; ++k) {
                joint[rj[k]].phi = u[k];
                joint[rj[k]].omega = u[n + k];
                apply_joint_state(rj[k]);
            }
        };

        // shared derivative f(t, u) -> du for both integrators: prescribe the moving base
        // and write u into the joints at the (sub-)step time, then solve the coupled
        // forward dynamics. du = [omega, q-ddot]. Setting time_ here threads the sub-step
        // time into any time-varying applied wrench / driven joint; coupled_step restores
        // time_ to t0 on exit, and step() advances it by dt.
        auto f = [&](value_t t, std::vector<value_t> const& u, std::vector<value_t>& du) {
            time_ = t;
            apply_driven_joints();
            write_u(u);
            auto const qdd = forward_dynamics(rj);
            for (size_t k = 0; k < n; ++k) {
                du[k] = u[n + k];   // dphi/dt = omega
                du[n + k] = qdd[k]; // domega/dt = q-ddot
            }
        };

        value_t const t0 = time_;
        if (integ_ == integrator_kind::abm2) {
            if (!abm_ || abm_->dim() != dim) abm_.emplace(dim); // (re)size the history
            abm_->step(f, u_mem, t0, dt);
        }
        else {
            rk4_integrator(dim).step(f, u_mem, t0, dt); // wraps the canonical rk4_step
        }

        time_ = t0; // restore the clock (step() advances it by dt)
        write_u(
            u_mem); // write the integrated state back into the joints + kinematic layers
    }
};

} // namespace hd::ga::pga

// value_t convenience aliases (inertia3dp, ...) live in ga_usr_types_mechanics.hpp, next
// to the other mechanics convenience types.
