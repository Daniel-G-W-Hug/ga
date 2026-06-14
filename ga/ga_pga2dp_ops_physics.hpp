#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_solver.hpp" // hd::ga::lu_decomp / lu_backsubs / det
#include "ga_pga2dp_ops.hpp"
#include "ga_usr_utilities.hpp" // hd::ga::rk4_step (shared RK4 integrator)
#include "ga_value_t.hpp"       // for value_t used in convenience type alias

#include <algorithm> // std::find, std::min, std::max
#include <array>
#include <cmath>      // std::abs
#include <functional> // std::function (time-varying applied wrench)
#include <limits>     // std::numeric_limits
#include <mdspan>
#include <stdexcept>
#include <string>
#include <unordered_map> // std::unordered_map (frame name -> index)
#include <utility>       // std::pair, std::move (assemble_mass_bias return)
#include <vector>


namespace hd::ga::pga {

// Forward declaration of the (optional) closed-loop layer in
// ga_pga2dp_ops_constraints.hpp. dynamic_system2dp grants it friendship so the
// constrained KKT / assembly solver can reuse the tree's private assembly seam
// (assemble_mass_bias, the joint screws, dof_joints/is_ancestor, apply_joint_state)
// WITHOUT widening the open-loop public API. Open-loop users never include the
// constraints header.
class closed_loop_system2dp;

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
//
//   I_cm = m * [  0    1         0        ]
//              [ -1    0         0        ]
//              [  0    0    (w^2+h^2)/12  ]
//
// I[0,1] =  m : total mass (Newton p = mv, e2 -> e23)
// I[1,0] = -m : total mass (Newton p = mv, e1 -> e31, opposite sign)
// I[2,2] = m*(w^2+h^2)/12 : moment of inertia about centroid (classical rectangle
// formula)
//
// Optional P_pivot parameter (default = body origin = cm):
// When P_pivot != O_b, the scalar parallel-axis (Steiner) correction is applied:
//
//   I_pivot[2,2] = J_cm + m*(Px² + Py²)
//
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


/////////////////////////////////////////////////////////////////////////////////////////
// static_frame2dp / static_system2dp: a tree of right-handed coordinate frames, each
// posed RELATIVE to its parent. Provides rigid coordinate transformations between frames.
/////////////////////////////////////////////////////////////////////////////////////////

struct pose2dp {
    vec2dp origin; // origin of frame in parent coordinates
    value_t phi;   // orientation of frame vs. parent coordinates, [phi]: rad
};


// Build the body->parent motor M = translate(origin) (x) rotate(phi) from a pose. The
// rotation about the parent origin is exp(0.5 * vec2dp(0,0,phi)); the translation by
// origin is exp(0.5 * vec2dp(-origin.y, origin.x, 0)). Parallel to the 3D
// motor_from_pose3dp (the 2D rotation generator is a scalar phi where 3D uses an
// axis*angle bivector). Inverse of pose2dp_from_motor. (static_system2dp::step_pos_trafo
// returns its regressive reverse, the parent->child transform.)
inline mvec2dp_u motor_from_pose2dp(pose2dp const& p)
{
    auto const M_rot = exp(0.5 * vec2dp(0.0, 0.0, p.phi));
    auto const M_tra = exp(0.5 * vec2dp(-p.origin.y, p.origin.x, 0.0));
    return rgpr(M_tra, M_rot); // rotate about origin, then translate -> body->parent
}

// Decode a body->parent motor M into a pose (origin, phi). The origin is where the body
// origin lands, unitized: move2dp(O, M); the orientation is the angle of the rotated e1
// axis (a vector, so translation-invariant): atan2 of move2dp(e1, M). The 2D analog of
// the 3D pose3dp_from_motor -- here the rotation read-off is a single atan2 (no rotor
// axis-angle extraction needed). Inverse of motor_from_pose2dp.
inline pose2dp pose2dp_from_motor(mvec2dp_u const& M)
{
    auto const o = move2dp(O_2dp, M);                  // body origin in parent coords
    auto const e1 = move2dp(vec2dp{1.0, 0.0, 0.0}, M); // rotated e1 direction
    return pose2dp{vec2dp{o.x / o.z, o.y / o.z, 1.0}, std::atan2(e1.y, e1.x)};
}


class static_frame2dp {

    // Static basis frame, a coordinate system (cs) w/o change over time.
    // Cs parent frame is assumed to be right-handed with e1 in x-direction, e2 rotated
    // by 90° anti-clockwise vs. e1 defining positive phi direction in direction of e12

    std::string name; // display name (default: 'I' = inertial frame)

    // pose vs. parent coordinate system (remember: ALLWAYS relative to parent system!)
    pose2dp pose;

  public:

    static_frame2dp(std::string frame_name = "I", vec2dp const& origin_in_parent = O_2dp,
                    value_t angle_vs_parent = value_t(0.0)) :
        name(frame_name), pose(origin_in_parent, angle_vs_parent)
    {
        if (pose.origin.z != 1.0) {
            throw std::runtime_error(
                std::string("static_frame2dp: Unitized projective point expected. "
                            "Provided origin.z == ") +
                std::to_string(pose.origin.z));
        }
    }

    std::string const& get_name() const { return name; }
    pose2dp get_pose() const { return pose; }

    // reposition the frame relative to its parent (origin is expected unitized, z = 1)
    void set_pose(vec2dp const& new_origin, value_t new_phi)
    {
        pose.origin = new_origin;
        pose.phi = new_phi;
    }
};


class static_system2dp {

    std::vector<static_frame2dp> vfr; // reference frames (index 0 == root)
    std::vector<size_t> parent_of;    // parent frame index per frame (root: itself)
    std::unordered_map<std::string, size_t> name_to_idx; // frame name -> index

  public:

    static_system2dp() = default; // create an empty system

    // sentinel for add_frame: select the default parent (the previously added frame)
    static size_t constexpr prev_frame = std::numeric_limits<size_t>::max();

    // Add a frame to the system. By default its parent is the previously added frame, so
    // a plain sequence of add_frame() calls builds a linear dependency chain
    //
    //   rf[0] (root) -> rf[1] -> rf[2] -> ...
    //
    // Pass an explicit parent_idx to branch the tree off an earlier frame instead. This
    // allows for example a robots with two arms both linked to a body frame that can
    // rotate vs. the world system. The coordinate transform between arbitrary systems in
    // the tree will work automatically based on the parent->child relationships).
    //
    // A frame's pose is always interpreted RELATIVE to its parent (see static_frame2dp).
    // The first frame added is the root; its parent is itself (a self-loop terminating
    // the upward walk in get_pos_trafo).
    void add_frame(static_frame2dp const& rf, size_t parent_idx = prev_frame)
    {
        size_t const new_idx = vfr.size();
        if (parent_idx == prev_frame) {
            parent_idx =
                (new_idx == 0) ? 0 : new_idx - 1; // default: previous (root: self)
        }
        else if (parent_idx >= new_idx) {
            throw std::runtime_error(
                std::string("static_system2dp: parent_idx must refer to an already added "
                            "frame in [0,") +
                std::to_string(new_idx) + std::string("), but provided parent_idx == ") +
                std::to_string(parent_idx));
        }
        if (name_to_idx.contains(rf.get_name())) {
            throw std::runtime_error(
                std::string("static_system2dp: frame name must be unique, but '") +
                rf.get_name() + std::string("' was already added"));
        }
        name_to_idx.emplace(rf.get_name(), new_idx);
        vfr.push_back(rf);
        parent_of.push_back(parent_idx);
    }

    // Look up a frame index by its name (throws if no such frame exists). Lets callers
    // address frames by name instead of raw index, e.g. get_pos_trafo("H0", "P").
    size_t index_of(std::string const& frame_name) const
    {
        auto const it = name_to_idx.find(frame_name);
        if (it == name_to_idx.end()) {
            throw std::runtime_error(std::string("static_system2dp: no frame named '") +
                                     frame_name + std::string("'"));
        }
        return it->second;
    }

    // Return a point transformation M for a point p_from (or vector) in system from_idx
    // to a point p_to (or vector) in system to_idx, by walking the frame tree: up from
    // `from` to the lowest common ancestor (LCA), then down to `to`.
    //
    //   Forward transformation:   from_idx is an ancestor of to_idx   (down only)
    //   Identity transformation:  from_idx == to_idx
    //   Backward transformation:  to_idx is an ancestor of from_idx   (up only)
    //   General (branching) case: from/to in different subtrees    (up to LCA, then down)
    //
    // A linear chain is the degenerate tree where each parent is the predecessor: then
    // the LCA is the shallower frame, one segment is empty, and this reduces to the
    // consecutive-index composition.
    //
    // To be used as "p_to = move2dp(p_from, M);" with M as returned from get_pos_trafo()
    //
    mvec2dp_u get_pos_trafo(size_t from_idx, size_t to_idx)
    {
        // check permissible range
        // (at least one static_frame2dp required for identity, otherwise at least two)
        if (vfr.empty() || from_idx >= vfr.size() || to_idx >= vfr.size()) {
            throw std::runtime_error(
                std::string("static_system2dp: Index range of static_frame2dp must be "
                            "within permissible limits [0,") +
                std::to_string(vfr.size()) + std::string("), but provided from_idx == ") +
                std::to_string(from_idx) + std::string(", to_idx == ") +
                std::to_string(to_idx));
        }

        // identity transformation (M is the pseudoscalar, the neutral element of rgpr())
        if (from_idx == to_idx) return I_2dp_mv_u; // unit pseudoscalar = identity trafo

        // ancestor chain of `to` (incl. to and the root); used to locate the LCA
        auto const to_chain = ancestor_chain(to_idx); // [to, parent(to), ..., root]

        // M_up: from -> LCA. Each child -> parent step is rrev(step_pos_trafo(child)).
        // A step further up is applied AFTER the accumulated motor -> multiplied on the
        // LEFT (regressive sandwich: move2dp(move2dp(p, A), B) == move2dp(p, rgpr(B,
        // A))).
        mvec2dp_u M_up = I_2dp_mv_u;
        size_t node = from_idx;
        while (std::find(to_chain.begin(), to_chain.end(), node) == to_chain.end()) {
            M_up = rgpr(rrev(step_pos_trafo(node)), M_up);
            node = parent_of[node];
        }
        size_t const lca = node;

        // M_down: LCA -> to. Each parent -> child step is step_pos_trafo(child). Walk the
        // prefix of to_chain above the LCA; the deepest child (to) is applied last and so
        // ends up on the far LEFT (rgpr is associative, so left-folding preserves order).
        mvec2dp_u M_down = I_2dp_mv_u;
        for (size_t j = 0; to_chain[j] != lca; ++j) {
            M_down = rgpr(M_down, step_pos_trafo(to_chain[j]));
        }

        // apply the up-segment first, then the down-segment
        return rgpr(M_down, M_up);
    }

    // Name-addressed convenience overload of get_pos_trafo (resolves names to indices).
    mvec2dp_u get_pos_trafo(std::string const& from_name, std::string const& to_name)
    {
        return get_pos_trafo(index_of(from_name), index_of(to_name));
    }

    // Single-step rigid coordinate change from parent(child_idx) to frame child_idx,
    // built from child_idx's pose -- stored RELATIVE to its parent (see static_frame2dp:
    // "pose vs. parent coordinate system"). The parent->child transform is the regressive
    // reverse of the body->parent motor built by motor_from_pose2dp.
    //
    // Translation and rotation are COMPOSED, not summed as generators: they do not
    // commute (Baker-Campbell-Hausdorff), so the single-exp form would be a screw about a
    // wrong centre. The composition reduces correctly to a pure translation (phi == 0) or
    // pure rotation (origin == parent origin).
    mvec2dp_u step_pos_trafo(size_t child_idx) const
    {
        return rrev(motor_from_pose2dp(vfr[child_idx].get_pose())); // parent -> child
    }

    // True if every frame's parent is its predecessor, i.e. the system is a plain linear
    // chain rather than a branching tree (root counts as linear).
    bool is_linear_chain() const
    {
        for (size_t i = 1; i < parent_of.size(); ++i) {
            if (parent_of[i] != i - 1) return false;
        }
        return true;
    }

    // at least two reference frames are required for transformations between them
    bool is_valid_for_transformations() const { return vfr.size() >= 2; }

    // read-only access (used by the fmt formatter and callers)
    size_t size() const { return vfr.size(); }
    bool empty() const { return vfr.empty(); }
    static_frame2dp const& frame(size_t idx) const { return vfr[idx]; }
    size_t parent(size_t idx) const { return parent_of[idx]; }

    // reposition frame idx relative to its parent (origin expected unitized, z = 1)
    void set_pose(size_t idx, vec2dp const& origin, value_t phi)
    {
        vfr[idx].set_pose(origin, phi);
    }

    // reposition frame idx from a pose2dp (parallel to the 3D set_pose(idx, pose3dp))
    void set_pose(size_t idx, pose2dp const& p) { vfr[idx].set_pose(p.origin, p.phi); }

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
// kinematic_system2dp: a static_system2dp augmented with a momentary kinematic state per
// frame (velocity & acceleration, linear & angular), all RELATIVE to the parent. Pure
// kinematics -- no forces (that is the future dynamic_system2dp). The initial pose of
// each frame comes from its static_frame2dp (inherited); the kinematic state is layered
// on top.
/////////////////////////////////////////////////////////////////////////////////////////

// A 2dp twist (instantaneous screw: angular + linear velocity, or their accelerations).
// It is stored in a vec2dp, but with a SPECIAL interpretation -- it is NOT an ordinary
// point/vector. In the right-handed plane omega > 0 is counter-clockwise in a
// right-handed system, i.e. in the positive direction of the bivector e12.
//
// The components encode the pga2dp motor generator twist2dp(.x, .y, .z) = (-v_y, v_x,
// omega) so that exp(0.5 * twist) is the motor and move2dp(twist, M) is its adjoint.
// Decode as omega = .z and v = (.y, -.x). The alias documents this intent at every
// signature and adds no overloads (it IS vec2dp). In pga3dp the corresponding twist is a
// genuine BiVec3dp (grade 2) -- see the 2D->3D notes.
using twist2dp = vec2dp;


// Momentary kinematic state of a frame RELATIVE to its parent (physical inputs). Linear
// quantities live in the parent frame (z = 0 directions); angular ones are scalars.
// Carries NO pose -- the pose is held by the frame's static_frame2dp.
struct kin_state2dp {
    vec2dp vel{0.0, 0.0, 0.0}; // linear velocity of the frame origin vs. parent [1/s]
    vec2dp acc{0.0, 0.0,
               0.0};    // linear acceleration of the frame origin vs. parent [1/s^2]
    value_t omega{0.0}; // angular velocity vs. parent [rad/s]
    value_t alpha{0.0}; // angular acceleration vs. parent [rad/s^2]
};


class kinematic_system2dp : public static_system2dp {

    // Relative velocity / acceleration twists per frame (see twist2dp for the encoding).
    // One entry per frame, kept in sync with the base frame list via the add_frame
    // overrides below.
    std::vector<twist2dp> rel_vtwist; // from (linear velocity, angular velocity)
    std::vector<twist2dp> rel_atwist; // from (linear acceleration, angular acceleration)

  public:

    // pack a (linear velocity, angular velocity) pair into a relative twist (see
    // twist2dp)
    static twist2dp to_twist(vec2dp const& v, value_t omega)
    {
        return twist2dp(-v.y, v.x,
                        omega); // 2D builder (the only dimension-specific piece)
    }

    // add a frame together with its momentary kinematic state (relative to the parent);
    // the frame's static_frame2dp still provides the initial pose
    void add_frame(static_frame2dp const& f, kin_state2dp const& k,
                   size_t parent_idx = prev_frame)
    {
        static_system2dp::add_frame(f, parent_idx);
        rel_vtwist.push_back(to_twist(k.vel, k.omega));
        rel_atwist.push_back(to_twist(k.acc, k.alpha));
    }

    // add a frame at rest (zero kinematic state) -- hides the base add_frame so the
    // parallel rel_vtwist vector always stays in sync
    void add_frame(static_frame2dp const& f, size_t parent_idx = prev_frame)
    {
        add_frame(f, kin_state2dp{}, parent_idx);
    }

    // set / replace the momentary kinematic state (relative to the parent) of an
    // already-added frame -- e.g. build a static tree first, then drive a joint
    void set_state(size_t idx, kin_state2dp const& k)
    {
        rel_vtwist[idx] = to_twist(k.vel, k.omega);
        rel_atwist[idx] = to_twist(k.acc, k.alpha);
    }

    void set_state(std::string const& frame_name, kin_state2dp const& k)
    {
        set_state(index_of(frame_name), k);
    }

    // set the relative velocity twist directly (the Lie-algebra-native setter; see
    // twist2dp). The physical alternative is set_state(...). Leaves the relative
    // acceleration unchanged.
    void set_twist(size_t idx, twist2dp const& B) { rel_vtwist[idx] = B; }

    void set_twist(std::string const& frame_name, twist2dp const& B)
    {
        rel_vtwist[index_of(frame_name)] = B;
    }

    // set the relative ACCELERATION twist directly (Lie-algebra-native; parallel to
    // set_twist). Used by the dynamics tier to set q-ddot per frame for the Newton-Euler
    // acceleration pass, and to zero it for the velocity-product (bias) pass.
    void set_accel_twist(size_t idx, twist2dp const& B) { rel_atwist[idx] = B; }

    // Advance the system by one time step dt (in-place evolution). Prescribed-motion
    // kinematics only -- no forces; the root frame (index 0) is left unchanged.
    //
    // Each non-root frame's relative pose is evolved on the motor manifold by its
    // body-frame relative twist:
    //
    //     P_new = P (x) exp(0.5 * B_rel * dt)        [P = body->parent motor]
    //
    // which is EXACT for a constant twist (N steps of dt == one step of N*dt). Its
    // relative velocity twist is then ramped by the relative acceleration (explicit
    // Euler):
    //
    //     B_rel += B_accel * dt
    void step(value_t dt)
    {
        for (size_t i = 1; i < size(); ++i) { // frame 0 is the root (left unchanged)
            // evolve the body->parent relative motor by the relative twist over dt
            auto const P = rrev(step_pos_trafo(i)); // T(origin) (x) R(phi)
            auto const P_new = rgpr(P, exp(0.5 * rel_vtwist[i] * dt));
            // decode the new relative pose (origin, phi) and store it
            set_pose(i, pose2dp_from_motor(P_new));
            // ramp the relative velocity by the (constant) relative acceleration
            rel_vtwist[i] = rel_vtwist[i] + rel_atwist[i] * dt;
        }
    }

    // World-frame velocity twist of frame idx: the sum along the root -> idx path of each
    // relative twist transported into the world frame by that frame's world motor (the
    // adjoint move2dp). Spatial twists add (Lie algebra), so V_i = V_parent + Ad(xi_i).
    twist2dp twist_world(size_t idx)
    {
        twist2dp V{0.0, 0.0, 0.0};
        for (size_t n = idx; parent(n) != n; n = parent(n)) {
            V = V + move2dp(rel_vtwist[n], get_pos_trafo(n, 0)); // Ad_{M_n}(xi_n)
        }
        return V;
    }

    twist2dp twist_world(std::string const& frame_name)
    {
        return twist_world(index_of(frame_name));
    }

    // Velocity field of a twist V at point X -- the PGA rate of change of a point:
    //
    //   Xdot = rcmt(V, X)        (ga_docu/3_ga_modelling_motion.tex, eq:rcmt_pga_world)
    //
    // In 2D PGA the twist is a vector (twist2dp) and rcmt(vec, vec) -> vec; the argument
    // ORDER matters: rcmt(V, X) == -rcmt(X, V).
    static vec2dp velocity_field(twist2dp const& V, vec2dp const& X)
    {
        return rcmt(V, X);
    }

    // Acceleration field at point X of a rigid body with velocity twist V and
    // acceleration twist A (5_ga_modelling_physics.tex, "Moving coordinate systems"):
    //
    //   a(X) = rcmt(A, X)              [frame/Euler (alpha x r) + origin acceleration]
    //        + rcmt(V, rcmt(V, X))     [centripetal:  rcmt(Omega, rcmt(Omega, r))]
    static vec2dp accel_field(twist2dp const& V, twist2dp const& A, vec2dp const& X)
    {
        return rcmt(A, X) + rcmt(V, rcmt(V, X));
    }

    // Velocity of a world-space point X rigidly attached to frame idx. (The
    // dimension-neutral rcmt form is a later refinement; this explicit field is exact.)
    vec2dp point_velocity(vec2dp const& X_world, size_t idx)
    {
        return velocity_field(twist_world(idx), X_world);
    }

    vec2dp point_velocity(vec2dp const& X_world, std::string const& frame_name)
    {
        return point_velocity(X_world, index_of(frame_name));
    }

    // Velocity of a point P that is rigidly attached to frame `src` (given in src's
    // coordinates), as MEASURED BY an observer riding frame `obs`, expressed in obs's
    // coordinates. This is the rotating-frame transport theorem: the relative velocity is
    //
    //   v_rel = velocity_field(V_src - V_obs, P_world)   (a free vector, in world)
    //
    // then rotated into obs. Subtracting the velocity the obs-frame imparts at P (its
    // origin motion plus its rotation) means an observer reads zero for any point fixed
    // in their own frame, and reads the absolute world velocity when obs is the inertial
    // root.
    vec2dp relative_point_velocity(vec2dp const& P_src, size_t src, size_t obs)
    {
        auto const P_world = move2dp(P_src, get_pos_trafo(src, 0));
        auto const v_rel_world =
            velocity_field(twist_world(src) - twist_world(obs), P_world);
        return move2dp(v_rel_world, get_pos_trafo(0, obs)); // express in obs coordinates
    }

    vec2dp relative_point_velocity(vec2dp const& P_src, std::string const& src_name,
                                   std::string const& obs_name)
    {
        return relative_point_velocity(P_src, index_of(src_name), index_of(obs_name));
    }

    // World-frame acceleration twist A = dV/dt of frame idx (angular acceleration in
    // A.z).
    twist2dp accel_twist_world(size_t idx) { return world_VA(idx).A; }

    twist2dp accel_twist_world(std::string const& frame_name)
    {
        return accel_twist_world(index_of(frame_name));
    }

    // Acceleration of a world-space point X rigidly attached to frame idx. Planar rigid-
    // body acceleration field: the linear field of the acceleration twist A, plus the
    // centripetal term -w^2 * (X - centre) carried by the velocity twist V = (a, b, w):
    //
    //   a(X) = field(A, X) + w * (V.x - w*X.x, V.y - w*X.y)
    //
    // (field(A,X) supplies the tangential alpha x r and linear-acceleration parts).
    vec2dp point_acceleration(vec2dp const& X_world, size_t idx)
    {
        auto const va = world_VA(idx);
        return accel_field(va.V, va.A, X_world);
    }

    vec2dp point_acceleration(vec2dp const& X_world, std::string const& frame_name)
    {
        return point_acceleration(X_world, index_of(frame_name));
    }

    // Acceleration of a point P rigidly attached to frame `src` (given in src's
    // coordinates), as MEASURED BY an observer riding frame `obs`, expressed in obs's
    // coordinates. Rotating-frame transport theorem for acceleration:
    //
    //   a_rel = a_abs - a_transport - 2*w_obs x v_rel        (all in world, then -> obs)
    //
    // where a_abs/a_transport are the accel fields of src/obs at P, v_rel the world
    // relative velocity, w_obs the observer's angular rate, and 2*w_obs x v_rel the
    // Coriolis term. (The centrifugal -w_obs^2 r contribution is contained in
    // -a_transport.) Reads zero for a point fixed in obs; equals the absolute world
    // acceleration when obs is the root.
    vec2dp relative_point_acceleration(vec2dp const& P_src, size_t src, size_t obs)
    {
        auto const P_world = move2dp(P_src, get_pos_trafo(src, 0));
        auto const vs = world_VA(src);
        auto const vo = world_VA(obs);
        auto const a_abs = accel_field(vs.V, vs.A, P_world); // world accel of P
        auto const a_tr = accel_field(vo.V, vo.A, P_world);  // obs-frame accel at P
        auto const v_rel =
            velocity_field(vs.V - vo.V, P_world);      // world relative velocity
        auto const coriolis = 2.0 * rcmt(vo.V, v_rel); // 2 * rcmt(Omega_obs, v_rel)
        auto const a_rel_world = a_abs - a_tr - coriolis;
        return move2dp(a_rel_world, get_pos_trafo(0, obs)); // express in obs coordinates
    }

    vec2dp relative_point_acceleration(vec2dp const& P_src, std::string const& src_name,
                                       std::string const& obs_name)
    {
        return relative_point_acceleration(P_src, index_of(src_name), index_of(obs_name));
    }

    // True if frame idx is (instantaneously) an inertial frame: it neither rotates nor
    // accelerates relative to the inertial root. It becomes NON-inertial as soon as any
    // of these is non-zero in the world frame: angular velocity, angular acceleration, or
    // the linear acceleration of its origin. A constant non-zero linear velocity stays
    // inertial. tol is the absolute zero-tolerance (default loose enough to survive the
    // rounding accumulated by the motor transport).
    bool is_inertial_frame(size_t idx, value_t tol = value_t(1e-9))
    {
        auto const O_world = move2dp(O_2dp, get_pos_trafo(idx, 0));
        auto const a_O = point_acceleration(O_world, idx); // accel of the frame origin
        return std::abs(twist_world(idx).z) < tol          // no world angular velocity
               &&
               std::abs(accel_twist_world(idx).z) < tol // no world angular acceleration
               && std::abs(a_O.x) < tol &&
               std::abs(a_O.y) < tol; // origin not accelerating
    }

    bool is_inertial_frame(std::string const& frame_name, value_t tol = value_t(1e-9))
    {
        return is_inertial_frame(index_of(frame_name), tol);
    }

    // True if EVERY frame in the system is inertial (whole-system check).
    bool is_inertial_system(value_t tol = value_t(1e-9))
    {
        for (size_t idx = 0; idx < size(); ++idx)
            if (!is_inertial_frame(idx, tol)) return false;
        return true;
    }

    // read-only access to the stored relative twists of a frame (see twist2dp for the
    // encoding; used by the fmt formatter). Decode as omega = .z, v = (.y, -.x).
    twist2dp relative_twist(size_t idx) const { return rel_vtwist[idx]; }

    twist2dp relative_accel_twist(size_t idx) const { return rel_atwist[idx]; }

    twist2dp relative_twist(std::string const& frame_name) const
    {
        return relative_twist(index_of(frame_name));
    }

    twist2dp relative_accel_twist(std::string const& frame_name) const
    {
        return relative_accel_twist(index_of(frame_name));
    }

  private:

    // World velocity & acceleration twists of frame idx, propagated root -> idx by the
    // recursive Newton-Euler relations (twists transported to world by the adjoint):
    //
    //   V_i = V_parent + Ad(xi_i)
    //   A_i = A_parent + Ad(xidot_i) + [V_i, Ad(xi_i)]   (Coriolis / centrifugal
    //   coupling). The se(2) twist Lie bracket [.,.] is the regressive commutator
    //   rcmt(.,.) on twists.
    struct world_va2dp {
        twist2dp V; // velocity twist
        twist2dp A; // acceleration twist
    };

    world_va2dp world_VA(size_t idx)
    {
        // build the path root -> idx
        std::vector<size_t> path;
        for (size_t n = idx;; n = parent(n)) {
            path.push_back(n);
            if (parent(n) == n) break; // reached the root (self-parent)
        }
        std::reverse(path.begin(), path.end());

        twist2dp V{0.0, 0.0, 0.0};
        twist2dp A{0.0, 0.0, 0.0};
        for (size_t n : path) {
            if (parent(n) == n) continue; // root contributes nothing
            auto const M = get_pos_trafo(n, 0);
            auto const zeta = move2dp(rel_vtwist[n], M);    // Ad(xi_n)     world rel. vel
            auto const zetadot = move2dp(rel_atwist[n], M); // Ad(xidot_n)  world rel. acc
            V = V + zeta;
            A = A + zetadot + rcmt(V, zeta);
        }
        return {V, A};
    }
};


/////////////////////////////////////////////////////////////////////////////////////////
// dynamic_system2dp: the forces/inertia tier on top of kinematic_system2dp. Each body
// carries an inertia map and a mass; a frame's relative acceleration twist is no longer
// PRESCRIBED (as in the kinematic layer) but COMPUTED from the applied wrench via the
// se(2) Euler equation
//
//   Omega_dot = I^-1[ W - rcmt(Omega, I(Omega)) ]            (= compute_omega_dot)
//
// and integrated with RK4. The Coriolis/centrifugal coupling sits entirely in the
// regressive commutator rcmt(Omega, I(Omega)) -- "twists do not commute" -- which is the
// geometric-algebra showcase of this tier (see TODO/kinematics_explanation.md).
//
// Milestone 1 scope: independent FREE rigid bodies under gravity (no joint coupling or
// constraints yet -- the revolute-joint / articulated tier follows). The pose is evolved
// on the motor manifold; the integration state is the Lie-algebra pair (B, Omega).
/////////////////////////////////////////////////////////////////////////////////////////

// Rigid-body inertial properties of a frame (body frame, about the body origin = cm).
struct body2dp {
    Inertia2dp<value_t> I;     // inertia map: body twist -> body momentum
    Inertia2dp<value_t> I_inv; // its inverse (cached)
    value_t mass{0.0};         // total mass (gravity + energy)
};

// Build a body2dp for a uniform rectangular plate (width w along e1, height h along e2)
// of total mass m, with the body origin at the centre of mass.
inline body2dp make_plate_body(value_t m, value_t w, value_t h)
{
    auto const I = get_plate_inertia(m, w, h); // about cm (default pivot = body origin)
    return body2dp{I, get_inertia_inverse(I), m};
}

// Joint type connecting a body to its parent (the reduced-coordinate degrees of freedom).
// Both 1-DOF kinds run through the SAME code -- only the screw generator differs (the PGA
// unification of rotation and translation): M(q) = rest (x) exp(1/2 q * screw), relative
// twist = q-dot * screw, Jacobian = velocity_field(screw, .).
//
//   free      : unconstrained 3-DOF rigid body (Milestone 1); state in the base layer.
//
//   revolute  : 1-DOF hinge; screw = a FINITE point Q_b (z = 1); q rotates about Q_b.
//
//   prismatic : 1-DOF slider; screw = IDEAL point/direction (z = 0); q translates it.
enum class joint2dp { free, revolute, prismatic };

// Per-frame joint state (parallel to the body[] list). Meaningful for 1-DOF joints; the
// free body ignores it and uses the base relative twist + pose directly.
struct joint_state2dp {
    joint2dp type{joint2dp::free};
    // joint screw generator (a twist2dp), in the body frame: a finite point (z = 1) for a
    // revolute hinge, an ideal point / direction (z = 0) for a prismatic slider. The same
    // field drives both via the same exponential -- see joint2dp above.
    vec2dp screw_b{0.0, 0.0, 1.0};
    mvec2dp_u rest;     // body->parent motor at q = 0 (the reference pose)
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


class dynamic_system2dp : public kinematic_system2dp {

    // The optional closed-loop layer (ga_pga2dp_ops_constraints.hpp) composes a
    // dynamic_system2dp as its spanning tree and reuses this class's private assembly
    // seam (assemble_mass_bias, the joint screws, dof_joints/is_ancestor,
    // apply_joint_state) to build the loop-closure constraint Jacobian and the
    // constrained dynamics on top -- without those internals becoming public API.
    friend class closed_loop_system2dp;

    std::vector<body2dp> body;         // per-frame inertial properties (index = frame)
    std::vector<joint_state2dp> joint; // per-frame joint state (index = frame)
    vec2dp grav{0.0, -9.81, 0.0}; // uniform gravity field [world frame, z = 0 direction]

    // Optional time-varying applied wrench per frame (world frame), folded into tau as an
    // external generalised force. Sparse map (frame index -> wrench function of time); a
    // frame without an entry has no applied wrench. The simulation clock `time_` is the
    // argument passed at each RK4 sub-step (set transiently by coupled_step). Only the
    // jointed (assemble_mass_bias) path consumes it; free bodies see gravity only.
    using wrench_fn = std::function<bivec2dp(value_t)>;
    std::unordered_map<size_t, wrench_fn> wrench_;
    value_t time_{0.0}; // simulation clock [s], advanced by step()

  public:

    dynamic_system2dp() = default;

    // add a frame WITHOUT inertia (e.g. the inertial root); keeps body[]/joint[] in sync
    // with the base frame list. Mirrors the two base add_frame overloads.
    void add_frame(static_frame2dp const& f, kin_state2dp const& k,
                   size_t parent_idx = prev_frame)
    {
        kinematic_system2dp::add_frame(f, k, parent_idx);
        body.push_back(body2dp{});
        joint.push_back(joint_state2dp{});
    }

    void add_frame(static_frame2dp const& f, size_t parent_idx = prev_frame)
    {
        kinematic_system2dp::add_frame(f, parent_idx);
        body.push_back(body2dp{});
        joint.push_back(joint_state2dp{});
    }

    // add a free (3-DOF) dynamic rigid body: pose + inertial properties + initial state
    void add_body(static_frame2dp const& f, body2dp const& b,
                  kin_state2dp const& k = kin_state2dp{}, size_t parent_idx = prev_frame)
    {
        add_frame(f, k, parent_idx);
        body.back() = b;
    }

    // add a revolute-jointed body: a 1-DOF hinge about the body-fixed pivot Q_b. The
    // frame `f` provides the REST pose (q = 0); the body then rotates about Q_b by q. The
    // joint screw is the pivot point itself (a finite point, z = 1) -- exp(1/2 q Q_b) is
    // a rotation. The effective hinge inertia (Steiner) emerges automatically from the
    // spatial Jacobian in the forward dynamics; nothing extra to precompute.
    void add_revolute_body(static_frame2dp const& f, body2dp const& b,
                           vec2dp const& pivot_b, value_t phi0 = 0.0,
                           value_t omega0 = 0.0, size_t parent_idx = prev_frame)
    {
        add_screw_joint(f, b, joint2dp::revolute, pivot_b, phi0, omega0, parent_idx);
    }

    // add a prismatic-jointed body: a 1-DOF slider along the body-fixed unit direction
    // `dir`. The joint screw is the translation generator (an ideal point, z = 0):
    // exp(1/2 s * screw) is a pure translation by s along dir. The machinery is IDENTICAL
    // to the revolute joint -- only the generator differs (the PGA unification).
    void add_prismatic_body(static_frame2dp const& f, body2dp const& b, vec2dp const& dir,
                            value_t s0 = 0.0, value_t v0 = 0.0,
                            size_t parent_idx = prev_frame)
    {
        // translation generator along dir: to_twist(dir, 0) = (-dir.y, dir.x, 0)
        add_screw_joint(f, b, joint2dp::prismatic, vec2dp{-dir.y, dir.x, 0.0}, s0, v0,
                        parent_idx);
    }

    void set_gravity(vec2dp const& g) { grav = g; } // g is a direction (z = 0)
    vec2dp gravity() const { return grav; }

    // Attach a time-varying applied wrench (world frame) to frame `idx`, folded into tau
    // as an external generalised force projected onto each supporting joint screw. The
    // wrench is a bivec2dp (force/moment); for a pure force F through a point P use
    // wdg(P, F). Evaluated at each RK4 sub-step time. Pass an empty function to clear.
    void set_applied_wrench(size_t idx, wrench_fn fn)
    {
        if (fn) wrench_[idx] = std::move(fn);
        else wrench_.erase(idx);
    }

    value_t time() const { return time_; }  // simulation clock [s]
    void set_time(value_t t) { time_ = t; } // reset / seed the clock

    value_t joint_phi(size_t idx) const { return joint[idx].phi; }     // revolute angle
    value_t joint_omega(size_t idx) const { return joint[idx].omega; } // revolute rate

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

    // Current angular acceleration of revolute joint `idx`, from the COUPLED joint-space
    // forward dynamics at the present state (no integration).
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
    // velocity-product bias left after step()). Call this before reading accelerations
    // for force / acceleration visualisation.
    void sync_accelerations()
    {
        auto const rj = dof_joints();
        auto const qdd = forward_dynamics(rj); // zeroes rel_atwist internally
        for (size_t k = 0; k < rj.size(); ++k)
            set_accel_twist(rj[k], qdd[k] * joint[rj[k]].screw_b);
    }

    // Advance the system by dt. The revolute joints form a COUPLED chain integrated
    // together in their reduced (joint-angle) coordinates via the joint-space forward
    // dynamics; free bodies are integrated independently (Milestone 1). RK4 throughout.
    void step(value_t dt)
    {
        auto const rj = dof_joints();
        if (!rj.empty()) coupled_step(rj, dt); // uses time_ for sub-step wrench eval
        for (size_t i = 1; i < size(); ++i)
            if (joint[i].type == joint2dp::free && body[i].mass > 0.0)
                step_free_body(i, dt);
        time_ += dt; // advance the simulation clock (coupled_step restores it to t0)
    }

    // --- energy diagnostics (inertial / world frame) ---------------------------------

    // total kinetic energy: sum over bodies of  1/2 <V_body, I(V_body)>, the spatial
    // inertia-map quadratic form (= 1/2 m|v_cm|^2 + 1/2 I_cm omega^2), evaluated with the
    // body twist V_body = the world twist pulled into the body frame where the inertia
    // map lives. Same dimension-agnostic form used by the mass matrix and the 3D lift.
    value_t kinetic_energy()
    {
        value_t ke = 0.0;
        for (size_t i = 1; i < size(); ++i) {
            auto const Minv = rrev(get_pos_trafo(i, 0));       // world -> body i
            twist2dp const Vb = move2dp(twist_world(i), Minv); // body twist of body i
            ke += 0.5 * spatial_dot(Vb, body[i].I(Vb));
        }
        return ke;
    }

    // total gravitational potential energy: sum over bodies of  -m (g . r_cm)
    value_t potential_energy()
    {
        value_t pe = 0.0;
        for (size_t i = 1; i < size(); ++i) {
            vec2dp const cm_w = move2dp(O_2dp, get_pos_trafo(i, 0));
            pe += -body[i].mass * (grav.x * cm_w.x + grav.y * cm_w.y);
            // joint-spring potential 1/2 k (q - q0)^2 (zero unless a spring is attached)
            value_t const dq = joint[i].phi - joint[i].q_rest;
            pe += 0.5 * joint[i].stiffness * dq * dq;
        }
        return pe;
    }

    value_t total_energy() { return kinetic_energy() + potential_energy(); }

    // Joint-space mass matrix M(q) (n*n row-major, n = number of revolute joints) at the
    // current configuration. A diagnostic / showcase quantity: the reduced inertia of the
    // articulated system, with the identity  1/2 * qdot^T M(q) qdot == kinetic_energy().
    std::vector<value_t> mass_matrix()
    {
        auto const rj = dof_joints();
        size_t const n = rj.size();
        std::vector<twist2dp> S(n); // world joint screws (unit rate)
        std::vector<mvec2dp_u> Minv(
            n); // world -> body i motor (pulls a twist into body i)
        for (size_t i = 0; i < n; ++i) {
            auto const M = get_pos_trafo(rj[i], 0); // body i -> world
            S[i] = move2dp(joint[rj[i]].screw_b, M);
            Minv[i] = rrev(M);
        }
        std::vector<value_t> Mmat(n * n, 0.0);
        for (size_t i = 0; i < n; ++i) {
            auto const& I = body[rj[i]].I; // inertia map about body i's cm (body frame)
            for (size_t j = 0; j < n; ++j) {
                if (!is_ancestor(rj[j], rj[i])) continue;
                twist2dp const xj = move2dp(S[j], Minv[i]); // joint-j screw in body i
                for (size_t k = 0; k < n; ++k) {
                    if (!is_ancestor(rj[k], rj[i])) continue;
                    twist2dp const xk = move2dp(S[k], Minv[i]);
                    // inertia-map quadratic form: carries mass + angular term uniformly
                    // (no S.z split), and is the form that lifts unchanged to 3D.
                    Mmat[j * n + k] += spatial_dot(xj, I(xk));
                }
            }
        }
        return Mmat;
    }

  private:

    // Spatial (reciprocal / Klein) pairing of a velocity twist with a momentum or
    // wrench bivector -> scalar. This is the dimension-agnostic inertia-map quadratic
    // form: with mom = I(xi), the value <xi, I(xi)> equals 2x the kinetic energy of a
    // body moving with body twist xi -- the mass AND the angular term are carried
    // uniformly by the inertia map I, with no separate m|v|^2 + I_cm*w^2 split. This is
    // the form that lifts UNCHANGED to 3D (I becomes the 6x6 Inertia3dp, the pairing the
    // BiVec3dp screw reciprocal product). In PGA2DP the pairing is the regressive product
    //
    //   spatial_dot(xi, mom) = -rwdg(xi, mom) = xi.y*mom.x - xi.x*mom.y + xi.z*mom.z.
    static value_t spatial_dot(twist2dp const& xi, bivec2dp const& mom)
    {
        return -value_t(rwdg(xi, mom));
    }

    // RK4-integrate one free rigid body (frame idx) over dt under gravity. The
    // integration state is the Lie-algebra pair (B, Omega): B is the relative generator
    // accumulated from the current relative pose M0 (so M(t) = M0 (x) exp(1/2 B)), Omega
    // the body twist. dB/dt = Omega; dOmega/dt = I^-1[ W_body - rcmt(Omega, I(Omega)) ].
    void step_free_body(size_t idx, value_t dt)
    {
        auto const M0 = rrev(step_pos_trafo(idx)); // current body -> parent motor
        auto const& I = body[idx].I;
        auto const& I_inv = body[idx].I_inv;
        value_t const m = body[idx].mass;

        // body-frame twist rate from the gravity wrench acting at the cm (= body origin)
        auto omega_dot = [&](vec2dp const& B, twist2dp const& Om) -> twist2dp {
            auto const M = rgpr(M0, exp(0.5 * B));             // pose at this stage
            auto const W_w = wdg(move2dp(O_2dp, M), m * grav); // gravity wrench (world)
            auto const W_b = move2dp(W_w, rrev(M)); // pulled into the body frame
            return compute_omega_dot(I_inv, W_b, Om, I);
        };

        // RK4 (shared rk4_step) on the Lie-algebra pair u = (B, Omega): dB/dt = Omega,
        // dOmega/dt = omega_dot(B, Omega). B starts at 0 (M(t) = M0 (x) exp(1/2 B)).
        std::array<twist2dp, 2> u_mem{twist2dp{}, relative_twist(idx)};
        std::array<twist2dp, 4> uh_mem{};
        std::array<twist2dp, 2> rhs_mem{};
        auto u = std::mdspan<twist2dp, std::dextents<size_t, 1>>(u_mem.data(), 2);
        auto uh = std::mdspan<twist2dp, std::dextents<size_t, 2>>(uh_mem.data(), 2, 2);
        auto const rhs =
            std::mdspan<twist2dp const, std::dextents<size_t, 1>>(rhs_mem.data(), 2);
        for (size_t s = 1; s <= 4; ++s) {
            rhs_mem[0] = u[1];                  // dB/dt = Omega
            rhs_mem[1] = omega_dot(u[0], u[1]); // dOmega/dt
            rk4_step(u, uh, rhs, dt, s);
        }

        // decode the evolved pose (origin, phi) and write pose + twist into the base
        // layer
        set_pose(idx, pose2dp_from_motor(rgpr(M0, exp(0.5 * u[0]))));
        set_twist(idx, u[1]);
    }

    // shared constructor for the 1-DOF screw joints (revolute / prismatic): the ONLY
    // difference between them is the screw generator passed in.
    void add_screw_joint(static_frame2dp const& f, body2dp const& b, joint2dp type,
                         vec2dp const& screw_b, value_t q0, value_t qdot0,
                         size_t parent_idx)
    {
        add_frame(f, parent_idx);
        size_t const idx = size() - 1;
        body[idx] = b;
        joint_state2dp js;
        js.type = type;
        js.screw_b = screw_b;
        js.rest = rrev(step_pos_trafo(idx)); // body->parent motor at the rest pose
        js.phi = q0;
        js.omega = qdot0;
        joint[idx] = js;
        apply_joint_state(idx); // write q0/qdot0 into the base pose + twist
    }

    // body->parent motor at generalised coordinate q: M(q) = rest (x) exp(1/2 q * screw).
    // The SAME exponential builds a rotation (revolute, screw = finite point) or a
    // translation (prismatic, screw = ideal point); the operators do not change.
    mvec2dp_u joint_motor(size_t idx, value_t q) const
    {
        return rgpr(joint[idx].rest, exp(0.5 * q * joint[idx].screw_b));
    }

    // Write the joint state (q, q-dot) into the base pose + relative twist, so all
    // kinematic queries (get_pos_trafo, twist_world, point_velocity, energy) stay valid.
    void apply_joint_state(size_t idx)
    {
        auto const M = joint_motor(idx, joint[idx].phi); // body->parent motor at q
        set_pose(idx, pose2dp_from_motor(M));
        set_twist(idx,
                  joint[idx].omega * joint[idx].screw_b); // rel twist = q-dot * screw
    }

    // Frame indices of all 1-DOF joints (revolute or prismatic) -- the generalised
    // coords.
    std::vector<size_t> dof_joints() const
    {
        std::vector<size_t> rj;
        for (size_t i = 1; i < size(); ++i)
            if (joint[i].type == joint2dp::revolute ||
                joint[i].type == joint2dp::prismatic)
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
    // 1-DOF joint chain `rj`, by virtual work over the bodies (all quantities GA-native):
    //
    //   M[j][k] = sum_i  spatial_dot( S_j^body_i , I_i( S_k^body_i ) )  over bodies i
    //             having BOTH joints j,k as ancestors  (the spatial inertia-map form);
    //   RHS[j]  = sum_i m_i ( vcm_i(S_j) . g  -  vcm_i(S_j) . b_cm_i )  over bodies i
    //   having
    //             joint j as ancestor  (gravity generalised force minus Coriolis bias).
    //
    // S_j = move2dp(Q_j, M_{j->world}) is the world joint screw; S_j^body_i =
    // move2dp(S_j, rrev(M_i)) transports it into body i's frame (where its inertia map
    // I_i lives); vcm_i(S) = velocity_field(S, cm_i) is the cm velocity per unit joint
    // rate (a spatial-Jacobian column); b_cm_i is the velocity-product
    // (Coriolis/centripetal) cm acceleration at q-ddot = 0. Pre: the joint state (phi,
    // omega) is already applied to the base pose + relative twist.
    //
    // SIDE EFFECT: this runs the bias pass, zeroing the chain's relative accel twists
    // (rel_atwist) so the world accel queries return only the q-ddot-independent part.
    //
    // Split out of forward_dynamics() as the reuse seam for the closed-loop layer
    // (ga_pga2dp_ops_constraints.hpp): the same spatial-Jacobian columns (velocity_field)
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
        // queries return only the q-ddot-independent Coriolis/centripetal part.
        for (size_t k = 0; k < n; ++k)
            set_accel_twist(rj[k], twist2dp{0.0, 0.0, 0.0});

        std::vector<vec2dp> S(n), cm(n), bcm(n);
        std::vector<mvec2dp_u> Minv(
            n); // world -> body i motor (pulls a twist into body i)
        std::vector<value_t> mass(n);
        for (size_t i = 0; i < n; ++i) {
            size_t const fi = rj[i];
            auto const M = get_pos_trafo(fi, 0);
            cm[i] = move2dp(O_2dp, M);
            bcm[i] = point_acceleration(cm[i], fi); // bias cm accel (rel_atwist = 0)
            S[i] = move2dp(joint[fi].screw_b, M);   // world joint screw (unit rate)
            Minv[i] = rrev(M);
            mass[i] = body[fi].mass;
        }

        std::vector<value_t> Mmat(n * n, 0.0), RHS(n, 0.0);
        for (size_t j = 0; j < n; ++j) {
            for (size_t i = 0; i < n; ++i) { // contribution of body i to coordinate j
                if (!is_ancestor(rj[j], rj[i])) continue;
                vec2dp const vj = velocity_field(S[j], cm[i]); // rcmt(S_j, cm_i)
                RHS[j] += mass[i] * (vj.x * grav.x + vj.y * grav.y) -
                          mass[i] * (vj.x * bcm[i].x + vj.y * bcm[i].y);
                auto const& I =
                    body[rj[i]].I; // inertia map about body i's cm (body frame)
                twist2dp const xj = move2dp(S[j], Minv[i]); // joint-j screw in body i
                for (size_t k = 0; k < n; ++k) {
                    if (!is_ancestor(rj[k], rj[i])) continue;
                    twist2dp const xk = move2dp(S[k], Minv[i]);
                    // inertia-map quadratic form: carries mass + angular term uniformly
                    // (no S.z split), and is the form that lifts unchanged to 3D.
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
            bivec2dp const W = fn(time_);
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
        std::vector<value_t> u_mem(2 * n), uh_mem(2 * 2 * n), rhs_mem(2 * n);
        for (size_t k = 0; k < n; ++k) {
            u_mem[k] = joint[rj[k]].phi;
            u_mem[n + k] = joint[rj[k]].omega;
        }
        auto u = std::mdspan<value_t, std::dextents<size_t, 1>>(u_mem.data(), 2 * n);
        auto uh = std::mdspan<value_t, std::dextents<size_t, 2>>(uh_mem.data(), 2, 2 * n);
        auto const rhs =
            std::mdspan<value_t const, std::dextents<size_t, 1>>(rhs_mem.data(), 2 * n);

        // apply the current u to the joint state, then solve forward dynamics -> rhs
        auto apply_u = [&] {
            for (size_t k = 0; k < n; ++k) {
                joint[rj[k]].phi = u[k];
                joint[rj[k]].omega = u[n + k];
                apply_joint_state(rj[k]);
            }
        };
        // thread the RK4 sub-step time into time_ so a time-varying applied wrench is
        // sampled at the correct stage time; restore the clock on exit (step() advances
        // it by dt). t_i + {0, dt/2, dt/2, dt} for the four stages.
        value_t const t0 = time_;
        for (size_t s = 1; s <= 4; ++s) {
            time_ = rk4_get_time(t0, dt, s - 1);
            apply_u();
            auto const qdd = forward_dynamics(rj);
            for (size_t k = 0; k < n; ++k) {
                rhs_mem[k] = u[n + k];   // dphi/dt = omega
                rhs_mem[n + k] = qdd[k]; // domega/dt = q-ddot
            }
            rk4_step(u, uh, rhs, dt, s);
        }
        time_ = t0;
        apply_u(); // write the integrated state back into the joints + kinematic layers
    }
};

} // namespace hd::ga::pga

// value_t convenience aliases (inertia2dp, ...) live in ga_usr_types_physics.hpp, next to
// the other physics convenience types.
