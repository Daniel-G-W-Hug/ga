#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_solver.hpp" // hd::ga::lu_decomp / lu_backsubs / det
#include "ga_pga2dp_ops.hpp"
#include "ga_value_t.hpp" // for value_t used in convenience type alias

#include <algorithm> // std::find, std::min, std::max
#include <array>
#include <cmath>  // std::abs
#include <limits> // std::numeric_limits
#include <mdspan>
#include <stdexcept>
#include <string>
#include <unordered_map> // std::unordered_map (frame name -> index)
#include <vector>


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
// I[2,2] = m*(w^2+h^2)/12 : moment of inertia about centroid (classical rectangle
// formula)
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

/////////////////////////////////////////////////////////////////////////////////////////
// static_frame2dp / static_system2dp: a tree of right-handed coordinate frames, each
// posed RELATIVE to its parent. Provides rigid coordinate transformations between frames.
/////////////////////////////////////////////////////////////////////////////////////////

struct pose2dp {
    vec2dp origin; // origin of frame in parent coordinates
    value_t phi;   // orientation of frame vs. parent coordinates, [phi]: rad
};

class static_frame2dp {

    // Static basis frame, a coordinate system (cs) w/o change over time.
    // Cs parent frame is assumed to be right-handed with e1 in x-direction, e2 rotated
    // by 90° anti-clockwise vs. e1 defining positive phi direction in direction of e12

    std::string name; // display name (default: 'I' = inertial frame)

    // pose vs. parent coordinate system
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
    //   rf[0] (root) -> rf[1] -> rf[2] -> ...
    // Pass an explicit parent_idx to branch the tree off an earlier frame instead. A
    // frame's pose is always interpreted RELATIVE to its parent (see static_frame2dp).
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
    // To be used as "p_to = move2dp(p_from, M);"
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
    // "pose vs. parent coordinate system").
    //
    // Translation and rotation are COMPOSED, not summed as generators: they do not
    // commute, so exp(-0.5*(t+r)) != exp(-0.5*r) (x) exp(-0.5*t) (Baker-Campbell-
    // Hausdorff), and the single-exp form is a screw about a wrong centre (it does not
    // even map the child origin to (0,0)). The composition reduces correctly to a pure
    // translation (rel.phi == 0) or pure rotation (rel.origin == parent origin).
    mvec2dp_u step_pos_trafo(size_t child_idx) const
    {
        auto const rel = vfr[child_idx].get_pose(); // pose relative to parent frame
        auto M_rot = exp(-0.5 * vec2dp(0.0, 0.0, rel.phi)); // rotate about parent origin
        auto M_tra =
            exp(-0.5 * vec2dp(-rel.origin.y, rel.origin.x, 0.0)); // translate by rel
        return rgpr(M_rot, M_tra);                                // parent -> child
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
// point/vector. In the right-handed plane (omega > 0 is counter-clockwise) the components
// encode the pga2dp motor generator
//     twist2dp(.x, .y, .z) = (-v_y, v_x, omega)
// so that exp(0.5 * twist) is the motor and move2dp(twist, M) is its adjoint. Decode as
// omega = .z and v = (.y, -.x). The alias documents this intent at every signature and
// adds no overloads (it IS vec2dp). In pga3dp the corresponding twist is a genuine
// BiVec3dp (grade 2) -- see the 2D->3D notes.
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

    // Advance the system by one time step dt (in-place evolution). Each non-root frame's
    // relative pose is evolved on the motor manifold by the body-frame relative twist,
    //   P_new = P (x) exp(0.5 * B_rel * dt)         [P = body->parent motor]
    // which is EXACT for a constant twist (N steps of dt == one step of N*dt), and its
    // relative velocity twist is ramped by the relative acceleration (explicit Euler):
    //   B_rel += B_accel * dt.
    // The root frame is left unchanged. Prescribed-motion kinematics only -- no forces.
    void step(value_t dt)
    {
        for (size_t i = 1; i < size(); ++i) { // frame 0 is the root (left unchanged)
            // evolve the body->parent relative motor by the relative twist over dt
            auto const P = rrev(step_pos_trafo(i));                  // T(origin) (x) R(phi)
            auto const P_new = rgpr(P, exp(0.5 * rel_vtwist[i] * dt));
            // decode the new relative pose (origin, phi) and store it
            auto const o = move2dp(O_2dp, P_new);                    // child origin in parent
            auto const e1 = move2dp(vec2dp{1.0, 0.0, 0.0}, P_new);   // rotated e1 direction
            set_pose(i, vec2dp{o.x / o.z, o.y / o.z, 1.0}, std::atan2(e1.y, e1.x));
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

    // Velocity field of a twist V = (a, b, w): the velocity of the point at world
    // position X moving rigidly with V, i.e. v(X) = v0 + w * (zhat x X) with v0 = (b,
    // -a).
    static vec2dp velocity_field(twist2dp const& V, vec2dp const& X)
    {
        return vec2dp(V.y - V.z * X.y, -V.x + V.z * X.x, 0.0);
    }

    // Acceleration field at point X of a rigid body with velocity twist V and acceleration
    // twist A: the linear field of A (tangential alpha x r + origin acceleration) plus the
    // centripetal term carried by V = (a, b, w):  a(X) = field(A, X) + w*(V.x - w*X.x, ...).
    static vec2dp accel_field(twist2dp const& V, twist2dp const& A, vec2dp const& X)
    {
        auto const cent = vec2dp(V.z * (V.x - V.z * X.x), V.z * (V.y - V.z * X.y), 0.0);
        return velocity_field(A, X) + cent;
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
    //   v_rel = velocity_field(V_src - V_obs, P_world)   (a free vector, in world)
    // then rotated into obs. Subtracting the velocity the obs-frame imparts at P (its
    // origin motion plus its rotation) means an observer reads zero for any point fixed in
    // their own frame, and reads the absolute world velocity when obs is the inertial root.
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
    //   a(X) = field(A, X) + w * (V.x - w*X.x, V.y - w*X.y)
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
    //   a_rel = a_abs - a_transport - 2*w_obs x v_rel        (all in world, then -> obs)
    // where a_abs/a_transport are the accel fields of src/obs at P, v_rel the world relative
    // velocity, w_obs the observer's angular rate, and 2*w_obs x v_rel the Coriolis term.
    // (The centrifugal -w_obs^2 r contribution is contained in -a_transport.) Reads zero for
    // a point fixed in obs; equals the absolute world acceleration when obs is the root.
    vec2dp relative_point_acceleration(vec2dp const& P_src, size_t src, size_t obs)
    {
        auto const P_world = move2dp(P_src, get_pos_trafo(src, 0));
        auto const vs = world_VA(src);
        auto const vo = world_VA(obs);
        auto const a_abs = accel_field(vs.V, vs.A, P_world);       // world accel of P
        auto const a_tr = accel_field(vo.V, vo.A, P_world);        // obs-frame accel at P
        auto const v_rel = velocity_field(vs.V - vo.V, P_world);   // world relative velocity
        auto const coriolis = 2.0 * vo.V.z * vec2dp(-v_rel.y, v_rel.x, 0.0);
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

    // se(2) twist Lie bracket in the twist2dp generator encoding (used for the Coriolis
    // coupling in the acceleration propagation; angular part of [.,.] is zero)
    static twist2dp bracket(twist2dp const& a, twist2dp const& b)
    {
        return vec2dp(b.z * a.y - a.z * b.y, a.z * b.x - b.z * a.x, 0.0);
    }

    // World velocity & acceleration twists of frame idx, propagated root -> idx by the
    // recursive Newton-Euler relations (twists transported to world by the adjoint):
    //   V_i = V_parent + Ad(xi_i)
    //   A_i = A_parent + Ad(xidot_i) + [V_i, Ad(xi_i)]   (Coriolis / centrifugal
    //   coupling)
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
            A = A + zetadot + bracket(V, zeta);
        }
        return {V, A};
    }
};

} // namespace hd::ga::pga

// value_t convenience aliases (inertia2dp, ...) live in ga_usr_types_physics.hpp, next to
// the other physics convenience types.
