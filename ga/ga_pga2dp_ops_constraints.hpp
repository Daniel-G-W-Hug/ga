#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Closed-loop (parallel mechanism) layer for PGA2DP. This is a SEPARATE, ADDITIVE
// capability on top of the open-chain dynamic_system2dp: a user who only needs open
// chains includes ga_pga2dp_ops_physics.hpp and pays no complexity. closed_loop_system2dp
// HAS-A dynamic_system2dp (composition) as its spanning tree and closes kinematic loops
// by registering loop-closure constraints between existing tree frames, then solving the
// resulting constrained problem. See TODO/closed_loop_system_consideration.md for the
// design and ga_docu/5_ga_modelling_physics.tex ("Reduced vs. maximal coordinates").
//
// Phase 1 (this file, position level): the point-coincidence constraint, the residual
// g(q), the constraint Jacobian G = dg/dq (reusing the spatial-Jacobian columns
// velocity_field(S_j, P) of the open-loop assembly), and assemble() -- solve g(q) = 0 by
// Newton for a consistent closed configuration. Velocity (Phase 2) and dynamic (Phase 3)
// closed loops follow.

#include "detail/ga_solver.hpp"      // hd::ga::lu_solve (shared LU)
#include "ga_pga2dp_ops_physics.hpp" // dynamic_system2dp (the reused spanning tree)
#include "ga_value_t.hpp"            // value_t

#include <cmath>     // std::abs
#include <cstddef>   // size_t
#include <stdexcept> // std::runtime_error
#include <string>
#include <vector>


namespace hd::ga::pga {

// Loop-closure constraint kinds. Phase 1 implements point coincidence only (a planar pin
// / a spatial spherical-joint attachment); point-on-line, distance / rigid-link length
// and parallel/perpendicular follow (all expressible as vanishing wdg/rwdg/dot
// expressions).
enum class constraint2dp { coincidence };

// Loop-closure constraint descriptor: anchor point anchor_a (given in frame_a's
// coordinates) must satisfy `type` against anchor point anchor_b (in frame_b's
// coordinates). For coincidence the two anchors must coincide in world coordinates -- the
// loop-closure condition for a four-bar pin or a Stewart-Gough leg pinned to the
// platform.
//
// A pure-data aggregate (all public fields): it auto-binds in ga_py on a plain bindgen
// regeneration (give it an fmt formatter so the generated __str__/__format__ work, and
// cover it with a dedicated ga_py test).
struct loop_constraint2dp {
    size_t frame_a;  // tree frame holding anchor a
    vec2dp anchor_a; // anchor point in frame_a coordinates (unitized, z = 1)
    size_t frame_b;  // tree frame holding anchor b
    vec2dp anchor_b; // anchor point in frame_b coordinates (unitized, z = 1)
    constraint2dp type{constraint2dp::coincidence};
};


// A closed-loop / parallel mechanism: an open-chain spanning tree (the reused
// dynamic_system2dp) plus loop-closure constraints between tree frames. Build the tree
// through the familiar forwarded API, then close loops with add_loop_constraint(...).
class closed_loop_system2dp {

    dynamic_system2dp tree_;                // the open-chain spanning tree (reused as-is)
    std::vector<loop_constraint2dp> loops_; // the extra closure edges

  public:

    closed_loop_system2dp() = default;

    // --- build the spanning tree (forward the familiar dynamic_system2dp builders)
    // ------

    void add_frame(static_frame2dp const& f,
                   size_t parent_idx = dynamic_system2dp::prev_frame)
    {
        tree_.add_frame(f, parent_idx);
    }

    void add_body(static_frame2dp const& f, body2dp const& b,
                  kin_state2dp const& k = kin_state2dp{},
                  size_t parent_idx = dynamic_system2dp::prev_frame)
    {
        tree_.add_body(f, b, k, parent_idx);
    }

    void add_revolute_body(static_frame2dp const& f, body2dp const& b,
                           vec2dp const& pivot_b, value_t phi0 = 0.0,
                           value_t omega0 = 0.0,
                           size_t parent_idx = dynamic_system2dp::prev_frame)
    {
        tree_.add_revolute_body(f, b, pivot_b, phi0, omega0, parent_idx);
    }

    void add_prismatic_body(static_frame2dp const& f, body2dp const& b, vec2dp const& dir,
                            value_t s0 = 0.0, value_t v0 = 0.0,
                            size_t parent_idx = dynamic_system2dp::prev_frame)
    {
        tree_.add_prismatic_body(f, b, dir, s0, v0, parent_idx);
    }

    size_t index_of(std::string const& frame_name) const
    {
        return tree_.index_of(frame_name);
    }

    // --- close loops
    // --------------------------------------------------------------------

    // register a loop-closure constraint between two existing tree frames
    void add_loop_constraint(loop_constraint2dp const& c) { loops_.push_back(c); }

    // drive a joint to a generalised coordinate q (set q, refresh the kinematic state).
    // For the driver / independent joints of a closed loop -- the dependent ones are
    // produced by assemble().
    void set_joint(size_t joint_frame, value_t q)
    {
        tree_.joint[joint_frame].phi = q;
        tree_.apply_joint_state(joint_frame);
    }

    value_t joint_phi(size_t joint_frame) const { return tree_.joint_phi(joint_frame); }

    size_t loop_count() const { return loops_.size(); }
    loop_constraint2dp const& loop(size_t c) const { return loops_[c]; }

    // access the underlying spanning tree for diagnostics shared with the open-loop tier
    // (get_pos_trafo, total_energy, joint_phi, ...) and for further configuration
    dynamic_system2dp& system() { return tree_; }
    dynamic_system2dp const& system() const { return tree_; }

    // --- constraint evaluation (position level)
    // -----------------------------------------

    // Residual g(q): for each coincidence constraint, the unitized world-coordinate
    // difference P_a - P_b (2 scalar equations per constraint in 2D). g == 0 <=> every
    // loop is closed. Length 2 * loop_count().
    std::vector<value_t> residual()
    {
        std::vector<value_t> g(2 * loops_.size());
        for (size_t c = 0; c < loops_.size(); ++c) {
            auto const& lc = loops_[c];
            vec2dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec2dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            g[2 * c + 0] = Pa.x - Pb.x;
            g[2 * c + 1] = Pa.y - Pb.y;
        }
        return g;
    }

    // infinity norm of the residual (the closure error ‖g‖)
    value_t residual_norm()
    {
        value_t n = 0.0;
        for (value_t const gi : residual())
            n = std::max(n, std::abs(gi));
        return n;
    }

    // Assemble a consistent closed configuration: hold the `driven` joints fixed and
    // solve the remaining (dependent) joint coordinates so that g(q) = 0, by Newton
    // iteration on the constraint Jacobian G (§3.4: G's columns are the relative partial
    // velocities velocity_field(S_j, P_a) - velocity_field(S_j, P_b), reused from the
    // open-loop spatial Jacobian). Returns the achieved residual norm ‖g‖. Throws if the
    // dependent Jacobian is singular or Newton does not converge within max_iter.
    //
    // For a 1-DOF four-bar (3 revolute joints, 2 coincidence equations) drive one joint;
    // the 2 dependent joints are then solved by the square Newton step.
    value_t assemble(std::vector<size_t> const& driven = {}, value_t tol = value_t(1e-12),
                     size_t max_iter = 50)
    {
        std::vector<size_t> const dep = dependent_joints(driven);

        for (size_t it = 0; it <= max_iter; ++it) {
            std::vector<value_t> const g = residual();
            value_t gnorm = 0.0;
            for (value_t const gi : g)
                gnorm = std::max(gnorm, std::abs(gi));
            if (gnorm < tol) return gnorm;
            if (it == max_iter) break;

            // Newton step: solve G_dep * delta = -g for the dependent joint increments
            std::vector<value_t> const G = constraint_jacobian(dep);
            std::vector<value_t> const delta = newton_step(G, g, dep.size());
            for (size_t k = 0; k < dep.size(); ++k) {
                tree_.joint[dep[k]].phi += delta[k];
                tree_.apply_joint_state(dep[k]);
            }
        }
        throw std::runtime_error(
            std::string(
                "closed_loop_system2dp::assemble: Newton did not converge within ") +
            std::to_string(max_iter) + std::string(" iterations (residual ") +
            std::to_string(residual_norm()) + std::string(")"));
    }

  private:

    // world-coordinate, unitized (z = 1) position of an anchor point of a tree frame
    vec2dp anchor_world(size_t frame, vec2dp const& anchor)
    {
        vec2dp const P = move2dp(anchor, tree_.get_pos_trafo(frame, 0));
        return vec2dp{P.x / P.z, P.y / P.z, 1.0};
    }

    // the dof joints NOT in `driven` -- the coordinates assemble() solves for
    std::vector<size_t> dependent_joints(std::vector<size_t> const& driven) const
    {
        std::vector<size_t> dep;
        for (size_t const j : tree_.dof_joints()) {
            bool is_driven = false;
            for (size_t const d : driven)
                if (d == j) {
                    is_driven = true;
                    break;
                }
            if (!is_driven) dep.push_back(j);
        }
        return dep;
    }

    // Constraint Jacobian G (m x ndep, row-major), m = 2 * loop_count(), columns indexed
    // by the dependent joints `dep`. Column k for constraint c is the relative partial
    // velocity of the two anchors w.r.t. joint dep[k]:
    //
    //     dP_a/dq_j - dP_b/dq_j = velocity_field(S_j, P_a) - velocity_field(S_j, P_b)
    //
    // (each term present only when joint j is an ancestor of the anchor's frame), where
    // S_j = move2dp(screw_b, M_{j->world}) is the world joint screw -- the SAME column
    // that builds the open-loop mass matrix and bias forces (assemble_mass_bias).
    std::vector<value_t> constraint_jacobian(std::vector<size_t> const& dep)
    {
        size_t const m = 2 * loops_.size();
        size_t const ndep = dep.size();
        std::vector<value_t> G(m * ndep, 0.0);

        // precompute the world joint screws of the dependent joints
        std::vector<twist2dp> S(ndep);
        for (size_t k = 0; k < ndep; ++k)
            S[k] = move2dp(tree_.joint[dep[k]].screw_b, tree_.get_pos_trafo(dep[k], 0));

        for (size_t c = 0; c < loops_.size(); ++c) {
            auto const& lc = loops_[c];
            vec2dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec2dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            for (size_t k = 0; k < ndep; ++k) {
                size_t const j = dep[k];
                vec2dp dPa{0.0, 0.0, 0.0};
                vec2dp dPb{0.0, 0.0, 0.0};
                if (tree_.is_ancestor(j, lc.frame_a))
                    dPa = dynamic_system2dp::velocity_field(S[k], Pa);
                if (tree_.is_ancestor(j, lc.frame_b))
                    dPb = dynamic_system2dp::velocity_field(S[k], Pb);
                G[(2 * c + 0) * ndep + k] = dPa.x - dPb.x;
                G[(2 * c + 1) * ndep + k] = dPa.y - dPb.y;
            }
        }
        return G;
    }

    // Solve G * delta = -g for the Newton increment (G is m x ndep, row-major). Three
    // regimes, all routed through the shared square LU solver (detail/ga_solver.hpp):
    //   ndep == m : square        -> solve directly
    //   ndep >  m : underdetermined -> minimum-norm  delta = G^T (G G^T)^-1 (-g)
    //   ndep <  m : overdetermined  -> least-squares delta = (G^T G)^-1 G^T (-g)
    // The four-bar is the square case; the other two keep assemble() well-posed for
    // redundant / over-constrained loops.
    static std::vector<value_t> newton_step(std::vector<value_t> const& G,
                                            std::vector<value_t> const& g, size_t ndep)
    {
        size_t const m = g.size();
        std::vector<value_t> rhs(m);
        for (size_t i = 0; i < m; ++i)
            rhs[i] = -g[i];

        if (ndep == m) {
            return hd::ga::lu_solve(G, rhs, m); // square G
        }
        if (ndep > m) {
            // minimum-norm: solve (G G^T) y = -g (m x m), then delta = G^T y
            std::vector<value_t> GGt(m * m, 0.0);
            for (size_t i = 0; i < m; ++i)
                for (size_t j = 0; j < m; ++j) {
                    value_t s = 0.0;
                    for (size_t k = 0; k < ndep; ++k)
                        s += G[i * ndep + k] * G[j * ndep + k];
                    GGt[i * m + j] = s;
                }
            std::vector<value_t> const y = hd::ga::lu_solve(GGt, rhs, m);
            std::vector<value_t> delta(ndep, 0.0);
            for (size_t k = 0; k < ndep; ++k) {
                value_t s = 0.0;
                for (size_t i = 0; i < m; ++i)
                    s += G[i * ndep + k] * y[i];
                delta[k] = s;
            }
            return delta;
        }
        // overdetermined: normal equations (G^T G) delta = G^T (-g) (ndep x ndep)
        std::vector<value_t> GtG(ndep * ndep, 0.0), Gtr(ndep, 0.0);
        for (size_t a = 0; a < ndep; ++a) {
            for (size_t i = 0; i < m; ++i)
                Gtr[a] += G[i * ndep + a] * rhs[i];
            for (size_t b = 0; b < ndep; ++b) {
                value_t s = 0.0;
                for (size_t i = 0; i < m; ++i)
                    s += G[i * ndep + a] * G[i * ndep + b];
                GtG[a * ndep + b] = s;
            }
        }
        return hd::ga::lu_solve(GtG, Gtr, ndep);
    }
};

} // namespace hd::ga::pga
