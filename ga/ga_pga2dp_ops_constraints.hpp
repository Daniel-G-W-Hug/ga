#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Closed-loop (parallel mechanism) layer for PGA2DP. A SEPARATE, ADDITIVE capability on
// top of the open-chain dynamic_system2dp: a user who only needs open chains includes
// ga_pga2dp_ops_mechanics.hpp and pays no complexity. closed_loop_system2dp HAS-A
// dynamic_system2dp (composition) as its spanning tree and closes kinematic loops by
// registering loop-closure constraints between existing tree frames, then solving the
// resulting constrained problem. (Background: ga_docu/5_ga_modelling_mechanics.tex,
// "Reduced vs. maximal coordinates".)
//
// Everything is built on the constraint residual g(q) (the closure error) and its
// Jacobian G = dg/dq, whose columns are the spatial-Jacobian partial velocities
// velocity_field(S_j, P) already produced by the open-loop assembly -- so the loop
// closure reuses the open-loop machinery rather than introducing new geometry. The class
// offers three tiers of solve:
//
//   - position. assemble() Newton-solves g(q) = 0 for a consistent closed configuration.
//
//   - kinematics. solve_velocities() / solve_accelerations() distribute a driver joint's
//     rate / acceleration to the dependent joints via G q-dot = 0 and
//     G q-ddot = -G-dot q-dot.
//
//   - dynamics. joint_accelerations() / step() solve the bordered acceleration-level KKT
//     system for the joint accelerations and the constraint (Lagrange) forces,
//     integrating by RK4 with post-step GGL projection (position + velocity) for
//     energy-clean drift control.
//
// The dimension-agnostic numeric kernels (lstsq_solve / kkt_solve, in
// detail/ga_solver.hpp) are shared verbatim with the 3D layer
// (ga_pga3dp_ops_constraints.hpp).

#include "detail/ga_solver.hpp"        // hd::ga::lstsq_solve / kkt_solve (shared solvers)
#include "ga_pga2dp_ops_mechanics.hpp" // dynamic_system2dp (the reused spanning tree)
#include "ga_usr_utilities.hpp"        // hd::ga::rk4_step (shared RK4 integrator)
#include "ga_value_t.hpp"              // value_t

#include <cmath>     // std::abs
#include <cstddef>   // size_t
#include <mdspan>    // RK4 state views (step)
#include <stdexcept> // std::runtime_error
#include <string>
#include <vector>


namespace hd::ga::pga {

// Loop-closure constraint kinds. Currently point coincidence only (a planar pin / a
// spatial spherical-joint attachment); point-on-line, distance / rigid-link length and
// parallel/perpendicular are natural extensions (all expressible as vanishing
// wdg/rwdg/dot expressions).
enum class constraint2dp { coincidence };

// Loop-closure constraint descriptor: anchor point anchor_a (given in frame_a's
// coordinates) must satisfy `type` against anchor point anchor_b (in frame_b's
// coordinates). For coincidence the two anchors must coincide in world coordinates -- the
// loop-closure condition for a four-bar pin or a Stewart-Gough leg pinned to the
// platform.
//
// A pure-data aggregate (all public fields) with an fmt formatter, so it is also exposed
// in the Python wrapper (ga_py).
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

    // set the generalised RATE q-dot of a joint (the driver / independent joints);
    // refresh its relative velocity twist. The dependent rates come from
    // solve_velocities().
    void set_joint_rate(size_t joint_frame, value_t qdot)
    {
        tree_.joint[joint_frame].omega = qdot;
        tree_.apply_joint_state(joint_frame);
    }

    value_t joint_rate(size_t joint_frame) const
    {
        return tree_.joint_omega(joint_frame);
    }

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
            std::vector<value_t> b(g.size());
            for (size_t i = 0; i < g.size(); ++i)
                b[i] = -g[i];
            std::vector<value_t> const delta = hd::ga::lstsq_solve(G, b, dep.size());
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

    // --- velocity / acceleration distribution (kinematic closed loop)
    // -------------------

    // Distribute the driver joint rates to the dependent joint rates enforcing the
    // velocity-level loop closure  G q-dot = 0  (the time derivative of g(q) = 0). The
    // `driven` joints keep their currently-set rate (set_joint_rate); partitioning
    // G = [G_drv | G_dep] gives
    //
    //     G_dep q-dot_dep = -G_drv q-dot_drv,
    //
    // solved through the shared Jacobian solve. The solved dependent rates are written
    // into the tree (so point_velocity / twist_world become consistent) and returned in
    // dependent-joint order. Pre: assemble()d configuration (g ~ 0).
    std::vector<value_t> solve_velocities(std::vector<size_t> const& driven = {})
    {
        std::vector<size_t> const dep = dependent_joints(driven);
        size_t const m = 2 * loops_.size();

        // rhs = -G_drv q-dot_drv  (driver columns times the prescribed driver rates)
        std::vector<value_t> const Gdrv = constraint_jacobian(driven);
        std::vector<value_t> rhs(m, 0.0);
        for (size_t i = 0; i < m; ++i)
            for (size_t k = 0; k < driven.size(); ++k)
                rhs[i] -= Gdrv[i * driven.size() + k] * tree_.joint[driven[k]].omega;

        std::vector<value_t> const Gdep = constraint_jacobian(dep);
        std::vector<value_t> const qdot_dep = hd::ga::lstsq_solve(Gdep, rhs, dep.size());
        for (size_t k = 0; k < dep.size(); ++k) {
            tree_.joint[dep[k]].omega = qdot_dep[k];
            tree_.apply_joint_state(dep[k]); // rel velocity twist = q-dot * screw
        }
        return qdot_dep;
    }

    // Distribute the driver joint accelerations to the dependent joint accelerations
    // enforcing the acceleration-level loop closure  G q-ddot + G-dot q-dot = 0  (the
    // second derivative of g(q) = 0):
    //
    //     G_dep q-ddot_dep = -G_drv q-ddot_drv - G-dot q-dot.
    //
    // The velocity-product term G-dot q-dot is NOT formed explicitly: it is the relative
    // acceleration of the two anchor points evaluated with all DEPENDENT q-ddot = 0 and
    // the DRIVEN q-ddot prescribed -- i.e. the point-acceleration field already provided
    // by the kinematic layer (the same bias-pass trick assemble_mass_bias uses). So the
    // rhs is just -(a(P_a) - a(P_b)) read off with that accel-twist configuration.
    // `driven_accels` are the prescribed q-ddot of the `driven` joints (parallel; missing
    // entries default to 0, the constant-rate driver case). The solved dependent
    // accelerations are written into the relative accel twists (so point_acceleration is
    // consistent) and returned. Pre: solve_velocities() has run (the bias term needs the
    // velocities).
    std::vector<value_t>
    solve_accelerations(std::vector<size_t> const& driven = {},
                        std::vector<value_t> const& driven_accels = {})
    {
        std::vector<size_t> const dep = dependent_joints(driven);
        size_t const m = 2 * loops_.size();

        // accel-twist configuration for the bias read-off: dependent q-ddot = 0; driven
        // q-ddot prescribed (so G_drv q-ddot_drv folds into the point-acceleration
        // field).
        for (size_t k = 0; k < dep.size(); ++k)
            tree_.set_accel_twist(dep[k], twist2dp{0.0, 0.0, 0.0});
        for (size_t k = 0; k < driven.size(); ++k) {
            value_t const qdd = (k < driven_accels.size()) ? driven_accels[k] : 0.0;
            tree_.set_accel_twist(driven[k], qdd * tree_.joint[driven[k]].screw_b);
        }

        // rhs = -(a(P_a) - a(P_b)) per constraint, with the above accel-twist config
        std::vector<value_t> rhs(m, 0.0);
        for (size_t c = 0; c < loops_.size(); ++c) {
            auto const& lc = loops_[c];
            vec2dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec2dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            vec2dp const aa = tree_.point_acceleration(Pa, lc.frame_a);
            vec2dp const ab = tree_.point_acceleration(Pb, lc.frame_b);
            rhs[2 * c + 0] = -(aa.x - ab.x);
            rhs[2 * c + 1] = -(aa.y - ab.y);
        }

        std::vector<value_t> const Gdep = constraint_jacobian(dep);
        std::vector<value_t> const qddot_dep = hd::ga::lstsq_solve(Gdep, rhs, dep.size());
        for (size_t k = 0; k < dep.size(); ++k)
            tree_.set_accel_twist(dep[k], qddot_dep[k] * tree_.joint[dep[k]].screw_b);
        return qddot_dep;
    }

    // --- dynamic closed loop (constrained forward dynamics + integration)
    // ---------------

    // Constrained forward dynamics at the current state: solve the acceleration-level KKT
    // system for the joint accelerations q-ddot and the Lagrange multipliers lambda
    // (the constraint / leg forces):
    //
    //     | M   Gᵀ | | q-ddot |   |  tau    |
    //     |        | |        | = |         |
    //     | G   0  | | lambda |   | -G-dot q-dot |
    //
    // where M(q), tau(q,q-dot) (gravity + Coriolis/centripetal bias) come from the OPEN-
    // loop assembly (dynamic_system2dp::assemble_mass_bias), G is the constraint
    // Jacobian, and -G-dot q-dot is the velocity-product term (constraint_bias). Returns
    // q-ddot for all dof joints (in dof_joints() order); writes lambda into `lambda_out`
    // if non-null. The bordered system is solved by the shared LU.
    std::vector<value_t> joint_accelerations(std::vector<value_t>* lambda_out = nullptr)
    {
        return kkt_dynamics(tree_.dof_joints(), lambda_out);
    }

    // Advance the closed-loop system by dt. The coupled joint chain is integrated
    // together in its joint coordinates by RK4 (shared rk4_step), with the constrained
    // KKT solve (joint_accelerations) supplying q-ddot at each sub-step. After the step
    // the state is projected back onto the constraint manifold -- position by a min-norm
    // Newton step (assemble), velocity by the projection q-dot <- q-dot - G⁺(G q-dot) --
    // which removes numerical drift WITHOUT injecting energy (GGL-style stabilisation),
    // preserving the energy-conservation property the reduced-coordinate design is tuned
    // for.
    void step(value_t dt)
    {
        auto const rj = tree_.dof_joints();
        size_t const n = rj.size();
        if (n == 0) return;

        std::vector<value_t> u_mem(2 * n), uh_mem(2 * 2 * n), rhs_mem(2 * n);
        for (size_t k = 0; k < n; ++k) {
            u_mem[k] = tree_.joint[rj[k]].phi;
            u_mem[n + k] = tree_.joint[rj[k]].omega;
        }
        auto u = std::mdspan<value_t, std::dextents<size_t, 1>>(u_mem.data(), 2 * n);
        auto uh = std::mdspan<value_t, std::dextents<size_t, 2>>(uh_mem.data(), 2, 2 * n);
        auto const rhs =
            std::mdspan<value_t const, std::dextents<size_t, 1>>(rhs_mem.data(), 2 * n);

        auto apply_u = [&] {
            for (size_t k = 0; k < n; ++k) {
                tree_.joint[rj[k]].phi = u[k];
                tree_.joint[rj[k]].omega = u[n + k];
                tree_.apply_joint_state(rj[k]);
            }
        };
        for (size_t s = 1; s <= 4; ++s) {
            apply_u();
            auto const qdd = kkt_dynamics(rj, nullptr);
            for (size_t k = 0; k < n; ++k) {
                rhs_mem[k] = u[n + k];   // dq/dt = q-dot
                rhs_mem[n + k] = qdd[k]; // dq-dot/dt = q-ddot
            }
            rk4_step(u, uh, rhs, dt, s);
        }
        apply_u();

        // stabilisation: project (q, q-dot) back onto the constraint manifold
        assemble(/*driven*/ {}); // position: min-norm Newton -> g ~ 0
        project_velocities(rj);  // velocity:  q-dot <- q-dot - G⁺(G q-dot)
    }

  private:

    // world-coordinate, unitized (z = 1) position of an anchor point of a tree frame
    vec2dp anchor_world(size_t frame, vec2dp const& anchor)
    {
        vec2dp const P = move2dp(anchor, tree_.get_pos_trafo(frame, 0));
        return vec2dp{P.x / P.z, P.y / P.z, 1.0};
    }

    // Solve the bordered acceleration-level KKT system for q-ddot (and lambda) at the
    // current state -- see joint_accelerations() for the system. `rj` is the joint list
    // (all dof joints for dynamics). Reuses the open-loop M / tau (assemble_mass_bias),
    // the constraint Jacobian G, and the velocity-product term G-dot q-dot
    // (constraint_bias). Returns q-ddot (length n); writes lambda (length m) if
    // requested.
    std::vector<value_t> kkt_dynamics(std::vector<size_t> const& rj,
                                      std::vector<value_t>* lambda_out)
    {
        size_t const n = rj.size();
        size_t const m = 2 * loops_.size();

        // M (n*n) and tau (n) from the open-loop assembly; this also runs the bias pass
        // (zeroes the chain's relative accel twists), which is exactly the q-ddot = 0
        // state constraint_bias() needs below.
        auto const mb = tree_.assemble_mass_bias(rj);
        std::vector<value_t> const& M = mb.first;
        std::vector<value_t> const& tau = mb.second;

        std::vector<value_t> const G = constraint_jacobian(rj);
        std::vector<value_t> const gd = constraint_bias(); // G-dot q-dot (m)

        // bottom rhs of the saddle-point system is -G-dot q-dot
        std::vector<value_t> gbias(m);
        for (size_t c = 0; c < m; ++c)
            gbias[c] = -gd[c];

        // solve | M  Gᵀ | |q-ddot|   | tau   |  through the shared bordered LU
        //       | G  0  | |lambda| = | gbias |
        return hd::ga::kkt_solve(M, G, tau, gbias, n, m, lambda_out);
    }

    // Velocity-product term G-dot q-dot (length m): the relative acceleration of the two
    // anchor points evaluated at q-ddot = 0. With every dof joint's relative accel twist
    // zeroed, point_acceleration returns exactly the Coriolis/centripetal part, so
    // a(P_a) - a(P_b) = (G-dot q-dot)_c. Pre: the velocities (relative velocity twists)
    // are current. (Idempotent with the bias pass assemble_mass_bias just ran.)
    std::vector<value_t> constraint_bias()
    {
        for (size_t const j : tree_.dof_joints())
            tree_.set_accel_twist(j, twist2dp{0.0, 0.0, 0.0});
        std::vector<value_t> gd(2 * loops_.size(), 0.0);
        for (size_t c = 0; c < loops_.size(); ++c) {
            auto const& lc = loops_[c];
            vec2dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec2dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            vec2dp const aa = tree_.point_acceleration(Pa, lc.frame_a);
            vec2dp const ab = tree_.point_acceleration(Pb, lc.frame_b);
            gd[2 * c + 0] = aa.x - ab.x;
            gd[2 * c + 1] = aa.y - ab.y;
        }
        return gd;
    }

    // Project the joint rates onto the constraint tangent space (enforce G q-dot = 0)
    // without changing the feasible component: q-dot <- q-dot - G⁺(G q-dot), where the
    // min-norm pseudo-inverse correction G⁺(G q-dot) = Gᵀ(G Gᵀ)^-1 (G q-dot) comes from
    // the shared solve. Removes velocity drift after a step with no energy injection.
    void project_velocities(std::vector<size_t> const& rj)
    {
        size_t const n = rj.size();
        size_t const m = 2 * loops_.size();
        std::vector<value_t> const G = constraint_jacobian(rj);

        std::vector<value_t> qd(n), Gv(m, 0.0);
        for (size_t k = 0; k < n; ++k)
            qd[k] = tree_.joint[rj[k]].omega;
        for (size_t i = 0; i < m; ++i)
            for (size_t k = 0; k < n; ++k)
                Gv[i] += G[i * n + k] * qd[k];

        std::vector<value_t> const dqd = hd::ga::lstsq_solve(G, Gv, n); // G⁺(G q-dot)
        for (size_t k = 0; k < n; ++k) {
            tree_.joint[rj[k]].omega -= dqd[k];
            tree_.apply_joint_state(rj[k]);
        }
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
};

} // namespace hd::ga::pga
