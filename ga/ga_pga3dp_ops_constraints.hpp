#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Closed-loop (parallel mechanism) layer for PGA3DP -- the 3D counterpart of
// ga_pga2dp_ops_constraints.hpp. Same classes, methods and algorithms; only the
// representation changes (vec3dp points with weight w, twist3dp = BiVec3dp screw axes,
// move3dp, mvec3dp_e motors, and 3 residual components per coincidence instead of 2). The
// dimension-agnostic numeric kernels (lstsq_solve / kkt_solve, in detail/ga_solver.hpp)
// are shared verbatim with the 2D layer. closed_loop_system3dp HAS-A dynamic_system3dp
// (composition) as its spanning tree and closes kinematic loops by registering
// loop-closure constraints between existing tree frames.
//
// As in 2D, everything is built on the constraint residual g(q) and its Jacobian G =
// dg/dq (columns = the spatial-Jacobian partials velocity_field(S_j, P) of the open-loop
// assembly), with three tiers of solve:
//
//   - position. assemble() Newton-solves g(q) = 0 for a consistent closed configuration.
//
//   - kinematics. solve_velocities() / solve_accelerations() distribute a driver joint's
//     rate / acceleration via G q-dot = 0 and G q-ddot = -G-dot q-dot.
//
//   - dynamics. joint_accelerations() / step() solve the bordered acceleration-level KKT
//     system for the joint accelerations and the constraint (Lagrange) forces (M / tau
//     from the open-loop assemble_mass_bias), integrating by RK4 with post-step GGL
//     projection (position + velocity).

#include "detail/ga_solver.hpp"        // hd::ga::lstsq_solve / kkt_solve (shared solvers)
#include "ga_pga3dp_ops_mechanics.hpp" // dynamic_system3dp (the reused spanning tree)
#include "ga_usr_utilities.hpp"        // hd::ga::rk4_step (shared RK4 integrator)
#include "ga_value_t.hpp"              // value_t

#include <cmath>     // std::abs
#include <cstddef>   // size_t
#include <mdspan>    // RK4 state views (step)
#include <stdexcept> // std::runtime_error
#include <string>
#include <vector>


namespace hd::ga::pga {

// Loop-closure constraint kinds. Currently point coincidence only (a spatial spherical /
// ball joint); point-on-line, distance / rigid-link length and parallel/perpendicular are
// natural extensions.
enum class constraint3dp { coincidence };

// Loop-closure constraint descriptor: anchor point anchor_a (given in frame_a's
// coordinates) must satisfy `type` against anchor point anchor_b (in frame_b's
// coordinates). For coincidence the two anchors must coincide in world coordinates -- the
// loop-closure condition for a Stewart-Gough leg pinned to the platform.
//
// A pure-data aggregate (all public fields) with an fmt formatter, so it is also exposed
// in the Python wrapper (ga_py).
struct loop_constraint3dp {
    size_t frame_a;  // tree frame holding anchor a
    vec3dp anchor_a; // anchor point in frame_a coordinates (unitized, w = 1)
    size_t frame_b;  // tree frame holding anchor b
    vec3dp anchor_b; // anchor point in frame_b coordinates (unitized, w = 1)
    constraint3dp type{constraint3dp::coincidence};
};


// A closed-loop / parallel mechanism: an open-chain spanning tree (the reused
// dynamic_system3dp) plus loop-closure constraints between tree frames. Build the tree
// through the familiar forwarded API, then close loops with add_loop_constraint(...).
class closed_loop_system3dp {

    dynamic_system3dp tree_;                // the open-chain spanning tree (reused as-is)
    std::vector<loop_constraint3dp> loops_; // the extra closure edges

  public:

    closed_loop_system3dp() = default;

    // --- build the spanning tree (forward the familiar dynamic_system3dp builders)
    // ------

    void add_frame(static_frame3dp const& f,
                   size_t parent_idx = dynamic_system3dp::prev_frame)
    {
        tree_.add_frame(f, parent_idx);
    }

    void add_body(static_frame3dp const& f, body3dp const& b,
                  kin_state3dp const& k = kin_state3dp{},
                  size_t parent_idx = dynamic_system3dp::prev_frame)
    {
        tree_.add_body(f, b, k, parent_idx);
    }

    // revolute hinge about the body-fixed axis line through `pivot_b` (a finite point,
    // w = 1) along unit direction `axis_b` (w = 0)
    void add_revolute_body(static_frame3dp const& f, body3dp const& b,
                           vec3dp const& pivot_b, vec3dp const& axis_b,
                           value_t phi0 = 0.0, value_t omega0 = 0.0,
                           size_t parent_idx = dynamic_system3dp::prev_frame)
    {
        tree_.add_revolute_body(f, b, pivot_b, axis_b, phi0, omega0, parent_idx);
    }

    void add_prismatic_body(static_frame3dp const& f, body3dp const& b, vec3dp const& dir,
                            value_t s0 = 0.0, value_t v0 = 0.0,
                            size_t parent_idx = dynamic_system3dp::prev_frame)
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
    void add_loop_constraint(loop_constraint3dp const& c) { loops_.push_back(c); }

    // drive a joint to a generalised coordinate q (set q, refresh the kinematic state)
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
    loop_constraint3dp const& loop(size_t c) const { return loops_[c]; }

    // access the underlying spanning tree for diagnostics shared with the open-loop tier
    dynamic_system3dp& system() { return tree_; }
    dynamic_system3dp const& system() const { return tree_; }

    // --- constraint evaluation (position level)
    // -----------------------------------------

    // Residual g(q): for each coincidence constraint, the unitized world-coordinate
    // difference P_a - P_b (3 scalar equations per constraint in 3D). g == 0 <=> every
    // loop is closed. Length 3 * loop_count().
    std::vector<value_t> residual()
    {
        std::vector<value_t> g(3 * loops_.size());
        for (size_t c = 0; c < loops_.size(); ++c) {
            auto const& lc = loops_[c];
            vec3dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec3dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            g[3 * c + 0] = Pa.x - Pb.x;
            g[3 * c + 1] = Pa.y - Pb.y;
            g[3 * c + 2] = Pa.z - Pb.z;
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
    // iteration on the constraint Jacobian G (the relative partial velocities
    // velocity_field(S_j, P_a) - velocity_field(S_j, P_b), reused from the open-loop
    // spatial Jacobian). Returns the achieved residual norm ‖g‖. Throws if Newton does
    // not converge within max_iter.
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
                "closed_loop_system3dp::assemble: Newton did not converge within ") +
            std::to_string(max_iter) + std::string(" iterations (residual ") +
            std::to_string(residual_norm()) + std::string(")"));
    }

    // --- velocity / acceleration distribution (kinematic closed loop)
    // -------------------

    // Distribute the driver joint rates to the dependent joint rates enforcing G q-dot =
    // 0:
    //     G_dep q-dot_dep = -G_drv q-dot_drv.
    // The solved dependent rates are written into the tree and returned. Pre:
    // assemble()d.
    std::vector<value_t> solve_velocities(std::vector<size_t> const& driven = {})
    {
        std::vector<size_t> const dep = dependent_joints(driven);
        size_t const m = 3 * loops_.size();

        std::vector<value_t> const Gdrv = constraint_jacobian(driven);
        std::vector<value_t> rhs(m, 0.0);
        for (size_t i = 0; i < m; ++i)
            for (size_t k = 0; k < driven.size(); ++k)
                rhs[i] -= Gdrv[i * driven.size() + k] * tree_.joint[driven[k]].omega;

        std::vector<value_t> const Gdep = constraint_jacobian(dep);
        std::vector<value_t> const qdot_dep = hd::ga::lstsq_solve(Gdep, rhs, dep.size());
        for (size_t k = 0; k < dep.size(); ++k) {
            tree_.joint[dep[k]].omega = qdot_dep[k];
            tree_.apply_joint_state(dep[k]);
        }
        return qdot_dep;
    }

    // Distribute the driver joint accelerations to the dependent ones enforcing
    //     G_dep q-ddot_dep = -G_drv q-ddot_drv - G-dot q-dot,
    // with the velocity-product term read off as the relative anchor acceleration at
    // dependent q-ddot = 0 (driven q-ddot prescribed) via point_acceleration. Pre:
    // solve_velocities() has run.
    std::vector<value_t>
    solve_accelerations(std::vector<size_t> const& driven = {},
                        std::vector<value_t> const& driven_accels = {})
    {
        std::vector<size_t> const dep = dependent_joints(driven);
        size_t const m = 3 * loops_.size();

        for (size_t k = 0; k < dep.size(); ++k)
            tree_.set_accel_twist(dep[k], twist3dp{});
        for (size_t k = 0; k < driven.size(); ++k) {
            value_t const qdd = (k < driven_accels.size()) ? driven_accels[k] : 0.0;
            tree_.set_accel_twist(driven[k], qdd * tree_.joint[driven[k]].screw_b);
        }

        std::vector<value_t> rhs(m, 0.0);
        for (size_t c = 0; c < loops_.size(); ++c) {
            auto const& lc = loops_[c];
            vec3dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec3dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            vec3dp const aa = tree_.point_acceleration(Pa, lc.frame_a);
            vec3dp const ab = tree_.point_acceleration(Pb, lc.frame_b);
            rhs[3 * c + 0] = -(aa.x - ab.x);
            rhs[3 * c + 1] = -(aa.y - ab.y);
            rhs[3 * c + 2] = -(aa.z - ab.z);
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
    // system for the joint accelerations q-ddot and the Lagrange multipliers lambda (the
    // constraint / leg forces). M(q), tau(q,q-dot) come from the open-loop assembly
    // (assemble_mass_bias), G from constraint_jacobian, and -G-dot q-dot from
    // constraint_bias. Returns q-ddot (dof_joints() order); writes lambda if requested.
    std::vector<value_t> joint_accelerations(std::vector<value_t>* lambda_out = nullptr)
    {
        return kkt_dynamics(tree_.dof_joints(), lambda_out);
    }

    // Advance the closed-loop system by dt: RK4 (shared rk4_step) with the constrained
    // KKT solve at each sub-step, then post-step GGL projection (position by min-norm
    // Newton, velocity by q-dot <- q-dot - G⁺(G q-dot)) for energy-clean drift control.
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

        assemble(/*driven*/ {}); // position: min-norm Newton -> g ~ 0
        project_velocities(rj);  // velocity:  q-dot <- q-dot - G⁺(G q-dot)
    }

  private:

    // world-coordinate, unitized (w = 1) position of an anchor point of a tree frame
    vec3dp anchor_world(size_t frame, vec3dp const& anchor)
    {
        return unitize(move3dp(anchor, tree_.get_pos_trafo(frame, 0)));
    }

    // Solve the bordered acceleration-level KKT system for q-ddot (and lambda) at the
    // current state -- see joint_accelerations(). `rj` is the joint list (all dof joints
    // for dynamics). Reuses the open-loop M / tau (assemble_mass_bias), G
    // (constraint_jacobian) and -G-dot q-dot (constraint_bias).
    std::vector<value_t> kkt_dynamics(std::vector<size_t> const& rj,
                                      std::vector<value_t>* lambda_out)
    {
        size_t const n = rj.size();
        size_t const m = 3 * loops_.size();

        // M (n*n) and tau (n) from the open-loop assembly; this also runs the bias pass
        // (zeroes the chain's relative accel twists), the q-ddot = 0 state
        // constraint_bias() needs below.
        auto const mb = tree_.assemble_mass_bias(rj);
        std::vector<value_t> const& M = mb.first;
        std::vector<value_t> const& tau = mb.second;

        std::vector<value_t> const G = constraint_jacobian(rj);
        std::vector<value_t> const gd = constraint_bias(); // G-dot q-dot (m)

        std::vector<value_t> gbias(m);
        for (size_t c = 0; c < m; ++c)
            gbias[c] = -gd[c];

        return hd::ga::kkt_solve(M, G, tau, gbias, n, m, lambda_out);
    }

    // Velocity-product term G-dot q-dot (length m): the relative acceleration of the two
    // anchor points at q-ddot = 0. With every dof joint's relative accel twist zeroed,
    // point_acceleration returns exactly the Coriolis/centripetal part.
    std::vector<value_t> constraint_bias()
    {
        for (size_t const j : tree_.dof_joints())
            tree_.set_accel_twist(j, twist3dp{});
        std::vector<value_t> gd(3 * loops_.size(), 0.0);
        for (size_t c = 0; c < loops_.size(); ++c) {
            auto const& lc = loops_[c];
            vec3dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec3dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            vec3dp const aa = tree_.point_acceleration(Pa, lc.frame_a);
            vec3dp const ab = tree_.point_acceleration(Pb, lc.frame_b);
            gd[3 * c + 0] = aa.x - ab.x;
            gd[3 * c + 1] = aa.y - ab.y;
            gd[3 * c + 2] = aa.z - ab.z;
        }
        return gd;
    }

    // Project the joint rates onto the constraint tangent space (enforce G q-dot = 0):
    // q-dot <- q-dot - G⁺(G q-dot). Removes velocity drift with no energy injection.
    void project_velocities(std::vector<size_t> const& rj)
    {
        size_t const n = rj.size();
        size_t const m = 3 * loops_.size();
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

    // the dof joints NOT in `driven` -- the coordinates assemble() / solvers solve for
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

    // Constraint Jacobian G (m x ndep, row-major), m = 3 * loop_count(), columns indexed
    // by the joints `dep`. Column k for constraint c is the relative partial velocity of
    // the two anchors w.r.t. joint dep[k]:
    //
    //     dP_a/dq_j - dP_b/dq_j = velocity_field(S_j, P_a) - velocity_field(S_j, P_b)
    //
    // (each term present only when joint j is an ancestor of the anchor's frame), where
    // S_j = move3dp(screw_b, M_{j->world}) is the world joint screw -- the SAME column
    // that builds the open-loop mass matrix and bias forces (assemble_mass_bias).
    std::vector<value_t> constraint_jacobian(std::vector<size_t> const& dep)
    {
        size_t const m = 3 * loops_.size();
        size_t const ndep = dep.size();
        std::vector<value_t> G(m * ndep, 0.0);

        std::vector<twist3dp> S(ndep);
        for (size_t k = 0; k < ndep; ++k)
            S[k] = move3dp(tree_.joint[dep[k]].screw_b, tree_.get_pos_trafo(dep[k], 0));

        for (size_t c = 0; c < loops_.size(); ++c) {
            auto const& lc = loops_[c];
            vec3dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec3dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            for (size_t k = 0; k < ndep; ++k) {
                size_t const j = dep[k];
                vec3dp dPa{};
                vec3dp dPb{};
                if (tree_.is_ancestor(j, lc.frame_a))
                    dPa = dynamic_system3dp::velocity_field(S[k], Pa);
                if (tree_.is_ancestor(j, lc.frame_b))
                    dPb = dynamic_system3dp::velocity_field(S[k], Pb);
                G[(3 * c + 0) * ndep + k] = dPa.x - dPb.x;
                G[(3 * c + 1) * ndep + k] = dPa.y - dPb.y;
                G[(3 * c + 2) * ndep + k] = dPa.z - dPb.z;
            }
        }
        return G;
    }
};

} // namespace hd::ga::pga
