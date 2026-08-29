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
//   - impacts. impact() applies the impulsive velocity jump a newly active (or violated)
//     constraint demands -- the same bordered system at velocity level, mass-weighted,
//     with a coefficient of restitution -- and reports the impulses.
//
// A loop closure is a JOINT between two tree frames (Featherstone, Rigid Body Dynamics
// Algorithms, ch. 8): its kind is the subspace of relative motion it removes, and the
// rows it contributes to g / G are that constraint-force subspace. The kinds provided
// here are the planar members of that table (see constraint2dp), and a constraint can be
// switched off and on at run time (loop_constraint2dp::active) -- a contact that opens
// and closes -- with the row count m following.
//
// The tree's generalised coordinates are dynamic_system2dp::coord -- one per screw of
// every dynamic joint, so a floating body (a free joint with mass) takes part in the
// loop closure like a revolute one; the layer indexes q, q-dot, q-ddot and the Jacobian
// columns by those coordinates.
//
// The dimension-agnostic numeric kernels (lstsq_solve / kkt_solve, in
// detail/ga_solver.hpp) are shared verbatim with the 3D layer
// (ga_pga3dp_ops_constraints.hpp).

#include "detail/ga_solver.hpp"        // hd::ga::lstsq_solve / kkt_solve (shared solvers)
#include "ga_pga2dp_ops_mechanics.hpp" // dynamic_system2dp (the reused spanning tree)
#include "ga_usr_utilities.hpp"        // hd::ga::rk4_step (shared RK4 integrator)
#include "ga_value_t.hpp"              // value_t

#include <cmath>     // std::abs, std::sqrt
#include <cstddef>   // size_t
#include <mdspan>    // RK4 state views (step)
#include <stdexcept> // std::runtime_error
#include <string>
#include <vector>


namespace hd::ga::pga {

// Loop-closure constraint kinds -- each a subspace of relative motion removed between
// the two anchors, i.e. the rows of g and G it contributes (Featherstone's constraint-
// force subspace T of the loop joint):
//
//   coincidence : the two anchor POINTS coincide (2 rows) -- a planar pin joint. Rotation
//                 about the pin stays free.
//
//   distance    : the two anchors keep the DISTANCE `length` (1 row) -- a rigid rod /
//                 a Stewart-Gough leg of fixed length; the multiplier is the rod force.
//
//   frame       : a WELD -- the two bodies are rigidly attached to each other (3 rows).
//                 Position: anchor_a and anchor_b coincide in the world, as for
//                 coincidence. Orientation: frame_a and frame_b have the same world
//                 heading (the anchors only locate the attachment point; they do not
//                 rotate the frames, so a weld between frames built with different
//                 headings must be closed by assemble() or built to match). Every
//                 relative dof of the loop is removed -- the zero-dof loop joint: a tool
//                 fixed in a hand, a foot standing flat, a rigid coupler.
//
// Multiplier convention (joint_accelerations): the KKT system is M q-ddot + G^T lambda
// = tau, so the constraint force on the system is -G^T lambda. Per kind, on anchor a:
// coincidence -- the force -lambda (x, y); distance -- the force -lambda n along the rod
// axis n from anchor b to anchor a, so lambda > 0 is TENSION; frame -- the force
// -lambda_F at the anchor and the couple -lambda_M. Anchor b receives the opposite.
//
// Point-on-line and parallel/perpendicular are natural further members (vanishing
// wdg/rwdg/dot expressions), added when a mechanism needs them.
enum class constraint2dp { coincidence, distance, frame };

// Loop-closure constraint descriptor: anchor point anchor_a (given in frame_a's
// coordinates) must satisfy `type` against anchor point anchor_b (in frame_b's
// coordinates). For coincidence the two anchors must coincide in world coordinates -- the
// loop-closure condition for a four-bar pin or a Stewart-Gough leg pinned to the
// platform; for distance they keep `length` apart; for frame the two frames coincide at
// the anchors.
//
// `active` switches the constraint on and off at run time (contact state: a foot that
// lifts, a hand that lets go). An inactive constraint contributes no rows -- the system
// behaves exactly as if it had not been added -- while keeping its index, so a repeating
// gait can toggle the same constraints without renumbering.
//
// A pure-data aggregate (all public fields) with an fmt formatter, so it is also exposed
// in the Python wrapper (ga_py). Adding a field here widens the bound constructor
// (CLAUDE.md, "Adding a FIELD to an already-bound pure-data struct").
struct loop_constraint2dp {
    size_t frame_a;  // tree frame holding anchor a
    vec2dp anchor_a; // anchor point in frame_a coordinates (unitized, z = 1)
    size_t frame_b;  // tree frame holding anchor b
    vec2dp anchor_b; // anchor point in frame_b coordinates (unitized, z = 1)
    constraint2dp type{constraint2dp::coincidence};
    value_t length{0.0}; // the rod length (distance kind only)
    bool active{true};   // contributes rows to g / G when true
};


// A closed-loop / parallel mechanism: an open-chain spanning tree (the reused
// dynamic_system2dp) plus loop-closure constraints between tree frames. Build the tree
// through the familiar forwarded API, then close loops with add_loop_constraint(...).
class closed_loop_system2dp {

    using coord = dynamic_system2dp::coord;

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

    // register a loop-closure constraint between two existing tree frames; returns its
    // index (stable for the life of the system, also while inactive)
    size_t add_loop_constraint(loop_constraint2dp const& c)
    {
        loops_.push_back(c);
        return loops_.size() - 1;
    }

    // switch constraint c on / off (contact state). Off: it contributes no rows and the
    // system behaves as if it were absent; on: it takes part again. The multipliers
    // returned by joint_accelerations() are those of the ACTIVE constraints, in index
    // order, each contributing rows_of(loop(c)) entries.
    void set_loop_active(size_t c, bool on) { loops_.at(c).active = on; }
    bool loop_active(size_t c) const { return loops_.at(c).active; }

    // move a constraint's anchors (a contact that closes at a NEW point each time -- a
    // foot touching down further along the ground). The kind and the frames stay; the
    // typical sequence is set_loop_anchors, then activate_loop_with_impact.
    void set_loop_anchors(size_t c, vec2dp const& anchor_a, vec2dp const& anchor_b)
    {
        loops_.at(c).anchor_a = anchor_a;
        loops_.at(c).anchor_b = anchor_b;
    }

    // drive a 1-dof joint to a generalised coordinate q (set q, refresh the kinematic
    // state). For the driver / independent joints of a closed loop -- the dependent ones
    // are produced by assemble(). A floating body is posed through
    // system().set_joint_motor().
    void set_joint(size_t joint_frame, value_t q) { tree_.set_joint(joint_frame, q); }

    value_t joint_phi(size_t joint_frame) const { return tree_.joint_phi(joint_frame); }

    // set the generalised RATE q-dot of a 1-dof joint (the driver / independent joints);
    // refresh its relative velocity twist. The dependent rates come from
    // solve_velocities().
    void set_joint_rate(size_t joint_frame, value_t qdot)
    {
        tree_.set_joint_rate(joint_frame, qdot);
    }

    value_t joint_rate(size_t joint_frame) const
    {
        return tree_.joint_omega(joint_frame);
    }

    size_t loop_count() const { return loops_.size(); } // registered, active or not

    // First row index of constraint c in g / G / lambda (active constraints in
    // registration order, rows_of() each) -- the row mapping a caller needs to
    // address one constraint's block, e.g. a commanded reaction lambda_d. For an
    // inactive constraint (no rows) it returns constraint_rows().
    size_t constraint_row_offset(size_t c) const
    {
        if (!loops_.at(c).active) return constraint_rows();
        size_t r = 0;
        for (size_t k = 0; k < c; ++k)
            if (loops_[k].active) r += rows_of(loops_[k]);
        return r;
    }
    size_t active_loop_count() const
    {
        size_t n = 0;
        for (auto const& lc : loops_)
            if (lc.active) ++n;
        return n;
    }
    loop_constraint2dp const& loop(size_t c) const { return loops_.at(c); }

    // rows a constraint of this kind contributes to g / G (its constraint subspace)
    static size_t rows_of(loop_constraint2dp const& lc)
    {
        switch (lc.type) {
            case constraint2dp::coincidence:
                return 2;
            case constraint2dp::distance:
                return 1;
            case constraint2dp::frame:
                return 3;
        }
        return 0;
    }

    // m: the total number of active constraint rows (the size of g and of lambda)
    size_t constraint_rows() const
    {
        size_t m = 0;
        for (auto const& lc : loops_)
            if (lc.active) m += rows_of(lc);
        return m;
    }

    // Rank of the active constraint Jacobian over the dof coordinates -- the number of
    // independent constraints at THIS configuration. Featherstone 8.10: it varies as the
    // mechanism moves, and mobility = n - rank. Below constraint_rows() the constraints
    // are dependent (a knee lock, an over-constrained loop) and the multipliers are only
    // determined up to that dependence; joint_accelerations() then returns their
    // minimum-norm values.
    size_t constraint_rank()
    {
        auto const rc = tree_.dof_coords();
        size_t const m = constraint_rows();
        if (m == 0 || rc.empty()) return 0;
        return hd::ga::matrix_rank(constraint_jacobian(rc), m, rc.size());
    }

    // The active constraint Jacobian G (m x n, row-major; m = constraint_rows(),
    // columns over dof_coords() in their order) -- the G of
    //
    //     M(q) q-ddot = RHS(q, q-dot) + tau - G^T lambda,        G q-ddot = bias
    //
    // (the lambda convention of joint_accelerations(): a loaded foot has lambda_y < 0).
    // Together with the tree's mass_bias() this is the whole constrained equation of
    // motion read off for an external law (computed torque commanding a desired
    // reaction G^T lambda_d). The row blocks follow the active constraints in
    // registration order, rows_of() each.
    std::vector<value_t> constraint_jacobian()
    {
        return constraint_jacobian(tree_.dof_coords());
    }

    // access the underlying spanning tree for diagnostics shared with the open-loop tier
    // (get_pos_trafo, total_energy, joint_phi, ...) and for further configuration
    dynamic_system2dp& system() { return tree_; }
    dynamic_system2dp const& system() const { return tree_; }

    // --- constraint evaluation (position level)
    // -----------------------------------------

    // Residual g(q), one block per ACTIVE constraint in index order: coincidence -- the
    // unitized world-coordinate difference P_a - P_b (2 rows); distance -- |P_a - P_b| -
    // length (1 row); frame -- P_a - P_b (2 rows) and the relative rotation angle of the
    // two frames, theta_a - theta_b (1 row). g == 0 <=> every active loop is closed.
    // Length constraint_rows().
    std::vector<value_t> residual()
    {
        std::vector<value_t> g(constraint_rows());
        size_t r = 0;
        for (auto const& lc : loops_) {
            if (!lc.active) continue;
            vec2dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec2dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            switch (lc.type) {
                case constraint2dp::coincidence:
                    g[r++] = Pa.x - Pb.x;
                    g[r++] = Pa.y - Pb.y;
                    break;
                case constraint2dp::distance:
                    g[r++] = separation(Pa, Pb) - lc.length;
                    break;
                case constraint2dp::frame:
                    g[r++] = Pa.x - Pb.x;
                    g[r++] = Pa.y - Pb.y;
                    g[r++] = relative_rotation(lc.frame_a, lc.frame_b);
                    break;
            }
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
    // solve the remaining (dependent) coordinates so that g(q) = 0, by Newton iteration
    // on the constraint Jacobian G (§3.4: G's columns are the relative partial
    // velocities velocity_field(S_j, P_a) - velocity_field(S_j, P_b), reused from the
    // open-loop spatial Jacobian). Each Newton increment is applied as a coordinate step
    // -- q += delta for a 1-dof joint, the retraction M <- M (x) rexp(1/2 delta . S) for
    // a floating body. Returns the achieved residual norm ‖g‖. Throws if Newton does not
    // converge within max_iter.
    //
    // For a 1-DOF four-bar (3 revolute joints, 2 coincidence equations) drive one joint;
    // the 2 dependent joints are then solved by the square Newton step.
    value_t assemble(std::vector<size_t> const& driven = {}, value_t tol = value_t(1e-12),
                     size_t max_iter = 50, value_t stall_tol = value_t(1e-3))
    {
        std::vector<coord> const dep = dependent_coords(driven);

        for (size_t it = 0; it <= max_iter; ++it) {
            std::vector<value_t> const g = residual();
            value_t gnorm = 0.0;
            for (value_t const gi : g)
                gnorm = std::max(gnorm, std::abs(gi));
            if (gnorm < tol) return gnorm;
            if (it == max_iter) break;

            // Newton step: solve G_dep * delta = -g for the dependent increments
            std::vector<value_t> const G = constraint_jacobian(dep);
            std::vector<value_t> b(g.size());
            for (size_t i = 0; i < g.size(); ++i)
                b[i] = -g[i];
            std::vector<value_t> const delta = hd::ga::lstsq_solve(G, b, dep.size());
            tree_.increment_coords(dep, delta);
        }
        // A rank-deficient closure (a straightened knee, an over-constrained loop)
        // leaves a component of g outside the constraint Jacobian's row space, and
        // the min-norm Newton STALLS at that least-squares floor rather than
        // converging. A sub-millimetre floor is the singular configuration's
        // answer, not a failure (measured: a walking biped drifting near knee lock
        // stalled at 3e-4 and the throw took the whole app down) -- return the
        // achieved residual. A large floor still throws.
        {
            std::vector<value_t> const g = residual();
            value_t gnorm = 0.0;
            for (value_t const gi : g)
                gnorm = std::max(gnorm, std::abs(gi));
            if (gnorm < stall_tol) return gnorm;
        }
        throw std::runtime_error(
            std::string(
                "closed_loop_system2dp::assemble: Newton did not converge within ") +
            std::to_string(max_iter) + std::string(" iterations (residual ") +
            std::to_string(residual_norm()) + std::string(")"));
    }

    // --- velocity / acceleration distribution (kinematic closed loop)
    // -------------------

    // Distribute the driver joint rates to the dependent coordinate rates enforcing the
    // velocity-level loop closure  G q-dot = 0  (the time derivative of g(q) = 0). The
    // `driven` joints keep their currently-set rate (set_joint_rate); partitioning
    // G = [G_drv | G_dep] gives
    //
    //     G_dep q-dot_dep = -G_drv q-dot_drv,
    //
    // solved through the shared Jacobian solve. The solved dependent rates are written
    // into the tree (so point_velocity / twist_world become consistent) and returned in
    // dependent-coordinate order. Pre: assemble()d configuration (g ~ 0).
    std::vector<value_t> solve_velocities(std::vector<size_t> const& driven = {})
    {
        std::vector<coord> const drv = tree_.coords_of(driven);
        std::vector<coord> const dep = dependent_coords(driven);
        size_t const m = constraint_rows();

        // rhs = -G_drv q-dot_drv  (driver columns times the prescribed driver rates)
        std::vector<value_t> const Gdrv = constraint_jacobian(drv);
        std::vector<value_t> rhs(m, 0.0);
        for (size_t i = 0; i < m; ++i)
            for (size_t k = 0; k < drv.size(); ++k)
                rhs[i] -= Gdrv[i * drv.size() + k] * tree_.coord_rate(drv[k]);

        std::vector<value_t> const Gdep = constraint_jacobian(dep);
        std::vector<value_t> const qdot_dep = hd::ga::lstsq_solve(Gdep, rhs, dep.size());
        write_rates(dep, qdot_dep);
        return qdot_dep;
    }

    // Distribute the driver joint accelerations to the dependent coordinate
    // accelerations enforcing the acceleration-level loop closure  G q-ddot + G-dot q-dot
    // = 0  (the second derivative of g(q) = 0):
    //
    //     G_dep q-ddot_dep = -G_drv q-ddot_drv - G-dot q-dot.
    //
    // The velocity-product term G-dot q-dot is NOT formed explicitly: it is the second
    // derivative of the residual evaluated with all DEPENDENT q-ddot = 0 and the DRIVEN
    // q-ddot prescribed -- i.e. the point-acceleration field already provided by the
    // kinematic layer (the same bias-pass trick assemble_mass_bias uses). So the rhs is
    // just -(the relative anchor accelerations, per kind) read off with that accel-twist
    // configuration. `driven_accels` are the prescribed q-ddot of the driven joints'
    // coordinates (parallel to coords_of(driven); missing entries default to 0, the
    // constant-rate driver case). The solved dependent accelerations are written into the
    // relative accel twists (so point_acceleration is consistent) and returned. Pre:
    // solve_velocities() has run (the bias term needs the velocities).
    std::vector<value_t>
    solve_accelerations(std::vector<size_t> const& driven = {},
                        std::vector<value_t> const& driven_accels = {})
    {
        std::vector<coord> const drv = tree_.coords_of(driven);
        std::vector<coord> const dep = dependent_coords(driven);

        // accel-twist configuration for the bias read-off: dependent q-ddot = 0; driven
        // q-ddot prescribed (so G_drv q-ddot_drv folds into the point-acceleration
        // field).
        for (auto const& c : dep)
            tree_.set_accel_twist(c.frame, twist2dp{0.0, 0.0, 0.0});
        for (auto const& c : drv)
            tree_.set_accel_twist(c.frame, twist2dp{0.0, 0.0, 0.0});
        for (size_t k = 0; k < drv.size(); ++k) {
            value_t const qdd = (k < driven_accels.size()) ? driven_accels[k] : 0.0;
            tree_.set_accel_twist(drv[k].frame, tree_.relative_accel_twist(drv[k].frame) +
                                                    qdd * tree_.screw_of(drv[k]));
        }

        // rhs = -(second derivative of g at that accel-twist configuration)
        std::vector<value_t> rhs = relative_accelerations();
        for (value_t& v : rhs)
            v = -v;

        std::vector<value_t> const Gdep = constraint_jacobian(dep);
        std::vector<value_t> const qddot_dep = hd::ga::lstsq_solve(Gdep, rhs, dep.size());
        for (size_t k = 0; k < dep.size(); ++k)
            tree_.set_accel_twist(dep[k].frame,
                                  tree_.relative_accel_twist(dep[k].frame) +
                                      qddot_dep[k] * tree_.screw_of(dep[k]));
        return qddot_dep;
    }

    // --- dynamic closed loop (constrained forward dynamics + integration)
    // ---------------

    // Constrained forward dynamics at the current state: solve the acceleration-level KKT
    // system for the coordinate accelerations q-ddot and the Lagrange multipliers lambda
    // (the constraint / leg forces):
    //
    //     | M   Gᵀ | | q-ddot |   |  tau    |
    //     |        | |        | = |         |
    //     | G   0  | | lambda |   | -G-dot q-dot |
    //
    // where M(q), tau(q,q-dot) (gravity + Coriolis/centripetal bias) come from the OPEN-
    // loop assembly (dynamic_system2dp::assemble_mass_bias), G is the constraint
    // Jacobian, and -G-dot q-dot is the velocity-product term (constraint_bias). Returns
    // q-ddot for all dof coordinates (in dof_coords() order); writes lambda into
    // `lambda_out` if non-null -- one block per ACTIVE constraint, rows_of() entries
    // each, with the sign convention stated at constraint2dp (the force on anchor a is
    // -lambda). The bordered system is solved by the shared LU.
    std::vector<value_t> joint_accelerations(std::vector<value_t>* lambda_out = nullptr)
    {
        return kkt_dynamics(tree_.dof_coords(), lambda_out);
    }

    // The CONSTRAINED analogue of dynamic_system2dp::sync_accelerations(): solve the
    // KKT accelerations and write them into the per-frame relative acceleration
    // twists, so point_acceleration / centre_of_mass_acceleration then report the
    // actual constrained dynamics rather than the velocity-product bias left after
    // step(). The base class's version cannot serve here: it runs the FREE forward
    // dynamics, which for a mechanism standing on two closed contacts is a different
    // problem with a different answer. Call it before reading any acceleration field
    // -- for visualisation, and for a controller that measures a task's response
    // rather than differentiating a Jacobian for it.
    std::vector<value_t> sync_accelerations()
    {
        auto const rc = tree_.dof_coords();
        auto const qdd = joint_accelerations();
        for (auto const& c : rc)
            tree_.set_accel_twist(c.frame, twist2dp{});
        for (size_t c = 0; c < rc.size(); ++c)
            tree_.set_accel_twist(rc[c].frame, tree_.relative_accel_twist(rc[c].frame) +
                                                   qdd[c] * tree_.screw_of(rc[c]));
        return qdd;
    }

    // Advance the closed-loop system by dt. The coupled coordinate state (read_state /
    // write_state / state_rates / commit_state of the tree -- 1-dof joints and floating
    // bodies alike) is integrated by RK4 (shared rk4_step), with the constrained KKT
    // solve (joint_accelerations) supplying q-ddot at each sub-step. After the step the
    // state is projected back onto the constraint manifold -- position by a min-norm
    // Newton step (assemble), velocity by the projection q-dot <- q-dot - G⁺(G q-dot) --
    // which removes numerical drift WITHOUT injecting energy (GGL-style stabilisation),
    // preserving the energy-conservation property the reduced-coordinate design is tuned
    // for. That projection is for ROUND-OFF only: a constraint switched on while the
    // velocities disagree with it is an impact, and impact() (the mass-weighted jump)
    // must be applied before the next step -- see activate_loop_with_impact().
    void step(value_t dt)
    {
        auto const rc = tree_.dof_coords();
        size_t const n = rc.size();
        if (n == 0) return;

        std::vector<value_t> u_mem, uh_mem(2 * 2 * n), rhs_mem(2 * n);
        tree_.read_state(rc, u_mem);
        auto u = std::mdspan<value_t, std::dextents<size_t, 1>>(u_mem.data(), 2 * n);
        auto uh = std::mdspan<value_t, std::dextents<size_t, 2>>(uh_mem.data(), 2, 2 * n);
        auto const rhs =
            std::mdspan<value_t const, std::dextents<size_t, 1>>(rhs_mem.data(), 2 * n);

        // thread the RK4 stage time into the tree's clock so a time-varying applied
        // wrench is sampled at the stage time (as dynamic_system's step does), and
        // advance it by dt at the end
        value_t const t0 = tree_.time();
        for (size_t s = 1; s <= 4; ++s) {
            tree_.set_time(rk4_get_time(t0, dt, s - 1));
            tree_.write_state(rc, u_mem);
            auto const qdd = kkt_dynamics(rc, nullptr);
            tree_.state_rates(rc, u_mem, qdd, rhs_mem);
            rk4_step(u, uh, rhs, dt, s);
        }
        tree_.set_time(t0 + dt);
        tree_.commit_state(rc, u_mem);

        // stabilisation: project (q, q-dot) back onto the constraint manifold
        if (constraint_rows() > 0) {
            assemble(/*driven*/ {}); // position: min-norm Newton -> g ~ 0
            project_velocities(rc);  // velocity:  q-dot <- q-dot - G⁺(G q-dot)
        }
    }


    // The IMPACT MAP: the impulsive velocity jump the active constraints demand when the
    // current rates violate G q-dot = 0 (a constraint that just appeared -- a foot
    // touching down -- or velocities set by hand). Featherstone sec. 11.7: over the
    // instant the configuration is frozen and momentum balance reads
    //
    //     M (q-dot+ - q-dot-) = -Gᵀ Λ,        G q-dot+ = -e G q-dot-
    //
    // with the impulse Λ per constraint row and the coefficient of restitution e (0:
    // plastic, the default -- the rows then stick; 1: elastic). Written for the jump
    // dq = q-dot+ - q-dot- this is the bordered system of joint_accelerations with a
    // zero top rhs and -(1 + e) G q-dot- below, so the shared rank-aware kkt_solve
    // does it (min-norm at a dependent row set); a consistent state gives dq = 0
    // EXACTLY, so touching down at zero closing velocity leaves the state untouched.
    // Λ follows the multiplier convention of joint_accelerations: the impulse on the
    // system is -Gᵀ Λ (a loaded foot has Λ_y < 0 like its λ_y). The kinetic energy
    // change is -1/2 dqᵀ M dq <= 0. The rates are written back; the position is
    // unchanged (an impact is instantaneous). Returns Λ, one block per active row set.
    //
    // This is distinct from the drift stabilisation of step(): that projection is the
    // Euclidean min-norm correction in coordinate space and is right for round-off, not
    // for physics -- a light foot and a heavy torso would be corrected alike.
    std::vector<value_t> impact(value_t e = 0.0)
    {
        auto const rc = tree_.dof_coords();
        size_t const n = rc.size();
        size_t const m = constraint_rows();
        if (n == 0 || m == 0) return {};

        auto const mb = tree_.assemble_mass_bias(rc); // M (n*n); the bias is unused
        std::vector<value_t> const& M = mb.first;
        std::vector<value_t> const G = constraint_jacobian(rc);

        std::vector<value_t> qd(n), Gv(m, 0.0);
        for (size_t k = 0; k < n; ++k)
            qd[k] = tree_.coord_rate(rc[k]);
        for (size_t i = 0; i < m; ++i)
            for (size_t k = 0; k < n; ++k)
                Gv[i] += G[i * n + k] * qd[k];

        std::vector<value_t> const zero(n, 0.0);
        std::vector<value_t> gb(m);
        for (size_t i = 0; i < m; ++i)
            gb[i] = -(1.0 + e) * Gv[i];

        std::vector<value_t> Lambda;
        std::vector<value_t> const dq = hd::ga::kkt_solve(M, G, zero, gb, n, m, &Lambda);
        for (size_t k = 0; k < n; ++k)
            qd[k] += dq[k];
        write_rates(rc, qd);
        return Lambda;
    }

    // touchdown in one call: switch constraint c on and resolve the velocities it now
    // demands (impact with restitution e). The anchors must already coincide (g ~ 0 for
    // that constraint) -- an impact fixes velocities, not positions.
    std::vector<value_t> activate_loop_with_impact(size_t c, value_t e = 0.0)
    {
        set_loop_active(c, true);
        return impact(e);
    }

  private:

    // world-coordinate, unitized (z = 1) position of an anchor point of a tree frame
    vec2dp anchor_world(size_t frame, vec2dp const& anchor)
    {
        vec2dp const P = move2dp(anchor, tree_.get_pos_trafo(frame, 0));
        return vec2dp{P.x / P.z, P.y / P.z, 1.0};
    }

    // Euclidean separation of two unitized points: the bulk norm of their difference (a
    // direction, z = 0)
    static value_t separation(vec2dp const& Pa, vec2dp const& Pb)
    {
        return value_t(bulk_nrm(vec2dp{Pa.x - Pb.x, Pa.y - Pb.y, 0.0}));
    }

    // unit direction from P_b to P_a (the rod axis of a distance constraint)
    static vec2dp rod_axis(vec2dp const& Pa, vec2dp const& Pb)
    {
        value_t const L = separation(Pa, Pb);
        return vec2dp{(Pa.x - Pb.x) / L, (Pa.y - Pb.y) / L, 0.0};
    }

    // Relative rotation angle of frame fa against frame fb, theta_a - theta_b, from the
    // world motor that carries frame b's orientation to frame a's, M_a (x) rrev(M_b).
    // The motor is rexp(1/2 B), so the angle is twice the log's rotational component;
    // in the plane that is the twist's z slot. Its time derivative is the difference of
    // the world angular velocities, which is what the frame kind's angular row of G
    // carries.
    //
    // Motors double-cover the motions: M and -M move everything identically, and a joint
    // passing through +-pi puts its motor on the far sheet, so the relative motor of two
    // frames with EQUAL orientation can arrive as -I. Its log would then report a full
    // turn (2 pi) as the residual, and the position projection would "close" it by
    // turning joints through whole revolutions -- measured as a sudden energy jump when a
    // welded arm swung through -pi. Normalising the sign so the identity-like component
    // (the pseudoscalar, the rgpr identity in the plane) is positive picks the near
    // sheet, and the angle comes out principal.
    value_t relative_rotation(size_t fa, size_t fb)
    {
        auto M = rgpr(tree_.get_pos_trafo(fa, 0), rrev(tree_.get_pos_trafo(fb, 0)));
        if (value_t(gr3(M)) < 0.0) M = -M;
        return 2.0 * rlog(M).z;
    }

    // write coordinate rates into the tree (per coordinate, then refresh each joint once)
    void write_rates(std::vector<coord> const& rc, std::vector<value_t> const& rates)
    {
        for (size_t k = 0; k < rc.size(); ++k)
            tree_.set_coord_rate(rc[k], rates[k]);
        for (auto const& c : rc)
            tree_.apply_joint_state(c.frame);
    }

    // Solve the bordered acceleration-level KKT system for q-ddot (and lambda) at the
    // current state -- see joint_accelerations() for the system. `rc` is the coordinate
    // list (all dof coordinates for dynamics). Reuses the open-loop M / tau
    // (assemble_mass_bias), the constraint Jacobian G, and the velocity-product term
    // G-dot q-dot (constraint_bias). Returns q-ddot (length n); writes lambda (length m)
    // if requested.
    std::vector<value_t> kkt_dynamics(std::vector<coord> const& rc,
                                      std::vector<value_t>* lambda_out)
    {
        size_t const n = rc.size();
        size_t const m = constraint_rows();

        // M (n*n) and tau (n) from the open-loop assembly; this also runs the bias pass
        // (zeroes the chain's relative accel twists), which is exactly the q-ddot = 0
        // state constraint_bias() needs below.
        auto const mb = tree_.assemble_mass_bias(rc);
        std::vector<value_t> const& M = mb.first;
        std::vector<value_t> const& tau = mb.second;

        std::vector<value_t> const G = constraint_jacobian(rc);
        std::vector<value_t> const gd = constraint_bias(); // G-dot q-dot (m)

        // bottom rhs of the saddle-point system is -G-dot q-dot
        std::vector<value_t> gbias(m);
        for (size_t c = 0; c < m; ++c)
            gbias[c] = -gd[c];

        // solve | M  Gᵀ | |q-ddot|   | tau   |  through the shared bordered LU
        //       | G  0  | |lambda| = | gbias |
        return hd::ga::kkt_solve(M, G, tau, gbias, n, m, lambda_out);
    }

    // Second time derivative of the residual at the CURRENT accel-twist configuration,
    // one block per active constraint (length m): for coincidence and the point rows of
    // frame the relative anchor acceleration a(P_a) - a(P_b); for distance the second
    // derivative of |d| - length,
    //
    //     n . (a_a - a_b) + (|v_rel|^2 - (n . v_rel)^2) / |d|      n = d / |d|
    //
    // (the transverse relative velocity's centripetal term); for the angular row of
    // frame the relative world angular acceleration. With every dof joint's relative
    // accel twist zeroed this is the velocity-product term G-dot q-dot (constraint_bias);
    // with the driven accelerations prescribed it is what solve_accelerations needs.
    std::vector<value_t> relative_accelerations()
    {
        std::vector<value_t> gd(constraint_rows(), 0.0);
        size_t r = 0;
        for (auto const& lc : loops_) {
            if (!lc.active) continue;
            vec2dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec2dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            vec2dp const aa = tree_.point_acceleration(Pa, lc.frame_a);
            vec2dp const ab = tree_.point_acceleration(Pb, lc.frame_b);
            switch (lc.type) {
                case constraint2dp::coincidence:
                    gd[r++] = aa.x - ab.x;
                    gd[r++] = aa.y - ab.y;
                    break;
                case constraint2dp::distance: {
                    vec2dp const n = rod_axis(Pa, Pb);
                    value_t const L = separation(Pa, Pb);
                    vec2dp const va = tree_.point_velocity(Pa, lc.frame_a);
                    vec2dp const vb = tree_.point_velocity(Pb, lc.frame_b);
                    value_t const vx = va.x - vb.x, vy = va.y - vb.y;
                    value_t const vn = n.x * vx + n.y * vy;
                    gd[r++] = n.x * (aa.x - ab.x) + n.y * (aa.y - ab.y) +
                              (vx * vx + vy * vy - vn * vn) / L;
                    break;
                }
                case constraint2dp::frame:
                    gd[r++] = aa.x - ab.x;
                    gd[r++] = aa.y - ab.y;
                    gd[r++] = tree_.accel_twist_world(lc.frame_a).z -
                              tree_.accel_twist_world(lc.frame_b).z;
                    break;
            }
        }
        return gd;
    }

    // Velocity-product term G-dot q-dot (length m): relative_accelerations() at
    // q-ddot = 0. With every dof joint's relative accel twist zeroed, point_acceleration
    // returns exactly the Coriolis/centripetal part. Pre: the velocities (relative
    // velocity twists) are current. (Idempotent with the bias pass assemble_mass_bias
    // just ran.)
    std::vector<value_t> constraint_bias()
    {
        for (size_t const j : tree_.dof_joints())
            tree_.set_accel_twist(j, twist2dp{0.0, 0.0, 0.0});
        return relative_accelerations();
    }

    // Project the coordinate rates onto the constraint tangent space (enforce G q-dot =
    // 0) without changing the feasible component: q-dot <- q-dot - G⁺(G q-dot), where the
    // min-norm pseudo-inverse correction G⁺(G q-dot) = Gᵀ(G Gᵀ)^-1 (G q-dot) comes from
    // the shared solve. Removes velocity drift after a step with no energy injection.
    void project_velocities(std::vector<coord> const& rc)
    {
        size_t const n = rc.size();
        size_t const m = constraint_rows();
        std::vector<value_t> const G = constraint_jacobian(rc);

        std::vector<value_t> qd(n), Gv(m, 0.0);
        for (size_t k = 0; k < n; ++k)
            qd[k] = tree_.coord_rate(rc[k]);
        for (size_t i = 0; i < m; ++i)
            for (size_t k = 0; k < n; ++k)
                Gv[i] += G[i * n + k] * qd[k];

        std::vector<value_t> const dqd = hd::ga::lstsq_solve(G, Gv, n); // G⁺(G q-dot)
        for (size_t k = 0; k < n; ++k)
            qd[k] -= dqd[k];
        write_rates(rc, qd);
    }

    // the dof coordinates NOT belonging to a `driven` joint -- the ones assemble() solves
    // for
    std::vector<coord> dependent_coords(std::vector<size_t> const& driven) const
    {
        std::vector<coord> dep;
        for (auto const& c : tree_.dof_coords()) {
            bool is_driven = false;
            for (size_t const d : driven)
                if (d == c.frame) {
                    is_driven = true;
                    break;
                }
            if (!is_driven) dep.push_back(c);
        }
        return dep;
    }

    // Constraint Jacobian G (m x ndep, row-major), m = constraint_rows(), columns
    // indexed by the coordinates `dep`, one row block per active constraint. The point
    // rows of column k are the relative partial velocity of the two anchors w.r.t.
    // coordinate dep[k]:
    //
    //     dP_a/dq_j - dP_b/dq_j = velocity_field(S_j, P_a) - velocity_field(S_j, P_b)
    //
    // (each term present only when the coordinate's joint is an ancestor of the
    // anchor's frame), where S_j is the coordinate's world screw -- the SAME column that
    // builds the open-loop mass matrix and bias forces (assemble_mass_bias). The
    // distance row is that relative velocity projected on the rod axis n; the angular
    // row of frame is the relative angular rate, the z slot of the world screw with the
    // same ancestor test.
    std::vector<value_t> constraint_jacobian(std::vector<coord> const& dep)
    {
        size_t const m = constraint_rows();
        size_t const ndep = dep.size();
        std::vector<value_t> G(m * ndep, 0.0);

        // precompute the world screws of the dependent coordinates
        std::vector<twist2dp> S(ndep);
        for (size_t k = 0; k < ndep; ++k)
            S[k] = tree_.world_screw(dep[k]);

        size_t r = 0;
        for (auto const& lc : loops_) {
            if (!lc.active) continue;
            vec2dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec2dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            vec2dp const n =
                (lc.type == constraint2dp::distance) ? rod_axis(Pa, Pb) : vec2dp{};
            for (size_t k = 0; k < ndep; ++k) {
                size_t const j = dep[k].frame;
                bool const on_a = tree_.is_ancestor(j, lc.frame_a);
                bool const on_b = tree_.is_ancestor(j, lc.frame_b);
                vec2dp dPa{0.0, 0.0, 0.0};
                vec2dp dPb{0.0, 0.0, 0.0};
                if (on_a) dPa = dynamic_system2dp::velocity_field(S[k], Pa);
                if (on_b) dPb = dynamic_system2dp::velocity_field(S[k], Pb);
                value_t const dx = dPa.x - dPb.x, dy = dPa.y - dPb.y;
                switch (lc.type) {
                    case constraint2dp::coincidence:
                        G[(r + 0) * ndep + k] = dx;
                        G[(r + 1) * ndep + k] = dy;
                        break;
                    case constraint2dp::distance:
                        G[r * ndep + k] = n.x * dx + n.y * dy;
                        break;
                    case constraint2dp::frame:
                        G[(r + 0) * ndep + k] = dx;
                        G[(r + 1) * ndep + k] = dy;
                        G[(r + 2) * ndep + k] =
                            (on_a ? S[k].z : 0.0) - (on_b ? S[k].z : 0.0);
                        break;
                }
            }
            r += rows_of(lc);
        }
        return G;
    }
};

} // namespace hd::ga::pga
