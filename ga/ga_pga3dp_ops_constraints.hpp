#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Closed-loop (parallel mechanism) layer for PGA3DP -- the 3D counterpart of
// ga_pga2dp_ops_constraints.hpp. Same classes, methods and algorithms; only the
// representation changes (vec3dp points with weight w, twist3dp = BiVec3dp screw axes,
// move3dp, mvec3dp_e motors, and 3 residual components per coincidence instead of 2, 6
// per frame instead of 3). The dimension-agnostic numeric kernels (lstsq_solve /
// kkt_solve, in detail/ga_solver.hpp) are shared verbatim with the 2D layer.
// closed_loop_system3dp HAS-A dynamic_system3dp (composition) as its spanning tree and
// closes kinematic loops by registering loop-closure constraints between existing tree
// frames.
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
//
// A loop closure is a JOINT between two tree frames (Featherstone ch. 8): its kind is
// the subspace of relative motion it removes, and its rows in g / G are that
// constraint-force subspace (see constraint3dp). A constraint can be switched off and
// on at run time (loop_constraint3dp::active), with the row count m following.
//
// The tree's generalised coordinates are dynamic_system3dp::coord -- one per screw of
// every dynamic joint, so a floating base (a free joint with mass) or a spherical joint
// takes part in the loop closure like a revolute one; the layer indexes q, q-dot, q-ddot
// and the Jacobian columns by those coordinates.

#include "detail/ga_solver.hpp"        // hd::ga::lstsq_solve / kkt_solve (shared solvers)
#include "ga_pga3dp_ops_mechanics.hpp" // dynamic_system3dp (the reused spanning tree)
#include "ga_usr_utilities.hpp"        // hd::ga::rk4_step (shared RK4 integrator)
#include "ga_value_t.hpp"              // value_t

#include <cmath>     // std::abs
#include <cstddef>   // size_t
#include <stdexcept> // std::runtime_error
#include <string>
#include <vector>


namespace hd::ga::pga {

// Loop-closure constraint kinds -- each a subspace of relative motion removed between
// the two anchors, i.e. the rows of g and G it contributes (Featherstone's constraint-
// force subspace T of the loop joint):
//
//   coincidence : the two anchor POINTS coincide (3 rows) -- a spherical / ball joint.
//                 Rotation about the point stays free.
//
//   distance    : the two anchors keep the DISTANCE `length` (1 row) -- a rigid rod, a
//                 Stewart-Gough leg of fixed length; the multiplier is the rod force.
//
//   frame       : a WELD -- the two bodies are rigidly attached to each other (6 rows).
//                 Position: anchor_a and anchor_b coincide in the world, as for
//                 coincidence. Orientation: frame_a and frame_b have the same world
//                 orientation (the anchors only locate the attachment point; they do
//                 not rotate the frames, so a weld between frames built with different
//                 orientations must be closed by assemble() or built to match). Every
//                 relative dof of the loop is removed -- the zero-dof loop joint: a tool
//                 fixed in a hand, a foot standing flat, a rigid coupler.
//
// Multiplier convention (joint_accelerations): the KKT system is M q-ddot + G^T lambda
// = tau, so the constraint force on the system is -G^T lambda. Per kind, on anchor a:
// coincidence -- the force -lambda (x, y, z); distance -- the force -lambda n along the
// rod axis n from anchor b to anchor a, so lambda > 0 is TENSION; frame -- the force
// -lambda_F (3) at the anchor and the couple -lambda_M (3, world axes). Anchor b
// receives the opposite.
//
// Point-on-line, point-on-plane and parallel/perpendicular are natural further members
// (vanishing wdg/rwdg/dot expressions), added when a mechanism needs them.
enum class constraint3dp { coincidence, distance, frame };

// Loop-closure constraint descriptor: anchor point anchor_a (given in frame_a's
// coordinates) must satisfy `type` against anchor point anchor_b (in frame_b's
// coordinates). For coincidence the two anchors must coincide in world coordinates -- the
// loop-closure condition for a Stewart-Gough leg pinned to the platform; for distance
// they keep `length` apart; for frame the two frames coincide at the anchors.
//
// `active` switches the constraint on and off at run time (contact state). An inactive
// constraint contributes no rows -- the system behaves exactly as if it had not been
// added -- while keeping its index.
//
// A pure-data aggregate (all public fields) with an fmt formatter, so it is also exposed
// in the Python wrapper (ga_py). Adding a field here widens the bound constructor
// (CLAUDE.md, "Adding a FIELD to an already-bound pure-data struct").
struct loop_constraint3dp {
    size_t frame_a;  // tree frame holding anchor a
    vec3dp anchor_a; // anchor point in frame_a coordinates (unitized, w = 1)
    size_t frame_b;  // tree frame holding anchor b
    vec3dp anchor_b; // anchor point in frame_b coordinates (unitized, w = 1)
    constraint3dp type{constraint3dp::coincidence};
    value_t length{0.0}; // the rod length (distance kind only)
    bool active{true};   // contributes rows to g / G when true
};


// A closed-loop / parallel mechanism: an open-chain spanning tree (the reused
// dynamic_system3dp) plus loop-closure constraints between tree frames. Build the tree
// through the familiar forwarded API, then close loops with add_loop_constraint(...).
class closed_loop_system3dp {

    using coord = dynamic_system3dp::coord;

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

    void add_spherical_body(static_frame3dp const& f, body3dp const& b,
                            vec3dp const& pivot_b,
                            size_t parent_idx = dynamic_system3dp::prev_frame)
    {
        tree_.add_spherical_body(f, b, pivot_b, parent_idx);
    }

    size_t index_of(std::string const& frame_name) const
    {
        return tree_.index_of(frame_name);
    }

    // --- close loops
    // --------------------------------------------------------------------

    // register a loop-closure constraint between two existing tree frames; returns its
    // index (stable for the life of the system, also while inactive)
    size_t add_loop_constraint(loop_constraint3dp const& c)
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

    // drive a 1-dof joint to a generalised coordinate q (set q, refresh the kinematic
    // state); a motor joint is posed through system().set_joint_motor()
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
    size_t active_loop_count() const
    {
        size_t n = 0;
        for (auto const& lc : loops_)
            if (lc.active) ++n;
        return n;
    }
    loop_constraint3dp const& loop(size_t c) const { return loops_.at(c); }

    // rows a constraint of this kind contributes to g / G (its constraint subspace)
    static size_t rows_of(loop_constraint3dp const& lc)
    {
        switch (lc.type) {
            case constraint3dp::coincidence:
                return 3;
            case constraint3dp::distance:
                return 1;
            case constraint3dp::frame:
                return 6;
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
    // independent constraints at THIS configuration (Featherstone 8.10: it varies as the
    // mechanism moves, and mobility = n - rank). Below constraint_rows() the multipliers
    // are only determined up to the dependence; joint_accelerations() then returns
    // their minimum-norm values.
    size_t constraint_rank()
    {
        auto const rc = tree_.dof_coords();
        size_t const m = constraint_rows();
        if (m == 0 || rc.empty()) return 0;
        return hd::ga::matrix_rank(constraint_jacobian(rc), m, rc.size());
    }

    // access the underlying spanning tree for diagnostics shared with the open-loop tier
    dynamic_system3dp& system() { return tree_; }
    dynamic_system3dp const& system() const { return tree_; }

    // --- constraint evaluation (position level)
    // -----------------------------------------

    // Residual g(q), one block per ACTIVE constraint in index order: coincidence -- the
    // unitized world-coordinate difference P_a - P_b (3 rows); distance -- |P_a - P_b| -
    // length (1 row); frame -- P_a - P_b (3 rows) and the relative rotation vector of the
    // two frames, axis * angle in world coordinates (3 rows). g == 0 <=> every active
    // loop is closed. Length constraint_rows().
    std::vector<value_t> residual()
    {
        std::vector<value_t> g(constraint_rows());
        size_t r = 0;
        for (auto const& lc : loops_) {
            if (!lc.active) continue;
            vec3dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec3dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            switch (lc.type) {
                case constraint3dp::coincidence:
                    g[r++] = Pa.x - Pb.x;
                    g[r++] = Pa.y - Pb.y;
                    g[r++] = Pa.z - Pb.z;
                    break;
                case constraint3dp::distance:
                    g[r++] = separation(Pa, Pb) - lc.length;
                    break;
                case constraint3dp::frame: {
                    g[r++] = Pa.x - Pb.x;
                    g[r++] = Pa.y - Pb.y;
                    g[r++] = Pa.z - Pb.z;
                    vec3dp const rot = relative_rotation(lc.frame_a, lc.frame_b);
                    g[r++] = rot.x;
                    g[r++] = rot.y;
                    g[r++] = rot.z;
                    break;
                }
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
    // on the constraint Jacobian G (the relative partial velocities velocity_field(S_j,
    // P_a) - velocity_field(S_j, P_b), reused from the open-loop spatial Jacobian). Each
    // Newton increment is applied as a coordinate step -- q += delta for a 1-dof joint,
    // the retraction M <- M (x) rexp(1/2 delta . S) for a motor joint. Returns the
    // achieved residual norm ‖g‖. Throws if Newton does not converge within max_iter.
    value_t assemble(std::vector<size_t> const& driven = {}, value_t tol = value_t(1e-12),
                     size_t max_iter = 50)
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
        throw std::runtime_error(
            std::string(
                "closed_loop_system3dp::assemble: Newton did not converge within ") +
            std::to_string(max_iter) + std::string(" iterations (residual ") +
            std::to_string(residual_norm()) + std::string(")"));
    }

    // --- velocity / acceleration distribution (kinematic closed loop)
    // -------------------

    // Distribute the driver joint rates to the dependent coordinate rates enforcing
    // G q-dot = 0:
    //     G_dep q-dot_dep = -G_drv q-dot_drv.
    // The solved dependent rates are written into the tree and returned. Pre:
    // assemble()d.
    std::vector<value_t> solve_velocities(std::vector<size_t> const& driven = {})
    {
        std::vector<coord> const drv = tree_.coords_of(driven);
        std::vector<coord> const dep = dependent_coords(driven);
        size_t const m = constraint_rows();

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

    // Distribute the driver joint accelerations to the dependent ones enforcing
    //     G_dep q-ddot_dep = -G_drv q-ddot_drv - G-dot q-dot,
    // with the velocity-product term read off as the second derivative of the residual
    // at dependent q-ddot = 0 (driven q-ddot prescribed) through the kinematic layer's
    // acceleration fields. `driven_accels` are the prescribed q-ddot of the driven
    // joints' coordinates (parallel to coords_of(driven); missing entries default to 0).
    // Pre: solve_velocities() has run.
    std::vector<value_t>
    solve_accelerations(std::vector<size_t> const& driven = {},
                        std::vector<value_t> const& driven_accels = {})
    {
        std::vector<coord> const drv = tree_.coords_of(driven);
        std::vector<coord> const dep = dependent_coords(driven);

        for (auto const& c : dep)
            tree_.set_accel_twist(c.frame, twist3dp{});
        for (auto const& c : drv)
            tree_.set_accel_twist(c.frame, twist3dp{});
        for (size_t k = 0; k < drv.size(); ++k) {
            value_t const qdd = (k < driven_accels.size()) ? driven_accels[k] : 0.0;
            tree_.set_accel_twist(drv[k].frame, tree_.relative_accel_twist(drv[k].frame) +
                                                    qdd * tree_.screw_of(drv[k]));
        }

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
    // (the constraint / leg forces). M(q), tau(q,q-dot) come from the open-loop assembly
    // (assemble_mass_bias), G from constraint_jacobian, and -G-dot q-dot from
    // constraint_bias. Returns q-ddot (dof_coords() order); writes lambda if requested --
    // one block per ACTIVE constraint, rows_of() entries each, with the sign convention
    // stated at constraint3dp (the force on anchor a is -lambda).
    std::vector<value_t> joint_accelerations(std::vector<value_t>* lambda_out = nullptr)
    {
        return kkt_dynamics(tree_.dof_coords(), lambda_out);
    }

    // Advance the closed-loop system by dt: RK4 (shared rk4_step) on the tree's
    // coordinate state with the constrained KKT solve at each sub-step, then post-step
    // GGL projection (position by min-norm Newton, velocity by q-dot <- q-dot - G⁺(G
    // q-dot)) for energy-clean drift control.
    void step(value_t dt)
    {
        auto const rc = tree_.dof_coords();
        size_t const n = rc.size();
        if (n == 0) return;

        std::vector<value_t> u;
        tree_.read_state(rc, u);
        auto f = [&](value_t, std::vector<value_t> const& uu, std::vector<value_t>& du) {
            tree_.write_state(rc, uu);
            auto const qdd = kkt_dynamics(rc, nullptr);
            tree_.state_rates(rc, uu, qdd, du);
        };
        rk4_integrator(2 * n).step(f, u, tree_.time(), dt);
        tree_.commit_state(rc, u);

        if (constraint_rows() > 0) {
            assemble(/*driven*/ {}); // position: min-norm Newton -> g ~ 0
            project_velocities(rc);  // velocity:  q-dot <- q-dot - G⁺(G q-dot)
        }
    }

  private:

    // world-coordinate, unitized (w = 1) position of an anchor point of a tree frame
    vec3dp anchor_world(size_t frame, vec3dp const& anchor)
    {
        return unitize(move3dp(anchor, tree_.get_pos_trafo(frame, 0)));
    }

    // Euclidean separation of two unitized points: the bulk norm of their difference (a
    // direction, w = 0)
    static value_t separation(vec3dp const& Pa, vec3dp const& Pb)
    {
        return value_t(bulk_nrm(vec3dp{Pa.x - Pb.x, Pa.y - Pb.y, Pa.z - Pb.z, 0.0}));
    }

    // unit direction from P_b to P_a (the rod axis of a distance constraint)
    static vec3dp rod_axis(vec3dp const& Pa, vec3dp const& Pb)
    {
        value_t const L = separation(Pa, Pb);
        return vec3dp{(Pa.x - Pb.x) / L, (Pa.y - Pb.y) / L, (Pa.z - Pb.z) / L, 0.0};
    }

    // Relative rotation of frame fa against frame fb as a world-frame rotation vector
    // (axis * angle, w = 0), from the world motor that carries frame b's orientation to
    // frame a's, M_a (x) rrev(M_b). The motor is rexp(1/2 B), so the rotation vector is
    // twice the log's weight (direction) part -- the translation lives in the bulk. Its
    // time derivative at closure is the difference of the world angular velocities,
    // which is what the frame kind's angular rows of G carry.
    //
    // Motors double-cover the motions: M and -M move everything identically, and a joint
    // passing through +-pi puts its motor on the far sheet, so the relative motor of two
    // frames with EQUAL orientation can arrive as -1. Its log would then report a full
    // turn as the residual and the position projection would "close" it by turning
    // joints through whole revolutions (measured in 2D as a sudden energy jump).
    // Normalising the sign so the identity-like component is positive picks the near
    // sheet -- and that component is the PSEUDOSCALAR, the rgpr identity (the scalar
    // part of the identity motor is zero here; testing gr0 instead broke the closure).
    vec3dp relative_rotation(size_t fa, size_t fb)
    {
        auto M = rgpr(tree_.get_pos_trafo(fa, 0), rrev(tree_.get_pos_trafo(fb, 0)));
        if (value_t(gr4(M)) < 0.0) M = -M;
        twist3dp const B = rlog(M);
        return vec3dp{2.0 * B.vx, 2.0 * B.vy, 2.0 * B.vz, 0.0};
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
    // current state -- see joint_accelerations(). `rc` is the coordinate list (all dof
    // coordinates for dynamics). Reuses the open-loop M / tau (assemble_mass_bias), G
    // (constraint_jacobian) and -G-dot q-dot (constraint_bias).
    std::vector<value_t> kkt_dynamics(std::vector<coord> const& rc,
                                      std::vector<value_t>* lambda_out)
    {
        size_t const n = rc.size();
        size_t const m = constraint_rows();

        // M (n*n) and tau (n) from the open-loop assembly; this also runs the bias pass
        // (zeroes the chain's relative accel twists), the q-ddot = 0 state
        // constraint_bias() needs below.
        auto const mb = tree_.assemble_mass_bias(rc);
        std::vector<value_t> const& M = mb.first;
        std::vector<value_t> const& tau = mb.second;

        std::vector<value_t> const G = constraint_jacobian(rc);
        std::vector<value_t> const gd = constraint_bias(); // G-dot q-dot (m)

        std::vector<value_t> gbias(m);
        for (size_t c = 0; c < m; ++c)
            gbias[c] = -gd[c];

        return hd::ga::kkt_solve(M, G, tau, gbias, n, m, lambda_out);
    }

    // Second time derivative of the residual at the CURRENT accel-twist configuration,
    // one block per active constraint (length m): the relative anchor acceleration for
    // the point rows; for distance
    //
    //     n . (a_a - a_b) + (|v_rel|^2 - (n . v_rel)^2) / |d|      n = d / |d|
    //
    // (the transverse relative velocity's centripetal term); for the angular rows of
    // frame the relative world angular acceleration (the weight slots of the world
    // acceleration twists). With every dof joint's relative accel twist zeroed this is
    // G-dot q-dot (constraint_bias); with the driven accelerations prescribed it is what
    // solve_accelerations needs.
    std::vector<value_t> relative_accelerations()
    {
        std::vector<value_t> gd(constraint_rows(), 0.0);
        size_t r = 0;
        for (auto const& lc : loops_) {
            if (!lc.active) continue;
            vec3dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec3dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            vec3dp const aa = tree_.point_acceleration(Pa, lc.frame_a);
            vec3dp const ab = tree_.point_acceleration(Pb, lc.frame_b);
            switch (lc.type) {
                case constraint3dp::coincidence:
                    gd[r++] = aa.x - ab.x;
                    gd[r++] = aa.y - ab.y;
                    gd[r++] = aa.z - ab.z;
                    break;
                case constraint3dp::distance: {
                    vec3dp const n = rod_axis(Pa, Pb);
                    value_t const L = separation(Pa, Pb);
                    vec3dp const va = tree_.point_velocity(Pa, lc.frame_a);
                    vec3dp const vb = tree_.point_velocity(Pb, lc.frame_b);
                    value_t const vx = va.x - vb.x, vy = va.y - vb.y, vz = va.z - vb.z;
                    value_t const vn = n.x * vx + n.y * vy + n.z * vz;
                    gd[r++] = n.x * (aa.x - ab.x) + n.y * (aa.y - ab.y) +
                              n.z * (aa.z - ab.z) +
                              (vx * vx + vy * vy + vz * vz - vn * vn) / L;
                    break;
                }
                case constraint3dp::frame: {
                    gd[r++] = aa.x - ab.x;
                    gd[r++] = aa.y - ab.y;
                    gd[r++] = aa.z - ab.z;
                    twist3dp const Aa = tree_.accel_twist_world(lc.frame_a);
                    twist3dp const Ab = tree_.accel_twist_world(lc.frame_b);
                    gd[r++] = Aa.vx - Ab.vx;
                    gd[r++] = Aa.vy - Ab.vy;
                    gd[r++] = Aa.vz - Ab.vz;
                    break;
                }
            }
        }
        return gd;
    }

    // Velocity-product term G-dot q-dot (length m): relative_accelerations() at
    // q-ddot = 0. With every dof joint's relative accel twist zeroed, point_acceleration
    // returns exactly the Coriolis/centripetal part.
    std::vector<value_t> constraint_bias()
    {
        for (size_t const j : tree_.dof_joints())
            tree_.set_accel_twist(j, twist3dp{});
        return relative_accelerations();
    }

    // Project the coordinate rates onto the constraint tangent space (enforce G q-dot =
    // 0): q-dot <- q-dot - G⁺(G q-dot). Removes velocity drift with no energy injection.
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

    // the dof coordinates NOT belonging to a `driven` joint -- the ones assemble() /
    // the solvers solve for
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
    // rows of frame are the relative angular rate, the weight slots (vx, vy, vz) of the
    // world screw with the same ancestor test.
    std::vector<value_t> constraint_jacobian(std::vector<coord> const& dep)
    {
        size_t const m = constraint_rows();
        size_t const ndep = dep.size();
        std::vector<value_t> G(m * ndep, 0.0);

        std::vector<twist3dp> S(ndep);
        for (size_t k = 0; k < ndep; ++k)
            S[k] = tree_.world_screw(dep[k]);

        size_t r = 0;
        for (auto const& lc : loops_) {
            if (!lc.active) continue;
            vec3dp const Pa = anchor_world(lc.frame_a, lc.anchor_a);
            vec3dp const Pb = anchor_world(lc.frame_b, lc.anchor_b);
            vec3dp const n =
                (lc.type == constraint3dp::distance) ? rod_axis(Pa, Pb) : vec3dp{};
            for (size_t k = 0; k < ndep; ++k) {
                size_t const j = dep[k].frame;
                bool const on_a = tree_.is_ancestor(j, lc.frame_a);
                bool const on_b = tree_.is_ancestor(j, lc.frame_b);
                vec3dp dPa{};
                vec3dp dPb{};
                if (on_a) dPa = dynamic_system3dp::velocity_field(S[k], Pa);
                if (on_b) dPb = dynamic_system3dp::velocity_field(S[k], Pb);
                value_t const dx = dPa.x - dPb.x, dy = dPa.y - dPb.y, dz = dPa.z - dPb.z;
                switch (lc.type) {
                    case constraint3dp::coincidence:
                        G[(r + 0) * ndep + k] = dx;
                        G[(r + 1) * ndep + k] = dy;
                        G[(r + 2) * ndep + k] = dz;
                        break;
                    case constraint3dp::distance:
                        G[r * ndep + k] = n.x * dx + n.y * dy + n.z * dz;
                        break;
                    case constraint3dp::frame: {
                        G[(r + 0) * ndep + k] = dx;
                        G[(r + 1) * ndep + k] = dy;
                        G[(r + 2) * ndep + k] = dz;
                        value_t const sa = on_a ? 1.0 : 0.0, sb = on_b ? 1.0 : 0.0;
                        G[(r + 3) * ndep + k] = (sa - sb) * S[k].vx;
                        G[(r + 4) * ndep + k] = (sa - sb) * S[k].vy;
                        G[(r + 5) * ndep + k] = (sa - sb) * S[k].vz;
                        break;
                    }
                }
            }
            r += rows_of(lc);
        }
        return G;
    }
};

} // namespace hd::ga::pga
