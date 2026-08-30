#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Unilateral ground contact for PGA3DP -- a SEPARATE, ADDITIVE layer on top of
// closed_loop_system3dp, for the feet of a walking or hopping mechanism (the one
// surface such a model touches). A contact is an INEQUALITY, in two halves
// (Featherstone, Rigid Body Dynamics Algorithms, ch. 11): the point may not go below
// the ground, h(P) >= 0, and the ground can only push, the normal reaction >= 0. Both
// are handled by EVENTS on the equality machinery of the closed-loop layer:
//
//   - touchdown: a contact point crossing the ground line downward is located by
//     bisection on the last step (the system is a value type, so a copy is stepped
//     from the pre-event state), the loop constraint is anchored at the crossing point
//     projected onto the ground and switched on WITH THE IMPACT MAP
//     (activate_loop_with_impact): the closing velocity is absorbed, plastically;
//
//   - lift-off: an active contact whose normal reaction turns to a pull of more than
//     release_threshold (0.5 N: a foot at zero load neither pulls nor pushes, and
//     round-off must not release it a thousand times) is released where it is
//     (Selig's criterion, lambda_y -> 0) -- provided it would SEPARATE: a copy of
//     the system with the constraint removed is stepped briefly and the contact is
//     released only if its point moves away from the ground (a contact that pulls
//     but would sink straight back is a RESTING contact and stays; without this test
//     a landed foot at zero load was released and re-landed every millisecond); it
//     left with a velocity and may land again at once (a hop, with its impact and
//     events). A contact released by hand (release()) that only sags lands again just
//     the same: it is a resting contact, and the separation test keeps it there until
//     something really lifts it.
//
// Two contact KINDS: a `point` contact is a coincidence loop (a ball joint, 3 rows: the
// foot may rotate about it), a `flat` contact is a frame loop (a weld to a fixed
// ground frame posed under the foot in its orientation, 6 rows: a flat foot whose
// ankle carries a moment). The flat kind reports its centre of pressure -- the moment
// about the two in-plane axes over the normal force, a 2-vector in the ground plane --
// and flags TIPPING when that leaves the foot rectangle (|cop_x| > half_length or
// |cop_y| > half_width); the transition to an edge contact is left to the caller.
//
// The ground is a PLANE (trivec3dp), so a slope is as natural as a floor: the signed
// height of a point is the incidence wdg(P, Pi) (a point on the plane has zero, above
// = the side the plane's normal points to, P0 -> P1 -> P2 right-handed), the landing
// point is project_onto(P, Pi), and the upward normal is the height's gradient. EVERY
// point of a mechanism that can reach the ground is registered as a contact -- feet,
// knees, the corners of a torso: a falling mechanism then lands on whatever touches
// first and lies on the ground, instead of sinking through it.

#include "ga_pga3dp_ops_constraints.hpp" // closed_loop_system3dp
#include "ga_usr_consts.hpp"             // O_3dp
#include "ga_usr_types.hpp"              // value_t type aliases
#include "ga_value_t.hpp"                // value_t

#include <cmath>   // std::abs, std::sqrt
#include <cstddef> // size_t
#include <string>
#include <vector>

namespace hd::ga::pga {

// provides unilateral ground-contact operations:
//
// - contact_kind3dp                     -> point (a pin) / flat (a weld: the flat foot)
// - ground_contact3dp                   -> the event layer: step(), heights, reactions,
//                                          the centre of pressure, the unilateral switch;
//                                          its nested spec (which frame point touches,
//                                          of which kind) and event (a touchdown /
//                                          lift-off record) -- nested so the binding
//                                          generator leaves them with the class

enum class contact_kind3dp { point, flat };

class ground_contact3dp {

  public:

    struct spec {
        size_t frame{0};                    // the frame carrying the contact point
        vec3dp point_b{0.0, 0.0, 0.0, 1.0}; // the point in that frame's coordinates
        contact_kind3dp kind{contact_kind3dp::point};
        value_t half_length{0.0}; // flat: the foot's half extents along its own x
        value_t half_width{0.0};  // and y (the centre of pressure may not leave the
                                  // rectangle -- tipping otherwise)
        bool unilateral{true};    // false: a BILATERAL pin that may pull (a clamped
                                  // foot, or a held-torso statics whose pins pull by
                                  // construction) -- it lands and can be released by
                                  // hand, but the lift-off rule leaves it alone
    };

    struct event {
        value_t t{};                  // simulation time of the event
        size_t contact{};             // which contact
        bool touchdown{};             // true: touchdown (with impact), false: lift-off
        vec3dp point{};               // the world point it happened at
        std::vector<value_t> impulse; // touchdown: the impact map's Lambda (all rows)
    };

    // the ground line and the tree frame the contacts are anchored to (the world root)
    ground_contact3dp(closed_loop_system3dp& cl, trivec3dp const& ground,
                      size_t ground_frame = 0) :
        cl_{&cl}, Pi_{unitize(ground)}, ground_frame_{ground_frame}
    {
        // the upward normal is the gradient of the height, which is linear in P
        n_ = vec3dp{height_of(vec3dp{1.0, 0.0, 0.0, 0.0}),
                    height_of(vec3dp{0.0, 1.0, 0.0, 0.0}),
                    height_of(vec3dp{0.0, 0.0, 1.0, 0.0}), 0.0};
        value_t const nn = std::sqrt(n_.x * n_.x + n_.y * n_.y + n_.z * n_.z);
        n_ = vec3dp{n_.x / nn, n_.y / nn, n_.z / nn, 0.0};
    }

    // the ground plane through three points; ABOVE the ground is the side the normal
    // of P0 -> P1 -> P2 (right-handed) points to (a floor z = z0: (0,0,z0), (1,0,z0),
    // (0,1,z0))
    static trivec3dp ground_plane(vec3dp const& P0, vec3dp const& P1, vec3dp const& P2)
    {
        return unitize(wdg(wdg(P0, P1), P2));
    }

    // signed height of a world point above the ground: the incidence of the point
    // with the plane (zero on it, positive above)
    value_t height(vec3dp const& P_world) const { return height_of(unitize(P_world)); }

    vec3dp const& normal() const { return n_; } // the upward unit normal
    trivec3dp const& ground() const { return Pi_; }

    // register a contact (its loop constraints are added to the system, inactive);
    // returns the contact index
    size_t add(spec const& spec)
    {
        contact c;
        c.spec = spec;
        vec3dp const P = point_world(spec.frame, spec.point_b);
        c.pin = cl_->add_loop_constraint(loop_constraint3dp{
            spec.frame, spec.point_b, ground_frame_, P, constraint3dp::coincidence});
        cl_->set_loop_active(c.pin, false);
        if (spec.kind == contact_kind3dp::flat) {
            std::string const name = "ground_contact_" + std::to_string(contacts_.size());
            cl_->add_frame(static_frame3dp(name, P), ground_frame_);
            c.ground = cl_->index_of(name);
            c.weld = cl_->add_loop_constraint(loop_constraint3dp{
                spec.frame, spec.point_b, c.ground, O_3dp, constraint3dp::frame});
            cl_->set_loop_active(c.weld, false);
        }
        contacts_.push_back(c);
        return contacts_.size() - 1;
    }

    // Advance the system by dt with the contact events resolved in order: each
    // touchdown is located by bisection on the sub-interval it happens in, applied,
    // and the remainder of dt is integrated on -- so a landing is exact in time rather
    // than late by up to one step. Lift-offs are checked at the end of the step.
    void step(value_t dt, value_t tol = value_t(1e-10))
    {
        value_t remaining = dt;
        int guard = 0;
        while (remaining > 0.0 && guard++ < 16) {
            closed_loop_system3dp const before = *cl_;
            std::vector<value_t> h0(contacts_.size());
            for (size_t c = 0; c < contacts_.size(); ++c)
                h0[c] = contact_height(c);
            cl_->step(remaining);
            // the first contact that crossed the ground downward in this interval
            size_t hit = contacts_.size();
            value_t tau = remaining;
            for (size_t c = 0; c < contacts_.size(); ++c) {
                // (h0 >= -1e-9: a foot released at the ground line sits at a height
                // of round-off size, possibly negative; it must still be able to land)
                if (contacts_[c].active || h0[c] < -1e-9 || contact_height(c) >= 0.0)
                    continue;
                // bisection on tau in (0, remaining]: not crossed at 0, crossed at
                // remaining; a copy of the pre-event state is stepped by tau
                value_t lo = 0.0, hi = remaining;
                for (int i = 0; i < 60 && hi - lo > tol; ++i) {
                    value_t const mid = 0.5 * (lo + hi);
                    closed_loop_system3dp trial = before;
                    trial.step(mid);
                    if (height(trial_point(trial, c)) < 0.0) hi = mid;
                    else lo = mid;
                }
                if (hi < tau) {
                    tau = hi;
                    hit = c;
                }
            }
            if (hit == contacts_.size()) break; // no event: the whole interval is done
            *cl_ = before;
            cl_->step(tau);
            touchdown(hit);
            remaining -= tau;
        }
        // lift-off: an active contact whose normal reaction pulls is released
        read_reactions();
        for (size_t c = 0; c < contacts_.size(); ++c)
            if (contacts_[c].active && contacts_[c].spec.unilateral &&
                contacts_[c].normal_force <= -release_threshold_ && would_separate(c))
                release(c);
    }

    // switch an ACTIVE contact's kind at run time (a pin becoming a flat foot when the
    // other leg lifts, and back): flat engages the weld with the impact map (the
    // foot's angular velocity is absorbed), point drops the angular row
    void set_kind(size_t idx, contact_kind3dp kind)
    {
        contact& c = contacts_[idx];
        if (c.weld == npos && kind == contact_kind3dp::flat) {
            // a point contact promoted to flat: create its weld on demand
            std::string const name = "ground_contact_" + std::to_string(idx);
            cl_->add_frame(static_frame3dp(name, contact_point(idx)), ground_frame_);
            c.ground = cl_->index_of(name);
            c.weld = cl_->add_loop_constraint(loop_constraint3dp{
                c.spec.frame, c.spec.point_b, c.ground, O_3dp, constraint3dp::frame});
            cl_->set_loop_active(c.weld, false);
        }
        c.spec.kind = kind;
        if (!c.active) return;
        if (kind == contact_kind3dp::flat) {
            pose_ground_frame(idx, contact_point(idx));
            cl_->set_loop_active(c.pin, false);
            cl_->activate_loop_with_impact(c.weld, 0.0);
        }
        else {
            cl_->set_loop_anchors(c.pin, c.spec.point_b, contact_point(idx));
            if (c.weld != npos) cl_->set_loop_active(c.weld, false);
            cl_->set_loop_active(c.pin, true);
        }
        read_reactions();
    }

    // engage a contact where its point is NOW, without an impact -- the initial
    // placement of a foot standing on the ground (the point must be on the ground to
    // round-off; a point above or below it would jump). Recorded as a touchdown event
    // with no impulse.
    void engage(size_t idx)
    {
        contact& c = contacts_[idx];
        if (c.active) return;
        vec3dp const at = unitize(project_onto(contact_point(idx), Pi_));
        if (c.spec.kind == contact_kind3dp::flat) {
            pose_ground_frame(idx, at);
            cl_->set_loop_active(c.weld, true);
        }
        else {
            cl_->set_loop_anchors(c.pin, c.spec.point_b, at);
            cl_->set_loop_active(c.pin, true);
        }
        c.active = true;
        events_.push_back({cl_->system().time(), idx, true, at, {}});
        read_reactions();
    }

    // release a contact by hand (a foot lifted by its muscles rather than pulled off)
    void release(size_t idx)
    {
        contact& c = contacts_[idx];
        if (!c.active) return;
        cl_->set_loop_active(c.pin, false);
        if (c.weld != npos) cl_->set_loop_active(c.weld, false);
        c.active = false;
        c.normal_force = 0.0;
        c.moment = vec3dp{0.0, 0.0, 0.0, 0.0};
        events_.push_back({cl_->system().time(), idx, false, contact_point(idx), {}});
    }

    // --- queries ---------------------------------------------------------------

    size_t count() const { return contacts_.size(); }
    bool active(size_t idx) const { return contacts_[idx].active; }
    contact_kind3dp kind(size_t idx) const { return contacts_[idx].spec.kind; }
    vec3dp contact_point(size_t idx) const
    {
        return point_world(contacts_[idx].spec.frame, contacts_[idx].spec.point_b);
    }
    value_t contact_height(size_t idx) const { return height(contact_point(idx)); }

    // the reaction the ground exerts on the contact point (world), from the last
    // multipliers read: the force, its normal component, and (flat) the moment
    // vector (world)
    vec3dp force(size_t idx) const { return contacts_[idx].force; }
    value_t normal_force(size_t idx) const { return contacts_[idx].normal_force; }
    vec3dp moment(size_t idx) const { return contacts_[idx].moment; }

    // flat: the centre of pressure in the FOOT's own axes (x, y of the contact frame),
    // the in-plane moment over the normal force (0 when unloaded), and whether it has
    // left the foot rectangle
    vec3dp cop(size_t idx) const
    {
        auto const& c = contacts_[idx];
        if (c.normal_force <= 0.0) return vec3dp{0.0, 0.0, 0.0, 0.0};
        // the moment about the foot's x axis shifts the pressure along y and vice
        // versa: cop = (M_y, -M_x) / N in the foot frame
        auto const M = cl_->system().get_pos_trafo(c.spec.frame, 0);
        vec3dp const mb = move3dp(c.moment, rrev(M)); // the moment in the foot frame
        return vec3dp{mb.y / c.normal_force, -mb.x / c.normal_force, 0.0, 0.0};
    }
    bool tipping(size_t idx) const
    {
        auto const& c = contacts_[idx];
        if (!c.active || c.spec.kind != contact_kind3dp::flat) return false;
        vec3dp const p = cop(idx);
        return std::abs(p.x) > c.spec.half_length || std::abs(p.y) > c.spec.half_width;
    }

    // The net wrench the ground exerts, summed over the ACTIVE contacts: each
    // contact's force line wdg(P_i, f_i), plus -- for a flat foot -- its ankle
    // couple. A couple is a weight-zero bivector (its moment is the same about EVERY
    // point), which is what makes it addable to a force line without a reference
    // point. The sum is ONE bivector carrying both resultants, the force as its
    // weight (att) and the moment as its bulk. Same expression as the 2D case.
    bivec3dp reaction_wrench() const
    {
        bivec3dp W{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        for (size_t i = 0; i < contacts_.size(); ++i) {
            if (!contacts_[i].active) continue;
            auto const& m = contacts_[i].moment;
            W += wdg(contact_point(i), contacts_[i].force) +
                 bivec3dp{0.0, 0.0, 0.0, m.x, m.y, m.z};
        }
        return W;
    }

    // The zero-moment point: the point of the ground where the net ground reaction
    // exerts no moment TANGENTIAL to it -- the MEET of the reaction wrench with the
    // ground, exactly as in 2D and by the same one-line expression:
    //
    //     ZMP = rwdg(reaction_wrench(), ground())
    //
    // The 3D definition is weaker than the 2D one and the meet tracks the difference
    // by itself. Two feet pushing and rubbing sum to a SCREW, not to a force line:
    // its moment cannot be made to vanish at any point, so there is no line of action
    // to intersect, and the classical definition accordingly asks only that the two
    // GROUND-TANGENTIAL components vanish -- the ground carries a torque about its
    // own normal through friction. The meet delivers precisely that, because
    //
    //     rwdg(couple along n, ground()) == 0        (a TANGENTIAL couple does not)
    //
    // so the regressive product annihilates exactly the one moment component the
    // definition permits to survive and no other. Measured on a frictional
    // two-contact stance: the full moment at the returned point is purely normal
    // (0, 0, -0.189 N m), its tangential part 1e-16, and 5 cm away 16.5 N m.
    // Nothing is projected out by hand and the same code stands in both dimensions;
    // in 2D there is no normal component to discard in the first place. Exact on a
    // SLOPE, since ground() is whatever plane was built. Tedrake's centre-of-pressure
    // formula sum(p_i N_i)/sum(N_i) is the flat-ground, frictionless special case.
    //
    // Returned unitized (w == 1) whenever the ground carries load. With the
    // resultant parallel to the ground there is no such point and the raw meet is
    // returned at infinity (w == 0) -- test .w before using .x/.y/.z.
    vec3dp zmp() const
    {
        vec3dp const X = rwdg(reaction_wrench(), Pi_);
        return (std::abs(X.w) > eps) ? unitize(X) : X;
    }

    // switch a contact between unilateral (lets go when it pulls) and bilateral (a pin
    // that may pull: released by hand only) at run time -- a foot about to be lifted
    // by its muscles is made unilateral so the lift-off rule lets it go
    void set_unilateral(size_t idx, bool on) { contacts_[idx].spec.unilateral = on; }
    bool unilateral(size_t idx) const { return contacts_[idx].spec.unilateral; }
    // the pull (a negative normal reaction, in force units) at which a unilateral
    // contact is released; 0 releases at the first tensile round-off
    void set_release_threshold(value_t f) { release_threshold_ = f; }
    value_t release_threshold() const { return release_threshold_; }

    // the loop-constraint indices behind a contact: the coincidence pin, and (flat)
    // the frame weld -- size_t(-1) while none has been created. A caller building
    // per-row quantities on the closed-loop layer (a commanded reaction lambda_d)
    // finds the contact's lambda block at constraint_row_offset(pin_loop(idx)).
    size_t pin_loop(size_t idx) const { return contacts_[idx].pin; }
    size_t weld_loop(size_t idx) const { return contacts_[idx].weld; }
    std::vector<event> const& events() const { return events_; }

    // refresh the reactions from the system's current state (step() does it)
    void read_reactions()
    {
        std::vector<value_t> lam;
        cl_->joint_accelerations(&lam);
        size_t r = 0;
        for (size_t k = 0; k < cl_->loop_count(); ++k) {
            if (!cl_->loop_active(k)) continue;
            bool const weld = cl_->loop(k).type == constraint3dp::frame;
            size_t const rows = weld ? 6 : 3;
            for (auto& c : contacts_) {
                if (k != c.pin && k != c.weld) continue;
                c.force = vec3dp{-lam[r], -lam[r + 1], -lam[r + 2], 0.0}; // on the foot
                c.normal_force = c.force.x * n_.x + c.force.y * n_.y + c.force.z * n_.z;
                c.moment = weld ? vec3dp{-lam[r + 3], -lam[r + 4], -lam[r + 5], 0.0}
                                : vec3dp{0.0, 0.0, 0.0, 0.0};
            }
            r += rows;
        }
    }

  private:

    static constexpr size_t npos = size_t(-1);

    struct contact {
        spec spec;
        size_t pin{npos};    // the coincidence loop
        size_t weld{npos};   // the frame loop (flat), npos if none yet
        size_t ground{npos}; // the fixed ground frame of the weld
        bool active{false};
        vec3dp force{0.0, 0.0, 0.0, 0.0};
        value_t normal_force{0.0};
        vec3dp moment{0.0, 0.0, 0.0, 0.0};
    };

    // the separation test of the lift-off rule: does the contact point move away from
    // the ground once its constraint is removed? (a short trial step on a copy)
    bool would_separate(size_t idx) const
    {
        contact const& c = contacts_[idx];
        closed_loop_system3dp trial = *cl_;
        trial.set_loop_active(c.pin, false);
        if (c.weld != npos) trial.set_loop_active(c.weld, false);
        value_t const h0 = height(trial_point(trial, idx));
        trial.step(separation_dt_);
        return height(trial_point(trial, idx)) > h0 + value_t(1e-9);
    }

    value_t height_of(vec3dp const& P) const { return value_t(wdg(P, Pi_)); }

    vec3dp point_world(size_t frame, vec3dp const& point_b) const
    {
        return unitize(move3dp(point_b, cl_->system().get_pos_trafo(frame, 0)));
    }
    vec3dp trial_point(closed_loop_system3dp& trial, size_t idx) const
    {
        auto const& s = contacts_[idx].spec;
        return unitize(move3dp(s.point_b, trial.system().get_pos_trafo(s.frame, 0)));
    }

    // the weld's ground frame: at the landing point, in the contact frame's orientation
    // (so the weld engages without a jump)
    void pose_ground_frame(size_t idx, vec3dp const& at)
    {
        contact const& c = contacts_[idx];
        auto M = cl_->system().get_pos_trafo(c.spec.frame, 0);
        if (value_t(gr4(M)) < 0.0) M = -M; // the near sheet of the double cover
        pose3dp p = pose3dp_from_motor(M);
        p.origin = at;
        cl_->system().set_pose(c.ground, p);
    }

    // touchdown of contact idx at the current state (its point on or just below the
    // ground): anchor at the projection onto the ground line, engage with the impact
    void touchdown(size_t idx)
    {
        contact& c = contacts_[idx];
        vec3dp const P = contact_point(idx);
        vec3dp const landing = unitize(project_onto(P, Pi_));
        std::vector<value_t> Lam;
        if (c.spec.kind == contact_kind3dp::flat) {
            pose_ground_frame(idx, landing);
            Lam = cl_->activate_loop_with_impact(c.weld, 0.0);
        }
        else {
            cl_->set_loop_anchors(c.pin, c.spec.point_b, landing);
            Lam = cl_->activate_loop_with_impact(c.pin, 0.0);
        }
        c.active = true;
        events_.push_back({cl_->system().time(), idx, true, landing, Lam});
    }

    closed_loop_system3dp* cl_;
    trivec3dp Pi_;
    vec3dp n_{0.0, 0.0, 1.0, 0.0};
    size_t ground_frame_;
    std::vector<contact> contacts_;
    std::vector<event> events_;
    value_t release_threshold_{value_t(0.5)};
    value_t separation_dt_{value_t(1e-4)}; // the trial step of would_separate
};

} // namespace hd::ga::pga
