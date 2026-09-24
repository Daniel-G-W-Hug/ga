#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Unilateral ground contact for PGA2DP -- a SEPARATE, ADDITIVE layer on top of
// closed_loop_system2dp, for the feet of a walking or hopping mechanism (the one
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
// Two contact KINDS: a `point` contact is a coincidence loop (a pin: the foot may
// rotate about it), a `flat` contact is a frame loop (a weld to a fixed ground frame
// posed under the foot at its heading: a flat foot whose ankle carries a moment).
// The flat kind reports its centre of pressure, moment / normal force, and flags
// TIPPING when that leaves the foot (|cop| > half_length); the transition to a point
// contact at the edge is left to the caller (a gait decides what a tipping foot does).
//
// A SEGMENT contact (spec.segment, a second point end_b) is a rigid segment standing on
// the ground line -- a sole, an edge, a rod -- and the layer runs its transitions: while
// only one end is down that end is a pin; when the other end lands, or both land within
// the simultaneity window, the segment is WELDED flat at its first point, at the exact
// touchdown and inside the same impact map; a flat segment whose centre of pressure lies
// beyond one of its ends ROLLS onto a pin at that end (a lift-off of the other end); and
// a segment lifts off by the rule of every contact. Two pins on one rigid segment would
// over-constrain it -- four rows on three degrees of freedom, exactly dependent while the
// segment is level and NEARLY dependent once it tilts by a hair, where the solve returns
// enormous reactions -- so the modes never hold more than three rows.
//
// The ground is a LINE (bivec2dp), so a slope is as natural as a floor: the signed
// height of a point is the incidence wdg(P, L) (a point on the line has zero, above =
// to the left of the line's direction P0 -> P1), the landing point is project_onto(P,
// L), and the upward normal is the height's gradient. EVERY point of a mechanism that
// can reach the ground is registered as a contact -- feet, knees, the corners of a
// torso: a falling mechanism then lands on whatever touches first and lies on the
// ground, instead of sinking through it.

#include "../ga_usr_consts.hpp"          // O_2dp
#include "../ga_usr_types.hpp"           // value_t type aliases
#include "../ga_value_t.hpp"             // value_t
#include "ga_pga2dp_ops_constraints.hpp" // closed_loop_system2dp

#include <algorithm> // std::min
#include <cmath>     // std::abs, std::sqrt
#include <cstddef>   // size_t
#include <stdexcept> // std::invalid_argument
#include <string>
#include <vector>

/////////////////////////////////////////////////////////////////////////////////////////
// provides unilateral ground-contact operations:
//
// - contact_kind2dp                     -> point (a pin) / flat (a weld: the flat foot);
//                                          a SEGMENT spec (a second point) switches
//                                          between the two by itself
// - ground_contact2dp                   -> the event layer: step(), heights, reactions,
//                                          the centre of pressure, the unilateral switch;
//                                          its nested spec (which frame point touches,
//                                          of which kind) and event (a touchdown /
//                                          lift-off record) -- nested so the binding
//                                          generator leaves them with the class
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga::pga {

enum class contact_kind2dp { point, flat };

class ground_contact2dp {

  public:

    struct spec {
        size_t frame{0};               // the frame carrying the contact point
        vec2dp point_b{0.0, 0.0, 1.0}; // the point in that frame's coordinates
        contact_kind2dp kind{contact_kind2dp::point};
        value_t half_length{0.0}; // flat: the foot's half length (the centre of
                                  // pressure may not leave +-half_length -- tipping)
        bool unilateral{true};    // false: a BILATERAL pin that may pull (a clamped
                                  // foot, or a held-torso statics whose pins pull by
                                  // construction) -- it lands and can be released by
                                  // hand, but the lift-off rule leaves it alone
        // a SEGMENT: the rigid segment from point_b to end_b (the same frame) standing on
        // the ground -- the layer switches it between a pin at either end and a weld at
        // point_b (see the header); `kind` and `half_length` are not used by it
        bool segment{false};
        vec2dp end_b{0.0, 0.0, 1.0};
    };

    struct event {
        value_t t{};                  // simulation time of the event
        size_t contact{};             // which contact
        bool touchdown{};             // true: touchdown (with impact), false: lift-off
        vec2dp point{};               // the world point it happened at
        std::vector<value_t> impulse; // touchdown: the impact map's Lambda (all rows)
        size_t end{0}; // a segment: the end it happened at (0 = point_b, 1 = end_b)
    };

    // the ground line and the tree frame the contacts are anchored to (the world root)
    ground_contact2dp(closed_loop_system2dp& cl, bivec2dp const& ground,
                      size_t ground_frame = 0) :
        cl_{&cl}, L_{unitize(ground)}, ground_frame_{ground_frame}
    {
        // the upward normal is the gradient of the height, which is linear in P
        n_ = vec2dp{height_of(vec2dp{1.0, 0.0, 0.0}), height_of(vec2dp{0.0, 1.0, 0.0}),
                    0.0};
        value_t const nn = std::sqrt(n_.x * n_.x + n_.y * n_.y);
        n_ = vec2dp{n_.x / nn, n_.y / nn, 0.0};
    }

    // A COPY OF THE LAYER, ATTACHED TO A COPY OF ITS SYSTEM -- a forecast run on a copied
    // state. The layer acts through a pointer to its system, so a plain copy would keep
    // stepping the ORIGINAL; that is why the plain copy is not public. `cl` must be a
    // copy of the system this layer was built on: the contacts' loop constraints and
    // ground frames live there, by index (a system with a different loop count is
    // refused). Functions registered on the system are copied as they are, so one that
    // captured the original system by reference still reads the original.
    ground_contact2dp(ground_contact2dp const& other, closed_loop_system2dp& cl) :
        ground_contact2dp(other)
    {
        if (cl.loop_count() != other.cl_->loop_count())
            throw std::invalid_argument(
                "ground_contact2dp: the system is not a copy of the layer's system");
        cl_ = &cl;
    }
    ground_contact2dp& operator=(ground_contact2dp const&) = delete;

    // the ground line through two points; ABOVE the ground is to the left of P0 -> P1
    // (a floor: P0 = (0, y0), P1 = (1, y0))
    static bivec2dp ground_line(vec2dp const& P0, vec2dp const& P1)
    {
        return unitize(wdg(P0, P1));
    }

    // signed height of a world point above the ground: the incidence of the point
    // with the line (zero on it, positive above)
    value_t height(vec2dp const& P_world) const { return height_of(unitize(P_world)); }

    vec2dp const& normal() const { return n_; } // the upward unit normal
    bivec2dp const& ground() const { return L_; }

    // register a contact (its loop constraints are added to the system, inactive);
    // returns the contact index
    size_t add(spec const& spec)
    {
        contact c;
        c.spec = spec;
        vec2dp const P = point_world(spec.frame, spec.point_b);
        c.pin = cl_->add_loop_constraint(loop_constraint2dp{
            spec.frame, spec.point_b, ground_frame_, P, constraint2dp::coincidence});
        cl_->set_loop_active(c.pin, false);
        if (spec.segment) c.spec.kind = contact_kind2dp::point; // its state, not a choice
        if (spec.kind == contact_kind2dp::flat || spec.segment) {
            std::string const name = "ground_contact_" + std::to_string(contacts_.size());
            cl_->add_frame(static_frame2dp(name, P, 0.0), ground_frame_);
            c.ground = cl_->index_of(name);
            c.weld = cl_->add_loop_constraint(loop_constraint2dp{
                spec.frame, spec.point_b, c.ground, O_2dp, constraint2dp::frame});
            cl_->set_loop_active(c.weld, false);
        }
        contacts_.push_back(c);
        return contacts_.size() - 1;
    }

    // Advance the system by dt with the contact events resolved in order: each
    // touchdown is located by bisection on the sub-interval it happens in, applied,
    // and the remainder of dt is integrated on -- so a landing is exact in time rather
    // than late by up to one step. Contacts reaching the ground within simultaneity() of
    // the first are engaged TOGETHER, in one impact map (see touchdown). Lift-offs are
    // checked at the end of the step.
    void step(value_t dt, value_t tol = value_t(1e-10))
    {
        value_t remaining = dt;
        int guard = 0;
        while (remaining > 0.0 && guard++ < 16) {
            closed_loop_system2dp const before = *cl_;
            // the points that may touch down: an inactive contact's point; a segment's
            // two ends while it is off the ground, its free end while one end is pinned
            std::vector<landing> cand;
            for (size_t c = 0; c < contacts_.size(); ++c) {
                contact const& k = contacts_[c];
                if (!k.active) {
                    cand.push_back({c, 0});
                    if (k.spec.segment) cand.push_back({c, 1});
                }
                else if (k.spec.segment && k.pinned_end >= 0)
                    cand.push_back({c, size_t(1 - k.pinned_end)});
            }
            std::vector<value_t> h0(cand.size());
            for (size_t i = 0; i < cand.size(); ++i)
                h0[i] = height(end_point(cand[i].idx, cand[i].end));
            cl_->step(remaining);
            // the first point that crossed the ground downward in this interval
            size_t hit = cand.size();
            value_t tau = remaining;
            std::vector<value_t> t_cross(cand.size(), value_t(-1.0));
            for (size_t i = 0; i < cand.size(); ++i) {
                // (h0 >= -1e-9: a foot released at the ground line sits at a height
                // of round-off size, possibly negative; it must still be able to land)
                if (h0[i] < -1e-9 || height(end_point(cand[i].idx, cand[i].end)) >= 0.0)
                    continue;
                // bisection on tau in (0, remaining]: not crossed at 0, crossed at
                // remaining; a copy of the pre-event state is stepped by tau
                value_t lo = 0.0, hi = remaining;
                for (int j = 0; j < 60 && hi - lo > tol; ++j) {
                    value_t const mid = 0.5 * (lo + hi);
                    closed_loop_system2dp trial = before;
                    trial.step(mid);
                    if (height(trial_point(trial, cand[i].idx, cand[i].end)) < 0.0)
                        hi = mid;
                    else lo = mid;
                }
                t_cross[i] = hi;
                if (hi < tau) {
                    tau = hi;
                    hit = i;
                }
            }
            if (hit == cand.size()) break; // no event: the whole interval is done
            // every point crossing within the simultaneity window of the first lands
            // with it (a negative window: the first one alone)
            std::vector<landing> lands{cand[hit]};
            for (size_t i = 0; i < cand.size(); ++i)
                if (i != hit && t_cross[i] >= 0.0 && t_cross[i] <= tau + simultaneity_)
                    lands.push_back(cand[i]);
            *cl_ = before;
            cl_->step(tau);
            touchdown(lands);
            remaining -= tau;
        }
        // lift-off: an active contact whose normal reaction pulls is released; a flat
        // segment whose centre of pressure lies beyond an end rolls onto that end
        //
        // THE ROLL IS INDEPENDENT OF `unilateral`, and the release below is not. That is
        // deliberate: `unilateral` says whether the ground may PULL at this contact, i.e.
        // whether the lift-off rule applies to it, while a segment's pin-or-weld is the
        // segment's own kinematic state -- the thing that makes it a segment rather than
        // two pins. A caller that switches a segment to bilateral for a while (to keep a
        // contact from vanishing while its load goes to zero) therefore still gets the
        // roll when the centre of pressure leaves an end, which is what keeps the
        // constraint set consistent with where the load actually is. If a segment must
        // stay welded flat regardless, do not model it as a contact: weld it with a
        // `frame` loop closure, which has no rules of its own.
        read_reactions();
        bool rolled = false;
        for (size_t c = 0; c < contacts_.size(); ++c) {
            // (whatever the sign of the TOTAL normal force: a segment lifted at one end
            // and pressed at the other pulls in total and must still roll -- measured, a
            // sole at N -4.5 N, one end -9.7 N, stayed welded when this required N > 0)
            if (contacts_[c].active && contacts_[c].spec.segment &&
                contacts_[c].spec.kind == contact_kind2dp::flat && roll_if_beyond(c)) {
                rolled = true;
                continue; // its reaction is re-read below before anything else acts
            }
            if (contacts_[c].active && contacts_[c].spec.unilateral &&
                contacts_[c].normal_force <= -release_threshold_ && would_separate(c))
                release(c);
        }
        if (rolled) read_reactions();
    }

    // switch an ACTIVE contact's kind at run time (a pin becoming a flat foot when the
    // other leg lifts, and back): flat engages the weld with the impact map (the
    // foot's angular velocity is absorbed), point drops the angular row
    void set_kind(size_t idx, contact_kind2dp kind)
    {
        contact& c = contacts_[idx];
        if (c.spec.segment)
            throw std::invalid_argument("ground_contact2dp::set_kind: a segment's kind "
                                        "follows its contact state");
        if (c.weld == npos && kind == contact_kind2dp::flat) {
            // a point contact promoted to flat: create its weld on demand
            std::string const name = "ground_contact_" + std::to_string(idx);
            cl_->add_frame(static_frame2dp(name, contact_point(idx), 0.0), ground_frame_);
            c.ground = cl_->index_of(name);
            c.weld = cl_->add_loop_constraint(loop_constraint2dp{
                c.spec.frame, c.spec.point_b, c.ground, O_2dp, constraint2dp::frame});
            cl_->set_loop_active(c.weld, false);
        }
        c.spec.kind = kind;
        if (!c.active) return;
        if (kind == contact_kind2dp::flat) {
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
        if (c.spec.segment) {
            // both ends on the ground: flat, welded at point_b; otherwise a pin at the
            // lower end
            value_t const h_0 = height(end_point(idx, 0)),
                          h_1 = height(end_point(idx, 1));
            if (std::abs(h_0) < 1e-6 && std::abs(h_1) < 1e-6) {
                vec2dp const at = unitize(project_onto(end_point(idx, 0), L_));
                pose_ground_frame(idx, at);
                cl_->set_loop_active(c.weld, true);
                c.spec.kind = contact_kind2dp::flat;
                c.pinned_end = -1;
                c.active = true;
                events_.push_back({cl_->system().time(), idx, true, at, {}, 0});
            }
            else {
                size_t const e = (h_1 < h_0) ? 1 : 0;
                vec2dp const at = unitize(project_onto(end_point(idx, e), L_));
                cl_->set_loop_anchors(c.pin, e == 1 ? c.spec.end_b : c.spec.point_b, at);
                cl_->set_loop_active(c.pin, true);
                c.pinned_end = int(e);
                c.active = true;
                events_.push_back({cl_->system().time(), idx, true, at, {}, e});
            }
            read_reactions();
            return;
        }
        vec2dp const at = unitize(project_onto(contact_point(idx), L_));
        if (c.spec.kind == contact_kind2dp::flat) {
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
        c.moment = 0.0;
        size_t const end = (c.spec.segment && c.pinned_end == 1) ? 1 : 0;
        events_.push_back(
            {cl_->system().time(), idx, false, contact_point(idx), {}, end});
        if (c.spec.segment) {
            c.pinned_end = -1;
            c.spec.kind = contact_kind2dp::point;
        }
    }

    // --- queries ---------------------------------------------------------------

    size_t count() const { return contacts_.size(); }
    bool active(size_t idx) const { return contacts_[idx].active; }
    contact_kind2dp kind(size_t idx) const { return contacts_[idx].spec.kind; }
    // the contact's point in the world: point_b -- for a segment the pinned end (point_b
    // while flat or off the ground), where its reaction acts
    vec2dp contact_point(size_t idx) const
    {
        auto const& c = contacts_[idx];
        return point_world(c.spec.frame, (c.spec.segment && c.pinned_end == 1)
                                             ? c.spec.end_b
                                             : c.spec.point_b);
    }
    // a segment's end in the world (0 = point_b, 1 = end_b; any contact's end 0 is its
    // point), and the end it is pinned at (-1: off the ground, or flat)
    vec2dp end_point(size_t idx, size_t end) const
    {
        auto const& s = contacts_[idx].spec;
        return point_world(s.frame, (s.segment && end == 1) ? s.end_b : s.point_b);
    }
    int pinned_end(size_t idx) const { return contacts_[idx].pinned_end; }
    value_t contact_height(size_t idx) const { return height(contact_point(idx)); }

    // SIMULTANEOUS TOUCHDOWNS. Contacts that reach the ground within this time of the
    // first one in a step are engaged together and resolved by ONE impact map. Resolved
    // one after the other, simultaneous impacts come out order-dependent -- a level rod
    // falling flat puts the pinned-end impulse, about a quarter of m v, on whichever end
    // is taken first, where the combined map gives each end half. Default 1 microsecond
    // (simultaneous to well within any control step); a negative window resolves every
    // contact on its own, in the order found.
    void set_simultaneity(value_t window) { simultaneity_ = window; }
    value_t simultaneity() const { return simultaneity_; }

    // the reaction the ground exerts on the contact point (world), from the last
    // multipliers read: the force, its normal component, and (flat) the moment
    vec2dp force(size_t idx) const { return contacts_[idx].force; }
    value_t normal_force(size_t idx) const { return contacts_[idx].normal_force; }
    value_t moment(size_t idx) const { return contacts_[idx].moment; }

    // flat: the centre of pressure along the foot, moment / normal force (0 when the
    // foot is unloaded), and whether it has left the foot
    value_t cop(size_t idx) const
    {
        auto const& c = contacts_[idx];
        return (c.normal_force > 0.0) ? c.moment / c.normal_force : 0.0;
    }
    bool tipping(size_t idx) const
    {
        auto const& c = contacts_[idx];
        // (a segment does not tip: it rolls onto its end by itself)
        return c.active && !c.spec.segment && c.spec.kind == contact_kind2dp::flat &&
               std::abs(cop(idx)) > c.spec.half_length;
    }

    // The net wrench the ground exerts, summed over the ACTIVE contacts: each
    // contact's force line wdg(P_i, f_i), plus -- for a flat foot -- its ankle
    // couple. A couple is the weight-zero bivector m*e12: its moment is m about
    // EVERY point, which is what makes it addable to a force line without a
    // reference point. The sum is ONE bivector carrying both resultants, the force
    // as its weight (att) and the moment as its bulk, so a caller needs no separate
    // force and torque accumulator. A null bivector while nothing touches.
    bivec2dp reaction_wrench() const
    {
        bivec2dp W{0.0, 0.0, 0.0};
        for (size_t i = 0; i < contacts_.size(); ++i) {
            if (!contacts_[i].active) continue;
            W += wdg(contact_point(i), contacts_[i].force) +
                 bivec2dp{0.0, 0.0, contacts_[i].moment};
        }
        return W;
    }

    // The zero-moment point: the point of the ground about which the net ground
    // reaction exerts no moment. IN THE PLANE a wrench simply IS its line of action
    // -- the moment about a point has the single component e12, so killing it is the
    // entire condition -- and the ZMP is therefore the MEET of that line with the
    // ground:
    //
    //     ZMP = rwdg(reaction_wrench(), ground())
    //
    // No moment balance is solved and no coordinate frame is chosen; it is exact for
    // an arbitrarily slanted resultant and on a SLOPE, since ground() is whatever
    // line was built. Tedrake's centre-of-pressure formula sum(p_i N_i)/sum(N_i) is
    // the flat-ground, vertical-force special case of this meet.
    //
    // The 3D header carries the SAME expression, though it earns it differently:
    // there the definition weakens to "no TANGENTIAL moment", and the meet tracks
    // that by annihilating the couple along the ground normal -- see its zmp().
    //
    // Returned unitized (w == 1) whenever the ground carries load. With no net
    // normal load the reaction line is parallel to the ground and the honest answer
    // is the point at infinity where they meet, so the raw meet is returned with
    // w == 0 -- test .z before using .x/.y.
    vec2dp zmp() const
    {
        vec2dp const X = rwdg(reaction_wrench(), L_);
        return (std::abs(X.z) > eps) ? unitize(X) : X;
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
        read_reactions(lam);
    }

    // the same, from the multipliers of a solve the caller already holds (the
    // lambda of joint_accelerations / sync_accelerations at the current state) --
    // no second solve
    void read_reactions(std::vector<value_t> const& lam)
    {
        size_t r = 0;
        for (size_t k = 0; k < cl_->loop_count(); ++k) {
            if (!cl_->loop_active(k)) continue;
            bool const weld = cl_->loop(k).type == constraint2dp::frame;
            size_t const rows = weld ? 3 : 2;
            for (auto& c : contacts_) {
                if (k != c.pin && k != c.weld) continue;
                c.force = vec2dp{-lam[r], -lam[r + 1], 0.0}; // -lambda: on the foot
                c.normal_force = c.force.x * n_.x + c.force.y * n_.y;
                c.moment = weld ? -lam[r + 2] : 0.0;
            }
            r += rows;
        }
    }

  private:

    // the plain copy still points at the original system: reached only through the
    // re-attaching copy above
    ground_contact2dp(ground_contact2dp const&) = default;

    static constexpr size_t npos = size_t(-1);

    struct contact {
        spec spec;
        size_t pin{npos};    // the coincidence loop
        size_t weld{npos};   // the frame loop (flat), npos if none yet
        size_t ground{npos}; // the fixed ground frame of the weld
        bool active{false};
        vec2dp force{0.0, 0.0, 0.0};
        value_t normal_force{0.0};
        value_t moment{0.0};
        int pinned_end{-1}; // a segment: the end its pin holds (-1: off, or flat)
    };
    // a point that may touch down: the contact, and (a segment) which end
    struct landing {
        size_t idx{0};
        size_t end{0};
    };

    // the separation test of the lift-off rule: does the contact point move away from
    // the ground once its constraint is removed? (a short trial step on a copy)
    bool would_separate(size_t idx) const
    {
        contact const& c = contacts_[idx];
        closed_loop_system2dp trial = *cl_;
        trial.set_loop_active(c.pin, false);
        if (c.weld != npos) trial.set_loop_active(c.weld, false);
        // the point that must rise: the contact's; a pinned segment's pinned end; a flat
        // segment's LOWER end (lifting one end only is a roll, not a lift-off)
        auto low = [&](closed_loop_system2dp& s) {
            if (!c.spec.segment) return height(trial_point(s, idx, 0));
            if (c.pinned_end >= 0)
                return height(trial_point(s, idx, size_t(c.pinned_end)));
            return std::min(height(trial_point(s, idx, 0)),
                            height(trial_point(s, idx, 1)));
        };
        value_t const h0 = low(trial);
        trial.step(separation_dt_);
        return low(trial) > h0 + value_t(1e-9);
    }

    // A flat segment rolls onto one end when the other end PULLS. The weld's reaction --
    // normal force N and moment m about point_b -- is the same as two normal forces at
    // the ends,
    //
    //     F_B = m / s_B,    F_A = N - F_B,
    //
    // s_B the segment's extent along the ground tangent t = (n.y, -n.x) (the direction in
    // which a force line through point_b + s t carries the moment s N about point_b, so
    // cop() = m / N is that s). An end pulling by more than release_threshold lets go
    // while the other end still pushes: the weld goes, the pin holds the pushing end
    // where it is (a subset of the weld's rows, which the velocities already satisfy --
    // no impact), and the lifted end's lift-off is recorded. The total N may be a pull
    // (lifted at one end, pressed at the other). Both ends pulling is not a roll: the
    // whole segment leaves by the lift-off rule. Returns whether it rolled.
    bool roll_if_beyond(size_t idx)
    {
        contact& c = contacts_[idx];
        vec2dp const A = end_point(idx, 0), B = end_point(idx, 1);
        value_t const sB = (B.x - A.x) * n_.y - (B.y - A.y) * n_.x;
        if (std::abs(sB) < 1e-9) return false; // a segment standing on end has no lever
        value_t const F_B = c.moment / sB, F_A = c.normal_force - F_B;
        int e = -1;
        if (F_A <= -release_threshold_ && F_B > -release_threshold_) e = 1;      // onto B
        else if (F_B <= -release_threshold_ && F_A > -release_threshold_) e = 0; // onto A
        if (e < 0) return false;
        size_t const pe = size_t(e), other = 1 - pe;
        cl_->set_loop_active(c.weld, false);
        cl_->set_loop_anchors(c.pin, pe == 1 ? c.spec.end_b : c.spec.point_b,
                              unitize(project_onto(end_point(idx, pe), L_)));
        cl_->set_loop_active(c.pin, true);
        c.spec.kind = contact_kind2dp::point;
        c.pinned_end = e;
        c.moment = 0.0;
        events_.push_back(
            {cl_->system().time(), idx, false, end_point(idx, other), {}, other});
        return true;
    }

    value_t height_of(vec2dp const& P) const { return -value_t(wdg(P, L_)); }

    vec2dp point_world(size_t frame, vec2dp const& point_b) const
    {
        return unitize(move2dp(point_b, cl_->system().get_pos_trafo(frame, 0)));
    }
    vec2dp trial_point(closed_loop_system2dp& trial, size_t idx, size_t end = 0) const
    {
        auto const& s = contacts_[idx].spec;
        vec2dp const& p_b = (s.segment && end == 1) ? s.end_b : s.point_b;
        return unitize(move2dp(p_b, trial.system().get_pos_trafo(s.frame, 0)));
    }

    // the weld's ground frame: at the landing point, with the contact frame's heading
    // (so the weld engages without a jump)
    void pose_ground_frame(size_t idx, vec2dp const& at)
    {
        contact const& c = contacts_[idx];
        auto M = cl_->system().get_pos_trafo(c.spec.frame, 0);
        if (value_t(gr3(M)) < 0.0) M = -M; // the near sheet of the double cover
        cl_->system().set_pose(c.ground, at, 2.0 * rlog(M).z);
    }

    // touchdown of the contacts idxs at the current state (their points on or just below
    // the ground; one landing within the simultaneity window may still sit above it by
    // that window's travel, which the next step's position projection closes): each is
    // anchored at its projection onto the ground line and switched on, and ONE impact map
    // resolves them together. Each event records that shared Lambda (all active rows).
    // A segment's end landing: a pin at that end while the segment is off the ground and
    // its other end is not landing with it; the WELD (at point_b) when the other end is
    // already pinned or lands in the same window
    void touchdown(std::vector<landing> const& lands)
    {
        std::vector<vec2dp> at;
        at.reserve(lands.size());
        for (size_t k = 0; k < lands.size(); ++k) {
            size_t const idx = lands[k].idx, end = lands[k].end;
            contact& c = contacts_[idx];
            vec2dp const landing = unitize(project_onto(end_point(idx, end), L_));
            at.push_back(landing);
            if (!c.spec.segment) {
                if (c.spec.kind == contact_kind2dp::flat) {
                    pose_ground_frame(idx, landing);
                    cl_->set_loop_active(c.weld, true);
                }
                else {
                    cl_->set_loop_anchors(c.pin, c.spec.point_b, landing);
                    cl_->set_loop_active(c.pin, true);
                }
                c.active = true;
                continue;
            }
            bool both_now = false;
            for (auto const& l : lands)
                both_now = both_now || (l.idx == idx && l.end != end);
            bool const other_pinned = c.active && c.pinned_end == int(1 - end);
            if (both_now || other_pinned) {
                if (c.spec.kind != contact_kind2dp::flat) { // weld once
                    cl_->set_loop_active(c.pin, false);
                    pose_ground_frame(idx, unitize(project_onto(end_point(idx, 0), L_)));
                    cl_->set_loop_active(c.weld, true);
                    c.spec.kind = contact_kind2dp::flat;
                    c.pinned_end = -1;
                }
            }
            else {
                cl_->set_loop_anchors(c.pin, end == 1 ? c.spec.end_b : c.spec.point_b,
                                      landing);
                cl_->set_loop_active(c.pin, true);
                c.pinned_end = int(end);
            }
            c.active = true;
        }
        std::vector<value_t> const Lam = cl_->impact(0.0);
        for (size_t k = 0; k < lands.size(); ++k)
            events_.push_back(
                {cl_->system().time(), lands[k].idx, true, at[k], Lam, lands[k].end});
    }

    closed_loop_system2dp* cl_;
    bivec2dp L_;
    vec2dp n_{0.0, 1.0, 0.0};
    size_t ground_frame_;
    std::vector<contact> contacts_;
    std::vector<event> events_;
    value_t release_threshold_{value_t(0.5)};
    value_t separation_dt_{value_t(1e-4)}; // the trial step of would_separate
    value_t simultaneity_{value_t(1e-6)};  // contacts landing together (set_simultaneity)
};

} // namespace hd::ga::pga
