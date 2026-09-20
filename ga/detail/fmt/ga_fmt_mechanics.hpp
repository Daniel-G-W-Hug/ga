#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_fmt_core.hpp"

// Pull in the mechanics inertia types so their fmt::formatter specialisations
// below are available wherever ga_pga.hpp is included.
#include "../../pga/ga_pga2dp_ops_mechanics.hpp"
#include "../../pga/ga_pga3dp_ops_mechanics.hpp"

// Pull in the (optional) closed-loop layers so the loop_constraint{2,3}dp formatters
// below are available wherever ga_pga.hpp is included (the constraints headers themselves
// stay free of any fmt dependency).
#include "../../pga/ga_pga2dp_ops_constraints.hpp"
#include "../../pga/ga_pga3dp_ops_constraints.hpp"

// Pull in the unilateral contact layers for the contact_kind{2,3}dp formatters below
// (same reasoning as the closed-loop layers above).
#include "../../pga/ga_pga2dp_ops_contact.hpp"
#include "../../pga/ga_pga3dp_ops_contact.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// Formatting support for PGA mechanics types (Inertia matrices)
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// the scoped enums of the mechanics / constraints / contact layers
//
// One switch per enum, in hd::ga::detail so the binding generator does not scan
// it. No default label: adding an enumerator makes the compiler report the
// unhandled case instead of silently printing the trailing "?".
//
// The formatters below inherit fmt::formatter<fmt::string_view>, so width, fill
// and alignment work on a joint kind as they do on any string, while a numeric
// spec is rejected by fmt's own parser. The composite formatters (joint_state,
// loop_constraint) print their kind member with a plain "{}".
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga::detail {

inline char const* enum_name(hd::ga::pga::joint2dp t)
{
    using hd::ga::pga::joint2dp;
    switch (t) {
        case joint2dp::free:
            return "free";
        case joint2dp::revolute:
            return "revolute";
        case joint2dp::prismatic:
            return "prismatic";
    }
    return "?";
}

inline char const* enum_name(hd::ga::pga::joint3dp t)
{
    using hd::ga::pga::joint3dp;
    switch (t) {
        case joint3dp::free:
            return "free";
        case joint3dp::revolute:
            return "revolute";
        case joint3dp::prismatic:
            return "prismatic";
        case joint3dp::helical:
            return "helical";
        case joint3dp::cylindrical:
            return "cylindrical";
        case joint3dp::spherical:
            return "spherical";
        case joint3dp::planar:
            return "planar";
    }
    return "?";
}

inline char const* enum_name(hd::ga::pga::constraint2dp t)
{
    using hd::ga::pga::constraint2dp;
    switch (t) {
        case constraint2dp::coincidence:
            return "coincidence";
        case constraint2dp::distance:
            return "distance";
        case constraint2dp::frame:
            return "frame";
    }
    return "?";
}

inline char const* enum_name(hd::ga::pga::constraint3dp t)
{
    using hd::ga::pga::constraint3dp;
    switch (t) {
        case constraint3dp::coincidence:
            return "coincidence";
        case constraint3dp::distance:
            return "distance";
        case constraint3dp::frame:
            return "frame";
    }
    return "?";
}

inline char const* enum_name(hd::ga::pga::contact_kind2dp t)
{
    using hd::ga::pga::contact_kind2dp;
    switch (t) {
        case contact_kind2dp::point:
            return "point";
        case contact_kind2dp::flat:
            return "flat";
    }
    return "?";
}

inline char const* enum_name(hd::ga::pga::contact_kind3dp t)
{
    using hd::ga::pga::contact_kind3dp;
    switch (t) {
        case contact_kind3dp::point:
            return "point";
        case contact_kind3dp::flat:
            return "flat";
    }
    return "?";
}

inline char const* enum_name(hd::ga::pga::integrator_kind t)
{
    using hd::ga::pga::integrator_kind;
    switch (t) {
        case integrator_kind::rk4:
            return "rk4";
        case integrator_kind::abm2:
            return "abm2";
    }
    return "?";
}

// shared body of the enum formatters below: print the enumerator's name, with the
// string formatter's own spec handling
template <typename E> struct enum_formatter : fmt::formatter<fmt::string_view> {

    template <typename FormatContext> auto format(E e, FormatContext& ctx) const
    {
        return fmt::formatter<fmt::string_view>::format(enum_name(e), ctx);
    }
};

} // namespace hd::ga::detail

template <>
struct fmt::formatter<hd::ga::pga::joint2dp>
    : hd::ga::detail::enum_formatter<hd::ga::pga::joint2dp> {};

template <>
struct fmt::formatter<hd::ga::pga::joint3dp>
    : hd::ga::detail::enum_formatter<hd::ga::pga::joint3dp> {};

template <>
struct fmt::formatter<hd::ga::pga::constraint2dp>
    : hd::ga::detail::enum_formatter<hd::ga::pga::constraint2dp> {};

template <>
struct fmt::formatter<hd::ga::pga::constraint3dp>
    : hd::ga::detail::enum_formatter<hd::ga::pga::constraint3dp> {};

template <>
struct fmt::formatter<hd::ga::pga::contact_kind2dp>
    : hd::ga::detail::enum_formatter<hd::ga::pga::contact_kind2dp> {};

template <>
struct fmt::formatter<hd::ga::pga::contact_kind3dp>
    : hd::ga::detail::enum_formatter<hd::ga::pga::contact_kind3dp> {};

template <>
struct fmt::formatter<hd::ga::pga::integrator_kind>
    : hd::ga::detail::enum_formatter<hd::ga::pga::integrator_kind> {};


////////////////////////////////////////////////////////////////////////////////
// Inertia2dp<T> - 3x3 inertia matrix for 2D projective GA
////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct fmt::formatter<hd::ga::pga::Inertia2dp<T>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::Inertia2dp<T> const& I, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        auto v = I.view();

        auto out = fmt::format_to(ctx.out(), "Inertia2dp(\n");
        for (size_t row = 0; row < 3; ++row) {
            out = fmt::format_to(out, "  ");
            for (size_t col = 0; col < 3; ++col) {
                out = fmt::format_to(out, "{}",
                                     this->nested(suppress_negative_zero(v[row, col])));
                if (col < 2) {
                    out = fmt::format_to(out, ", ");
                    ctx.advance_to(out);
                }
            }
            if (row < 2) {
                out = fmt::format_to(out, ",\n");
            }
            else {
                out = fmt::format_to(out, "\n");
            }
        }
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// static_frame2dp - a coordinate frame posed relative to its parent
//
// A single nested format spec (e.g. "{:.3f}") flows into *both* the vec2dp origin
// and the value_t angle. fmt::nested_formatter<T> can only nest one type T, so we
// capture the raw spec once in parse() and forward it to each sub-object's own
// formatter (vec2dp's and value_t's), each of which applies it to its scalars.
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::static_frame2dp> {

    // the raw nested spec captured between ':' and '}', e.g. ".3f"
    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it; // fmt expects the iterator left at the closing '}'
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::static_frame2dp const& frame, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;

        // rebuild a child format string "{:<spec>}" and forward it verbatim to
        // vec2dp's and value_t's own formatters
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out =
            fmt::format_to(ctx.out(), "static_frame2dp(name = '{}'", frame.get_name());
        out = fmt::format_to(out, ", origin = ");
        out = fmt::format_to(out, fmt::runtime(child), frame.get_pose().origin);
        out = fmt::format_to(out, ", phi = ");
        out = fmt::format_to(out, fmt::runtime(child),
                             suppress_negative_zero(frame.get_pose().phi));
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// static_system2dp - a tree of static_frame2dp frames
//
// Same spec-forwarding pattern as static_frame2dp, one level up: capture the
// nested spec once and forward it to each static_frame2dp's own formatter,
// which in turn cascades it down to vec2dp / value_t.
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::static_system2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::static_system2dp const& sys, FormatContext& ctx) const
    {
        // forward the captured spec to each static_frame2dp's own formatter
        auto const child = fmt::format("{{:{}}}", spec_); // -> "{:<spec>}"

        auto out = fmt::format_to(ctx.out(), "static_system2dp(");
        if (sys.empty()) return fmt::format_to(out, ")");

        out = fmt::format_to(out, "\n");
        for (size_t idx = 0; idx < sys.size(); ++idx) {
            // show the parent index (root prints its own index, a self-loop)
            out = fmt::format_to(out, "    rf[{}] (parent {}) = ", idx, sys.parent(idx));
            out = fmt::format_to(out, fmt::runtime(child), sys.frame(idx));
            out = fmt::format_to(out, "\n");
        }
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// static_frame3dp - a coordinate frame posed relative to its parent
//
// Same spec-forwarding pattern as static_frame2dp: both sub-objects are vec3dp
// here (origin point and axis*angle rotation), so the captured spec goes to two
// vec3dp formatters instead of a vec2dp and a value_t.
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::static_frame3dp> {

    // the raw nested spec captured between ':' and '}', e.g. ".3f"
    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it; // fmt expects the iterator left at the closing '}'
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::static_frame3dp const& frame, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out =
            fmt::format_to(ctx.out(), "static_frame3dp(name = '{}'", frame.get_name());
        out = fmt::format_to(out, ", origin = ");
        out = fmt::format_to(out, fmt::runtime(child), frame.get_pose().origin);
        out = fmt::format_to(out, ", rot = ");
        out = fmt::format_to(out, fmt::runtime(child), frame.get_pose().rot);
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// static_system3dp - a tree of static_frame3dp frames
//
// Same spec-forwarding pattern as static_system2dp, one level up: capture the
// nested spec once and forward it to each static_frame3dp's own formatter,
// which in turn cascades it down to vec3dp.
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::static_system3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::static_system3dp const& sys, FormatContext& ctx) const
    {
        // forward the captured spec to each static_frame3dp's own formatter
        auto const child = fmt::format("{{:{}}}", spec_); // -> "{:<spec>}"

        auto out = fmt::format_to(ctx.out(), "static_system3dp(");
        if (sys.empty()) return fmt::format_to(out, ")");

        out = fmt::format_to(out, "\n");
        for (size_t idx = 0; idx < sys.size(); ++idx) {
            // show the parent index (root prints its own index, a self-loop)
            out = fmt::format_to(out, "    rf[{}] (parent {}) = ", idx, sys.parent(idx));
            out = fmt::format_to(out, fmt::runtime(child), sys.frame(idx));
            out = fmt::format_to(out, "\n");
        }
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// pose2dp - a frame's pose vs. its parent (origin point + orientation angle)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::pose2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::pose2dp const& p, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "pose2dp(origin = ");
        out = fmt::format_to(out, fmt::runtime(child), p.origin);
        out = fmt::format_to(out, ", phi = ");
        out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(p.phi));
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// kin_state2dp - momentary kinematic state (velocity / acceleration) of a frame
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::kin_state2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::kin_state2dp const& k, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "kin_state2dp(vel = ");
        out = fmt::format_to(out, fmt::runtime(child), k.vel);
        out = fmt::format_to(out, ", acc = ");
        out = fmt::format_to(out, fmt::runtime(child), k.acc);
        out = fmt::format_to(out, ", omega = ");
        out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(k.omega));
        out = fmt::format_to(out, ", alpha = ");
        out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(k.alpha));
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// joint_range2dp / joint_drive2dp - a joint's own specification: where it may
// go, and what its actuator can deliver. Both print their NEUTRAL value (inf
// bounds, an ideal drive) as inf, which is what an unrestricted joint carries.
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::joint_range2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::joint_range2dp const& r, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "joint_range2dp(lo = ");
        out = fmt::format_to(out, fmt::runtime(child), r.lo);
        out = fmt::format_to(out, ", hi = ");
        out = fmt::format_to(out, fmt::runtime(child), r.hi);
        out = fmt::format_to(out, ", k_stop = ");
        out = fmt::format_to(out, fmt::runtime(child), r.k_stop);
        out = fmt::format_to(out, ", c_stop = ");
        out = fmt::format_to(out, fmt::runtime(child), r.c_stop);
        return fmt::format_to(out, ")");
    }
};


template <> struct fmt::formatter<hd::ga::pga::joint_drive2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::joint_drive2dp const& d_, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "joint_drive2dp(tau_max = ");
        out = fmt::format_to(out, fmt::runtime(child), d_.tau_max);
        out = fmt::format_to(out, ", qd_max = ");
        out = fmt::format_to(out, fmt::runtime(child), d_.qd_max);
        out = fmt::format_to(out, ", qdd_max = ");
        out = fmt::format_to(out, fmt::runtime(child), d_.qdd_max);
        out = fmt::format_to(out, ", armature = ");
        out = fmt::format_to(out, fmt::runtime(child), d_.armature);
        return fmt::format_to(out, ", actuated = {})", d_.actuated);
    }
};


////////////////////////////////////////////////////////////////////////////////
// joint2dp / joint_state2dp - the reduced-coordinate joint of a body vs. its
// parent (joint kind + screw generator + rest motor + generalised coord/rate)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::joint_state2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::joint_state2dp const& j, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "joint_state2dp(type = {}", j.type);
        out = fmt::format_to(out, ", screw_b = ");
        out = fmt::format_to(out, fmt::runtime(child), j.screw_b);
        out = fmt::format_to(out, ", rest = ");
        out = fmt::format_to(out, fmt::runtime(child), j.rest);
        out = fmt::format_to(out, ", phi = ");
        out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(j.phi));
        out = fmt::format_to(out, ", omega = ");
        out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(j.omega));
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// kinematic_system2dp - a frame tree plus per-frame momentary kinematic state
//
// Prints each frame's pose (inherited static_frame2dp) and its parent, then the
// stored RELATIVE velocity twist decoded back to (omega, v) -- the configured
// kinematic state of the system.
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::kinematic_system2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::kinematic_system2dp const& sys, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "kinematic_system2dp(");
        if (sys.empty()) return fmt::format_to(out, ")");

        out = fmt::format_to(out, "\n");
        for (size_t idx = 0; idx < sys.size(); ++idx) {
            out = fmt::format_to(out, "    rf[{}] (parent {}) = ", idx, sys.parent(idx));
            out = fmt::format_to(out, fmt::runtime(child), sys.frame(idx));

            // decode the relative twists vec2dp(-vy, vx, w) -> (omega, v) and (alpha, a)
            auto const tw = sys.relative_twist(idx);
            auto const aw = sys.relative_accel_twist(idx);
            out = fmt::format_to(out, "  | rel omega = ");
            out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(tw.z));
            out = fmt::format_to(out, ", v = (");
            out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(tw.y));
            out = fmt::format_to(out, ", ");
            out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(-tw.x));
            out = fmt::format_to(out, "); alpha = ");
            out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(aw.z));
            out = fmt::format_to(out, ", a = (");
            out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(aw.y));
            out = fmt::format_to(out, ", ");
            out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(-aw.x));
            out = fmt::format_to(out, ")\n");
        }
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// kinematic_system3dp - a frame tree plus per-frame momentary kinematic state
//
// Prints each frame's pose (inherited static_frame3dp) and its parent, then the
// stored RELATIVE velocity twist decoded back to (omega, v) -- the configured
// kinematic state of the system. The twist packs the angular velocity into the
// weight slots and the linear one into the bulk (kinematic_system3dp::to_twist),
// so the decode is omega = (vx,vy,vz), v = (mx,my,mz) -- no sign juggling, unlike
// the planar vec2dp(-vy, vx, w).
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::kinematic_system3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::kinematic_system3dp const& sys, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "kinematic_system3dp(");
        if (sys.empty()) return fmt::format_to(out, ")");

        // print one (x, y, z) triple through the captured spec
        auto const triple = [&](auto o, auto x, auto y, auto z) {
            o = fmt::format_to(o, "(");
            o = fmt::format_to(o, fmt::runtime(child), suppress_negative_zero(x));
            o = fmt::format_to(o, ", ");
            o = fmt::format_to(o, fmt::runtime(child), suppress_negative_zero(y));
            o = fmt::format_to(o, ", ");
            o = fmt::format_to(o, fmt::runtime(child), suppress_negative_zero(z));
            return fmt::format_to(o, ")");
        };

        out = fmt::format_to(out, "\n");
        for (size_t idx = 0; idx < sys.size(); ++idx) {
            out = fmt::format_to(out, "    rf[{}] (parent {}) = ", idx, sys.parent(idx));
            out = fmt::format_to(out, fmt::runtime(child), sys.frame(idx));

            auto const tw = sys.relative_twist(idx);
            auto const aw = sys.relative_accel_twist(idx);
            out = fmt::format_to(out, "  | rel omega = ");
            out = triple(out, tw.vx, tw.vy, tw.vz);
            out = fmt::format_to(out, ", v = ");
            out = triple(out, tw.mx, tw.my, tw.mz);
            out = fmt::format_to(out, "; alpha = ");
            out = triple(out, aw.vx, aw.vy, aw.vz);
            out = fmt::format_to(out, ", a = ");
            out = triple(out, aw.mx, aw.my, aw.mz);
            out = fmt::format_to(out, "\n");
        }
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// joint_range3dp / joint_drive3dp - a joint's own specification: where it may
// go, and what its actuator can deliver. Both print their NEUTRAL value (inf
// bounds, an ideal drive) as inf, which is what an unrestricted joint carries.
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::joint_range3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::joint_range3dp const& r, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "joint_range3dp(lo = ");
        out = fmt::format_to(out, fmt::runtime(child), r.lo);
        out = fmt::format_to(out, ", hi = ");
        out = fmt::format_to(out, fmt::runtime(child), r.hi);
        out = fmt::format_to(out, ", k_stop = ");
        out = fmt::format_to(out, fmt::runtime(child), r.k_stop);
        out = fmt::format_to(out, ", c_stop = ");
        out = fmt::format_to(out, fmt::runtime(child), r.c_stop);
        return fmt::format_to(out, ")");
    }
};


template <> struct fmt::formatter<hd::ga::pga::joint_drive3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::joint_drive3dp const& d_, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "joint_drive3dp(tau_max = ");
        out = fmt::format_to(out, fmt::runtime(child), d_.tau_max);
        out = fmt::format_to(out, ", qd_max = ");
        out = fmt::format_to(out, fmt::runtime(child), d_.qd_max);
        out = fmt::format_to(out, ", qdd_max = ");
        out = fmt::format_to(out, fmt::runtime(child), d_.qdd_max);
        out = fmt::format_to(out, ", armature = ");
        out = fmt::format_to(out, fmt::runtime(child), d_.armature);
        return fmt::format_to(out, ", actuated = {})", d_.actuated);
    }
};


////////////////////////////////////////////////////////////////////////////////
// joint3dp / joint_state3dp - the reduced-coordinate joint of a body vs. its
// parent in 3D (joint kind + screw-axis line + rest motor + generalised
// coord/rate)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::joint_state3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::joint_state3dp const& j, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "joint_state3dp(type = {}", j.type);
        out = fmt::format_to(out, ", screw_b = ");
        out = fmt::format_to(out, fmt::runtime(child), j.screw_b);
        out = fmt::format_to(out, ", rest = ");
        out = fmt::format_to(out, fmt::runtime(child), j.rest);
        out = fmt::format_to(out, ", phi = ");
        out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(j.phi));
        out = fmt::format_to(out, ", omega = ");
        out = fmt::format_to(out, fmt::runtime(child), suppress_negative_zero(j.omega));
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// screw_axis2dp - the screw axis of a planar motion (centre, angle, translation)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::screw_axis2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::screw_axis2dp const& s, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);
        auto out = fmt::format_to(ctx.out(), "screw_axis2dp(centre = ");
        out = fmt::format_to(out, fmt::runtime(child), s.centre);
        out = fmt::format_to(out, ", angle = ");
        out = fmt::format_to(out, fmt::runtime(child), s.angle);
        out = fmt::format_to(out, ", translation = ");
        out = fmt::format_to(out, fmt::runtime(child), s.translation);
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// screw_axis3dp - the screw axis of a spatial motion (unit line, angle, distance)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::screw_axis3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::screw_axis3dp const& s, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);
        auto out = fmt::format_to(ctx.out(), "screw_axis3dp(axis = ");
        out = fmt::format_to(out, fmt::runtime(child), s.axis);
        out = fmt::format_to(out, ", angle = ");
        out = fmt::format_to(out, fmt::runtime(child), s.angle);
        out = fmt::format_to(out, ", dist = ");
        out = fmt::format_to(out, fmt::runtime(child), s.dist);
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// screw_system{2,3}dp - span and Lie closure dimensions of a set of joint screws
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::screw_system2dp> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto it = ctx.begin();
        while (it != ctx.end() && *it != '}')
            ++it;
        return it;
    }
    template <typename FormatContext>
    auto format(hd::ga::pga::screw_system2dp const& s, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "screw_system2dp(span = {}, closure = {})",
                              s.span, s.closure);
    }
};

template <> struct fmt::formatter<hd::ga::pga::screw_system3dp> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto it = ctx.begin();
        while (it != ctx.end() && *it != '}')
            ++it;
        return it;
    }
    template <typename FormatContext>
    auto format(hd::ga::pga::screw_system3dp const& s, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "screw_system3dp(span = {}, closure = {})",
                              s.span, s.closure);
    }
};

////////////////////////////////////////////////////////////////////////////////
// loop_constraint2dp - a closed-loop point-coincidence constraint between two
// tree frames (anchor_a in frame_a must coincide with anchor_b in frame_b)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::loop_constraint2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::loop_constraint2dp const& c, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out =
            fmt::format_to(ctx.out(), "loop_constraint2dp(frame_a = {}", c.frame_a);
        out = fmt::format_to(out, ", anchor_a = ");
        out = fmt::format_to(out, fmt::runtime(child), c.anchor_a);
        out = fmt::format_to(out, ", frame_b = {}", c.frame_b);
        out = fmt::format_to(out, ", anchor_b = ");
        out = fmt::format_to(out, fmt::runtime(child), c.anchor_b);
        out = fmt::format_to(out, ", type = {}", c.type);
        out = fmt::format_to(out, ", length = ");
        out = fmt::format_to(out, fmt::runtime(child), c.length);
        return fmt::format_to(out, ", active = {})", c.active);
    }
};


////////////////////////////////////////////////////////////////////////////////
// loop_constraint3dp - a closed-loop point-coincidence constraint between two
// tree frames (anchor_a in frame_a must coincide with anchor_b in frame_b)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::loop_constraint3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::loop_constraint3dp const& c, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out =
            fmt::format_to(ctx.out(), "loop_constraint3dp(frame_a = {}", c.frame_a);
        out = fmt::format_to(out, ", anchor_a = ");
        out = fmt::format_to(out, fmt::runtime(child), c.anchor_a);
        out = fmt::format_to(out, ", frame_b = {}", c.frame_b);
        out = fmt::format_to(out, ", anchor_b = ");
        out = fmt::format_to(out, fmt::runtime(child), c.anchor_b);
        out = fmt::format_to(out, ", type = {}", c.type);
        out = fmt::format_to(out, ", length = ");
        out = fmt::format_to(out, fmt::runtime(child), c.length);
        return fmt::format_to(out, ", active = {})", c.active);
    }
};


////////////////////////////////////////////////////////////////////////////////
// Inertia3dp<T> - 6x6 inertia matrix for 3D projective GA
////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct fmt::formatter<hd::ga::pga::Inertia3dp<T>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::Inertia3dp<T> const& I, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        auto v = I.view();

        auto out = fmt::format_to(ctx.out(), "Inertia3dp(\n");
        for (size_t row = 0; row < 6; ++row) {
            out = fmt::format_to(out, "  ");
            for (size_t col = 0; col < 6; ++col) {
                out = fmt::format_to(out, "{}",
                                     this->nested(suppress_negative_zero(v[row, col])));
                if (col < 5) {
                    out = fmt::format_to(out, ", ");
                    ctx.advance_to(out);
                }
            }
            if (row < 5) {
                out = fmt::format_to(out, ",\n");
            }
            else {
                out = fmt::format_to(out, "\n");
            }
        }
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// body2dp - a body's inertial record (inertia map, its cached inverse, total mass)
template <> struct fmt::formatter<hd::ga::pga::body2dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::body2dp const& b, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "body2dp(mass = ");
        out = fmt::format_to(out, fmt::runtime(child), b.mass);
        out = fmt::format_to(out, ", I = ");
        out = fmt::format_to(out, fmt::runtime(child), b.I);
        out = fmt::format_to(out, ", I_inv = ");
        out = fmt::format_to(out, fmt::runtime(child), b.I_inv);
        return fmt::format_to(out, ")");
    }
};

// body3dp - a body's inertial record (inertia map, its cached inverse, total mass)
template <> struct fmt::formatter<hd::ga::pga::body3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::body3dp const& b, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "body3dp(mass = ");
        out = fmt::format_to(out, fmt::runtime(child), b.mass);
        out = fmt::format_to(out, ", I = ");
        out = fmt::format_to(out, fmt::runtime(child), b.I);
        out = fmt::format_to(out, ", I_inv = ");
        out = fmt::format_to(out, fmt::runtime(child), b.I_inv);
        return fmt::format_to(out, ")");
    }
};

// pose3dp - a frame's pose vs. its parent (origin point + axis*angle rotation)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::pose3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::pose3dp const& p, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "pose3dp(origin = ");
        out = fmt::format_to(out, fmt::runtime(child), p.origin);
        out = fmt::format_to(out, ", rot = ");
        out = fmt::format_to(out, fmt::runtime(child), p.rot);
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// kin_state3dp - momentary kinematic state (velocity / acceleration) of a frame
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::kin_state3dp> {

    fmt::string_view spec_{};

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        auto const begin = ctx.begin();
        auto it = begin;
        while (it != ctx.end() && *it != '}')
            ++it;
        spec_ = fmt::string_view(begin, static_cast<size_t>(it - begin));
        return it;
    }

    template <typename FormatContext>
    auto format(hd::ga::pga::kin_state3dp const& k, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "kin_state3dp(vel = ");
        out = fmt::format_to(out, fmt::runtime(child), k.vel);
        out = fmt::format_to(out, ", acc = ");
        out = fmt::format_to(out, fmt::runtime(child), k.acc);
        out = fmt::format_to(out, ", omega = ");
        out = fmt::format_to(out, fmt::runtime(child), k.omega);
        out = fmt::format_to(out, ", alpha = ");
        out = fmt::format_to(out, fmt::runtime(child), k.alpha);
        return fmt::format_to(out, ")");
    }
};
