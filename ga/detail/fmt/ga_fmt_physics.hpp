#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_fmt_core.hpp"

// Pull in the physics inertia types so their fmt::formatter specialisations
// below are available wherever ga_pga.hpp is included.
#include "../../ga_pga2dp_ops_physics.hpp"
#include "../../ga_pga3dp_ops_physics.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// Formatting support for PGA Physics types (Inertia matrices)
/////////////////////////////////////////////////////////////////////////////////////////


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
