#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "fmt/format.h"

// Pull in the magnetostatics functionality
#include "../../ga_usr_magnetostatics.hpp"

////////////////////////////////////////////////////////////////////////////////
// Formatters for the magnetostatics data types (ga_usr_magnetostatics.hpp), which
// includes this header itself at its end -- the types must be complete before the
// specializations below.
//
// The format spec is forwarded to the members, so "{:.3f}" etc. work as they do
// for the GA types.
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::line_current> {

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
    auto format(hd::ga::line_current const& x, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "line_current(");
        out = fmt::format_to(out, "p = ");
        out = fmt::format_to(out, fmt::runtime(child), x.p);
        out = fmt::format_to(out, ", I = ");
        out = fmt::format_to(out, fmt::runtime(child), x.I);
        return fmt::format_to(out, ")");
    }
};

template <> struct fmt::formatter<hd::ga::contour_point> {

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
    auto format(hd::ga::contour_point const& x, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "contour_point(");
        out = fmt::format_to(out, "P = ");
        out = fmt::format_to(out, fmt::runtime(child), x.P);
        out = fmt::format_to(out, ", n = ");
        out = fmt::format_to(out, fmt::runtime(child), x.n);
        out = fmt::format_to(out, ", ds = ");
        out = fmt::format_to(out, fmt::runtime(child), x.ds);
        return fmt::format_to(out, ")");
    }
};

template <> struct fmt::formatter<hd::ga::current_segment> {

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
    auto format(hd::ga::current_segment const& x, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "current_segment(");
        out = fmt::format_to(out, "a = ");
        out = fmt::format_to(out, fmt::runtime(child), x.a);
        out = fmt::format_to(out, ", b = ");
        out = fmt::format_to(out, fmt::runtime(child), x.b);
        out = fmt::format_to(out, ", I = ");
        out = fmt::format_to(out, fmt::runtime(child), x.I);
        return fmt::format_to(out, ")");
    }
};

template <> struct fmt::formatter<hd::ga::surface_point> {

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
    auto format(hd::ga::surface_point const& x, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "surface_point(");
        out = fmt::format_to(out, "P = ");
        out = fmt::format_to(out, fmt::runtime(child), x.P);
        out = fmt::format_to(out, ", n = ");
        out = fmt::format_to(out, fmt::runtime(child), x.n);
        out = fmt::format_to(out, ", dA = ");
        out = fmt::format_to(out, fmt::runtime(child), x.dA);
        return fmt::format_to(out, ")");
    }
};
