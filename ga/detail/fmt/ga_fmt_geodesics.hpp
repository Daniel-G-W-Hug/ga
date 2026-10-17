#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "fmt/format.h"

// Pull in the geodesics functionality
#include "../../ga_usr_geodesics.hpp"

////////////////////////////////////////////////////////////////////////////////
// Formatters for the geodetic types (ga_usr_geodesics.hpp), which includes this
// header itself -- the types must be complete before the specializations below.
//
// The format spec is forwarded to the value_t members, so "{:.3f}" etc. work as
// they do for the GA types.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// ellipsoid - reference ellipsoid given by its two radii
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::ellipsoid> {

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
    auto format(hd::ga::ellipsoid const& e, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "ellipsoid(r_equator = ");
        out = fmt::format_to(out, fmt::runtime(child), e.r_equator);
        out = fmt::format_to(out, ", r_pole = ");
        out = fmt::format_to(out, fmt::runtime(child), e.r_pole);
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// geo_pos_dms - position as provided by the user (two strings + elevation)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::geo_pos_dms> {

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
    auto format(hd::ga::geo_pos_dms const& p, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(
            ctx.out(), "geo_pos_dms(lat = {}, lon = {}, height = ", p.lat, p.lon);
        out = fmt::format_to(out, fmt::runtime(child), p.height);
        return fmt::format_to(out, " m)");
    }
};


////////////////////////////////////////////////////////////////////////////////
// geo_pos - position ready to calculate with (angles in rad + total radius)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::geo_pos> {

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
    auto format(hd::ga::geo_pos const& p, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "geo_pos(lat = ");
        out = fmt::format_to(out, fmt::runtime(child), p.lat);
        out = fmt::format_to(out, " rad, lon = ");
        out = fmt::format_to(out, fmt::runtime(child), p.lon);
        out = fmt::format_to(out, " rad, height = ");
        out = fmt::format_to(out, fmt::runtime(child), p.height);
        return fmt::format_to(out, " m)");
    }
};


////////////////////////////////////////////////////////////////////////////////
// enu_frame - the local east/north/up directions, expressed in ECEF
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::enu_frame> {

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
    auto format(hd::ga::pga::enu_frame const& f, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);

        auto out = fmt::format_to(ctx.out(), "enu_frame(east = ");
        out = fmt::format_to(out, fmt::runtime(child), f.east);
        out = fmt::format_to(out, ", north = ");
        out = fmt::format_to(out, fmt::runtime(child), f.north);
        out = fmt::format_to(out, ", up = ");
        out = fmt::format_to(out, fmt::runtime(child), f.up);
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// geo_pos_dms2dp - meridian-section position as the user provides it
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::geo_pos_dms2dp> {

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
    auto format(hd::ga::geo_pos_dms2dp const& p, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);
        auto out = fmt::format_to(ctx.out(), "geo_pos_dms2dp(lat = {}, height = ", p.lat);
        out = fmt::format_to(out, fmt::runtime(child), p.height);
        return fmt::format_to(out, " m)");
    }
};


////////////////////////////////////////////////////////////////////////////////
// geo_pos2dp - meridian-section position, ready to calculate with
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::geo_pos2dp> {

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
    auto format(hd::ga::geo_pos2dp const& p, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);
        auto out = fmt::format_to(ctx.out(), "geo_pos2dp(lat = ");
        out = fmt::format_to(out, fmt::runtime(child), p.lat);
        out = fmt::format_to(out, " rad, height = ");
        out = fmt::format_to(out, fmt::runtime(child), p.height);
        return fmt::format_to(out, " m)");
    }
};


////////////////////////////////////////////////////////////////////////////////
// un_frame - the local (up, north) directions of the meridian section
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::un_frame> {

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
    auto format(hd::ga::pga::un_frame const& p, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);
        auto out = fmt::format_to(ctx.out(), "un_frame(up = ");
        out = fmt::format_to(out, fmt::runtime(child), p.up);
        out = fmt::format_to(out, ", north = ");
        out = fmt::format_to(out, fmt::runtime(child), p.north);
        return fmt::format_to(out, ")");
    }
};


////////////////////////////////////////////////////////////////////////////////
// azimuthal - which projection onto the tangent plane (printed by name)
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::azimuthal> : fmt::formatter<string_view> {

    template <typename FormatContext>
    auto format(hd::ga::pga::azimuthal k, FormatContext& ctx) const
    {
        string_view name = "unknown";
        switch (k) {
            case hd::ga::pga::azimuthal::orthographic:
                name = "orthographic";
                break;
            case hd::ga::pga::azimuthal::gnomonic:
                name = "gnomonic";
                break;
            case hd::ga::pga::azimuthal::stereographic:
                name = "stereographic";
                break;
            case hd::ga::pga::azimuthal::perspective:
                name = "perspective";
                break;
        }
        return fmt::formatter<string_view>::format(name, ctx);
    }
};


////////////////////////////////////////////////////////////////////////////////
// azimuthal_options - the per-kind settings of a projection onto the tangent plane
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::azimuthal_options> {

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
    auto format(hd::ga::pga::azimuthal_options const& o, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);
        auto out = fmt::format_to(ctx.out(), "azimuthal_options(height = ");
        out = fmt::format_to(out, fmt::runtime(child), o.height);
        out = fmt::format_to(out, " m, max_angle = ");
        out = fmt::format_to(out, fmt::runtime(child), o.max_angle);
        return fmt::format_to(out, " rad)");
    }
};


////////////////////////////////////////////////////////////////////////////////
// azimuthal_map - one projection onto the tangent plane at one station
////////////////////////////////////////////////////////////////////////////////

template <> struct fmt::formatter<hd::ga::pga::azimuthal_map> {

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
    auto format(hd::ga::pga::azimuthal_map const& m, FormatContext& ctx) const
    {
        auto const child = fmt::format("{{:{}}}", spec_);
        auto out = fmt::format_to(ctx.out(), "azimuthal_map({}, station = ", m.kind);
        out = fmt::format_to(out, fmt::runtime(child), m.station);
        out = fmt::format_to(out, ", C = ");
        out = fmt::format_to(out, fmt::runtime(child), m.C);
        out = fmt::format_to(out, ", limit = ");
        out = fmt::format_to(out, fmt::runtime(child), m.limit);
        return fmt::format_to(out, ")");
    }
};
