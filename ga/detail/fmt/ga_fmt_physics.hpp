#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_fmt_core.hpp"

// Include physics headers to get the HD_GA_PGA*_HAS_PHYSICS_OPS macros
// These headers are already included before ga_fmt_support.hpp in ga_pga.hpp
#include "../../ga_pga2dp_ops_physics.hpp"
#include "../../ga_pga3dp_ops_physics.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// Formatting support for PGA Physics types (Inertia matrices)
// These formatters are conditionally compiled based on availability of hd_solver.hpp
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Inertia2dp<T> - 3x3 inertia matrix for 2D projective GA
////////////////////////////////////////////////////////////////////////////////
#if HD_GA_PGA2DP_HAS_PHYSICS_OPS

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

#endif // HD_GA_PGA2DP_HAS_PHYSICS_OPS


////////////////////////////////////////////////////////////////////////////////
// Inertia3dp<T> - 6x6 inertia matrix for 3D projective GA
////////////////////////////////////////////////////////////////////////////////
#if HD_GA_PGA3DP_HAS_PHYSICS_OPS

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

#endif // HD_GA_PGA3DP_HAS_PHYSICS_OPS
