#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_fmt_core.hpp"

// hint: inclusion sequence is important
// (scalar and vector types need to come before multivector types for mv ctors to work)

#include "../type_t/ga_scalar_t.hpp"

#include "../type_t/ga_vec4_t.hpp"

#include "../type_t/ga_bvec6_t.hpp"

#include "../type_t/ga_mvec16_t.hpp"
#include "../type_t/ga_mvec2_t.hpp"
#include "../type_t/ga_mvec8_t.hpp"

// specific convenience types for pga
#include "../type_t/ga_type4ds.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// Formatting support for Space-Time Geometric Algebra (STA) types
// Includes formatters for scalars, vectors, bi-, tri-, quadvectors (ps) and multivectors
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag> for PGA types: Scalar4ds<T>, PScalar4ds<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::scalar4ds_tag> ||
             std::is_same_v<Tag, hd::ga::pscalar4ds_tag>)
struct fmt::formatter<hd::ga::Scalar_t<T, Tag>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::Scalar_t<T, Tag>& v, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                     hd::ga::Scalar_t<T, hd::ga::scalar4ds_tag>>) {
            auto out = fmt::format_to(ctx.out(), "Scalar4ds(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(T(v))));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar4ds_tag>>) {
            auto out = fmt::format_to(ctx.out(), "PScalar4ds(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(T(v))));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(T(v))));
            return fmt::format_to(out, ")");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec4_t<T, Tag> includes Vec4ds<T>, TriVec4ds<T> (STA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::vec4ds_tag> ||
             std::is_same_v<Tag, hd::ga::trivec4ds_tag>)
struct fmt::formatter<hd::ga::Vec4_t<T, Tag>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::Vec4_t<T, Tag>& v, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
                                     hd::ga::Vec4_t<T, hd::ga::vec4ds_tag>>) {
            auto out = fmt::format_to(ctx.out(), "Vec4ds(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.x)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.y)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.z)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.w)));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
                                          hd::ga::Vec4_t<T, hd::ga::trivec4ds_tag>>) {
            auto out = fmt::format_to(ctx.out(), "TriVec4ds(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.x)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.y)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.z)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.w)));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.x)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.y)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.z)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.w)));
            return fmt::format_to(out, ")");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T, Tag> - includes BiVec4ds<T> (STA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::bivec4ds_tag>
struct fmt::formatter<hd::ga::BVec6_t<T, Tag>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::BVec6_t<T, Tag>& v, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        // Format as BiVec4ds(vx, vy, vz, mx, my, mz)
        auto out = fmt::format_to(ctx.out(), "BiVec4ds(");
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.vx)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.vy)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.vz)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.mx)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.my)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.mz)));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T, Tag> includes DualNumber4ds<T> (STA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::dual_number4ds_tag>)
struct fmt::formatter<hd::ga::MVec2_t<T, Tag>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec2_t<T, Tag>& v, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        if constexpr (std::is_same_v<hd::ga::MVec2_t<T, Tag>,
                                     hd::ga::MVec2_t<T, hd::ga::dual_number4ds_tag>>) {
            auto out = fmt::format_to(ctx.out(), "DualNum4ds(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c0)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c1)));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c0)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c1)));
            return fmt::format_to(out, ")");
        }
    }
};


////////////////////////////////////////////////////////////////////////////////
// MVec8_t<T, Tag> includes MVec4ds_E<T> and MVec4ds_U<T> (STA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::mvec4ds_e_tag> ||
             std::is_same_v<Tag, hd::ga::mvec4ds_u_tag>)
struct fmt::formatter<hd::ga::MVec8_t<T, Tag>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec8_t<T, Tag>& v, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                     hd::ga::MVec8_t<T, hd::ga::mvec4ds_e_tag>>) {
            // Format as MVec4ds_E(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec4ds_E(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c0)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c1)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c2)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c3)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c4)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c5)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c6)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c7)));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                          hd::ga::MVec8_t<T, hd::ga::mvec4ds_u_tag>>) {
            // Format as MVec4ds_U(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec4ds_U(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c0)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c1)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c2)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c3)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c4)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c5)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c6)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c7)));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c0)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c1)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c2)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c3)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c4)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c5)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c6)));
            out = fmt::format_to(out, ",");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c7)));
            return fmt::format_to(out, ")");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> includes MVec4ds<T> (STA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::mvec4ds_tag>
struct fmt::formatter<hd::ga::MVec16_t<T, Tag>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec16_t<T, Tag>& v, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        // Format as MVec4ds(c0, c1, c2, ..., c15)
        auto out = fmt::format_to(ctx.out(), "MVec4ds(");
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c0)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c1)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c2)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c3)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c4)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c5)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c6)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c7)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c8)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c9)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c10)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c11)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c12)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c13)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c14)));
        out = fmt::format_to(out, ",");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(v.c15)));
        return fmt::format_to(out, ")");
    }
};
