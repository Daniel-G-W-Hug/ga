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
#include "../type_t/ga_type2dc.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// Formatting support for Conformal Geometric Algebra (CGA, 2dc) types
// Includes formatters for scalars, vectors, bi-, tri-, quadvectors (ps) and multivectors
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag> for PGA types: Scalar2dc<T>, PScalar2dc<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::scalar2dc_tag> ||
             std::is_same_v<Tag, hd::ga::pscalar2dc_tag>)
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
                                     hd::ga::Scalar_t<T, hd::ga::scalar2dc_tag>>) {
            auto out = fmt::format_to(ctx.out(), "Scalar2dc(");
            out = fmt::format_to(out, "{}", this->nested(suppress_negative_zero(T(v))));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar2dc_tag>>) {
            auto out = fmt::format_to(ctx.out(), "PScalar2dc(");
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
// Vec4_t<T, Tag> includes Vec2dc<T>, TriVec2dc<T> (CGA 2dc types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::vec2dc_tag> ||
             std::is_same_v<Tag, hd::ga::trivec2dc_tag>)
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
                                     hd::ga::Vec4_t<T, hd::ga::vec2dc_tag>>) {
            auto out = fmt::format_to(ctx.out(), "Vec2dc(");
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
                                          hd::ga::Vec4_t<T, hd::ga::trivec2dc_tag>>) {
            auto out = fmt::format_to(ctx.out(), "TriVec2dc(");
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
// BVec6_t<T, Tag> - includes BiVec2dc<T> (CGA 2dc types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::bivec2dc_tag>
struct fmt::formatter<hd::ga::BVec6_t<T, Tag>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::BVec6_t<T, Tag>& v, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        // Format as BiVec2dc(vx, vy, vz, mx, my, mz)
        auto out = fmt::format_to(ctx.out(), "BiVec2dc(");
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
// MVec2_t<T, Tag> includes DualNumber2dc<T> (CGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::dual_number2dc_tag>)
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
                                     hd::ga::MVec2_t<T, hd::ga::dual_number2dc_tag>>) {
            auto out = fmt::format_to(ctx.out(), "DualNum2dc(");
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
// MVec8_t<T, Tag> includes MVec2dc_E<T> and MVec2dc_U<T> (CGA 2dc types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::mvec2dc_e_tag> ||
             std::is_same_v<Tag, hd::ga::mvec2dc_u_tag>)
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
                                     hd::ga::MVec8_t<T, hd::ga::mvec2dc_e_tag>>) {
            // Format as MVec2dc_E(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec2dc_E(");
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
                                          hd::ga::MVec8_t<T, hd::ga::mvec2dc_u_tag>>) {
            // Format as MVec2dc_U(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec2dc_U(");
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
// MVec16_t<T, Tag> includes MVec2dc<T> (CGA 2dc types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::mvec2dc_tag>
struct fmt::formatter<hd::ga::MVec16_t<T, Tag>> : fmt::nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec16_t<T, Tag>& v, FormatContext& ctx) const
    {
        using hd::ga::detail::suppress_negative_zero;
        // Format as MVec2dc(c0, c1, c2, ..., c15)
        auto out = fmt::format_to(ctx.out(), "MVec2dc(");
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
