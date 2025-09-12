#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_fmt_core.hpp"

// basic types used by all algebras
#include "../type_t/ga_bvec6_t.hpp"
#include "../type_t/ga_mvec16_t.hpp"
#include "../type_t/ga_mvec2_t.hpp"
#include "../type_t/ga_mvec4_t.hpp"
#include "../type_t/ga_mvec8_t.hpp"
#include "../type_t/ga_scalar_t.hpp"
#include "../type_t/ga_vec3_t.hpp"
#include "../type_t/ga_vec4_t.hpp"

// specific convenience types for pga
#include "../type_t/ga_type2dp.hpp"
#include "../type_t/ga_type3dp.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// Formatting support for Projective Geometric Algebra (PGA) types
// Includes formatters for PGA scalars, vectors, multivectors, and geometric objects
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag> for PGA types: Scalar2dp<T>, PScalar2dp<T>, Scalar3dp<T>,
// PScalar3dp<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::scalar2dp_tag> ||
             std::is_same_v<Tag, hd::ga::pscalar2dp_tag> ||
             std::is_same_v<Tag, hd::ga::scalar3dp_tag> ||
             std::is_same_v<Tag, hd::ga::pscalar3dp_tag>)
struct fmt::formatter<hd::ga::Scalar_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::Scalar_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                     hd::ga::Scalar_t<T, hd::ga::scalar2dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "Scalar2dp(");
            out = fmt::format_to(out, "{}", this->nested(T(v)));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar2dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "PScalar2dp(");
            out = fmt::format_to(out, "{}", this->nested(T(v)));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::scalar3dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "Scalar3dp(");
            out = fmt::format_to(out, "{}", this->nested(T(v)));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar3dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "PScalar3dp(");
            out = fmt::format_to(out, "{}", this->nested(T(v)));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(T(v)));
            return fmt::format_to(out, ")");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec3_t<T, Tag> includes Vec2dp<T>, BiVec2dp<T> (PGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::vec2dp_tag> ||
             std::is_same_v<Tag, hd::ga::bivec2dp_tag>)
struct fmt::formatter<hd::ga::Vec3_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::Vec3_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::Vec3_t<T, Tag>,
                                     hd::ga::Vec3_t<T, hd::ga::vec2dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "Vec2dp(");
            out = fmt::format_to(out, "{}", this->nested(v.x));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.y));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.z));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::Vec3_t<T, Tag>,
                                          hd::ga::Vec3_t<T, hd::ga::bivec2dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "BiVec2dp(");
            out = fmt::format_to(out, "{}", this->nested(v.x));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.y));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.z));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(v.x));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.y));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.z));
            return fmt::format_to(out, ")");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec4_t<T, Tag> includes Vec3dp<T>, TriVec3dp<T> (PGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::vec3dp_tag> ||
             std::is_same_v<Tag, hd::ga::trivec3dp_tag>)
struct fmt::formatter<hd::ga::Vec4_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::Vec4_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
                                     hd::ga::Vec4_t<T, hd::ga::vec3dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "Vec3dp(");
            out = fmt::format_to(out, "{}", this->nested(v.x));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.y));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.z));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.w));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
                                          hd::ga::Vec4_t<T, hd::ga::trivec3dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "TriVec3dp(");
            out = fmt::format_to(out, "{}", this->nested(v.x));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.y));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.z));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.w));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(v.x));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.y));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.z));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.w));
            return fmt::format_to(out, ")");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T, Tag> - includes BiVec3dp<T> (PGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::bivec3dp_tag>
struct fmt::formatter<hd::ga::BVec6_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::BVec6_t<T, Tag>& v, FormatContext& ctx) const
    {
        // Format as BiVec3dp(vx, vy, vz, mx, my, mz)
        auto out = fmt::format_to(ctx.out(), "BiVec3dp(");
        out = fmt::format_to(out, "{}", this->nested(v.vx));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.vy));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.vz));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.mx));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.my));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.mz));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T, Tag> includes DualNumber2dp<T>, DualNumber3dp<T> (PGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::dual_number2dp_tag> ||
             std::is_same_v<Tag, hd::ga::dual_number3dp_tag>)
struct fmt::formatter<hd::ga::MVec2_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec2_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::MVec2_t<T, Tag>,
                                     hd::ga::MVec2_t<T, hd::ga::dual_number2dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "DualNum2dp(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<
                               hd::ga::MVec2_t<T, Tag>,
                               hd::ga::MVec2_t<T, hd::ga::dual_number3dp_tag>>) {
            auto out = fmt::format_to(ctx.out(), "DualNum3dp(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            return fmt::format_to(out, ")");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec4_t<T, Tag> includes MVec2dp_E<T>, MVec2dp_U<T> (PGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::mvec2dp_e_tag> ||
             std::is_same_v<Tag, hd::ga::mvec2dp_u_tag>)
struct fmt::formatter<hd::ga::MVec4_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec4_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::MVec4_t<T, Tag>,
                                     hd::ga::MVec4_t<T, hd::ga::mvec2dp_e_tag>>) {
            auto out = fmt::format_to(ctx.out(), "MVec2dp_E(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c2));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c3));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::MVec4_t<T, Tag>,
                                          hd::ga::MVec4_t<T, hd::ga::mvec2dp_u_tag>>) {
            auto out = fmt::format_to(ctx.out(), "MVec2dp_U(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c2));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c3));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c2));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c3));
            return fmt::format_to(out, ")");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec8_t<T, Tag> includes MVec2dp<T>, MVec3dp_E<T> and MVec3dp_U<T> (PGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::mvec2dp_tag> ||
             std::is_same_v<Tag, hd::ga::mvec3dp_e_tag> ||
             std::is_same_v<Tag, hd::ga::mvec3dp_u_tag>)
struct fmt::formatter<hd::ga::MVec8_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec8_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                     hd::ga::MVec8_t<T, hd::ga::mvec2dp_tag>>) {
            // Format as MVec2dp(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec2dp(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c2));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c3));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c4));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c5));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c6));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c7));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                          hd::ga::MVec8_t<T, hd::ga::mvec3dp_e_tag>>) {
            // Format as MVec3dp_E(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec3dp_E(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c2));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c3));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c4));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c5));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c6));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c7));
            return fmt::format_to(out, ")");
        }
        else if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                          hd::ga::MVec8_t<T, hd::ga::mvec3dp_u_tag>>) {
            // Format as MVec3dp_U(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec3dp_U(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c2));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c3));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c4));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c5));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c6));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c7));
            return fmt::format_to(out, ")");
        }
        else {
            auto out = fmt::format_to(ctx.out(), "(");
            out = fmt::format_to(out, "{}", this->nested(v.c0));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c1));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c2));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c3));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c4));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c5));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c6));
            out = fmt::format_to(out, ", ");
            ctx.advance_to(out);
            out = fmt::format_to(out, "{}", this->nested(v.c7));
            return fmt::format_to(out, ")");
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> includes MVec3dp<T> (PGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::mvec3dp_tag>
struct fmt::formatter<hd::ga::MVec16_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec16_t<T, Tag>& v, FormatContext& ctx) const
    {
        // Format as MVec3dp(c0, c1, c2, ..., c15)
        auto out = fmt::format_to(ctx.out(), "MVec3dp(");
        out = fmt::format_to(out, "{}", this->nested(v.c0));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c1));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c2));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c3));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c4));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c5));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c6));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c7));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c8));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c9));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c10));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c11));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c12));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c13));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c14));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c15));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// PGA Geometric Objects - 2dp Projective Geometry
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Vector2d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Vector2d<T>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::pga::Vector2d<T>& v, FormatContext& ctx) const
    {
        auto out = fmt::format_to(ctx.out(), "Vector2d(");
        out = fmt::format_to(out, "{}", this->nested(v.x));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.y));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// Point2d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Point2d<T>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::pga::Point2d<T>& v, FormatContext& ctx) const
    {
        auto out = fmt::format_to(ctx.out(), "Point2d(");
        out = fmt::format_to(out, "{}", this->nested(v.x));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.y));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// Point2dp<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Point2dp<T>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::pga::Point2dp<T>& v, FormatContext& ctx) const
    {
        auto out = fmt::format_to(ctx.out(), "Point2dp(");
        out = fmt::format_to(out, "{}", this->nested(v.x));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.y));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.z));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// Line2d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Line2d<T>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::pga::Line2d<T>& v, FormatContext& ctx) const
    {
        auto out = fmt::format_to(ctx.out(), "Line2d(");
        out = fmt::format_to(out, "{}", this->nested(v.x));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.y));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.z));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// PGA Geometric Objects - 3dp Projective Geometry
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Vector3d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Vector3d<T>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::pga::Vector3d<T>& v, FormatContext& ctx) const
    {
        auto out = fmt::format_to(ctx.out(), "Vector3d(");
        out = fmt::format_to(out, "{}", this->nested(v.x));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.y));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.z));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// Point3d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Point3d<T>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::pga::Point3d<T>& p, FormatContext& ctx) const
    {
        auto out = fmt::format_to(ctx.out(), "Point3d(");
        out = fmt::format_to(out, "{}", this->nested(p.x));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(p.y));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(p.z));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// Point3dp<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Point3dp<T>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::pga::Point3dp<T>& p, FormatContext& ctx) const
    {
        auto out = fmt::format_to(ctx.out(), "Point3dp(");
        out = fmt::format_to(out, "{}", this->nested(p.x));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(p.y));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(p.z));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(p.w));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// Line3d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Line3d<T>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::pga::Line3d<T>& l, FormatContext& ctx) const
    {
        auto out = fmt::format_to(ctx.out(), "Line3d(");
        out = fmt::format_to(out, "{}", this->nested(l.vx));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(l.vy));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(l.vz));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(l.mx));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(l.my));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(l.mz));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// Plane3d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Plane3d<T>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::pga::Plane3d<T>& p, FormatContext& ctx) const
    {
        auto out = fmt::format_to(ctx.out(), "Plane3d(");
        out = fmt::format_to(out, "{}", this->nested(p.x));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(p.y));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(p.z));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(p.w));
        return fmt::format_to(out, ")");
    }
};