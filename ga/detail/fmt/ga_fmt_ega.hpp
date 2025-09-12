#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_fmt_core.hpp"

#include "../type_t/ga_bvec6_t.hpp"
#include "../type_t/ga_mvec16_t.hpp"
#include "../type_t/ga_mvec2_t.hpp"
#include "../type_t/ga_mvec4_t.hpp"
#include "../type_t/ga_mvec8_t.hpp"
#include "../type_t/ga_scalar_t.hpp"
#include "../type_t/ga_vec2_t.hpp"
#include "../type_t/ga_vec3_t.hpp"
#include "../type_t/ga_vec4_t.hpp"

/////////////////////////////////////////////////////////////////////////////////////////
// Formatting support for Euclidean Geometric Algebra (EGA) types
// Includes formatters for scalars, vectors, bivectors, and multivectors
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag>: Scalar2d<T>, PScalar2d<T>, Scalar3d<T>, PScalar3d<T>, etc.
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::Scalar_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::Scalar_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                     hd::ga::Scalar_t<T, hd::ga::scalar2d_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar2d({})", this->nested(T(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar2d_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar2d({})", this->nested(T(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::scalar3d_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar3d({})", this->nested(T(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar3d_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar3d({})", this->nested(T(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::scalar4d_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar4d({})", this->nested(T(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar4d_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar4d({})", this->nested(T(v)));
        }
        else {
            return fmt::format_to(ctx.out(), "({})", this->nested(T(v)));
        }
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec2_t<T, Tag> includes Vec2d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::vec2d_tag>
struct fmt::formatter<hd::ga::Vec2_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::Vec2_t<T, Tag>& v, FormatContext& ctx) const
    {
        // Format as Vec2d(x, y)
        auto out = fmt::format_to(ctx.out(), "Vec2d(");

        // Use nested formatter for x component
        out = fmt::format_to(out, "{}", this->nested(v.x));

        out = fmt::format_to(out, ", ");

        // Update context and use nested formatter for y component
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.y));

        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec3_t<T, Tag> includes Vec3d<T>, BiVec3d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::vec3d_tag> ||
             std::is_same_v<Tag, hd::ga::bivec3d_tag>)
struct fmt::formatter<hd::ga::Vec3_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::Vec3_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::Vec3_t<T, Tag>,
                                     hd::ga::Vec3_t<T, hd::ga::vec3d_tag>>) {
            // Format as Vec3d(x, y, z)
            auto out = fmt::format_to(ctx.out(), "Vec3d(");
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
                                          hd::ga::Vec3_t<T, hd::ga::bivec3d_tag>>) {
            // Format as BiVec3d(x, y, z)
            auto out = fmt::format_to(ctx.out(), "BiVec3d(");
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
// Vec4_t<T, Tag> includes Vec4d<T>, TriVec4d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::vec4d_tag> ||
             std::is_same_v<Tag, hd::ga::trivec4d_tag>)
struct fmt::formatter<hd::ga::Vec4_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::Vec4_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
                                     hd::ga::Vec4_t<T, hd::ga::vec4d_tag>>) {
            // Format as Vec4d(x, y, z, w)
            auto out = fmt::format_to(ctx.out(), "Vec4d(");
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
                                          hd::ga::Vec4_t<T, hd::ga::trivec4d_tag>>) {
            // Format as TriVec4d(x, y, z, w)
            auto out = fmt::format_to(ctx.out(), "TriVec4d(");
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
// BVec6_t<T, Tag> - includes BiVec4d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::bivec4d_tag>
struct fmt::formatter<hd::ga::BVec6_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::BVec6_t<T, Tag>& v, FormatContext& ctx) const
    {
        // Format as BiVec4d(vx, vy, vz, mx, my, mz)
        auto out = fmt::format_to(ctx.out(), "BiVec4d(");
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
// MVec2_t<T, Tag> includes MVec2d_E<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::mvec2d_e_tag>
struct fmt::formatter<hd::ga::MVec2_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec2_t<T, Tag>& v, FormatContext& ctx) const
    {
        // Format as MVec2d_E(c0, c1)
        auto out = fmt::format_to(ctx.out(), "MVec2d_E(");
        out = fmt::format_to(out, "{}", this->nested(v.c0));
        out = fmt::format_to(out, ", ");
        ctx.advance_to(out);
        out = fmt::format_to(out, "{}", this->nested(v.c1));
        return fmt::format_to(out, ")");
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec4_t<T, Tag> includes MVec2d<T>, MVec3d_E<T>, MVec3d_U<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::mvec2d_tag> ||
             std::is_same_v<Tag, hd::ga::mvec3d_e_tag> ||
             std::is_same_v<Tag, hd::ga::mvec3d_u_tag>)
struct fmt::formatter<hd::ga::MVec4_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec4_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::MVec4_t<T, Tag>,
                                     hd::ga::MVec4_t<T, hd::ga::mvec2d_tag>>) {
            // Format as MVec2d(c0, c1, c2, c3)
            auto out = fmt::format_to(ctx.out(), "MVec2d(");
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
                                          hd::ga::MVec4_t<T, hd::ga::mvec3d_e_tag>>) {
            // Format as MVec3d_E(c0, c1, c2, c3)
            auto out = fmt::format_to(ctx.out(), "MVec3d_E(");
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
                                          hd::ga::MVec4_t<T, hd::ga::mvec3d_u_tag>>) {
            // Format as MVec3d_U(c0, c1, c2, c3)
            auto out = fmt::format_to(ctx.out(), "MVec3d_U(");
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
// MVec8_t<T, Tag> includes MVec3d<T>, MVec4d_E<T>, MVec4d_U<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::mvec3d_tag> ||
             std::is_same_v<Tag, hd::ga::mvec4d_e_tag> ||
             std::is_same_v<Tag, hd::ga::mvec4d_u_tag>)
struct fmt::formatter<hd::ga::MVec8_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec8_t<T, Tag>& v, FormatContext& ctx) const
    {
        if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                     hd::ga::MVec8_t<T, hd::ga::mvec3d_tag>>) {
            // Format as MVec3d(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec3d(");
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
                                          hd::ga::MVec8_t<T, hd::ga::mvec4d_e_tag>>) {
            // Format as MVec4d_E(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec4d_E(");
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
                                          hd::ga::MVec8_t<T, hd::ga::mvec4d_u_tag>>) {
            // Format as MVec4d_U(c0, c1, c2, c3, c4, c5, c6, c7)
            auto out = fmt::format_to(ctx.out(), "MVec4d_U(");
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
// MVec16_t<T, Tag> includes MVec4d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::mvec4d_tag>
struct fmt::formatter<hd::ga::MVec16_t<T, Tag>> : nested_formatter<T> {
    constexpr auto parse(format_parse_context& ctx)
    {
        return fmt::nested_formatter<T>::parse(ctx);
    }

    template <typename FormatContext>
    auto format(const hd::ga::MVec16_t<T, Tag>& v, FormatContext& ctx) const
    {
        // Format as MVec4d(c0, c1, c2, ..., c15)
        auto out = fmt::format_to(ctx.out(), "MVec4d(");
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