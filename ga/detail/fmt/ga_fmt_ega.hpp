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
struct fmt::formatter<hd::ga::Scalar_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Scalar_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                     hd::ga::Scalar_t<T, hd::ga::scalar2d_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar2d({})", nested(double(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar2d_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar2d({})", nested(double(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::scalar3d_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar3d({})", nested(double(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar3d_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar3d({})", nested(double(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::scalar4d_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar4d({})", nested(double(v)));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar4d_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar4d({})", nested(double(v)));
        }
        else {
            return fmt::format_to(ctx.out(), "({})", nested(double(v)));
        }
#else
        return fmt::format_to(ctx.out(), "({})", nested(double(v)));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec2_t<T, Tag> includes Vec2d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::vec2d_tag>
struct fmt::formatter<hd::ga::Vec2_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec2_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        return fmt::format_to(ctx.out(), "Vec2d({}, {})", nested(v.x), nested(v.y));
#else
        return fmt::format_to(ctx.out(), "({}, {})", nested(v.x), nested(v.y));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec3_t<T, Tag> includes Vec3d<T>, BiVec3d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::vec3d_tag> ||
             std::is_same_v<Tag, hd::ga::bivec3d_tag>)
struct fmt::formatter<hd::ga::Vec3_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec3_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::Vec3_t<T, Tag>,
                                     hd::ga::Vec3_t<T, hd::ga::vec3d_tag>>) {
            return fmt::format_to(ctx.out(), "Vec3d({}, {}, {})", nested(v.x),
                                  nested(v.y), nested(v.z));
        }
        else if constexpr (std::is_same_v<hd::ga::Vec3_t<T, Tag>,
                                          hd::ga::Vec3_t<T, hd::ga::bivec3d_tag>>) {
            return fmt::format_to(ctx.out(), "BiVec3d({}, {}, {})", nested(v.x),
                                  nested(v.y), nested(v.z));
        }
        else {
            return fmt::format_to(ctx.out(), "({}, {}, {})", nested(v.x), nested(v.y),
                                  nested(v.z));
        }
#else
        return fmt::format_to(ctx.out(), "({}, {}, {})", nested(v.x), nested(v.y),
                              nested(v.z));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec4_t<T, Tag> includes Vec4d<T>, TriVec4d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::vec4d_tag> ||
             std::is_same_v<Tag, hd::ga::trivec4d_tag>)
struct fmt::formatter<hd::ga::Vec4_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec4_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
                                     hd::ga::Vec4_t<T, hd::ga::vec4d_tag>>) {
            return fmt::format_to(ctx.out(), "Vec4d({}, {}, {}, {})", nested(v.x),
                                  nested(v.y), nested(v.z), nested(v.w));
        }
        else if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
                                          hd::ga::Vec4_t<T, hd::ga::trivec4d_tag>>) {
            return fmt::format_to(ctx.out(), "TriVec4d({}, {}, {}, {})", nested(v.x),
                                  nested(v.y), nested(v.z), nested(v.w));
        }
        else {
            return fmt::format_to(ctx.out(), "({}, {}, {}, {})", nested(v.x), nested(v.y),
                                  nested(v.z), nested(v.w));
        }
#else
        return fmt::format_to(ctx.out(), "({}, {}, {}, {})", nested(v.x), nested(v.y),
                              nested(v.z), nested(v.w));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T, Tag> - includes BiVec4d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::bivec4d_tag>
struct fmt::formatter<hd::ga::BVec6_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::BVec6_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        return fmt::format_to(ctx.out(), "BiVec4d({}, {}, {}, {}, {}, {})", nested(v.vx),
                              nested(v.vy), nested(v.vz), nested(v.mx), nested(v.my),
                              nested(v.mz));
#else
        return fmt::format_to(ctx.out(), "({}, {}, {}, {}, {}, {})", nested(v.vx),
                              nested(v.vy), nested(v.vz), nested(v.mx), nested(v.my),
                              nested(v.mz));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T, Tag> includes MVec2d_E<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::mvec2d_e_tag>
struct fmt::formatter<hd::ga::MVec2_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec2_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        return fmt::format_to(ctx.out(), "MVec2d_E({}, {})", nested(v.c0), nested(v.c1));
#else
        return fmt::format_to(ctx.out(), "({}, {})", nested(v.c0), nested(v.c1));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec4_t<T, Tag> includes MVec2d<T>, MVec3d_E<T>, MVec3d_U<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::mvec2d_tag> ||
             std::is_same_v<Tag, hd::ga::mvec3d_e_tag> ||
             std::is_same_v<Tag, hd::ga::mvec3d_u_tag>)
struct fmt::formatter<hd::ga::MVec4_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec4_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::MVec4_t<T, Tag>,
                                     hd::ga::MVec4_t<T, hd::ga::mvec2d_tag>>) {
            return fmt::format_to(ctx.out(), "MVec2d({}, {}, {}, {})", nested(v.c0),
                                  nested(v.c1), nested(v.c2), nested(v.c3));
        }
        else if constexpr (std::is_same_v<hd::ga::MVec4_t<T, Tag>,
                                          hd::ga::MVec4_t<T, hd::ga::mvec3d_e_tag>>) {
            return fmt::format_to(ctx.out(), "MVec3d_E({}, {}, {}, {})", nested(v.c0),
                                  nested(v.c1), nested(v.c2), nested(v.c3));
        }
        else if constexpr (std::is_same_v<hd::ga::MVec4_t<T, Tag>,
                                          hd::ga::MVec4_t<T, hd::ga::mvec3d_u_tag>>) {
            return fmt::format_to(ctx.out(), "MVec3d_U({}, {}, {}, {})", nested(v.c0),
                                  nested(v.c1), nested(v.c2), nested(v.c3));
        }
        else {
            return fmt::format_to(ctx.out(), "({}, {}, {}, {})", nested(v.c0),
                                  nested(v.c1), nested(v.c2), nested(v.c3));
        }
#else
        return fmt::format_to(ctx.out(), "({}, {}, {}, {})", nested(v.c0), nested(v.c1),
                              nested(v.c2), nested(v.c3));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec8_t<T, Tag> includes MVec3d<T>, MVec4d_E<T>, MVec4d_U<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::is_same_v<Tag, hd::ga::mvec3d_tag> ||
             std::is_same_v<Tag, hd::ga::mvec4d_e_tag> ||
             std::is_same_v<Tag, hd::ga::mvec4d_u_tag>)
struct fmt::formatter<hd::ga::MVec8_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec8_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                     hd::ga::MVec8_t<T, hd::ga::mvec3d_tag>>) {
            return fmt::format_to(ctx.out(), "MVec3d({}, {}, {}, {}, {}, {}, {}, {})",
                                  nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3),
                                  nested(v.c4), nested(v.c5), nested(v.c6), nested(v.c7));
        }
        else if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                          hd::ga::MVec8_t<T, hd::ga::mvec4d_e_tag>>) {
            return fmt::format_to(ctx.out(), "MVec4d_E({}, {}, {}, {}, {}, {}, {}, {})",
                                  nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3),
                                  nested(v.c4), nested(v.c5), nested(v.c6), nested(v.c7));
        }
        else if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                          hd::ga::MVec8_t<T, hd::ga::mvec4d_u_tag>>) {
            return fmt::format_to(ctx.out(), "MVec4d_U({}, {}, {}, {}, {}, {}, {}, {})",
                                  nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3),
                                  nested(v.c4), nested(v.c5), nested(v.c6), nested(v.c7));
        }
        else {
            return fmt::format_to(ctx.out(), "({}, {}, {}, {}, {}, {}, {}, {})",
                                  nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3),
                                  nested(v.c4), nested(v.c5), nested(v.c6), nested(v.c7));
        }
#else
        return fmt::format_to(ctx.out(), "({}, {}, {}, {}, {}, {}, {}, {})", nested(v.c0),
                              nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
                              nested(v.c5), nested(v.c6), nested(v.c7));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> includes MVec4d<T> (EGA types only)
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires std::is_same_v<Tag, hd::ga::mvec4d_tag>
struct fmt::formatter<hd::ga::MVec16_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec16_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        return fmt::format_to(
            ctx.out(),
            "MVec4d({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})",
            nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
            nested(v.c5), nested(v.c6), nested(v.c7), nested(v.c8), nested(v.c9),
            nested(v.c10), nested(v.c11), nested(v.c12), nested(v.c13), nested(v.c14),
            nested(v.c15));
#else
        return fmt::format_to(
            ctx.out(), "({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})",
            nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
            nested(v.c5), nested(v.c6), nested(v.c7), nested(v.c8), nested(v.c9),
            nested(v.c10), nested(v.c11), nested(v.c12), nested(v.c13), nested(v.c14),
            nested(v.c15));
#endif
    }
};