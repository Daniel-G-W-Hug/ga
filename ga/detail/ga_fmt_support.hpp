#pragma once

// author: Daniel Hug, 2024 & 2025

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "type_t/ga_type_tags.hpp"
#include <type_traits>

#include "type_t/ga_scalar_t.hpp"
#include "type_t/ga_vec2_t.hpp"
#include "type_t/ga_vec3_t.hpp"
#include "type_t/ga_vec4_t.hpp"

#include "type_t/ga_bvec6_t.hpp"

#include "type_t/ga_mvec16_t.hpp"
#include "type_t/ga_mvec2_t.hpp"
#include "type_t/ga_mvec4_t.hpp"
#include "type_t/ga_mvec8_t.hpp"

#include "ga_pga_2dp_objects.hpp"

////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library for all types defined in ga library
////////////////////////////////////////////////////////////////////////////////
// Usage example:
//
// std::vector<BiVec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
// BiVec3d p{1.0, 2.0, 3.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));

////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag>: Scalar2d<T>, PScalar2d<T>, Scalar3d<T>, PScalar3d<T>, etc.
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::Scalar_t<T, Tag>> : formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Scalar_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                     hd::ga::Scalar_t<T, hd::ga::scalar2d_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar2d({})", double(v));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar2d_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar2d({})", double(v));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::scalar2dp_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar2dp({})", double(v));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar2dp_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar2d({})", double(v));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::scalar3d_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar3d({})", double(v));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar3d_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar3d({})", double(v));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::scalar3dp_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar3dp({})", double(v));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar3dp_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar3dp({})", double(v));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::scalar4d_tag>>) {
            return fmt::format_to(ctx.out(), "Scalar4d({})", double(v));
        }
        else if constexpr (std::is_same_v<hd::ga::Scalar_t<T, Tag>,
                                          hd::ga::Scalar_t<T, hd::ga::pscalar4d_tag>>) {
            return fmt::format_to(ctx.out(), "PScalar4d({})", double(v));
        }
        else {
            return fmt::format_to(ctx.out(), "({})", double(v));
        }
#else
        return fmt::format_to(ctx.out(), "({})", double(v));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec2_t<T, Tag> includes Vec2d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::Vec2_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec2_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::Vec2_t<T, Tag>,
                                     hd::ga::Vec2_t<T, hd::ga::vec2d_tag>>) {
            return fmt::format_to(ctx.out(), "Vec2d({}, {})", nested(v.x), nested(v.y));
        }
        else {
            return fmt::format_to(ctx.out(), "({}, {})", nested(v.x), nested(v.y));
        }
#else
        return fmt::format_to(ctx.out(), "({}, {})", nested(v.x), nested(v.y));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec3_t<T, Tag> includes Vec3d<T>, BiVec3d<T>, Vec2dp<T>, BiVec2dp<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
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
        else if constexpr (std::is_same_v<hd::ga::Vec3_t<T, Tag>,
                                          hd::ga::Vec3_t<T, hd::ga::vec2dp_tag>>) {
            return fmt::format_to(ctx.out(), "Vec2dp({}, {}, {})", nested(v.x),
                                  nested(v.y), nested(v.z));
        }
        else if constexpr (std::is_same_v<hd::ga::Vec3_t<T, Tag>,
                                          hd::ga::Vec3_t<T, hd::ga::bivec2dp_tag>>) {
            return fmt::format_to(ctx.out(), "BiVec2dp({}, {}, {})", nested(v.x),
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
// Vec4_t<T, Tag> includes Vec3dp<T>, TriVec3dp<T>, Vec4d<T>, TriVec4d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::Vec4_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec4_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
                                     hd::ga::Vec4_t<T, hd::ga::vec3dp_tag>>) {
            return fmt::format_to(ctx.out(), "Vec3dp({}, {}, {}, {})", nested(v.x),
                                  nested(v.y), nested(v.z), nested(v.w));
        }
        else if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
                                          hd::ga::Vec4_t<T, hd::ga::trivec3dp_tag>>) {
            return fmt::format_to(ctx.out(), "TriVec3dp({}, {}, {}, {})", nested(v.x),
                                  nested(v.y), nested(v.z), nested(v.w));
        }
        else if constexpr (std::is_same_v<hd::ga::Vec4_t<T, Tag>,
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
// BVec6_t<T, Tag> - includes BiVec3dp<T>, BiVec4d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::BVec6_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::BVec6_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::BVec6_t<T, Tag>,
                                     hd::ga::BVec6_t<T, hd::ga::bivec3dp_tag>>) {
            return fmt::format_to(ctx.out(), "BiVec3dp({}, {}, {}, {}, {}, {})",
                                  nested(v.vx), nested(v.vy), nested(v.vz), nested(v.mx),
                                  nested(v.my), nested(v.mz));
        }
        else if constexpr (std::is_same_v<hd::ga::BVec6_t<T, Tag>,
                                          hd::ga::BVec6_t<T, hd::ga::bivec4d_tag>>) {
            return fmt::format_to(ctx.out(), "BiVec4d({}, {}, {}, {}, {}, {})",
                                  nested(v.vx), nested(v.vy), nested(v.vz), nested(v.mx),
                                  nested(v.my), nested(v.mz));
        }
        else {
            return fmt::format_to(ctx.out(), "({}, {}, {}, {}, {}, {})", nested(v.vx),
                                  nested(v.vy), nested(v.vz), nested(v.mx), nested(v.my),
                                  nested(v.mz));
        }
#else
        return fmt::format_to(ctx.out(), "({}, {}, {}, {}, {}, {})", nested(v.vx),
                              nested(v.vy), nested(v.vz), nested(v.mx), nested(v.my),
                              nested(v.mz));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec2_t<T, Tag> includes MVec2d_E<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::MVec2_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec2_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::MVec2_t<T, Tag>,
                                     hd::ga::MVec2_t<T, hd::ga::mvec2d_e_tag>>) {
            return fmt::format_to(ctx.out(), "MVec2d_E({}, {})", nested(v.c0),
                                  nested(v.c1));
        }
        else {
            return fmt::format_to(ctx.out(), "({}, {})", nested(v.c0), nested(v.c1));
        }
#else
        return fmt::format_to(ctx.out(), "({}, {})", nested(v.c0), nested(v.c1));
#endif
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec4_t<T, Tag> includes MVec2d<T>, MVec2dp_E<T>, MVec2dp_U<T>,
//                          MVec3d_E<T>, MVec3d_U<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
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
                                          hd::ga::MVec4_t<T, hd::ga::mvec2dp_e_tag>>) {
            return fmt::format_to(ctx.out(), "MVec2dp_E({}, {}, {}, {})", nested(v.c0),
                                  nested(v.c1), nested(v.c2), nested(v.c3));
        }
        else if constexpr (std::is_same_v<hd::ga::MVec4_t<T, Tag>,
                                          hd::ga::MVec4_t<T, hd::ga::mvec2dp_u_tag>>) {
            return fmt::format_to(ctx.out(), "MVec2dp_U({}, {}, {}, {})", nested(v.c0),
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
// MVec8_t<T, Tag> includes MVec3d<T>, MVec3dp_E<T> and MVec3dp_U<T>,
//                          MVec4d_E<T>, MVec4d_U<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
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
                                          hd::ga::MVec8_t<T, hd::ga::mvec3dp_e_tag>>) {
            return fmt::format_to(ctx.out(), "MVec3dp_E({}, {}, {}, {}, {}, {}, {}, {})",
                                  nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3),
                                  nested(v.c4), nested(v.c5), nested(v.c6), nested(v.c7));
        }
        else if constexpr (std::is_same_v<hd::ga::MVec8_t<T, Tag>,
                                          hd::ga::MVec8_t<T, hd::ga::mvec3dp_u_tag>>) {
            return fmt::format_to(ctx.out(), "MVec3dp_U({}, {}, {}, {}, {}, {}, {}, {})",
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
// MVec16_t<T, Tag> includes MVec3dp<T>, MVec4d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::MVec16_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec16_t<T, Tag>& v, FormatContext& ctx) const
    {
#if defined(_HD_GA_PRINT_WITH_TYPE_INFO)
        if constexpr (std::is_same_v<hd::ga::MVec16_t<T, Tag>,
                                     hd::ga::MVec16_t<T, hd::ga::mvec3dp_tag>>) {
            return fmt::format_to(
                ctx.out(),
                "MVec3dp({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})",
                nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
                nested(v.c5), nested(v.c6), nested(v.c7), nested(v.c8), nested(v.c9),
                nested(v.c10), nested(v.c11), nested(v.c12), nested(v.c13), nested(v.c14),
                nested(v.c15));
        }
        else if constexpr (std::is_same_v<hd::ga::MVec16_t<T, Tag>,
                                          hd::ga::MVec16_t<T, hd::ga::mvec4d_tag>>) {
            return fmt::format_to(
                ctx.out(),
                "MVec4d({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})",
                nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
                nested(v.c5), nested(v.c6), nested(v.c7), nested(v.c8), nested(v.c9),
                nested(v.c10), nested(v.c11), nested(v.c12), nested(v.c13), nested(v.c14),
                nested(v.c15));
        }
        else {
            return fmt::format_to(
                ctx.out(),
                "({}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {})",
                nested(v.c0), nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
                nested(v.c5), nested(v.c6), nested(v.c7), nested(v.c8), nested(v.c9),
                nested(v.c10), nested(v.c11), nested(v.c12), nested(v.c13), nested(v.c14),
                nested(v.c15));
        }
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

////////////////////////////////////////////////////////////////////////////////
// Vector2d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Vector2d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::pga::Vector2d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "Vector2d({}, {})", nested(v.x), nested(v.y));
    }
};

////////////////////////////////////////////////////////////////////////////////
// Point2d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Point2d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::pga::Point2d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "Point2d({}, {})", nested(v.x), nested(v.y));
    }
};

////////////////////////////////////////////////////////////////////////////////
// Point2dp<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Point2dp<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::pga::Point2dp<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "Point2dp({}, {}, {})", nested(v.x), nested(v.y),
                              nested(v.z));
    }
};

////////////////////////////////////////////////////////////////////////////////
// Line2dp<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::pga::Line2dp<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::pga::Line2dp<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "Line2dp({}, {}, {})", nested(v.x), nested(v.y),
                              nested(v.z));
    }
};