#pragma once

// author: Daniel Hug, 2024

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "type_t/ga_scalar_t.hpp"
#include "type_t/ga_vec2_t.hpp"
#include "type_t/ga_vec3_t.hpp"
#include "type_t/ga_vec4_t.hpp"

#include "type_t/ga_bvec6_t.hpp"

#include "type_t/ga_mvec16_t.hpp"
#include "type_t/ga_mvec2_t.hpp"
#include "type_t/ga_mvec4_t.hpp"
#include "type_t/ga_mvec8_t.hpp"

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
// Scalar_t<T, Tag>: Scalar2d<T>, PScalar2d<T>, PScalar3d<T>, PScalar4d<T>, etc.
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::Scalar_t<T, Tag>> : formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Scalar_t<T, Tag>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({})", double(v));
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
        return fmt::format_to(ctx.out(), "({},{})", nested(v.x), nested(v.y));
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec3_t<T, Tag> includes Vec3d<T> and BiVec3d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::Vec3_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec3_t<T, Tag>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{})", nested(v.x), nested(v.y),
                              nested(v.z));
    }
};

////////////////////////////////////////////////////////////////////////////////
// Point2dp<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::Point2dp<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Point2dp<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{})", nested(v.x), nested(v.y),
                              nested(v.z));
    }
};

////////////////////////////////////////////////////////////////////////////////
// Vec4_t<T, Tag> includes Vec4d<T> and TriVec4d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::Vec4_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Vec4_t<T, Tag>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{})", nested(v.x), nested(v.y),
                              nested(v.z), nested(v.w));
    }
};

////////////////////////////////////////////////////////////////////////////////
// BVec6_t<T, Tag> - includes BiVec4d<T>, BiVec3dp<T>, Line3dp<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::BVec6_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::BVec6_t<T, Tag>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{},{},{})", nested(v.vx),
                              nested(v.vy), nested(v.vz), nested(v.mx), nested(v.my),
                              nested(v.mz));
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
        return fmt::format_to(ctx.out(), "({},{})", nested(v.c0), nested(v.c1));
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec4_t<T, Tag> includes MVec2d<T>, MVec3d_E<T> and MVec3d_U<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::MVec4_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec4_t<T, Tag>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{})", nested(v.c0), nested(v.c1),
                              nested(v.c2), nested(v.c3));
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec8_t<T, Tag> includes MVec3d<T>, MVec4d_E<T> and MVec4d_U<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::MVec8_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec8_t<T, Tag>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{},{},{},{},{})", nested(v.c0),
                              nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
                              nested(v.c5), nested(v.c6), nested(v.c7));
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec16_t<T, Tag> includes MVec4d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::MVec16_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec16_t<T, Tag>& v, FormatContext& ctx) const
    {
        return fmt::format_to(
            ctx.out(), "({},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{})", nested(v.c0),
            nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4), nested(v.c5),
            nested(v.c6), nested(v.c7), nested(v.c8), nested(v.c9), nested(v.c10),
            nested(v.c11), nested(v.c12), nested(v.c13), nested(v.c14), nested(v.c15));
    }
};