#pragma once

// author: Daniel Hug, 2024

////////////////////////////////////////////////////////////////////////////////
// printing support via fmt library for all types defined in ga library
////////////////////////////////////////////////////////////////////////////////

// Usage example:
//
// std::vector<BiVec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
// BiVec3d p{1.0, 2.0, 3.0};
// fmt::print(" p = {}\n", p);
// fmt::print(" vp1 = {}\n", fmt::join(vp1, ", "));

#include "fmt/format.h"
#include "fmt/ranges.h" // support printing of (nested) containers & tuples

#include "ga_scalar_t.hpp"
#include "ga_vec2_t.hpp"
#include "ga_vec3_t.hpp"
#include "ga_vec4_t.hpp"

#include "ga_mvec2d.hpp"
#include "ga_mvec2d_e.hpp"

#include "ga_mvec3d_e.hpp"
#include "ga_mvec3d_u.hpp"

#include "ga_mvec16_t.hpp"
#include "ga_mvec6_t.hpp"
#include "ga_mvec8_t.hpp"

////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag> includes Scalar<T>, PScalar2d<T>, PScalar3d<T>, PScalar4d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::Scalar_t<T, Tag>> : formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::Scalar_t<T, Tag>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "{}", double(v));
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
// MVec2d<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::MVec2d<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec2d<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{})", nested(v.c0), nested(v.c1),
                              nested(v.c2), nested(v.c3));
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec2d_E<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::MVec2d_E<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec2d_E<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{})", nested(v.c0), nested(v.c1));
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
// MVec3d_E<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::MVec3d_E<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec3d_E<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{})", nested(v.c0), nested(v.c1),
                              nested(v.c2), nested(v.c3));
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec3d_U<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T>
struct fmt::formatter<hd::ga::MVec3d_U<T>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec3d_U<T>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{})", nested(v.c0), nested(v.c1),
                              nested(v.c2), nested(v.c3));
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
// MVec6_t<T, Tag> - includes BiVec4d<T>, BiVec3dp<T>, Line3dp<T>
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
struct fmt::formatter<hd::ga::MVec6_t<T, Tag>> : nested_formatter<double> {
    template <typename FormatContext>
    auto format(const hd::ga::MVec6_t<T, Tag>& v, FormatContext& ctx) const
    {
        return fmt::format_to(ctx.out(), "({},{},{},{},{},{})", nested(v.c0),
                              nested(v.c1), nested(v.c2), nested(v.c3), nested(v.c4),
                              nested(v.c5));
    }
};

////////////////////////////////////////////////////////////////////////////////
// MVec8_t<T, Tag> includes MVec3d<T>
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
// MVec16_t<T, Tag> includes MVec3d<T>
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