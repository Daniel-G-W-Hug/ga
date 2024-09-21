#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_type_2d.hpp"
#include "type_t/ga_type_2dp.hpp"
#include "type_t/ga_vec2_t.hpp"
#include "type_t/ga_vec3_t.hpp"


namespace hd::ga {

// Vector2d: 2d vector of projective algebra storing only components x, y explicitly.
// The z component is assumed to be z = 0
template <typename T> struct Vector2d : public Vec2d<T> {

    using Vec2d<T>::Vec2d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec2d<T>::x;
    using Vec2d<T>::y;

    Vector2d(Vec2d<T> const& v) : Vec2d<T>(v) {};
};

// Point2d: 2d point of projective algebra storing only components x, y explicitly.
// The z component is assumed to by z=1
template <typename T> struct Point2d : public Vec2d<T> {

    using Vec2d<T>::Vec2d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec2d<T>::x;
    using Vec2d<T>::y;

    Point2d(Vec2d<T> const& v) : Vec2d<T>(v) {};
};

// Point2dp: 2d point of projective algebra storing all three components x, y, z
//           explicitly
//
// a Point2dp is a Vec2dp, thus all operations defined for Vec2dp
// work directly for Point2dp - only deviations will be specified
template <typename T> struct Point2dp : public Vec2dp<T> {

    using Vec2dp<T>::Vec2dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec2dp<T>::x;
    using Vec2dp<T>::y;
    using Vec2dp<T>::z;

    Point2dp(Vec2dp<T> const& v) : Vec2dp<T>(v) {};
    Point2dp(Point2d<T> const& p) : Vec2dp<T>(p.x, p.y, T(1.0)) {};
    Point2dp(Vec2d<T> const& v) : Vec2dp<T>(v.x, v.y, T(1.0)) {};

    Point2dp& unitize()
    {
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (std::abs(z) < std::numeric_limits<T>::epsilon()) {
            throw std::runtime_error("z-component too small for unitization" +
                                     std::to_string(z) + "\n");
        }
#endif
        x /= z;
        y /= z;
        z = T(1.0);
        return (*this);
    }
};

template <typename T>
    requires std::floating_point<T>
Point2dp<T> unitize(Point2dp<T> const& p)
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(p.z) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("z-component too small for unitization" +
                                 std::to_string(p.z) + "\n");
    }
#endif
    T inv = T(1.0) / p.z;
    return Point2dp<T>(p.x * inv, p.y * inv, T(1.0));
}


// Line2dp: 2d line of projective algebra storing all three components
//          explicitly as components x, y, z of a BiVec2dp
//
// a Line2dp is a BiVec2dp, thus all operations defined for BiVec2dp
// work directly for Line2dp - only deviations will be specified
template <typename T> struct Line2dp : public BiVec2dp<T> {

    using BiVec2dp<T>::BiVec2dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using BiVec2dp<T>::x;
    using BiVec2dp<T>::y;
    using BiVec2dp<T>::z;

    Line2dp() = default;
    Line2dp(T x, T y, T z) : BiVec2dp<T>(x, y, z) {};
    Line2dp(BiVec2dp<T> const& b) : BiVec2dp<T>(b) {};
    Line2dp(Point2d<T> const& p, Point2d<T> const& q) :
        // Line2dp = wdg(p,q), but wdg() cannot be used here to avoid circular dependency
        BiVec2dp<T>(p.y - q.y, q.x - p.x, p.x * q.y - p.y * q.x) {};
    Line2dp(Point2d<T> const& p, Vec2d<T> const& v) :
        // Line2dp constructed from a point and a direction vector
        BiVec2dp<T>(-v.y, v.x, p.x * v.y - p.y * v.x) {};

    Line2dp& unitize()
    {
        // unitization for a 2d bivector means std::sqrt(x^2 + y^2) = 1
        T wn = std::sqrt(x * x + y * y);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (wn < std::numeric_limits<T>::epsilon()) {
            throw std::runtime_error("bivector weight norm too small for unitization" +
                                     std::to_string(wn) + "\n");
        }
#endif
        T inv = T(1.0) / wn;
        x *= inv;
        y *= inv;
        z *= inv;
        return (*this);
    }
};

template <typename T>
    requires std::floating_point<T>
Line2dp<T> unitize(Line2dp<T> const& l)
{
    // unitization for a 2d bivector means std::sqrt(x^2 + y^2) = 1
    T wn = std::sqrt(l.x * l.x + l.y * l.y);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (wn < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector weight norm too small for unitization" +
                                 std::to_string(wn) + "\n");
    }
#endif
    T inv = T(1.0) / wn;
    return Line2dp<T>(l.x * inv, l.y * inv, l.z * inv);
}

} // namespace hd::ga