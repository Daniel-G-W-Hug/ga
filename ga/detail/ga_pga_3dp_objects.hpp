#pragma once

// author: Daniel Hug, 2024 & 2025

#include "type_t/ga_type_3d.hpp"
#include "type_t/ga_type_3dp.hpp"
#include "type_t/ga_vec2_t.hpp"
#include "type_t/ga_vec3_t.hpp"


namespace hd::ga::pga {

////////////////////////////////////////////////////////////////////////////////
// bulk
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> bulk(Vec3dp<T> const& v)
{
    return Vec3dp<T>(v.x, v.y, v.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> bulk(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), B.mx, B.my, B.mz);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> bulk(TriVec3dp<T> const& t)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), t.w);
}


////////////////////////////////////////////////////////////////////////////////
// weight
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> weight(Vec3dp<T> const& v)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), v.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> weight(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(B.vx, B.vy, B.vz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> weight(TriVec3dp<T> const& t)
{
    return TriVec3dp<T>(t.x, t.y, t.z, T(0.0));
}


////////////////////////////////////////////////////////////////////////////////
// bulk norm
////////////////////////////////////////////////////////////////////////////////

// return squared bulk norm of scalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(Scalar3dp<T> s)
{
    // ||bulk(s)||^2 = dot(s,s) = s^(T) G s = gr0(s*rev(s)) = s*s
    // using rev(s) = (-1)^[k(k-1)/2] s for a k-blade: 0-blade => rev(s) = s
    return Scalar3dp<T>(T(s) * T(s));
}

// return bulk norm of scalar
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(Scalar3dp<T> s)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(s)));
}

// return squared bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(Vec3dp<T> const& v)
{
    // ||bulk(v)||^2 = dot(v,v) = v^(T) G v = gr0(v*rev(v))
    //               = v.x * v.x + v.y * v.y + v.z * v.z
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    return Scalar3dp<T>(v.x * v.x + v.y * v.y + v.z * v.z);
}

// return bulk norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(Vec3dp<T> const& v)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(v)));
}

// return squared bulk magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(BiVec3dp<T> const& B)
{
    // ||bulk(B)||^2 = dot(B,B) = B^(T) G B = gr0(B*rev(B))
    //               = B.mx * B.mx + B.my * B.my + B.mz * B.mz
    // using rev(B) = (-1)^[k(k-1)/2] B for a k-blade: 2-blade => rev(B) = -B
    return Scalar3dp<T>(B.mx * B.mx + B.my * B.my + B.mz * B.mz);
}

// return magnitude of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(BiVec3dp<T> const& B)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(B)));
}

// return squared bulk norm of trivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm_sq(TriVec3dp<T> const& t)
{
    // ||bulk(t)||^2 = dot(t,t) = t^(T) G t = gr0(t*rev(t)) = t.w * t.w
    // using rev(t) = (-1)^[k(k-1)/2] t for a k-blade: 3-blade => rev(t) = -t
    return Scalar3dp<T>(t.w * t.w);
}

// return bulk norm of trivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_nrm(TriVec3dp<T> const& t)
{
    return Scalar3dp<T>(std::sqrt(bulk_nrm_sq(t)));
}

////////////////////////////////////////////////////////////////////////////////
// weight norm
////////////////////////////////////////////////////////////////////////////////

// return squared weight norm of vector
// ||weight(v)||^2 = lcmpl( dot(rcmpl(v), rcmpl(v))) ) = rdot(v,v) = v.w * v.w
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(Vec3dp<T> const& v)
{
    return PScalar3dp<T>(v.w * v.w);
}

// return weigth norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(Vec3dp<T> const& v)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(v)));
}

// return squared weight norm of bivector
// ||weight(B)||^2 = lcmpl( dot(rcmpl(B),cmpl(B))) ) = rdot(B, B)
//                 = B.vx * B.vx + B.vy * B.vy + B.vz * B.vz
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(BiVec3dp<T> const& B)
{
    return PScalar3dp<T>(B.vx * B.vx + B.vy * B.vy + B.vz * B.vz);
}

// return weight norm of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(BiVec3dp<T> const& B)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(B)));
}

// return squared weight norm of trivector
// ||weight(t)||^2 = lcmpl( dot(rcmpl(t),cmpl(t))) ) = rdot(t, t)
//                 = t.x * t.x + t.y * t.y + t.z * t.z
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm_sq(TriVec3dp<T> const& t)
{
    return PScalar3dp<T>(t.x * t.x + t.y * t.y + t.z * t.z);
}

// return weight norm of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_nrm(TriVec3dp<T> const& t)
{
    return PScalar3dp<T>(std::sqrt(weight_nrm_sq(t)));
}


////////////////////////////////////////////////////////////////////////////////
// geometric norm: (perpendicular) distance to the origin (distance = c0/c1)
//
// returns a dual number, for correct handling of objects at infinity (c1 == 0)
//
////////////////////////////////////////////////////////////////////////////////

// provide the distance of the point from the origin
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(Vec3dp<T> const& v)
{
    return DualNum3dp<T>(bulk_nrm_sq(v), weight_nrm_sq(v));
}

// return geometric norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(Vec3dp<T> const& v)
{
    return DualNum3dp<T>(bulk_nrm(v), weight_nrm(v));
}

// provide the perpendicular distance of the line to the origin
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(BiVec3dp<T> const& B)
{
    return DualNum3dp<T>(bulk_nrm_sq(B), weight_nrm_sq(B));
}

// return geometric norm of bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(BiVec3dp<T> const& B)
{
    return DualNum3dp<T>(bulk_nrm(B), weight_nrm(B));
}

// provide the perpendicular distance of the plane to the origin
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm_sq(TriVec3dp<T> const& t)
{
    return DualNum3dp<T>(bulk_nrm_sq(t), weight_nrm_sq(t));
}

// return geometric norm of vector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr DualNum3dp<T> geom_nrm(TriVec3dp<T> const& t)
{
    return DualNum3dp<T>(bulk_nrm(t), weight_nrm(t));
}

////////////////////////////////////////////////////////////////////////////////
// unitization operations
////////////////////////////////////////////////////////////////////////////////

// return a vector unitized to v.w == 1.0  (implies weight_nrm(v) = 1.0)
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> unitize(Vec3dp<T> const& v)
{
    T n = v.w; // v.w == sign(v.w)*weight_nrm(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(n) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("Vector weight_nrm too small for unitization " +
                                 std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return Vec3dp<T>(v.x * inv, v.y * inv, v.z * inv, T(1.0));
}

// return a bivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> unitize(BiVec3dp<T> const& B)
{
    T n = weight_nrm(B);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(n) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("Bivector weight_norm too small for unitization " +
                                 std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * B;
}

// return a trivector unitized to weight_nrm == 1.0
template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> unitize(TriVec3dp<T> const& t)
{
    T n = weight_nrm(t);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(n) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("Trivector weight_norm too small for unitization " +
                                 std::to_string(n) + "\n");
    }
#endif
    T inv = T(1.0) / n; // for multiplication with inverse of norm
    return inv * t;
}


////////////////////////////////////////////////////////////////////////////////
// (right) bulk_dual (=complement operation applied to the bulk)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> bulk_dual(Scalar3dp<T> s)
{
    return PScalar3dp<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> bulk_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(v.x, v.y, v.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> bulk_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(-B.mx, -B.my, -B.mz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> bulk_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(T(0.0), T(0.0), T(0.0), -t.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> bulk_dual([[maybe_unused]] PScalar3dp<T>)
{
    return Scalar3dp<T>(0.0);
}


////////////////////////////////////////////////////////////////////////////////
// (right) weight_dual (=complement operation applied to the weight)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3dp<T> weight_dual([[maybe_unused]] Scalar3dp<T>)
{
    return PScalar3dp<T>(0.0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> weight_dual(Vec3dp<T> const& v)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), v.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> weight_dual(BiVec3dp<T> const& B)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> weight_dual(TriVec3dp<T> const& t)
{
    return Vec3dp<T>(-t.x, -t.y, -t.z, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> weight_dual(PScalar3dp<T> ps)
{
    return Scalar3dp<T>(T(ps));
}

////////////////////////////////////////////////////////////////////////////////
// convenient type aliases
////////////////////////////////////////////////////////////////////////////////

// Vector3d: 3d vector of projective algebra storing only components x, y, z explicitly.
// The w component is assumed to be w = 0.0
template <typename T>
    requires(std::floating_point<T>)
struct Vector3d : public Vec3d<T> {

    using Vec3d<T>::Vec3d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec3d<T>::x;
    using Vec3d<T>::y;
    using Vec3d<T>::z;

    Vector3d(Vec3d<T> const& v) : Vec3d<T>(v) {};
};

// Point3d: 3d point of projective algebra storing only components x, y, z explicitly.
// The w component is assumed to by w = 1.0
template <typename T>
    requires(std::floating_point<T>)
struct Point3d : public Vec3d<T> {

    using Vec3d<T>::Vec3d; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec3d<T>::x;
    using Vec3d<T>::y;
    using Vec3d<T>::z;

    Point3d(Vec3d<T> const& v) : Vec3d<T>(v) {};
    Point3d(T x, T y, T z) : Vec3d<T>(x, y, z) {};
};

// Point3dp: 3d point of projective algebra storing all four components x, y, z, w
//           explicitly
//
// a Point3dp is a Vec3dp, thus all operations defined for Vec3dp
// work directly for Point3dp - only deviations will be specified
template <typename T>
    requires(std::floating_point<T>)
struct Point3dp : public Vec3dp<T> {

    using Vec3dp<T>::Vec3dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using Vec3dp<T>::x;
    using Vec3dp<T>::y;
    using Vec3dp<T>::z;
    using Vec3dp<T>::w;

    Point3dp(Vec3dp<T> const& v) : Vec3dp<T>(v) {};
    Point3dp(Point3d<T> const& p) : Vec3dp<T>(p.x, p.y, p.z, T(1.0)) {};
    Point3dp(Vec3d<T> const& v) : Vec3dp<T>(v.x, v.y, v.z, T(1.0)) {};

    Point3dp& unitize()
    {
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (std::abs(z) < std::numeric_limits<T>::epsilon()) {
            throw std::runtime_error("z-component too small for unitization " +
                                     std::to_string(z) + "\n");
        }
#endif
        x /= w;
        y /= w;
        z /= w;
        w = T(1.0);
        return (*this);
    }
};

template <typename T>
    requires std::floating_point<T>
inline constexpr Point3dp<T> unitize(Point3dp<T> const& p)
{
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(p.w) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("w-component too small for unitization " +
                                 std::to_string(p.w) + "\n");
    }
#endif
    T inv = T(1.0) / p.w;
    return Point3dp<T>(p.x * inv, p.y * inv, p.z * inv, 1.0);
}


// Line3d: 3d line of projective algebra storing all six components
//          explicitly as components vx, vy, vz, mx, my, mz of a BiVec3dp
//
// a Line3d is a BiVec3dp, thus all operations defined for BiVec3dp
// work directly for Line3d - only deviations will be specified
template <typename T>
    requires(std::floating_point<T>)
struct Line3d : public BiVec3dp<T> {

    using BiVec3dp<T>::BiVec3dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using BiVec3dp<T>::vx;
    using BiVec3dp<T>::vy;
    using BiVec3dp<T>::vz;
    using BiVec3dp<T>::mx;
    using BiVec3dp<T>::my;
    using BiVec3dp<T>::mz;

    Line3d() = default;
    Line3d(T vx, T vy, T vz, T mx, T my, T mz) : BiVec3dp<T>(vx, vy, vz, mx, my, mz) {};
    Line3d(BiVec3dp<T> const& b) : BiVec3dp<T>(b) {};
    Line3d(Vec3d<T> const& dir, BiVec3d<T> const& mom) :
        // direction vector dir and moment bivector mom must match, i.e. be perpendicular
        // to each other (direction << moment) == 0
        BiVec3dp<T>(dir.x, dir.y, dir.z, mom.x, mom.y, mom.z) {};

    Line3d& unitize()
    {
        // unitization for a 3d bivector means std::sqrt(x^2 + y^2) = 1
        T wn = weight_nrm(*this);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
        if (wn < std::numeric_limits<T>::epsilon()) {
            throw std::runtime_error("bivector weight norm too small for unitization " +
                                     std::to_string(wn) + "\n");
        }
#endif
        T inv = T(1.0) / wn;
        vx *= inv;
        vy *= inv;
        vz *= inv;
        mx *= inv;
        my *= inv;
        mz *= inv;
        return (*this);
    }
};

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Line3d<T> unitize(Line3d<T> const& l)
{
    // unitization for a 3d bivector means std::sqrt(x^2 + y^2 + z^2) = 1
    // i.e. unitization of the direction vector of the line
    T wn = weight_nrm(l);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (wn < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector weight norm too small for unitization " +
                                 std::to_string(wn) + "\n");
    }
#endif
    T inv = T(1.0) / wn;
    return Line3d<T>(l.vx * inv, l.vy * inv, l.vz * inv, l.mx * inv, l.my * inv,
                     l.mz * inv);
}

// Plane3d: 3d plane of projective algebra (a trivector in the modeling 4d space)
template <typename T>
    requires(std::floating_point<T>)
struct Plane3d : public TriVec3dp<T> {

    using TriVec3dp<T>::TriVec3dp; // inherit base class ctors

    // When using partial template specialization, the compiler needs to know whether
    // we want to treat symbols as class members.
    using TriVec3dp<T>::x;
    using TriVec3dp<T>::y;
    using TriVec3dp<T>::z;
    using TriVec3dp<T>::w;

    Plane3d(TriVec3dp<T> const& t) : TriVec3dp<T>(t) {};
    Plane3d(Line3d<T> const& l, Point3d<T> const& p) : TriVec3dp<T>(join(l, p)) {};
    Plane3d(Point3d<T> const& p1, Point3d<T> const& p2, Point3d<T> const& p3) :
        TriVec3dp<T>(join(join(p1, p2), p3)) {};
};

////////////////////////////////////////////////////////////////////////////////
// attitude operations: att = rwdg( u, rcmpl(e4_3dp) ) = rwdg(u, horizon_3dp)
//
// (the attitude is the intersection of the object with the horizon)
// the result of att(object_with_grade_k) is an object with grade k-1
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> att(Vec3dp<T> const& v)
{
    return Scalar3dp<T>(v.w);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3dp<T> att(Point3dp<T> const& p)
{
    return Scalar3dp<T>(p.w);
}

// return the attitude (i.e. the direction vector) of the line
template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> att(BiVec3dp<T> const& B)
{
    return Vec3dp<T>(B.vx, B.vy, B.vz, T(0.0));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> att(Line3d<T> const& l)
{
    return Vec3dp<T>(l.vx, l.vy, l.vz, T(0.0));
}

// return the attitude (i.e. the normal vector) of the plane
template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> att(TriVec3dp<T> const& t)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), t.x, t.y, t.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3dp<T> att(Plane3d<T> const& l)
{
    return BiVec3dp<T>(T(0.0), T(0.0), T(0.0), l.x, l.y, l.z);
}


template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> att(PScalar3dp<T> ps)
{
    return TriVec3dp<T>(T(0.0), T(0.0), T(0.0), ps);
}

} // namespace hd::ga::pga