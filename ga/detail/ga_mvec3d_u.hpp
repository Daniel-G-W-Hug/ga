#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_type_3d.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec3d_U<T> M = (c0 * e1 + c1 * e2 + c2 * e3) + c4 e1^e2^e3
//
// with the term in brackets being the vector c0 * e1 + c1 * e2 + c2 * e3
// and the term c4 e1^e2^e3 being the trivector
//
// This is the definition of a multivector in the uneven subalgebra of G<3,0,0>,
// i.e. it models only multivectors with even grades 1 and 3
// This subalgebra is used to store intermediate results for rotations.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.

////////////////////////////////////////////////////////////////////////////////
// use MVec4_t including its ctors and add specific ctors for MVec4_t<T, Tag>
// by using partial template specialization for the Tag=mvec3d_u_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec4_t<T, mvec3d_u_tag> : public MVec4_t<T, default_tag> {

    using MVec4_t<T, default_tag>::MVec4_t; // inherit base class ctors

    // assign a vector part exclusively (other grades = 0)
    MVec4_t(Vec3d<T> const& v) : MVec4_t(v.x, v.y, v.z, T(0.0)) {}

    // assign a scalar part exclusively (other grades = 0)
    MVec4_t(PScalar3d<T> ps) : MVec4_t(T(0.0), T(0.0), T(0.0), T(ps)) {}

    // assign vector and pseudoscalar parts
    MVec4_t(Vec3d<T> const& v, PScalar3d<T> ps) : MVec4_t(v.x, v.y, v.z, T(ps)) {}
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec4_t<T, mvec3d_u_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of the uneven multivector
//
// grade 1: gr1() - vector
// grade 2: gr3() - pseudoscalar

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> gr1(MVec3d_U<T> const& M)
{
    return Vec3d<T>(M.c0, M.c1, M.c2);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> gr3(MVec3d_U<T> const& M)
{
    return PScalar3d<T>(M.c3);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine vectors and trivectors to uneven grade multivectors
////////////////////////////////////////////////////////////////////////////////

// vector + pseudoscalar => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator+(Vec3d<T> const& v,
                                                              PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(v, ps);
}

// pseudoscalar + vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator+(PScalar3d<T> ps,
                                                              Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(v, ps);
}

// uneven grade mulivector + pseudoscalar => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator+(MVec3d_U<T> const& M,
                                                              PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(M.c0, M.c1, M.c2, M.c3 + ps);
}

// pseudoscalar + uneven grade vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator+(PScalar3d<T> ps,
                                                              MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(M.c0, M.c1, M.c2, M.c3 + ps);
}

// uneven grade mulivector + vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator+(MVec3d_U<T> const& M,
                                                              Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(M.c0 + v.x, M.c1 + v.y, M.c2 + v.z, M.c3);
}

// vector + uneven grade vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator+(Vec3d<T> const& v,
                                                              MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(M.c0 + v.x, M.c1 + v.y, M.c2 + v.z, M.c3);
}

} // namespace hd::ga