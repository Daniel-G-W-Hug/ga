#pragma once

// author: Daniel Hug, 2024 & 2025

#include "type_t/ga_type_3dp.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dp_u_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec3dp_u_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign the vector part exclusively (other grades = 0)
    MVec8_t(Vec3dp<T> const& v) :
        MVec8_t(v.x, v.y, v.z, v.w, T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign the trivector part exclusively (other grades = 0)
    MVec8_t(TriVec3dp<T> const& t) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), t.x, t.y, t.z, t.w)
    {
    }

    // assign from a vector and a trivector directly
    MVec8_t(Vec3dp<T> const& v, TriVec3dp<T> const& t) :
        MVec8_t(v.x, v.y, v.z, v.w, t.x, t.y, t.z, t.w)
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec8_t<T, mvec3dp_u_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a 3dp uneven multivector
//
// grade 1: gr1() - vector
// grade 3: gr3() - trivector

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3dp<T> gr1(MVec3dp_U<T> const& M)
{
    return Vec3dp<T>(M.c0, M.c1, M.c2, M.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr TriVec3dp<T> gr3(MVec3dp_U<T> const& M)
{
    return TriVec3dp<T>(M.c4, M.c5, M.c6, M.c7);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine vectors and trivectors to uneven grade multivectors
////////////////////////////////////////////////////////////////////////////////

// vector + trivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(Vec3dp<T> const& v,
                                                               TriVec3dp<T> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(v, t);
}

// trivector + vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator+(TriVec3dp<T> const& t,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(v, t);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine vectors and trivectors to uneven grade multivectors
////////////////////////////////////////////////////////////////////////////////

// vector + trivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(Vec3dp<T> const& v,
                                                               TriVec3dp<T> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(v, -t);
}

// trivector + vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3dp_U<std::common_type_t<T, U>> operator-(TriVec3dp<T> const& t,
                                                               Vec3dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_U<ctype>(-v, t);
}

} // namespace hd::ga