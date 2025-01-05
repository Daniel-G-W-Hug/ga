#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_type_3d.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec3d_E<T> M = c0 + (c1 e2^e3 + c2 e3^e1 + c3 e1^e2)
//
// with the term in brackets being the bivector c1 e2^e3 + c2 e3^e1 + c3 e1^e2
// that models a plane in 3D defined by the bivector coordinates (c1, c2, c3)
//
// This is the definition of a multivector in the even subalgebra of G<3,0,0>,
// i.e. it models only multivectors with even grades 0 and 2
// which are quaternions.
// This subalgebra is closed under addition and multiplication.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of a fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.

template <typename T> struct MVec4_t<T, mvec3d_e_tag> : public MVec4_t<T, default_tag> {

    using MVec4_t<T, default_tag>::MVec4_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec4_t(Scalar3d<T> s) : MVec4_t(T(s), T(0.0), T(0.0), T(0.0)) {}

    // assign a bivector part exclusively (other grades = 0)
    MVec4_t(BiVec3d<T> const& B) : MVec4_t(T(0.0), B.x, B.y, B.z) {}

    // assign scalar and bivector parts
    MVec4_t(Scalar3d<T> s, BiVec3d<T> const& B) : MVec4_t(T(s), B.x, B.y, B.z) {}
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec4_t<T, mvec3d_e_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of the even multivector
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> gr0(MVec3d_E<T> const& M)
{
    return Scalar3d<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> gr2(MVec3d_E<T> const& M)
{
    return BiVec3d<T>(M.c1, M.c2, M.c3);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar + bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator+(Scalar3d<T> s,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(s, B);
}

// bivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator+(BiVec3d<T> const& B,
                                                              Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(s, B);
}

// scalar + even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator+(Scalar3d<T> s,
                                                              MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(s + M.c0, M.c1, M.c2, M.c3);
}

// even grade multivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator+(MVec3d_E<T> const& M,
                                                              Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(s + M.c0, M.c1, M.c2, M.c3);
}

// bivector + even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator+(BiVec3d<T> const& B,
                                                              MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(M.c0, B.x + M.c1, B.y + M.c2, B.z + M.c3);
}

// even grade multivector + bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator+(MVec3d_E<T> const& M,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(M.c0, B.x + M.c1, B.y + M.c2, B.z + M.c3);
}

} // namespace hd::ga