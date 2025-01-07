#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_type_2dp.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec2dp_E<T> M = c0 + (c1 * e2^e3 + c2 * e3^e1 + c3 * e1^e2)
//
// with the term in brackets being the bivector
// c1 * e2^e3 + c2 * e3^e1 + c3 * e1^e2 that models a plane in 3D defined by
// the bivector coordinates (c1, c2, c3)
//
// This is the definition of a multivector in the even subalgebra of G<2,0,1>,
// i.e. it models only multivectors with even grades 0 and 2
// which are quaternions.
// This subalgebra is closed under addition and multiplication.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of a fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.

template <typename T> struct MVec4_t<T, mvec2dp_e_tag> : public MVec4_t<T, default_tag> {

    using MVec4_t<T, default_tag>::MVec4_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec4_t(Scalar2dp<T> s) : MVec4_t(T(s), T(0.0), T(0.0), T(0.0)) {}

    // assign a bivector part exclusively (other grades = 0)
    MVec4_t(BiVec2dp<T> const& B) : MVec4_t(T(0.0), B.x, B.y, B.z) {}

    // assign scalar and bivector parts
    MVec4_t(Scalar2dp<T> s, BiVec2dp<T> const& B) : MVec4_t(T(s), B.x, B.y, B.z) {}
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec4_t<T, mvec2dp_e_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of the even multivector
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2dp<T> gr0(MVec2dp_E<T> const& M)
{
    return Scalar2dp<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec2dp<T> gr2(MVec2dp_E<T> const& M)
{
    return BiVec2dp<T>(M.c1, M.c2, M.c3);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar + bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator+(Scalar2dp<T> s,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(s, B);
}

// bivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator+(BiVec2dp<T> const& B,
                                                               Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(s, B);
}

// scalar + even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator+(Scalar2dp<T> s,
                                                               MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(s + M.c0, M.c1, M.c2, M.c3);
}

// even grade multivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator+(MVec2dp_E<T> const& M,
                                                               Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(s + M.c0, M.c1, M.c2, M.c3);
}

// bivector + even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator+(BiVec2dp<T> const& B,
                                                               MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(M.c0, B.x + M.c1, B.y + M.c2, B.z + M.c3);
}

// even grade multivector + bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator+(MVec2dp_E<T> const& M,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(M.c0, B.x + M.c1, B.y + M.c2, B.z + M.c3);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar - bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator-(Scalar2dp<T> s,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(s, -B);
}

// bivector - scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator-(BiVec2dp<T> const& B,
                                                               Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(-s, B);
}

// scalar - even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator-(Scalar2dp<T> s,
                                                               MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(s - M.c0, -M.c1, -M.c2, -M.c3);
}

// even grade multivector - scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator-(MVec2dp_E<T> const& M,
                                                               Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(M.c0 - s, M.c1, M.c2, M.c3);
}

// bivector - even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator-(BiVec2dp<T> const& B,
                                                               MVec2dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(-M.c0, B.x - M.c1, B.y - M.c2, B.z - M.c3);
}

// even grade multivector - bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec2dp_E<std::common_type_t<T, U>> operator-(MVec2dp_E<T> const& M,
                                                               BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(M.c0, M.c1 - B.x, M.c2 - B.y, M.c3 - B.z);
}

} // namespace hd::ga