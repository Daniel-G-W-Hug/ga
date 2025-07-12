#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type3dp.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec3dp_e_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec3dp_e_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec8_t(Scalar3dp<T> s) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec8_t(BiVec3dp<T> const& B) :
        MVec8_t(T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec8_t(PScalar3dp<T> ps) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign the scalar and the bivector part only, assuming the pseudoscalar part as 0.0
    MVec8_t(Scalar3dp<T> s, BiVec3dp<T> const& B) :
        MVec8_t(T(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(0.0))
    {
    }

    // assign the bivector and the pseudoscalar part only, assuming the scalar part as 0.0
    MVec8_t(BiVec3dp<T> const& B, PScalar3dp<T> ps) :
        MVec8_t(T(0.0), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(ps))
    {
    }

    // assign the scalar and the pseudoscalar part only, assuming the bivector part as 0.0
    MVec8_t(Scalar3dp<T> s, PScalar3dp<T> ps) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign all three parts, the scalar, the bivector and the pseudoscalar explicitly
    MVec8_t(Scalar3dp<T> s, BiVec3dp<T> const& B, PScalar3dp<T> ps) :
        MVec8_t(T(s), B.vx, B.vy, B.vz, B.mx, B.my, B.mz, T(ps))
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec8_t<T, mvec3dp_e_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of an even 3dp multivector
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector
// grade 3: gr4() - quadvector (= pseudoscalar in 3dp)

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar3dp<T> gr0(MVec3dp_E<T> const& M)
{
    return Scalar3dp<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec3dp<T> gr2(MVec3dp_E<T> const& M)
{
    return BiVec3dp<T>(M.c1, M.c2, M.c3, M.c4, M.c5, M.c6);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar3dp<T> gr4(MVec3dp_E<T> const& M)
{
    return PScalar3dp<T>(M.c7);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar + bivector => even grade multivector (in 3dp imply a zero pseudoscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(Scalar3dp<T> s,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s, B);
}

// bivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(BiVec3dp<T> const& B,
                                                        Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s, B);
}

// pseudoscalar + bivector => even grade multivector (in 3dp imply a zero scalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(PScalar3dp<T> ps,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(B, ps);
}

// bivector + pseudoscalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(BiVec3dp<T> const& B,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(B, ps);
}

// pseudoscalar + scalar => even grade multivector (in 3dp imply a zero bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s, ps);
}

// scalar + pseudoscalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s, ps);
}

// scalar + even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(Scalar3dp<T> s,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s + M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7);
}

// even grade multivector + scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(MVec3dp_E<T> const& M,
                                                        Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s + M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7);
}

// bivector + even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(BiVec3dp<T> const& B,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(M.c0, B.vx + M.c1, B.vy + M.c2, B.vz + M.c3, B.mx + M.c4,
                            B.my + M.c5, B.mz + M.c6, M.c7);
}

// even grade multivector + bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(MVec3dp_E<T> const& M,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(M.c0, B.vx + M.c1, B.vy + M.c2, B.vz + M.c3, B.mx + M.c4,
                            B.my + M.c5, B.mz + M.c6, M.c7);
}

// pseudoscalar + even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(PScalar3dp<T> ps,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, ps + M.c7);
}

// even grade multivector + pseudoscalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator+(MVec3dp_E<T> const& M,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, ps + M.c7);
}


////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars and bivectors to even grade multivectors
////////////////////////////////////////////////////////////////////////////////

// scalar - bivector => even grade multivector (in 3dp imply a zero pseudoscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(Scalar3dp<T> s,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s, -B);
}

// bivector - scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(BiVec3dp<T> const& B,
                                                        Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(-s, B);
}

// pseudoscalar - bivector => even grade multivector (in 3dp imply a zero scalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(PScalar3dp<T> ps,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(-B, ps);
}

// bivector - pseudoscalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(BiVec3dp<T> const& B,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(B, -ps);
}

// pseudoscalar - scalar => even grade multivector (in 3dp imply a zero bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(PScalar3dp<T> ps, Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(-s, ps);
}

// scalar - pseudoscalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(Scalar3dp<T> s, PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s, -ps);
}

// scalar - even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(Scalar3dp<T> s,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(s - M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, -M.c7);
}

// even grade multivector - scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(MVec3dp_E<T> const& M,
                                                        Scalar3dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(-s + M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, M.c7);
}

// bivector - even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(BiVec3dp<T> const& B,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(-M.c0, B.vx - M.c1, B.vy - M.c2, B.vz - M.c3, B.mx - M.c4,
                            B.my - M.c5, B.mz - M.c6, -M.c7);
}

// even grade multivector - bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(MVec3dp_E<T> const& M,
                                                        BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(M.c0, -B.vx + M.c1, -B.vy + M.c2, -B.vz + M.c3, -B.mx + M.c4,
                            -B.my + M.c5, -B.mz + M.c6, M.c7);
}

// pseudoscalar - even grade mulivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(PScalar3dp<T> ps,
                                                        MVec3dp_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(-M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, ps - M.c7);
}

// even grade multivector - pseudoscalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec3dp_E<std::common_type_t<T, U>> operator-(MVec3dp_E<T> const& M,
                                                        PScalar3dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3dp_E<ctype>(M.c0, M.c1, M.c2, M.c3, M.c4, M.c5, M.c6, -ps + M.c7);
}

} // namespace hd::ga