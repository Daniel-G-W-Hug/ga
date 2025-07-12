#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_type2dp.hpp" // type aliases for scalars, vector, pseudoscalar, etc.

#include "ga_mvec2dp_e.hpp" // even grade multivector
#include "ga_mvec2dp_u.hpp" // uneven (or odd) grade multivector


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec8_t including its ctors and add specific ctors for MVec8_t<T, Tag>
// by using partial template specialization for the Tag=mvec2dp_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec8_t<T, mvec2dp_tag> : public MVec8_t<T, default_tag> {

    using MVec8_t<T, default_tag>::MVec8_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec8_t(Scalar2dp<T> s) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a vector part exclusively (other grades = 0)
    MVec8_t(Vec2dp<T> const& v) :
        MVec8_t(T(0.0), v.x, v.y, v.z, T(0.0), T(0.0), T(0.0), T(0.0))
    {
    }

    // assign a bivector part exclusively (other grades = 0)
    MVec8_t(BiVec2dp<T> const& B) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), B.x, B.y, B.z, T(0.0))
    {
    }

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec8_t(PScalar2dp<T> ps) :
        MVec8_t(T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign from a scalar and a bivector
    // (other grades = 0)
    MVec8_t(Scalar2dp<T> s, BiVec2dp<T> const& B) :
        MVec8_t(T(s), T(0.0), T(0.0), T(0.0), B.x, B.y, B.z, T(0.0))
    {
    }

    // assign from the even subalgebra
    MVec8_t(MVec2dp_E<T> const& M) :
        MVec8_t(M.c0, T(0.0), T(0.0), T(0.0), M.c1, M.c2, M.c3, T(0.0))
    {
    }

    // assign a vector and a pseudoscalar
    MVec8_t(Vec2dp<T> const& v, PScalar2dp<T> ps) :
        MVec8_t(T(0.0), v.x, v.y, v.z, T(0.0), T(0.0), T(0.0), T(ps))
    {
    }

    // assign from the odd subalgebra
    MVec8_t(MVec2dp_U<T> const& M) :
        MVec8_t(T(0.0), M.c0, M.c1, M.c2, T(0.0), T(0.0), T(0.0), M.c3)
    {
    }

    // assign a full multivector
    MVec8_t(Scalar2dp<T> s, Vec2dp<T> const& v, BiVec2dp<T> const& B, PScalar2dp<T> ps) :
        MVec8_t(T(s), v.x, v.y, v.z, B.x, B.y, B.z, T(ps))
    {
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec8_t<T, mvec2dp_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a 3d multivector
//
// grade 0: gr0() - scalar
// grade 1: gr1() - vector
// grade 2: gr2() - bivector
// grade 3: gr3() - trivector (= pseudoscalar in 3d)

template <typename T>
    requires(std::floating_point<T>)
constexpr Scalar2dp<T> gr0(MVec2dp<T> const& M)
{
    return Scalar2dp<T>(M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr Vec2dp<T> gr1(MVec2dp<T> const& M)
{
    return Vec2dp<T>(M.c1, M.c2, M.c3);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr BiVec2dp<T> gr2(MVec2dp<T> const& M)
{
    return BiVec2dp<T>(M.c4, M.c5, M.c6);
}

template <typename T>
    requires(std::floating_point<T>)
constexpr PScalar2dp<T> gr3(MVec2dp<T> const& M)
{
    return PScalar2dp<T>(M.c7);
}

// return the grades of the basic types

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] Scalar2dp<T>)
{
    return 0;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] Vec2dp<T> const&)
{
    return 1;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] BiVec2dp<T> const&)
{
    return 2;
}

template <typename T>
    requires(std::floating_point<T>)
constexpr size_t gr([[maybe_unused]] PScalar2dp<T>)
{
    return 3;
}


////////////////////////////////////////////////////////////////////////////////
// addition operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec2dp_e.hpp and mvec2dp_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar + vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(s, v, BiVec2dp<ctype>(0.0, 0.0, 0.0), PScalar2dp<ctype>(0.0));
}

// vector + scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(s, v, BiVec2dp<ctype>(0.0, 0.0, 0.0), PScalar2dp<ctype>(0.0));
}

// scalar + pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(s, Vec2dp<ctype>(0.0, 0.0, 0.0), BiVec2dp<ctype>(0.0, 0.0, 0.0),
                          ps);
}

// pseudoscalar + scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(s, Vec2dp<ctype>(0.0, 0.0, 0.0), BiVec2dp<ctype>(0.0, 0.0, 0.0),
                          ps);
}

// vector + bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(Vec2dp<T> const& v,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(Scalar2dp<ctype>(0.0), v, B, PScalar2dp<ctype>(0.0));
}

// bivector + vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(BiVec2dp<T> const& B,
                                                      Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(Scalar2dp<ctype>(0.0), v, B, PScalar2dp<ctype>(0.0));
}

// pseudoscalar + bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(PScalar2dp<T> ps,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(Scalar2dp<ctype>(0.0), Vec2dp<ctype>(0.0, 0.0, 0.0), B, ps);
}

// bivector + pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(BiVec2dp<T> const& B,
                                                      PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(Scalar2dp<ctype>(0.0), Vec2dp<ctype>(0.0, 0.0, 0.0), B, ps);
}

// mulitivector + scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(MVec2dp<T> const& M, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(gr0(M) + s, gr1(M), gr2(M), gr3(M));
}

// mulitivector + vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(MVec2dp<T> const& M,
                                                      Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(gr0(M), gr1(M) + v, gr2(M), gr3(M));
}

// mulitivector + bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(MVec2dp<T> const& M,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(gr0(M), gr1(M), gr2(M) + B, gr3(M));
}

// mulitivector + pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator+(MVec2dp<T> const& M,
                                                      PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) + ps);
}


////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars, pseudoscalar and vectors to multivectors
// (only remaining combinations not covered in mvec2dp_e.hpp and mvec2dp_u.hpp)
////////////////////////////////////////////////////////////////////////////////

// scalar - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(Scalar2dp<T> s, Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(s, -v, BiVec2dp<ctype>(0.0, 0.0, 0.0), PScalar2dp<ctype>(0.0));
}

// vector - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(Vec2dp<T> const& v, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(-s, v, BiVec2dp<ctype>(0.0, 0.0, 0.0), PScalar2dp<ctype>(0.0));
}

// scalar - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(Scalar2dp<T> s, PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(s, Vec2dp<ctype>(0.0, 0.0, 0.0), BiVec2dp<ctype>(0.0, 0.0, 0.0),
                          -ps);
}

// pseudoscalar - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(PScalar2dp<T> ps, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(-s, Vec2dp<ctype>(0.0, 0.0, 0.0),
                          BiVec2dp<ctype>(0.0, 0.0, 0.0), ps);
}

// vector - bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(Vec2dp<T> const& v,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(Scalar2dp<ctype>(0.0), v, -B, PScalar2dp<ctype>(0.0));
}

// bivector - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(BiVec2dp<T> const& B,
                                                      Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(Scalar2dp<ctype>(0.0), -v, B, PScalar2dp<ctype>(0.0));
}

// pseudoscalar - bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(PScalar2dp<T> ps,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(Scalar2dp<ctype>(0.0), Vec2dp<ctype>(0.0, 0.0, 0.0), -B, ps);
}

// bivector - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(BiVec2dp<T> const& B,
                                                      PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(Scalar2dp<ctype>(0.0), Vec2dp<ctype>(0.0, 0.0, 0.0), B, -ps);
}

// mulitivector - scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(MVec2dp<T> const& M, Scalar2dp<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(gr0(M) - s, gr1(M), gr2(M), gr3(M));
}

// mulitivector - vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(MVec2dp<T> const& M,
                                                      Vec2dp<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(gr0(M), gr1(M) - v, gr2(M), gr3(M));
}

// mulitivector - bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(MVec2dp<T> const& M,
                                                      BiVec2dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(gr0(M), gr1(M), gr2(M) - B, gr3(M));
}

// mulitivector - pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2dp<std::common_type_t<T, U>> operator-(MVec2dp<T> const& M,
                                                      PScalar2dp<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp<ctype>(gr0(M), gr1(M), gr2(M), gr3(M) - ps);
}

} // namespace hd::ga