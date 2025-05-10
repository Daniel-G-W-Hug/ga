#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <algorithm> // std::max, std::clamp
#include <cmath>     // std::abs
#include <concepts>  // std::floating_point<T>
#include <iostream>  // std::cout, std::ostream
#include <limits>    // std::numeric_limits
#include <numbers>   // math constants like pi
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

#include "detail/ga_mvec3d.hpp" // inclusion of multivector imports all component types

namespace hd::ga::ega {

////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for k = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> gr_inv(Scalar3d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> gr_inv(Vec3d<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> gr_inv(BiVec3d<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> gr_inv(PScalar3d<T> ps)
{
    // grade 3: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> gr_inv(MVec3d_E<T> const& M)
{
    // grade 0 and 2: no sign change
    return M;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> gr_inv(MVec3d_U<T> const& M)
{
    // grade 1 and 3: sign reversal
    return -M;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> gr_inv(MVec3d<T> const& M)
{
    // grade 0 and 2: no sign change
    // grade 1 and 3: sign reversal
    return MVec3d<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)), gr_inv(gr3(M)));
}

////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> rev(Scalar3d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> rev(Vec3d<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> rev(BiVec3d<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> rev(PScalar3d<T> ps)
{
    // grade 3: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> rev(MVec3d_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec3d_E<T>(rev(gr0(M)), rev(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> rev(MVec3d_U<T> const& M)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec3d_U<T>(rev(gr1(M)), rev(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> rev(MVec3d<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    return MVec3d<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)), rev(gr3(M)));
}

////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation:
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> conj(Scalar3d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> conj(Vec3d<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> conj(BiVec3d<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> conj(PScalar3d<T> ps)
{
    // grade 3: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> conj(MVec3d_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec3d_E<T>(conj(gr0(M)), conj(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> conj(MVec3d_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec3d_U<T>(conj(gr1(M)), conj(gr3(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> conj(MVec3d<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    return MVec3d<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)), conj(gr3(M)));
}


////////////////////////////////////////////////////////////////////////////////
// dot product (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1, e3^2 = +1
//
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(MVec3d<T> const& A,
                                                        MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 +
                           A.c4 * B.c4 + A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(PScalar3d<T> ps1,
                                                        PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps1) * ctype(ps2));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(BiVec3d<T> const& B1,
                                                        BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(B1.x * B2.x + B1.y * B2.y + B1.z * B2.z);
}

// coordinate free definition: dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(Vec3d<T> const& v1,
                                                        Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(Scalar3d<T> s1, Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// wedge product (= outer product)
////////////////////////////////////////////////////////////////////////////////

// wedge product extended to fully populated multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& A,
                                                      MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype c3 = A.c0 * B.c3 + A.c3 * B.c0;
    ctype c4 = A.c0 * B.c4 + A.c4 * B.c0 + A.c2 * B.c3 - A.c3 * B.c2;
    ctype c5 = A.c0 * B.c5 + A.c5 * B.c0 + A.c3 * B.c1 - A.c1 * B.c3;
    ctype c6 = A.c0 * B.c6 + A.c6 * B.c0 + A.c1 * B.c2 - A.c2 * B.c1;
    ctype c7 = A.c0 * B.c7 + A.c7 * B.c0 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 +
               A.c4 * B.c1 + A.c5 * B.c2 + A.c6 * B.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// wedge product between a multivector M and a bivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& M,
                                                      BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = 0.0;
    ctype c1 = 0.0;
    ctype c2 = 0.0;
    ctype c3 = 0.0;
    ctype c4 = M.c0 * B.x;
    ctype c5 = M.c0 * B.y;
    ctype c6 = M.c0 * B.z;
    ctype c7 = M.c1 * B.x + M.c2 * B.y + M.c3 * B.z;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// wedge product between a bivector B and a multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& B,
                                                      MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = 0;
    ctype c1 = 0;
    ctype c2 = 0;
    ctype c3 = 0;
    ctype c4 = B.x * M.c0;
    ctype c5 = B.y * M.c0;
    ctype c6 = B.z * M.c0;
    ctype c7 = B.x * M.c1 + B.y * M.c2 + B.z * M.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// wedge product between a multivector M and a vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<T> const& M,
                                                      Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = 0.0;
    ctype c1 = M.c0 * v.x;
    ctype c2 = M.c0 * v.y;
    ctype c3 = M.c0 * v.z;
    ctype c4 = M.c2 * v.z - M.c3 * v.y;
    ctype c5 = -M.c1 * v.z + M.c3 * v.x;
    ctype c6 = M.c1 * v.y - M.c2 * v.x;
    ctype c7 = M.c4 * v.x + M.c5 * v.y + M.c6 * v.z;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// wedge product between a vector v and a multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v,
                                                      MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = 0;
    ctype c1 = v.x * M.c0;
    ctype c2 = v.y * M.c0;
    ctype c3 = v.z * M.c0;
    ctype c4 = v.y * M.c3 - v.z * M.c2;
    ctype c5 = -v.x * M.c3 + v.z * M.c1;
    ctype c6 = v.x * M.c2 - v.y * M.c1;
    ctype c7 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// wedge product between multivector M and scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> wdg(MVec3d<U> const& M, Scalar3d<T> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// wedge product between a scalar s and a multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// wedge product between two pseudoscalars
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3d<T>,
                                                        [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// wedge product between a pseudoscalar and a bivector
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
wdg([[maybe_unused]] PScalar3d<T>, [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// wedge product between a bivector and a pseudoscalar
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
wdg([[maybe_unused]] BiVec3d<T> const&, [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// wedge product between a pseudoscalar and a vector
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar3d<T>,
                                                        [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// wedge product between a vector and a pseudoscalar
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec3d<T> const&,
                                                        [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> wdg(PScalar3d<T> ps, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) * ctype(s));
}

// wedge product with a pscalar (returns a scaled pscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(s) * ctype(ps));
}

// wedge product between two bivectors
// => returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
wdg([[maybe_unused]] BiVec3d<T> const&, [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// wedge product between a bivector B and a vector v
// => returns a trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& B,
                                                         Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(B.x * v.x + B.y * v.y + B.z * v.z);
}

// wedge product between a vector v and a bivector B
// => returns a trivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v,
                                                         BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(v.x * B.x + v.y * B.y + v.z * B.z);
}

// wedge product with one scalar (returns a scaled bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> wdg(BiVec3d<T> const& B, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// wedge product with one scalar (returns a scaled bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

// wedge product between two vectors (returns a bivector in 3d)
// coordinate-free definition: wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v1,
                                                       Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z,
                          v1.x * v2.y - v1.y * v2.x);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// wedge product between two scalars (returns a scalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s1, Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive wedge product (= outer product for complements)
// as defined by E. Lengyel in "Projective geometric algebra illuminated"
// independent of the geometric product, just depending on the outer product (wdg)
// as well as the complement and thus the pseudoscalar of the space
// (in this definition is does NOT connect directly to the geometric product,
// but to the outer product exclusively)
//
//        rwdg(ul, ur) = cmpl(wdg(cmpl(ul),cmpl(ur))) = cmpl(cmpl(ul) ^ cmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

// regressive wedge product extended to a fully populated multivectors
// => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d<T> const& A,
                                                       MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;

    ctype c0 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
               A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    ctype c1 = A.c1 * B.c7 + A.c5 * B.c6 - A.c6 * B.c5 + A.c7 * B.c1;
    ctype c2 = A.c2 * B.c7 - A.c4 * B.c6 + A.c6 * B.c4 + A.c7 * B.c2;
    ctype c3 = A.c3 * B.c7 + A.c4 * B.c5 - A.c5 * B.c4 + A.c7 * B.c3;
    ctype c4 = A.c4 * B.c7 + A.c7 * B.c4;
    ctype c5 = A.c5 * B.c7 + A.c7 * B.c5;
    ctype c6 = A.c6 * B.c7 + A.c7 * B.c6;
    ctype c7 = A.c7 * B.c7;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// regressive wedge product between to pseudoscalars
// => returns a scaled pseudoscalar (=trivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps1,
                                                          PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps1) * ctype(ps2));
}

// regressive wedge product between a pseudoscalar ps (=trivector) and a bivector B
// => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps,
                                                        BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * BiVec3d<ctype>(B.x, B.y, B.z);
}

// regressive wedge product between a bivector B and a pseudoscalar ps (=trivector)
// => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> rwdg(BiVec3d<T> const& B,
                                                        PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(B.x, B.y, B.z) * ctype(ps);
}

// regressive wedge product between a pseudoscalar ps (=trivector) and a vector v
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec3d<ctype>(v.x, v.y, v.z);
}

// regressive wedge product between a vector v and a pseudoscalar ps (=trivector)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> rwdg(Vec3d<T> const& v, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v.x, v.y, v.z) * ctype(ps);
}


// regressive wedge product between a pseudoscalar ps (=trivector) and a scalar s
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> rwdg(PScalar3d<T> ps, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps) * ctype(s));
}

// regressive wedge product between a scalar s and a pseudoscalar ps (=trivector)
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> rwdg(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s) * ctype(ps));
}

// regressive wedge product between a bivector B1 and a bivector B2
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> rwdg(BiVec3d<T> const& B1,
                                                      BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(B1.y * B2.z - B1.z * B2.y, B1.z * B2.x - B1.x * B2.z,
                        B1.x * B2.y - B1.y * B2.x);
}

// regressive wedge product between a bivector B and a vector v
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> rwdg(BiVec3d<T> const& B,
                                                         Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(B.x * v.x + B.y * v.y + B.z * v.z);
}

// regressive wedge product between a vector v and a bivector B
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> rwdg(Vec3d<T> const& v,
                                                         BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v.x * B.x + v.y * B.y + v.z * B.z);
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The resulting object is a lies in B and is perpendicular to A
//
// L. Dorst: The contraction A << B of an a-blade A onto a b-blade B is a sub-blade
// of B of grade b-a which is perpendicular to A, and linear in both arguments
////////////////////////////////////////////////////////////////////////////////

// left contraction (A << B) extended to full multivectors
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d<T> const& A,
                                                             MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c4 * B.c4 +
               A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7;
    ctype c1 = A.c0 * B.c1 + A.c2 * B.c6 - A.c3 * B.c5 + A.c4 * B.c7;
    ctype c2 = A.c0 * B.c2 - A.c1 * B.c6 + A.c3 * B.c4 + A.c5 * B.c7;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c5 - A.c2 * B.c4 + A.c6 * B.c7;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c7;
    ctype c5 = A.c0 * B.c5 + A.c2 * B.c7;
    ctype c6 = A.c0 * B.c6 + A.c3 * B.c7;
    ctype c7 = A.c0 * B.c7;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// left contraction (B << M) of a bivector B onto a multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B,
                                                             MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = B.x * M.c4 + B.y * M.c5 + B.z * M.c6;
    ctype c1 = B.x * M.c7;
    ctype c2 = B.y * M.c7;
    ctype c3 = B.z * M.c7;
    ctype c4 = 0.0;
    ctype c5 = 0.0;
    ctype c6 = 0.0;
    ctype c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// left contraction (v << M) of a vector v onto a multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v,
                                                             MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3;
    ctype c1 = v.y * M.c6 - v.z * M.c5;
    ctype c2 = -v.x * M.c6 + v.z * M.c4;
    ctype c3 = v.x * M.c5 - v.y * M.c4;
    ctype c4 = v.x * M.c7;
    ctype c5 = v.y * M.c7;
    ctype c6 = v.z * M.c7;
    ctype c7 = 0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// left contraction (s << M) of a scalar s with a mulivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                             MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// left contraction (s << B) of a scalar s with an even grade mulivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                               MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// left contraction (s << B) of a scalar s with an uneven grade mulivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                               MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// left contraction (ps1 << ps2) - pseudoscalar ps1 contracted onto pseudoscalar ps2
// returns a scalar (the product of the pseudoscalars)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator<<(PScalar3d<T> ps1,
                                                               PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps1) * ctype(ps2));
}

// left contraction - pseudoscalar contracted onto bivector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3d<T>, [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// left contraction of a bivector B onto a pseudoscalar ps (=trivector)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B,
                                                            PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(B.x, B.y, B.z) * ctype(ps);
}

// left contraction -  pseudoscalar contracted onto vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3d<T>, [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// left contraction of a vector v onto a pseudoscalar ps (=trivector)
// => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v,
                                                              PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(v.x, v.y, v.z) * ctype(ps);
}

// left contraction - pseudoscalar contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] PScalar3d<T>, [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// left contraction (s << ps) - scalar s contracted onto pseudoscalar ps
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                                PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(s) * ctype(ps));
}

// left contraction (B1 << B2) - bivector B1 contracted onto B2
// (=identical to scalar product dot(B1,B2) of two bivectors; identical to (B1 >> B2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B1,
                                                               BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(B1.x * B2.x + B1.y * B2.y + B1.z * B2.z);
}

// left contraction - bivector contracted onto vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3d<T> const&, [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// left contraction of a vector v onto a bivector B ("v contracted onto B")
// => returns a vector
// (originally implemented as dot-product)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v,
                                                            BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v.y * B.z - v.z * B.y, -v.x * B.z + v.z * B.x,
                        v.x * B.y - v.y * B.x);
}

// left contraction - bivector contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] BiVec3d<T> const&, [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// left contraction (s << B) - scalar s contracted onto bivector B
// identical with scalar multiplication of the bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

// left contraction (v1 << v2) - vector v1 contracted onto v2
// (=identical to scalar product dot(v1,v2) of two vectors; identical to (v1 >> v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator<<(Vec3d<T> const& v1,
                                                               Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// left contraction - vector contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator<<([[maybe_unused]] Vec3d<T> const&, [[maybe_unused]] Scalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// left contraction (s << v) - scalar s contracted onto vector v
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                            Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// left contraction (s1 << s2) - scalar s1 contracted onto scalar s2
// returns a scalar (the product of the scalars)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s1,
                                                               Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B: "A contracted by B"
//
// The resulting object lies in A and is perpendicular to B
////////////////////////////////////////////////////////////////////////////////

// right contraction (A >> B) extended to a full multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& A,
                                                             MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c4 * B.c4 +
               A.c5 * B.c5 + A.c6 * B.c6 + A.c7 * B.c7;
    ctype c1 = A.c1 * B.c0 + A.c5 * B.c3 - A.c6 * B.c2 + A.c7 * B.c4;
    ctype c2 = A.c2 * B.c0 - A.c4 * B.c3 + A.c6 * B.c1 + A.c7 * B.c5;
    ctype c3 = A.c3 * B.c0 + A.c4 * B.c2 - A.c5 * B.c1 + A.c7 * B.c6;
    ctype c4 = A.c4 * B.c0 + A.c7 * B.c1;
    ctype c5 = A.c5 * B.c0 + A.c7 * B.c2;
    ctype c6 = A.c6 * B.c0 + A.c7 * B.c3;
    ctype c7 = A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// right contraction (M >> B) of a multivector M contracted by a bivector B
// right contraction (A >> B) extended to a full multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& M,
                                                             BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c4 * B.x + M.c5 * B.y + M.c6 * B.z;
    ctype c1 = M.c7 * B.x;
    ctype c2 = M.c7 * B.y;
    ctype c3 = M.c7 * B.z;
    ctype c4 = 0.0;
    ctype c5 = 0.0;
    ctype c6 = 0.0;
    ctype c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// right contraction (M >> v) of a multivector M contracted by a vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& M,
                                                             Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z;
    ctype c1 = M.c5 * v.z - M.c6 * v.y;
    ctype c2 = -M.c4 * v.z + M.c6 * v.x;
    ctype c3 = M.c4 * v.y - M.c5 * v.x;
    ctype c4 = M.c7 * v.x;
    ctype c5 = M.c7 * v.y;
    ctype c6 = M.c7 * v.z;
    ctype c7 = 0.0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// right contraction (M >> s) of a multivector by a scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& M,
                                                             Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// right contraction (A >> s) of even grade multivector A by scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator>>(MVec3d_E<T> const& M,
                                                               Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// right contraction (A >> s) of uneven grade multivector A by scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator>>(MVec3d_U<T> const& M,
                                                               Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// right contraction - pseudoscalar contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps1,
                                                               PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(ps1) * ctype(ps2));
}

// right contraction of pseudoscalar ps (=trivector) by bivector B "ps contracted by B"
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps,
                                                            BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec3d<ctype>(B.x, B.y, B.z);
}

// right contraction - bivector contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] BiVec3d<T> const&, [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// right contraction of pseudoscalar ps (=trivector) by vector v (" ps contracted by v")
// => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps,
                                                              Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * BiVec3d<ctype>(v.x, v.y, v.z);
}

// right contraction -  vector contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3d<T> const&, [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// right contraction (ps >> s) - pseudoscalar ps contracted by scalar s
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator>>(PScalar3d<T> ps,
                                                                Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) * ctype(s));
}

// right contraction - scalar contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3d<T>, [[maybe_unused]] PScalar3d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// right contraction (B1 >> B2) - bivector B1 contracted by bivector B2
// (=identical to scalar product dot(B1,B2) of two bivectors; identical to (B1 >> B2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B1,
                                                               BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(B1.x * B2.x + B1.y * B2.y + B1.z * B2.z);
}

// right contraction of bivector B by vector v ("B contracted by v")
// => returns a vector
// (originally implemented as dot-product)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B,
                                                            Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(B.y * v.z - B.z * v.y, -B.x * v.z + B.z * v.x,
                        B.x * v.y - B.y * v.x);
}

// right contraction - vector contracted by a bivector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Vec3d<T> const&, [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// right contraction (B >> s) - bivector B contracted by scalar s
// identical with scalar multiplication of the bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B,
                                                              Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// right contraction - scalar contracted by a bivector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3d<T>, [[maybe_unused]] BiVec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// right contraction (v1 >> v2) - vector v1 contracted by vector v2
// (=identical to scalar product dot(v1,v2) of two vectors; identical to (v1 << v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator>>(Vec3d<T> const& v1,
                                                               Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// right contraction (v >> s) - vector v contracted by scalar s
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator>>(Vec3d<T> const& v,
                                                            Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// right contraction - scalar contracted by a vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>>
operator>>([[maybe_unused]] Scalar3d<T>, [[maybe_unused]] Vec3d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(0.0);
}

// right contraction (s1 >> s2) - scalar s1 contracted by scalar s2
// returns a scalar (the product of the scalars)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator>>(Scalar3d<T> s1,
                                                               Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// classical cross product
////////////////////////////////////////////////////////////////////////////////

// cross-product between two vectors (only defined in 3d)
//  => returns a 3d vector
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> cross(Vec3d<T> const& v1,
                                                       Vec3d<U> const& v2)
{
    return Vec3d<std::common_type_t<T, U>>(
        v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}
//     double cross product identity:
//     a x (b x c) = -dot(a, wdg(b,c)) = -a.(b^c)


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> cmt(MVec3d<T> const& A,
                                                      MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(0.0, -A.c2 * B.c6 + A.c3 * B.c5 - A.c5 * B.c3 + A.c6 * B.c2,
                         A.c1 * B.c6 - A.c3 * B.c4 + A.c4 * B.c3 - A.c6 * B.c1,
                         -A.c1 * B.c5 + A.c2 * B.c4 - A.c4 * B.c2 + A.c5 * B.c1,
                         A.c2 * B.c3 - A.c3 * B.c2 - A.c5 * B.c6 + A.c6 * B.c5,
                         -A.c1 * B.c3 + A.c3 * B.c1 + A.c4 * B.c6 - A.c6 * B.c4,
                         A.c1 * B.c2 - A.c2 * B.c1 - A.c4 * B.c5 + A.c5 * B.c4, 0.0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> cmt(BiVec3d<T> const& B1,
                                                       BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(-B1.y * B2.z + B1.z * B2.y, B1.x * B2.z - B1.z * B2.x,
                          -B1.x * B2.y + B1.y * B2.x);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> cmt(BiVec3d<T> const& B,
                                                     Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(-B.y * v.z + B.z * v.y, B.x * v.z - B.z * v.x,
                        -B.x * v.y + B.y * v.x);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> cmt(Vec3d<T> const& v,
                                                     BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(-v.y * B.z + v.z * B.y, v.x * B.z - v.z * B.x,
                        -v.x * B.y + v.y * B.x);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> cmt(Vec3d<T> const& v1,
                                                       Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v1.y * v2.z - v1.z * v2.y, -v1.x * v2.z + v1.z * v2.x,
                        v1.x * v2.y - v1.y * v2.x);
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product A*B for fully populated 3d multivector
// Expensive! - Don't use if you don't have to! (64x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors:
// a * b = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector*vector = scalar + bivector)
//       = (b >> a) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector*vector = scalar + bivector)
//       = (a << b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector*vector = scalar + bivector)
// A * b = (a << A) + wdg(A,b) = gr1(Ab) + gr3(Ab)  (bivector*vector = vector + trivector)
// a * B = (B >> a) + wdg(a,B) = gr1(aB) + gr3(aB)  (vector*bivector = vector + trivector)
//
// multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                            MVec3d<U> const& B)
{
    // geometric product of two fully populated 3d multivectors
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c4 -
               A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 - A.c4 * B.c7 -
               A.c5 * B.c3 + A.c6 * B.c2 - A.c7 * B.c4;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4 + A.c4 * B.c3 -
               A.c5 * B.c7 - A.c6 * B.c1 - A.c7 * B.c5;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0 - A.c4 * B.c2 +
               A.c5 * B.c1 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c7 + A.c2 * B.c3 - A.c3 * B.c2 + A.c4 * B.c0 -
               A.c5 * B.c6 + A.c6 * B.c5 + A.c7 * B.c1;
    ctype c5 = A.c0 * B.c5 - A.c1 * B.c3 + A.c2 * B.c7 + A.c3 * B.c1 + A.c4 * B.c6 +
               A.c5 * B.c0 - A.c6 * B.c4 + A.c7 * B.c2;
    ctype c6 = A.c0 * B.c6 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c7 - A.c4 * B.c5 +
               A.c5 * B.c4 + A.c6 * B.c0 + A.c7 * B.c3;
    ctype c7 = A.c0 * B.c7 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 + A.c4 * B.c1 +
               A.c5 * B.c2 + A.c6 * B.c3 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of a multivector with an even grade multivector
// multivector * even grade multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                            MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c4 * B.c1 - A.c5 * B.c2 - A.c6 * B.c3;
    ctype c1 = A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2 - A.c7 * B.c1;
    ctype c2 = A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1 - A.c7 * B.c2;
    ctype c3 = -A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0 - A.c7 * B.c3;
    ctype c4 = A.c0 * B.c1 + A.c4 * B.c0 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype c5 = A.c0 * B.c2 + A.c4 * B.c3 + A.c5 * B.c0 - A.c6 * B.c1;
    ctype c6 = A.c0 * B.c3 - A.c4 * B.c2 + A.c5 * B.c1 + A.c6 * B.c0;
    ctype c7 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product of an even grade multivector with a multivector
// even grade multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                            MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 - A.c1 * B.c4 - A.c2 * B.c5 - A.c3 * B.c6;
    ctype c1 = A.c0 * B.c1 - A.c1 * B.c7 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c3 - A.c2 * B.c7 - A.c3 * B.c1;
    ctype c3 = A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 - A.c3 * B.c7;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5;
    ctype c5 = A.c0 * B.c5 + A.c1 * B.c6 + A.c2 * B.c0 - A.c3 * B.c4;
    ctype c6 = A.c0 * B.c6 - A.c1 * B.c5 + A.c2 * B.c4 + A.c3 * B.c0;
    ctype c7 = A.c0 * B.c7 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// multivector * uneven multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                            MVec3d_U<U> const& B)
{
    // geometric product of two fully populated 3d multivectors
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2 - A.c7 * B.c3;
    ctype c1 = A.c0 * B.c0 - A.c4 * B.c3 - A.c5 * B.c2 + A.c6 * B.c1;
    ctype c2 = A.c0 * B.c1 + A.c4 * B.c2 - A.c5 * B.c3 - A.c6 * B.c0;
    ctype c3 = A.c0 * B.c2 - A.c4 * B.c1 + A.c5 * B.c0 - A.c6 * B.c3;
    ctype c4 = A.c1 * B.c3 + A.c2 * B.c2 - A.c3 * B.c1 + A.c7 * B.c0;
    ctype c5 = -A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0 + A.c7 * B.c1;
    ctype c6 = A.c1 * B.c1 - A.c2 * B.c0 + A.c3 * B.c3 + A.c7 * B.c2;
    ctype c7 = A.c0 * B.c3 + A.c4 * B.c0 + A.c5 * B.c1 + A.c6 * B.c2;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// uneven multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                            MVec3d<U> const& B)
{
    // geometric product of two fully populated 3d multivectors
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 - A.c3 * B.c7;
    ctype c1 = A.c0 * B.c0 - A.c1 * B.c6 + A.c2 * B.c5 - A.c3 * B.c4;
    ctype c2 = A.c0 * B.c6 + A.c1 * B.c0 - A.c2 * B.c4 - A.c3 * B.c5;
    ctype c3 = -A.c0 * B.c5 + A.c1 * B.c4 + A.c2 * B.c0 - A.c3 * B.c6;
    ctype c4 = A.c0 * B.c7 + A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c1;
    ctype c5 = -A.c0 * B.c3 + A.c1 * B.c7 + A.c2 * B.c1 + A.c3 * B.c2;
    ctype c6 = A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c7 + A.c3 * B.c3;
    ctype c7 = A.c0 * B.c4 + A.c1 * B.c5 + A.c2 * B.c6 + A.c3 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product A * B of a multivector A with the trivector B from the right
// multivector * trivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& A,
                                                            PScalar3d<U> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(-A.c7, -A.c4, -A.c5, -A.c6, A.c1, A.c2, A.c3, A.c0) * ctype(B);
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
// to the multivector B
// trivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                            MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) * MVec3d<ctype>(-B.c7, -B.c4, -B.c5, -B.c6, B.c1, B.c2, B.c3, B.c0);
}

// geometric product B * M of bivector B with multivector M
// bivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& M,
                                                            BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = -M.c4 * B.x - M.c5 * B.y - M.c6 * B.z;
    ctype c1 = -M.c2 * B.z + M.c3 * B.y - M.c7 * B.x;
    ctype c2 = M.c1 * B.z - M.c3 * B.x - M.c7 * B.y;
    ctype c3 = -M.c1 * B.y + M.c2 * B.x - M.c7 * B.z;
    ctype c4 = M.c0 * B.x - M.c5 * B.z + M.c6 * B.y;
    ctype c5 = M.c0 * B.y + M.c4 * B.z - M.c6 * B.x;
    ctype c6 = M.c0 * B.z - M.c4 * B.y + M.c5 * B.x;
    ctype c7 = M.c1 * B.x + M.c2 * B.y + M.c3 * B.z;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product B * M of a bivector B with a mulitvector M
// bivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(BiVec3d<U> const& B,
                                                            MVec3d<T> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = -B.x * M.c4 - B.y * M.c5 - B.z * M.c6;
    ctype c1 = -B.x * M.c7 - B.y * M.c3 + B.z * M.c2;
    ctype c2 = B.x * M.c3 - B.y * M.c7 - B.z * M.c1;
    ctype c3 = -B.x * M.c2 + B.y * M.c1 - B.z * M.c7;
    ctype c4 = B.x * M.c0 - B.y * M.c6 + B.z * M.c5;
    ctype c5 = B.x * M.c6 + B.y * M.c0 - B.z * M.c4;
    ctype c6 = -B.x * M.c5 + B.y * M.c4 + B.z * M.c0;
    ctype c7 = B.x * M.c1 + B.y * M.c2 + B.z * M.c3;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}


// geometric product M * v of a mulitvector M with a vector v
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& M,
                                                            Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c1 * v.x + M.c2 * v.y + M.c3 * v.z;
    ctype c1 = M.c0 * v.x - M.c5 * v.z + M.c6 * v.y;
    ctype c2 = M.c0 * v.y + M.c4 * v.z - M.c6 * v.x;
    ctype c3 = M.c0 * v.z - M.c4 * v.y + M.c5 * v.x;
    ctype c4 = M.c2 * v.z - M.c3 * v.y + M.c7 * v.x;
    ctype c5 = -M.c1 * v.z + M.c3 * v.x + M.c7 * v.y;
    ctype c6 = M.c1 * v.y - M.c2 * v.x + M.c7 * v.z;
    ctype c7 = M.c4 * v.x + M.c5 * v.y + M.c6 * v.z;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product v * M of a vector v with a multivector M
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                            MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = v.x * M.c1 + v.y * M.c2 + v.z * M.c3;
    ctype c1 = v.x * M.c0 - v.y * M.c6 + v.z * M.c5;
    ctype c2 = v.x * M.c6 + v.y * M.c0 - v.z * M.c4;
    ctype c3 = -v.x * M.c5 + v.y * M.c4 + v.z * M.c0;
    ctype c4 = v.x * M.c7 + v.y * M.c3 - v.z * M.c2;
    ctype c5 = -v.x * M.c3 + v.y * M.c7 + v.z * M.c1;
    ctype c6 = v.x * M.c2 - v.y * M.c1 + v.z * M.c7;
    ctype c7 = v.x * M.c4 + v.y * M.c5 + v.z * M.c6;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// geometric product A * s of mulitvector A with scalar s
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(MVec3d<T> const& M,
                                                            Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// geometric product s * B of scalar s with multivector B
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                            MVec3d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// geometric product A * B for two multivectors from the even subalgebra (3d case)
// even grade multivector * even grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(
        Scalar3d<ctype>(A.c0 * B.c0 - A.c1 * B.c1 - A.c2 * B.c2 - A.c3 * B.c3),
        BiVec3d<ctype>(A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2,
                       A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1,
                       A.c0 * B.c3 - A.c1 * B.c2 + A.c2 * B.c1 + A.c3 * B.c0));
}


// geometric product A * B for two multivectors from the uneven subalgebra (3d case)
// uneven grade multivector * uneven grade multivector = even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(
        Scalar3d<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c3),
        BiVec3d<ctype>(A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0,
                       -A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 + A.c3 * B.c1,
                       A.c0 * B.c1 - A.c1 * B.c0 + A.c2 * B.c3 + A.c3 * B.c2));
}

// geometric product A * B of a multivector A from the even subalgebra (3d case)
// with a multivector B of the uneven subalgebra
// even grade multivector * uneven grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(
        Vec3d<ctype>(A.c0 * B.c0 - A.c1 * B.c3 - A.c2 * B.c2 + A.c3 * B.c1,
                     A.c0 * B.c1 + A.c1 * B.c2 - A.c2 * B.c3 - A.c3 * B.c0,
                     A.c0 * B.c2 - A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar3d<ctype>(A.c0 * B.c3 + A.c1 * B.c0 + A.c2 * B.c1 + A.c3 * B.c2));
}

// geometric product A * B of a multivector A from the uneven subalgebra (3d case)
// with a multivector B of the even subalgebra
// uneven grade multivector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(
        Vec3d<ctype>(A.c0 * B.c0 - A.c1 * B.c3 + A.c2 * B.c2 - A.c3 * B.c1,
                     A.c0 * B.c3 + A.c1 * B.c0 - A.c2 * B.c1 - A.c3 * B.c2,
                     -A.c0 * B.c2 + A.c1 * B.c1 + A.c2 * B.c0 - A.c3 * B.c3),
        PScalar3d<ctype>(A.c0 * B.c1 + A.c1 * B.c2 + A.c2 * B.c3 + A.c3 * B.c0));
}

// geometric product M * ps of an even multivector M multiplied from the right
// by the trivector ps (=3d pseudoscalar)
// even grade multivector * trivector => uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<U> const& M,
                                                              PScalar3d<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-M.c1, -M.c2, -M.c3), PScalar3d<ctype>(M.c0)) *
           ctype(ps);
}

// geometric product ps * M of a trivector ps (=3d pseudoscalar) multiplied from the left
// to the even grade multivector M
// trivector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(PScalar3d<T> ps,
                                                              MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) *
           MVec3d_U<ctype>(Vec3d<ctype>(-M.c1, -M.c2, -M.c3), PScalar3d<ctype>(M.c0));
}

// geometric product M * ps of an uneven grade multivector M multiplied from the right
// by the trivector ps (=3d pseudoscalar)
// uneven grade multivector * trivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<U> const& M,
                                                              PScalar3d<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(-M.c3), BiVec3d<ctype>(M.c0, M.c1, M.c2)) *
           ctype(ps);
}

// geometric product ps * M of a trivector ps (=3d pseudoscalar) multiplied from the left
// to the uneven grade multivector M
// trivector * uneven grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(PScalar3d<T> ps,
                                                              MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) *
           MVec3d_E<ctype>(Scalar3d<ctype>(-M.c3), BiVec3d<ctype>(M.c0, M.c1, M.c2));
}

// geometric product M * B of a multivector M from the even subalgebra (3d case)
// with a bivector B from the right
// even grade multivector * bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& M,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(-M.c1 * B.x - M.c2 * B.y - M.c3 * B.z),
                           BiVec3d<ctype>(M.c0 * B.x - M.c2 * B.z + M.c3 * B.y,
                                          M.c0 * B.y + M.c1 * B.z - M.c3 * B.x,
                                          M.c0 * B.z - M.c1 * B.y + M.c2 * B.x));
}

// geometric product B * M of a bivector A with an even grade multivector M
// from the right
// bivector * even grade multivector => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                              MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(-B.x * M.c1 - B.y * M.c2 - B.z * M.c3),
                           BiVec3d<ctype>(B.x * M.c0 - B.y * M.c3 + B.z * M.c2,
                                          B.x * M.c3 + B.y * M.c0 - B.z * M.c1,
                                          -B.x * M.c2 + B.y * M.c1 + B.z * M.c0));
}

// geometric product M * B of an uneven grade multivector M with a bivector B
// uneven grade multivector * bivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& M,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-M.c1 * B.z + M.c2 * B.y - M.c3 * B.x,
                                        M.c0 * B.z - M.c2 * B.x - M.c3 * B.y,
                                        -M.c0 * B.y + M.c1 * B.x - M.c3 * B.z),
                           PScalar3d<ctype>(M.c0 * B.x + M.c1 * B.y + M.c2 * B.z));
}

// geometric product B * M of a bivector B with an uneven grade multivector M
// bivector * uneven grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                              MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-B.x * M.c3 - B.y * M.c2 + B.z * M.c1,
                                        B.x * M.c2 - B.y * M.c3 - B.z * M.c0,
                                        -B.x * M.c1 + B.y * M.c0 - B.z * M.c3),
                           PScalar3d<ctype>(B.x * M.c0 + B.y * M.c1 + B.z * M.c2));
}

// geometric product M * v of an even grade multivector M with a vector v
// from the right
// even grade multivector * vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& M,
                                                              Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(M.c0 * v.x - M.c2 * v.z + M.c3 * v.y,
                                        M.c0 * v.y + M.c1 * v.z - M.c3 * v.x,
                                        M.c0 * v.z - M.c1 * v.y + M.c2 * v.x),
                           PScalar3d<ctype>(M.c1 * v.x + M.c2 * v.y + M.c3 * v.z));
}

// geometric product v * M of a vector v and an even grade multivector M from the left
// vector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                              MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(v.x * M.c0 - v.y * M.c3 + v.z * M.c2,
                                        v.x * M.c3 + v.y * M.c0 - v.z * M.c1,
                                        -v.x * M.c2 + v.y * M.c1 + v.z * M.c0),
                           PScalar3d<ctype>(v.x * M.c1 + v.y * M.c2 + v.z * M.c3));
}

// geometric product M * v of an uneven grade multivector M with a vector v
// from the right
// uneven grade multivector * vector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& M,
                                                              Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(M.c0 * v.x + M.c1 * v.y + M.c2 * v.z),
                           BiVec3d<ctype>(M.c1 * v.z - M.c2 * v.y + M.c3 * v.x,
                                          -M.c0 * v.z + M.c2 * v.x + M.c3 * v.y,
                                          M.c0 * v.y - M.c1 * v.x + M.c3 * v.z));
}

// geometric product v * M of a vector v and an uneven grade multivector B from the left
// vector * uneven grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                              MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(v.x * M.c0 + v.y * M.c1 + v.z * M.c2),
                           BiVec3d<ctype>(v.x * M.c3 + v.y * M.c2 - v.z * M.c1,
                                          -v.x * M.c2 + v.y * M.c3 + v.z * M.c0,
                                          v.x * M.c1 - v.y * M.c0 + v.z * M.c3));
}

// geometric product M * s of even mulitvector M with scalar s
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& M,
                                                              Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// geometric product s * M of scalar s with even multivector M
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                              MVec3d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// geometric product M * s of uneven multivector M with scalar s
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& M,
                                                              Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// geometric product s * M of scalar s with eneven multivector M
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                              MVec3d_U<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// geometric product ps1 * ps2 of two trivectors
// trivector * trivector => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator*(PScalar3d<T> ps1,
                                                              PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(-ctype(ps1) * ctype(ps2)); // trivectors in 3d square to -1
}

// geometric product B * ps of bivector B multiplied by a trivector ps from the right
// bivector * trivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                           PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(-B.x, -B.y, -B.z) * ctype(ps);
}

// geometric product ps * B of a trivector ps (=3d pseudoscalar) multiplied from the left
// to the bivector B
// trivector * bivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> ps,
                                                           BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec3d<ctype>(-B.x, -B.y, -B.z);
}

// geometric product v * ps of a vector a multiplied with a trivector B from the right
// vector * trivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                             PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(v.x, v.y, v.z) * ctype(ps);
}

// geometric product ps * v of a trivector A (=3d pseudoscalar) multiplied from the left
// to the vector v
// trivector * vector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(PScalar3d<T> ps,
                                                             Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * BiVec3d<ctype>(v.x, v.y, v.z);
}

// geometric product ps * s of a trivector ps and a scalar s
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator*(PScalar3d<T> ps,
                                                               Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) *
                            ctype(s)); // scalar multiplication with a trivector
}

// geometric product s * ps of a scalar s and a trivector ps
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                               PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(s) *
                            ctype(ps)); // scalar multiplication with a trivector
}

// geometric product A * B between two bivectors
// A * B = gr0(A * B) + gr2(A * B)
//
// the full geometric bivector product only exists in >= 4d spaces:
// A * B = dot(A,B) + cmt(A,B) + wdg(A,B) = gr0(A * B) + gr2(A * B) + gr4(A * B)
// In 3D we don't have gr4(A * B) and thus only the terms up to grade 3 remain.
// The bivector product AxB = cmt(A,B) = 0.5*(A*B-B*A)is called the commutator product.
//                            cmt(B,A) = -cmt(A,B)
//
// A * B = -dot(A,B) + cmt(A,B) + wdg(A,B)  (in 4d and higher dimensional spaces)
// A * B = -dot(A,B) + cmt(A,B)             (in 3d for the case with two bivectors)
//
// => bivector * bivector = scalar + bivector = even grade multivector (in 3d)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B1,
                                                              BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(-dot(B1, B2), cmt(B1, B2));
}

// geometric product B * v for a bivector B and a vector v
// B * v = (v << B) + wdg(B,v) = gr1(B * v) + gr3(B * v)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => bivector * vector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                              Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>((v << B), wdg(B, v));
}

// geometric product v * B for a vector v and a bivector B
// v * B = (B >> v) + wdg(v,B) = gr1(v * B) + gr3(v * B)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => vector * bivector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>((B >> v), wdg(v, B));
}


// (geometric) product B * s of bivector B multiplied by a scalar s from the right
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                             Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return B * ctype(s);
}

// (geometric) product s * B of a scalar s multiplied from the left to the bivector B
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                             BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

// geometric product a * b of two vectors
// a * b = (b >> a) + wdg(a,b) = gr0(a * b) + gr2(a * b)
//
// HINT: if a full 3d multivector is required as result, it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// vector * vector => even grade multivector (= scalar + bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                              Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>((b >> a), wdg(a, b));
}

// (geometric) product v * s of a vector v multiplied with a scalar s from the right
// vector * scalar => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                           Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// (geometric) product s * v of a scalar s multiplied from the left to the vector v
// scalar * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                           Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// geometric product s1 * s2 of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s1,
                                                              Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses of scalars, blades and multivectors
// w.r.t. the geometric product:
// for k-blades: A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ... (from reversion)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> inv(Scalar3d<T> s)
{
    T sq_n = nrm_sq(s);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("scalar norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n;

    return Scalar3d<T>(rev(s) * inv);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> inv(Vec3d<T> const& v)
{
    T sq_n = nrm_sq(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("vector norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n;
    return Vec3d<T>(rev(v) * inv);
}

// return the multiplicative inverse of the bivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> inv(BiVec3d<T> const& B)
{
    T sq_n = nrm_sq(B);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("bivector norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n;
    return BiVec3d<T>(rev(B) * inv); // minus sign due to reversion
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> inv(PScalar3d<T> ps)
{
    T sq_n = nrm_sq(ps);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("trivector norm too small for inversion " +
                                 std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector

    return PScalar3d<T>(rev(ps) * inv); // minus sign due to reversion
}

// return the multiplicative inverse of the quaternion (inv(z) = 1/nrm_sq(z)*rev(z))
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> inv(MVec3d_E<T> const& E)
{
    T sq_n = nrm_sq(E);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error(
            "norm of even grade multivector too small for inversion " +
            std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n;
    return MVec3d_E<T>(rev(E) * inv);
}

// return the multiplicative inverse of the uneven grade multivector
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> inv(MVec3d_U<T> const& U)
{
    T sq_n = nrm_sq(U);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (sq_n < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error(
            "norm of uneven grade multivector too small for inversion " +
            std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n;
    return MVec3d_U<T>(rev(U) * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord algebras",
// Hitzer, Sangwine, 2016
// left and a right inverse are the same (see paper of Hitzer, Sangwine)
template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> inv(MVec3d<T> const& M)
{
    T m_conjm = gr0(M * conj(M) * gr_inv(M) * rev(M));

#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(m_conjm) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("multivector norm too small for inversion " +
                                 std::to_string(m_conjm) + "\n");
    }
#endif
    T inv = T(1.0) / m_conjm;
    return MVec3d<T>(conj(M) * gr_inv(M) * rev(M) * inv);
}


////////////////////////////////////////////////////////////////////////////////
// angle operations
////////////////////////////////////////////////////////////////////////////////

// return the angle between of two vectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    ctype nrm_prod = nrm(v1) * nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle " +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}


// unsuccessful try to extend angle range to -pi ... pi,
// because orientation is not defined uniquely:
//
// // return the angle between of two vectors
// // range of angle: -pi <= angle <= pi
// template <typename T, typename U>
//     requires(std::floating_point<T> && std::floating_point<U>)
// inline constexpr std::common_type_t<T, U> angle(Vec3d<T> const& v1, Vec3d<U>
// const& v2)
// {
//     using ctype = std::common_type_t<T, U>;
//     using std::numbers::pi;

// ctype nrm_prod = nrm(v1) * nrm(v2);
// #if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
// if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
//     throw std::runtime_error("vector norm product too small for calculation "
//                              "of angle " +
//                              std::to_string(nrm_prod) + "\n");
// }
// #endif

// auto cos_angle = std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0));
// auto sin_angle = std::clamp(ctype(nrm(wdg(v1, v2))) / nrm_prod, ctype(-1.0),
// ctype(1.0));
//     // wdg() does contain magnitude, but no unique value of orientation
//     // so we chose one arbitrarily => but would deliver only pos. angles!

//     fmt::println("   c = {: .4f}, s = {: .4f}, wdg = {: .4f}, nrm_wdg = {: .4f}",
//                  cos_angle, sin_angle, wdg(v1, v2), nrm(wdg(v1, v2)));

//     if (cos_angle >= 0.0) {
//         // quadrant I or IV
//         return std::asin(sin_angle);
//     }
//     else if (cos_angle < 0.0 && sin_angle >= 0.0) {
//         // quadrant II
//         return pi - std::asin(sin_angle);
//     }
//     else {
//         // cos_angle < 0.0 && sin_angle < 0.0)
//         // quadrant III
//         return -pi - std::asin(sin_angle);
//     }
// }

// return the angle between two bivectors
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(BiVec3d<T> const& v1,
                                                BiVec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v1) * nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle " +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}


// return the angle between of a vector and a bivector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec3d<T> const& v, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(v) * nrm(B);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error("norm product too small for calculation of angle " +
                                 std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    // return std::acos(
    //     std::clamp(ctype(nrm(dot(v, B))) / nrm_prod, ctype(-1.0), ctype(1.0)));
    return std::acos(
        std::clamp(ctype(nrm((v << B))) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between of a bivector and a vector
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(BiVec3d<T> const& B, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype nrm_prod = nrm(B) * nrm(v);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error("norm product too small for calculation of angle " +
                                 std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    // return std::acos(
    //     std::clamp(ctype(nrm(dot(B, v))) / nrm_prod, ctype(-1.0), ctype(1.0)));
    return std::acos(
        std::clamp(ctype(nrm((B >> v))) / nrm_prod, ctype(-1.0), ctype(1.0)));
}


////////////////////////////////////////////////////////////////////////////////
// 3d rotation operations
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////
// exponential function with bivector as argument for setup of quaternions
// as geometric multivector with a scalar and a bivector part
// MVec3d_E<T> M = c0 + (c1 e2^e3 + c2 e3^e1 + c3 e1^e2)
//
// quaternion: q = a + b I with I being the bivector in brackets above
//             representing a plane in the algebra G^3
//
//             a rotation in 3D is represented by the plane and the
//             size of the rotation, the later is given by the angle
//             theta, which is the magnitude of the bivector
//
// Inputs:
//         - an arbitray bivector representing the oriented plane of rotation
//           (does not need to be normalized)
//         - a rotation angle
// Output:
//         - a rotor representing the rotation
//
// HINT:     For a rotation around an axis n (n = normalize(Vec3d<T>))
//           use the bivector B = n*I_3d  => B = Vec3d<T> * PScalar3d<T>
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> exp(BiVec3d<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(std::cos(theta)),
                           normalize(I) * std::sin(theta));
}

//////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//       1.) an arbitray bivector representing the oriented plane of rotation
//           (does not need to be normalized, defines what is a posive rotation angle)
//       2.) a rotation angle in that plane
// Output:
//           a rotor representing the requested rotation,
//           for applying the sandwich product as in rotate(v,rotor)
//
//////////////////////////////////////////////////////////////////////////////////////////
//
// for a rotation about an axis n (n = normalized vector) choose the ansatz n*B = I_3d
// and multiply both sides with n from the left (remember n*n = |n|^2 = 1)
//
// => choose: B = n*I_3d
//
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> rotor(BiVec3d<T> const& I, U theta)
{
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec3d_E<ctype>(Scalar3d<ctype>(std::cos(half_angle)),
                           normalize(I) * std::sin(half_angle));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> rotate(Vec3d<T> const& v,
                                                        MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3d_E<ctype> rr = rev(rotor);
    // MVec3d_U<ctype> tmp = rotor * v;
    // MVec3d_U<ctype> res = tmp * rr;

    // trivector part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the pseudoscalar part
    // which will be zero anyway
    // return Vec3d<ctype>(gr1<ctype>(rotor * v * rev(rotor)));

    return Vec3d<ctype>(gr1(rotor * v * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> rotate_opt1(Vec3d<T> const& v,
                                                             MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    // optimized version:
    MVec3d_E<ctype> rr = rev(rotor);
    MVec3d_U<ctype> tmp = rotor * v;
    // formular from operator*(MVec3d_U<T>, MVec3d_E<U>) - only vector part
    return Vec3d<ctype>(tmp.c0 * rr.c0 - tmp.c1 * rr.c3 + tmp.c2 * rr.c2 - tmp.c3 * rr.c1,
                        tmp.c0 * rr.c3 + tmp.c1 * rr.c0 - tmp.c2 * rr.c1 - tmp.c3 * rr.c2,
                        -tmp.c0 * rr.c2 + tmp.c1 * rr.c1 + tmp.c2 * rr.c0 -
                            tmp.c3 * rr.c3);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> rotate_opt2(Vec3d<T> const& v,
                                                             MVec3d_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;
    ctype k1 = R.c0 * v.x - R.c2 * v.z + R.c3 * v.y;
    ctype k2 = R.c0 * v.y + R.c1 * v.z - R.c3 * v.x;
    ctype k3 = R.c0 * v.z - R.c1 * v.y + R.c2 * v.x;
    ctype k4 = R.c1 * v.x + R.c2 * v.y + R.c3 * v.z;
    return Vec3d<ctype>(k1 * R.c0 + k2 * R.c3 - k3 * R.c2 + k4 * R.c1,
                        -k1 * R.c3 + k2 * R.c0 + k3 * R.c1 + k4 * R.c2,
                        k1 * R.c2 - k2 * R.c1 + k3 * R.c0 + k4 * R.c3);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> rotate(BiVec3d<T> const& B,
                                                          MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // scalar part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the scalar part
    // which will be zero anyway

    return BiVec3d<ctype>(gr2(rotor * B * rev(rotor)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> rotate_opt1(BiVec3d<T> const& B,
                                                               MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // optimized version:
    MVec3d_E<ctype> rr = rev(rotor);
    MVec3d_E<ctype> tmp = rotor * B;
    // formular from operator*(MVec3d_E<T>, MVec3d_E<U>) - only bivector part
    return BiVec3d<ctype>(
        tmp.c0 * rr.c1 + tmp.c1 * rr.c0 - tmp.c2 * rr.c3 + tmp.c3 * rr.c2,
        tmp.c0 * rr.c2 + tmp.c1 * rr.c3 + tmp.c2 * rr.c0 - tmp.c3 * rr.c1,
        tmp.c0 * rr.c3 - tmp.c1 * rr.c2 + tmp.c2 * rr.c1 + tmp.c3 * rr.c0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> rotate_opt2(BiVec3d<T> const& B,
                                                               MVec3d_E<U> const& R)
{
    using ctype = std::common_type_t<T, U>;

    ctype k1 = R.c1 * B.x + R.c2 * B.y + R.c3 * B.z;
    ctype k2 = R.c0 * B.x - R.c2 * B.z + R.c3 * B.y;
    ctype k3 = R.c0 * B.y + R.c1 * B.z - R.c3 * B.x;
    ctype k4 = R.c0 * B.z - R.c1 * B.y + R.c2 * B.x;
    return BiVec3d<ctype>(k1 * R.c1 + k2 * R.c0 + k3 * R.c3 - k4 * R.c2,
                          k1 * R.c2 - k2 * R.c3 + k3 * R.c0 + k4 * R.c1,
                          k1 * R.c3 + k2 * R.c2 - k3 * R.c1 + k4 * R.c0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> rotate(MVec3d<T> const& M,
                                                         MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(rotor * M * rev(rotor));
}


////////////////////////////////////////////////////////////////////////////////
// 3d complement operations
// (the concept of complement is defined w.r.t. the outer product)
// it depends on the pseudoscalar of the space
// here implemented for I_3d = e1^e2^e3 = e123
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B of the blade u as subspace of R^3 then
// compl(M) represents the subspace orthorgonal to B
// the complement exchanges basis vectors (exclusively, no impact on magnitude),
// which are in the k-blade u with the basis vectors which are NOT contained in the
// k-blade u and are needed to fill the space completely to the corresponding pseudoscalar
//
// left complement:  lcmpl(u) ^ u  = I_3d = e1^e2^e3  =>  lcmpl(u) = I_3d * rev(u)
// right complement: u ^ rcmpl(u)  = I_3d = e1^e2^e3  =>  rcmpl(u) = rev(v) * I_3d
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined lcmpl(u) == rcmpl(u) == cmpl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// complements have different signs

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> cmpl(Scalar3d<T> s)
{
    // u ^ cmpl(u) = e1^e2^e3
    // u = s 1:
    //     s ^ cmpl(u) = e1^e2^e3 => cmpl(u) = rev(s) * I_3d = s e1^e2^e3
    return PScalar3d<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> cmpl(Vec3d<T> const& v)
{
    // u ^ compl(u) = e1^e2^e3
    // u = v.x e1 + v.y e2 + v.z e3:
    //     u ^ cmpl(u) = e1^e2^e3 => cmpl(u) = rev(v) * I_3d
    //     e1 => cmpl(u) = v.x e23
    //     e2 => cmpl(u) = v.y e31
    //     e3 => cmpl(u) = v.z e12
    return BiVec3d<T>(v.x, v.y, v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> cmpl(BiVec3d<T> const& B)
{
    // u ^ compl(u) = e1^e2^e3
    // u = B.x e23 + B.y e31 + B.z e12:
    //     u ^ cmpl(u) = e1^e2^e3 => cmpl(u) = rev(B) * I_3d
    //     e23 => cmpl(u) = B.x e1
    //     e31 => cmpl(u) = B.y e2
    //     e12 => cmpl(u) = B.z e3
    return Vec3d<T>(B.x, B.y, B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> cmpl(PScalar3d<T> ps)
{
    // u ^ compl(u) = e1^e2^e3
    // u = ps e1^e2^e3:
    //     u ^ cmpl(u) = e1^e2^e3 => cmpl(u) = rev(ps) * I_3d = ps 1
    return Scalar3d<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> cmpl(MVec3d_E<T> const& M)
{
    // use the component complements directly
    return MVec3d_U<T>(cmpl(gr2(M)), cmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> cmpl(MVec3d_U<T> const& M)
{
    // use the component complements directly
    return MVec3d_E<T>(cmpl(gr3(M)), cmpl(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> cmpl(MVec3d<T> const& M)
{
    // use the component complements directly
    return MVec3d<T>(cmpl(gr3(M)), cmpl(gr2(M)), cmpl(gr1(M)), cmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// duality (as defined in Lengyel, "PGA illuminated")
// is defined w.r.t. the outer product
////////////////////////////////////////////////////////////////////////////////
//
// if M represents the subspace B as subspace of R^3 then
// dual(M) represents the subspace orthorgonal to B
//
// dual(A) = cmpl(A) in spaces of uneven dimension
//         = rcmpl(A) in spaces of even dimension (right dual)
//
// this dual satisfies (right) dual(A) = rev(A) * I_n
//
// for the left cmpl/dual the (left) dual: ldual(L) = I_n * rev(L)
//
// -> derived from the defining equation of the left and right complements


template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> dual(Scalar3d<T> s)
{
    return cmpl(s);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> dual(Vec3d<T> const& v)
{
    return cmpl(v);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> dual(BiVec3d<T> const& B)
{
    return cmpl(B);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> dual(PScalar3d<T> ps)
{
    return cmpl(ps);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> dual(MVec3d_E<T> const& M)
{
    return MVec3d_U<T>(cmpl(gr2(M)), cmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> dual(MVec3d_U<T> const& M)
{
    return MVec3d_E<T>(cmpl(gr3(M)), cmpl(gr1(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> dual(MVec3d<T> const& M)
{
    return MVec3d<T>(cmpl(gr3(M)), cmpl(gr2(M)), cmpl(gr1(M)), cmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> and BiVec3d<T> projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
// v_parallel = dot(v1, v2)) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v1,
                                                              Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2);
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v1,
                                                             Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v1 - project_onto(v1, v2));

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec3d<ctype>(gr1(wdg(v1, v2) * inv(v2)));
}

// projection of a vector v onto a bivector B
// v_parallel = gr1((B >> v) * inv(B))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v,
                                                              BiVec3d<U> const& B)
{
    return gr1((B >> v) * inv(B));
}

// rejection of vector v1 from a bivector v2
// u_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v,
                                                             BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v - project_onto(v, B));

    // works, but is more effort compared to solution via projection and vector difference
    // return gr1(wdg(v, B) * inv(B));
}

////////////////////////////////////////////////////////////////////////////////
// reflections
////////////////////////////////////////////////////////////////////////////////
// Macdonald p. 129:
//
// Reflect a j-blade u in the k-dimensional subspace B is the blade
//
//   u_reflected = (-1)^[j*(k+1)]*B*u*inv(B)
//
////////////////////////////////////////////////////////////////////////////////

// reflect a vector v on a hyperplane B orthogonal to vector nB
//
// hyperplane: a n-1 dimensional subspace in a space of dimension n (a plane in 3d space)
// orthogonal to vector nB: the hyperplane B is dual to nB (i.e. a two dimensional
// subspace)
//
// HINT: choose nB = dual(B) (for normalized nB)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reflect_on(Vec3d<T> const& v,
                                                            Vec3d<U> const& nB)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(gr1(-nB * v * inv(nB)));
}

// reflect a vector v in an arbitrary bivector B, i.e. a plane
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reflect_on(Vec3d<T> const& v,
                                                            BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(gr1(-B * v * inv(B)));
}

// reflect a bivector UB in an arbitrary bivector B (both modelling planes)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> reflect_on(BiVec3d<T> const& UB,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(gr2(B * UB * inv(B)));
}

// reflect a vector v on another vector b
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reflect_on_vec(Vec3d<T> const& v,
                                                                Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(gr1(b * v * inv(b)));
}

////////////////////////////////////////////////////////////////////////////////
// Gram-Schmidt-Orthogonalization part 1: 2d plane embedded in 3d space
////////////////////////////////////////////////////////////////////////////////
//
// input:  two linear independent vectors u and v in 3d defining a plane
// output: two orthogonal vectors with the first one being u and the second one a vector
// perpendicular to u in the orientation of v, both forming an orthogonal system
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
std::vector<Vec3d<std::common_type_t<T, U>>> gs_orthogonal(Vec3d<T> const& u,
                                                           Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec3d<ctype>> basis;
    basis.push_back(u);
    basis.emplace_back(reject_from(v, u));
    return basis;
}

// input:  two linear independent vectors u and v in 3d defining a plane
// output: two orthonormal vectors with the first one being normalize(u) and
// the second one a normalized vector perpendicular to u in the orientation of v,
// both forming an orthogonal system
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
std::vector<Vec3d<std::common_type_t<T, U>>> gs_orthonormal(Vec3d<T> const& u,
                                                            Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec3d<ctype>> basis;
    Vec3d<ctype> u_unitized{normalize(u)};
    basis.push_back(u_unitized);
    basis.emplace_back(normalize(reject_from(v, u_unitized)));
    return basis;
}

////////////////////////////////////////////////////////////////////////////////
// Gram-Schmidt-Orthogonalization part 2: 3d space
////////////////////////////////////////////////////////////////////////////////
//
// input:  three linear independent vectors u, v and w in 3d
// output: three orthogonal vectors with the first one being u and the second and third
// being  perpendicular to u and the plane spanned by u and v respectively.
// All three from an orthogonal system
template <typename U, typename V, typename W>
    requires(std::floating_point<U> && std::floating_point<V> && std::floating_point<W>)
inline constexpr std::vector<Vec3d<std::common_type_t<U, V, W>>>
gs_orthogonal(Vec3d<U> const& u, Vec3d<V> const& v, Vec3d<W> const& w)
{
    using ctype = std::common_type_t<U, V, W>;
    std::vector<Vec3d<ctype>> basis;
    basis.reserve(3);
    basis.push_back(u);
    basis.emplace_back(reject_from(v, u));
    basis.emplace_back(reject_from(w, wdg(u, v)));
    return basis;
}

// input:  three linear independent vectors u, v and w in 3d
// output: three orthonormal vectors with the first one being normalized u and the second
// and third being normalized and perpendicular to u and the plane spanned by u and v
// respectively. All three from an orthogonal system
template <typename U, typename V, typename W>
    requires(std::floating_point<U> && std::floating_point<V> && std::floating_point<W>)
inline constexpr std::vector<Vec3d<std::common_type_t<U, V, W>>>
gs_orthonormal(Vec3d<U> const& u, Vec3d<V> const& v, Vec3d<W> const& w)
{
    using ctype = std::common_type_t<U, V, W>;
    std::vector<Vec3d<ctype>> basis;
    basis.reserve(3);
    Vec3d<ctype> u_unitized{normalize(u)};
    basis.push_back(u_unitized);
    basis.emplace_back(normalize(reject_from(v, u_unitized)));
    basis.emplace_back(normalize(reject_from(w, wdg(u, v))));
    return basis;
}

} // namespace hd::ga::ega