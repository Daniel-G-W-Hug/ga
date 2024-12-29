#pragma once

// author: Daniel Hug, 2024

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
inline constexpr MVec3d_E<T> gr_inv(MVec3d_E<T> const& E)
{
    // grade 0 and 2: no sign change
    return E;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> gr_inv(MVec3d_U<T> const& U)
{
    // grade 1 and 3: sign reversal
    return -U;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> gr_inv(MVec3d<T> const& M)
{
    // grade 0 and 2: no sign change
    // grade 1 and 3: sign reversal
    return MVec3d<T>(M.c0, -M.c1, -M.c2, -M.c3, M.c4, M.c5, M.c6, -M.c7);
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
inline constexpr MVec3d_E<T> rev(MVec3d_E<T> const& E)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec3d_E<T>(E.c0, -E.c1, -E.c2, -E.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> rev(MVec3d_U<T> const& U)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec3d_U<T>(U.c0, U.c1, U.c2, -U.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> rev(MVec3d<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    return MVec3d<T>(M.c0, M.c1, M.c2, M.c3, -M.c4, -M.c5, -M.c6, -M.c7);
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
inline constexpr MVec3d_E<T> conj(MVec3d_E<T> const& E)
{
    // grade 0: no sign change
    // grade 2: sign change
    return MVec3d_E<T>(E.c0, -E.c1, -E.c2, -E.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> conj(MVec3d_U<T> const& U)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec3d_U<T>(-U.c0, -U.c1, -U.c2, U.c3);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> conj(MVec3d<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    return MVec3d<T>(M.c0, -M.c1, -M.c2, -M.c3, -M.c4, -M.c5, -M.c6, M.c7);
}


////////////////////////////////////////////////////////////////////////////////
// scalar product (= dot product defined for equal grades exclusively)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(Scalar3d<T> s1, Scalar3d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s1) * ctype(s2));
}

// return the dot product of two vectors in G<3,0,0>
// coordinate free definition: dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(Vec3d<T> const& v1,
                                                        Vec3d<U> const& v2)
{
    // definition: dot(v1, v2) = (v1)^T g_12 v2 using the metric g_12
    // definition: dot(A, B) = gr0(A*B) using the geometric product
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1, e3^2 = +1
    // as diagonal elements of g_12
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

// return dot product of two bivectors A and B (= a scalar)
// dot(A,B) = gr0(A * B)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(BiVec3d<T> const& B1,
                                                        BiVec3d<U> const& B2)
{
    // definition: dot(A, B) = gr0(A*B)
    // -> only the symmetric (i.e. scalar) part remains
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1, e3^2 = +1
    // and dot(e23, e23) = -1, dot(e31,e31) = -1, dot(e12,e12) = -1
    // and all other dot(exy,ezw)=0
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(-B1.x * B2.x - B1.y * B2.y - B1.z * B2.z);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(PScalar3d<T> ps1,
                                                        PScalar3d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(-ctype(ps1) * ctype(ps2));
}

// scalar product dot(a,b) (nrm_sq(a,b) = dot(a, rev(b)))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> dot(MVec3d<T> const& A,
                                                        MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 -
                           A.c4 * B.c4 - A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7);
}


////////////////////////////////////////////////////////////////////////////////
// wedge product (= outer product)
////////////////////////////////////////////////////////////////////////////////

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * Vec3d<ctype>(v.x, v.y, v.z);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> wdg(Vec3d<T> const& v, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v.x, v.y, v.z) * ctype(s);
}

// wedge product with a pscalar (returns a scaled pscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> wdg(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(s) * ctype(ps));
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> wdg(PScalar3d<T> ps, Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) * ctype(s));
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

// wedge product extended to a fully populated multivectors
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

// regressive wedge product between a scalar s and a pseudoscalar ps (=trivector)
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> rwdg(Scalar3d<T> s, PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(ctype(s) * ctype(ps));
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

// regressive wedge product extended to a fully populated multivectors
// => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> rwdg(MVec3d<T> const& A,
                                                       MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;

    ctype c0 = A.c0 * B.c7 + A.c7 * B.c0 + A.c1 * B.c4 + A.c2 * B.c5 + A.c3 * B.c6 +
               A.c4 * B.c1 + A.c5 * B.c2 + A.c6 * B.c3;
    ctype c1 = A.c1 * B.c7 + A.c7 * B.c1 - A.c6 * B.c5 + A.c5 * B.c6;
    ctype c2 = A.c2 * B.c7 + A.c7 * B.c2 - A.c4 * B.c6 + A.c6 * B.c4;
    ctype c3 = A.c3 * B.c7 + A.c7 * B.c3 - A.c5 * B.c4 + A.c4 * B.c5;
    ctype c4 = A.c4 * B.c7 + A.c7 * B.c4;
    ctype c5 = A.c5 * B.c7 + A.c7 * B.c5;
    ctype c6 = A.c6 * B.c7 + A.c7 * B.c6;
    ctype c7 = A.c7 * B.c7;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}


////////////////////////////////////////////////////////////////////////////////
// left contractions
////////////////////////////////////////////////////////////////////////////////

// left contraction (A << B) - "A contracted onto B" == "A removed from B"
//                             "A taken out of B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator<<(MVec3d<T> const& A,
                                                             MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c4 -
               A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype c1 = A.c0 * B.c1 - A.c4 * B.c7 - A.c2 * B.c6 + A.c3 * B.c5;
    ctype c2 = A.c0 * B.c2 - A.c5 * B.c7 - A.c3 * B.c4 + A.c1 * B.c6;
    ctype c3 = A.c0 * B.c3 - A.c6 * B.c7 - A.c1 * B.c5 + A.c2 * B.c4;
    ctype c4 = A.c0 * B.c4 + A.c1 * B.c7;
    ctype c5 = A.c0 * B.c5 + A.c2 * B.c7;
    ctype c6 = A.c0 * B.c6 + A.c3 * B.c7;
    ctype c7 = A.c0 * B.c7;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// left contraction (s << v) - scalar s contracted onto vector v
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                            Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * Vec3d<ctype>(v.x, v.y, v.z);
}

// left contraction (s << B) - scalar s contracted onto bivector B
// identical with scalar multiplication of the bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator<<(Scalar3d<T> s,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * BiVec3d<ctype>(B.x, B.y, B.z);
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

// left contraction (B1 << B2) - bivector B1 contracted onto B2
// (=identical to scalar product dot(B1,B2) of two bivectors; identical to (B1 >> B2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B1,
                                                               BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(-B1.x * B2.x - B1.y * B2.y - B1.z * B2.z);
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
    return Vec3d<ctype>(-v.y * B.z + v.z * B.y, -v.z * B.x + v.x * B.z,
                        -v.x * B.y + v.y * B.x);
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

// left contraction of a bivector B onto a pseudoscalar ps (=trivector)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator<<(BiVec3d<T> const& B,
                                                            PScalar3d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(-B.x, -B.y, -B.z) * ctype(ps);
}


////////////////////////////////////////////////////////////////////////////////
// right contractions
////////////////////////////////////////////////////////////////////////////////

// right contraction (A >> B) - "A contracted by B" == "B removed from A"
//                             "B taken out of A"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> operator>>(MVec3d<T> const& A,
                                                             MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c4 -
               A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype c1 = A.c1 * B.c0 - A.c7 * B.c4 - A.c5 * B.c3 + A.c6 * B.c2;
    ctype c2 = A.c2 * B.c0 - A.c7 * B.c5 - A.c6 * B.c1 + A.c4 * B.c3;
    ctype c3 = A.c3 * B.c0 - A.c7 * B.c6 - A.c4 * B.c2 + A.c5 * B.c1;
    ctype c4 = A.c4 * B.c0 + A.c7 * B.c1;
    ctype c5 = A.c5 * B.c0 + A.c7 * B.c2;
    ctype c6 = A.c6 * B.c0 + A.c7 * B.c3;
    ctype c7 = A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// right contraction (v >> s) - vector v contracted by scalar s
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator>>(Vec3d<U> const& v,
                                                            Scalar3d<T> s)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v.x, v.y, v.z) * ctype(s);
}

// right contraction (B >> s) - bivector B contracted by scalar s
// identical with scalar multiplication of the bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator>>(BiVec3d<U> const& B,
                                                              Scalar3d<T> s)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(B.x, B.y, B.z) * ctype(s);
}

// right contraction (ps >> s) - pseudoscalar ps contracted by scalar s
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr PScalar3d<std::common_type_t<T, U>> operator>>(PScalar3d<U> ps,
                                                                Scalar3d<T> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar3d<ctype>(ctype(ps) * ctype(s));
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

// right contraction (B1 >> B2) - bivector B1 contracted by bivector B2
// (=identical to scalar product dot(B1,B2) of two bivectors; identical to (B1 >> B2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar3d<std::common_type_t<T, U>> operator>>(BiVec3d<T> const& B1,
                                                               BiVec3d<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar3d<ctype>(-B1.x * B2.x - B1.y * B2.y - B1.z * B2.z);
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
    return Vec3d<ctype>(-B.y * v.z + B.z * v.y, -B.z * v.x + B.x * v.z,
                        -B.x * v.y + B.y * v.x);
}

// right contraction of pseudoscalar ps (=trivector) by vector v (" ps contracted by v")
// => returns a bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator>>(PScalar3d<U> ps,
                                                              Vec3d<T> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * BiVec3d<ctype>(v.x, v.y, v.z);
}

// right contraction of pseudoscalar ps (=trivector) by bivector B "ps contracted by B"
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator>>(PScalar3d<U> ps,
                                                            BiVec3d<T> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec3d<ctype>(-B.x, -B.y, -B.z);
}

////////////////////////////////////////////////////////////////////////////////
// alternative multivector products (in use instead of contractions)
////////////////////////////////////////////////////////////////////////////////

// inner product (as defined by Hestens; like fat_dot below, but without scalar parts)
//
// inner(A,B) := sum_(r!=0,s!=0)^n gr( gr(A)_r * gr(B)_s )_|s-r|
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<T> inner(MVec3d<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c4 - A.c5 * B.c5 -
               A.c6 * B.c6 - A.c7 * B.c7;
    ctype c1 = -A.c2 * B.c6 + A.c3 * B.c5 - A.c5 * B.c3 + A.c6 * B.c2 - A.c4 * B.c7 -
               A.c7 * B.c4;
    ctype c2 = -A.c3 * B.c4 + A.c1 * B.c6 - A.c6 * B.c1 + A.c4 * B.c3 - A.c5 * B.c7 -
               A.c7 * B.c5;
    ctype c3 = -A.c1 * B.c5 + A.c2 * B.c4 - A.c4 * B.c2 + A.c5 * B.c1 - A.c6 * B.c7 -
               A.c7 * B.c6;
    ctype c4 = A.c1 * B.c7 + A.c7 * B.c1;
    ctype c5 = A.c2 * B.c7 + A.c7 * B.c2;
    ctype c6 = A.c3 * B.c7 + A.c7 * B.c3;

    ctype c7 = 0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

// fat_dot product (including the scalar parts)
//
// inner(A,B) := sum_(r=0,s=0)^n gr( gr(A)_r * gr(B)_s )_|s-r|
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<T> fat_dot(MVec3d<T> const& A, MVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3 - A.c4 * B.c4 -
               A.c5 * B.c5 - A.c6 * B.c6 - A.c7 * B.c7;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c6 + A.c3 * B.c5 - A.c5 * B.c3 +
               A.c6 * B.c2 - A.c4 * B.c7 - A.c7 * B.c4;
    ctype c2 = A.c0 * B.c2 + A.c2 * B.c0 - A.c3 * B.c4 + A.c1 * B.c6 - A.c6 * B.c1 +
               A.c4 * B.c3 - A.c5 * B.c7 - A.c7 * B.c5;
    ctype c3 = A.c0 * B.c3 + A.c3 * B.c0 - A.c1 * B.c5 + A.c2 * B.c4 - A.c4 * B.c2 +
               A.c5 * B.c1 - A.c6 * B.c7 - A.c7 * B.c6;
    ctype c4 = A.c0 * B.c4 + A.c4 * B.c0 + A.c1 * B.c7 + A.c7 * B.c1;
    ctype c5 = A.c0 * B.c5 + A.c5 * B.c0 + A.c2 * B.c7 + A.c7 * B.c2;
    ctype c6 = A.c0 * B.c6 + A.c6 * B.c0 + A.c3 * B.c7 + A.c7 * B.c3;
    ctype c7 = A.c0 * B.c7 + A.c7 * B.c0;
    return MVec3d<ctype>(c0, c1, c2, c3, c4, c5, c6, c7);
}

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


// return commutator product cmt(A,B) of two bivectors A and B (= a bivector)
// cmt(A,B) = 0.5*(AB-BA) = gr2(A * B)
// the commutator product is antisymmetric, i.e. it is zero when a bivector is
// multiplied by itself, i.e. in that case only the dot product remains
// as the symmetric part
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> cmt(BiVec3d<T> const& A,
                                                       BiVec3d<U> const& B)
{
    // this implementation is only valid in an orthonormal basis
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(A.z * B.y - A.y * B.z, A.x * B.z - A.z * B.x,
                          A.y * B.x - A.x * B.y);
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product A*B for fully populated 3d multivector
// Expensive! - Don't use if you don't have to! (64x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors:
// a * b = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)  (vector*vector = scalar + bivector)
// A * b = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)  (bivector*vector = vector + trivector)
// a * B = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)  (vector*bivector = vector + trivector)
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

// geometric product A * B of a multivector A from the even subalgebra (3d case)
// with a bivector B from the right
// even grade multivector * bivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(-A.c1 * B.x - A.c2 * B.y - A.c3 * B.z),
                           BiVec3d<ctype>(A.c0 * B.x - A.c2 * B.z + A.c3 * B.y,
                                          A.c0 * B.y + A.c1 * B.z - A.c3 * B.x,
                                          A.c0 * B.z - A.c1 * B.y + A.c2 * B.x));
}

// geometric product A * B of a bivector A with an even grade multivector B
// from the right
// bivector * even grade multivector => returns an even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& A,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(-A.x * B.c1 - A.y * B.c2 - A.z * B.c3),
                           BiVec3d<ctype>(A.x * B.c0 - A.y * B.c3 + A.z * B.c2,
                                          A.x * B.c3 + A.y * B.c0 - A.z * B.c1,
                                          -A.x * B.c2 + A.y * B.c1 + A.z * B.c0));
}

// geometric product A * b of an even grade multivector A with a vector b
// from the right
// even grade multivector * vector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<T> const& A,
                                                              Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(A.c0 * b.x - A.c2 * b.z + A.c3 * b.y,
                                        A.c0 * b.y + A.c1 * b.z - A.c3 * b.x,
                                        A.c0 * b.z - A.c1 * b.y + A.c2 * b.x),
                           PScalar3d<ctype>(A.c1 * b.x + A.c2 * b.y + A.c3 * b.z));
}

// geometric product a * B of a vector a and an even grade multivector B from the left
// vector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(a.x * B.c0 - a.y * B.c3 + a.z * B.c2,
                                        a.x * B.c3 + a.y * B.c0 - a.z * B.c1,
                                        -a.x * B.c2 + a.y * B.c1 + a.z * B.c0),
                           PScalar3d<ctype>(a.x * B.c1 + a.y * B.c2 + a.z * B.c3));
}

// geometric product A * B of an uneven grade multivector A with a bivector B
// uneven grade multivector * bivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_U<T> const& A,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-A.c1 * B.z + A.c2 * B.y - A.c3 * B.x,
                                        A.c0 * B.z - A.c2 * B.x - A.c3 * B.y,
                                        -A.c0 * B.y + A.c1 * B.x - A.c3 * B.z),
                           PScalar3d<ctype>(A.c0 * B.x + A.c1 * B.y + A.c2 * B.z));
}

// geometric product A * B of a bivector A with an uneven grade multivector B
// bivector * uneven grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(BiVec3d<T> const& A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-A.x * B.c3 - A.y * B.c2 + A.z * B.c1,
                                        +A.x * B.c2 - A.y * B.c3 - A.z * B.c0,
                                        -A.x * B.c1 + A.y * B.c0 - A.z * B.c3),
                           PScalar3d<ctype>(A.x * B.c0 + A.y * B.c1 + A.z * B.c2));
}

// geometric product A * B between two bivectors
// A * B = gr0(A * B) + gr2(A * B)
//
// the full geometric bivector product only exists in >= 4d spaces:
// A * B = dot(A,B) + cmt(A,B) + wdg(A,B) = gr0(A * B) + gr2(A * B) + gr4(A * B)
// In 3D we don't have gr4(A * B) and thus only the terms up to grade 3 remain.
// The bivector product AxB = cmt(A,B) = 0.5*(ab-ba)is called the commutator product.
//
// A * B = dot(A,B) + cmt(A,B) + wdg(A,B)  (in 4d and higher dimensional spaces)
// A * B = dot(A,B) + cmt(A,B)             (in 3d)
//
// => bivector * bivector = scalar + bivector = even grade multivector (in 3d)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(BiVec3d<T> const& a,
                                                              BiVec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(dot(a, b)), cmt(a, b));
}

// geometric product A*b for a bivector and a vector
// A * b = dot(A,b) + wdg(A,b) = gr1(A * b) + gr3(A * b)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => bivector * vector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(BiVec3d<T> const& A,
                                                              Vec3d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    // return MVec3d_U<ctype>(dot(A, b), wdg(A, b));
    return MVec3d_U<ctype>((A >> b), wdg(A, b));
}

// geometric product a * B for a vector a and a bivector B
// a * B = dot(a,B) + wdg(a,B) = gr1(a * B) + gr3(a * B)
//
// HINT: if a full 3d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
//
// => vector * bivector = vector + trivector (= uneven multivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(Vec3d<T> const& a,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // return MVec3d_U<ctype>(dot(a, B), wdg(a, B));
    return MVec3d_U<ctype>((a << B), wdg(a, B));
}

// geometric product a * b of two vectors
// a * b = dot(a,b) + wdg(a,b) = gr0(a * b) + gr2(a * b)
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
    return MVec3d_E<ctype>(Scalar3d<ctype>(dot(a, b)), wdg(a, b));
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
// to the even grade multivector B
// trivector * even grade multivector => uneven grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                              MVec3d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) *
           MVec3d_U<ctype>(Vec3d<ctype>(-B.c1, -B.c2, -B.c3), PScalar3d<ctype>(B.c0));
}

// geometric product A * B of an even multivector A multiplied from the right
// by the trivector B (=3d pseudoscalar)
// even grade multivector * trivector => uneven multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_U<std::common_type_t<T, U>> operator*(MVec3d_E<U> const& A,
                                                              PScalar3d<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_U<ctype>(Vec3d<ctype>(-A.c1, -A.c2, -A.c3), PScalar3d<ctype>(A.c0)) *
           ctype(B);
}

// geometric product A * B of a trivector A (=3d pseudoscalar) multiplied from the left
// to the uneven grade multivector B
// trivector * uneven grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(PScalar3d<T> A,
                                                              MVec3d_U<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(A) *
           MVec3d_E<ctype>(Scalar3d<ctype>(-B.c3), BiVec3d<ctype>(B.c0, B.c1, B.c2));
}

// geometric product A * B of an uneven grade multivector A multiplied from the right
// by the trivector B (=3d pseudoscalar)
// uneven grade multivector * trivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d_E<std::common_type_t<T, U>> operator*(MVec3d_U<U> const& A,
                                                              PScalar3d<T> B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d_E<ctype>(Scalar3d<ctype>(-A.c3), BiVec3d<ctype>(A.c0, A.c1, A.c2)) *
           ctype(B);
}


// (geometric) product s * B of a scalar s multiplied from the left to the bivector B
// scalar * bivector => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                             BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * BiVec3d<ctype>(B.x, B.y, B.z);
}

// (geometric) product B * s of bivector B multiplied by a scalar s from the right
// bivector * scalar => bivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr BiVec3d<std::common_type_t<T, U>> operator*(BiVec3d<T> const& B,
                                                             Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return BiVec3d<ctype>(B.x, B.y, B.z) * ctype(s);
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

// (geometric) product s * v of a scalar s multiplied from the left to the vector v
// scalar * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(Scalar3d<T> s,
                                                           Vec3d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * Vec3d<ctype>(v.x, v.y, v.z);
}

// (geometric) product v * s of a vector v multiplied with a scalar s from the right
// vector * scalar => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> operator*(Vec3d<T> const& v,
                                                           Scalar3d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Vec3d<ctype>(v.x, v.y, v.z) * ctype(s);
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
        throw std::runtime_error("scalar norm too small for inversion" +
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
        throw std::runtime_error("vector norm too small for inversion" +
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
        throw std::runtime_error("bivector norm too small for inversion" +
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
        throw std::runtime_error("bivector norm too small for inversion" +
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
            "norm of even grade multivector too small for inversion" +
            std::to_string(sq_n) + "\n");
    }
#endif
    T inv = T(1.0) / sq_n;
    return MVec3d_E<T>(rev(E) * inv);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> inv(MVec3d<T> const& v)
{
    T m_conjm = gr0(v * conj(v));

#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (std::abs(m_conjm) < std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("multivector norm too small for inversion " +
                                 std::to_string(m_conjm) + "\n");
        // example: MVec2D(1,1,1,1) is not invertible
    }
#endif
    T inv = T(1.0) / m_conjm;
    return MVec3d<T>(conj(v) * inv);
    // ATTENTION: there is a left and a right inverse (see paper of Hitzer, Sangwine)
}


////////////////////////////////////////////////////////////////////////////////
// angle operations for vectors and multivectors of the even subalgebra
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
            "vector norm product too small for calculation of angle" +
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
//                              "of angle" +
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
            "vector norm product too small for calculation of angle" +
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
        throw std::runtime_error("norm product too small for calculation of angle" +
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
        throw std::runtime_error("norm product too small for calculation of angle" +
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
inline constexpr BiVec3d<std::common_type_t<T, U>> rotate(BiVec3d<T> const& v,
                                                          MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;

    // MVec3d_E<ctype> rr = rev(rotor);
    // MVec3d_E<ctype> tmp = rotor * v;
    // MVec3d_E<ctype> res = tmp * rr;

    // scalar part of res is 0 due to symmetric product  rotor * v * rev(rotor)
    //
    // optimization potential for sandwich product by replacing the second product
    // with a specific operation that skips the calculation of the scalar part
    // which will be zero anyway
    // return BiVec3d<ctype>(gr2<ctype>(rotor * v * rev(rotor)));

    // optimized version:
    MVec3d_E<ctype> rr = rev(rotor);
    MVec3d_E<ctype> tmp = rotor * v;
    // formular from operator*(MVec3d_E<T>, MVec3d_E<U>) - only bivector part
    return BiVec3d<ctype>(
        tmp.c0 * rr.c1 + tmp.c1 * rr.c0 - tmp.c2 * rr.c3 + tmp.c3 * rr.c2,
        tmp.c0 * rr.c2 + tmp.c1 * rr.c3 + tmp.c2 * rr.c0 - tmp.c3 * rr.c1,
        tmp.c0 * rr.c3 - tmp.c1 * rr.c2 + tmp.c2 * rr.c1 + tmp.c3 * rr.c0);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr MVec3d<std::common_type_t<T, U>> rotate(MVec3d<T> const& v,
                                                         MVec3d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec3d<ctype>(rotor * v * rev(rotor));
}

////////////////////////////////////////////////////////////////////////////////
// 3d duality operations
// (the concept of dual is defined w.r.t. the geometric product)
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B as subspace of R^3 then
// dual(M) represents the subspace orthorgonal to B
// => return the dual(M) of the multivector M

#if defined(_HD_GA_HESTENES_DORAN_LASENBY_DUALITY)
////////////////////////////////////////////////////////////////////////////////
// duality (as defined by Hestenes or by Doran, Lasenby in "GA for physicists"):
// (same subspace as for Macdonld's definition, but other sign)
// (will have influence on all formulae concerning duality)
////////////////////////////////////////////////////////////////////////////////
//
// dual(A) = I*A
//
// I_3d * 1 = e1^e2^e3 * 1 = e1^e2^e3
//
// I_3d * e1 = e1^e2^e3 * e1 = e_1231 =  e_1123 = e_23         = e_2^e_3
// I_3d * e2 = e1^e2^e3 * e2 = e_1232 = -e_1223 = -e_13 = e_31 = e_3^e_1
// I_3d * e3 = e1^e2^e3 * e3 = e_1233 =  e_12                  = e_1^e_2
//
// I_3d * e2^e3 = e1^e2^e3 * e2^e3 = e_12323 = -e_1            = -e1
// I_3d * e3^e1 = e1^e2^e3 * e3^e1 = e_12331 = -e_11233 = -e_2 = -e2
// I_3d * e1^e2 = e1^e2^e3 * e1^e2 = e_12312 = -e_11223 = -e_3 = -e3
//
// I_3d * e1^e2^e3 = e1^e2^e3 * e1^e2^e3 = e_123123 = -1

// using this duality definition, following identities hold for the cross-product:
//
//    a x b  = -dual(a ^ b) = -I_3d*(a ^ b) = -(a ^ b)*I_3d =  (a ^ b)*rev(I_3d)
//    a ^ b  =  dual(a x b) =  I_3d*(a x b) =  (a x b)*I_3d = -(a x b)*rev(I_3d)
//
//    the last 3 columns in both lines are valid independent of the duality definition!
//    the sign issue only occurs due to the duality definition itself

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> dual(PScalar3d<T> ps)
{
    // dual(A) = I*A
    // e123 * (ps * e123) = -ps
    return Scalar3d<T>(-T(ps));
}

// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars, this can only be avoided, when the scalar
// type is uniquely defined for the corresponding algebra
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> dual(Scalar3d<T> s)
{
    // dual(A) = I*A
    // e123 * (s) = s * e123
    return PScalar3d<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> dual(Vec3d<T> const& v)
{
    // dual(A) = I*A
    // e123 * (v.x * e1  + v.y * e2  + v.z * e3)
    //      =  v.x * e23 + v.y * e31 + v.z * e12
    return BiVec3d<T>(v.x, v.y, v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> dual(BiVec3d<T> const& B)
{
    // dual(A) = I*A
    // e123 * (  b.x * e23 + b.y * e31 + b.z * e12)
    //      =  - b.x * e1  - b.y * e2  - b.z * e3
    return Vec3d<T>(-B.x, -B.y, -B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> dual(MVec3d_E<T> const& M)
{
    // dual(A) = I*A
    // e123 * (s + b.x * e23 + b.y * e31 + b.z * e12)
    //      =    - b.x * e1  - b.y * e2  - b.z * e3 + s * e123
    return MVec3d_U<T>(-M.c1, -M.c2, -M.c3, M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> dual(MVec3d_U<T> const& M)
{
    // dual(A) = I*A
    // e123 * (      v.x * e1  + v.y * e2  + v.z * e3 + ps * e123)
    //      =  -ps + v.x * e23 + v.y * e31 + v.z * e12
    return MVec3d_E<T>(-M.c3, M.c0, M.c1, M.c2);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> dual(MVec3d<T> const& M)
{
    // dual(A) = I*A
    // e123 * (  s + v.x * e1  + v.y * e2  + v.z * e3
    //             + b.x * e23 + b.y * e31 + b.z * e12 + ps * e123)
    //      =  -ps - b.x * e1  - b.y * e2  - b.z * e3
    //             + v.x * e23 + v.y * e31 + v.z * e12 +  s * e123
    return MVec3d<T>(-M.c7, -M.c4, -M.c5, -M.c6, M.c1, M.c2, M.c3, M.c0);
}

#else
////////////////////////////////////////////////////////////////////////////////
// duality (as defined in Macdonald, "Linear and geometric algebra"):
////////////////////////////////////////////////////////////////////////////////
//
// dual(A) = A/I = A*I^(-1) = A*rev(I)
//
// 1 * rev(I_3d) = 1 * e3^e2^e1 = e3^e2^e1 = e_321 = -e_123 = -e1^e2^e3
//
// e1 * rev(I_3d) = e1 * e3^e2^e1 = e_1321 =  e_1132 = -e_23 = -e2^e3
// e2 * rev(I_3d) = e2 * e3^e2^e1 = e_2321 = -e_2231 = -e_31 = -e3^e1
// e3 * rev(I_3d) = e3 * e3^e2^e1 = e_3321 = -e_3312 = -e_12 = -e1^e2
//
// e2^e3 * rev(I_3d) = e2^e3 * e3^e2^e1 = e_23321 = e_1           = e1
// e3^e1 * rev(I_3d) = e3^e1 * e3^e2^e1 = e_31321 = e_33112 = e_2 = e2
// e1^e2 * rev(I_3d) = e1^e2 * e3^e2^e1 = e_12321 = e_11223 = e_3 = e3
//
// e1^e2^e3 * rev(I_3d) = e1^e2^e3 * e3^e2^e1 = e_123321= 1

// using this duality definition, following duality properties hold
// (A. Macdonald, p. 110):
//
// a) dual(aA) = a dual(A)
// b) dual(A + B) = dual(A) + dual(B)
// c) dual(dual(A)) = (-1)^(n*(n-1)/2) A   (with n as dimension of the (sub)space)
// d) |dual(B)| = |B|
// e) if B is a j-blade then dual(B) is an (n-j)-blade
// f) if A is a j-vector then dual(A) is an (n-j)-vector
//    (remember: a j-vector is a sum of j-blades, which are outer products)

// using this duality definition, following identities hold for the cross-product:
//
//    a x b  =  dual(a ^ b) =  (a ^ b)*rev(I_3d) = -I_3d*(a ^ b) = -(a ^ b)*I_3d
//    a ^ b  = -dual(a x b) = -(a x b)*rev(I_3d) =  I_3d*(a x b) =  (a x b)*I_3d
//
//    the last 3 columns in both lines are valid independent of the duality definition!
//    the sign issue only occurs due to the duality definition itself

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar3d<T> dual(PScalar3d<T> ps)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    // (ps * e123) * e321 = ps
    return Scalar3d<T>(T(ps));
}

// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars, this can only be avoided, when the scalar
// type is uniquely defined for the corresponding algebra
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> dual(Scalar3d<T> s)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    // (s) * e321 = -s * e123
    return PScalar3d<T>(-T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> dual(Vec3d<T> const& v)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   ( v.x * e1  + v.y * e2  + v.z * e3) * e321
    // =  -v.x * e23 - v.y * e31 - v.z * e12
    return BiVec3d<T>(-v.x, -v.y, -v.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec3d<T> dual(BiVec3d<T> const& B)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (b.x * e23 + b.y * e31 + b.z * e12) * e321
    // =  b.x * e1  + b.y * e2  + b.z * e3
    return Vec3d<T>(B.x, B.y, B.z);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_U<T> dual(MVec3d_E<T> const& M)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (s + b.x * e23 + b.y * e31 + b.z * e12 ) * e321
    // =      b.x * e1  + b.y * e2  + b.z * e3 - s * e123
    return MVec3d_U<T>(M.c1, M.c2, M.c3, -M.c0);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d_E<T> dual(MVec3d_U<T> const& M)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //   (     v.x * e1  + v.y * e2  + v.z * e3 + ps * e123) * e321
    // =  ps - v.x * e23 - v.y * e31 - v.z * e12
    return MVec3d_E<T>(M.c3, -M.c0, -M.c1, -M.c2);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec3d<T> dual(MVec3d<T> const& M)
{
    // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //  (  s + v.x * e1  + v.y * e2  + v.z * e3
    //       + b.x * e23 + b.y * e31 + b.z * e12 + ps * e123) * e321
    // =  ps + b.x * e1  + b.y * e2  + b.z * e3
    //       - v.x * e23 - v.y * e31 - v.z * e12 -  s * e123
    return MVec3d<T>(M.c7, M.c4, M.c5, M.c6, -M.c1, -M.c2, -M.c3, -M.c0);
}
#endif


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
// left complement:  l_cmpl(u) ^ u  = I_3d = e1^e2^e3
// right complement: u ^ r_cmpl(u)  = I_3d = e1^e2^e3
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_cmpl(u), r_cmpl(u) => cmpl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// complements have different signs

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar3d<T> cmpl(Scalar3d<T> s)
{
    // u ^ cmpl(u) = e1^e2^e3
    // u = s 1:
    //     s ^ cmpl(u) = e1^e2^e3 => cmpl(u) = s e1^e2^e3
    return PScalar3d<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr BiVec3d<T> cmpl(Vec3d<T> const& v)
{
    // u ^ compl(u) = e1^e2^e3
    // u = v.x e1 + v.y e2 + v.z e3:
    //     u ^ cmpl(u) = e1^e2^e3
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
    //     u ^ cmpl(u) = e1^e2^e3
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
    //     u ^ cmpl(u) = e1^e2^e3 => cmpl(u) = ps 1
    return Scalar3d<T>(T(ps));
}


////////////////////////////////////////////////////////////////////////////////
// Vec3d<T> and BiVec3d<T> projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of a vector v1 onto vector v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v1,
                                                              Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * Vec3d<ctype>(inv(v2));
}

// projection of v1 onto v2 (v2 must already be normalized to nrm(v2) == 1)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
project_onto_normalized(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * v2; // v2 is already its own reverse
}

// projection of a vector v1 onto a bivector v2
// v_parallel = dot(v1,v2) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> project_onto(Vec3d<T> const& v,
                                                              BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // // Vec3d<ctype> a = dot(v, B);
    // Vec3d<ctype> a = (v << B);
    // BiVec3d<ctype> Bi = inv(B);
    // // use the formular equivalent to the geometric product to save computational cost
    // // a * Bi = dot(a,Bi) + wdg(a,Bi)
    // // v_parallel = gr1(a * Bi) = dot(a,Bi)
    // // return Vec3d<ctype>(dot(a, Bi));
    // return Vec3d<ctype>(a << Bi);

    return Vec3d<ctype>((v << inv(B)) << B);
    // return Vec3d<ctype>(gr1((v << inv(B)) * B));
}

// projection of a vector v1 onto a normalized bivector v2
// u_parallel = gr1(dot(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
project_onto_normalized(Vec3d<T> const& v, BiVec3d<U> const& B)
{
    // requires B to be normalized

    using ctype = std::common_type_t<T, U>;
    // // Vec3d<ctype> a = dot(v, B);
    // Vec3d<ctype> a = (v << B);
    // // up to the sign B already is it's own inverse
    // BiVec3d<ctype> Bi = -B;
    // // use the formular equivalent to the geometric product to save computational cost
    // // a * Bi = dot(a,Bi) + wdg(a,Bi)
    // // v_parallel = gr1(a * Bi) = dot(a,Bi)
    // // return Vec3d<ctype>(dot(a, Bi));
    // return Vec3d<ctype>(a << Bi);

    // up to the sign B already is it's own inverse
    return Vec3d<ctype>((v << (-B)) << B);
}

// rejection of vector v1 from a vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v1,
                                                             Vec3d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    BiVec3d<ctype> B = wdg(v1, v2);
    Vec3d<ctype> v2_inv = inv(v2);
    // use the formular equivalent to the geometric product to save computational cost
    // B * b_inv = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(B * b_inv) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    // return Vec3d<ctype>(dot(B, v2_inv));
    return Vec3d<ctype>(B >> v2_inv);
}

// rejection of vector v1 from a normalized vector v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
reject_from_normalized(Vec3d<T> const& v1, Vec3d<U> const& v2)
{
    // requires v2 to be normalized
    using ctype = std::common_type_t<T, U>;
    BiVec3d<ctype> B = wdg(v1, v2);
    Vec3d<ctype> v2_inv = v2; // v2 is its own inverse, if normalized
    // use the formular equivalent to the geometric product to save computational cost
    // B * b_inv = dot(B,b_inv) + wdg(A,bi)
    // v_perp = gr1(B * b_inv) = dot(B,b_inv)
    // (the trivector part is zero, because v2 is part of the bivector in the product)
    // return Vec3d<ctype>(dot(B, v2_inv));
    return Vec3d<ctype>(B >> v2_inv);
}

// rejection of vector v1 from a bivector v2
// u_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reject_from(Vec3d<T> const& v,
                                                             BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    PScalar3d<ctype> A = wdg(v, B);
    BiVec3d<ctype> Bi = inv(B);
    // trivector * bivector = vector
    return A * Bi;
}

// rejection of vector v1 from a normalized bivector v2
// u_perp = wdg(v1,v2) * inv(v2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>>
reject_from_normalized(Vec3d<T> const& v, BiVec3d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    PScalar3d<ctype> a = wdg(v, B);
    // up to the sign v2 already is it's own inverse
    BiVec3d<ctype> Bi = -B;
    // trivector * bivector = vector (derived from full geometric product to save
    // costs)
    return Vec3d<ctype>(-a * Bi.x, -a * Bi.y, -a * Bi.z);
}

// reflect a vector v on a hyperplane B orthogonal to vector nB
//
// hyperplane: a n-1 dimensional subspace in a space of dimension n (a line in 2d space)
// orthogonal to vector b: the hyperplane is dual to b (i.e. a one dimensional subspace)
//
// HINT: choose b * B = I_3d  =>  B = b * I_3d  (for normalized b)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Vec3d<std::common_type_t<T, U>> reflect_on_hyp(Vec3d<T> const& v,
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
    basis.emplace_back(normalize(reject_from_normalized(v, u_unitized)));
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
    basis.emplace_back(normalize(reject_from_normalized(v, u_unitized)));
    basis.emplace_back(normalize(reject_from(w, wdg(u, v))));
    return basis;
}

} // namespace hd::ga::ega