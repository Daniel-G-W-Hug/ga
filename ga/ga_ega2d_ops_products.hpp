#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_ega2d_ops_basics.hpp" // ega2d ops basics


namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega2d product operations:
//
// - dot()              -> dot product
// - wdg()              -> wedge product
// - rwdg()             -> regressive wedge product
// - operator<<()       -> left contraction
// - operator>>()       -> right contraction
// - cmt()              -> commutator product
// - operator*()        -> geometric product
// - inv()              -> inversion operation
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot product (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume e1^2 = +1, e2^2 = +1
//
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> dot(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> dot(PScalar2d<T> ps1, PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// return dot-product of two vectors in G<2,0,0>
// dot(v1,v2) = nrm(v1)*nrm(v2)*cos(angle)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> dot(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    // definition: dot(v1, v2) = (v1)^T g_12 v2 with the metric g_12
    // this assumes an orthonormal basis with e1^2 = +1, e2^2 = +1
    // as diagonal elements of g_12
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.x + v1.y * v2.y);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> dot(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// wedge product (= outer product)
////////////////////////////////////////////////////////////////////////////////

// wedge product between two fully populated multivectors A and B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0;
    ctype c2 = A.c0 * B.c2 + A.c2 * B.c0;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product between multivectors A and even grade multivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& A, MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c1 * B.c0;
    ctype c2 = A.c2 * B.c0;
    ctype c3 = A.c0 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product between even grade multivector A and multivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c0 * B.c1;
    ctype c2 = A.c0 * B.c2;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product between multivector M and pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& M, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(M.c0 * ctype(ps));
}

// wedge product between pseudoscalar ps and multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(PScalar2d<T> ps, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * M.c0);
}

// wedge product M ^ v between multivector M and vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& M, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = 0.0;
    ctype c1 = M.c0 * v.x;
    ctype c2 = M.c0 * v.y;
    ctype c3 = M.c1 * v.y - M.c2 * v.x;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product v ^ M between vector v and multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = 0.0;
    ctype c1 = v.x * M.c0;
    ctype c2 = v.y * M.c0;
    ctype c3 = v.x * M.c2 - v.y * M.c1;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// wedge product between multivector M and scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// wedge product between scalar s and multivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// wedge product between two even grade multivectors A and B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& A,
                                                 MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// wedge product between even grade multivector M and pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& M, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(M.c0 * ctype(ps));
}

// wedge product between pseudoscalar ps and even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(PScalar2d<T> ps, MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * M.c0);
}

// wedge product between even grade multivector M and vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& M, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(M.c0 * v.x, M.c0 * v.y);
}

// wedge product between vector v and even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> wdg(Vec2d<T> const& v, MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.x * M.c0, v.y * M.c0);
}

// wedge product between even grade multivector M and scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> wdg(MVec2d_E<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// wedge product between scalar s and even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> wdg(Scalar2d<T> s, MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// wedge product of two pseudoscalars
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2d<T>,
                                                 [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// wedge product of a vector with a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> wdg([[maybe_unused]] Vec2d<T> const&,
                                                 [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// wedge product of a pseudoscalar with a vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> wdg([[maybe_unused]] PScalar2d<T>,
                                                 [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * ctype(s));
}

// wedge product with a pscalar (returns a scaled pscalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(s) * ctype(ps));
}

// wedge product (returns a bivector)
// wdg(v1,v2) = |v1| |v2| sin(theta)
// where theta: -pi <= theta <= pi (different to definition of angle for dot product!)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    return PScalar2d<std::common_type_t<T, U>>(v1.x * v2.y - v1.y * v2.x);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> wdg(Vec2d<T> const& v, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// wedge product with one scalar (returns a scaled vector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// wedge product between two scalars (returns a scalar)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> wdg(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// regressive wedge product (= outer product between complements)
// as defined by E. Lengyel in "Projective geometric algebra illuminated"
// independent of the geometric product (as is the case for the dual), but just depending
// on the outer product (wdg) as well as the complement of the multivector
// (i.e. this definition is does NOT connect directly to the geometric product, but to the
// outer product, exclusively)
//
// rwdg(ul, ur) = lcmpl(wdg(rcmpl(ul),rcmpl(ur))) = lcmpl(rcmpl(ul)^rcmpl(ur))
//              = rcmpl(wdg(lcmpl(ul),lcmpl(ur))) = rcmpl(lcmpl(ul)^lcmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////

// regressive wedge product extended to a fully populated multivectors
// => returns a multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> rwdg(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    ctype c1 = A.c1 * B.c3 + A.c3 * B.c1;
    ctype c2 = A.c2 * B.c3 + A.c3 * B.c2;
    ctype c3 = A.c3 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// regressive wedge product between to pseudoscalars
// => returns a scaled pseudoscalar (=bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps1, PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// regressive wedge product between a pseudoscalar ps (=bivector) and a vector v
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2d<ctype>(v.x, v.y);
}

// regressive wedge product between a vector v and a pseudoscalar ps (=bivector)
// => returns a vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> rwdg(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.x, v.y) * ctype(ps);
}

// regressive wedge product between a pseudoscalar ps (=bivector) and a scalar s
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps) * ctype(s));
}

// regressive wedge product between a scalar s and a pseudoscalar ps (=bivector)
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s) * ctype(ps));
}

// regressive wedge product between a vector v1 and a vector v2
// => returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(v1.x * v2.y - v1.y * v2.x);
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Vec2d<T> const&,
                                                  [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2d<T>,
                                                  [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>{0.0};
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> rwdg([[maybe_unused]] Scalar2d<T>,
                                                  [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>{0.0};
}


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B
//
// The resulting object is A lies in B and is perpendicular to A
//
// L. Dorst: The contraction A << B of an a-blade A onto a b-blade B is a sub-blade
// of B of grade b-a which is perpendicular to A, and linear in both arguments
////////////////////////////////////////////////////////////////////////////////

// left contraction (A << B) - "A contracted onto B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& A,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype c1 = A.c0 * B.c1 + A.c2 * B.c3;
    ctype c2 = A.c0 * B.c2 - A.c1 * B.c3;
    ctype c3 = A.c0 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (A << B) - "A contracted onto B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& A,
                                                      MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c3 * B.c1;
    ctype c1 = A.c2 * B.c1;
    ctype c2 = -A.c1 * B.c1;
    ctype c3 = A.c0 * B.c1;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (A << B) - "A contracted onto B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& A,
                                                      MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c3;
    ctype c1 = A.c0 * B.c1;
    ctype c2 = A.c0 * B.c2;
    ctype c3 = A.c0 * B.c3;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (M << ps) of mulivector A onto pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<U> const& M, PScalar2d<T> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c3 * ps;
    ctype c1 = M.c2 * ps;
    ctype c2 = -M.c1 * ps;
    ctype c3 = M.c0 * ps;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (ps << M) of a pseudoscalar ps onto mulivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(PScalar2d<T> ps,
                                                        MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ps * M.c3);
}

// left contraction (M << v) of multivector M onto vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& M,
                                                      Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c1 * v.x + M.c2 * v.y;
    ctype c1 = M.c0 * v.x;
    ctype c2 = M.c0 * v.y;
    ctype c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (v << M) of vector v onto multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v,
                                                      MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = v.x * M.c1 + v.y * M.c2;
    ctype c1 = v.y * M.c3;
    ctype c2 = -v.x * M.c3;
    ctype c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// left contraction (M << s) of multivector M onto scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c0 * s);
}

// left contraction (s << M) of scalar s onto multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;

    return ctype(s) * M;
}

// left contraction (A << B) - "A contracted onto B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& A,
                                                        MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1;
    ctype c1 = A.c0 * B.c1;
    return MVec2d_E<ctype>(c0, c1);
}

// left contraction (M << ps) of even grade mulivector M with pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& M,
                                                        PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c1 * ps;
    ctype c1 = M.c0 * ps;
    return MVec2d_E<ctype>(c0, c1);
}

// left contraction (ps << M) of pseudoscalar ps with an even grade mulivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(PScalar2d<T> ps,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ps * M.c1);
}

// left contraction (M << v) - "M contracted onto v"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& M,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(M.c0 * v.x, M.c0 * v.y);
}

// left contraction (v << M) - "v contracted onto M"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.y * M.c1, -v.x * M.c1);
}

// left contraction (M << s) of even grade mulivector M with scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(MVec2d_E<T> const& M,
                                                        Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c0 * s);
}

// left contraction (s << M) of a scalar s with an even grade mulivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator<<(Scalar2d<T> s,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * M;
}

// left contraction - pseudoscalar ps1 contracted onto pseudoscalar ps2
// returns a scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(PScalar2d<T> ps1,
                                                        PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ps1 * ps2);
}

// left contraction -  pseudoscalar contracted onto vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar2d<T>,
                                                        [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// left contraction (v << ps) - vector v contracted onto pseudoscalar ps
// (identical with the geometric product v * ps for this case)
// (identical with cmt(ps, v))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.y, -v.x) * ctype(ps);
}

// left contraction - pseudoscalar contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<([[maybe_unused]] PScalar2d<T>,
                                                        [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// left contraction (s << ps) - scalar s contracted onto pseudoscalar ps
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(s) * ctype(ps));
}

// left contraction (v1 << v2) - vector v1 contracted onto vector v2
// (=identical to scalar product dot(v1,v2) of two vectors and identical to (v1 >> v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(Vec2d<T> const& v1,
                                                        Vec2d<U> const& v2)
{
    return dot(v1, v2);
}

// left contraction - vector contracted onto scalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<([[maybe_unused]] Vec2d<T> const&,
                                                        [[maybe_unused]] Scalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// left contraction (s << v) - scalar s contracted onto vector v
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// left contraction (s1 << s2) - scalar s1 contracted onto scalar s2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator<<(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B
//
// The resulting object lies in A and is perpendicular to B
////////////////////////////////////////////////////////////////////////////////

// right contraction (A >> B)  - "A contracted by B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<T> operator>>(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 + A.c3 * B.c3;
    ctype c1 = A.c1 * B.c0 - A.c3 * B.c2;
    ctype c2 = A.c2 * B.c0 + A.c3 * B.c1;
    ctype c3 = A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (A >> B)  - "A contracted by B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<T> operator>>(MVec2d<T> const& A, MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c3 * B.c1;
    ctype c1 = A.c1 * B.c0;
    ctype c2 = A.c2 * B.c0;
    ctype c3 = A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (A >> B)  - "A contracted by B"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<T> operator>>(MVec2d_E<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c3;
    ctype c1 = -A.c1 * B.c2;
    ctype c2 = A.c1 * B.c1;
    ctype c3 = A.c1 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}
// right contraction (M >> ps) of a multivector M by a pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(MVec2d<T> const& M,
                                                        PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c3 * ps);
}

// right contraction (ps >> M) of a pseudoscalar ps by a multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = ps * M.c3;
    ctype c1 = -ps * M.c2;
    ctype c2 = ps * M.c1;
    ctype c3 = ps * M.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (M >> v) of multivector A contracted by vector v
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<T> operator>>(MVec2d<T> const& M, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = M.c1 * v.x + M.c2 * v.y;
    ctype c1 = -M.c3 * v.y;
    ctype c2 = M.c3 * v.x;
    ctype c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (v >> M) of vector v contracted by multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<T> operator>>(Vec2d<T> const& v, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = v.x * M.c1 + v.y * M.c2;
    ctype c1 = v.x * M.c0;
    ctype c2 = v.y * M.c0;
    ctype c3 = 0.0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// right contraction (M >> s) of multivector M by scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator>>(MVec2d<T> const& M, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// right contraction (s >> M) of scalar s by multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Scalar2d<T> s, MVec2d<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(s * M.c0);
}

// right contraction (A >> B) of an even grade multivector A by even grade multivector B
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& A,
                                                        MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1;
    ctype c1 = A.c1 * B.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// right contraction (M >> ps) of an even grade multivector M by a pseudoscalar ps
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& M,
                                                        PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(M.c1 * ps);
}

// right contraction (ps >> M) of pseudoscalar ps by even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    ctype c0 = ps * M.c1;
    ctype c1 = ps * M.c0;
    return MVec2d_E<ctype>(c0, c1);
}

// left contraction (M >> v) - "M contracted by v"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& M,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-M.c1 * v.y, M.c1 * v.x);
}

// left contraction (v >> M) - "v contracted by M"
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(Vec2d<T> const& v,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.x * M.c0, v.y * M.c0);
}

// right contraction (M >> s) of even grade multivector M by scalar s
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator>>(MVec2d_E<T> const& M,
                                                        Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return M * ctype(s);
}

// right contraction (s >> M) of scalar s by even grade multivector M
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Scalar2d<T> s,
                                                        MVec2d_E<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(s * M.c0);
}

// right contraction - pseudoscalar s1 contracted by pseudoscalar s2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps1,
                                                        PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(ps1) * ctype(ps2));
}

// right contraction (ps >> v) - pseudoscalar ps contracted by vector v
// (identical with the geometric product ps * v for this case)
// (identical to cmt(v,ps))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2d<ctype>(-v.y, v.x);
}

// right contraction -  vector contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>([[maybe_unused]] Vec2d<T> const&,
                                                        [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// right contraction (ps >> s) - pseudoscalar ps contracted by scalar s
// identical with scalar multiplication of the pseudoscalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> operator>>(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * ctype(s));
}

// right contraction - scalar contracted by a pseudoscalar
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar2d<T>,
                                                        [[maybe_unused]] PScalar2d<U>)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// right contraction (v1 >> v2) - vector v2 taken out of a vector v1
// (=identical to scalar product dot(v1,v2) of two vectors and identical to (v1 << v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Vec2d<T> const& v1,
                                                        Vec2d<U> const& v2)
{
    return dot(v1, v2);
}

// right contraction (v >> s) - vector v contracted by scalar s
// identical with scalar multiplication of the vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator>>(Vec2d<U> const& v, Scalar2d<T> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// right contraction - scalar contracted by a vector
// returns 0
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>([[maybe_unused]] Scalar2d<T>,
                                                        [[maybe_unused]] Vec2d<U> const&)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(0.0);
}

// left contraction (s1 >> s2) - scalar s1 contracted by scalar s2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator>>(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> cmt(MVec2d<T> const& A, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(0.0, -A.c2 * B.c3 + A.c3 * B.c2, A.c1 * B.c3 - A.c3 * B.c1,
                         A.c1 * B.c2 - A.c2 * B.c1);
}

// cmt(B,v) = -cmt(v,B)
// identical to (v << ps)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> cmt(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(ctype(ps) * v.y, -ctype(ps) * v.x);
}

// cmt(v,B) = -cmt(B,v)
// identical to (ps >> v)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> cmt(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-v.y * ctype(ps), v.x * ctype(ps));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> cmt(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(v1.x * v2.y - v1.y * v2.x);
}


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////

// geometric product A * B for fully populated 2d multivectors
// Expensive! - Don't use if you don't have to! (16x mul_add)
//
// Use equivalent formulae instead for not fully populated multivectors, e.g.:
// a*b = dot(a,b) + wdg(a,b) = gr0(a*b) + gr2(a*b) (vector * vector = scalar + bivector)
//
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                     MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    // geometric product A * B of two fully populated multivectors
    ctype c0 = A.c0 * B.c0 + A.c1 * B.c1 + A.c2 * B.c2 - A.c3 * B.c3;
    ctype c1 = A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2;
    ctype c2 = A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1;
    ctype c3 = A.c0 * B.c3 + A.c1 * B.c2 - A.c2 * B.c1 + A.c3 * B.c0;
    return MVec2d<ctype>(c0, c1, c2, c3);
}

// geometric product A * B for a multivector A with an even grade multivector B
// multivector * even grade multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                     MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c0 * B.c0 - A.c3 * B.c1, A.c1 * B.c0 - A.c2 * B.c1,
                         A.c1 * B.c1 + A.c2 * B.c0, A.c0 * B.c1 + A.c3 * B.c0);
}

// geometric product A * B for an even grade multivector A with a full 2d multivector B
// from the right
// even grade multivector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                     MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c0 * B.c0 - A.c1 * B.c3, A.c0 * B.c1 + A.c1 * B.c2,
                         A.c0 * B.c2 - A.c1 * B.c1, A.c0 * B.c3 + A.c1 * B.c0);
}

// geometric product A * ps of a multivector A with the pseudoscalar ps (=bivector)
// multivector * pseudoscalar (=bivector) => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(-A.c3, -A.c2, A.c1, A.c0) * ctype(ps);
}

// geometric product ps * B of pseudoscalar ps (=bivector) with the multivector B
// 2d pseudoscalar (=bivector) * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec2d<ctype>(-B.c3, B.c2, -B.c1, B.c0);
}

// geometric product A * v for a multivector A with a vector v
// multivector * vector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A,
                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(A.c1 * v.x + A.c2 * v.y, A.c0 * v.x + A.c3 * v.y,
                         A.c0 * v.y - A.c3 * v.x, A.c1 * v.y - A.c2 * v.x);
}

// geometric product v * B for a vector v with a full 2d multivector B
// vector * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v,
                                                     MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(v.x * B.c1 + v.y * B.c2, v.x * B.c0 - v.y * B.c3,
                         v.x * B.c3 + v.y * B.c0, v.x * B.c2 - v.y * B.c1);
}

// geometric product A * s of a multivector A multiplied from the right by the scalar s
// multivector * scalar => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(MVec2d<T> const& A, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return A * ctype(s);
}

// geometric product s * B of a scalar s multiplied to the multivector B from left
// scalar * multivector => multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s, MVec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

// geometric product A * B for two multivectors from the even subalgebra (2d case)
// even grade multivector * even grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                       MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(A.c0 * B.c0 - A.c1 * B.c1, A.c0 * B.c1 + A.c1 * B.c0);
}

// geometric product A * ps of an even grade multivector A multiplied from the right
// with the 2d pseudoscalar (=bivector) ps
// even grade multivector * 2d pseudoscalar (=bivector) => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                       PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(-A.c1, A.c0) * ctype(ps);
}

// geometric product ps * B of a bivector A multiplied from the left
// to the multivector from the even subalgebra B (MVec2d_E)
// bivector * even grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(PScalar2d<T> ps,
                                                       MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * MVec2d_E<ctype>(-B.c1, B.c0);
}

// geometric product A * v for a multivector from the even subalgebra A with a vector v
// even grade multivector * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                    Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(A.c0 * v.x + A.c1 * v.y, A.c0 * v.y - A.c1 * v.x);
}

// geometric product v * B of a vector v with an even grade multivector B
// vector * even grade multivector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v,
                                                    MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v.x * B.c0 - v.y * B.c1, v.x * B.c1 + v.y * B.c0);
}

// geometric product A * s of an even grade multivector A multiplied with a scalar s
// even grade multivector * scalar => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(MVec2d_E<T> const& A,
                                                       Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return A * ctype(s);
}

// geometric product s * B of a scalar with an even grade multivector B
// scalar * even grade multivector => even grade multivector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(Scalar2d<T> s,
                                                       MVec2d_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * B;
}

// geometric product ps * ps of two pseudoscalars
// pseudoscalar * pseudoscalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps1, PScalar2d<U> ps2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(-ctype(ps1) * ctype(ps2));
}

// geometric product ps * v of a pseudoscalar s with the vector b
// pseudoscalar (=bivector) * vector => vector
// this multiplication turns the vector by -90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(ps) * Vec2d<ctype>(v.y, -v.x);
}

// geometric product v * ps of the vector v multiplied with the pseudoscalar ps
// vector * pseudoscalar (=bivector) => vector
// this multiplication turns the vector by +90° in the plane e1^e2
// (positive angle is from e1 towards e2)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-v.y, v.x) * ctype(ps);
}

// geometric product ps * s of pseudscalar ps and scalar s
// pseudoscalar (=bivector) * scalar => pseudoscalar (=bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> operator*(PScalar2d<T> ps, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(ps) * ctype(s));
}

// geometric product s * ps of scalar s and pseudoscalar ps
// scalar * pseudoscalar (=bivector) => pseudoscalar (=bivector)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr PScalar2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s, PScalar2d<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return PScalar2d<ctype>(ctype(s) * ctype(ps));
}

// geometric product ab for two vectors (returns a multivector of the even subalgebra)
// v1*v2 = dot(v1,v2) + wdg(v1,v2) = gr0(v1*v2) + gr2(v1*v2)
// => vector * vector = scalar + bivector
//
// HINT: if a full 2d multivector is required as result it must be converted explicitly,
//       since C++ does not allow overloading on different return types
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> operator*(Vec2d<T> const& v1,
                                                       Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(dot(v1, v2), wdg(v1, v2));
}

// geometric product v * s of vector v and scalar s
// vector * scalar => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(Vec2d<T> const& v, Scalar2d<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return v * ctype(s);
}

// geometric product s * v of scalar s and vector v
// scalar * vector => vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s, Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(s) * v;
}

// geometric product s1 * s2 of two scalars
// scalar * scalar => scalar
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Scalar2d<std::common_type_t<T, U>> operator*(Scalar2d<T> s1, Scalar2d<U> s2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar2d<ctype>(ctype(s1) * ctype(s2));
}


////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses of scalars, blades and multivectors
// w.r.t. the geometric product:
// for k-blades: A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
////////////////////////////////////////////////////////////////////////////////
// HINT: inv() cannot be constexpr due to the checks for division by zero
//       which might throw

template <typename T>
    requires(std::floating_point<T>)
inline Scalar2d<T> inv(Scalar2d<T> s)
{
    T sq_n = nrm_sq(s);
    hd::ga::detail::check_normalization<T>(sq_n, "scalar");
    T inv = T(1.0) / sq_n;

    return Scalar2d<T>(rev(s) * inv);
}

template <typename T>
    requires(std::floating_point<T>)
inline Vec2d<T> inv(Vec2d<T> const& v)
{
    T sq_n = nrm_sq(v);
    hd::ga::detail::check_normalization<T>(sq_n, "vector");
    T inv = T(1.0) / sq_n;
    return Vec2d<T>(rev(v) * inv);
}

template <typename T>
    requires(std::floating_point<T>)
inline PScalar2d<T> inv(PScalar2d<T> ps)
{
    T sq_n = nrm_sq(ps);
    hd::ga::detail::check_normalization<T>(sq_n, "pseudoscalar");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector

    return PScalar2d<T>(rev(ps) * inv);
}

// return the multiplicative inverse of the complex number (inv(z) = 1/nrm_sq(z)*rev(z))
// with rev(z) being the complex conjugate
template <typename T>
    requires(std::floating_point<T>)
inline MVec2d_E<T> inv(MVec2d_E<T> const& E)
{
    T sq_n = nrm_sq(E);
    hd::ga::detail::check_normalization<T>(sq_n, "even grade multivector");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return MVec2d_E<T>(rev(E) * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord algebras",
// Hitzer, Sangwine, 2016
// left and a right inverse are the same (see paper of Hitzer, Sangwine)
template <typename T>
    requires(std::floating_point<T>)
inline MVec2d<T> inv(MVec2d<T> const& M)
{
    T m_conjm = T(gr0(M * conj(M)));
    hd::ga::detail::check_normalization<T>(std::abs(m_conjm), "multivector");
    T inv = T(1.0) / m_conjm; // inverse of squared norm for a vector
    return MVec2d<T>(conj(M) * inv);
}

} // namespace hd::ga::ega