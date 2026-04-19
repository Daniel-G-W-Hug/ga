#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_foundation.hpp"    // ga library headers and infrastructure
#include "detail/type_t/ga_mvec2d.hpp" // 2D multivector types (includes component types)

namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega2d basic operations:
//
// - gr_inv()                      -> grade inversion
// - rev()                         -> reversion
// - conj()                        -> conjugation
//
// - l_cmpl(), r_cmpl()            -> left and right complements
//
// - nrm_sq(), nrm()               -> norm
//
// - normalize()                   -> normalization (nrm scaled to 1.0)
//
// - l_dual(), r_dual()            -> left and right duals
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for k = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2d<T> gr_inv(Scalar2d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2d<T> gr_inv(Vec2d<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2d<T> gr_inv(PScalar2d<T> ps)
{
    // grade 2: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d_E<T> gr_inv(MVec2d_E<T> const& M)
{
    // grade 0 and 2: no sign change
    return M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d<T> gr_inv(MVec2d<T> const& M)
{
    // grade 0 and 2: no sign change
    // grade 1: sign reversal
    return MVec2d<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)));
}


////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2d<T> rev(Scalar2d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2d<T> rev(Vec2d<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2d<T> rev(PScalar2d<T> ps)
{
    // grade 2: sign reversal
    return -ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d_E<T> rev(MVec2d_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign reversal
    return MVec2d_E<T>(rev(gr0(M)), rev(gr2(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d<T> rev(MVec2d<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign reversal
    return MVec2d<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)));
}


////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation:
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for k = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2d<T> conj(Scalar2d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2d<T> conj(Vec2d<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2d<T> conj(PScalar2d<T> ps)
{
    // grade 2: sign reversal
    return -ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d_E<T> conj(MVec2d_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign reversal
    return MVec2d_E<T>(conj(gr0(M)), conj(gr2(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d<T> conj(MVec2d<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    return MVec2d<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)));
}


////////////////////////////////////////////////////////////////////////////////
// 2d complement operations (impact on basis vector exclusively!)
// (the concept of complement is defined w.r.t. the outer product)
// it depends on the pseudoscalar of the space
// here implemented for I_2d = e1^e2 = e12
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B of the blade u as subspace of R^2 then
// compl(M) represents the subspace orthorgonal to B
// the complement exchanges basis vectors (exclusively, no impact on magnitude),
// which are in the k-blade u with the basis vectors which are NOT contained in the
// k-blade u and are needed to fill the space completely to the corresponding pseudoscalar
//
// left complement:  l_cmpl(u) ^ u  = I_2d = e1^e2  =>  l_cmpl(u) = I_2d * rev(u)
// right complement: u ^ r_cmpl(u)  = I_2d = e1^e2  =>  r_cmpl(u) = rev(u) * I_2d
//
// (get formula on the rhs by multiplying with inv(u) from the right or left respectively)
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_cmpl(u) = r_cmpl(u) = cmpl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// complements have different signs

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2d<T> l_cmpl(Scalar2d<T> s)
{
    // l_cmpl(u) ^ u = e1^e2
    // u = s 1:
    //     l_cmpl(u) ^ u = e1^e2 => l_cmpl(u) = I_2d * rev(s) = s e1^e2
    return PScalar2d<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2d<T> l_cmpl(Vec2d<T> const& v)
{
    // l_cmpl(u) ^ u = e1^e2
    // u = v.x e1 + v.y e2:
    //     l_cmpl(u) ^ u = e1^e2 => cmpl(u) = I_2d * rev(v)
    //     l_cmpl(u) = e12 * (v.x * e1 + v.y * e2) = -v.x * e2 + v.y * e1
    //                                            b=  v.y * e1 - v.x * e2
    //     => l_cmpl(u).x =  v.y e1
    //     => l_cmpl(u).y = -v.x e2
    return Vec2d<T>(v.y, -v.x);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2d<T> l_cmpl(PScalar2d<T> ps)
{
    // l_cmpl(u) ^ u = e1^e2
    // u = ps e1^e2:
    //     l_cmpl(u) ^ u = e1^e2 => l_cmpl(u) = I_2d * rev(ps) = ps 1
    return Scalar2d<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d_E<T> l_cmpl(MVec2d_E<T> const& M)
{
    // use the component complements directly
    return MVec2d_E<T>(l_cmpl(gr2(M)), l_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d<T> l_cmpl(MVec2d<T> const& M)
{
    // use the component complements directly
    return MVec2d<T>(l_cmpl(gr2(M)), l_cmpl(gr1(M)), l_cmpl(gr0(M)));
}


template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2d<T> r_cmpl(Scalar2d<T> s)
{
    // u ^ r_cmpl(u) = e1^e2
    // u = s 1:
    //     u ^ r_cmpl(u) = e1^e2 => r_cmpl(u) = rev(s) * I_2d = s e1^e2
    return PScalar2d<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2d<T> r_cmpl(Vec2d<T> const& v)
{
    // u ^ r_cmpl(u) = e1^e2
    // u = v.x e1 + v.y e2:
    //     u ^ r_cmpl(u) = e1^e2 => r_cmpl(u) = rev(v) * I_2d
    //     r_cmpl(u) = (v.x * e1 + v.y * e2) * e12 =   v.x * e2 - v.y * e1
    //                                             = - v.y * e1 + v.x * e2
    //     => r_cmpl(u).x = -v.y e1
    //     => r_cmpl(u).y =  v.x e2

    return Vec2d<T>(-v.y, v.x);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2d<T> r_cmpl(PScalar2d<T> ps)
{
    // u ^ r_cmpl(u) = e1^e2
    // u = ps e1^e2:
    //     u ^ r_cmpl(u) = e1^e2 => r_cmpl(u) = rev(ps) * I_2d = ps 1
    return Scalar2d<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d_E<T> r_cmpl(MVec2d_E<T> const& M)
{
    // use the component complements directly
    return MVec2d_E<T>(r_cmpl(gr2(M)), r_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d<T> r_cmpl(MVec2d<T> const& M)
{
    // use the component complements directly
    return MVec2d<T>(r_cmpl(gr2(M)), r_cmpl(gr1(M)), r_cmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// nrm_sq(), nrm()
//
// nrm_sq(M)    -> return squared norm
// nrm(M)       -> return norm (Euclidean)
//
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(Scalar2d<T> s)
{
    return T(s) * T(s);
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm(Scalar2d<T> s)
{
    return std::sqrt(nrm_sq(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(Vec2d<T> const& v)
{
    return v.x * v.x + v.y * v.y;
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm(Vec2d<T> const& v)
{
    return std::sqrt(nrm_sq(v));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(PScalar2d<T> ps)
{
    return T(ps) * T(ps);
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm(PScalar2d<T> ps)
{
    return std::sqrt(nrm_sq(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(MVec2d_E<T> const& M)
{
    return nrm_sq(gr0(M)) + nrm_sq(gr2(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm(MVec2d_E<T> const& M)
{
    return std::sqrt(nrm_sq(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(MVec2d<T> const& M)
{
    return nrm_sq(gr0(M)) + nrm_sq(gr1(M)) + nrm_sq(gr2(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm(MVec2d<T> const& M)
{
    return std::sqrt(nrm_sq(M));
}


////////////////////////////////////////////////////////////////////////////////
// normalization operations:
// return an object normalized to nrm(u) == 1.0
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
inline Scalar2d<T> normalize(Scalar2d<T> s)
{
    T m = nrm(s);
    hd::ga::detail::check_normalization<T>(m, "scalar (2d)");
    T inv = T(1.0) / m;
    return Scalar2d<T>(inv * T(s));
}

template <typename T>
    requires(numeric_type<T>)
inline Vec2d<T> normalize(Vec2d<T> const& v)
{
    T m = nrm(v);
    hd::ga::detail::check_normalization<T>(m, "vector (2d)");
    T inv = T(1.0) / m;
    return inv * v;
}

template <typename T>
    requires(numeric_type<T>)
inline PScalar2d<T> normalize(PScalar2d<T> ps)
{
    T m = nrm(ps);
    hd::ga::detail::check_normalization<T>(m, "pseudoscalar (2d)");
    T inv = T(1.0) / m;
    return PScalar2d<T>(inv * T(ps));
}

template <typename T>
    requires(numeric_type<T>)
inline MVec2d_E<T> normalize(MVec2d_E<T> const& M)
{
    T m = nrm(M);
    hd::ga::detail::check_normalization<T>(m, "even-grade multivector (2d)");
    T inv = T(1.0) / m;
    return inv * M;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec2d<T> normalize(MVec2d<T> const& M)
{
    T m = nrm(M);
    hd::ga::detail::check_normalization<T>(m, "multivector (2d)");
    T inv = T(1.0) / m;
    return inv * M;
}

////////////////////////////////////////////////////////////////////////////////
// duality is defined w.r.t. the outer product,
// because the complement is defined w.r.t. the outer product.
// (as defined in Lengyel, "PGA illuminated")
//
// The only difference between the dual and the complement is the fact that
// for the dual the argument is multiplied with the extended metric before
// taking the complement. However, in ega the metric is the identity matrix.
// Thus the dual and the complement are identical.
////////////////////////////////////////////////////////////////////////////////
//
// if M represents the subspace B as subspace of R^2 then
// l_dual(M) and r_dual(M) represent a subspace orthorgonal to B
//
// r_dual(A) = r_cmpl(A) in spaces of even dimension
// l_dual(A) = l_cmpl(A) in spaces of even dimension
//
// r_dual(A) = l_dual(A) = cmpl(A) in spaces of odd dimension
//
// the right dual satisfies r_dual(A) = rev(A) * I_n
// the left dual satisfies  l_dual(A) = I_n * rev(A)
//
// -> derived from the defining equation of the left and right complements
////////////////////////////////////////////////////////////////////////////////


// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars, this can only be avoided, when the scalar
// type is uniquely defined for the corresponding algebra
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2d<T> r_dual(Scalar2d<T> s)
{
    return r_cmpl(s);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2d<T> r_dual(Vec2d<T> const& v)
{
    return r_cmpl(v);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2d<T> r_dual(PScalar2d<T> ps)
{
    return r_cmpl(ps);
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d_E<T> r_dual(MVec2d_E<T> const& M)
{
    return MVec2d_E<T>(r_cmpl(gr2(M)), r_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d<T> r_dual(MVec2d<T> const& M)
{
    return MVec2d<T>(r_cmpl(gr2(M)), r_cmpl(gr1(M)), r_cmpl(gr0(M)));
}


// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars, this can only be avoided, when the scalar
// type is uniquely defined for the corresponding algebra
template <typename T>
    requires(numeric_type<T>)
constexpr PScalar2d<T> l_dual(Scalar2d<T> s)
{
    return l_cmpl(s);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec2d<T> l_dual(Vec2d<T> const& v)
{
    return l_cmpl(v);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar2d<T> l_dual(PScalar2d<T> ps)
{
    return l_cmpl(ps);
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d_E<T> l_dual(MVec2d_E<T> const& M)
{
    return MVec2d_E<T>(l_cmpl(gr2(M)), l_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec2d<T> l_dual(MVec2d<T> const& M)
{
    return MVec2d<T>(l_cmpl(gr2(M)), l_cmpl(gr1(M)), l_cmpl(gr0(M)));
}

} // namespace hd::ga::ega