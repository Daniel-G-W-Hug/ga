#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "detail/ga_foundation.hpp"        // ga library headers and infrastructure
#include "detail/type_t/ga_dualnum4ds.hpp" // DualNum4ds
#include "detail/type_t/ga_mvec4ds.hpp"    // 4ds multivector types
#include "detail/type_t/ga_type4ds.hpp"    // Point4ds, Vector3d, Point3d, Line3d, Plane3d

#include "detail/ga_error_handling.hpp"


namespace hd::ga::sta {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga4ds basic operations:
//
// - gr_inv()                      -> grade inversion
// - rev()                         -> reversion
// - rrev()                        -> regressive reversion
// - conj()                        -> conjugation
//
// - l_cmpl(), r_cmpl()            -> left and right complement
// - nrm_sq(), nrm()               -> return norm
//
// - normalize()                   -> return normalized object (nrm scaled to 1.0)
//
// - r_dual()                      -> return right dual
// - l_dual()                      -> return right dual
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for r = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> gr_inv(Scalar4ds<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> gr_inv(Vec4ds<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> gr_inv(BiVec4ds<T> const& B)
{
    // grade 2: no sign change
    return B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> gr_inv(TriVec4ds<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> gr_inv(PScalar4ds<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> gr_inv(MVec4ds_E<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    return M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> gr_inv(MVec4ds_U<T> const& M)
{
    // grade 1 and 3: sign reversal
    return -M;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> gr_inv(MVec4ds<T> const& M)
{
    // grade 0, 2, and 4: no sign change
    // grade 1 and 3: sign reversal
    return MVec4ds<T>(gr_inv(gr0(M)), gr_inv(gr1(M)), gr_inv(gr2(M)), gr_inv(gr3(M)),
                      gr_inv(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// reversion operation: reverses the factors in a blade
// rev(A_r) = (-1)^(r*(r-1)/2) A_r
// pattern for r = 0, 1, 2, 3, ...: + + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> rev(Scalar4ds<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> rev(Vec4ds<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> rev(BiVec4ds<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> rev(TriVec4ds<T> const& t)
{
    // grade 3: sign reversal
    return -t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> rev(PScalar4ds<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> rev(MVec4ds_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec4ds_E<T>(rev(gr0(M)), rev(gr2(M)), rev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> rev(MVec4ds_U<T> const& M)
{
    // grade 1: no sign change
    // grade 3: sign change
    return MVec4ds_U<T>(rev(gr1(M)), rev(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> rev(MVec4ds<T> const& M)
{
    // grade 0: no sign change
    // grade 1: no sign change
    // grade 2: sign change
    // grade 3: sign change
    // grade 4: no sign change
    return MVec4ds<T>(rev(gr0(M)), rev(gr1(M)), rev(gr2(M)), rev(gr3(M)), rev(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// regressive reversion operation: reverse of complement + backtrafo
// rrev(A_r) = l_cmpl((-1)^(r*(r-1)/2) r_cmpl(A_r)) = (-1)^((n-r)*((n-r)-1)/2) A_r
// pattern for n=4, r = 0, 1, 2, 3, ...: + - - + + - - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> rrev(Scalar4ds<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> rrev(Vec4ds<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> rrev(BiVec4ds<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> rrev(TriVec4ds<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> rrev(PScalar4ds<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> rrev(MVec4ds_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign reversal
    return MVec4ds_E<T>(rrev(gr0(M)), rrev(gr2(M)), rrev(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> rrev(MVec4ds_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec4ds_U<T>(rrev(gr1(M)), rrev(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> rrev(MVec4ds<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    // grade 4: no sign change
    return MVec4ds<T>(rrev(gr0(M)), rrev(gr1(M)), rrev(gr2(M)), rrev(gr3(M)),
                      rrev(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// Clifford conjugation (identical sign pattern to regressive reverse for pga4ds):
// conj(A_r) = (-1)^(r*(r+1)/2) A_r
// pattern for r = 0, 1, 2, 3, ...: + - - + + - - + + ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> conj(Scalar4ds<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> conj(Vec4ds<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> conj(BiVec4ds<T> const& B)
{
    // grade 2: sign reversal
    return -B;
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> conj(TriVec4ds<T> const& t)
{
    // grade 3: no sign change
    return t;
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> conj(PScalar4ds<T> ps)
{
    // grade 4: no sign change
    return ps;
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> conj(MVec4ds_E<T> const& M)
{
    // grade 0, 4: no sign change
    // grade 2: sign change
    return MVec4ds_E<T>(conj(gr0(M)), conj(gr2(M)), conj(gr4(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> conj(MVec4ds_U<T> const& M)
{
    // grade 1: sign reversal
    // grade 3: no sign change
    return MVec4ds_U<T>(conj(gr1(M)), conj(gr3(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> conj(MVec4ds<T> const& M)
{
    // grade 0: no sign change
    // grade 1: sign reversal
    // grade 2: sign reversal
    // grade 3: no sign change
    // grade 4: no sign change
    return MVec4ds<T>(conj(gr0(M)), conj(gr1(M)), conj(gr2(M)), conj(gr3(M)),
                      conj(gr4(M)));
}


////////////////////////////////////////////////////////////////////////////////
// 4ds complement operations
// (the concept of complement is defined w.r.t. the outer product)
////////////////////////////////////////////////////////////////////////////////

// if M represents the subspace B of the blade u as subspace of R^2 then
// compl(M) represents the subspace orthorgonal to B
// the complement exchanges basis vectors which are in the k-blade u with
// the basis vectors which are NOT contained in the k-blade u
// and are needed to fill the space completely to the corresponding pseudoscalar
//
// left complement:  l_cmpl(u) ^ u  = I_4ds = g0^g1^g2^g3
// right complement: u ^ r_cmpl(u)  = I_4ds = g0^g1^g2^g3
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined l_compl(u), r_compl(u) => compl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// comploments have different signs
//
// complement operation with g0^g1^g2^g3 as the pseudoscalar

////////////////////////////////////////////////////////////////////////////////
// right complements
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> r_cmpl(Scalar4ds<T> s)
{
    // u ^ r_cmpl(u) = g0^g1^g2^g3
    // u = 1:
    //     1 ^ r_cmpl(u) = g0^g1^g2^g3 => r_cmpl(u) = s g0^g1^g2^g3
    return PScalar4ds<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> r_cmpl(Vec4ds<T> const& v)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = v.x e1 + v.y e2 + v.z e3 + v.w e4:
    //     u ^ r_cmpl(u) = e1^e2^e3^e4 =>
    //     u = e1 => r_cmpl(u) = e423
    //     u = e2 => r_cmpl(u) = e431
    //     u = e3 => r_cmpl(u) = e412
    //     u = e4 => r_cmpl(u) = e321
    return TriVec4ds<T>(v.x, v.y, v.z, v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> r_cmpl(BiVec4ds<T> const& B)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = B.vx e41 + B.vy e42 + B.vz e43 + B.mx e23 + B.my e31 + B.mz e12:
    //     u ^ r_cmpl(u) = e1^e2^e3^e4 =>
    //     u = e41 => r_cmpl(u) = -e23
    //     u = e42 => r_cmpl(u) = -e31
    //     u = e43 => r_cmpl(u) = -e12
    //     u = e23 => r_cmpl(u) = -e41
    //     u = e31 => r_cmpl(u) = -e42
    //     u = e12 => r_cmpl(u) = -e43
    return BiVec4ds<T>(-B.mx, -B.my, -B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> r_cmpl(TriVec4ds<T> const& t)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = t.x e423 + t.y e431 + t.z e412 + t.w e321:
    //     u ^ r_cmpl(u) = e1^e2^e3^e4 =>
    //     u = e423 => r_cmpl(u) = -e1
    //     u = e431 => r_cmpl(u) = -e2
    //     u = e412 => r_cmpl(u) = -e3
    //     u = e321 => r_cmpl(u) = -e4
    return Vec4ds<T>(-t.x, -t.y, -t.z, -t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> r_cmpl(PScalar4ds<T> ps)
{
    // u ^ r_compl(u) = e1^e2^e3^e4
    // u = e1^e2^e3^e4:
    //     e1^e2^e3^e4 ^ r_cmpl(u) = e1^e2^e3^e4 => r_cmpl(u) = ps 1
    return Scalar4ds<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> r_cmpl(MVec4ds_U<T> const& M)
{
    // use the component complements directly
    return MVec4ds_U<T>(r_cmpl(gr3(M)), r_cmpl(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> r_cmpl(MVec4ds_E<T> const& M)
{
    // use the component complements directly
    return MVec4ds_E<T>(r_cmpl(gr4(M)), r_cmpl(gr2(M)), r_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> r_cmpl(MVec4ds<T> const& M)
{
    // use the component complements directly
    return MVec4ds<T>(r_cmpl(gr4(M)), r_cmpl(gr3(M)), r_cmpl(gr2(M)), r_cmpl(gr1(M)),
                      r_cmpl(gr0(M)));
}

////////////////////////////////////////////////////////////////////////////////
// left complements
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> l_cmpl(Scalar4ds<T> s)
{
    // l_cmpl(u) ^ u = g0^g1^g2^g3
    // u = 1:
    //     l_cmpl(u) ^ 1 = g0^g1^g2^g3 => l_cmpl(u) = s g0^g1^g2^g3
    return PScalar4ds<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> l_cmpl(Vec4ds<T> const& v)
{
    // l_cmpl(u) ^ u = g0^g1^g2^g3
    // u = v.x g1 + v.y g2 + v.z g3 + v.w g0:
    //     l_cmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e1 => l_cmpl(u) = -e423
    //     u = e2 => l_cmpl(u) = -e431
    //     u = e3 => l_cmpl(u) = -e412
    //     u = e4 => l_cmpl(u) = -e321
    return TriVec4ds<T>(-v.x, -v.y, -v.z, -v.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> l_cmpl(BiVec4ds<T> const& B)
{
    // l_cmpl(u) ^ u = e1^e2^e3^e4
    // u = B.vx e41 + B.vy e42 + B.vz e43 + B.mx e23 + B.my e31 + B.mz e12:
    //     l_cmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e41 => l_cmpl(u) = -e23
    //     u = e42 => l_cmpl(u) = -e31
    //     u = e43 => l_cmpl(u) = -e12
    //     u = e23 => l_cmpl(u) = -e41
    //     u = e31 => l_cmpl(u) = -e42
    //     u = e12 => l_cmpl(u) = -e43
    return BiVec4ds<T>(-B.mx, -B.my, -B.mz, -B.vx, -B.vy, -B.vz);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> l_cmpl(TriVec4ds<T> const& t)
{
    // l_cmpl(u) ^ u = e1^e2^e3^e4
    // u = t.x e423 + t.y e431 + t.z e412 + t.w e321:
    //     l_cmpl(u) ^ u = e1^e2^e3^e4 =>
    //     u = e423 => l_cmpl(u) = e1
    //     u = e431 => l_cmpl(u) = e2
    //     u = e412 => l_cmpl(u) = e3
    //     u = e321 => l_cmpl(u) = e4
    return Vec4ds<T>(t.x, t.y, t.z, t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> l_cmpl(PScalar4ds<T> ps)
{
    // l_cmpl(u) ^ u = e1^e2^e3^e4
    // u = e1^e2^e3^e4:
    //     l_cmpl(u) ^ e1^e2^e3^e4 = e1^e2^e3^e4 => l_cmpl(u) = ps 1
    return Scalar4ds<T>(T(ps));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> l_cmpl(MVec4ds_U<T> const& M)
{
    // use the component complements directly
    return MVec4ds_U<T>(l_cmpl(gr3(M)), l_cmpl(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> l_cmpl(MVec4ds_E<T> const& M)
{
    // use the component complements directly
    return MVec4ds_E<T>(l_cmpl(gr4(M)), l_cmpl(gr2(M)), l_cmpl(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> l_cmpl(MVec4ds<T> const& M)
{
    // use the component complements directly
    return MVec4ds<T>(l_cmpl(gr4(M)), l_cmpl(gr3(M)), l_cmpl(gr2(M)), l_cmpl(gr1(M)),
                      l_cmpl(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// nrm(u): return the bulk norm of u (a scalar value)
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(Scalar4ds<T> s)
{
    return T(s) * T(s);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(Scalar4ds<T> s)
{
    return Scalar4ds<T>(std::sqrt(nrm_sq(s)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(Vec4ds<T> const& v)
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(Vec4ds<T> const& v)
{
    return Scalar4ds<T>(std::sqrt(nrm_sq(v)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(BiVec4ds<T> const& B)
{
    return B.mx * B.mx + B.my * B.my + B.mz * B.mz;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(BiVec4ds<T> const& B)
{
    return Scalar4ds<T>(std::sqrt(nrm_sq(B)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(TriVec4ds<T> const& t)
{
    return t.w * t.w;
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(TriVec4ds<T> const& t)
{
    return Scalar4ds<T>(std::sqrt(nrm_sq(t)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq([[maybe_unused]] PScalar4ds<T>)
{
    return T(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm([[maybe_unused]] PScalar4ds<T>)
{
    return Scalar4ds<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(MVec4ds_E<T> const& M)
{
    return nrm_sq(gr0(M)) + nrm_sq(gr2(M)) + nrm_sq(gr4(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(MVec4ds_E<T> const& M)
{
    return Scalar4ds<T>(std::sqrt(nrm_sq(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(MVec4ds_U<T> const& M)
{
    return nrm_sq(gr1(M)) + nrm_sq(gr3(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(MVec4ds_U<T> const& M)
{
    return Scalar4ds<T>(std::sqrt(nrm_sq(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr T nrm_sq(MVec4ds<T> const& M)
{
    return nrm_sq(gr0(M)) + nrm_sq(gr1(M)) + nrm_sq(gr2(M)) + nrm_sq(gr3(M)) +
           nrm_sq(gr4(M));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> nrm(MVec4ds<T> const& M)
{
    return Scalar4ds<T>(std::sqrt(nrm_sq(M)));
}


////////////////////////////////////////////////////////////////////////////////
// normalization operations:
// return an object normalized to nrm(u) == 1.0
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
inline Scalar4ds<T> normalize(Scalar4ds<T> s)
{
    T m = to_val(nrm(s));
    hd::ga::detail::check_normalization<T>(m, "scalar (4ds)");
    T scale = T(1.0) / m;
    return Scalar4ds<T>(scale * T(s));
}

template <typename T>
    requires(numeric_type<T>)
inline Vec4ds<T> normalize(Vec4ds<T> const& v)
{
    T n = to_val(nrm(v));
    hd::ga::detail::check_normalization<T>(std::abs(n), "vector (4ds)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * v;
}

template <typename T>
    requires(numeric_type<T>)
inline BiVec4ds<T> normalize(BiVec4ds<T> const& B)
{
    T n = to_val(nrm(B));
    hd::ga::detail::check_normalization<T>(n, "bivector (4ds)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * B;
}

template <typename T>
    requires(numeric_type<T>)
inline TriVec4ds<T> normalize(TriVec4ds<T> const& t)
{
    T n = to_val(nrm(t));
    hd::ga::detail::check_normalization<T>(n, "trivector (4ds)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * t;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> normalize(MVec4ds_E<T> const& M)
{
    T n = to_val(nrm(M));
    hd::ga::detail::check_normalization<T>(n, "even-grade multivector (4ds)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * M;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_U<T> normalize(MVec4ds_U<T> const& M)
{
    T n = to_val(nrm(M));
    hd::ga::detail::check_normalization<T>(n, "odd-grade multivector (4ds)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * M;
}

template <typename T>
    requires(numeric_type<T>)
inline MVec4ds<T> normalize(MVec4ds<T> const& M)
{
    T n = to_val(nrm(M));
    hd::ga::detail::check_normalization<T>(n, "multivector (4ds)");
    T scale = T(1.0) / n; // for multiplication with inverse of norm
    return scale * M;
}


////////////////////////////////////////////////////////////////////////////////
// r_dual(A) = r_cmpl(A) = r_cmpl( metric * A )
//
// -> right complement operation applied to the bulk
//
////////////////////////////////////////////////////////////////////////////////
// duality (as defined in Lengyel, "PGA illuminated")
// is defined w.r.t. the outer product, because the complement is defined w.r.t.
// the outer product.
//
// The only difference between the dual and the complement is the fact that
// for the dual the argument is multiplied with the extended metric before
// taking the complement.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> r_dual(Scalar4ds<T> s)
{
    return PScalar4ds<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> r_dual(Vec4ds<T> const& v)
{
    return TriVec4ds<T>(v.x, v.y, v.z, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> r_dual(BiVec4ds<T> const& B)
{
    return BiVec4ds<T>(-B.mx, -B.my, -B.mz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> r_dual(TriVec4ds<T> const& t)
{
    return Vec4ds<T>(T(0.0), T(0.0), T(0.0), -t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> r_dual([[maybe_unused]] PScalar4ds<T>)
{
    return Scalar4ds<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> r_dual(MVec4ds_E<T> const& M)
{
    return MVec4ds_E<T>(r_dual(gr4(M)), r_dual(gr2(M)), r_dual(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> r_dual(MVec4ds_U<T> const& M)
{
    return MVec4ds_U<T>(r_dual(gr3(M)), r_dual(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> r_dual(MVec4ds<T> const& M)
{
    return MVec4ds<T>(r_dual(gr4(M)), r_dual(gr3(M)), r_dual(gr2(M)), r_dual(gr1(M)),
                      r_dual(gr0(M)));
}


////////////////////////////////////////////////////////////////////////////////
// l_dual(A) = l_cmpl(A) = l_cmpl( metric * A )
//
// -> left complement operation
//
////////////////////////////////////////////////////////////////////////////////
// duality (as defined in Lengyel, "PGA illuminated")
// is defined w.r.t. the outer product, because the complement is defined w.r.t.
// the outer product.
//
// The only difference between the dual and the complement is the fact that
// for the dual the argument is multiplied with the extended metric before
// taking the complement.
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> l_dual(Scalar4ds<T> s)
{
    return PScalar4ds<T>(T(s));
}

template <typename T>
    requires(numeric_type<T>)
constexpr TriVec4ds<T> l_dual(Vec4ds<T> const& v)
{
    return TriVec4ds<T>(-v.x, -v.y, -v.z, T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr BiVec4ds<T> l_dual(BiVec4ds<T> const& B)
{
    return BiVec4ds<T>(-B.mx, -B.my, -B.mz, T(0.0), T(0.0), T(0.0));
}

template <typename T>
    requires(numeric_type<T>)
constexpr Vec4ds<T> l_dual(TriVec4ds<T> const& t)
{
    return Vec4ds<T>(T(0.0), T(0.0), T(0.0), t.w);
}

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> l_dual([[maybe_unused]] PScalar4ds<T>)
{
    return Scalar4ds<T>(0.0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_E<T> l_dual(MVec4ds_E<T> const& M)
{
    return MVec4ds_E<T>(l_dual(gr4(M)), l_dual(gr2(M)), l_dual(gr0(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds_U<T> l_dual(MVec4ds_U<T> const& M)
{
    return MVec4ds_U<T>(l_dual(gr3(M)), l_dual(gr1(M)));
}

template <typename T>
    requires(numeric_type<T>)
constexpr MVec4ds<T> l_dual(MVec4ds<T> const& M)
{
    return MVec4ds<T>(l_dual(gr4(M)), l_dual(gr3(M)), l_dual(gr2(M)), l_dual(gr1(M)),
                      l_dual(gr0(M)));
}


} // namespace hd::ga::sta