#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "detail/ga_foundation.hpp"    // ga library headers and infrastructure
#include "detail/type_t/ga_mvec3d.hpp" // 3D multivector types (includes component types)

namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega3d basic operations:
//
// - gr_inv()                     -> grade inversion
// - rev()                        -> reversion
// - conj()                       -> conjugation
// - cmpl()                       -> complement
// - dual()                       -> dual
//
/////////////////////////////////////////////////////////////////////////////////////////


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
    //     u ^ cmpl(u) = e1^e2^e3 => cmpl(u) = rev(u) * I_3d
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
// the outer product.
//
// The only difference between the dual and the complement is the fact that
// for the dual the argument is multiplied with the extended metric before
// taking the complement. However, in ega the metric is the identity matrix.
// Thus the dual and the complement are identical.
////////////////////////////////////////////////////////////////////////////////
//
// if M represents the subspace B as subspace of R^2 then
// left_dual(M) and right_dual(M) represent a subspace orthorgonal to B
//
// right_dual(A) = rcmpl(A) in spaces of even dimension
// left_dual(A)  = lcmpl(A) in spaces of even dimension
//
// right_dual(A) = left_dual(A) = cmpl(A) in spaces of odd dimension
//
// the right dual satisfies right_dual(A) = rev(A) * I_n
// the left dual satisfies   left_dual(A) = I_n * rev(A)
//
// -> derived from the defining equation of the left and right complements
////////////////////////////////////////////////////////////////////////////////


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

} // namespace hd::ga::ega