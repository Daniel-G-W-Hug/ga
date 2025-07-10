#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "detail/ga_foundation.hpp"    // ga library headers and infrastructure
#include "detail/type_t/ga_mvec2d.hpp" // 2D multivector types (includes component types)

namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega2d basic operations:
//
// - gr_inv()                     -> grade inversion
// - rev()                        -> reversion
// - conj()                       -> conjugation
// - lcmpl(), rcmpl()             -> left and right complements
// - left_dual(), right_dual()    -> left and right duals
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// grade inversion operation: reverse the sign of odd blades
// gr_inv(A_r) = (-1)^r A_r
// pattern for k = 0, 1, 2, 3, ...: + - + - + - ...
////////////////////////////////////////////////////////////////////////////////

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> gr_inv(Scalar2d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> gr_inv(Vec2d<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> gr_inv(PScalar2d<T> ps)
{
    // grade 2: no sign change
    return ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> gr_inv(MVec2d_E<T> const& M)
{
    // grade 0 and 2: no sign change
    return M;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> gr_inv(MVec2d<T> const& M)
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
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> rev(Scalar2d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> rev(Vec2d<T> const& v)
{
    // grade 1: no sign change
    return v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> rev(PScalar2d<T> ps)
{
    // grade 2: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> rev(MVec2d_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign reversal
    return MVec2d_E<T>(rev(gr0(M)), rev(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> rev(MVec2d<T> const& M)
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
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> conj(Scalar2d<T> s)
{
    // grade 0: no sign change
    return s;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> conj(Vec2d<T> const& v)
{
    // grade 1: sign reversal
    return -v;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> conj(PScalar2d<T> ps)
{
    // grade 2: sign reversal
    return -ps;
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> conj(MVec2d_E<T> const& M)
{
    // grade 0: no sign change
    // grade 2: sign reversal
    return MVec2d_E<T>(conj(gr0(M)), conj(gr2(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> conj(MVec2d<T> const& M)
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
// left complement:  lcmpl(u) ^ u  = I_2d = e1^e2  =>  lcmpl(u) = I_2d * rev(u)
// right complement: u ^ rcmpl(u)  = I_2d = e1^e2  =>  rcmpl(u) = rev(u) * I_2d
//
// (get formula on the rhs by multiplying with inv(u) from the right or left respectively)
//
// in spaces of odd dimension right and left complements are identical and thus there
// is only one complement operation defined lcmpl(u) = rcmpl(u) = cmpl(u)
//
// in spaces of even dimension and when the grade of the k-vector is odd left and right
// complements have different signs

template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> lcmpl(Scalar2d<T> s)
{
    // lcmpl(u) ^ u = e1^e2
    // u = s 1:
    //     lcmpl(u) ^ u = e1^e2 => lcmpl(u) = I_2d * rev(s) = s e1^e2
    return PScalar2d<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> lcmpl(Vec2d<T> const& v)
{
    // lcmpl(u) ^ u = e1^e2
    // u = v.x e1 + v.y e2:
    //     lcmpl(u) ^ u = e1^e2 => cmpl(u) = I_2d * rev(v)
    //     lcmpl(u) = e12 * (v.x * e1 + v.y * e2) = -v.x * e2 + v.y * e1
    //                                            =  v.y * e1 - v.x * e2
    //     => lcmpl(u).x =  v.y e1
    //     => lcmpl(u).y = -v.x e2
    return Vec2d<T>(v.y, -v.x);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> lcmpl(PScalar2d<T> ps)
{
    // lcmpl(u) ^ u = e1^e2
    // u = ps e1^e2:
    //     lcmpl(u) ^ u = e1^e2 => lcmpl(u) = I_2d * rev(ps) = ps 1
    return Scalar2d<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> lcmpl(MVec2d_E<T> const& M)
{
    // use the component complements directly
    return MVec2d_E<T>(lcmpl(gr2(M)), lcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> lcmpl(MVec2d<T> const& M)
{
    // use the component complements directly
    return MVec2d<T>(lcmpl(gr2(M)), lcmpl(gr1(M)), lcmpl(gr0(M)));
}


template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> rcmpl(Scalar2d<T> s)
{
    // u ^ rcmpl(u) = e1^e2
    // u = s 1:
    //     u ^ rcmpl(u) = e1^e2 => rcmpl(u) = rev(s) * I_2d = s e1^e2
    return PScalar2d<T>(T(s));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> rcmpl(Vec2d<T> const& v)
{
    // u ^ rcmpl(u) = e1^e2
    // u = v.x e1 + v.y e2:
    //     u ^ rcmpl(u) = e1^e2 => rcmpl(u) = rev(v) * I_2d
    //     rcmpl(u) = (v.x * e1 + v.y * e2) * e12 =   v.x * e2 - v.y * e1
    //                                            = - v.y * e1 + v.x * e2
    //     => rcmpl(u).x = -v.y e1
    //     => rcmpl(u).y =  v.x e2

    return Vec2d<T>(-v.y, v.x);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> rcmpl(PScalar2d<T> ps)
{
    // u ^ rcmpl(u) = e1^e2
    // u = ps e1^e2:
    //     u ^ rcmpl(u) = e1^e2 => rcmpl(u) = rev(ps) * I_2d = ps 1
    return Scalar2d<T>(T(ps));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> rcmpl(MVec2d_E<T> const& M)
{
    // use the component complements directly
    return MVec2d_E<T>(rcmpl(gr2(M)), rcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> rcmpl(MVec2d<T> const& M)
{
    // use the component complements directly
    return MVec2d<T>(rcmpl(gr2(M)), rcmpl(gr1(M)), rcmpl(gr0(M)));
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


// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars, this can only be avoided, when the scalar
// type is uniquely defined for the corresponding algebra
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> right_dual(Scalar2d<T> s)
{
    return rcmpl(s);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> right_dual(Vec2d<T> const& v)
{
    return rcmpl(v);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> right_dual(PScalar2d<T> ps)
{
    return rcmpl(ps);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> right_dual(MVec2d_E<T> const& M)
{
    return MVec2d_E<T>(rcmpl(gr2(M)), rcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> right_dual(MVec2d<T> const& M)
{
    return MVec2d<T>(rcmpl(gr2(M)), rcmpl(gr1(M)), rcmpl(gr0(M)));
}


// this one is problematic for overloading, because 2d and 3d case
// transform scalars to different pseudoscalars, this can only be avoided, when the scalar
// type is uniquely defined for the corresponding algebra
template <typename T>
    requires(std::floating_point<T>)
inline constexpr PScalar2d<T> left_dual(Scalar2d<T> s)
{
    return lcmpl(s);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Vec2d<T> left_dual(Vec2d<T> const& v)
{
    return lcmpl(v);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr Scalar2d<T> left_dual(PScalar2d<T> ps)
{
    return lcmpl(ps);
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d_E<T> left_dual(MVec2d_E<T> const& M)
{
    return MVec2d_E<T>(lcmpl(gr2(M)), lcmpl(gr0(M)));
}

template <typename T>
    requires(std::floating_point<T>)
inline constexpr MVec2d<T> left_dual(MVec2d<T> const& M)
{
    return MVec2d<T>(lcmpl(gr2(M)), lcmpl(gr1(M)), lcmpl(gr0(M)));
}

} // namespace hd::ga::ega