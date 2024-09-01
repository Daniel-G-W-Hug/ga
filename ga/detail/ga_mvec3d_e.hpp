#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_mvec4_t.hpp"

#include "type_t/ga_type_0d.hpp"
#include "type_t/ga_type_3d.hpp"


namespace hd::ga {

template <typename T> using MVec3d_E = MVec4_t<T, mvec3d_e_tag>;

////////////////////////////////////////////////////////////////////////////////
// MVec3d_E<T> M = c0 + (c1 e2^e3 + c2 e3^e1 + c3 e1^e2)
//
// with the term in brackets being the bivector c1 e2^e3 + c2 e3^e1 + c3 e1^e2
// that models a plane in 3D defined by the bivector coordinates (c1, c2, c3)
//
// This is the definition of a multivector in the even subalgebra of G<3,0,0>,
// i.e. it models only multivectors with even grades 0 and 2
// which are quaternions.
// This subalgebra is closed under addition and multiplication.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of a fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.

template <typename T> struct MVec4_t<T, mvec3d_e_tag> : public MVec4_t<T, default_tag> {

    using MVec4_t<T, default_tag>::MVec4_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec4_t(Scalar<T> s) : MVec4_t(T(s), T(0.0), T(0.0), T(0.0)) {}

    // assign a bivector part exclusively (other grades = 0)
    MVec4_t(BiVec3d<T> b) : MVec4_t(T(0.0), b.x, b.y, b.z) {}

    // assign scalar and bivector parts
    MVec4_t(Scalar<T> s, BiVec3d<T> b) : MVec4_t(T(s), b.x, b.y, b.z) {}
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec4_t<T, mvec3d_e_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of the even multivector
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector

template <typename T> inline constexpr Scalar<T> gr0(MVec3d_E<T> const& v)
{
    return Scalar<T>(v.c0);
}

template <typename T> inline constexpr BiVec3d<T> gr2(MVec3d_E<T> const& v)
{
    return BiVec3d<T>(v.c1, v.c2, v.c3);
}

} // namespace hd::ga