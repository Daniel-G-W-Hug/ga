#pragma once

// author: Daniel Hug, 2024

#include "type_t/ga_type_2d.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// MVec2d_E<T> M = c0 + c1 I (with I being the bivector of the plane e1^e2)
//
// This is the definition of a multivector in the even subalgebra of G<2,0,0>,
// i.e. it models only multivectors with even grades 0 and 2 in the plane e1^e2
// (complex numbers).
// This subalgebra is closed under addition and multiplication.
////////////////////////////////////////////////////////////////////////////////

// This is defined in order to limit memory requirements and computational effort
// for these sepecific multivectors vs. usage of a fully populated multivectors.
// At the same time this enables easy integration with fully populated
// multivectors, if required by the application.

////////////////////////////////////////////////////////////////////////////////
// use MVec2_t including its ctors and add specific ctors for MVec2_t<T, Tag>
// by using partial template specialization for the Tag=mvec2d_tag
////////////////////////////////////////////////////////////////////////////////

template <typename T> struct MVec2_t<T, mvec2d_e_tag> : public MVec2_t<T, default_tag> {

    using MVec2_t<T, default_tag>::MVec2_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    MVec2_t(Scalar2d<T> s) : MVec2_t(s, T(0.0)) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    MVec2_t(PScalar2d<T> ps) : MVec2_t(T(0.0), ps) {}

    // assign a geometric product resulting from a product of two vectors
    // via dot(v1,v2) and wdg(v1,v2) directly (other grades = 0)
    // (less expensive compared to full geometric product)
    MVec2_t(Scalar2d<T> s, PScalar2d<T> ps) : MVec2_t(T(s), T(ps)) {}
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization MVec2_t<T, mvec2d_e_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector from the even subalgebra
//
// grade 0: gr0() - scalar
// grade 2: gr2() - bivector (= pseudoscalar in 2d)

template <typename T> inline constexpr Scalar2d<T> gr0(MVec2d_E<T> const& v)
{
    return Scalar2d<T>(v.c0);
}

template <typename T> inline constexpr PScalar2d<T> gr2(MVec2d_E<T> const& v)
{
    return PScalar2d<T>(v.c1);
}

} // namespace hd::ga