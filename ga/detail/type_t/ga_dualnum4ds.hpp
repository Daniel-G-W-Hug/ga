#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_type4ds.hpp"


namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// DualNum4ds<T> M = c0 1 + c1 I_4ds (dual number with s 1 + ps I_4ds,
//                                    where I_4ds = e1234, e1234^2 = 0)
/////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
// use MVec2_t including its ctors and add specific ctors for MVec2_t<T, Tag>
// by using partial template specialization for the Tag=dual_number4ds_tag
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
struct MVec2_t<T, dual_number4ds_tag> : public MVec2_t<T, default_tag> {

    using MVec2_t<T, default_tag>::MVec2_t; // inherit base class ctors

    // assign a scalar part exclusively (other grades = 0)
    constexpr MVec2_t(Scalar4ds<T> s) : MVec2_t(T(s), T(0.0)) {}

    // assign a pseudoscalar part exclusively (other grades = 0)
    constexpr MVec2_t(PScalar4ds<T> ps) : MVec2_t(T(0.0), T(ps)) {}

    // assign from a scalar and a pseudoscalar part
    constexpr MVec2_t(Scalar4ds<T> s, PScalar4ds<T> ps) : MVec2_t(T(s), T(ps)) {}

    // Override compound assignment operators to return correct derived type
    // This ensures GCC+doctest can properly deduce the tag type without needing cross-tag
    // comparisons
    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator+=(MVec2_t<U, dual_number4ds_tag> const& v) noexcept
    {
        this->c0 += v.c0;
        this->c1 += v.c1;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator-=(MVec2_t<U, dual_number4ds_tag> const& v) noexcept
    {
        this->c0 -= v.c0;
        this->c1 -= v.c1;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator*=(U s) noexcept
    {
        this->c0 *= s;
        this->c1 *= s;
        return *this;
    }

    template <typename U>
        requires(numeric_type<U>)
    MVec2_t& operator/=(U s) noexcept(!detail::extended_testing_enabled())
    {
        detail::check_division_by_zero<T, U>(s, "multivector division 2 comp.");
        this->c0 /= s;
        this->c1 /= s;
        return *this;
    }
};

////////////////////////////////////////////////////////////////////////////////
// define grade operations for partial specialization
// DualNum4ds<T> = MVec2_t<T, dual_number4ds_tag>
////////////////////////////////////////////////////////////////////////////////

// returning various grades of a multivector from the even subalgebra
//
// grade 0: gr0() - scalar
// grade 4: gr4() - pseudoscalar in pga4ds

template <typename T>
    requires(numeric_type<T>)
constexpr Scalar4ds<T> gr0(DualNum4ds<T> const& M)
{
    return Scalar4ds<T>(M.c0);
}

template <typename T>
    requires(numeric_type<T>)
constexpr PScalar4ds<T> gr4(DualNum4ds<T> const& M)
{
    return PScalar4ds<T>(M.c1);
}

////////////////////////////////////////////////////////////////////////////////
// addition operations to combine dual numbers
////////////////////////////////////////////////////////////////////////////////

// scalar + pseudoscalar => dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator+(Scalar4ds<T> s, PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(s, ps);
}

// pseudoscalar + scalar => dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator+(PScalar4ds<T> ps, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(s, ps);
}

// scalar + dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator+(Scalar4ds<T> s,
                                                         DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(T(s) + M.c0, M.c1);
}

// dual number + scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator+(DualNum4ds<T> const& M,
                                                         Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(U(s) + M.c0, M.c1);
}

// pseudosscalar + dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator+(PScalar4ds<T> ps,
                                                         DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(M.c0, T(ps) + M.c1);
}

// dual number + pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator+(DualNum4ds<T> const& M,
                                                         PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(M.c0, U(ps) + M.c1);
}

////////////////////////////////////////////////////////////////////////////////
// subtraction operations to combine scalars and dual numbers
////////////////////////////////////////////////////////////////////////////////

// scalar - pseudoscalar => dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator-(Scalar4ds<T> s, PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(s, -ps);
}

// pseudoscalar - scalar => dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator-(PScalar4ds<T> ps, Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(-s, ps);
}

// scalar - dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator-(Scalar4ds<T> s,
                                                         DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(T(s) - M.c0, -M.c1);
}

// dual number - scalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator-(DualNum4ds<T> const& M,
                                                         Scalar4ds<U> s)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(M.c0 - U(s), M.c1);
}

// pseudosscalar - dual number
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator-(PScalar4ds<T> ps,
                                                         DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(-M.c0, T(ps) - M.c1);
}

// dual number - pseudoscalar
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr DualNum4ds<std::common_type_t<T, U>> operator-(DualNum4ds<T> const& M,
                                                         PScalar4ds<U> ps)
{
    using ctype = std::common_type_t<T, U>;
    return DualNum4ds<ctype>(M.c0, M.c1 - U(ps));
}


////////////////////////////////////////////////////////////////////////////////
// arithmetic between DualNum4ds and the richer multivector types
//
// In 4ds the dual number occupies grades 0+4, both of which lie within the
// even subalgebra (grades 0+2+4). Combinations with bivec or mvec_e
// therefore stay in the even subalgebra and produce MVec4ds_E. All other
// combinations break out of the even subalgebra and produce the full MVec4ds.
//
// Layouts (for reference):
//   MVec4ds_E (8 comps): c0=scalar, c1..c6=bivec(vx,vy,vz,mx,my,mz), c7=pscalar
//   MVec4ds_U (8 comps): c0..c3=vec(x,y,z,w), c4..c7=trivec(x,y,z,w)
//   MVec4ds  (16 comps): c0=scalar, c1..c4=vec, c5..c10=bivec, c11..c14=trivec,
//   c15=pscalar
////////////////////////////////////////////////////////////////////////////////

// dual number + bivector  =>  even multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator+(DualNum4ds<T> const& M,
                                                        BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(M.c0), B, PScalar4ds<ctype>(M.c1));
}

// bivector + dual number  =>  even multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator+(BiVec4ds<T> const& B,
                                                        DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(M.c0), B, PScalar4ds<ctype>(M.c1));
}

// dual number + even multivector  =>  even multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator+(DualNum4ds<T> const& M,
                                                        MVec4ds_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds_E<ctype>(M.c0 + E.c0, E.c1, E.c2, E.c3, E.c4, E.c5, E.c6, E.c7 + M.c1);
}

// even multivector + dual number  =>  even multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator+(MVec4ds_E<T> const& E,
                                                        DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds_E<ctype>(E.c0 + M.c0, E.c1, E.c2, E.c3, E.c4, E.c5, E.c6, E.c7 + M.c1);
}

// dual number + vector  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(DualNum4ds<T> const& M,
                                                      Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(M.c0), v, BiVec4ds<ctype>{},
                          TriVec4ds<ctype>{}, PScalar4ds<ctype>(M.c1));
}

// vector + dual number  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(Vec4ds<T> const& v,
                                                      DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(M.c0), v, BiVec4ds<ctype>{},
                          TriVec4ds<ctype>{}, PScalar4ds<ctype>(M.c1));
}

// dual number + trivector  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(DualNum4ds<T> const& M,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(M.c0), Vec4ds<ctype>{}, BiVec4ds<ctype>{}, t,
                          PScalar4ds<ctype>(M.c1));
}

// trivector + dual number  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(TriVec4ds<T> const& t,
                                                      DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(M.c0), Vec4ds<ctype>{}, BiVec4ds<ctype>{}, t,
                          PScalar4ds<ctype>(M.c1));
}

// dual number + odd multivector  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(DualNum4ds<T> const& M,
                                                      MVec4ds_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(M.c0), Vec4ds<ctype>(O.c0, O.c1, O.c2, O.c3),
                          BiVec4ds<ctype>{}, TriVec4ds<ctype>(O.c4, O.c5, O.c6, O.c7),
                          PScalar4ds<ctype>(M.c1));
}

// odd multivector + dual number  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(MVec4ds_U<T> const& O,
                                                      DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(M.c0), Vec4ds<ctype>(O.c0, O.c1, O.c2, O.c3),
                          BiVec4ds<ctype>{}, TriVec4ds<ctype>(O.c4, O.c5, O.c6, O.c7),
                          PScalar4ds<ctype>(M.c1));
}

// dual number + multivector  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(DualNum4ds<T> const& M,
                                                      MVec4ds<U> const& A)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(A.c0 + M.c0, A.c1, A.c2, A.c3, A.c4, A.c5, A.c6, A.c7, A.c8,
                          A.c9, A.c10, A.c11, A.c12, A.c13, A.c14, A.c15 + M.c1);
}

// multivector + dual number  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator+(MVec4ds<T> const& A,
                                                      DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(A.c0 + M.c0, A.c1, A.c2, A.c3, A.c4, A.c5, A.c6, A.c7, A.c8,
                          A.c9, A.c10, A.c11, A.c12, A.c13, A.c14, A.c15 + M.c1);
}


// ===== subtraction =====

// dual number - bivector  =>  even multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator-(DualNum4ds<T> const& M,
                                                        BiVec4ds<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(M.c0), -B, PScalar4ds<ctype>(M.c1));
}

// bivector - dual number  =>  even multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator-(BiVec4ds<T> const& B,
                                                        DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds_E<ctype>(Scalar4ds<ctype>(-M.c0), B, PScalar4ds<ctype>(-M.c1));
}

// dual number - even multivector  =>  even multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator-(DualNum4ds<T> const& M,
                                                        MVec4ds_E<U> const& E)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds_E<ctype>(M.c0 - E.c0, -E.c1, -E.c2, -E.c3, -E.c4, -E.c5, -E.c6,
                            M.c1 - E.c7);
}

// even multivector - dual number  =>  even multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds_E<std::common_type_t<T, U>> operator-(MVec4ds_E<T> const& E,
                                                        DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds_E<ctype>(E.c0 - M.c0, E.c1, E.c2, E.c3, E.c4, E.c5, E.c6, E.c7 - M.c1);
}

// dual number - vector  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(DualNum4ds<T> const& M,
                                                      Vec4ds<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(M.c0), -v, BiVec4ds<ctype>{},
                          TriVec4ds<ctype>{}, PScalar4ds<ctype>(M.c1));
}

// vector - dual number  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(Vec4ds<T> const& v,
                                                      DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(-M.c0), v, BiVec4ds<ctype>{},
                          TriVec4ds<ctype>{}, PScalar4ds<ctype>(-M.c1));
}

// dual number - trivector  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(DualNum4ds<T> const& M,
                                                      TriVec4ds<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(M.c0), Vec4ds<ctype>{}, BiVec4ds<ctype>{}, -t,
                          PScalar4ds<ctype>(M.c1));
}

// trivector - dual number  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(TriVec4ds<T> const& t,
                                                      DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(-M.c0), Vec4ds<ctype>{}, BiVec4ds<ctype>{}, t,
                          PScalar4ds<ctype>(-M.c1));
}

// dual number - odd multivector  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(DualNum4ds<T> const& M,
                                                      MVec4ds_U<U> const& O)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(M.c0),
                          Vec4ds<ctype>(-O.c0, -O.c1, -O.c2, -O.c3), BiVec4ds<ctype>{},
                          TriVec4ds<ctype>(-O.c4, -O.c5, -O.c6, -O.c7),
                          PScalar4ds<ctype>(M.c1));
}

// odd multivector - dual number  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(MVec4ds_U<T> const& O,
                                                      DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(Scalar4ds<ctype>(-M.c0), Vec4ds<ctype>(O.c0, O.c1, O.c2, O.c3),
                          BiVec4ds<ctype>{}, TriVec4ds<ctype>(O.c4, O.c5, O.c6, O.c7),
                          PScalar4ds<ctype>(-M.c1));
}

// dual number - multivector  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(DualNum4ds<T> const& M,
                                                      MVec4ds<U> const& A)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(M.c0 - A.c0, -A.c1, -A.c2, -A.c3, -A.c4, -A.c5, -A.c6, -A.c7,
                          -A.c8, -A.c9, -A.c10, -A.c11, -A.c12, -A.c13, -A.c14,
                          M.c1 - A.c15);
}

// multivector - dual number  =>  full multivector
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec4ds<std::common_type_t<T, U>> operator-(MVec4ds<T> const& A,
                                                      DualNum4ds<U> const& M)
{
    using ctype = std::common_type_t<T, U>;
    return MVec4ds<ctype>(A.c0 - M.c0, A.c1, A.c2, A.c3, A.c4, A.c5, A.c6, A.c7, A.c8,
                          A.c9, A.c10, A.c11, A.c12, A.c13, A.c14, A.c15 - M.c1);
}

} // namespace hd::ga