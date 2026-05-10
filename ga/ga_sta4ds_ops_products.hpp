#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_sta4ds_ops_basics.hpp"


namespace hd::ga::sta {

/////////////////////////////////////////////////////////////////////////////////////////
// provides pga4ds product operations:
//
// - dot()                   -> dot product
// - rdot()                  -> regressive dot product
// - wdg()                   -> wedge product (join as convenience interface)
// - rwdg()                  -> regressive wedge product (meet as convenience interface)
// - operator<<()            -> left contraction
// - operator>>()            -> right contraction
// - cmt()                   -> commutator product (= asymmetric part of gpr)
// - rcmt()                  -> regressive commutator product (=asymmetric part of rgpr)
// - operator*()             -> geometric product (= gpr)
// - rgpr()                  -> regressive geometric product
// - inv()                   -> inversion operation (w.r.t. geometric product)
//
// Expansions are the regressive versions of the corresponding contractions.
//
// - l_expand4ds()           -> left expansion
// - r_expand4ds()           -> right expansion
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// dot products for 4ds (=inner product defined for equal grades exclusively)
//
// dot(v1,v2) = v1^T * g_{ij} * v2 is the scalar product with g_{ij} as the metric
//
// here we assume g0^2 = -1, g1^2 = +1, g2^2 = +1, g3^2 = +1
//
// Via the exomorphism G(a ^ b) = G(a) ^ G(b) the extended metric G is defined as a 2^n by
// 2^n matrix via the outer product, which in turn allows for definition of an extended
// dot product for all grades of a full multi-vector from scalar to pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// regressive dot products for 4ds (= defined for equal grades exclusively)
//
// rdot(v1,v2) = l_cmpl( dot(r_cmpl(v1),r_cmpl(v2)) )
//
// returns a pseudoscalar
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// wedge products (= outer product) and join operations
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// regressive wedge product (= outer product for complements) and meet operations
// as defined by E. Lengyel in "Projective geometric algebra illuminated"
// independent of the geometric product, just depending on the outer product (wdg)
// as well as the complement and thus the pseudoscalar of the space
// (this definition does NOT connect directly to the geometric product,
// but to the outer product exclusively)
//
//        rwdg(ul, ur) = l_cmpl(wdg(r_cmpl(ul),r_cmpl(ur))) = l_cmpl(r_cmpl(ul) ^
//        r_cmpl(ur))
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// left contractions A << B: "A contracted onto B"
//
// The resulting object is a lies in B and is perpendicular to A
//
// Implements the left contraction as per "PGA Illuminated", E. Lengyel:
//
// operator<<(a,b) = l_contract(a,b) = rwdg( l_dual(a), b )
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// right contractions A >> B: "A contracted by B"
//
// The resulting object lies in A and is perpendicular to B
//
// Implements the right bulk contraction as per "PGA Illuminated", E. Lengyel:
//
// operator>>(a,b) = r_contract(a,b) = rwdg( a, r_dual(b) )
//
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// commutator product (the asymmetric part of the geometric product)
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// regressive commutator product
// (the asymmetric part of the regressive geometric product)
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// geometric products
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// regressive geometric products
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// multiplicative inverses of scalars, blades and multivectors
// w.r.t. the geometric product:
// for k-blades: A^(-1) = rev(A)/|A|^2 = (-1)^(k*(k-1)/2)*A/|A|^2
// pattern for k = 0, 1, 2, 3, ...: + + - - + + - - ... (from reversion)
////////////////////////////////////////////////////////////////////////////////
// HINT: inv() cannot be constexpr due to the checks for division by zero
//       which might throw

template <typename T>
    requires(numeric_type<T>)
inline Scalar4ds<T> inv(Scalar4ds<T> s)
{
    T sq_n = bulk_nrm_sq(s);
    hd::ga::detail::check_normalization<T>(sq_n, "scalar");
    T inv = T(1.0) / sq_n;

    return Scalar4ds<T>(rev(s) * inv);
}

template <typename T>
    requires(numeric_type<T>)
inline Vec4ds<T> inv(Vec4ds<T> const& v)
{
    // v^(-1) = rev(v)/|v|^2 = v/dot(v,v) = v/bulk_sq_nrm(v)
    // using rev(v) = (-1)^[k(k-1)/2] v for a k-blade: 1-blade => rev(v) = v
    T sq_n = bulk_nrm_sq(v);
    hd::ga::detail::check_normalization<T>(sq_n, "vector");
    T inv = T(1.0) / sq_n; // inverse of squared norm for a vector
    return Vec4ds<T>(v.x * inv, v.y * inv, v.z * inv, v.w * inv);
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline BiVec4ds<T> inv(BiVec4ds<T> const& B)
{
    auto bc = B * conj(B);
    auto bcmap = gr0(bc) + gr2(bc) - gr4(bc);
    // fmt::println("B={}", B);
    // fmt::println("bc={}", bc);
    // fmt::println("bcmap={}", bcmap);
    // fmt::println("bc*bcmap={}", bc * bcmap);
    T sq_n = T(gr0(bc * bcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "bivector");
    return gr2(conj(B) * bcmap) / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline TriVec4ds<T> inv(TriVec4ds<T> const& t)
{
    auto tc = t * conj(t);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    // fmt::println("t={}", t);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "trivector");
    return gr3(conj(t) * tcmap) / sq_n;
}

// TODO: invert the pseudoscalar

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_E<T> inv(MVec4ds_E<T> const& E)
{
    auto tc = E * conj(E);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    // fmt::println("E={}", E);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "even-grade multivector");
    return conj(E) * tcmap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds_U<T> inv(MVec4ds_U<T> const& U)
{
    auto tc = U * conj(U);
    auto tcmap = gr0(tc) + gr2(tc) - gr4(tc);
    // fmt::println("E={}", E);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "odd-grade multivector");
    return conj(U) * tcmap / sq_n;
}

// formula from "Multivector and multivector matrix inverses in real Cliﬀord
// algebras", Hitzer, Sangwine, 2016 left and a right inverse are the same (see paper
// of Hitzer, Sangwine)
template <typename T>
    requires(numeric_type<T>)
inline MVec4ds<T> inv(MVec4ds<T> const& M)
{
    auto tc = M * conj(M);
    auto tcmap = gr0(tc) + gr1(tc) + gr2(tc) - gr3(tc) - gr4(tc);
    // fmt::println("M={}", M);
    // fmt::println("tc={}", tc);
    // fmt::println("tcmap={}", tcmap);
    // fmt::println("tc*tcmap={}", tc * tcmap);
    T sq_n = T(gr0(tc * tcmap));
    hd::ga::detail::check_normalization<T>(sq_n, "multivector");
    return conj(M) * tcmap / sq_n;
}

} // namespace hd::ga::sta