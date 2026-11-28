#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// THE TRANSWEDGE DECOMPOSITION of the geometric product (Lengyel 2025), checked from an
// algebra's own primitives on every pair of basis blades:
//
//     a * b = sum_k (-1)^(k(k-1)/2) a twdg_k b,
//
//     a twdg_k b = sum_{c of grade k} wdg( rwdg(l_cmpl(c), a), rwdg(b, r_dual(c)) )
//
// (ga_docu.pdf, "Transwedge products"). The metric enters only through the dual -- in
// PGA the bulk dual. Where an algebra has a regressive geometric product, its
// decomposition is the complement of the decomposition of the complements,
//
//     rgpr(a, b) = r_cmpl( l_cmpl(a) * l_cmpl(b) )   with * decomposed as above.
//
// An algebra passes its operations as a struct of static functions: lc / rc (the left
// and the right complement, cmpl in an odd dimension), rd (the right dual), w / rw (the
// wedge and the regressive wedge), gp (the geometric product) and, if it has one, rgp.
// The full multivector stores its components grade by grade, scalar first.
/////////////////////////////////////////////////////////////////////////////////////////

#include "ga/ga_value_t.hpp"

#include <cstddef>
#include <utility>

namespace transwedge {

using hd::ga::value_t;

constexpr std::size_t binomial(std::size_t n, std::size_t k)
{
    std::size_t r = 1;
    for (std::size_t j = 1; j <= k; ++j)
        r = r * (n - k + j) / j;
    return r;
}

// the grade of the i-th component of a full multivector in n dimensions
constexpr std::size_t grade_of(std::size_t i, std::size_t n)
{
    std::size_t start = 0;
    for (std::size_t g = 0; g <= n; ++g) {
        start += binomial(n, g);
        if (i < start) return g;
    }
    return n;
}

// the i-th basis blade as a full multivector
template <typename M, std::size_t... I>
M basis_blade_impl(std::size_t i, std::index_sequence<I...>)
{
    return M((I == i ? value_t(1.0) : value_t(0.0))...);
}

template <typename M, std::size_t N> M basis_blade(std::size_t i)
{
    return basis_blade_impl<M>(i, std::make_index_sequence<N>{});
}

// sum_k sign(k) a twdg_k b. The two flags are the falsifiers: flip uses the sign pattern
// (-1)^(k(k+1)/2), cmpl_for_dual the right complement in place of the right dual (which
// only an algebra with the identity metric cannot tell apart)
template <typename M, std::size_t n, typename Ops>
M twdg_sum(M const& a, M const& b, bool flip = false, bool cmpl_for_dual = false)
{
    constexpr std::size_t N = std::size_t(1) << n;
    M sum{};
    for (std::size_t k = 0; k <= n; ++k) {
        M part{};
        for (std::size_t i = 0; i < N; ++i) {
            if (grade_of(i, n) != k) continue;
            M const c = basis_blade<M, N>(i);
            M const dc = cmpl_for_dual ? M(Ops::rc(c)) : M(Ops::rd(c));
            part = part + M(Ops::w(M(Ops::rw(M(Ops::lc(c)), a)), M(Ops::rw(b, dc))));
        }
        std::size_t const e = flip ? k * (k + 1) / 2 : (k == 0 ? 0 : k * (k - 1) / 2);
        sum = (e % 2 == 0) ? M(sum + part) : M(sum - part);
    }
    return sum;
}

// the number of basis pairs on which the decomposition differs from the product
template <typename M, std::size_t n, typename Ops>
std::size_t gpr_mismatches(bool flip = false, bool cmpl_for_dual = false)
{
    constexpr std::size_t N = std::size_t(1) << n;
    std::size_t bad = 0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            M const a = basis_blade<M, N>(i);
            M const b = basis_blade<M, N>(j);
            if (!(twdg_sum<M, n, Ops>(a, b, flip, cmpl_for_dual) == M(Ops::gp(a, b))))
                ++bad;
        }
    }
    return bad;
}

template <typename M, std::size_t n, typename Ops> std::size_t rgpr_mismatches()
{
    constexpr std::size_t N = std::size_t(1) << n;
    std::size_t bad = 0;
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            M const a = basis_blade<M, N>(i);
            M const b = basis_blade<M, N>(j);
            M const lca = Ops::lc(a);
            M const lcb = Ops::lc(b);
            if (!(M(Ops::rc(twdg_sum<M, n, Ops>(lca, lcb))) == M(Ops::rgp(a, b)))) ++bad;
        }
    }
    return bad;
}

} // namespace transwedge
