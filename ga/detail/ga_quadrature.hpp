#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <array>   // std::array
#include <cmath>   // std::abs, std::cos
#include <cstddef> // size_t

#include "../ga_usr_utilities.hpp" // pi
#include "../ga_value_t.hpp"       // value_t

/////////////////////////////////////////////////////////////////////////////////////////
// Gauss-Legendre quadrature nodes, shared by the headers that integrate along a curve or
// over a surface (the geodesic integrals, the closed-surface rules of the
// magnetostatics). Domain-agnostic: no GA here.
//
// provides in namespace hd::ga::detail:
//
// - gauss_legendre12_t, gauss_legendre12(): the 12-point rule on [-1, 1], computed once
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga::detail {

// 12-point Gauss-Legendre nodes and weights on [-1, 1], by Newton's method on the
// Legendre polynomial P_12; computed once
struct gauss_legendre12_t {
    std::array<value_t, 12> x{};
    std::array<value_t, 12> w{};
};

inline gauss_legendre12_t const& gauss_legendre12()
{
    static gauss_legendre12_t const g = [] {
        gauss_legendre12_t r{};
        int const n = 12;
        for (int i = 0; i < n / 2; ++i) {
            value_t z = std::cos(pi * (i + 0.75) / (n + 0.5));
            value_t dp = 1.0;
            for (int it = 0; it < 100; ++it) {
                value_t p1 = 1.0, p2 = 0.0;
                for (int j = 0; j < n; ++j) {
                    value_t const p3 = p2;
                    p2 = p1;
                    p1 = ((2.0 * j + 1.0) * z * p2 - j * p3) / (j + 1.0);
                }
                dp = n * (z * p1 - p2) / (z * z - 1.0); // P_12'(z)
                value_t const dz = p1 / dp;
                z -= dz;
                if (std::abs(dz) < 1.0e-16) break;
            }
            r.x[size_t(i)] = -z;
            r.x[size_t(n - 1 - i)] = z;
            r.w[size_t(i)] = r.w[size_t(n - 1 - i)] = 2.0 / ((1.0 - z * z) * dp * dp);
        }
        return r;
    }();
    return g;
}

} // namespace hd::ga::detail
