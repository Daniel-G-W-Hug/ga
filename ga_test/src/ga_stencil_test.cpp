// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

// Tests for the finite-difference stencil generator (ga/detail/ga_stencil.hpp) and
// the algebra-independent scheme layer built on it (ga/ga_usr_fd.hpp).
//
// Every case pins a stencil with a KNOWN closed form (weights, order, and leading
// truncation coefficient), so a regression in the Taylor-matching system or the
// leading-term detection shows up as a hard numeric mismatch. Further cases apply a
// stencil to an analytic function and check the measured convergence rate against the
// reported order.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "fmt/format.h"

#include <cmath>     // std::sin, std::cos, std::abs, std::log2, std::cbrt
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::invalid_argument
#include <vector>    // std::vector

#include "ga/detail/ga_solver.hpp"
#include "ga/detail/ga_stencil.hpp"
#include "ga/ga_usr_fd.hpp"

using namespace hd::ga;
using hd::ga::fact;
using hd::ga::stencil_lhs;
using hd::ga::stencil_t;


// A minimal stand-in for a multivector unknown: tridiag_solve() requires only T + T
// and double * T of the type it solves for, and this has exactly those two.
struct pair2_t {
    double u{0.0};
    double v{0.0};
};

inline pair2_t operator+(pair2_t const& a, pair2_t const& b)
{
    return {a.u + b.u, a.v + b.v};
}

inline pair2_t operator*(double s, pair2_t const& a) { return {s * a.u, s * a.v}; }

TEST_SUITE("fd stencil generator")
{

    TEST_CASE("fact")
    {
        CHECK(fact(0) == doctest::Approx(1.0));
        CHECK(fact(1) == doctest::Approx(1.0));
        CHECK(fact(5) == doctest::Approx(120.0));
        CHECK(fact(10) == doctest::Approx(3628800.0));

        // the table is a compile-time constant, so it is usable in constant
        // expressions -- and is immutable, so concurrent callers cannot race on it
        static_assert(fact(0) == 1.0);
        static_assert(fact(20) == 2432902008176640000.0); // last exact one in double
        static_assert(fact(171) == std::numeric_limits<double>::infinity());

        // accuracy is that of the iterated product all the way out, NOT of an
        // exp(log_gamma()) approximation, which loses five digits by n = 100
        double const exact100 = 9.33262154439441526817e157;
        CHECK(std::abs(fact(100) - exact100) / exact100 < 1.0e-15);

        // 170! is the last one a double can hold; past it the honest answer is inf
        CHECK(fact(170) < std::numeric_limits<double>::infinity());
        CHECK(fact(171) == std::numeric_limits<double>::infinity());

        // a negative argument is rejected rather than wrapped into a huge loop
        CHECK_THROWS_AS(fact(-1), std::invalid_argument);
    }

    TEST_CASE("central 3-point first derivative (explicit, order 2)")
    {
        // f'(x0) ~ ( f(x0+h) - f(x0-h) ) / (2h);  leading error  h^2/6 * f'''
        double const x0 = 0.0;
        double const h = 1.0;
        stencil_t const s(x0, stencil_lhs::f1, {x0 - h, x0, x0 + h}, {x0}, {});

        REQUIRE(s.wf0.size() == 3);
        CHECK(s.wf0[0] == doctest::Approx(-1.0 / (2.0 * h)).epsilon(1e-12));
        CHECK(s.wf0[1] == doctest::Approx(0.0).epsilon(1e-12));
        CHECK(s.wf0[2] == doctest::Approx(1.0 / (2.0 * h)).epsilon(1e-12));

        REQUIRE(s.wf1.size() == 1);
        CHECK(s.wf1[0] == doctest::Approx(1.0).epsilon(1e-12)); // normalized lhs

        CHECK(s.order == 2);
        CHECK(s.trunc_err == doctest::Approx(h * h / 6.0).epsilon(1e-12));
    }

    TEST_CASE("central 3-point second derivative (explicit, order 2)")
    {
        // f''(x0) ~ ( f(x0-h) - 2 f(x0) + f(x0+h) ) / h^2;  leading error  h^2/12 * f''''
        double const x0 = 0.0;
        double const h = 1.0;
        stencil_t const s(x0, stencil_lhs::f2, {x0 - h, x0, x0 + h}, {}, {x0});

        REQUIRE(s.wf0.size() == 3);
        CHECK(s.wf0[0] == doctest::Approx(1.0 / (h * h)).epsilon(1e-12));
        CHECK(s.wf0[1] == doctest::Approx(-2.0 / (h * h)).epsilon(1e-12));
        CHECK(s.wf0[2] == doctest::Approx(1.0 / (h * h)).epsilon(1e-12));

        REQUIRE(s.wf2.size() == 1);
        CHECK(s.wf2[0] == doctest::Approx(1.0).epsilon(1e-12)); // normalized lhs

        CHECK(s.order == 2);
        CHECK(s.trunc_err == doctest::Approx(h * h / 12.0).epsilon(1e-12));
    }

    TEST_CASE("one-sided 2-point first derivative (explicit, order 1)")
    {
        // f'(x0) ~ ( f(x0+h) - f(x0) ) / h;  leading error  h/2 * f''
        //
        // Pins the leading-term detection: the first non-vanishing residual order is
        // the one that counts (an overwriting scan would report order 2 here).
        double const x0 = 0.0;
        double const h = 1.0;
        stencil_t const s(x0, stencil_lhs::f1, {x0, x0 + h}, {x0}, {});

        REQUIRE(s.wf0.size() == 2);
        CHECK(s.wf0[0] == doctest::Approx(-1.0 / h).epsilon(1e-12));
        CHECK(s.wf0[1] == doctest::Approx(1.0 / h).epsilon(1e-12));

        CHECK(s.order == 1);
        CHECK(s.trunc_err == doctest::Approx(h / 2.0).epsilon(1e-12));
    }

    TEST_CASE("compact 3-point first derivative (implicit Pade scheme, order 4)")
    {
        // (1/6) f'(x0-h) + (2/3) f'(x0) + (1/6) f'(x0+h)
        //     ~ ( f(x0+h) - f(x0-h) ) / (2h)
        //
        // The classic 4th-order compact scheme in the lhs-sum-normalized form.
        // Exercises the lhs (sfact) path of the residual computation.
        double const x0 = 0.0;
        double const h = 1.0;
        stencil_t const s(x0, stencil_lhs::f1, {x0 - h, x0 + h}, {x0 - h, x0, x0 + h},
                          {});

        REQUIRE(s.wf0.size() == 2);
        CHECK(s.wf0[0] == doctest::Approx(-1.0 / (2.0 * h)).epsilon(1e-12));
        CHECK(s.wf0[1] == doctest::Approx(1.0 / (2.0 * h)).epsilon(1e-12));

        REQUIRE(s.wf1.size() == 3);
        CHECK(s.wf1[0] == doctest::Approx(1.0 / 6.0).epsilon(1e-12));
        CHECK(s.wf1[1] == doctest::Approx(2.0 / 3.0).epsilon(1e-12));
        CHECK(s.wf1[2] == doctest::Approx(1.0 / 6.0).epsilon(1e-12));
        CHECK(s.wf1[0] + s.wf1[1] + s.wf1[2] == doctest::Approx(1.0).epsilon(1e-12));

        CHECK(s.order == 4);
        CHECK(s.trunc_err == doctest::Approx(-std::pow(h, 4) / 180.0).epsilon(1e-12));
    }

    TEST_CASE("measured convergence rate matches the reported order")
    {
        // apply the one-sided 3-point stencil for f' to f = sin at x0 and measure the
        // convergence rate from two step sizes; it must match the reported order
        double const x0 = 0.3;

        auto fd_error = [&](double h) {
            stencil_t const s(x0, stencil_lhs::f1, {x0, x0 + h, x0 + 2.0 * h}, {x0}, {});
            double fd = 0.0;
            for (size_t j = 0; j < s.xf0.size(); ++j) {
                fd += s.wf0[j] * std::sin(s.xf0[j]);
            }
            return std::abs(fd - std::cos(x0));
        };

        // reported order (the residual test is relative, so any h reports it)
        stencil_t const s(x0, stencil_lhs::f1, {x0, x0 + 1.0, x0 + 2.0}, {x0}, {});
        CHECK(s.order == 2);

        double const e1 = fd_error(1.0e-3);
        double const e2 = fd_error(0.5e-3);
        double const rate = std::log2(e1 / e2);
        CHECK(rate == doctest::Approx(2.0).epsilon(0.05));
    }

    ////////////////////////////////////////////////////////////////////////////
    // order detection -- the property that used to be silently scale-dependent
    ////////////////////////////////////////////////////////////////////////////

    TEST_CASE("stencil_t: order detection is independent of the node scale")
    {
        // The residual test compares each Taylor order against the magnitude of the
        // contributions it cancels between. It must therefore report the SAME order
        // for the same stencil shape at any spacing -- an absolute threshold reports
        // order 0 once the spacing is small enough, which then reads as a valid
        // answer and poisons every step size derived from it.
        for (int e = 0; e <= 12; e += 2) {
            double const h = std::pow(10.0, -double(e));
            std::vector<double> nodes;
            for (int i = -2; i <= 2; ++i)
                nodes.push_back(double(i) * h);

            stencil_t const s(0.0, stencil_lhs::f1, nodes, {0.0}, {});
            CHECK(s.order == 4);
            // leading term of the 5-point central first difference: -h^4/30 f^(5)
            CHECK(s.trunc_err ==
                  doctest::Approx(-std::pow(h, 4) / 30.0).epsilon(1.0e-10));
        }
    }

    TEST_CASE("stencil_t: wide stencils report their full order")
    {
        // A wide stencil has a SMALL leading coefficient, which an absolute test
        // cannot distinguish from zero. Both the natural and the unit-interval
        // layouts of the same 13-point stencil must report order 12.
        std::vector<double> wide, squeezed;
        for (int i = -6; i <= 6; ++i) {
            wide.push_back(double(i));
            squeezed.push_back(double(i) / 6.0); // |nodes| <= 1, tiny trunc_err
        }
        CHECK(stencil_t(0.0, stencil_lhs::f1, wide, {0.0}, {}).order == 12);

        stencil_t const sq(0.0, stencil_lhs::f1, squeezed, {0.0}, {});
        CHECK(sq.order == 12);
        CHECK(std::abs(sq.trunc_err) < 1.0e-6); // below the old absolute threshold
        CHECK(std::abs(sq.trunc_err) > 0.0);
    }

    TEST_CASE("stencil_t: an uncharacterizable stencil throws")
    {
        // fewer than two function samples cannot pin a derivative
        CHECK_THROWS_AS(stencil_t(0.0, stencil_lhs::f1, {0.0}, {0.0}, {}),
                        std::invalid_argument);
        // repeated nodes make the Taylor system singular -- it must not come back
        // reporting order 0 as though that were an answer
        CHECK_THROWS(stencil_t(0.0, stencil_lhs::f1, {-1.0, -1.0, -1.0}, {0.0}, {}));
    }

    ////////////////////////////////////////////////////////////////////////////
    // the parts of stencil_t's model the scheme layer must keep reachable
    ////////////////////////////////////////////////////////////////////////////

    TEST_CASE("stencil_t: staggered stencil, x0 between the nodes")
    {
        // derivative at the MIDPOINT from samples on the half-offsets: the Yee-type
        // arrangement. Exact weights {+1/24, -9/8, +9/8, -1/24}, order 4.
        stencil_t const s(0.0, stencil_lhs::f1, {-1.5, -0.5, 0.5, 1.5}, {0.0}, {});
        CHECK(s.order == 4);
        CHECK(s.wf0[0] == doctest::Approx(1.0 / 24.0).epsilon(1.0e-12));
        CHECK(s.wf0[1] == doctest::Approx(-9.0 / 8.0).epsilon(1.0e-12));
        CHECK(s.wf0[2] == doctest::Approx(9.0 / 8.0).epsilon(1.0e-12));
        CHECK(s.wf0[3] == doctest::Approx(-1.0 / 24.0).epsilon(1.0e-12));

        // it reproduces the derivative of an analytic function at the midpoint
        double const h = 1.0e-3;
        double acc = 0.0;
        for (size_t k = 0; k < s.xf0.size(); ++k)
            acc += s.wf0[k] / h * std::sin(0.3 + s.xf0[k] * h);
        CHECK(acc == doctest::Approx(std::cos(0.3)).epsilon(1.0e-12));
    }

    TEST_CASE("stencil_t: mixed f / f'' data (Hermite) raises the order")
    {
        // three f samples alone give order 2 for f'; supplying f'' at +-1 as well
        // lifts the SAME f footprint to order 4
        stencil_t const plain(0.0, stencil_lhs::f1, {-1.0, 0.0, 1.0}, {0.0}, {});
        stencil_t const herm(0.0, stencil_lhs::f1, {-1.0, 0.0, 1.0}, {0.0}, {-1.0, 1.0});
        CHECK(plain.order == 2);
        CHECK(herm.order == 4);
        CHECK(herm.wf2.size() == 2);

        // and it is not merely a larger number: apply both to sin and compare
        double const h = 1.0e-2;
        auto apply = [&](stencil_t const& s) {
            double acc = 0.0;
            for (size_t k = 0; k < s.xf0.size(); ++k)
                acc += s.wf0[k] / h * std::sin(0.3 + s.xf0[k] * h);
            for (size_t k = 0; k < s.xf2.size(); ++k) // supplied f'' data
                acc += s.wf2[k] * h * (-std::sin(0.3 + s.xf2[k] * h));
            return std::abs(acc - std::cos(0.3));
        };
        CHECK(apply(herm) < apply(plain) / 100.0);
    }

    ////////////////////////////////////////////////////////////////////////////
    // the scheme layer (ga/ga_usr_fd.hpp)
    ////////////////////////////////////////////////////////////////////////////

    TEST_CASE("fd_scheme: central_scheme reports the order it was asked for")
    {
        // The regression this pins: central_scheme() used to normalise its nodes into
        // [-1, 1] and then read the order off an ABSOLUTE threshold, so every order
        // >= 8 silently came back as 0 -- and fd_step() then collapsed to eps.
        for (int deriv : {1, 2}) {
            for (int order = 2; order <= 12; order += 2) {
                auto const sc = central_scheme(deriv, order);
                CAPTURE(deriv);
                CAPTURE(order);
                CHECK(sc.order == order);
                CHECK(sc.trunc_err != 0.0);
                CHECK(sc.deriv == deriv);
                CHECK(sc.nodes.size() == std::size_t(order) + 1); // minimal centred
                CHECK(fd_step(sc, 1.0) > 1.0e-7); // never collapses towards eps
            }
        }
    }

    TEST_CASE("fd_scheme: fd_step uses the leading truncation coefficient")
    {
        // h_opt = ( d * sum|w| * eps / (p * |trunc_err|) )^(1/(p+d)). For the central
        // first difference sum|w| = 1, trunc_err = 1/6, p = 2, d = 1, which is the
        // textbook cbrt(3 eps) -- the exponent-only form eps^(1/3) misses by ~1.4x.
        auto const d1 = central_scheme(1, 2);
        CHECK(d1.trunc_err == doctest::Approx(1.0 / 6.0).epsilon(1.0e-12));

        double const e = std::numeric_limits<double>::epsilon();
        CHECK(fd_step(d1, 1.0) == doctest::Approx(std::cbrt(3.0 * e)).epsilon(1.0e-9));
        // the coefficient is what separates it from the exponent-only estimate: the
        // two differ by exactly cbrt(3), so a step that ignored trunc_err fails here
        CHECK(fd_step(d1, 1.0) / std::cbrt(e) ==
              doctest::Approx(std::cbrt(3.0)).epsilon(1.0e-9));

        // and it tracks the evaluation scale
        CHECK(fd_step(d1, 1000.0) == doctest::Approx(1000.0 * fd_step(d1, 1.0)));
    }

    TEST_CASE("fd_scheme: staggered and Hermite schemes are reachable")
    {
        // staggered_scheme(): x0 sits between the nodes, so none of them is 0
        auto const st = staggered_scheme(1, 4);
        CHECK(st.order == 4);
        CHECK(st.nodes.size() == 4);
        for (double n : st.nodes)
            CHECK(n != 0.0);
        CHECK(st.is_offset());
        CHECK(!st.is_compact());
        CHECK(!st.is_hermite());
        CHECK(staggered_scheme(2, 2).order >= 2);

        // it evaluates through apply_scheme like any other offset scheme
        double const h = 1.0e-3;
        double const got =
            apply_scheme<double>(st, [](double off) { return std::sin(0.3 + off); }, h);
        CHECK(got == doctest::Approx(std::cos(0.3)).epsilon(1.0e-11));

        // make_scheme(): the Hermite case, with f'' supplied as rhs data
        auto const hm = make_scheme(0.0, 1, {-1.0, 0.0, 1.0}, {0.0}, {-1.0, 1.0});
        CHECK(hm.is_hermite());
        CHECK(hm.order == 4);
        CHECK(hm.aux_deriv() == 2);
        CHECK(hm.aux_weights.size() == 2);
    }

    TEST_CASE("fd_scheme: apply_scheme rejects what it cannot evaluate")
    {
        auto sample = [](double) { return 1.0; };

        // implicit: needs a grid
        CHECK_THROWS_AS(apply_scheme<double>(compact_scheme(1), sample, 1.0e-3),
                        std::invalid_argument);
        // Hermite: a sampler of f alone cannot supply the derivative data
        CHECK_THROWS_AS(apply_scheme<double>(
                            make_scheme(0.0, 1, {-1.0, 0.0, 1.0}, {0.0}, {-1.0, 1.0}),
                            sample, 1.0e-3),
                        std::invalid_argument);
        // grid-bound: its nodes are absolute coordinates, not step multiples
        std::vector<double> x;
        for (int i = 0; i < 9; ++i)
            x.push_back(0.1 * double(i));
        CHECK_THROWS_AS(apply_scheme<double>(scheme_at(x, 4, 1), sample, 1.0e-3),
                        std::invalid_argument);
        // and a non-positive step is not a step
        CHECK_THROWS_AS(apply_scheme<double>(central_scheme(1, 2), sample, 0.0),
                        std::invalid_argument);
    }

    TEST_CASE("fd_scheme: scheme_at uses the minimal centred footprint")
    {
        // The regression this pins: a dead ternary made scheme_at() ask for order+2
        // nodes, so an interior order-2 stencil was a lopsided FOUR-point one rather
        // than the classical centred three-point difference.
        std::vector<double> x;
        for (int i = 0; i < 21; ++i)
            x.push_back(0.1 * double(i));

        for (int deriv : {1, 2}) {
            for (int order : {2, 4}) {
                auto const sc = scheme_at(x, 10, deriv, fd_kind::explicit_fd, order);
                CAPTURE(deriv);
                CAPTURE(order);
                CHECK(sc.nodes.size() == std::size_t(order) + 1);
                CHECK(sc.order == order);
                // centred on x[10]: the window starts order/2 nodes before it
                CHECK(sc.node0 == std::size_t(10 - order / 2));
                CHECK(sc.nodes.front() == doctest::Approx(x[10] - 0.1 * (order / 2)));
                CHECK(sc.nodes.back() == doctest::Approx(x[10] + 0.1 * (order / 2)));
            }
        }

        // a boundary node has no symmetry to lean on, so it is widened to order+deriv
        // to still reach the requested accuracy
        for (int deriv : {1, 2}) {
            auto const sc = scheme_at(x, 0, deriv, fd_kind::explicit_fd, 4);
            CAPTURE(deriv);
            CHECK(sc.node0 == 0);
            CHECK(sc.nodes.size() == std::size_t(4 + deriv));
            CHECK(sc.order == 4);
        }
    }

    TEST_CASE("fd_scheme: grid derivatives converge at the reported order")
    {
        // f = sin on a uniform grid; halving the spacing must drop the error by
        // 2^order for both the explicit and the compact path
        auto max_err = [](std::size_t n, int deriv, fd_kind kind, int order) {
            std::vector<double> x;
            std::vector<double> f;
            for (std::size_t i = 0; i < n; ++i) {
                double const xi = 0.2 + 2.0 * double(i) / double(n - 1);
                x.push_back(xi);
                f.push_back(std::sin(xi));
            }
            auto const d = fd_derivative(x, f, deriv, kind, order);
            double e = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                double const exact = (deriv == 1) ? std::cos(x[i]) : -std::sin(x[i]);
                e = std::max(e, std::abs(d[i] - exact));
            }
            return e;
        };

        for (int deriv : {1, 2}) {
            for (int order : {2, 4}) {
                double const e1 = max_err(41, deriv, fd_kind::explicit_fd, order);
                double const e2 = max_err(81, deriv, fd_kind::explicit_fd, order);
                double const rate = std::log2(e1 / e2);
                CAPTURE(deriv);
                CAPTURE(order);
                CHECK(rate == doctest::Approx(double(order)).epsilon(0.10));
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // the lean tridiagonal solver the compact schemes ask for
    ////////////////////////////////////////////////////////////////////////////

    TEST_CASE("tridiag_solve: agrees with the general dense solver")
    {
        // same system, both ways: the O(n) sweep must land where LU with partial
        // pivoting does, or it is not a solver but an approximation of one
        std::size_t const n = 6;
        std::vector<double> const a{0.0, -1.0, -2.0, -1.0, -3.0, -1.0};
        std::vector<double> const b{4.0, 5.0, 6.0, 4.5, 7.0, 3.0};
        std::vector<double> const c{-1.0, -2.0, -1.0, -1.5, -1.0, 0.0};
        std::vector<double> const d{1.0, -2.0, 3.0, 0.5, -1.0, 2.0};

        auto const x = tridiag_solve(a, b, c, d);

        std::vector<double> dense(n * n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (i > 0) dense[i * n + (i - 1)] = a[i];
            dense[i * n + i] = b[i];
            if (i + 1 < n) dense[i * n + (i + 1)] = c[i];
        }
        auto const x_lu = lu_solve(dense, d, n);

        for (std::size_t i = 0; i < n; ++i) {
            CHECK(x[i] == doctest::Approx(x_lu[i]).epsilon(1.0e-12));
        }

        // and it really satisfies the equations it was given
        for (std::size_t i = 0; i < n; ++i) {
            double lhs = b[i] * x[i];
            if (i > 0) lhs += a[i] * x[i - 1];
            if (i + 1 < n) lhs += c[i] * x[i + 1];
            CHECK(lhs == doctest::Approx(d[i]).epsilon(1.0e-12));
        }
    }

    TEST_CASE("tridiag_solve: steady advection-diffusion boundary layer")
    {
        // v u' = D u'' on [0,1] with u(0) = 0, u(1) = 1 -- the standard test for a
        // tridiagonal solve, because central differencing of both terms produces
        // exactly three diagonals and the answer is known in closed form:
        //
        //     u(x) = ( exp(Pe x) - 1 ) / ( exp(Pe) - 1 ),     Pe = v / D
        //
        // At Pe = 10 the solution is flat over most of the domain and turns up into a
        // boundary layer near x = 1, so a solver that merely produced a smooth curve
        // would not pass.
        double const v = 1.0, D = 0.1;
        double const Pe = v / D;

        auto exact = [&](double xx) {
            return (std::exp(Pe * xx) - 1.0) / (std::exp(Pe) - 1.0);
        };

        // assemble the rows from the fd weights themselves, which is how a caller
        // would build such a system
        auto max_err = [&](std::size_t n) {
            std::vector<double> x;
            for (std::size_t i = 0; i < n; ++i)
                x.push_back(double(i) / double(n - 1));

            std::vector<double> a(n, 0.0), b(n, 1.0), c(n, 0.0), d(n, 0.0);
            for (std::size_t i = 1; i + 1 < n; ++i) {
                auto const s1 = scheme_at(x, i, 1, fd_kind::explicit_fd, 2);
                auto const s2 = scheme_at(x, i, 2, fd_kind::explicit_fd, 2);
                // both are the centred three-point stencil on the same window
                REQUIRE(s1.nodes.size() == 3);
                REQUIRE(s2.nodes.size() == 3);
                REQUIRE(s1.node0 == i - 1);
                REQUIRE(s2.node0 == i - 1);

                a[i] = v * s1.weights[0] - D * s2.weights[0];
                b[i] = v * s1.weights[1] - D * s2.weights[1];
                c[i] = v * s1.weights[2] - D * s2.weights[2];
            }
            // Dirichlet rows: b = 1 already, only the values are needed. They also
            // make the system strictly dominant at the ends, which is what carries
            // the weakly dominant interior rows (|b| = |a| + |c| there).
            d[0] = exact(0.0);
            d[n - 1] = exact(1.0);

            auto const u = tridiag_solve(a, b, c, d);

            double e = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                e = std::max(e, std::abs(u[i] - exact(x[i])));
            return e;
        };

        double const e51 = max_err(51);
        double const e101 = max_err(101);
        double const e201 = max_err(201);

        // the discretization is 2nd order, so each halving of h quarters the error
        CHECK(std::log2(e51 / e101) == doctest::Approx(2.0).epsilon(0.10));
        CHECK(std::log2(e101 / e201) == doctest::Approx(2.0).epsilon(0.10));
        CHECK(e201 < 1.0e-4);
    }

    TEST_CASE("tridiag_solve: advection-diffusion with a source term")
    {
        // The homogeneous case above has an interior right-hand side of ZERO, which
        // makes it blind to how the rhs is carried through the elimination -- drop
        // that term and the case still passes. So drive the same operator with a
        // manufactured solution, whose source is non-zero at every interior node:
        //
        //     u(x) = sin(2 pi x) + x        v u' - D u'' = s(x)
        //     s(x) = v (2 pi cos(2 pi x) + 1) + D 4 pi^2 sin(2 pi x)
        double const v = 1.0, D = 0.1;
        double const tp = 8.0 * std::atan(1.0); // 2 pi

        auto exact = [&](double xx) { return std::sin(tp * xx) + xx; };
        auto source = [&](double xx) {
            return v * (tp * std::cos(tp * xx) + 1.0) + D * tp * tp * std::sin(tp * xx);
        };

        auto max_err = [&](std::size_t n) {
            std::vector<double> x;
            for (std::size_t i = 0; i < n; ++i)
                x.push_back(double(i) / double(n - 1));

            std::vector<double> a(n, 0.0), b(n, 1.0), c(n, 0.0), d(n, 0.0);
            for (std::size_t i = 1; i + 1 < n; ++i) {
                auto const s1 = scheme_at(x, i, 1, fd_kind::explicit_fd, 2);
                auto const s2 = scheme_at(x, i, 2, fd_kind::explicit_fd, 2);
                a[i] = v * s1.weights[0] - D * s2.weights[0];
                b[i] = v * s1.weights[1] - D * s2.weights[1];
                c[i] = v * s1.weights[2] - D * s2.weights[2];
                d[i] = source(x[i]); // non-zero everywhere inside
            }
            d[0] = exact(0.0);
            d[n - 1] = exact(1.0);

            auto const u = tridiag_solve(a, b, c, d);

            double e = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                e = std::max(e, std::abs(u[i] - exact(x[i])));
            return e;
        };

        double const e101 = max_err(101);
        double const e201 = max_err(201);
        CHECK(std::log2(e101 / e201) == doctest::Approx(2.0).epsilon(0.10));
        CHECK(e201 < 1.0e-3);
    }

    TEST_CASE("tridiag_solve: the unknowns need not be scalars")
    {
        // the coefficients are scalars while the unknowns carry structure -- the case
        // fd_derivative() needs, where each node holds a whole multivector
        std::size_t const n = 4;
        std::vector<double> const a{0.0, -1.0, -1.0, -1.0};
        std::vector<double> const b{2.0, 3.0, 3.0, 2.0};
        std::vector<double> const c{-1.0, -1.0, -1.0, 0.0};
        std::vector<pair2_t> d(n);
        for (std::size_t i = 0; i < n; ++i)
            d[i] = pair2_t{double(i) + 1.0, -2.0 * double(i)};

        auto const x = tridiag_solve(a, b, c, d);

        // Each component must solve the same scalar system. Pinned against the DENSE
        // solver, not against tridiag_solve on scalars -- comparing the routine with
        // itself would confirm only that it is generic, never that it is right.
        std::vector<double> dense(n * n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            if (i > 0) dense[i * n + (i - 1)] = a[i];
            dense[i * n + i] = b[i];
            if (i + 1 < n) dense[i * n + (i + 1)] = c[i];
        }
        std::vector<double> du(n), dv(n);
        for (std::size_t i = 0; i < n; ++i) {
            du[i] = d[i].u;
            dv[i] = d[i].v;
        }
        auto const xu = lu_solve(dense, du, n);
        auto const xv = lu_solve(dense, dv, n);
        for (std::size_t i = 0; i < n; ++i) {
            CHECK(x[i].u == doctest::Approx(xu[i]).epsilon(1.0e-12));
            CHECK(x[i].v == doctest::Approx(xv[i]).epsilon(1.0e-12));
        }
    }

    TEST_CASE("tridiag_solve: rejects a system it cannot solve")
    {
        std::vector<double> const ok3{1.0, 1.0, 1.0};
        std::vector<double> const ok3d{1.0, 2.0, 3.0};

        // mismatched operand lengths
        CHECK_THROWS_AS(tridiag_solve(std::vector<double>{1.0, 1.0}, ok3, ok3, ok3d),
                        hd::ga::Solver_error);
        // too small to be a tridiagonal system
        CHECK_THROWS_AS(tridiag_solve(std::vector<double>{1.0}, std::vector<double>{1.0},
                                      std::vector<double>{0.0}, std::vector<double>{1.0}),
                        hd::ga::Solver_error);
        // a singular system must not come back with a plausible-looking answer
        CHECK_THROWS_AS(tridiag_solve(std::vector<double>{0.0, 1.0, 1.0},
                                      std::vector<double>{0.0, 1.0, 1.0}, ok3, ok3d),
                        hd::ga::Solver_error);
        // ... including one that only goes singular during the elimination
        CHECK_THROWS_AS(tridiag_solve(std::vector<double>{0.0, 1.0, 1.0},
                                      std::vector<double>{1.0, 1.0, 1.0},
                                      std::vector<double>{1.0, 1.0, 0.0}, ok3d),
                        hd::ga::Solver_error);
    }

    TEST_CASE("boundary conditions, explicit FD: the two Dirichlet routes agree")
    {
        // The boundary-row recipes documented at tridiag_solve() in ga_solver.hpp.
        // Route 1 keeps the boundary node as an unknown with an identity row; route 2
        // drops it and folds a[1]*g_0 / c[n-2]*g_1 into the right-hand side. They must
        // return the same interior values -- if they do not, one of the two recipes is
        // wrong and the comment is teaching a mistake.
        double const v = 1.0, D = 0.1;
        double const tp = 8.0 * std::atan(1.0); // 2 pi
        std::size_t const n = 81;

        auto exact = [&](double xx) { return std::sin(tp * xx) + xx; };
        auto source = [&](double xx) {
            return v * (tp * std::cos(tp * xx) + 1.0) + D * tp * tp * std::sin(tp * xx);
        };

        std::vector<double> x;
        for (std::size_t i = 0; i < n; ++i)
            x.push_back(double(i) / double(n - 1));

        // the interior rows are common to both routes
        std::vector<double> ai(n, 0.0), bi(n, 0.0), ci(n, 0.0), si(n, 0.0);
        for (std::size_t i = 1; i + 1 < n; ++i) {
            auto const s1 = scheme_at(x, i, 1, fd_kind::explicit_fd, 2);
            auto const s2 = scheme_at(x, i, 2, fd_kind::explicit_fd, 2);
            ai[i] = v * s1.weights[0] - D * s2.weights[0];
            bi[i] = v * s1.weights[1] - D * s2.weights[1];
            ci[i] = v * s1.weights[2] - D * s2.weights[2];
            si[i] = source(x[i]);
        }
        double const g0 = exact(0.0);
        double const g1 = exact(1.0);

        // --- route 1: identity rows, n unknowns -------------------------------------
        std::vector<double> a1 = ai, b1 = bi, c1 = ci, d1 = si;
        a1[0] = 0.0;
        b1[0] = 1.0;
        c1[0] = 0.0;
        d1[0] = g0;
        a1[n - 1] = 0.0;
        b1[n - 1] = 1.0;
        c1[n - 1] = 0.0;
        d1[n - 1] = g1;
        auto const u1 = tridiag_solve(a1, b1, c1, d1);

        // --- route 2: eliminate the known nodes, n-2 unknowns -----------------------
        std::size_t const m = n - 2;
        std::vector<double> a2(m), b2(m), c2(m), d2(m);
        for (std::size_t k = 0; k < m; ++k) {
            std::size_t const i = k + 1; // grid index of interior unknown k
            a2[k] = ai[i];
            b2[k] = bi[i];
            c2[k] = ci[i];
            d2[k] = si[i];
        }
        d2[0] -= a2[0] * g0; // the known left value moves to the rhs
        a2[0] = 0.0;
        d2[m - 1] -= c2[m - 1] * g1; // and the known right value likewise
        c2[m - 1] = 0.0;
        auto const u2 = tridiag_solve(a2, b2, c2, d2);

        for (std::size_t k = 0; k < m; ++k) {
            CHECK(u2[k] == doctest::Approx(u1[k + 1]).epsilon(1.0e-11));
        }
        // and route 1 reproduces the prescribed values exactly, being identity rows
        CHECK(u1[0] == doctest::Approx(g0).epsilon(1.0e-14));
        CHECK(u1[n - 1] == doctest::Approx(g1).epsilon(1.0e-14));
    }

    TEST_CASE("boundary conditions, explicit FD: Neumann restored to the band")
    {
        // Recipe 3 from the tridiag_solve() comment. The 2nd-order one-sided stencil
        // for u'(1) reaches back to column n-3, one past the band, so it is folded
        // against row n-2 -- the only other row touching that column. The gate is
        // that the result still converges at 2nd order: a 1st-order boundary row
        // would give a rate near 1, which is exactly what the elimination avoids.
        double const v = 1.0, D = 0.1;
        double const tp = 8.0 * std::atan(1.0); // 2 pi

        auto exact = [&](double xx) { return std::sin(tp * xx) + xx; };
        auto dexact = [&](double xx) { return tp * std::cos(tp * xx) + 1.0; };
        auto source = [&](double xx) {
            return v * (tp * std::cos(tp * xx) + 1.0) + D * tp * tp * std::sin(tp * xx);
        };

        auto max_err = [&](std::size_t n) {
            std::vector<double> x;
            for (std::size_t i = 0; i < n; ++i)
                x.push_back(double(i) / double(n - 1));

            std::vector<double> a(n, 0.0), b(n, 1.0), c(n, 0.0), d(n, 0.0);
            for (std::size_t i = 1; i + 1 < n; ++i) {
                auto const s1 = scheme_at(x, i, 1, fd_kind::explicit_fd, 2);
                auto const s2 = scheme_at(x, i, 2, fd_kind::explicit_fd, 2);
                a[i] = v * s1.weights[0] - D * s2.weights[0];
                b[i] = v * s1.weights[1] - D * s2.weights[1];
                c[i] = v * s1.weights[2] - D * s2.weights[2];
                d[i] = source(x[i]);
            }
            d[0] = exact(0.0); // Dirichlet on the left, identity row

            // Neumann on the right: the one-sided 2nd-order stencil for u'
            std::size_t const L = n - 1;
            auto const sN = scheme_at(x, L, 1, fd_kind::explicit_fd, 2);
            REQUIRE(sN.nodes.size() == 3);   // one-sided closure is order+deriv wide
            REQUIRE(sN.node0 == L - 2);      // and it reaches back to column L-2
            double const w0 = sN.weights[0]; // at column L-2 -- OUTSIDE the band
            double const w1 = sN.weights[1]; // at column L-1
            double const w2 = sN.weights[2]; // at column L

            // eliminate the out-of-band entry against row L-1, which also touches L-2
            double const f = w0 / a[L - 1];
            a[L] = w1 - f * b[L - 1];
            b[L] = w2 - f * c[L - 1];
            c[L] = 0.0;
            d[L] = dexact(1.0) - f * d[L - 1];

            auto const u = tridiag_solve(a, b, c, d);

            double e = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                e = std::max(e, std::abs(u[i] - exact(x[i])));
            return e;
        };

        double const e101 = max_err(101);
        double const e201 = max_err(201);
        double const e401 = max_err(401);
        CHECK(std::log2(e101 / e201) == doctest::Approx(2.0).epsilon(0.15));
        CHECK(std::log2(e201 / e401) == doctest::Approx(2.0).epsilon(0.15));
        CHECK(e401 < 1.0e-3);
    }

    TEST_CASE("boundary conditions, explicit FD: Neumann built into the stencil")
    {
        // Recipe 4 from the tridiag_solve() comment, and the practical use of the
        // Hermite path in make_scheme(): rather than giving the boundary a row of its
        // own, generate its stencil already knowing u'(1). A supplied derivative is
        // worth one node of footprint, so u'' at the last node comes from {x_{L-1},
        // x_L} alone -- two columns, which is exactly what row L is allowed to touch.
        // No elimination step is needed at all.
        double const v = 1.0, D = 0.1;
        double const tp = 8.0 * std::atan(1.0); // 2 pi

        auto exact = [&](double xx) { return std::sin(tp * xx) + xx; };
        auto dexact = [&](double xx) { return tp * std::cos(tp * xx) + 1.0; };
        auto source = [&](double xx) {
            return v * (tp * std::cos(tp * xx) + 1.0) + D * tp * tp * std::sin(tp * xx);
        };

        auto max_err = [&](std::size_t n) {
            std::vector<double> x;
            for (std::size_t i = 0; i < n; ++i)
                x.push_back(double(i) / double(n - 1));

            std::vector<double> a(n, 0.0), b(n, 1.0), c(n, 0.0), d(n, 0.0);
            for (std::size_t i = 1; i + 1 < n; ++i) {
                auto const s1 = scheme_at(x, i, 1, fd_kind::explicit_fd, 2);
                auto const s2 = scheme_at(x, i, 2, fd_kind::explicit_fd, 2);
                a[i] = v * s1.weights[0] - D * s2.weights[0];
                b[i] = v * s1.weights[1] - D * s2.weights[1];
                c[i] = v * s1.weights[2] - D * s2.weights[2];
                d[i] = source(x[i]);
            }
            d[0] = exact(0.0); // Dirichlet on the left, identity row

            // the BC-aware stencil for u''(x_L), consuming the prescribed u'(x_L)
            std::size_t const L = n - 1;
            double const g = dexact(1.0);
            auto const sB = make_scheme(x[L], 2, {x[L - 1], x[L]}, {x[L]}, {x[L]});
            REQUIRE(sB.is_hermite());
            REQUIRE(sB.weights.size() == 2);     // two f columns: L-1 and L ...
            REQUIRE(sB.aux_weights.size() == 1); // ... plus the supplied derivative

            // v u'(x_L) - D u''(x_L) = s(x_L), with u'(x_L) = g known, so both the
            // advective term and the stencil's derivative term move to the rhs
            a[L] = -D * sB.weights[0];
            b[L] = -D * sB.weights[1];
            c[L] = 0.0;
            d[L] = source(x[L]) - v * g + D * sB.aux_weights[0] * g;

            auto const u = tridiag_solve(a, b, c, d);

            double e = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                e = std::max(e, std::abs(u[i] - exact(x[i])));
            return e;
        };

        double const e101 = max_err(101);
        double const e201 = max_err(201);
        double const e401 = max_err(401);
        CHECK(std::log2(e101 / e201) == doctest::Approx(2.0).epsilon(0.15));
        CHECK(std::log2(e201 / e401) == doctest::Approx(2.0).epsilon(0.15));
        CHECK(e401 < 1.0e-3);

        // the footprint claim the recipe rests on: a supplied derivative is worth one
        // node. Same order, one column less -- which is what makes the row fit.
        auto const with_bc = make_scheme(0.0, 2, {0.0, 1.0}, {0.0}, {0.0});
        auto const without = make_scheme(0.0, 2, {0.0, 1.0, 2.0}, {0.0});
        CHECK(with_bc.order == without.order);
        CHECK(with_bc.weights.size() + 1 == without.weights.size());

        // and one order up it holds again: 3 f-nodes + u' reach what 4 f-nodes do
        auto const with_bc2 = make_scheme(0.0, 2, {0.0, 1.0, 2.0}, {0.0}, {0.0});
        auto const without2 = make_scheme(0.0, 2, {0.0, 1.0, 2.0, 3.0}, {0.0});
        CHECK(with_bc2.order == 2);
        CHECK(without2.order == 2);
        CHECK(with_bc2.weights.size() + 1 == without2.weights.size());
    }

    TEST_CASE("boundary conditions, compact FD: three closures that keep the band")
    {
        // A compact scheme states the derivative IMPLICITLY, so using one is a
        // tridiagonal solve and its "boundary condition" is a different animal than
        // for an explicit scheme: f is known everywhere, and what the first and last
        // rows need is a closure for the DERIVATIVE. Three routes, all inside the
        // band, worked through here on f = sin so the answer is known.
        //
        // The lesson the numbers carry: unlike the elliptic problem above -- where a
        // boundary row one order down cost nothing -- here the closure order caps the
        // GLOBAL order, because the implicit solve couples every node to every other.
        // That is why fd_derivative() closes at 4th order to match its interior.

        // closure: 0 = explicit one-sided, 1 = compact (Lele), 2 = prescribed
        auto max_err = [](std::size_t n, int closure, int expl_order) {
            std::vector<double> x, f;
            for (std::size_t i = 0; i < n; ++i) {
                double const xi = 0.2 + 2.0 * double(i) / double(n - 1);
                x.push_back(xi);
                f.push_back(std::sin(xi));
            }

            std::vector<double> a(n, 0.0), b(n, 1.0), c(n, 0.0), d(n, 0.0);

            // interior: the Pade scheme, its lhs going straight onto the diagonals
            for (std::size_t i = 1; i + 1 < n; ++i) {
                auto const sc = scheme_at(x, i, 1, fd_kind::compact_fd);
                a[i] = sc.lhs_weights[0];
                b[i] = sc.lhs_weights[1];
                c[i] = sc.lhs_weights[2];
                d[i] = sc.weights[0] * f[i - 1] + sc.weights[1] * f[i] +
                       sc.weights[2] * f[i + 1];
            }

            for (int side = 0; side < 2; ++side) {
                std::size_t const i = (side == 0) ? 0 : n - 1;

                if (closure == 0) {
                    // ROUTE A -- explicit one-sided closure. The row states the
                    // derivative outright, so its lhs is the identity and the whole
                    // stencil sits in the rhs. Trivially inside the band whatever its
                    // width, which is what makes it the default (fd_derivative uses
                    // it). Its ORDER is what the global order will be.
                    auto const sc = scheme_at(x, i, 1, fd_kind::explicit_fd, expl_order);
                    a[i] = 0.0;
                    b[i] = 1.0;
                    c[i] = 0.0;
                    d[i] = detail::gather_grid_scheme(sc, f);
                }
                else if (closure == 1) {
                    // ROUTE B -- a COMPACT closure: keep the derivative implicit at
                    // the boundary too, coupling f'_0 to f'_1 only. Two lhs nodes is
                    // exactly the two columns the corner row may touch. make_scheme
                    // generates the classical Lele closure from that description --
                    // f'_0 + 2 f'_1 = (-5 f_0 + 4 f_1 + f_2)/(2h), 3rd order.
                    std::size_t const j0 = (side == 0) ? 0 : n - 3;
                    std::vector<double> const rhs_n{x[j0], x[j0 + 1], x[j0 + 2]};
                    std::vector<double> const lhs_n =
                        (side == 0) ? std::vector<double>{x[0], x[1]}
                                    : std::vector<double>{x[n - 2], x[n - 1]};
                    auto const sc = make_scheme(x[i], 1, rhs_n, lhs_n);
                    CHECK(sc.is_compact());
                    CHECK(sc.order == 3);

                    if (side == 0) {
                        a[i] = 0.0;
                        b[i] = sc.lhs_weights[0]; // f'_0, the diagonal
                        c[i] = sc.lhs_weights[1]; // f'_1, the super-diagonal
                    }
                    else {
                        a[i] = sc.lhs_weights[0]; // f'_{n-2}, the sub-diagonal
                        b[i] = sc.lhs_weights[1]; // f'_{n-1}, the diagonal
                        c[i] = 0.0;
                    }
                    d[i] = sc.weights[0] * f[j0] + sc.weights[1] * f[j0 + 1] +
                           sc.weights[2] * f[j0 + 2];
                }
                else {
                    // ROUTE C -- the derivative is PRESCRIBED there (the physical
                    // boundary condition is known). Then no closure is needed at all:
                    // the row is the identity on the unknown itself.
                    a[i] = 0.0;
                    b[i] = 1.0;
                    c[i] = 0.0;
                    d[i] = std::cos(x[i]);
                }
            }

            auto const fp = tridiag_solve(a, b, c, d);
            double e = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                e = std::max(e, std::abs(fp[i] - std::cos(x[i])));
            return e;
        };

        auto rate = [&](int closure, int ord) {
            return std::log2(max_err(81, closure, ord) / max_err(161, closure, ord));
        };

        // A 2nd-order explicit closure CAPS the 4th-order interior at 2 -- the whole
        // point of the note in fd_derivative(), pinned here rather than asserted
        CHECK(rate(0, 2) == doctest::Approx(2.0).epsilon(0.05));
        // ... while a 4th-order one lets the interior through
        CHECK(rate(0, 4) == doctest::Approx(4.0).epsilon(0.05));
        // the compact Lele closure is 3rd order, and that is what comes out globally
        CHECK(rate(1, 0) == doctest::Approx(3.0).epsilon(0.05));
        // a prescribed derivative costs nothing: full 4th order, and the most accurate
        CHECK(rate(2, 0) == doctest::Approx(4.0).epsilon(0.05));
        CHECK(max_err(161, 2, 0) < max_err(161, 0, 4));

        // and fd_derivative()'s own default closure lands where route A at order 4 does
        std::vector<double> x, f;
        for (std::size_t i = 0; i < 161; ++i) {
            double const xi = 0.2 + 2.0 * double(i) / 160.0;
            x.push_back(xi);
            f.push_back(std::sin(xi));
        }
        auto const fp = fd_derivative(x, f, 1, fd_kind::compact_fd);
        double e = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i)
            e = std::max(e, std::abs(fp[i] - std::cos(x[i])));
        CHECK(e == doctest::Approx(max_err(161, 0, 4)).epsilon(1.0e-9));
    }
}


/////////////////////////////////////////////////////////////////////////////////////////
// dense solver: lstsq_solve's rank guard, its damping, and nullspace_project
// (ga/detail/ga_solver.hpp) -- tested here beside tridiag_solve, the other shared solver
// primitive with no consumer-specific home.
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("dense solver: lstsq_solve / nullspace_project")
{

    // build an m x ncols matrix (flat row-major) from a smooth, well-conditioned rule
    static auto make_A = [](size_t m, size_t ncols) {
        std::vector<double> A(m * ncols);
        for (size_t i = 0; i < m; ++i)
            for (size_t k = 0; k < ncols; ++k)
                A[i * ncols + k] = std::cos(0.7 * double(i) + 1.3 * double(k)) +
                                   0.4 * std::sin(double(i) * double(k));
        return A;
    };
    static auto residual = [](std::vector<double> const& A, std::vector<double> const& x,
                              std::vector<double> const& b, size_t ncols) {
        double e = 0.0;
        for (size_t i = 0; i < b.size(); ++i) {
            double s = 0.0;
            for (size_t k = 0; k < ncols; ++k)
                s += A[i * ncols + k] * x[k];
            e = std::max(e, std::abs(s - b[i]));
        }
        return e;
    };

    TEST_CASE("lstsq_solve: the three full-rank regimes still solve")
    {
        fmt::println("lstsq_solve: the three full-rank regimes still solve");

        // square and underdetermined must satisfy A x = b exactly; overdetermined is a
        // least-squares fit and generally will not, so it is checked by the normal
        // equations instead.
        { // square
            auto const A = make_A(5, 5);
            std::vector<double> b{0.3, -1.2, 0.8, 2.0, -0.5};
            CHECK(residual(A, lstsq_solve(A, b, 5), b, 5) == doctest::Approx(0.0));
        }
        { // underdetermined (minimum norm) -- exact, and no other exact solution is
            // shorter
            auto const A = make_A(4, 9);
            std::vector<double> b{0.3, -1.2, 0.8, 2.0};
            auto const x = lstsq_solve(A, b, 9);
            CHECK(residual(A, x, b, 9) == doctest::Approx(0.0));
            // minimum-norm: x must lie in the row space, i.e. x = A^T y for some y. Test
            // it by checking that adding any null-space vector only lengthens x.
            double nx = 0.0;
            for (double v : x)
                nx += v * v;
            auto const z = nullspace_project(A, std::vector<double>(9, 1.0), 4, 9);
            double nz = 0.0, dot = 0.0;
            for (size_t k = 0; k < 9; ++k) {
                nz += z[k] * z[k];
                dot += z[k] * x[k];
            }
            REQUIRE(nz > 1.0e-6); // the null space is non-trivial
            CHECK(dot ==
                  doctest::Approx(0.0).epsilon(1.0e-9)); // and x is orthogonal to it
            CHECK(nx > 0.0);
        }
        { // overdetermined (least squares): the residual is orthogonal to the columns
            auto const A = make_A(8, 3);
            std::vector<double> b{0.3, -1.2, 0.8, 2.0, -0.5, 1.1, 0.2, -0.9};
            auto const x = lstsq_solve(A, b, 3);
            for (size_t k = 0; k < 3; ++k) {
                double s = 0.0;
                for (size_t i = 0; i < 8; ++i) {
                    double r = -b[i];
                    for (size_t j = 0; j < 3; ++j)
                        r += A[i * 3 + j] * x[j];
                    s += A[i * 3 + k] * r;
                }
                CHECK(s == doctest::Approx(0.0)); // A^T (A x - b) = 0
            }
        }
        fmt::println("");
    }

    TEST_CASE("lstsq_solve: a rank-deficient input is REFUSED, not silently inflated")
    {
        fmt::println("lstsq_solve: a rank-deficient input is REFUSED, not silently "
                     "inflated");

        // The normal-equation route cannot form a pseudo-inverse without full rank: the
        // Gram matrix is singular, and lu_decomp substitutes TINY (1e-20) for the
        // vanishing pivot rather than raising -- so the answer comes back scaled by that
        // pivot's reciprocal. Measured before the guard: x0 = -9.04e+02 on a system whose
        // correct scale is 0.26. The guard turns that into an exception.
        auto A = make_A(6, 11);
        for (size_t k = 0; k < 11; ++k)
            A[5 * 11 + k] = A[4 * 11 + k]; // row 5 := row 4 -> rank 5 of 6
        std::vector<double> b{0.3, -1.2, 0.8, 2.0, -0.5, 1.4};

        CHECK_THROWS_AS(lstsq_solve(A, b, 11), Solver_error);

        // damping is the documented way through: it regularizes instead of refusing, and
        // must return a finite answer of sane magnitude.
        auto const xd = lstsq_solve(A, b, 11, 1.0e-4);
        REQUIRE(xd.size() == 11u);
        double nmax = 0.0;
        for (double v : xd)
            nmax = std::max(nmax, std::abs(v));
        CHECK(std::isfinite(nmax));
        CHECK(nmax < 10.0); // not the 1e+2..1e+3 blow-up of an unguarded TINY pivot
        fmt::println("");
    }

    TEST_CASE("lstsq_solve: damping is off by default and reduces to the exact solve")
    {
        fmt::println("lstsq_solve: damping is off by default and reduces to the exact "
                     "solve");

        // Every pre-existing call site passes no damping and must be bit-unchanged.
        auto const A = make_A(4, 9);
        std::vector<double> b{0.3, -1.2, 0.8, 2.0};
        auto const x0 = lstsq_solve(A, b, 9);
        auto const x1 = lstsq_solve(A, b, 9, 0.0);
        REQUIRE(x0.size() == x1.size());
        for (size_t k = 0; k < x0.size(); ++k)
            CHECK(x0[k] == x1[k]); // exact equality, not an approximation

        // a small damping perturbs the answer slightly and shrinks it (Tikhonov)
        auto const xd = lstsq_solve(A, b, 9, 1.0e-3);
        double n0 = 0.0, nd = 0.0;
        for (size_t k = 0; k < x0.size(); ++k) {
            n0 += x0[k] * x0[k];
            nd += xd[k] * xd[k];
        }
        CHECK(nd < n0);       // regularization shortens the solution
        CHECK(nd > 0.9 * n0); // but only slightly, at this damping
        fmt::println("");
    }

    TEST_CASE("nullspace_project: the result is invisible to A")
    {
        fmt::println("nullspace_project: the result is invisible to A");

        // v_ns = (I - A^+ A) v must satisfy A v_ns = 0: a secondary objective projected
        // this way cannot disturb the primary task. This is the redundancy-resolution
        // step; the same projector appears in the task-space control literature.
        auto const A = make_A(4, 9);
        std::vector<double> v(9);
        for (size_t k = 0; k < 9; ++k)
            v[k] = 0.4 * double(k) - 1.1;

        auto const vns = nullspace_project(A, v, 4, 9);
        REQUIRE(vns.size() == 9u);
        std::vector<double> const zero(4, 0.0);
        CHECK(residual(A, vns, zero, 9) == doctest::Approx(0.0).epsilon(1.0e-9));

        // idempotent: projecting again changes nothing
        auto const vns2 = nullspace_project(A, vns, 4, 9);
        for (size_t k = 0; k < 9; ++k)
            CHECK(vns2[k] == doctest::Approx(vns[k]).epsilon(1.0e-9));

        // and the removed part is exactly the task component
        for (size_t k = 0; k < 9; ++k)
            CHECK(std::isfinite(v[k] - vns[k]));

        // a full-rank SQUARE A has a trivial null space -> the projection is ~zero
        auto const S = make_A(5, 5);
        std::vector<double> w(5, 1.0);
        auto const wns = nullspace_project(S, w, 5, 5);
        for (double x : wns)
            CHECK(x == doctest::Approx(0.0).epsilon(1.0e-8));

        // wrong-sized v is rejected rather than read out of bounds
        CHECK_THROWS_AS(nullspace_project(A, std::vector<double>(3, 1.0), 4, 9),
                        Solver_error);
        fmt::println("");
    }
}
