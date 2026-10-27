// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

/////////////////////////////////////////////////////////////////////////////////////////
// The preconditions of the dense solver and the finite-difference builders: each entry
// point refuses a malformed input with its own message instead of returning a silently
// wrong result. The message is pinned, so a different error thrown for another reason
// cannot pass for the one under test.
/////////////////////////////////////////////////////////////////////////////////////////

#include <array>
#include <mdspan>

TEST_SUITE("numerics: preconditions")
{
    TEST_CASE("dense solver: malformed shapes are refused by every entry point")
    {
        // lu_decomp / lu_backsubs: the matrix must be square and the permutation and the
        // right-hand side must match it
        std::array<double, 6> a23{1, 2, 3, 4, 5, 6};
        std::array<int, 2> perm{};
        std::mdspan<double, std::dextents<size_t, 2>> m23(a23.data(), 2, 3);
        std::mdspan<int, std::dextents<size_t, 1>> p2(perm.data(), 2);
        CHECK_THROWS_WITH_AS(lu_decomp(m23, p2), doctest::Contains("non-square"),
                             Solver_error);

        std::array<double, 4> a22{4, 1, 1, 3};
        std::array<double, 3> b3{1, 2, 3};
        std::mdspan<double const, std::dextents<size_t, 2>> c22(a22.data(), 2, 2);
        std::mdspan<int const, std::dextents<size_t, 1>> cp2(perm.data(), 2);
        std::mdspan<double, std::dextents<size_t, 1>> mb3(b3.data(), 3);
        CHECK_THROWS_WITH_AS(lu_backsubs(c22, cp2, mb3),
                             doctest::Contains("right-hand side size"), Solver_error);

        // det: only a square matrix has one
        std::mdspan<double, std::dextents<size_t, 2>> d23(a23.data(), 2, 3);
        CHECK_THROWS_AS((void)det(d23), std::invalid_argument);

        // the vector interfaces: A must hold rows * ncols entries
        std::vector<double> const A5(5, 1.0), b2{1.0, 2.0};
        CHECK_THROWS_WITH_AS((void)minnorm_solve(A5, b2, 3),
                             doctest::Contains("b.size() * ncols"), Solver_error);
        CHECK_THROWS_WITH_AS((void)detail::qr_decomp(A5, 2, 3),
                             doctest::Contains("does not match m x n"), Solver_error);

        std::vector<double> const lo3(3, -1.0), hi3(3, 1.0), lo2(2, -1.0);
        CHECK_THROWS_WITH_AS((void)bvls_solve(A5, b2, 3, lo3, hi3),
                             doctest::Contains("b.size() * ncols"), Solver_error);
        std::vector<double> const A6(6, 1.0);
        CHECK_THROWS_WITH_AS((void)bvls_solve(A6, b2, 3, lo2, hi3),
                             doctest::Contains("lo and hi must have ncols"),
                             Solver_error);

        std::vector<double> const none{};
        std::vector<double> x3(3, 0.0), x2(2, 0.0);
        CHECK_THROWS_WITH_AS(
            (void)qp_ls_solve(A5, b2, 3, none, none, none, none, lo3, hi3, x3),
            doctest::Contains("rows * ncols"), Solver_error);
        CHECK_THROWS_WITH_AS(
            (void)qp_ls_solve(A6, b2, 3, none, none, none, none, lo3, hi3, x2),
            doctest::Contains("must have ncols entries"), Solver_error);
    }

    TEST_CASE("dense solver: a zero row is singular, and LU says so")
    {
        // a row of zeros cannot be pivoted: lu_solve must refuse, not divide by zero
        std::vector<double> const A{1.0, 2.0, 0.0, 0.0};
        std::vector<double> const b{1.0, 1.0};
        CHECK_THROWS_WITH_AS((void)lu_solve(A, b, 2), doctest::Contains("singular"),
                             Solver_error);
    }

    TEST_CASE("dense solver: kkt multipliers need a factorization made for them")
    {
        // M = I (2x2), one constraint row G = [1 1]: x1 + x2 = g
        std::vector<double> const M{1.0, 0.0, 0.0, 1.0};
        std::vector<double> const G{1.0, 1.0};
        std::vector<double> const f{0.0, 0.0}, g{2.0};

        auto const with = kkt_factorize(M, G, 2, 1);
        std::vector<double> lambda;
        auto const x = with.solve(f, g, &lambda);
        CHECK(x[0] == doctest::Approx(1.0));
        CHECK(x[1] == doctest::Approx(1.0));
        REQUIRE(lambda.size() == 1);

        // made without multipliers it still solves for x, and refuses lambda
        auto const without = kkt_factorize(M, G, 2, 1, false);
        auto const x_only = without.solve(f, g);
        CHECK(x_only[0] == doctest::Approx(1.0));
        CHECK(x_only[1] == doctest::Approx(1.0));
        CHECK_THROWS_WITH_AS((void)without.solve(f, g, &lambda),
                             doctest::Contains("multipliers were asked for"),
                             Solver_error);
    }

    TEST_CASE("dense solver: damped lstsq_solve is the Tikhonov solution")
    {
        // lambda^2 = damping * trace(Gram) / n is added to the Gram matrix's diagonal
        double const damping = 0.5;

        // overdetermined, A^T A = I: x = A^T b / (1 + lambda^2), lambda^2 = damping
        std::vector<double> const Ao{1.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        std::vector<double> const bo{1.0, 2.0, 3.0};
        auto const xo = lstsq_solve(Ao, bo, 2, damping);
        REQUIRE(xo.size() == 2);
        CHECK(xo[0] == doctest::Approx(1.0 / 1.5));
        CHECK(xo[1] == doctest::Approx(2.0 / 1.5));

        // underdetermined, A A^T = I: x = A^T b / (1 + lambda^2)
        std::vector<double> const Au{1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
        std::vector<double> const bu{1.0, 2.0};
        auto const xu = lstsq_solve(Au, bu, 3, damping);
        REQUIRE(xu.size() == 3);
        CHECK(xu[0] == doctest::Approx(1.0 / 1.5));
        CHECK(xu[1] == doctest::Approx(2.0 / 1.5));
        CHECK(xu[2] == doctest::Approx(0.0));

        // rank deficient: the damping IS the answer to the rank drop -- the two equal
        // columns share the load equally, (3 + 3d) t + 3 t = 3 with lambda^2 = 3 d
        std::vector<double> const Ar{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
        std::vector<double> const br{1.0, 1.0, 1.0};
        auto const xr = lstsq_solve(Ar, br, 2, 1.0);
        CHECK(xr[0] == doctest::Approx(1.0 / 3.0));
        CHECK(xr[1] == doctest::Approx(1.0 / 3.0));
    }

    TEST_CASE("fd schemes: the builders refuse what they cannot construct")
    {
        CHECK_THROWS_WITH_AS((void)make_scheme(0.0, 3, {-1.0, 0.0, 1.0}, {0.0}),
                             doctest::Contains("deriv must be 1 or 2"),
                             std::invalid_argument);
        CHECK_THROWS_WITH_AS((void)make_scheme(0.0, 1, {-1.0, 0.0, 1.0}, {}),
                             doctest::Contains("at least one lhs node"),
                             std::invalid_argument);

        // a centred scheme only reaches even orders
        CHECK_THROWS_AS((void)central_scheme(1, 3), std::invalid_argument);
        CHECK_THROWS_AS((void)central_scheme(2, 0), std::invalid_argument);
        CHECK_THROWS_AS((void)staggered_scheme(1, 1), std::invalid_argument);
        CHECK_THROWS_AS((void)staggered_scheme(1, 5), std::invalid_argument);

        // a compact scheme couples its neighbours' derivative values, so it needs two
        std::vector<double> const x{0.0, 0.1, 0.25, 0.45, 0.7};
        CHECK_THROWS_WITH_AS((void)scheme_at(x, 0, 1, fd_kind::compact_fd),
                             doctest::Contains("interior node"), std::invalid_argument);
        CHECK_THROWS_WITH_AS((void)scheme_at(x, x.size() - 1, 1, fd_kind::compact_fd),
                             doctest::Contains("interior node"), std::invalid_argument);
        CHECK_NOTHROW((void)scheme_at(x, 2, 1, fd_kind::compact_fd));
    }

    TEST_CASE("dense solver: a damped lstsq_solve of a SQUARE system is damped too")
    {
        // a square A takes the normal-equation route like an overdetermined one,
        // (A^T A + l^2 I) x = A^T b with l^2 = damping * trace(A^T A) / n. Until
        // 2026-10-27 this case pinned the opposite -- the damping ignored, the plain
        // solution returned -- which is the behaviour that was changed
        std::vector<double> const A{4.0, 1.0, 1.0, 3.0};
        std::vector<double> const b{1.0, 2.0};
        auto const x = lstsq_solve(A, b, 2, 0.5);
        double const l2 = 0.5 * (17.0 + 10.0) / 2.0; // trace(A^T A) = 17 + 10
        std::vector<double> const G{17.0 + l2, 7.0, 7.0, 10.0 + l2};
        auto const ref = lu_solve(G, std::vector<double>{6.0, 7.0}, 2); // A^T b
        CHECK(x[0] == doctest::Approx(ref[0]));
        CHECK(x[1] == doctest::Approx(ref[1]));
        auto const plain = lu_solve(A, b, 2);
        CHECK(std::abs(x[0] - plain[0]) > 0.01); // and it is not the plain solution
    }
}
