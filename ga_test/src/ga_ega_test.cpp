// Copyright 2024-2025, Daniel Hug. All rights reserved.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <chrono>
#include <iostream>
#include <tuple>
#include <vector>

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_ega.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)

TEST_SUITE("Euclidean Geometric Algebra (EGA)")
{

    TEST_CASE("G<2,0,0>: ega2d")
    {
        fmt::println("G<2,0,0>: ega2d");
        // 2d euklidean geometric algebra
        const algebra<2> alg;
        CHECK(alg.p() == 2);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 2);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 4);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 3); // == dim_space + 1
        fmt::println("   ega2d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   ega2d: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }


    ////////////////////////////////////////////////////////////////////////////////
    // Vec2d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2d: default init")
    {
        fmt::println("Vec2d: default init");
        vec2d v;
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.x) < eps);
        CHECK(std::abs(v.y) < eps);
    }
    TEST_CASE("Vec2d: with curly braced intializer")
    {
        fmt::println("Vec2d: with curly braced intializer");
        vec2d v{0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.x) < eps);
        CHECK(std::abs(v.y) < eps);
    }
    TEST_CASE("Vec2d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("Vec2d: cp ctor & cp assign incl. type deduction");
        vec2d v1{1.0, 2.0}; // init with double (type deduction)
        vec2d v2{v1};       // cp ctor
        vec2d v3 = v2;      // cp assign
        vec2d v4 = -v2;     // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(std::abs(v1.x - 1.0) < eps);
        CHECK(std::abs(v1.y - 2.0) < eps);
        CHECK(std::abs(v2.x - 1.0) < eps);
        CHECK(std::abs(v2.y - 2.0) < eps);
        CHECK(std::abs(v3.x - 1.0) < eps);
        CHECK(std::abs(v3.y - 2.0) < eps);
        CHECK(v4 == -v2);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }


    TEST_CASE("Vec2d: fmt & cout printing")
    {
        fmt::println("Vec2d: fmt & cout printing:\n");

        vec2d pf{1.0f, 2.00001f};
        vec2d pd{1.0, 2.00001};

        std::cout << "       cout: pf = " << pf << std::endl;
        fmt::println("       fmt:  pf = {}", pf);
        fmt::println("       fmt:  pf = {:.4f}", pf);
        fmt::println("");

        std::cout << "       cout: pd = " << pd << std::endl;
        fmt::println("       fmt:  pd = {}", pd);
        fmt::println("       fmt:  pd = {:.4f}", pd);
        fmt::println("");

        std::vector<Vec2d<double>> vp1{{1.0, 1.0}, {1.5, 2.0}};
        fmt::println("       fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("       fmt: vp1 = {:g}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("Vec2d: comparison float")
    {
        fmt::println("Vec2d: comparison float");

        Vec2d<float> v1f{1.0f, 2.0f};
        Vec2d<float> v2f{2.0f, 4.0f};
        Vec2d<float> v3f{1.0f, 2.0000001f};
        Vec2d<float> v4f{v1f};

        // fmt::println("   v1f = {}", v1f);
        // fmt::println("   v2f = {}", v2f);
        // fmt::println("   v3f = {}", v3f);
        // fmt::println("   v4f = {}", v4f);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<float>::epsilon());

        CHECK(v1f == v4f);           // comparison (equality)
        CHECK(v1f != v2f);           // comparison (inequality)
        CHECK(nrm(v1f) < nrm(v2f));  // comparison (less than)
        CHECK(nrm(v2f) >= nrm(v1f)); // comparison (greater than or equal)
        CHECK(v3f == v1f);           // comparison (eqality)
    }

    TEST_CASE("Vec2d: comparison double")
    {
        fmt::println("Vec2d: comparison double");

        Vec2d<double> v1d{1.0, 2.0};
        Vec2d<double> v2d{2.0, 4.0};
        Vec2d<double> v3d{1.0, 2.0000000000000001};
        Vec2d<double> v4d{v1d};

        // fmt::println("   v1d = {}", v1d);
        // fmt::println("   v2d = {}", v2d);
        // fmt::println("   v3d = {}", v3d);
        // fmt::println("   v4d = {}", v4d);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<double>::epsilon());

        CHECK(v1d == v4d);           // comparison (equality)
        CHECK(v1d != v2d);           // comparison (inequality)
        CHECK(nrm(v1d) < nrm(v2d));  // comparison norm
        CHECK(nrm(v2d) >= nrm(v1d)); // comparison norm
        CHECK(v3d == v1d);           // comparison (eqality)
    }

    TEST_CASE("Vec2d: vector space and linearity tests")
    {
        fmt::println("Vec2d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        vec2d p0;
        vec2d p1{1.0, 2.0};
        vec2d p2{2.0, 4.0};
        vec2d p3{3.0, 6.0};
        vec2d p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t. multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t. addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
        CHECK(p4 == -p1);
    }

    TEST_CASE("Vec2d: inner product properties")
    {
        fmt::println("Vec2d: inner product properties");

        double a = 2.35;
        vec2d u{1.0, 2.0};
        vec2d v{-0.5, 3.0};
        vec2d w{3.0, 6.0};

        CHECK(dot(a * u, v) == a * dot(u, v));
        CHECK(dot(u + v, w) == dot(u, w) + dot(v, w));
        CHECK(dot(u, v) == dot(v, u));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec2d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2d: operations - norm, inverse, dot")
    {
        fmt::println("Vec2d: operations - norm, inverse, dot");

        scalar2d s1{3.2};

        vec2d v1{2.0, 1.0};
        vec2d v2{normalize(v1)};

        vec2d v3{2.0, 6.0};
        vec2d v4{inv(v3)};

        pscalar2d ps1{-4.7};
        mvec2d_e mve1{s1, ps1};
        mvec2d mv1{s1, v1, ps1};

        // fmt::println("v1 = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2 = normalize(v1) = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));

        CHECK(std::abs(nrm_sq(v1) - 5.0) < eps);
        CHECK(std::abs(nrm_sq(v2) - 1.0) < eps);
        CHECK(std::abs(dot(v4, v3) - 1.0) < eps);

        auto m = vec2d{13.0, 5.0};
        auto prd = m * inv(m);
        CHECK(std::abs(gr0(prd) - 1.0) < eps);
        CHECK(std::abs(gr2(prd) - 0.0) < eps);

        // check inverses - scalar
        // fmt::println("");
        // fmt::println("s1 * inv(s1) = {}", s1 * inv(s1)); // s
        CHECK(std::abs(nrm(s1 * inv(s1)) - 1) < eps);
        CHECK(std::abs(inv(s1) - rev(s1) / nrm_sq(s1)) < eps);

        // check inverses - vector
        // fmt::println("v1 * inv(v1) = {}", v1 * inv(v1)); // mv_e
        CHECK(std::abs(nrm(gr0(v1 * inv(v1))) - 1) < eps);
        CHECK(std::abs(nrm(gr2(v1 * inv(v1))) - 0) < eps);
        CHECK(std::abs(nrm(inv(v1) - rev(v1) / nrm_sq(v1))) < eps);

        // check inverses - pseudoscalar
        // fmt::println("ps1 * inv(ps1) = {}", ps1 * inv(ps1)); // s
        CHECK(std::abs(nrm(ps1 * inv(ps1)) - 1) < eps);
        CHECK(std::abs(nrm(inv(ps1) - rev(ps1) / nrm_sq(ps1))) < eps);

        // check inverses - even grade multivector
        // fmt::println("mve1 * inv(mve1) = {}", mve1 * inv(mve1)); // mv_e
        CHECK(std::abs(nrm(gr0(mve1 * inv(mve1))) - 1) < eps);
        CHECK(std::abs(nrm(gr2(mve1 * inv(mve1))) - 0) < eps);
        CHECK(std::abs(nrm(inv(mve1) - rev(mve1) / nrm_sq(mve1))) < eps);

        // check inverses - multivector
        // fmt::println("mv1 * inv(mv1) = {}", mv1 * inv(mv1)); // mv
        CHECK(std::abs(nrm(gr0(mv1 * inv(mv1))) - 1) < eps);
        CHECK(std::abs(nrm(gr1(mv1 * inv(mv1))) - 0) < eps);
        CHECK(std::abs(nrm(gr2(mv1 * inv(mv1))) - 0) < eps);
        CHECK(std::abs(nrm(gr0(inv(mv1) * mv1)) - 1) < eps); // left and right inverse
                                                             // are equal
        // fmt::println("");
    }

    TEST_CASE("Vec2d: operations - angle")
    {
        fmt::println("Vec2d: operations - angle");

        std::vector<std::tuple<double, Vec2d<double>>> v1;
        std::vector<std::tuple<double, Vec2d<double>>> v2;
        std::vector<std::tuple<double, Vec2d<double>>> v3;

        for (int i = -12; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2d<double>(std::cos(phi), std::sin(phi));
            v1.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_2d, c));
        }
        // fmt::println("");

        for (int i = -12; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2d<double>(std::cos(phi + pi / 2), std::sin(phi + pi / 2));
            v2.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e2_2d, c));
        }
        // fmt::println("");

        for (int i = -12; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2d<double>(std::cos(phi + pi / 4), std::sin(phi + pi / 4));
            v3.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_2d + e2_2d, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v1) {
            CHECK(std::abs(phi - angle(e1_2d, c)) < eps);
        }
        for (auto const& [phi, c] : v2) {
            CHECK(std::abs(phi - angle(e2_2d, c)) < eps);
        }
        auto ref_vec = normalize(e1_2d + e2_2d);
        for (auto const& [phi, c] : v3) {
            CHECK(std::abs(phi - angle(ref_vec, c)) < eps);
        }

        // auto v = Vec2d<double>{1.0, 0.0};
        // // auto v = Vec2d<double>{1.0, 1.0};
        // for (auto const& [phi, c] : v1) {
        //     auto u1 = v * c;
        //     auto u2 = c * v;
        //     fmt::println("   phi={: .4f}, phi={:> 4.0f}°, c={: .3f},"
        //                  "  u1={: .3f}, u2={: .3f}",
        //                  phi, phi * 180 / pi, c, u1, u2);
        // }
        // fmt::println("");
    }

    TEST_CASE("Vec2d: operations - wedge")
    {
        fmt::println("Vec2d: operations - wedge");

        std::vector<std::tuple<double, Vec2d<double>>> v;

        for (int i = -12; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2d<double>(std::cos(phi), std::sin(phi));
            v.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_2d, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v) {
            CHECK(std::abs(wdg(e1_2d, c) - sin(angle(e1_2d, c))) < eps);
        }
    }

    TEST_CASE("Vec2d: operations - project / reject / reflect")
    {
        fmt::println("Vec2d: operations - project / reject / reflect");

        vec2d v1{1.5, 4.0};
        vec2d v2{3.0, 1.5};

        vec2d v3{project_onto(v1, v2)};
        vec2d v4{reject_from(v1, v2)};
        vec2d v5{v3 + v4};

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        // fmt::println("v2u = {: .4f}, nrm(v2) = {: .4f}", v2u, nrm(v2u));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .4f},"
        //              " nrm(v3) = {: .4f}",
        //              v3, nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .4f},"
        //              " nrm(v4) = {: .4f}",
        //              v4, nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .4f},"
        //              " nrm(v5) = {: .4f}",
        //              v5, nrm(v5));
        // fmt::println("");

        vec2d v{1.0, 3.0};
        vec2d b{e1_2d + e2_2d};

        // fmt::println("v  = {}", v);
        // fmt::println("b  = {}", b);
        // fmt::println("reflect_on_vec(v,b)  = {}", reflect_on_vec(v, b));
        // fmt::println("reflect_on(v,b)  = {}", reflect_on(v, b));
        // fmt::println("");


        CHECK(v3 + v4 == v5);
        CHECK(std::abs(dot(v3, v4)) < eps);
        CHECK(v5 == v1);

        CHECK(inv(v2) * wdg(v2, v1) == wdg(v1, v2) * inv(v2));

        CHECK(reflect_on_vec(v, b).x == 3);
        CHECK(reflect_on_vec(v, b).y == 1);
        CHECK(reflect_on(v, b).x == -3);
        CHECK(reflect_on(v, b).y == -1);

        // checking time required
        //
        // auto start = std::chrono::system_clock::now();
        // for (size_t i = 0; i < 10000000; ++i) {
        //     vec2d v = reject_from(v1, v2);
        // }
        // auto end = std::chrono::system_clock::now();
        // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end -
        // start); fmt::println("The measurement took {}", elapsed);

        // point reflected on a vector
        vec2d p{4, 1};
        CHECK(reflect_on_vec(p, x_axis_2d) == vec2d{4, -1});

        // coordinate axis reflected on perpendicular axis yield their negatives
        CHECK(reflect_on_vec(y_axis_2d, x_axis_2d) == -y_axis_2d);
        CHECK(reflect_on_vec(x_axis_2d, y_axis_2d) == -x_axis_2d);

        // coordinate axis reflected on itself remains itself (identity)
        CHECK(reflect_on_vec(x_axis_2d, x_axis_2d) == x_axis_2d);
        CHECK(reflect_on_vec(y_axis_2d, y_axis_2d) == y_axis_2d);

        // point reflected on a hyperplane that the vector is a normal to
        // the hyperplane can be created by taking the dual (or the rcmpl) of the normal
        CHECK(reflect_on(p, dual(x_axis_2d)) == vec2d{4, -1});

        // coordinate axis reflected on perpendicular axis yield their negatives
        CHECK(reflect_on(y_axis_2d, dual(x_axis_2d)) == -y_axis_2d);
        CHECK(reflect_on(x_axis_2d, dual(y_axis_2d)) == -x_axis_2d);

        // coordinate axis reflected on itself remains itself (identity)
        CHECK(reflect_on(x_axis_2d, dual(x_axis_2d)) == x_axis_2d);
        CHECK(reflect_on(y_axis_2d, dual(y_axis_2d)) == y_axis_2d);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec2d<T> Gram-Schmidt-Orthogonalization
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2d: Gram-Schmidt-Orthogonalization")
    {
        fmt::println("Vec2d: Gram-Schmidt-Orthogonalization");

        vec2d v1{2.0, 2.0};
        vec2d v2{-1.5, 2.5};
        vec2d v3{1.5, -2.5};

        std::vector<vec2d> og12 = gs_orthogonal(v1, v2);
        std::vector<vec2d> og13 = gs_orthogonal(v1, v3);

        std::vector<vec2d> on12 = gs_orthonormal(v1, v2);
        std::vector<vec2d> on13 = gs_orthonormal(v1, v3);

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}, angle = {: .2f}°", v1, nrm(v1),
        //              rad2deg(angle(e1_2d, v1)));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v2, nrm(v2),
        //              rad2deg(angle(e1_2d, v2)));
        // fmt::println("v3  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v3, nrm(v3),
        //              rad2deg(angle(e1_2d, v3)));
        // fmt::println("");
        // fmt::println("og12[0]  = {: .4f}, nrm(og12[0]) = {: .4f}, angle = {: .2f}°",
        //              og12[0], nrm(og12[0]), rad2deg(angle(e1_2d, og12[0])));
        // fmt::println("og12[1]  = {: .4f}, nrm(og12[1]) = {: .4f}, angle = {: .2f}°",
        //              og12[1], nrm(og12[1]), rad2deg(angle(e1_2d, og12[1])));
        // fmt::println("on12[0]  = {: .4f}, nrm(on12[0]) = {: .4f}, angle = {: .2f}°",
        //              on12[0], nrm(on12[0]), rad2deg(angle(e1_2d, on12[0])));
        // fmt::println("on12[1]  = {: .4f}, nrm(on12[1]) = {: .4f}, angle = {: .2f}°",
        //              on12[1], nrm(on12[1]), rad2deg(angle(e1_2d, on12[1])));
        // fmt::println("");
        // fmt::println("og13[0]  = {: .4f}, nrm(og13[0]) = {: .4f}, angle = {: .2f}°",
        //              og13[0], nrm(og13[0]), rad2deg(angle(e1_2d, og13[0])));
        // fmt::println("og13[1]  = {: .4f}, nrm(og13[1]) = {: .4f}, angle = {: .2f}°",
        //              og13[1], nrm(og13[1]), rad2deg(angle(e1_2d, og13[1])));
        // fmt::println("on13[0]  = {: .4f}, nrm(on13[0]) = {: .4f}, angle = {: .2f}°",
        //              on13[0], nrm(on13[0]), rad2deg(angle(e1_2d, on13[0])));
        // fmt::println("on13[1]  = {: .4f}, nrm(on13[1]) = {: .4f}, angle = {: .2f}°",
        //              on13[1], nrm(on13[1]), rad2deg(angle(e1_2d, on13[1])));
        // fmt::println("");
        // fmt::println(
        //     "sign(wdg(on12[0],on12[1])/I_2d) = {} (+1: right-handed, -1: left-handed)",
        //     sign(wdg(on12[0], on12[1]) / I_2d));
        // fmt::println(
        //     "sign(wdg(on13[0],on13[1])/I_2d) = {} (+1: right-handed, -1: left-handed)",
        //     sign(wdg(on13[0], on13[1]) / I_2d));
        // fmt::println("");

        CHECK(std::abs(dot(og12[0], og12[1])) < eps);
        CHECK(std::abs(dot(on12[0], on12[1])) < eps);
        CHECK(std::abs(nrm(on12[0]) - 1.0) < eps);
        CHECK(std::abs(nrm(on12[1]) - 1.0) < eps);
        CHECK(std::abs(dot(og13[0], og13[1])) < eps);
        CHECK(std::abs(dot(on13[0], on13[1])) < eps);
        CHECK(std::abs(nrm(on13[0]) - 1.0) < eps);
        CHECK(std::abs(nrm(on13[1]) - 1.0) < eps);
    }

    TEST_CASE("Vec2d: vector ratio")
    {
        fmt::println("Vec2d: vector ratio");

        auto v1 = vec2d{1, 0};
        auto v2 = vec2d{2, 2};
        auto v3 = vec2d{0, 2};

        // x : v3 = v2 : v1 => x * inv(v3) = v2 * inv(v1)
        //                  => x = v2 * inv(v1) * v3
        //
        // should be same angle between v3 and x as between v1 and v2
        // should be same factor of norms between x and v3 as between v2 and v1

        auto Rs = v2 * inv(v1);
        auto x = Rs * v3;

        // fmt::println("");
        // fmt::println("v1 = {}", v1);
        // fmt::println("v2 = {}", v2);
        // fmt::println("v3 = {}", v3);
        // fmt::println("x  = {}", x);
        // fmt::println("");
        // fmt::println("Rs = {}", Rs);
        // fmt::println("");

        CHECK(std::abs(angle(v1, v2) - angle(v3, x)) < eps);
        CHECK(std::abs(nrm(v2) / nrm(v1) - nrm(x) / nrm(v3)) < eps);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec2d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2d: default init")
    {
        fmt::println("MVec2d: default init");
        // default initialization
        mvec2d v;
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.c0) < eps);
        CHECK(std::abs(v.c1) < eps);
        CHECK(std::abs(v.c2) < eps);
        CHECK(std::abs(v.c3) < eps);
    }
    TEST_CASE("MVec2d: with curly braced intializer")
    {
        fmt::println("MVec2d: with curly braced intializer");
        // default initialization
        mvec2d v{0.0, 1.0, 2.0, 3.0};
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.c0 - 0.0) < eps);
        CHECK(std::abs(v.c1 - 1.0) < eps);
        CHECK(std::abs(v.c2 - 2.0) < eps);
        CHECK(std::abs(v.c3 - 3.0) < eps);
    }

    TEST_CASE("MVec2d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("MVec2d: cp ctor & cp assign incl. type deduction");
        // default initialization
        mvec2d v1{0.0, 1.0, 2.0, 3.0}; // init with double (type deduction)
        mvec2d v2{v1};                 // cp ctor
        mvec2d v3 = v2;                // cp assign
        mvec2d v4 = -v3;               // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(std::abs(v2.c0 - 0.0) < eps);
        CHECK(std::abs(v2.c1 - 1.0) < eps);
        CHECK(std::abs(v2.c2 - 2.0) < eps);
        CHECK(std::abs(v2.c3 - 3.0) < eps);
        CHECK(std::abs(v3.c0 - 0.0) < eps);
        CHECK(std::abs(v3.c1 - 1.0) < eps);
        CHECK(std::abs(v3.c2 - 2.0) < eps);
        CHECK(std::abs(v3.c3 - 3.0) < eps);
        CHECK(v4 == -v3);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }

    TEST_CASE("MVec2d: defining basic types and ctor checks")
    {
        fmt::println("MVec2d: defining basic types and ctor checks");

        auto s = scalar2d(5.0);
        auto v = vec2d{1.0, 2.0};
        auto ps = pscalar2d(-5.0);

        auto mv1 = mvec2d{s};
        auto mv2 = mvec2d{v};
        auto mv3 = mvec2d{ps};
        auto mv4 = mvec2d{s, ps};
        auto mv5 = mvec2d_e{s};
        auto mv6 = mvec2d_e{ps};
        auto mv7 = mvec2d_e{s, ps};
        auto mv8 = mvec2d{mv7};
        auto mv9 = mvec2d{s, v, ps};

        // fmt::println("   mv1 = {}", mv1);
        // fmt::println("   mv2 = {}", mv2);
        // fmt::println("   mv3 = {}", mv3);
        // fmt::println("   mv4 = {}", mv4);
        // fmt::println("   mv5 = {}", mv5);
        // fmt::println("   mv6 = {}", mv6);
        // fmt::println("   mv7 = {}", mv7);
        // fmt::println("   mv8 = {}", mv8);
        // fmt::println("   mv9 = {}", mv9);

        CHECK(gr0(mv1) == s);
        CHECK(gr1(mv2) == v);
        CHECK(gr2(mv3) == ps);

        CHECK(gr0(mv4) == s);
        CHECK(gr2(mv4) == ps);

        CHECK(gr0(mv5) == s);
        CHECK(gr2(mv6) == ps);

        CHECK(gr0(mv7) == s);
        CHECK(gr2(mv7) == ps);

        CHECK(gr0(mv8) == s);
        CHECK(gr2(mv8) == ps);

        CHECK(gr0(mv9) == s);
        CHECK(gr1(mv9) == v);
        CHECK(gr2(mv9) == ps);

        CHECK(gr(s) == 0);
        CHECK(gr(v) == 1);
        CHECK(gr(ps) == 2);
    }

    TEST_CASE("MVec2d: fmt & cout printing")
    {
        fmt::println("MVec2d: fmt & cout printing:\n");

        mvec2d pf{1.0f, 2.00001f, 0.0f, 3.0f};
        mvec2d pd{1.0, 2.00001, 0.0, 3.0};

        std::cout << "   cout: pf = " << pf << std::endl;
        fmt::println("    fmt: pf = {}", pf);
        fmt::println("    fmt: pf = {:.4f}", pf);
        fmt::println("");

        std::cout << "   cout: pd = " << pd << std::endl;
        fmt::println("    fmt: pd = {}", pd);
        fmt::println("    fmt: pd = {:.4f}", pd);
        fmt::println("");

        std::vector<MVec2d<double>> vp1{{1.0, 1.0, 1.0, 2.0}, {0.5, 1.5, 2.0, 2.5}};
        fmt::println("    fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("    fmt: vp1 = {:g}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("MVec2d: vector space and linearity tests")
    {
        fmt::println("MVec2d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        mvec2d p0;
        mvec2d p1{0.0, 1.0, 2.0, 3.0};
        mvec2d p2{0.0, 2.0, 4.0, 6.0};
        mvec2d p3{0.0, 3.0, 6.0, 9.0};
        mvec2d p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t. multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t. addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
        CHECK(p4 == -p1);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec2d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2d: wedge product - basic properties")
    {
        fmt::println("MVec2d: wedge product - basic properties");

        vec2d v1{1.0, 2.0};
        vec2d v2{0.5, 3.0};
        vec2d v3{-2.0, 6.0};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar2d{sd};
        auto t = scalar2d{st};

        CHECK(wdg(v1, v1) == pscalar2d{});                    // wdg=0 for collin. vectors
        CHECK(wdg(v1, v2) == -wdg(v2, v1));                   // anticommutative for vect.
        CHECK(wdg(wdg(v1, v2), v3) == wdg(v1, wdg(v2, v3)));  // wdg is associative
        CHECK(wdg(v1, v2 + v3) == wdg(v1, v2) + wdg(v1, v3)); // wdg distributes over add.
        CHECK(wdg(v1 + v2, v3) == wdg(v1, v3) + wdg(v2, v3)); // wdg distributes over add.
        CHECK(wdg(sd * v1, v2) == wdg(v1, sd * v2)); // scalars can be factored out of wdg
        CHECK(wdg(sd * v1, v2) == sd * wdg(v1, v2)); // scalars can be factored out of wdg
        CHECK(wdg(s, t) == wdg(t, s));   // wdg between scalars equivalent to scalar mult.
        CHECK(wdg(s, v1) == wdg(v1, s)); // wdg between scalar and vector
        CHECK(wdg(s, v1) == sd * v1);    // wdg between scalar and vector
    }

    TEST_CASE("MVec2d: geometric product - basic properties")
    {
        fmt::println("MVec2d: geometric product - basic properties");

        vec2d v1{1.0, 2.0};
        vec2d v2{0.5, 3.0};
        vec2d v3{-2.0, 6.0};
        pscalar2d B1{-4};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar2d{sd};
        auto t = scalar2d{st};

        // nrm_sq(v1) = dot(v1,v1)  for every vector
        CHECK(nrm_sq(v1) == dot(v1, v1));

        CHECK(v1 * v2 == dot(v1, v2) + wdg(v1, v2)); // valid for vectors only
        CHECK(v1 * v2 == (v2 >> v1) + wdg(v1, v2));  // contraction = dot for same grades

        // valid also for gr(B1) != gr(v1)
        CHECK(v1 * B1 == gr1((B1 >> v1) + wdg(v1, B1)));
        CHECK(B1 * v1 == gr1((v1 << B1) + wdg(B1, v1)));

        // only valid for vectors (!):
        // dot = gpr symmetric part
        CHECK(dot(v1, v2) == gr0(0.5 * (v1 * v2 + v2 * v1)));
        // wdg = gpr anti-symmetric part
        CHECK(wdg(v1, v2) == gr2(0.5 * (v1 * v2 - v2 * v1)));

        // mathematical characteristics
        CHECK((v1 * v2) * v3 == v1 * (v2 * v3));    // gpr is associative
        CHECK(v1 * (v2 + v3) == v1 * v2 + v1 * v3); // gpr distributes over addition
        CHECK((v1 + v2) * v3 == v1 * v3 + v2 * v3); // wdg distributes over addition
        CHECK((sd * v1) * v2 == v1 * (sd * v2));    // scalars can be factored out of gpr
        CHECK((sd * v1) * v2 == sd * (v1 * v2));    // scalars can be factored out of gpr
        CHECK(s * t == t * s);    // gpr between scalars equivalent to scalar mult.
        CHECK(s * v1 == v1 * s);  // gpr between scalar and vector
        CHECK(s * v1 == sd * v1); // gpr between scalar and vector
    }

    TEST_CASE("MVec2d: geometric product tests")
    {
        fmt::println("MVec2d: geometric product tests");

        auto v1 = vec2d{1.0, 2.0};
        auto v2 = vec2d{0.5, 3.0};
        auto d12 = dot(v1, v2);
        auto w12 = wdg(v1, v2);

        mvec2d mv1{v1};
        mvec2d mv2{v2};
        auto wdp = 0.5 * (v1 * v2 + v2 * v1);
        auto wdp_mv12 = 0.5 * (mv1 * mv2 + mv2 * mv1);
        auto wdm = 0.5 * (v1 * v2 - v2 * v1);
        auto wdm_mv12 = 0.5 * (mv1 * mv2 - mv2 * mv1);

        // fmt::println("");
        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   dot(v1,v2) = {}", d12);
        // fmt::println("   wdg(v1,v2) = {}", w12);
        // fmt::println("");
        // fmt::println("   mv1 = {}", mv1);
        // fmt::println("   mv2 = {}", mv2);
        // fmt::println("   wdp      = 0.5*(v1 * v2 + v2 * v1)     = {}", wdp);
        // fmt::println("   gr0(wdp)                               = {}", gr0(wdp));
        // fmt::println("   wdp_mv12 = 0.5*(mv1 * mv2 + mv2 * mv1) = {}", wdp_mv12);
        // fmt::println("   wdm      = 0.5*(v1 * v2 - v2 * v1)     = {}", wdm);
        // fmt::println("   gr2(wdm)                               = {}", gr2(wdm));
        // fmt::println("   wdm_mv12 = 0.5*(mv1 * mv2 - mv2 * mv1) = {}", wdm_mv12);
        // fmt::println("");
        // fmt::println("   gr0(wdp_mv12) = {}", gr0(wdp_mv12));
        // fmt::println("   gr1(wdp_mv12) = {}", gr1(wdp_mv12));
        // fmt::println("   gr2(wdp_mv12) = {}", gr2(wdp_mv12));
        // fmt::println("");
        // fmt::println("   gr0(wdm_mv12) = {}", gr0(wdm_mv12));
        // fmt::println("   gr1(wdm_mv12) = {}", gr1(wdm_mv12));
        // fmt::println("   gr2(wdm_mv12) = {}", gr2(wdm_mv12));
        // fmt::println("");

        CHECK(d12 == gr0(0.5 * (v1 * v2 + v2 * v1)));
        CHECK(d12 == gr0(0.5 * (mv1 * mv2 + mv2 * mv1)));
        CHECK(w12 == gr2(0.5 * (v1 * v2 - v2 * v1)));
        CHECK(w12 == gr2(0.5 * (mv1 * mv2 - mv2 * mv1)));
        // redundant checks (just do avoid unused variable warnings)
        CHECK(d12 == gr0(wdp_mv12));
        CHECK(w12 == gr2(wdm_mv12));
        CHECK(d12 == gr0(wdp));
        CHECK(w12 == gr2(wdm));
    }

    TEST_CASE("MVec2d: geometric product - combinatorial tests")
    {
        fmt::println("MVec2d: geometric product - combinatorial tests");

        scalar2d s1{3.0};
        scalar2d s2{-1.5};
        vec2d v1{1.0, 2.0};
        vec2d v2{0.5, 3.0};
        pscalar2d ps1{-5.0};
        pscalar2d ps2{2.0};

        // fmt::println("s1  = {}", s1);
        // fmt::println("s2  = {}", s2);
        // fmt::println("v1  = {}", v1);
        // fmt::println("v2  = {}", v2);
        // fmt::println("ps1 = {}", ps1);
        // fmt::println("ps2 = {}", ps2);

        // checks all ctor combinations and equivalence of simplified products with
        // geometric product of a fully populated multivector
        //
        // assumption: the geometric product of a fully populated multivector serves as
        //             reference and is correctly implemented
        //             results are relative to full geometric product

        // MVec2d_E * MVec2d
        CHECK(mvec2d{s1, ps1} * mvec2d{s1, v1, ps1} ==
              mvec2d{mvec2d_e{s1, ps1} * mvec2d{s1, v1, ps1}});

        // MVec2d * MVec2d_E
        CHECK(mvec2d{s1, v1, ps1} * mvec2d{s1, ps1} ==
              mvec2d{mvec2d{s1, v1, ps1} * mvec2d_e{s1, ps1}});

        // MVec2d_E * MVec2d_E
        CHECK(mvec2d{s1, ps1} * mvec2d{s2, ps2} ==
              mvec2d{mvec2d_e{s1, ps1} * mvec2d_e{s2, ps2}});

        // MVec2d_E * Vec2d
        CHECK(mvec2d{s1, ps1} * mvec2d{v2} == mvec2d{mvec2d_e{s1, ps1} * v2});

        // Vec2d * MVec2d_E
        CHECK(mvec2d{v1} * mvec2d{s2, ps2} == mvec2d{v1 * mvec2d_e{s2, ps2}});

        // Vec2d * Vec2d
        CHECK(mvec2d{v1} * mvec2d{v2} == mvec2d{v1 * v2});

        // PScalar2d * MVec2d
        CHECK(mvec2d{ps1} * mvec2d{s2, v2, ps2} == mvec2d{ps1 * mvec2d{s2, v2, ps2}});

        // MVec2d * Pscalar2d
        CHECK(mvec2d{s1, v1, ps1} * mvec2d{ps2} == mvec2d{mvec2d{s1, v1, ps1} * ps2});

        // PScalar2d * MVec2d_E
        CHECK(mvec2d{ps1} * mvec2d{s2, ps2} == mvec2d{ps1 * mvec2d_e{s2, ps2}});

        // MVec2d_E * Pscalar2d
        CHECK(mvec2d{s1, ps1} * mvec2d{ps2} == mvec2d{mvec2d_e{s1, ps1} * ps2});

        // PScalar2d * Vec2d
        CHECK(mvec2d{ps1} * mvec2d{v2} == mvec2d{ps1 * v2});

        // Vec2d * Pscalar2d
        CHECK(mvec2d{v1} * mvec2d{ps2} == mvec2d{v1 * ps2});

        // Pscalar2d * Pscalar2d
        CHECK(mvec2d{ps1} * mvec2d{ps2} == mvec2d{ps1 * ps2});

        // Scalar * Scalar
        CHECK(mvec2d{s1} * mvec2d{s2} == mvec2d{s1 * s2});
    }

    TEST_CASE("MVec2d geometric product tests - recovering vectors from the "
              "geometric product")
    {
        fmt::println("MVec2d: geometric product tests - recovering vectors from the "
                     "geometric product");

        // Two multivectors mv1 and mv2 formed from vectors v1 and v2.
        // (gr0(mv1)==0 && gr1(mv1) != 0 && gr2(mv1)==0 &&
        //  gr0(mv2)==0 && gr1(mv2) != 0 && gr2(mv2)==0 )
        //
        // They are multiplied by the geometric product to form a multivector C
        // C = mv1(v1) * mv2(v2) = mv1 * mv2
        //
        // C contains a scalar part and a bivector part exclusively,
        // the remaining components are zero.
        // (gr0(C) != 0 && gr1(C)==0 && gr2(C) !=0)
        //
        // The scalar part of C represents the parts of v1 and v2
        // that are parallel to each other.
        // The bivector part of C represents the parts of v1 and v2
        // that are perpendicular to each other.
        //
        // multiply C from the right with inv(v2) recovers v1
        // multiply C from the left the the inv(v1) recovers v2

        vec2d a{1.0, 2.0};
        vec2d b{0.5, 3.0};
        mvec2d C{scalar2d(dot(a, b)), pscalar2d(wdg(a, b))};
        mvec2d gpr_right = C * mvec2d{inv(b)};
        mvec2d gpr_left = mvec2d{inv(a)} * C;

        // fmt::println("   a  = {}", a);
        // fmt::println("   b  = {}", b);
        // fmt::println("   C = a * b = {}", C);
        // fmt::println("");
        // fmt::println("   C * bi = gpr_right = {}", gpr_right);
        // fmt::println("   ai * C = gpr_left  = {}", gpr_left);
        // fmt::println("   gr1(gpr_right) = a = {}", gr1(gpr_right));
        // fmt::println("   gr1(gpr_left)  = b = {}", gr1(gpr_left));

        CHECK(a == gr1(gpr_right));
        CHECK(b == gr1(gpr_left));
    }

    TEST_CASE("MVec2d: geometric product tests - equivalence tests")
    {
        fmt::println("MVec2d: geometric product tests - equivalence tests");

        vec2d a{1.0, 2.0};
        vec2d b{0.5, 3.0};
        mvec2d mva{a};
        mvec2d mvb{b};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        mvec2d ab = a * b;
        mvec2d abm = mva * mvb;
        mvec2d abd{scalar2d(dot_ab), wdg_ab};

        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   ab  = a * b                      = {}", ab);
        // fmt::println("   abm = mva * mvb                  = {}", abm);
        // fmt::println("   abd = MVec2d(dot(a,b), wdg(a,b)) = {}", abd);

        CHECK(ab == abm);
        CHECK(ab == abd);
    }

    TEST_CASE("MVec2d: assignment tests")
    {
        fmt::println("MVec2d: assignment tests");

        vec2d v1{1.0, 2.0};
        vec2d v2{0.5, 3.0};

        mvec2d mv1{0.0, 1.0, 2.0, 0.0};
        mvec2d mv2{0.0, 0.5, 3.0, 0.0};
        mvec2d mv3{v1};
        mvec2d mv4 = v2;

        mvec2d mv5(scalar2d(5.0));
        mvec2d mv6(pscalar2d(6.0));

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("");
        // fmt::println("   mv1 = {}", mv1);
        // fmt::println("   mv2 = {}", mv2);
        // fmt::println("   mv3 = {}", mv3);
        // fmt::println("   mv4 = {}", mv4);
        // fmt::println("   mv5 = {}", mv5);
        // fmt::println("   mv6 = {}", mv6);
        // fmt::println("");
        // fmt::println("   gr1(mv1) = {}", gr1(mv1));
        // fmt::println("   gr1(mv2) = {}", gr1(mv2));
        // fmt::println("   gr1(mv3) = {}", gr1(mv3));
        // fmt::println("   gr1(mv3) = {}", gr1(mv4));

        CHECK(gr1(mv1) == v1);
        CHECK(gr1(mv2) == v2);
        CHECK(gr1(mv3) == v1);
        CHECK(gr1(mv4) == v2);
        CHECK(mv1 == mv3);
        CHECK(mv4 == mv2);
    }

    TEST_CASE("MVec2d: modeling complex numbers")
    {
        fmt::println("MVec2d: modeling complex numbers");

        vec2d v1{1.0, -1.0};
        mvec2d v1m{v1}; // full 2d multivector

        // multiplying with e1 from the left should make it a complex number
        // i.e. a multivector with a scalar (=Re) and a bivector part (=Im)
        // (for test purposes here, the even subalgebra would be sufficient)
        auto vc = e1_2d * v1;
        auto vcm = e1m_2d * v1m; // full gpr

        // multiplying with I2 from the right should rotate by +90°
        auto vr = vc * I_2d;
        auto vrm = vcm * Im_2d; // full gpr

        // multiplying with I2 from the left should rotate by -90°
        auto vl = I_2d * vc;
        auto vlm = Im_2d * vcm; // full gpr

        // defining a complex number in all three forms
        auto u = vec2d{1.0, 0.0};
        auto v = vec2d(std::cos(pi / 6.0), std::sin(pi / 6.0)); // unit vec +30%
        auto angle_uv = angle(u, v);

        auto uv = u * v; // complex number with real part and bivector part
        auto a = gr0(uv);
        auto b = gr2(uv);
        auto r = sqrt(a * a + b * b);

        // fmt::println("   I_2d          = {}", I_2d);
        // fmt::println("   Im_2d         = {}", Im_2d);
        // fmt::println("   I_2d * I_2d   = {}", I_2d * I_2d);
        // fmt::println("   Im_2d * Im_2d = {}", Im_2d * Im_2d);
        // fmt::println("");
        // fmt::println("   e1_2d  = {}", e1_2d);
        // fmt::println("   e1m_2d = {}", e1m_2d);
        // fmt::println("   e2_2d  = {}", e2_2d);
        // fmt::println("   e2m_2d = {}", e2m_2d);
        // fmt::println("");
        // fmt::println("   vc   = {}", vc);
        // fmt::println("   vcm  = {}", vcm);
        // fmt::println("   vr   = {}", vr);
        // fmt::println("   vrm  = {}", vrm);
        // fmt::println("   vl   = {}", vl);
        // fmt::println("   vlm  = {}", vlm);
        // fmt::println("");
        // fmt::println("   v1         = {}", v1);
        // fmt::println("   v1 * I_2d  = {}", v1 * I_2d);
        // fmt::println("   I_2d * v1  = {}", I_2d * v1);
        // fmt::println("");
        // fmt::println("   u        = {}", u);
        // fmt::println("   v        = {}", v);
        // fmt::println("   angle_uv = {:.3}°", angle_uv * 180 / pi);
        // fmt::println("");
        // fmt::println("   uv                  = {}", uv);
        // fmt::println("   a = gr0(uv)         = {}", a);
        // fmt::println("   b = gr2(uv)         = {}", b);
        // fmt::println("   r = sqrt(a^2 + b^2) = {}", r);
        // fmt::println("   r exp(angle_uv) = {}", r *
        // exp(pscalar2d(angle_uv))); HINT:declaring angle_uv a PScalar2d
        // makes it a bivector angle, i.e. a multiple of the bivector I_2d ATTENTION:
        // if you don't declare it as such, the normal exponential function
        //            will be called, resulting in a scalar result!

        CHECK(std::abs(angle_uv - pi / 6.0) < eps);
        CHECK(std::abs(r - 0.5 * std::sqrt(2.0)) < eps);
        CHECK(gr0(vc) == gr0(vcm));
        CHECK(gr2(vc) == gr2(vcm));
        CHECK(gr0(vr) == gr0(vrm));
        CHECK(gr2(vr) == gr2(vrm));
        CHECK(gr0(vl) == gr0(vlm));
        CHECK(gr2(vl) == gr2(vlm));
        CHECK(v1.x == (v1 * I_2d).y); // rotation +90°
        CHECK(v1.y == -(v1 * I_2d).x);
        CHECK(v1.x == -(I_2d * v1).y); // rotation -90°
        CHECK(v1.y == (I_2d * v1).x);
        CHECK(std::abs(nrm(vc) - std::sqrt(2.0)) < eps);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec2d_E<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2d_E: modelling complex numbers - basics")
    {
        fmt::println("MVec2d_E: modelling complex numbers - basics");

        // defining a complex number in all three forms as multivector
        auto u = vec2d{1.0, 0.0};
        auto v = vec2d(std::cos(pi / 6.0), std::sin(pi / 6.0)); // unit vec +30°

        auto angle_uv = angle(u, v);

        auto uv = u * v; // complex number with real part and bivector part
        auto v2 = exp(I_2d, angle_uv);
        auto re = gr0(uv);
        auto im = gr2(uv);
        auto r = sqrt(re * re + im * im);

        mvec2d_e a{1.0, 0.0};
        mvec2d_e b{1.0, 1.0};
        mvec2d_e c{a + b};
        mvec2d_e d{a - b};
        mvec2d_e e{2.0 * b};
        mvec2d_e f{b * 2.0};
        mvec2d_e g{-e};
        mvec2d_e h{0.0, 1.0};
        mvec2d_e as = a * a;
        mvec2d_e hs = h * h;

        mvec2d_e j = b * c;
        auto k = I_2d;
        mvec2d_e l = exp(I_2d, pi / 2);
        mvec2d_e m = Im_2d_E;
        mvec2d n = Im_2d;
        // fmt::println("   Multivector form of complex numbers:");
        // fmt::println("   u                     = {}", u);
        // fmt::println("   v                     = {}", v);
        // fmt::println("   angle(u,v)            = {:.3}°", angle_uv * 180 / pi);
        // fmt::println("   uv = u * v            = {}", uv);
        // fmt::println("   re = gr0(uv)          = {}", re);
        // fmt::println("   im = gr2(uv)          = {}", im);
        // fmt::println("   r = sqrt(re^2 + im^2) = {}", r);
        // fmt::println("");
        // fmt::println("   Using the even subalgebra only (std form of complex
        // numbers):");
        // // declaring angle_uv a PScalar2d makes it a bivector angle,
        // // i.e. a multiple of the bivector I_2d
        // // ATTENTION: if you don't declare it as such, the normal exponential
        // function
        // //            will be called, resulting in a scalar result!
        // fmt::println("   v2=exp(angle_uv) = {}", v2);
        // fmt::println("");
        // fmt::println("   a         = {}", a);
        // fmt::println("   b         = {}", b);
        // fmt::println("   c = a+b   = {}", c);
        // fmt::println("   d = a-b   = {}", d);
        // fmt::println("   e = 2.0*b = {}", e);
        // fmt::println("   f = b*2.0 = {}", f);
        // fmt::println("   g = -e    = {}", g);
        // fmt::println("");
        // fmt::println("   h =           = {}", h);
        // fmt::println("   as = a * a    = {}", as);
        // fmt::println("   hs = h * h    = {}", hs);
        // fmt::println("   b * h         = {}", b * h);
        // fmt::println("   h * b         = {}", h * b);
        // fmt::println("");
        // fmt::println("   j = b * c     = {}", j);
        // fmt::println("");
        // fmt::println("   k = I_2d                         = {}", k);
        // fmt::println("   l = exp(pscalar2d(pi/2)) = {:.3f}", l);
        // fmt::println("   m = Im_2d_E                      = {}", m);
        // fmt::println("   n = Im_2d                        = {}", n);

        CHECK(std::abs(r - 0.5 * std::sqrt(2.0)) < eps);
        CHECK(c == a + b);
        CHECK(d == a - b);
        CHECK(e == 2.0 * b);
        CHECK(f == b * 2.0);
        CHECK(g == -e);
        CHECK(as == a);
        CHECK(hs == mvec2d_e(-1.0, 0.0));
        CHECK(j == b * c);
        CHECK(k == I_2d);
        CHECK(v.x == v2.c0);
        CHECK(v.y == v2.c1);
        CHECK(b * h ==
              h * b); // the 2d pseudoscalar commutes commutes with complex numbers
        CHECK(l == m);
        CHECK(n == Im_2d);
        CHECK(rev(b + c) == rev(b) + rev(c));
        CHECK(rev(b * c) == rev(b) * rev(c));
        CHECK(nrm(b * c) == nrm(b) * nrm(c));
        CHECK(b * c == c * b);

        CHECK(std::abs(nrm_sq(mvec2d_e{1.0, 1.0}) - 2.0) < eps);
        CHECK(std::abs(nrm(mvec2d_e{1.0, 1.0}) - std::sqrt(2.0)) < eps);
        CHECK(rev(mvec2d_e{1.0, 1.0}) == mvec2d_e{1.0, -1.0});
        CHECK(std::abs(nrm(mvec2d_e{scalar2d(1.0), pscalar2d(1.0)}) - std::sqrt(2.0)) <
              eps);

        CHECK(mvec2d_e{-1.0, 1.0} * inv(mvec2d_e{-1.0, 1.0}) == mvec2d_e{1.0, 0.0});
        CHECK(std::abs(gr0(mvec2d_e{-1.0, 1.0} * rev(mvec2d_e{-1.0, 1.0})) -
                       nrm_sq(mvec2d_e{-1.0, 1.0})) < eps);
        CHECK(std::abs(gr2(mvec2d_e{-1.0, 1.0} * rev(mvec2d_e{-1.0, 1.0}))) < eps);

        CHECK(std::abs(angle_to_re(mvec2d_e{1.0, 0.0}) - 0.0) < eps);
        CHECK(std::abs(angle_to_re(mvec2d_e{1.0, 1.0}) - pi / 4.0) < eps);
        CHECK(std::abs(angle_to_re(mvec2d_e{0.0, 1.0}) - pi / 2.0) < eps);
        CHECK(std::abs(angle_to_re(mvec2d_e{-1.0, 1.0}) - pi * 3.0 / 4.0) < eps);
        CHECK(std::abs(angle_to_re(mvec2d_e{-1.0, 0.0}) - pi) < eps);
        CHECK(std::abs(angle_to_re(mvec2d_e{1.0, -1.0}) - (-pi / 4.0)) < eps);
        CHECK(std::abs(angle_to_re(mvec2d_e{0.0, -1.0}) - (-pi / 2.0)) < eps);
        CHECK(std::abs(angle_to_re(mvec2d_e{-1.0, -1.0}) - (-pi * 3.0 / 4.0)) < eps);

        CHECK(vec2d{1.0, 0.0} * vec2d{1.1, 1.1} ==
              rev(vec2d{1.1, 1.1} * vec2d{1.0, 0.0}));
        CHECK(exp(I_2d, pi / 4) == rev(exp(I_2d, -pi / 4)));
        CHECK(exp(I_2d, -angle_uv) * u == u * exp(I_2d, angle_uv)); // 2d rotation direct
        CHECK(exp(I_2d, -angle_uv) * u == v);
        CHECK(rotate(u, rotor(I_2d, angle_uv)) == v); // 2d rotation with double product
                                                      // completely as in the 3d case
                                                      // more effort computationally,
                                                      // but independent of dimension of
                                                      // space
    }

    TEST_CASE("MVec2d_E: modelling complex numbers - products")
    {
        fmt::println("MVec2d_E: modelling complex numbers - products");

        // std::vector<std::tuple<double, MVec2d_E<double>>> c_v;
        // for (int i = -12; i <= 12; ++i) {
        //     double phi = i * pi / 12;
        //     mvec2d_e c = exp(pscalar2d(phi));
        //     c_v.push_back(std::make_tuple(phi, c));
        //     fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
        //                  " angle={: .4f}",
        //                  i, phi, phi * 180 / pi, c, angle(c));
        // }
        // fmt::println("");

        // auto v = Vec2d<double>{1.0, 0.0};
        // // auto v = Vec2d<double>{1.0, 1.0};
        // for (auto const& [phi, c] : c_v) {
        //     auto u1 = v * c;
        //     auto u2 = c * v;
        //     fmt::println("   phi={: .4f}, phi={:> 4.0f}°, c={: .3f},"
        //                  "  u1={: .3f}, u2={: .3f}",
        //                  phi, phi * 180 / pi, c, u1, u2);
        // }
        // fmt::println("");


        CHECK(mvec2d_e{2.0, 3.0} * mvec2d{-1.0, 1.5, -2.0, -3.0} ==
              mvec2d{2.0, 0.0, 0.0, 3.0} * mvec2d{-1.0, 1.5, -2.0, -3.0});
        CHECK(mvec2d_e{2.0, 3.0} * vec2d{1.5, -2.0} ==
              gr1(mvec2d{2.0, 0.0, 0.0, 3.0} * mvec2d{0.0, 1.5, -2.0, 0.0}));

        CHECK(gr0(vec2d{1.5, -2.0} * vec2d{2.0, 3.0}) ==
              gr0(mvec2d{0.0, 1.5, -2.0, 0.0} * mvec2d{0.0, 2.0, 3.0, 0.0}));
        CHECK(gr2(vec2d{1.5, -2.0} * vec2d{2.0, 3.0}) ==
              gr2(mvec2d{0.0, 1.5, -2.0, 0.0} * mvec2d{0.0, 2.0, 3.0, 0.0}));

        // multiply from left
        CHECK(pscalar2d(1.5) * mvec2d{-1.0, 1.5, -2.0, -3.0} ==
              mvec2d{0.0, 0.0, 0.0, 1.5} * mvec2d{-1.0, 1.5, -2.0, -3.0});

        CHECK(mvec2d(pscalar2d(1.5) * mvec2d_e{-1.0, -3.0}) ==
              mvec2d{0.0, 0.0, 0.0, 1.5} * mvec2d{-1.0, 0.0, 0.0, -3.0});

        CHECK(mvec2d(pscalar2d(1.5) * vec2d{-1.0, -3.0}) ==
              mvec2d{0.0, 0.0, 0.0, 1.5} * mvec2d{0.0, -1.0, -3.0, 0.0});

        // multiply from right
        CHECK(mvec2d{-1.0, 1.5, -2.0, -3.0} * pscalar2d(1.5) ==
              mvec2d{-1.0, 1.5, -2.0, -3.0} * mvec2d{0.0, 0.0, 0.0, 1.5});

        CHECK(mvec2d_e{-1.0, -3.0} * mvec2d(pscalar2d(1.5)) ==
              mvec2d{-1.0, 0.0, 0.0, -3.0} * mvec2d{0.0, 0.0, 0.0, 1.5});

        CHECK(mvec2d(vec2d{-1.0, -3.0} * pscalar2d(1.5)) ==
              mvec2d{0.0, -1.0, -3.0, 0.0} * mvec2d{0.0, 0.0, 0.0, 1.5});

        // two bivectors
        CHECK(mvec2d(scalar2d(pscalar2d(1.5) * pscalar2d(3.0))) ==
              mvec2d{0.0, 0.0, 0.0, 1.5} * mvec2d{0.0, 0.0, 0.0, 3.0});

        // mvec2d_e tests multiply from left
        CHECK(mvec2d_e{-1.0, -3.0} * mvec2d(-1.0, 1.5, -2.0, -3.0) ==
              mvec2d{-1.0, 0.0, 0.0, -3.0} * mvec2d{-1.0, 1.5, -2.0, -3.0});

        CHECK(mvec2d(mvec2d_e{-1.0, -3.0} * vec2d(1.5, -2.0)) ==
              mvec2d{-1.0, 0.0, 0.0, -3.0} * mvec2d{0.0, 1.5, -2.0, 0.0});

        // mvec2d_e tests multiply from right
        CHECK(mvec2d(-1.0, 1.5, -2.0, -3.0) * mvec2d_e{-1.0, -3.0} ==
              mvec2d{-1.0, 1.5, -2.0, -3.0} * mvec2d{-1.0, 0.0, 0.0, -3.0});

        CHECK(mvec2d(vec2d(1.5, -2.0) * mvec2d_e{-1.0, -3.0}) ==
              mvec2d{0.0, 1.5, -2.0, 0.0} * mvec2d{-1.0, 0.0, 0.0, -3.0});

        // multiply two mvec2d_e
        CHECK(mvec2d(mvec2d_e(-3.0, 2.0) * mvec2d_e{-1.0, -3.0}) ==
              mvec2d{-3.0, 0.0, 0.0, 2.0} * mvec2d{-1.0, 0.0, 0.0, -3.0});

        auto m = mvec2d_e{13.0, 5.0};
        auto prd = m * inv(m);
        CHECK(std::abs(gr0(prd) - 1.0) < eps);
        CHECK(std::abs(gr2(prd) - 0.0) < eps);
    }


    TEST_CASE("MVec2d: dualization")
    {
        fmt::println("MVec2d: dualization");

        vec2d v{1.0, 2.0};                    // 2d vector
        mvec2d vm{10.0, 1.0, 2.0, 30.0};      // full 2d multivector
        mvec2d vm2{-7.0, 3.0, -42.0, 5.0};    // full 2d multivector
        mvec2d vm_even{10.0, 0.0, 0.0, 30.0}; // full 2d multivector - even content
        mvec2d_e vm_E{10.0, 30.0};            // even grade 2d multivector

        ////////////////////////////////////////////////////////////////////////////////
        // duality (as defined in Lengyel, "PGA illuminated"):
        ////////////////////////////////////////////////////////////////////////////////
        //
        // dual(A) = cmpl(A) in spaces of uneven dimension
        //         = rcmpl(A) in spaces of even dimension
        //

        auto vm_dual_manual = rev(vm) * Im_2d;
        auto vm_dual = dual(vm);

        auto vm_dual_even_manual = rev(vm_even) * Im_2d;
        auto vm_dual_even = dual(vm_even);

        auto vm_dual_manual_E = rev(vm_E) * Im_2d_E;
        auto vm_dual_E = dual(vm_E);

        auto v_dual_manual = rev(v) * I_2d;
        auto v_dual = dual(v);

        // fmt::println("   I_2d               = {}", I_2d);
        // fmt::println("   Im_2d              = {}", Im_2d);
        // fmt::println("   Im_2d_E            = {}", Im_2d_E);
        // fmt::println("");
        // fmt::println("   vm                 = {}", vm);
        // fmt::println("   vm*rev(Im_2d)      = {}", vm_dual_manual);
        // fmt::println("   dual(vm)           = {}", vm_dual);
        // fmt::println("");
        // fmt::println("   vm_even            = {}", vm_even);
        // fmt::println("   vm_even*rev(Im_2d) = {}", vm_dual_even_manual);
        // fmt::println("   dual(vm_even)      = {}", vm_dual_even);
        // fmt::println("");
        // fmt::println("   vm_E               = {}", vm_E);
        // fmt::println("   vm_E*rev(Im_2d_E)  = {}", vm_dual_manual_E);
        // fmt::println("   dual(vm_E)         = {}", vm_dual_E);
        // fmt::println("");
        // fmt::println("   v                  = {}", v);
        // fmt::println("   v*rev(I_2d)        = {}", v_dual_manual);
        // fmt::println("   dual(v)            = {}", v_dual);

        CHECK(vm_dual == vm_dual_manual);
        CHECK(vm_dual_even == vm_dual_even_manual);
        CHECK(vm_dual_E == vm_dual_manual_E);
        CHECK(v_dual == v_dual_manual);
        CHECK(dual(scalar2d(5)) == pscalar2d(5));
        CHECK(dual(pscalar2d(5)) == scalar2d(5));
        CHECK(dual(I_2d) == 1);

        // dual properties (A. Macdonald, p. 110):
        //
        // a) dual(aA) = a dual(A)
        // b) dual(A + B) = dual(A) + dual(B)
        // c) dual(dual(A)) = (-1)^(n*(n-1)/2) A   (with n as dimension of the (sub)space)
        // d) |dual(B)| = |B|
        // e) if B is a j-blade then dual(B) is an (n-j)-blade
        // f) if A is a j-vector then dual(A) is an (n-j)-vector
        //    (remember: a j-vector is a sum of j-blades, which are outer products)

        CHECK(dual(3.0 * vm) == 3.0 * dual(vm));
        CHECK(dual(vm + vm2) == dual(vm) + dual(vm2));
        CHECK(dual(I_2d) == 1);
        CHECK(dual(v) == (vec2d(-v.y, v.x))); // = rcmpl(v)

        // inner and outer products are in G^n are dual (A. Macdonald, p. 111):
        //
        // dual(dot(A,B)) = wdg(A,dual(B)), dual(wdg(A,B)) = dot(A,dual(B))

        auto a = vec2d{5.0, 1.0};
        auto b = vec2d{3.0, 3.0};

        // fmt::println("   a                 = {}", a);
        // fmt::println("   b                 = {}", b);
        // fmt::println("   dual(b)           = {}", dual(b));
        // fmt::println("   dot(a, b)         = {}", dot(a, b));
        // fmt::println("   wdg(a, b)         = {}", wdg(a, b));
        // fmt::println("");
        // fmt::println("   dual(dot(a, b))   = {}", dual(dot(a, b)));
        // fmt::println("   wdg(a, dual(b))   = {}", wdg(a, dual(b)));
        // fmt::println("");
        // fmt::println("   dual(wdg(a, b))   = {}", dual(wdg(a, b)));
        // fmt::println("   dot(a, dual(b))   = {}", dot(a, dual(b)));

        CHECK(dual(scalar2d(dot(a, b))) == wdg(a, dual(b)));
        CHECK(wdg(a, dual(b)) == dot(a, b) * I_2d);
    }

    TEST_CASE("MVec2d: product tests")
    {
        fmt::println("MVec2d: product tests");
        // tests based on "The inner products of geometric algebra", Leo Dorst

        auto s1 = scalar2d{2.0};
        auto v1 = vec2d{1.0, -3.0};
        auto ps1 = pscalar2d{-2.0};
        auto mv1 = mvec2d{s1, v1, ps1};

        auto s2 = scalar2d{-1.0};
        auto v2 = vec2d{2.0, 1.5};
        auto ps2 = pscalar2d{3.0};
        auto mv2 = mvec2d{s2, v2, ps2};

        auto s3 = scalar2d{-10.0};
        auto v3 = vec2d{-2.0, 4.0};
        auto ps3 = pscalar2d{-12.0};
        auto mv3 = mvec2d{s3, v3, ps3};

        auto V1 = mvec2d{v1};
        auto V2 = mvec2d{v2};
        auto V3 = mvec2d{v3};
        auto B = mvec2d{wdg(v1, v2)};

        // show values
        // fmt::println("\n   mv1      = {}, nrm_sq = {}, nrm = {}", mv1, nrm_sq(mv1),
        //              nrm(mv1));
        // fmt::println("   gr0(mv1) = {}, nrm_sq = {}, nrm = {}", gr0(mv1), nrm_sq(s1),
        //              nrm(s1));
        // fmt::println("   gr1(mv1) = {}, nrm_sq = {}, nrm = {}", gr1(mv1), nrm_sq(v1),
        //              nrm(v1));
        // fmt::println("   gr2(mv1) = {}, nrm_sq = {}, nrm = {}", gr2(mv1), nrm_sq(ps1),
        //              nrm(ps1));
        // fmt::println("");
        // fmt::println("   mv2      = {}, nrm_sq = {}, nrm = {}", mv2, nrm_sq(mv2),
        //              nrm(mv2));
        // fmt::println("   gr0(mv2) = {}, nrm_sq = {}, nrm = {}", gr0(mv2), nrm_sq(s2),
        //              nrm(s2));
        // fmt::println("   gr1(mv2) = {}, nrm_sq = {}, nrm = {}", gr1(mv2), nrm_sq(v2),
        //              nrm(v2));
        // fmt::println("   gr2(mv2) = {}, nrm_sq = {}, nrm = {}", gr2(mv2), nrm_sq(ps2),
        //              nrm(ps2));
        // fmt::println("");

        // scalar product
        // fmt::println("   gr0(gr0(mv1)*gr0(mv2))={}, gr0(gr0(mv1)*gr1(mv2))={}, "
        //              "gr0(gr0(mv1)*gr2(mv2))={}",
        //              gr0(mvec2d{gr0(mv1) * gr0(mv2)}), gr0(mvec2d{gr0(mv1) *
        //              gr1(mv2)}), gr0(mvec2d{gr0(mv1) * gr2(mv2)}));
        // fmt::println("   gr0(gr1(mv1)*gr0(mv2))={}, gr0(gr1(mv1)*gr1(mv2))={}, "
        //              "gr0(gr1(mv1)*gr2(mv2))={}",
        //              gr0(mvec2d{gr1(mv1) * gr0(mv2)}), gr0(mvec2d{gr1(mv1) *
        //              gr1(mv2)}), gr0(mvec2d{gr1(mv1) * gr2(mv2)}));
        // fmt::println("   gr0(gr2(mv1)*gr0(mv2))={}, gr0(gr2(mv1)*gr1(mv2))={}, "
        //              "gr0(gr2(mv1)*gr2(mv2))={}",
        //              gr0(mvec2d{gr2(mv1) * gr0(mv2)}), gr0(mvec2d{gr2(mv1) *
        //              gr1(mv2)}), gr0(mvec2d{gr2(mv1) * gr2(mv2)}));
        // fmt::println("");
        CHECK(nrm_sq(mv1) == nrm_sq(s1) + nrm_sq(v1) + nrm_sq(ps1));
        CHECK(nrm_sq(mv2) == nrm_sq(s2) + nrm_sq(v2) + nrm_sq(ps2));

        CHECK(dot(mv1, mv2) == dot(s1, s2) + dot(v1, v2) + dot(ps1, ps2));
        CHECK(dot(mv1, mv1) == nrm_sq(mv1));
        CHECK(dot(mv2, mv2) == nrm_sq(mv2));

        // wedge product
        CHECK(wdg(v1, v2) == gr2(wdg(mvec2d(v1), mvec2d(v2))));
        // check full permissible range of arguments, even the ones delivering 0)
        CHECK(wdg(v1, ps1) == scalar2d(0.0));
        CHECK(wdg(ps1, v1) == scalar2d(0.0));
        CHECK(wdg(ps1, ps1) == scalar2d(0.0));

        // contractions and comparison with other products
        // fmt::println("   V1      = {}", V1);
        // fmt::println("   V2      = {}", V2);
        // fmt::println("   B       = {}", B);
        // fmt::println("");
        // fmt::println("   V1 << B = {}", V1 << B);
        // fmt::println("   V2 << B = {}", V2 << B);
        // fmt::println("   B << V1 = {}", B << V1);
        // fmt::println("   B << V2 = {}", B << V2);
        // fmt::println("");
        // fmt::println("   B >> V1 = {}", B >> V1);
        // fmt::println("   B >> V2 = {}", B >> V2);
        // fmt::println("   V1 >> B = {}", V1 >> B);
        // fmt::println("   V2 >> B = {}", V2 >> B);
        // fmt::println("");
        // fmt::println("   V1 << V1      = {}", V1 << V1);
        // fmt::println("   V1 >> V1      = {}", V1 >> V1);
        // fmt::println("   V1 >> rev(V1) = {}", V1 >> rev(V1));
        // fmt::println("");
        // fmt::println("   B << B        = {}", B << B);
        // fmt::println("   B >> B        = {}", B >> B);
        // fmt::println("   B >> rev(B)   = {}", B >> rev(B));
        // fmt::println("");
        // fmt::println("   dot(V1,V2)    = {}", dot(V1, V2));
        // fmt::println("   dot(V2,V1)    = {}", dot(V2, V1));
        // fmt::println("   V1 << V2      = {}", V1 << V2);
        // fmt::println("   V2 << V1      = {}", V2 << V1);
        // fmt::println("   V1 >> V2      = {}", V1 >> V2);
        // fmt::println("   V2 >> V2      = {}", V2 >> V1);
        // fmt::println("");
        // fmt::println("   V1 << Im_2d   = {}", V1 << Im_2d);
        // fmt::println("   Im_2d >> V1   = {}", Im_2d >> V2);
        // fmt::println("   B << Im2d     = {}", B << Im_2d);
        // fmt::println("   Im2d >> B     = {}", Im_2d >> B);
        // fmt::println("");

        // connection between products (2.2.5)
        CHECK(dot(wdg(mv1, mv2), mv3) == dot(mv1, mv2 << mv3));
        CHECK(dot(mv3, wdg(mv2, mv1)) == dot(mv3 >> mv2, mv1));
        CHECK(rev(mv1 >> mv2) == rev(mv2) << rev(mv1));
        CHECK((V1 >> B) == mvec2d{0.0, 0.0, 0.0, 0.0});
        CHECK(dot(V1, B) == 0);

        // understanding the contraction (2.3.1)
        // fmt::println("   s1                = {}", s1);
        // fmt::println("   mv2               = {}", mv2);
        // fmt::println("   mvec2d(s1) << mv2 = {}", mvec2d{s1} << mv2);
        // fmt::println("   value_t(s1) * mv2 = {}", value_t(s1) * mv2);
        // fmt::println("");
        CHECK((mvec2d{s1} << mv2) == value_t(s1) * mv2);

        CHECK((mvec2d{s2, v2, pscalar2d{0.0}} << mvec2d{s1}) ==
              mvec2d{s1 * s2, 0.0, 0.0, 0.0});
        CHECK((mvec2d{v2} << mvec2d{s1}) == mvec2d{0.0, 0.0, 0.0, 0.0});

        CHECK((wdg(v2, v3) >> v1) == wdg(v1 << v2, v3) + wdg(gr_inv(v2), v1 << v3));
        // CHECK((wdg(v1, v2) << v3) == (v1 << scalar2d(v2 << v3)));

        // contractions - check full permissible range of arguments, even the ones
        // delivering 0 as a result
        CHECK((ps1 << s1) == 0.0);
        CHECK((ps1 << v1) == 0.0);
        CHECK((v1 << s1) == 0.0);
        //
        CHECK((s1 >> ps1) == 0.0);
        CHECK((v1 >> ps1) == 0.0);
        CHECK((s1 >> v1) == 0.0);

        // 2.3.2
        CHECK(V1 * mv2 == (mv2 >> V1) + wdg(V1, mv2));
        CHECK(mv2 * V1 == (V1 << mv2) + wdg(mv2, V1));
        CHECK((mv2 >> V1) == -(V1 << gr_inv(mv2))); // replacing the right contraction
        CHECK(V1 * mv2 == -(V1 << gr_inv(mv2)) + wdg(V1, mv2)); // yep, replacing works

        // express contraction and wedge product through the geometric product
        CHECK((mv2 >> V1) == 0.5 * (V1 * mv2 - gr_inv(mv2) * V1));
        CHECK(wdg(V1, mv2) == 0.5 * (V1 * mv2 + gr_inv(mv2) * V1));
        CHECK(wdg(mv2, V1) == 0.5 * (mv2 * V1 + V1 * gr_inv(mv2)));

        // cross-check direct implementation of rwdg by comparing with wdg
        CHECK(rwdg(mv1, mv2) == lcmpl(wdg(rcmpl(mv1), rcmpl(mv2))));
    }

    TEST_CASE("MVec2d: simple applications & complements")
    {
        fmt::println("MVec2d: simple applications & complements");

        auto s1 = scalar2d{2.0};
        auto v1 = vec2d{1.0, -3.0};
        auto ps1 = pscalar2d{-2.0};
        auto mv1 = mvec2d{s1, v1, ps1};

        // simple projections
        auto v = vec2d{4.0, 3.5};
        auto u = vec2d{-1.0, 2.3};

        auto v_in_u = project_onto(v, u);
        auto v_perp_u = reject_from(v, u);

        // fmt::println("   v         = {:.3f}", v);
        // fmt::println("   u         = {:.3f}", u);
        // fmt::println("");
        // fmt::println("   v_in_u    = {:.3f}", v_in_u);
        // fmt::println("   v_perp_u  = {:.3f}", v_perp_u);
        // fmt::println("");

        // x^B = 0 for every point x in B
        CHECK(nrm_sq(wdg(v_in_u, u)) < eps);

        // v_perp_u should be proportional to the normal vector n of u
        // n = dual(u) or n = cmpl(u)
        // thus, wdg(v_perp_B, n) == 0 is required
        CHECK(nrm_sq(wdg(v_perp_u, dual(u))) < eps);
        CHECK(nrm_sq(wdg(v_perp_u, rcmpl(u))) < eps);
        CHECK(nrm_sq(wdg(lcmpl(u), v_perp_u)) < eps);

        // v_in_B and v_perp_B should be perpendicular to each other
        CHECK(nrm_sq(dot(v_in_u, v_perp_u)) < eps);

        // v should be the sum of v_in_B and v_perp_B
        CHECK(v == v_in_u + v_perp_u);

        // check complements
        // fmt::println("");
        // auto s = scalar2d(1);
        // auto e1 = e1_2d;
        // auto e2 = e2_2d;
        // auto ps = I_2d;
        // fmt::println("   u:  s, e1, e2, ps \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              s, e1, e2, ps);
        // fmt::println("");
        // fmt::println("   rcmpl(u):  rcmpl(s), rcmpl(e1), rcmpl(e2), rcmpl(ps) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              rcmpl(s), rcmpl(e1), rcmpl(e2), rcmpl(ps));
        // fmt::println("");
        // fmt::println(
        //     "   rev(u)*I_2d:  rev(s)*I_2d, rev(e1)*I_2d, rev(e2)*I_2d, rev(ps)*I_2d \n"
        //     "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //     rev(s) * I_2d, rev(e1) * I_2d, rev(e2) * I_2d, rev(ps) * I_2d);
        // fmt::println("");
        // fmt::println("   lcmpl(u):  lcmpl(s), lcmpl(e1), lcmpl(e2), lcmpl(ps) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              lcmpl(s), lcmpl(e1), lcmpl(e2), lcmpl(ps));
        // fmt::println("");
        // fmt::println(
        //     "   I_2d*rev(u):  I_2d*rev(s), I_2d*rev(e1), I_2d*rev(e2), I_2d*rev(ps) \n"
        //     "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //     I_2d * rev(s), I_2d * rev(e1), I_2d * rev(e2), I_2d * rev(ps));
        // fmt::println("");

        CHECK(rcmpl(scalar2d(1.0)) == I_2d);
        CHECK(rcmpl(scalar2d(1.0)) == rev(scalar2d(1.0)) * I_2d);
        CHECK(lcmpl(scalar2d(1.0)) == I_2d);
        CHECK(lcmpl(scalar2d(1.0)) == I_2d * rev(scalar2d(1.0)));
        CHECK(rcmpl(e1_2d) == e2_2d);
        CHECK(rcmpl(e1_2d) == rev(e1_2d) * I_2d);
        CHECK(lcmpl(e1_2d) == -e2_2d);
        CHECK(lcmpl(e1_2d) == I_2d * rev(e1_2d));
        CHECK(rcmpl(e2_2d) == -e1_2d);
        CHECK(rcmpl(e2_2d) == rev(e2_2d) * I_2d);
        CHECK(lcmpl(e2_2d) == e1_2d);
        CHECK(lcmpl(e2_2d) == I_2d * rev(e2_2d));
        CHECK(rcmpl(I_2d) == scalar2d(1.0));
        CHECK(rcmpl(I_2d) == rev(I_2d) * I_2d);
        CHECK(lcmpl(I_2d) == scalar2d(1.0));
        CHECK(lcmpl(I_2d) == I_2d * rev(I_2d));
        //
        CHECK(lcmpl(rcmpl(mv1)) == mv1);
        CHECK(lcmpl(rcmpl(mvec2d_e(s1, ps1))) == mvec2d_e(s1, ps1));
        //
        CHECK(wdg(scalar2d(5), rcmpl(scalar2d(5))) / nrm_sq(scalar2d(5)) == I_2d);
        CHECK(wdg(lcmpl(scalar2d(5)), scalar2d(5)) / nrm_sq(scalar2d(5)) == I_2d);
        CHECK(wdg(v, rcmpl(v)) / nrm_sq(v) == I_2d);
        CHECK(wdg(lcmpl(v), v) / nrm_sq(v) == I_2d);
        CHECK(wdg(pscalar2d(3), rcmpl(pscalar2d(3))) / nrm_sq(pscalar2d(3)) == I_2d);
        CHECK(wdg(lcmpl(pscalar2d(3)), pscalar2d(3)) / nrm_sq(pscalar2d(3)) == I_2d);

        CHECK((rev(v) * I_2d) == rcmpl(v));
        CHECK((I_2d * rev(v)) == lcmpl(v));
    }

    TEST_CASE("G<3,0,0>: ega3d")
    {
        fmt::println("");
        fmt::println("G<3,0,0>: ega3d");
        // 3d euklidean geometric algebra
        const algebra<3> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 3);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 8);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 4); // == dim_space + 1
        fmt::println("   ega3d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   ega3d: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3d: default init")
    {
        fmt::println("Vec3d: default init");
        vec3d v;
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.x) < eps);
        CHECK(std::abs(v.y) < eps);
        CHECK(std::abs(v.z) < eps);
    }
    TEST_CASE("Vec3d: with curly braced intializer")
    {
        fmt::println("Vec3d: with curly braced intializer");
        vec3d v{0.0, 0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.x) < eps);
        CHECK(std::abs(v.y) < eps);
        CHECK(std::abs(v.z) < eps);
    }
    TEST_CASE("Vec3d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("Vec3d: cp ctor & cp assign incl. type deduction");
        vec3d v1{1.0, 2.0, 3.0}; // init with double (type deduction)
        vec3d v2{v1};            // cp ctor
        vec3d v3 = v2;           // cp assign
        vec3d v4 = -v2;          // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(std::abs(v1.x - 1.0) < eps);
        CHECK(std::abs(v1.y - 2.0) < eps);
        CHECK(std::abs(v1.z - 3.0) < eps);
        CHECK(std::abs(v2.x - 1.0) < eps);
        CHECK(std::abs(v2.y - 2.0) < eps);
        CHECK(std::abs(v2.z - 3.0) < eps);
        CHECK(std::abs(v3.x - 1.0) < eps);
        CHECK(std::abs(v3.y - 2.0) < eps);
        CHECK(std::abs(v3.z - 3.0) < eps);
        CHECK(v4 == -v2);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }


    TEST_CASE("Vec3d: fmt & cout printing")
    {
        fmt::println("Vec3d: fmt & cout printing:\n");

        vec3d pf{1.0f, 2.00001f, 3.0f};
        vec3d pd{1.0, 2.00001, 3.0};

        // std::cout << "       cout: pf = " << pf << std::endl;
        fmt::println("       fmt:  pf = {}", pf);
        fmt::println("       fmt:  pf = {:.4f}", pf);
        fmt::println("");

        // std::cout << "       cout: pd = " << pd << std::endl;
        fmt::println("       fmt:  pd = {}", pd);
        fmt::println("       fmt:  pd = {:.4f}", pd);
        fmt::println("");

        std::vector<Vec3d<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
        fmt::println("       fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("       fmt: vp1 = {:g}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("Vec3d: comparison float")
    {
        fmt::println("Vec3d: comparison float");

        Vec3d<float> v1f{1.0f, 2.0f, 3.0f};
        Vec3d<float> v2f{2.0f, 4.0f, 3.0f};
        Vec3d<float> v3f{1.0f, 2.0000001f, 3.0f};
        Vec3d<float> v4f{v1f};

        // fmt::println("   v1f = {}", v1f);
        // fmt::println("   v2f = {}", v2f);
        // fmt::println("   v3f = {}", v3f);
        // fmt::println("   v4f = {}", v4f);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<float>::epsilon());

        CHECK(v1f == v4f);           // comparison (equality)
        CHECK(v1f != v2f);           // comparison (inequality)
        CHECK(nrm(v1f) < nrm(v2f));  // comparison (less than)
        CHECK(nrm(v2f) >= nrm(v1f)); // comparison (greater than or equal)
        CHECK(v3f == v1f);           // comparison (eqality)
    }

    TEST_CASE("Vec3d: comparison double")
    {
        fmt::println("Vec3d: comparison double");

        Vec3d<double> v1d{1.0, 2.0, 3.0};
        Vec3d<double> v2d{2.0, 4.0, 3.0};
        Vec3d<double> v3d{1.0, 2.0000000000000001, 3.0};
        Vec3d<double> v4d{v1d};

        // fmt::println("   v1d = {}", v1d);
        // fmt::println("   v2d = {}", v2d);
        // fmt::println("   v3d = {}", v3d);
        // fmt::println("   v4d = {}", v4d);

        // fmt::println("    fmt: eps = {}", std::numeric_limits<double>::epsilon());

        CHECK(v1d == v4d);           // comparison (equality)
        CHECK(v1d != v2d);           // comparison (inequality)
        CHECK(nrm(v1d) < nrm(v2d));  // comparison norm
        CHECK(nrm(v2d) >= nrm(v1d)); // comparison norm
        CHECK(v3d == v1d);           // comparison (eqality)
    }

    TEST_CASE("Vec3d: vector space and linearity tests")
    {
        fmt::println("Vec3d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        vec3d p0;
        vec3d p1{1.0, 2.0, 3.0};
        vec3d p2{2.0, 4.0, 6.0};
        vec3d p3{3.0, 6.0, 9.0};
        vec3d p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t. multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t.addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
        CHECK(p4 == -p1);
    }

    TEST_CASE("Vec3d: inner product properties")
    {
        fmt::println("Vec3d: inner product properties");

        double a = 2.35;
        vec3d u{1.0, 2.0, 1.0};
        vec3d v{-0.5, 3.0, 0.5};
        vec3d w{3.0, 6.0, -3.0};

        CHECK(dot(a * u, v) == a * dot(u, v));
        CHECK(dot(u + v, w) == dot(u, w) + dot(v, w));
        CHECK(dot(u, v) == dot(v, u));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3d: operations - norm, inverse, dot")
    {
        fmt::println("Vec3d: operations - norm, inverse, dot");

        scalar3d s1{3.2};

        vec3d v1{2.0, 1.0, 2.0};
        vec3d v2{normalize(v1)};

        vec3d v3{2.0, 6.0, -4.0};
        vec3d v4{inv(v3)};

        bivec3d b1{-2.3, 1.2, 4.5};
        pscalar3d ps1{-4.7};
        mvec3d_e mve1{s1, b1};
        mvec3d_u mvu1{v1, ps1};
        mvec3d mv1{s1, v1, b1, ps1};

        // fmt::println("v1 = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2 = normalize(v1) = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        // fmt::println("v3 = {: .4f}, nrm(v1) = {: .4f}", v3, nrm(v3));
        // fmt::println(
        //     "v4 = inv(v3) = {: .4f}, nrm(v3) = {: .4f}, nrm(v3)*nrm(v4) = {: .4f}", v4,
        //     nrm(v4), nrm(v3) * nrm(v4));

        CHECK(std::abs(nrm_sq(v1) - 9.0) < eps);
        CHECK(std::abs(nrm_sq(v2) - 1.0) < eps);
        CHECK(std::abs(dot(v4, v3) - 1.0) < eps);

        // check inverses - scalar
        // fmt::println("");
        // fmt::println("s1 * inv(s1) = {}", s1 * inv(s1)); // s
        CHECK(std::abs(nrm(s1 * inv(s1)) - 1) < eps);
        CHECK(std::abs(inv(s1) - rev(s1) / nrm_sq(s1)) < eps);

        // check inverses - vector
        // fmt::println("v1 * inv(v1) = {}", v1 * inv(v1)); // mv_e
        CHECK(std::abs(nrm(gr0(v1 * inv(v1))) - 1) < eps);
        CHECK(std::abs(nrm(gr2(v1 * inv(v1))) - 0) < eps);
        CHECK(std::abs(nrm(inv(v1) - rev(v1) / nrm_sq(v1))) < eps);

        // check inverses - bivector
        // fmt::println("b1 * inv(b1) = {}", b1 * inv(b1)); // mv_e
        CHECK(std::abs(nrm(gr0(b1 * inv(b1))) - 1) < eps);
        CHECK(std::abs(nrm(gr2(b1 * inv(b1))) - 0) < eps);
        CHECK(std::abs(nrm(inv(b1) - rev(b1) / nrm_sq(b1))) < eps);

        // check inverses - pseudoscalar
        // fmt::println("ps1 * inv(ps1) = {}", ps1 * inv(ps1)); // s
        CHECK(std::abs(nrm(ps1 * inv(ps1)) - 1) < eps);
        CHECK(std::abs(nrm(inv(ps1) - rev(ps1) / nrm_sq(ps1))) < eps);

        // check inverses - even grade multivector
        // fmt::println("mve1 * inv(mve1) = {}", mve1 * inv(mve1)); // mv_e
        CHECK(std::abs(nrm(gr0(mve1 * inv(mve1))) - 1) < eps);
        CHECK(std::abs(nrm(gr2(mve1 * inv(mve1))) - 0) < eps);
        CHECK(std::abs(nrm(inv(mve1) - rev(mve1) / nrm_sq(mve1))) < eps);

        // check inverses - uneven grade multivector
        // fmt::println("mvu1 * inv(mvu1) = {}", mvu1 * inv(mvu1)); // mv_e
        CHECK(std::abs(nrm(gr0(mvu1 * inv(mvu1))) - 1) < eps);
        CHECK(std::abs(nrm(gr2(mvu1 * inv(mvu1))) - 0) < eps);
        CHECK(std::abs(nrm(inv(mvu1) - rev(mvu1) / nrm_sq(mvu1))) < eps);

        // check inverses - multivector
        // fmt::println("mv1 * inv(mv1) = {}", mv1 * inv(mv1)); // mv
        CHECK(std::abs(nrm(gr0(mv1 * inv(mv1))) - 1) < eps);
        CHECK(std::abs(nrm(gr1(mv1 * inv(mv1))) - 0) < eps);
        CHECK(std::abs(nrm(gr2(mv1 * inv(mv1))) - 0) < eps);
        CHECK(std::abs(nrm(gr0(inv(mv1) * mv1)) - 1) < eps); // left and right inverse
                                                             // are equal

        // fmt::println("");
    }

    TEST_CASE("Vec3d: operations - angle I")
    {
        fmt::println("Vec3d: operations - angle");

        vec3d v1{1.0, 0.0, 0.0};
        vec3d v2{normalize(vec3d(1.0, 1.0, 0.0))};
        vec3d v3{0.0, 1.0, 0.0};
        vec3d v4{normalize(vec3d(-1.0, 1.0, 0.0))};
        vec3d v5{-1.0, 0.0, 0.0};
        vec3d v6{normalize(vec3d(-1.0, -1.0, 0.0))};
        vec3d v7{0.0, -1.0, 0.0};
        vec3d v8{normalize(vec3d(1.0, -1.0, 0.0))};

        // fmt::println("v1 = {: .4f}, nrm(v1) = {:.8f}, "
        //              "angle(v1,v1) = {:.8f}, {:.8f}",
        //              v1, nrm(v1), angle(v1, v1), angle(v1, v1) / pi);
        // fmt::println("v2 = {: .4f}, nrm(v2) = {:.8f}, "
        //              "angle(v1,v2) = {:.8f}, {:.8f}",
        //              v2, nrm(v2), angle(v1, v2), angle(v1, v2) / pi);
        // fmt::println("v3 = {: .4f}, nrm(v3) = {:.8f}, "
        //              "angle(v1,v3) = {:.8f}, {:.8f} ",
        //              v3, nrm(v3), angle(v1, v3), angle(v1, v3) / pi);
        // fmt::println("v4 = {: .4f}, nrm(v4) = {:.8f}, "
        //              "angle(v1,v4) = {:.8f}, {:.8f} ",
        //              v4, nrm(v4), angle(v1, v4), angle(v1, v4) / pi);
        // fmt::println("v5 = {: .4f}, nrm(v5) = {:.8f}, "
        //              "angle(v1,v5) = {:.8f}, {:.8f} ",
        //              v5, nrm(v5), angle(v1, v5), angle(v1, v5) / pi);
        // fmt::println("v6 = {: .4f}, nrm(v6) = {:.8f}, "
        //              "angle(v1,v6) = {:.8f}, {:.8f} ",
        //              v6, nrm(v6), angle(v1, v6), angle(v1, v6) / pi);
        // fmt::println("v7 = {: .4f}, nrm(v7) = {:.8f}, "
        //              "angle(v1,v7) = {:.8f}, {:.8f} ",
        //              v7, nrm(v7), angle(v1, v7), angle(v1, v7) / pi);
        // fmt::println("v8 = {: .4f}, nrm(v8) = {:.8f}, "
        //              "angle(v1,v8) = {:.8f}, {:.8f} ",
        //              v8, nrm(v8), angle(v1, v8), angle(v1, v8) / pi);

        CHECK(std::abs(angle(v1, v1) - 0.0) < eps);
        CHECK(std::abs(angle(v1, v2) - pi * 0.25) < eps);
        CHECK(std::abs(angle(v1, v3) - pi * 0.5) < eps);
        CHECK(std::abs(angle(v1, v4) - pi * 0.75) < eps);
        CHECK(std::abs(angle(v1, v5) - pi) < eps);

        // just to suppress unused variable warnings
        CHECK(v6 == normalize(vec3d(-1.0, -1.0, 0.0)));
        CHECK(v8 == normalize(vec3d(1.0, -1.0, 0.0)));
    }

    TEST_CASE("Vec3d: operations - angle II")
    {
        fmt::println("Vec3d: operations - angle II");

        std::vector<std::tuple<double, Vec3d<double>>> v1;
        std::vector<std::tuple<double, Vec3d<double>>> v2;
        std::vector<std::tuple<double, Vec3d<double>>> v3;

        // only positive angles are easy to implement vs. the 2d case

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3d<double>(std::cos(phi), std::sin(phi), 0.0);
            v1.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_3d, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3d<double>(std::cos(phi + pi / 2), std::sin(phi + pi / 2), 0.0);
            v2.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e2_3d, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec3d<double>(std::cos(phi + pi / 4), std::sin(phi + pi / 4), 0.0);
            v3.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_3d + e2_3d, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v1) {
            CHECK(std::abs(phi - angle(e1_3d, c)) < eps);
        }
        for (auto const& [phi, c] : v2) {
            CHECK(std::abs(phi - angle(e2_3d, c)) < eps);
        }
        auto ref_vec = normalize(e1_3d + e2_3d);
        for (auto const& [phi, c] : v3) {
            CHECK(std::abs(phi - angle(ref_vec, c)) < eps);
        }
    }

    TEST_CASE("Vec3d: operations - wedge")
    {
        fmt::println("Vec3d: operations - wedge");

        vec3d v1{1.0, 0.0, 0.0};
        vec3d v2{normalize(vec3d(1.0, 1.0, 0.0))};
        vec3d v3{0.0, 1.0, 0.0};
        vec3d v4{normalize(vec3d(-1.0, 1.0, 0.0))};
        vec3d v5{-1.0, 0.0, 0.0};
        vec3d v6{normalize(vec3d(-1.0, -1.0, 0.0))};
        vec3d v7{0.0, -1.0, 0.0};
        vec3d v8{normalize(vec3d(1.0, -1.0, 0.0))};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3d{sd};
        auto t = scalar3d{st};

        // fmt::println("v1 = {: .4f}, wdg(v1,v1) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v1, wdg(v1, v1), angle(v1, v1));
        // fmt::println("v2 = {: .4f}, wdg(v1,v2) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v2, wdg(v1, v2), angle(v1, v2));
        // fmt::println("v3 = {: .4f}, wdg(v1,v3) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v3, wdg(v1, v3), angle(v1, v3));
        // fmt::println("v4 = {: .4f}, wdg(v1,v4) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v4, wdg(v1, v4), angle(v1, v4));
        // fmt::println("v5 = {: .4f}, wdg(v1,v5) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v5, wdg(v1, v5), angle(v1, v5));
        // fmt::println("v6 = {: .4f}, wdg(v1,v6) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v6, wdg(v1, v6), angle(v1, v6));
        // fmt::println("v7 = {: .4f}, wdg(v1,v7) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v7, wdg(v1, v7), angle(v1, v7));
        // fmt::println("v8 = {: .4f}, wdg(v1,v8) = {: .4f}, "
        //              "angle = {: .4f}",
        //              v8, wdg(v1, v8), angle(v1, v8));

        CHECK(std::abs(nrm(wdg(v1, v1)) - sin(angle(v1, v1))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v2)) - sin(angle(v1, v2))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v3)) - sin(angle(v1, v3))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v4)) - sin(angle(v1, v4))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v5)) - sin(angle(v1, v5))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v6)) - sin(angle(v1, v6))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v7)) - sin(angle(v1, v7))) < eps);
        CHECK(std::abs(nrm(wdg(v1, v8)) - sin(angle(v1, v8))) < eps);
    }

    TEST_CASE("Vec3d: operations - project / reject / reflect (vector - vector)")
    {
        fmt::println("Vec3d: operations - project / reject / reflect (vector - vector)");

        vec3d v1{5.0, 1.0, 1.0};
        vec3d v2{2.0, 2.0, 1.0};

        vec3d v3{project_onto(v1, v2)};
        vec3d v4{reject_from(v1, v2)};
        vec3d v5{v3 + v4};


        // fmt::println("");
        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        // fmt::println("v2u = {: .4f}, nrm(v2) = {: .4f}", v2u, nrm(v2u));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .4f}, nrm(v3) = {: .4f}", v3,
        //              nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .4f}, nrm(v4) = {: .4f}", v4,
        //              nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .4f}, nrm(v5) = {: .4f}", v5,
        //              nrm(v5));
        // fmt::println("");
        // // this helps to understand, why the dot-product is sufficient
        // auto w = wdg(v1, v2);
        // auto i = inv(v2);
        // fmt::println("wdg(v1,v2)         = {: .4f}", w);
        // fmt::println("inv(v2)            = {: .4f}", i);
        // fmt::println("wdg(v1,v2)*inv(v2) = {: .4f}", w * i);
        // fmt::println("");

        vec3d v{4.0, 1.0, 1.0};
        vec3d b{e2_3d};
        auto B = bivec3d{e12_3d};

        // auto UB = bivec3d{e23_3d + e12_3d};
        // fmt::println("v   = {}", v);
        // fmt::println("b   = {}", b);
        // fmt::println("B   = {}", B);
        // fmt::println("UB  = {}", UB);
        // fmt::println("reflect_on_vec(v,b)  = {}", reflect_on_vec(v, b));
        // fmt::println("reflect_on(v,e3_3d)  = {}", reflect_on(v, e3_3d));
        // fmt::println("reflect_on(v,B)  = {}", reflect_on(v, B));
        // fmt::println("reflect_on(UB,B) = {}", reflect_on(UB, B));
        // fmt::println("");

        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);

        // just to suppress unused variable warnings
        CHECK(b == e2_3d);

        CHECK(reflect_on(v, e3_3d) == reflect_on(v, B));

        // checking time required
        //
        // auto start = std::chrono::system_clock::now();
        // for (size_t i = 0; i < 10000000; ++i) {
        //     vec3d v = reject_from(v1, v2);
        // }
        // auto end = std::chrono::system_clock::now();
        // auto elapsed =
        // std::chrono::duration_cast<std::chrono::milliseconds>(end
        // -
        // start); fmt::println("The measurement took {}", elapsed);

        // point reflected on a vector
        vec3d p{4, 1, 0};
        CHECK(reflect_on_vec(p, x_axis_3d) == vec3d{4, -1, 0});

        // coordinate axis reflected on perpendicular axis yield their negatives
        CHECK(reflect_on_vec(y_axis_3d, x_axis_3d) == -y_axis_3d);
        CHECK(reflect_on_vec(z_axis_3d, x_axis_3d) == -z_axis_3d);
        CHECK(reflect_on_vec(x_axis_3d, y_axis_3d) == -x_axis_3d);
        CHECK(reflect_on_vec(z_axis_3d, y_axis_3d) == -z_axis_3d);
        CHECK(reflect_on_vec(x_axis_3d, z_axis_3d) == -x_axis_3d);
        CHECK(reflect_on_vec(y_axis_3d, z_axis_3d) == -y_axis_3d);

        // coordinate axis reflected on itself remains itself (identity)
        CHECK(reflect_on_vec(x_axis_3d, x_axis_3d) == x_axis_3d);
        CHECK(reflect_on_vec(y_axis_3d, y_axis_3d) == y_axis_3d);
        CHECK(reflect_on_vec(z_axis_3d, z_axis_3d) == z_axis_3d);

        // point reflected on a hyperplane that the vector is a normal to
        // the hyperplane can be created by taking the dual (or the rcmpl) of the normal
        CHECK(reflect_on(p, dual(y_axis_3d)) == vec3d{4, -1, 0});
        // alternatively the plane can be used directly (represented by a bivector)
        CHECK(reflect_on(p, zx_plane_3d) == vec3d{4, -1, 0});

        // coordinate axis reflected on perpendicular base planes yield their negatives
        CHECK(reflect_on(x_axis_3d, yz_plane_3d) == -x_axis_3d);
        CHECK(reflect_on(x_axis_3d, dual(x_axis_3d)) == -x_axis_3d);
        CHECK(reflect_on(y_axis_3d, zx_plane_3d) == -y_axis_3d);
        CHECK(reflect_on(y_axis_3d, dual(y_axis_3d)) == -y_axis_3d);
        CHECK(reflect_on(z_axis_3d, xy_plane_3d) == -z_axis_3d);
        CHECK(reflect_on(z_axis_3d, dual(z_axis_3d)) == -z_axis_3d);

        // a coordinate plane reflected on itself remains itself (identity)
        CHECK(reflect_on(yz_plane_3d, yz_plane_3d) == yz_plane_3d);
        CHECK(reflect_on(zx_plane_3d, zx_plane_3d) == zx_plane_3d);
        CHECK(reflect_on(xy_plane_3d, xy_plane_3d) == xy_plane_3d);

        // reflect planes on planes directly
        CHECK(reflect_on(e23_3d + e12_3d, e12_3d) == -e23_3d + e12_3d);
    }

    TEST_CASE("Vec3d: operations - project / reject / reflect (vector - bivector)")
    {
        fmt::println(
            "Vec3d: operations - project / reject / reflect (vector - bivector)");

        vec3d v1{5.0, 3.0, 1.0};
        bivec3d v2 = wdg(vec3d{0.0, 0.0, 2.0}, vec3d{2.0, 0.0, 0.0});

        vec3d v3{project_onto(v1, v2)};
        vec3d v4{reject_from(v1, v2)};
        vec3d v5{v3 + v4};

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .4f}, nrm(v3) = {: .4f}", v3,
        //              nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .4f}, nrm(v4) = {: .4f}", v4,
        //              nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .4f}, nrm(v5) = {: .4f}", v5,
        //              nrm(v5));

        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);
        CHECK(project_onto(v1, e31_3d) == vec3d{5, 0, 1});
        CHECK(reject_from(v1, e31_3d) == vec3d{0, 3, 0});

        // from Macdonald, "Linear and Geometric Algebra", Exercise 7.14, p. 129
        CHECK(reflect_on(wdg(e1_3d + e3_3d, e2_3d), e12_3d) == wdg(e1_3d - e3_3d, e2_3d));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3d<T> Gram-Schmidt-Orthogonalization
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3d: Gram-Schmidt-Orthogonalization")
    {
        fmt::println("Vec3d: Gram-Schmidt-Orthogonalization - part 1 (plane)");

        vec3d v1{2.0, 2.0, 0.0};
        vec3d v2{-1.5, 2.5, 1.0};
        vec3d v3{1.5, -2.5, -1.0};

        std::vector<vec3d> og12 = gs_orthogonal(v1, v2);
        std::vector<vec3d> og13 = gs_orthogonal(v1, v3);

        std::vector<vec3d> on12 = gs_orthonormal(v1, v2);
        std::vector<vec3d> on13 = gs_orthonormal(v1, v3);

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}, angle = {: .2f}°", v1, nrm(v1),
        //              rad2deg(angle(e1_3d, v1)));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v2, nrm(v2),
        //              rad2deg(angle(e1_3d, v2)));
        // fmt::println("v3  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v3, nrm(v3),
        //              rad2deg(angle(e1_3d, v3)));
        // fmt::println("");
        // fmt::println("og12[0]  = {: .4f}, nrm(og12[0]) = {: .4f}, angle = {: .2f}°",
        //              og12[0], nrm(og12[0]), rad2deg(angle(e1_3d, og12[0])));
        // fmt::println("og12[1]  = {: .4f}, nrm(og12[1]) = {: .4f}, angle = {: .2f}°",
        //              og12[1], nrm(og12[1]), rad2deg(angle(e1_3d, og12[1])));
        // fmt::println("on12[0]  = {: .4f}, nrm(on12[0]) = {: .4f}, angle = {: .2f}°",
        //              on12[0], nrm(on12[0]), rad2deg(angle(e1_3d, on12[0])));
        // fmt::println("on12[1]  = {: .4f}, nrm(on12[1]) = {: .4f}, angle = {: .2f}°",
        //              on12[1], nrm(on12[1]), rad2deg(angle(e1_3d, on12[1])));
        // fmt::println("");
        // fmt::println("og13[0]  = {: .4f}, nrm(og13[0]) = {: .4f}, angle = {: .2f}°",
        //              og13[0], nrm(og13[0]), rad2deg(angle(e1_3d, og13[0])));
        // fmt::println("og13[1]  = {: .4f}, nrm(og13[1]) = {: .4f}, angle = {: .2f}°",
        //              og13[1], nrm(og13[1]), rad2deg(angle(e1_3d, og13[1])));
        // fmt::println("on13[0]  = {: .4f}, nrm(on13[0]) = {: .4f}, angle = {: .2f}°",
        //              on13[0], nrm(on13[0]), rad2deg(angle(e1_3d, on13[0])));
        // fmt::println("on13[1]  = {: .4f}, nrm(on13[1]) = {: .4f}, angle = {: .2f}°",
        //              on13[1], nrm(on13[1]), rad2deg(angle(e1_3d, on13[1])));
        // fmt::println("");

        CHECK(std::abs(dot(og12[0], og12[1])) < eps);
        CHECK(std::abs(dot(og13[0], og13[1])) < eps);

        CHECK(std::abs(dot(on12[0], on12[1])) < eps);
        CHECK(std::abs(nrm(on12[0]) - 1.0) < eps);
        CHECK(std::abs(nrm(on12[1]) - 1.0) < eps);

        CHECK(std::abs(dot(on13[0], on13[1])) < eps);
        CHECK(std::abs(nrm(on13[0]) - 1.0) < eps);
        CHECK(std::abs(nrm(on13[1]) - 1.0) < eps);

        fmt::println("Vec3d: Gram-Schmidt-Orthogonalization - part 2 (space)");

        v1 = vec3d{2.0, 2.0, 0.0};
        v2 = vec3d{-1.5, 2.5, 1.0};
        v3 = vec3d{1.5, -2.5, -6.0}; // all three vectors must be linear independent

        std::vector<vec3d> og123 = gs_orthogonal(v1, v2, v3);
        std::vector<vec3d> on123 = gs_orthonormal(v1, v2, v3);

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}, angle = {: .2f}°", v1, nrm(v1),
        //              rad2deg(angle(e1_3d, v1)));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v2, nrm(v2),
        //              rad2deg(angle(e1_3d, v2)));
        // fmt::println("v3  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v3, nrm(v3),
        //              rad2deg(angle(e1_3d, v3)));
        // fmt::println("");
        // fmt::println("og123[0]  = {: .4f}, nrm(og123[0]) = {: .4f}, angle = {: .2f}°",
        //              og123[0], nrm(og123[0]), rad2deg(angle(e1_3d, og123[0])));
        // fmt::println("og123[1]  = {: .4f}, nrm(og123[1]) = {: .4f}, angle = {: .2f}°",
        //              og123[1], nrm(og123[1]), rad2deg(angle(e1_3d, og123[1])));
        // fmt::println("og123[2]  = {: .4f}, nrm(og123[2]) = {: .4f}, angle = {: .2f}°",
        //              og123[2], nrm(og123[2]), rad2deg(angle(e1_3d, og123[2])));
        // fmt::println("");
        // fmt::println("on123[0]  = {: .4f}, nrm(on123[0]) = {: .4f}, angle = {: .2f}°",
        //              on123[0], nrm(on123[0]), rad2deg(angle(e1_3d, on123[0])));
        // fmt::println("on123[1]  = {: .4f}, nrm(on123[1]) = {: .4f}, angle = {: .2f}°",
        //              on123[1], nrm(on123[1]), rad2deg(angle(e1_3d, on123[1])));
        // fmt::println("on123[2]  = {: .4f}, nrm(on123[2]) = {: .4f}, angle = {: .2f}°",
        //              on123[2], nrm(on123[2]), rad2deg(angle(e1_3d, on123[2])));
        // fmt::println("");
        // fmt::println("sign(wdg(on123[0], wdg(on123[1], on123[2]))/I_3d) = {} "
        //              "(+1: right-handed, -1: left-handed)",
        //              sign(wdg(on123[0], wdg(on123[1], on123[2])) / I_3d));
        // fmt::println("");

        CHECK(std::abs(dot(og123[0], og123[1])) < eps);
        CHECK(std::abs(dot(og123[1], og123[2])) < eps);
        CHECK(std::abs(dot(og123[0], og123[2])) < eps);

        CHECK(std::abs(dot(on123[0], on123[1])) < eps);
        CHECK(std::abs(dot(on123[1], on123[2])) < eps);
        CHECK(std::abs(dot(on123[0], on123[2])) < eps);

        CHECK(std::abs(nrm(on123[0]) - 1.0) < eps);
        CHECK(std::abs(nrm(on123[1]) - 1.0) < eps);
        CHECK(std::abs(nrm(on123[2]) - 1.0) < eps);
    }

    TEST_CASE("Vec3d: cross-product")
    {
        fmt::println("Vec3d: cross-product");

        vec3d u{1.0, 1.0, 0.0};
        vec3d v{0.0, 1.0, 1.0};
        vec3d w{1.0, 1.0, 1.0};

        vec3d u_cross_v = cross(u, v);
        bivec3d u_wdg_v = wdg(u, v);

        CHECK(u_cross_v == dual(u_wdg_v));
        CHECK(u_wdg_v == dual(u_cross_v));

        // definitions using the pseudoscalars directly are valid independent of the
        // duality definitions

        CHECK(u_cross_v == u_wdg_v * rev(I_3d));
        CHECK(u_cross_v == rev(u_wdg_v) * I_3d);
        CHECK(u_cross_v == -I_3d * u_wdg_v);
        CHECK(u_cross_v == -u_wdg_v * I_3d);

        CHECK(wdg(u, v) == u_cross_v * I_3d);
        CHECK(wdg(u, v) == I_3d * u_cross_v);
        CHECK(wdg(u, v) == -u_cross_v * rev(I_3d));

        // double cross product identity
        CHECK(cross(u, cross(v, w)) == (u << wdg(v, w)));
        CHECK(cross(u, cross(v, w)) == (dot(u, w) * v - dot(u, v) * w));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec3d<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3d: default init")
    {
        fmt::println("MVec3d: default init");
        // default initialization
        mvec3d v;
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.c0) < eps);
        CHECK(std::abs(v.c1) < eps);
        CHECK(std::abs(v.c2) < eps);
        CHECK(std::abs(v.c3) < eps);
        CHECK(std::abs(v.c4) < eps);
        CHECK(std::abs(v.c5) < eps);
        CHECK(std::abs(v.c6) < eps);
        CHECK(std::abs(v.c7) < eps);
    }
    TEST_CASE("MVec3d: with curly braced intializer")
    {
        fmt::println("MVec3d: with curly braced intializer");
        // default initialization
        mvec3d v{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.c0 - 0.0) < eps);
        CHECK(std::abs(v.c1 - 1.0) < eps);
        CHECK(std::abs(v.c2 - 2.0) < eps);
        CHECK(std::abs(v.c3 - 3.0) < eps);
        CHECK(std::abs(v.c4 - 23.0) < eps);
        CHECK(std::abs(v.c5 - 31.0) < eps);
        CHECK(std::abs(v.c6 - 12.0) < eps);
        CHECK(std::abs(v.c7 - 123.0) < eps);
    }

    TEST_CASE("MVec3d: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("MVec3d: cp ctor & cp assign incl. type deduction");
        // default initialization
        mvec3d v1{0.0,  1.0,  2.0,  3.0,
                  23.0, 31.0, 12.0, 123.0}; // init with double (type deduction)
        mvec3d v2{v1};                      // cp ctor
        mvec3d v3 = v2;                     // cp assign
        mvec3d v4 = -v3;                    // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(std::abs(v2.c0 - 0.0) < eps);
        CHECK(std::abs(v2.c1 - 1.0) < eps);
        CHECK(std::abs(v2.c2 - 2.0) < eps);
        CHECK(std::abs(v2.c3 - 3.0) < eps);
        CHECK(std::abs(v2.c4 - 23.0) < eps);
        CHECK(std::abs(v2.c5 - 31.0) < eps);
        CHECK(std::abs(v2.c6 - 12.0) < eps);
        CHECK(std::abs(v2.c7 - 123.0) < eps);

        CHECK(std::abs(v3.c0 - 0.0) < eps);
        CHECK(std::abs(v3.c1 - 1.0) < eps);
        CHECK(std::abs(v3.c2 - 2.0) < eps);
        CHECK(std::abs(v3.c3 - 3.0) < eps);
        CHECK(std::abs(v3.c4 - 23.0) < eps);
        CHECK(std::abs(v3.c5 - 31.0) < eps);
        CHECK(std::abs(v3.c6 - 12.0) < eps);
        CHECK(std::abs(v3.c7 - 123.0) < eps);

        CHECK(v4 == -v3);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }

    TEST_CASE("MVec3d: defining basic types and ctor checks")
    {
        fmt::println("MVec3d: defining basic types and ctor checks");

        auto s = scalar3d(5.0);
        auto v = vec3d{1.0, 2.0, 1.0};
        auto B = bivec3d{-1.0, 2.0, 1.0};
        auto ps = pscalar3d(-5.0);

        auto mv0 = mvec3d{s};
        auto mv1 = mvec3d{v};
        auto mv2 = mvec3d{B};
        auto mv3 = mvec3d{ps};

        auto mv5a = mvec3d{s, B};
        auto mv5b = mvec3d_e{s};
        auto mv5c = mvec3d_e{B};
        auto mv5d = mvec3d_e{s, B};
        auto mv5e = mvec3d{mv5d};

        auto mv6a = mvec3d{v, ps};
        auto mv6b = mvec3d_u{v};
        auto mv6c = mvec3d_u{ps};
        auto mv6d = mvec3d_u{v, ps};
        auto mv6e = mvec3d{mv6d};

        auto mv7 = mvec3d{s, v, B, ps};

        // fmt::println("   mv0  = {}", mv0);
        // fmt::println("   mv1  = {}", mv1);
        // fmt::println("   mv2  = {}", mv2);
        // fmt::println("   mv3  = {}", mv3);
        // fmt::println("   mv5a = {}", mv5a);
        // fmt::println("   mv5b = {}", mv5b);
        // fmt::println("   mv5c = {}", mv5c);
        // fmt::println("   mv5d = {}", mv5d);
        // fmt::println("   mv5e = {}", mv5e);
        // fmt::println("   mv6a = {}", mv6a);
        // fmt::println("   mv6b = {}", mv6b);
        // fmt::println("   mv6c = {}", mv6c);
        // fmt::println("   mv6d = {}", mv6d);
        // fmt::println("   mv6e = {}", mv6e);
        // fmt::println("   mv7  = {}", mv7);

        CHECK(gr0(mv0) == s);
        CHECK(gr1(mv1) == v);
        CHECK(gr2(mv2) == B);
        CHECK(gr3(mv3) == ps);

        CHECK(gr0(mv5a) == s);
        CHECK(gr1(mv5a) == vec3d{});
        CHECK(gr2(mv5a) == B);
        CHECK(gr3(mv5a) == pscalar3d{});

        CHECK(gr0(mv5b) == s);
        CHECK(gr2(mv5b) == bivec3d{});

        CHECK(gr0(mv5c) == scalar3d{});
        CHECK(gr2(mv5c) == B);

        CHECK(gr0(mv5d) == s);
        CHECK(gr2(mv5d) == B);

        CHECK(gr0(mv5e) == s);
        CHECK(gr1(mv5e) == vec3d{});
        CHECK(gr2(mv5e) == B);
        CHECK(gr3(mv5e) == pscalar3d{});

        CHECK(gr0(mv6a) == scalar3d{});
        CHECK(gr1(mv6a) == v);
        CHECK(gr2(mv6a) == bivec3d{});
        CHECK(gr3(mv6a) == ps);

        CHECK(gr1(mv6b) == v);
        CHECK(gr3(mv6b) == pscalar3d{});

        CHECK(gr1(mv6c) == vec3d{});
        CHECK(gr3(mv6c) == ps);

        CHECK(gr1(mv6d) == v);
        CHECK(gr3(mv6d) == ps);

        CHECK(gr0(mv6e) == scalar3d{});
        CHECK(gr1(mv6e) == v);
        CHECK(gr2(mv6e) == bivec3d{});
        CHECK(gr3(mv6e) == ps);

        CHECK(gr0(mv7) == s);
        CHECK(gr1(mv7) == v);
        CHECK(gr2(mv7) == B);
        CHECK(gr3(mv7) == ps);

        CHECK(gr(s) == 0);
        CHECK(gr(v) == 1);
        CHECK(gr(B) == 2);
        CHECK(gr(ps) == 3);
    }

    TEST_CASE("MVec3d: fmt & cout printing")
    {
        fmt::println("MVec3d: fmt & cout printing:\n");

        mvec3d pf{1.0f, 2.00001f, 0.0f, 3.0f, 1.0f, 2.00001f, 0.0f, 3.0f};
        mvec3d pd{1.0, 2.00001, 0.0, 3.0, 1.0, 2.00001, 0.0, 3.0};

        std::cout << "    cout: pf = " << pf << std::endl;
        fmt::println("    fmt:  pf = {}", pf);
        fmt::println("    fmt:  pf = {:.4f}", pf);
        fmt::println("");

        std::cout << "    cout: pd = " << pd << std::endl;
        fmt::println("    fmt:  pd = {}", pd);
        fmt::println("    fmt:  pd = {:.4f}", pd);
        fmt::println("");

        std::vector<MVec3d<double>> vp1{{1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0},
                                        {0.5, 1.5, 2.0, 2.5, 1.0, 1.0, 1.0, 2.0}};
        fmt::println("    fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("    fmt: vp1 = {:g}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("MVec3d: vector space and linearity tests")
    {
        fmt::println("MVec3d: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        mvec3d p0;
        mvec3d p1{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
        mvec3d p2{0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0};
        mvec3d p3{0.0, 3.0, 6.0, 9.0, 0.0, 3.0, 6.0, 9.0};
        mvec3d p4 = -p1; // assignment using unary minus
        double s = 2.35;
        double t = -1.3;

        CHECK(p1 + p1 == p2); // addition is defined

        // vector addition
        CHECK(p2 + p1 == p1 + p2);               // addition is commutative
        CHECK((p1 + p2) + p3 == p1 + (p2 + p3)); // addition is associative
        CHECK(p1 + p0 == p1);                    // zero is the additive identity
        CHECK(p1 * 0.0 == p0); // scalar multplication with null creates the null vector

        // scalar multiplication
        CHECK(p1 * 1.0 == p1);                   // 1.0 is the multiplicative identity
        CHECK((s * t) * p1 == s * (t * p1));     // is associative w.r.t.multiplication
        CHECK(s * (p1 + p2) == s * p1 + s * p2); // scalar multiplication distributes
        CHECK((p1 + p2) * s == p1 * s + p2 * s); // over vector addition
        CHECK((s + t) * p1 == s * p1 + t * p1);  // and is associative w.r.t. addition

        // additional tests
        CHECK(p1 + (-p1) == p0); // there is an inverse element with respect to addition
        CHECK(p1 + p2 == p3);    // component wise addition
        CHECK(p1 * 2.0 == p2);   // component wise multiplication
        CHECK(p4 == -p1);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec3d<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3d: wedge product - basic properties")
    {
        fmt::println("MVec3d: wedge product - basic properties");

        vec3d v1{1.0, 2.0, 1.0};
        vec3d v2{0.5, 3.0, 2.0};
        vec3d v3{-2.0, 6.0, 3.0};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3d{sd};
        auto t = scalar3d{st};

        CHECK(wdg(v1, v1) == bivec3d{});                      // wdg=0 for collin. vectors
        CHECK(wdg(v1, v2) == -wdg(v2, v1));                   // anticommutative for vect.
        CHECK(wdg(wdg(v1, v2), v3) == wdg(v1, wdg(v2, v3)));  // wdg is associative
        CHECK(wdg(v1, v2 + v3) == wdg(v1, v2) + wdg(v1, v3)); // wdg distributes over add.
        CHECK(wdg(v1 + v2, v3) == wdg(v1, v3) + wdg(v2, v3)); // wdg distributes over add.
        CHECK(wdg(sd * v1, v2) == wdg(v1, sd * v2)); // scalars can be factored out of wdg
        CHECK(wdg(sd * v1, v2) == sd * wdg(v1, v2)); // scalars can be factored out of wdg
        CHECK(wdg(s, t) == wdg(t, s));   // wdg between scalars equivalent to scalar mult.
        CHECK(wdg(s, v1) == wdg(v1, s)); // wdg between scalar and vector
        CHECK(wdg(s, v1) == sd * v1);    // wdg between scalar and vector
    }

    TEST_CASE("MVec3d: geometric product - basic properties")
    {
        fmt::println("MVec3d: geometric product - basic properties");

        vec3d v1{1.0, 2.0, 1.0};
        vec3d v2{0.5, 3.0, 2.0};
        vec3d v3{-2.0, 6.0, 3.0};
        bivec3d B1{-4, 2, 1};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3d{sd};
        auto t = scalar3d{st};

        // nrm_sq(v1) = dot(v1,v1)  for every vector
        CHECK(nrm_sq(v1) == dot(v1, v1));

        CHECK(v1 * v2 == dot(v1, v2) + wdg(v1, v2)); // valid for vectors only
        CHECK(v1 * v2 == (v2 >> v1) + wdg(v1, v2));  // contraction = dot for same grades

        // valid also for gr(B1) != gr(v1)
        CHECK(v1 * B1 == (B1 >> v1) + wdg(v1, B1));
        CHECK(B1 * v1 == (v1 << B1) + wdg(B1, v1));

        // only valid for vectors (!):
        // dot = gpr symmetric part
        CHECK(dot(v1, v2) == gr0(0.5 * (v1 * v2 + v2 * v1)));
        // wdg = gpr anti-symmetric part
        CHECK(wdg(v1, v2) == gr2(0.5 * (v1 * v2 - v2 * v1)));

        // mathematical characteristics
        CHECK((v1 * v2) * v3 == v1 * (v2 * v3));    // gpr is associative
        CHECK(v1 * (v2 + v3) == v1 * v2 + v1 * v3); // gpr distributes over addition
        CHECK((v1 + v2) * v3 == v1 * v3 + v2 * v3); // wdg distributes over addition
        CHECK((sd * v1) * v2 == v1 * (sd * v2));    // scalars can be factored out of gpr
        CHECK((sd * v1) * v2 == sd * (v1 * v2));    // scalars can be factored out of gpr
        CHECK(s * t == t * s);    // gpr between scalars equivalent to scalar mult.
        CHECK(s * v1 == v1 * s);  // gpr between scalar and vector
        CHECK(s * v1 == sd * v1); // gpr between scalar and vector
    }

    TEST_CASE("MVec3d: geometric product tests - vec * vec")
    {
        fmt::println("MVec3d: geometric product tests - vec * vec");

        // ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)
        //
        // dot(a,b) = 0.5*(ab + ba)   (symmetric part)
        // wdg(a,b) = 0.5*(ab - ba)   (antisymmetric part)

        vec3d a{1.0, 2.0, 3.0};
        vec3d b{0.5, 3.0, -2.0};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        mvec3d mva{a};
        mvec3d mvb{b};
        auto mvab = mva * mvb;
        auto mvab_sym = 0.5 * (mva * mvb + mvb * mva);
        auto mvab_asym = 0.5 * (mva * mvb - mvb * mva);

        // fmt::println("   a = {}", a);
        // fmt::println("   b = {}", b);
        // fmt::println("   dot(a,b) = {}", dot_ab);
        // fmt::println("   wdg(a,b) = {}", wdg_ab);
        // fmt::println("");
        // fmt::println("   mva  = {}", mva);
        // fmt::println("   mvb  = {}", mvb);
        // fmt::println("   mvab = {}", mvab);
        // fmt::println("   mvab_sym  = 0.5*(mva * mvb + mvb * mva) = {}", mvab_sym);
        // fmt::println("   mvab_asym = 0.5*(mva * mvb - mvb * mva) = {}", mvab_asym);
        // fmt::println("");
        // fmt::println("   gr0(mvab) = {}", gr0(mvab));
        // fmt::println("   gr1(mvab) = {}", gr1(mvab));
        // fmt::println("   gr2(mvab) = {}", gr2(mvab));
        // fmt::println("   gr3(mvab) = {}", gr3(mvab));

        CHECK(dot_ab == gr0(mvab));
        CHECK(dot_ab == gr0(mvab_sym));
        CHECK(wdg_ab == gr2(mvab));
        CHECK(wdg_ab == gr2(mvab_asym));
    }

    TEST_CASE("MVec3d: geometric product tests - bivec * vec")
    {
        fmt::println("MVec3d: geometric product tests - bivec * vec");

        // Ab = (b << A) + wdg(A,b) = gr1(Ab) + gr3(Ab)
        //
        // (b << A) = 0.5*(Ab - Aa)   (antisymmetric part)
        // wdg(A,b) = 0.5*(Ab + Aa)   (symmetric part)

        bivec3d A{1.0, 2.0, 3.0};
        vec3d b{0.5, 3.0, -2.0};
        // auto dot_ab = dot(a, b);
        auto b_onto_A = (b << A);
        auto wdg_Ab = wdg(A, b);

        mvec3d mva{A};
        mvec3d mvb{b};
        auto mvab = mva * mvb;
        auto mvab_sym = 0.5 * (mva * mvb + mvb * mva);
        auto mvab_asym = 0.5 * (mva * mvb - mvb * mva);

        // fmt::println("");
        // fmt::println("   A = {}", A);
        // fmt::println("   b = {}", b);
        // fmt::println("   (b << A) = {}, gr1(A*b) = {}", b_onto_A, gr1(A * b));
        // fmt::println("   wdg(A,b) = {}, gr3(A*b) = {}", wdg_Ab, gr3(A * b));
        // fmt::println("");
        // fmt::println("   mva  = {}", mva);
        // fmt::println("   mvb  = {}", mvb);
        // fmt::println("   mvab = {}", mvab);
        // fmt::println("   mvab_sym  = 0.5*(mva * mvb + mvb * mva) = {}", mvab_sym);
        // fmt::println("   mvab_asym = 0.5*(mva * mvb - mvb * mva) = {}", mvab_asym);
        // fmt::println("");
        // fmt::println("   gr0(mvab) = {}", gr0(mvab));
        // fmt::println("   gr1(mvab) = {}", gr1(mvab));
        // fmt::println("   gr2(mvab) = {}", gr2(mvab));
        // fmt::println("   gr3(mvab) = {}", gr3(mvab));
        // fmt::println("");

        CHECK(b_onto_A == gr1(mvab));
        CHECK(b_onto_A == gr1(mvab_asym));
        CHECK(wdg_Ab == gr3(mvab));
        CHECK(wdg_Ab == gr3(mvab_sym));
    }

    TEST_CASE("MVec3d: geometric product tests - vec * bivec")
    {
        fmt::println("MVec3d: geometric product tests - vec * bivec");

        // a*B = (B >> a) + wdg(a,B) = gr1(aB) + gr3(aB)
        //
        // (B >> a) = 0.5*(aB - Ba)   (antisymmetric part)
        // wdg(a,B) = 0.5*(aB + Ba)   (symmetric part)

        vec3d a{1.0, 2.0, 3.0};
        bivec3d B{0.5, 3.0, -2.0};
        // auto dot_ab = dot(a, b);
        auto B_by_a = (B >> a);
        auto wdg_aB = wdg(a, B);

        mvec3d mva{a};
        mvec3d mvb{B};
        auto mvab = mva * mvb;
        auto mvab_sym = 0.5 * (mva * mvb + mvb * mva);
        auto mvab_asym = 0.5 * (mva * mvb - mvb * mva);

        // fmt::println("");
        // fmt::println("   a = {}", a);
        // fmt::println("   B = {}", B);
        // fmt::println("   (B >> a) = {}", B_by_a);
        // fmt::println("   wdg(a,b) = {}", wdg_aB);
        // fmt::println("");
        // fmt::println("   mva  = {}", mva);
        // fmt::println("   mvb  = {}", mvb);
        // fmt::println("   mvab = {}", mvab);
        // fmt::println("   mvab_sym  = 0.5*(mva * mvb + mvb * mva) = {}", mvab_sym);
        // fmt::println("   mvab_asym = 0.5*(mva * mvb - mvb * mva) = {}", mvab_asym);
        // fmt::println("");
        // fmt::println("   gr0(mvab) = {}", gr0(mvab));
        // fmt::println("   gr1(mvab) = {}", gr1(mvab));
        // fmt::println("   gr2(mvab) = {}", gr2(mvab));
        // fmt::println("   gr3(mvab) = {}", gr3(mvab));
        // fmt::println("");

        CHECK(B_by_a == gr1(mvab));
        CHECK(B_by_a == gr1(mvab_asym));
        CHECK(wdg_aB == gr3(mvab));
        CHECK(wdg_aB == gr3(mvab_sym));
    }

    TEST_CASE("MVec3d: geometric product - combinatorial tests")
    {
        fmt::println("MVec3d: geometric product - combinatorial tests");

        scalar3d s1{3.0};
        scalar3d s2{-1.5};
        vec3d v1{1.0, 2.0, 3.0};
        vec3d v2{0.5, 3.0, -2.0};
        bivec3d b1{-1.0, -3.0, -5.0};
        bivec3d b2{2.0, 0.5, -4.0};
        pscalar3d ps1{-5.0};
        pscalar3d ps2{2.0};

        // checks all ctor combinations and equivalence of simplified products with
        // geometric product of a fully populated multivector
        //
        // assumption: the geometric product of a fully populated multivector serves as
        //             reference and is correctly implemented
        //             results are relative to full geometric product

        // MVec3d_E * MVec3d
        CHECK(mvec3d{s1, b1} * mvec3d{s1, v1, b1, ps1} ==
              mvec3d{mvec3d_e{s1, b1} * mvec3d{s1, v1, b1, ps1}});

        // MVec3d * MVec3d_E
        CHECK(mvec3d{s1, v1, b1, ps1} * mvec3d{s1, b1} ==
              mvec3d{mvec3d{s1, v1, b1, ps1} * mvec3d_e{s1, b1}});

        // MVec3d_E * MVec3d_E
        CHECK(mvec3d{s1, b1} * mvec3d{s2, b2} ==
              mvec3d{mvec3d_e{s1, b1} * mvec3d_e{s2, b2}});

        // MVec3d_U * MVec3d_U
        CHECK(mvec3d{v1, ps1} * mvec3d{v2, ps2} ==
              mvec3d{mvec3d_u{v1, ps1} * mvec3d_u{v2, ps2}});

        // MVec3d_E * MVec3d_U
        CHECK(mvec3d{s1, b1} * mvec3d{v2, ps2} ==
              mvec3d{mvec3d_e{s1, b1} * mvec3d_u{v2, ps2}});

        // MVec3d_U * MVec3d_E
        CHECK(mvec3d{v1, ps1} * mvec3d{s2, b2} ==
              mvec3d{mvec3d_u{v1, ps1} * mvec3d_e{s2, b2}});

        // MVec3d_E * BiVec3d
        CHECK(mvec3d{s1, b1} * mvec3d{b2} == mvec3d{mvec3d_e{s1, b1} * mvec3d{b2}});

        // BiVec3d * MVec3d_E
        CHECK(mvec3d{b1} * mvec3d{s2, b2} == mvec3d{mvec3d{b1} * mvec3d_e{s2, b2}});

        // MVec3d_E * Vec3d
        CHECK(mvec3d{s1, b1} * mvec3d{v2} == mvec3d{mvec3d_e{s1, b1} * mvec3d{v2}});

        // Vec3d * MVec3d_E
        CHECK(mvec3d{v1} * mvec3d{s2, b2} == mvec3d{mvec3d{v1} * mvec3d_e{s2, b2}});

        // MVec3d_U * BiVec3d
        CHECK(mvec3d{v1, ps1} * mvec3d{b2} == mvec3d{mvec3d_u{v1, ps1} * mvec3d{b2}});

        // BiVec3d * MVec3d_U
        CHECK(mvec3d{b1} * mvec3d{v2, ps2} == mvec3d{mvec3d{b1} * mvec3d_u{v2, ps2}});

        // BiVec3d * BiVec3d
        CHECK(mvec3d{b1} * mvec3d{b2} == mvec3d{b1 * b2});

        // BiVec3d * Vec3d
        CHECK(mvec3d{b1} * mvec3d{v1} == mvec3d{b1 * v1});

        // Vec3d * BiVec3d
        CHECK(mvec3d{v1} * mvec3d{b1} == mvec3d{v1 * b1});

        // Vec3d * Vec3d
        CHECK(mvec3d{v1} * mvec3d{v2} == mvec3d{v1 * v2});

        // PScalar3d * MVec3d
        CHECK(mvec3d{ps1} * mvec3d{s2, v2, b2, ps2} ==
              mvec3d{ps1 * mvec3d{s2, v2, b2, ps2}});

        // MVec3d * Pscalar3d
        CHECK(mvec3d{s1, v1, b1, ps1} * mvec3d{ps2} ==
              mvec3d{mvec3d{s1, v1, b1, ps1} * ps2});

        // PScalar3d * MVec3d_E
        CHECK(mvec3d{ps1} * mvec3d{s2, b2} == mvec3d{ps1 * mvec3d_e{s2, b2}});

        // MVec3d_E * Pscalar3d
        CHECK(mvec3d{s1, b1} * mvec3d{ps2} == mvec3d{mvec3d_e{s1, b1} * ps2});

        // PScalar3d * MVec3d_U
        CHECK(mvec3d{ps1} * mvec3d{v2, ps2} == mvec3d{ps1 * mvec3d_u{v2, ps2}});

        // MVec3d_U * Pscalar3d
        CHECK(mvec3d{v1, ps1} * mvec3d{ps2} == mvec3d{mvec3d_u{v1, ps1} * ps2});

        // PScalar3d * BiVec3d
        CHECK(mvec3d{ps1} * mvec3d{b2} == mvec3d{ps1 * b2});

        // BiVec3d * Pscalar3d
        CHECK(mvec3d{b1} * mvec3d{ps2} == mvec3d{b1 * ps2});

        // PScalar3d * Vec3d
        CHECK(mvec3d{ps1} * mvec3d{v2} == mvec3d{ps1 * v2});

        // Vec3d * Pscalar3d
        CHECK(mvec3d{v1} * mvec3d{ps2} == mvec3d{v1 * ps2});

        // Pscalar3d * Pscalar3d
        CHECK(mvec3d{ps1} * mvec3d{ps2} == mvec3d{ps1 * ps2});

        // Scalar * Scalar
        CHECK(mvec3d{s1} * mvec3d{s2} == mvec3d{s1 * s2});
    }

    TEST_CASE("MVec3d: geometric product tests - recovering vectors from the"
              "geometric product")
    {
        fmt::println("MVec3d: geometric product tests - recovering vectors from the "
                     "geometric product");

        // Two multivectors mv1 and mv2 formed from vectors v1 and v2.
        // (gr0(mv1)==0 && gr1(mv1) != 0 && gr2(mv1)==0 &&
        //  gr0(mv2)==0 && gr1(mv2) != 0 && gr2(mv2)==0 )
        //
        // They are multiplied by the geometric product to form a multivector C
        // C = mv1(v1) * mv2(v2) = mv1 * mv2
        //
        // C contains a scalar part and a bivector part exclusively,
        // the remaining components are zero.
        // (gr0(C) != 0 && gr1(C)==0 && gr2(C) !=0)
        //
        // The scalar part of C represents the parts of v1 and v2
        // that are parallel to each other.
        // The bivector part of C represents the parts of v1 and v2
        // that are perpendicular to each other.
        //
        // multiply C from the right with inv(v2) recovers v1
        // multiply C from the left the the inv(v1) recovers v2

        vec3d a{1.0, 2.0, 3.0};
        vec3d b{0.5, 3.0, -4.0};
        mvec3d mva{a};
        mvec3d mvb{b};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);
        mvec3d C = a * b;
        mvec3d Cm = mva * mvb;
        mvec3d Cd{scalar3d(dot_ab), wdg_ab};

        mvec3d gpr_right = C * mvec3d{inv(b)};
        mvec3d gpr_left = mvec3d{inv(a)} * C;

        // fmt::println("   a                           = {}", a);
        // fmt::println("   b                           = {}", b);
        // fmt::println("   C  = a * b                  = {}", C);
        // fmt::println("   Cm = mva * mvb              = {}", Cm);
        // fmt::println("   Cd = mv(dot(a,b), wdg(a,b)) = {}", Cd);
        // fmt::println("");
        // fmt::println("   C * bi = gpr_right = {}", gpr_right);
        // fmt::println("   ai * C = gpr_left  = {}", gpr_left);
        // fmt::println("   gr1(gpr_right) = a = {}", gr1(gpr_right));
        // fmt::println("   gr1(gpr_left)  = b = {}", gr1(gpr_left));

        CHECK(C == Cm);
        CHECK(C == Cd);
        CHECK(a == gr1(gpr_right));
        CHECK(b == gr1(gpr_left));

        auto m = mvec2d{13.0, -27.0, 3.0, 5.0};
        auto prd = m * inv(m);
        CHECK(std::abs(gr0(prd) - 1.0) < eps);
        CHECK(nrm(gr1(prd)) < eps);
        CHECK(std::abs(gr2(prd) - 0.0) < eps);
    }

    TEST_CASE("MVec3d: geometric product tests - equivalence tests")
    {
        fmt::println("MVec3d: geometric product tests - equivalence tests");

        vec3d a{1.0, 2.0, 3.0};
        vec3d b{0.5, 3.0, -4.0};
        mvec3d mva{a};
        mvec3d mvb{b};

        bivec3d A{1.0, 2.0, 3.0};
        bivec3d B{0.5, 3.0, -4.0};
        mvec3d mvA{A};
        mvec3d mvB{B};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        auto b_onto_A = (b << A);
        auto wdg_Ab = wdg(A, b);

        auto B_by_a = (B >> a);
        auto wdg_aB = wdg(a, B);

        mvec3d_e ab = a * b;
        mvec3d abm = mva * mvb;
        mvec3d abd{scalar3d{dot_ab}, wdg_ab};

        mvec3d_u Ab = A * b;
        mvec3d Abm = mvA * mvb;
        mvec3d Abd{b_onto_A, wdg_Ab};

        mvec3d_u aB = a * B;
        mvec3d aBm = mva * mvB;
        mvec3d aBd{B_by_a, wdg_aB};

        // fmt::println("");
        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   ab  = mvec3d_e(a * b)            = {}", ab);
        // fmt::println("   abm = mva * mvb                  = {}", abm);
        // fmt::println("   abd = mvec3d(dot(a,b), wdg(a,b)) = {}", abd);
        // fmt::println("");
        // fmt::println("   A                                = {}", A);
        // fmt::println("   mvA                              = {}", mvA);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   Ab  = mvec3d_u(A * b)            = {}", Ab);
        // fmt::println("   Abm = mvA * mvb                  = {}", Abm);
        // fmt::println("   Abd = mvec3d((b << A), wdg(A,b)) = {}", Abd);
        // fmt::println("");
        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   B                                = {}", B);
        // fmt::println("   mvB                              = {}", mvB);
        // fmt::println("   aB  = mvec3d_u(a * B)            = {}", aB);
        // fmt::println("   aBm = mva * mvB                  = {}", aBm);
        // fmt::println("   aBd = mvec3d((B >> a), wdg(a,B)) = {}", aBd);
        // fmt::println("");

        CHECK(gr0(ab) == gr0(abm));
        CHECK(gr1(abm) == vec3d{});
        CHECK(gr2(ab) == gr2(abm));
        CHECK(gr3(abm) == pscalar3d{0.0});

        CHECK(gr0(ab) == gr0(abd));
        CHECK(gr1(abd) == vec3d{});
        CHECK(gr2(ab) == gr2(abd));
        CHECK(gr3(abd) == pscalar3d{0.0});

        CHECK(gr0(Abm) == 0);
        CHECK(gr1(Ab) == gr1(Abm));
        CHECK(gr2(Abm) == bivec3d{});
        CHECK(gr3(Ab) == gr3(Abm));

        CHECK(gr0(Abd) == 0);
        CHECK(gr1(Ab) == gr1(Abd));
        CHECK(gr2(Abd) == bivec3d{});
        CHECK(gr3(Ab) == gr3(Abd));

        CHECK(gr0(aBm) == 0);
        CHECK(gr1(aB) == gr1(aBm));
        CHECK(gr2(aBm) == bivec3d{});
        CHECK(gr3(aB) == gr3(aBm));

        CHECK(gr0(aBd) == 0);
        CHECK(gr1(aB) == gr1(aBd));
        CHECK(gr2(aBd) == bivec3d{});
        CHECK(gr3(aB) == gr3(aBd));
    }

    TEST_CASE("MVec3d: assignment tests")
    {
        fmt::println("MVec3d: assignment tests");

        vec3d v1{1.0, 2.0, 3.0};
        vec3d v2{0.5, 1.0, 1.5};
        vec3d v3{0.5, 1.0, -4.5};
        bivec3d b1{1.0, 2.0, 3.0};

        mvec3d mv1{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        mvec3d mv2{0.0, 0.5, 1.0, 1.5, 11.5, 15.5, 6.0, 61.5};
        mvec3d mv3{mv1};
        mvec3d mv4 = mv2;

        mvec3d mv5(scalar3d(5.0));
        mvec3d mv6(pscalar3d(6.0));
        mvec3d mv7{v1};
        mvec3d mv8{b1};
        mvec3d mv9{scalar3d{dot(v1, v3)}, wdg(v1, v3)};

        mvec3d mv10{v1, pscalar3d(10.0)};
        // This must not compile! Implict conversion to vec3d possible
        // possible solution: explicitly deleted constructor for MVec3d
        // mvec3d mv11{b1, pscalar3d_t(10.0)};

        // this does not compile (which is fine, a base cannot convert to derived)
        // mvec3d mv12{scalar3d(10.0), v1};

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("");
        // fmt::println("   mv1 = {}", mv1);
        // fmt::println("   mv2 = {}", mv2);
        // fmt::println("   mv3 = {}", mv3);
        // fmt::println("   mv4 = {}", mv4);
        // fmt::println("   mv5 = {}", mv5);
        // fmt::println("   mv6 = {}", mv6);
        // fmt::println("");
        // fmt::println("   gr1(mv1) = {}", gr1(mv1));
        // fmt::println("   gr1(mv2) = {}", gr1(mv2));
        // fmt::println("   gr1(mv3) = {}", gr1(mv3));
        // fmt::println("   gr1(mv3) = {}", gr1(mv4));
        // fmt::println("");
        // fmt::println("   v1 = {}", v1);
        // fmt::println("   mv7 = v1 = {}", mv7);
        // fmt::println("   b1 = {}", b1);
        // fmt::println("   mv8 = b1 = {}", mv8);
        // fmt::println("");
        // fmt::println("   mv9 = {}", mv9);
        // fmt::println("   mv10 = {}", mv10);

        CHECK(gr1(mv1) == v1);
        CHECK(gr1(mv2) == v2);
        CHECK(gr1(mv3) == v1);
        CHECK(gr1(mv4) == v2);
        CHECK(gr0(mv5) == 5.0);
        CHECK(gr3(mv6) == 6.0);
        CHECK(mv1 == mv3);
        CHECK(mv4 == mv2);
        CHECK(gr1(mv7) == v1);
        CHECK(gr2(mv8) == b1);
        CHECK(gr0(mv9) == dot(v1, v3));
        CHECK(gr2(mv9) == wdg(v1, v3));
        CHECK(gr1(mv10) == v1);
        CHECK(gr3(mv10) == 10.0);
    }

    TEST_CASE("MVec3d: bivector product properties")
    {
        fmt::println("MVec3d: bivector product properties");

        bivec3d b1{1.0, 2.0, 3.0};
        mvec3d mb1{b1};
        bivec3d b2{-3.0, 1.0, 2.0};
        mvec3d mb2{b2};

        auto gpr12_m = mb1 * mb2;
        auto gpr21_m = mb2 * mb1;
        auto gpr12_m_sym = 0.5 * (gpr12_m + gpr21_m);
        auto gpr12_m_asym = 0.5 * (gpr12_m - gpr21_m);

        auto gpr12_d = b1 * b2;
        auto gpr21_d = b2 * b1;
        auto gpr12_d_sym = 0.5 * (gpr12_d + gpr21_d);
        auto gpr12_d_asym = 0.5 * (gpr12_d - gpr21_d);

        // fmt::println("   b1  = {}", b1);
        // fmt::println("   mb1 = {}", mb1);
        // fmt::println("   b2  = {}", b2);
        // fmt::println("   mb2 = {}", mb2);
        // fmt::println("");
        // fmt::println("   dot(b1, b2) = {}", dot(b1, b2));
        // fmt::println("   cmt(b1, b2) = {}", cmt(b1, b2));
        // fmt::println("   dot(b2, b1) = {}", dot(b2, b1));
        // fmt::println("   cmt(b2, b1) = {}", cmt(b2, b1));
        // fmt::println("");
        // fmt::println("   gpr12_m = mb1 * mb2 = {}", gpr12_m);
        // fmt::println("   gpr21_m = mb2 * mb1 = {}", gpr21_m);
        // fmt::println("   gpr12_m_sym  = 0.5*(gpr12_d + gpr21_d) = {}", gpr12_m_sym);
        // fmt::println("   gpr12_m_asym = 0.5*(gpr12_m - gpr21_m) = {}", gpr12_m_asym);
        // fmt::println("");
        // fmt::println("   gpr12_d = b1 * b2 = {} ", gpr12_d);
        // fmt::println("   gpr21_d = b2 * b1) = {} ", gpr21_d);
        // fmt::println("   gpr12_d_sym  = 0.5*(gpr12_d + gpr21_d) = {}", gpr12_d_sym);
        // fmt::println("   gpr12_d_asym = 0.5*(gpr12_d - gpr21_d) = {}", gpr12_d_asym);
        // fmt::println("");

        CHECK(gr2(mb1) == b1);

        // just to silence unused variable warnings
        CHECK(gpr12_m_sym == 0.5 * (gpr12_m + gpr21_m));
        CHECK(gpr12_m_asym == 0.5 * (gpr12_m - gpr21_m));
        CHECK(gpr12_d_sym == 0.5 * (gpr12_d + gpr21_d));
        CHECK(gpr12_d_asym == 0.5 * (gpr12_d - gpr21_d));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec3d_E<T> and MVec3d_U<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec3d_E/_U: modelling even and uneven parts of 3d algebra - basics")
    {
        fmt::println(
            "MVec3d_E/_U: modelling even and uneven parts of 3d algebra - basics");

        // defining a complex number in all three forms as multivector
        auto u = normalize(vec3d{1.0, 0.0, 0.0});
        auto v =
            normalize(vec3d(std::cos(pi / 12), std::sin(pi / 12), 0.0)); // unit vec +15%
        auto angle_uv = angle(u, v);
        auto B = wdg(u, v); // normalized bivector describing the plane spanned by u and v

        auto my_exp = exp(-B, angle_uv);
        auto my_rot = rotor(B, 2.0 * angle_uv);

        // definition of rotor used here: B = u^v
        // => B determines the meaning of the positive sign of the rotation
        //
        auto R_m =
            mvec3d(exp(-B, angle_uv)); // Rotor formed by u and v (normalized bivector)
        auto Rr_m = mvec3d(rev(R_m));  // and its reverse

        auto c = vec3d{1.0, 1.0, 1.0};
        auto c_m = mvec3d{c};

        auto c_tmp_m = R_m * c_m;
        auto c_rot_m = c_tmp_m * Rr_m;

        auto R = exp(-B, angle_uv); // Rotor formed by u and v (normalized bivector)
        auto Rr = rev(R);           // and its reverse

        auto c_tmp_l = R * c;
        auto c_rot_u_l = c_tmp_l * Rr;
        auto c_rot_l = gr1(c_rot_u_l);
        // due to symmetry of R and Rr the gr3(c_rot) part will be zero
        // and thus can be assumed to be zero for further computations

        auto c_tmp_r = c * Rr;
        auto c_rot_u_r = R * c_tmp_r;
        auto c_rot_r = gr1(c_rot_u_r);
        // due to symmetry of R and Rr the gr3(c_rot) part will be zero
        // and thus can be assumed to be zero for further computations

        auto angle_c_c_rot = angle(c, c_rot_l); // not that easy in 3D!
        // (angle in plane of both vectors is not the angle in the plane
        // represented by the bivector!)
        // => requires projection of vectors onto plane and then taking
        // the angle between the projected vectors to be correct (bivector angle!)

        auto c_proj = project_onto(c, B);
        auto c_rot_proj = project_onto(c_rot_l, B);
        auto angle_proj = angle(c_proj, c_rot_proj);

        // fmt::println("   u                     = {: .3}", u);
        // fmt::println("   v                     = {: .3}", v);
        // fmt::println("   B = u^v = wdg(u,v)    = {: .3}", B);
        // fmt::println("   angle(u,v)            = {: .3}°", rad2deg(angle_uv));
        // fmt::println("   sin(angle(u,v))       = {: .3}", std::sin(angle_uv));
        // fmt::println("");
        // fmt::println("   c                     = {: .3}", c);
        // fmt::println("");
        // fmt::println("Implemented as full multivector operation:");
        // fmt::println("   R_m  = mvec3d(exp(-B,angle_uv))  = {: .3}", R_m);
        // fmt::println("   Rr_m = rev(R_m)                  = {: .3}", Rr_m);
        // fmt::println("   Rr_m*R_m                         = {: .3}", Rr_m * R_m);
        // fmt::println("   c_m                              = {: .3}", c_m);
        // fmt::println("   c_tmp_m = R_m*c_m                = {: .3}", c_tmp_m);
        // fmt::println("   c_rot_m = c_tmp_m*Rr_m           = {: .3}", c_rot_m);
        // fmt::println("   gr1(c_rot_m)                     = {: .3}", gr1(c_rot_m));
        // fmt::println("");
        // fmt::println("Implemented as reduced grade multivector operation:");
        // fmt::println("   R  = exp(-B,angle_uv)            = {: .3}", R);
        // fmt::println("   Rr = rev(R)                      = {: .3}", Rr);
        // fmt::println("   my_exp = exp(-B, angle_uv)       = {: .3}", my_exp);
        // fmt::println("   my_rot = rotor(B, 2*angle_uv)    = {: .3}", my_rot);
        // fmt::println("");
        // fmt::println("Left multiplication of rotor first:");
        // fmt::println("   c_tmp_l = R*c            = {: .3}", c_tmp_l);
        // fmt::println("   c_rot_u_l = c_tmp_l*Rr   = {: .3}", c_rot_u_l);
        // fmt::println("   c_rot_l = gr1(c_rot_u_l) = {: .3}", c_rot_l);
        // fmt::println("");
        // fmt::println("Right multiplication of rotor first:");
        // fmt::println("   c_tmp_r = c*Rr           = {: .3}", c_tmp_r);
        // fmt::println("   c_rot_u_r = R*c_tmp_r    = {: .3}", c_rot_u_r);
        // fmt::println("   c_rot_r = gr1(c_rot_u_r) = {: .3}", c_rot_r);
        // fmt::println("");
        // fmt::println("   angle(c, c_rot_l) = {: .3}°", rad2deg(angle_c_c_rot));
        // fmt::println("   angle(projected)  = {: .3}°", rad2deg(angle_proj));
        // fmt::println("");
        // fmt::println("direct calclulation:");
        // fmt::println("   c_rot = rotate(c,R)          = {: .3}", rotate(c, R));

        CHECK(nrm(rotate(c, R)) == nrm(c));
        CHECK(gr1(c_rot_m) == rotate(c, R));
        CHECK(rotate(c, R) == rotate_opt1(c, R));
        CHECK(rotate(c, R) == rotate_opt2(c, R));
        // n I_3d approach:
        CHECK(rotate(vec3d{1.0, 0.0, 0.0}, rotor(e3_3d * I_3d, pi / 4)) ==
              normalize(vec3d{1.0, 1.0, 0.0}));
        // using a bivector directly:
        CHECK(rotate(vec3d{1.0, 0.0, 0.0}, rotor(e12_3d, pi / 4)) ==
              normalize(vec3d{1.0, 1.0, 0.0}));

        // direct rotation of a bivector
        CHECK(rotate(bivec3d{0.0, 0.0, 1.0}, rotor(e23_3d, pi / 2)) == -e31_3d);

        // example see Macdonald "Linear and Geometric Algebra", Exercise 7.12, p. 127
        auto Bv =
            wdg(e2_3d, e1_3d + std::sqrt(3.0) * e3_3d); // bivector describing the plane
        CHECK(std::abs(nrm(Bv) - 2.0) < eps);
        CHECK(rotate(Bv, rotor(e31_3d, pi / 3)) == -2.0 * e12_3d);
        CHECK(rotate(Bv, rotor(e31_3d, pi / 3)) ==
              rotate_opt1(Bv, rotor(e31_3d, pi / 3)));
        CHECK(rotate(Bv, rotor(e31_3d, pi / 3)) ==
              rotate_opt2(Bv, rotor(e31_3d, pi / 3)));

        // just to silence unused variable warnings
        CHECK(my_exp == exp(-B, angle_uv));
        CHECK(my_rot == rotor(B, 2.0 * angle_uv));
        CHECK(c_rot_r == gr1(c_rot_u_r));
        CHECK(angle_c_c_rot == angle(c, c_rot_l));
        CHECK(angle_proj == angle(c_proj, c_rot_proj));
    }

    TEST_CASE("MVec3d: dualization")
    {
        fmt::println("MVec3d: dualization");

        vec3d v{1.0, 2.0, 3.0};                                   // 3d vector
        bivec3d B{10.0, 20.0, 30.0};                              // 3d bivector
        mvec3d vm{100.0, 1.0, 2.0, 3.0, 10.0, 20.0, 30.0, 300.0}; // full 3d multivector
        mvec3d vm2{-20.0, 3.0,  7.0,    -4.0,
                   -1.0,  40.0, -330.0, -70.0}; // full 3d multivector

        // full 3d multivector - even content
        mvec3d vm_even{100.0, 0.0, 0.0, 0.0, 10.0, 20.0, 30.0, 0.0};
        // even grade 3d multivector
        mvec3d_e vm_E{100.0, 10.0, 20.0, 30.0};

        // full 3d multivector - uneven content
        mvec3d vm_uneven{0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 300.0};
        // uneven grade 3d multivector
        mvec3d_u vm_U{1.0, 2.0, 3.0, 300.0};

        ////////////////////////////////////////////////////////////////////////////////
        // duality (as defined in Lengyel, "PGA illuminated"):
        ////////////////////////////////////////////////////////////////////////////////
        //
        // dual(A) = cmpl(A) in spaces of uneven dimension
        //         = rcmpl(A) in spaces of even dimension
        //

        auto vm_dual_manual = rev(vm) * Im_3d;
        auto vm_dual = dual(vm);

        auto vm_dual_even_manual = rev(vm_even) * Im_3d;
        auto vm_dual_even = dual(vm_even);

        auto vm_dual_uneven_manual = rev(vm_uneven) * Im_3d;
        auto vm_dual_uneven = dual(vm_uneven);

        // result is uneven, naming chosen for consistency
        auto vm_dual_manual_E = rev(vm_E) * I_3d;
        auto vm_dual_E = dual(vm_E);

        // result is even, naming chosen for consistency
        auto vm_dual_manual_U = rev(vm_U) * Im_3d_U;
        auto vm_dual_U = dual(vm_U);

        auto v_dual_manual = rev(v) * I_3d;
        auto v_dual = dual(v);

        auto B_dual_manual = rev(B) * I_3d;
        auto B_dual = dual(B);

        // fmt::println("");
        // fmt::println("   I_3d                   = {}", I_3d);
        // fmt::println("   Im_3d                  = {}", Im_3d);
        // fmt::println("   Im_3d_U                = {}", Im_3d_U);
        // fmt::println("");
        // fmt::println("   v                      = {}", v);
        // fmt::println("   B                      = {}", B);
        // fmt::println("");
        // fmt::println("   vm                     = {}", vm);
        // fmt::println("   rev(vm) * Im_3d        = {}", vm_dual_manual);
        // fmt::println("   dual(vm)               = {}", vm_dual);
        // fmt::println("");
        // fmt::println("   vm_even                = {}", vm_even);
        // fmt::println("   rev(vm_even) * Im_3d   = {}", vm_dual_even_manual);
        // fmt::println("   dual(vm_even)          = {}", vm_dual_even);
        // fmt::println("");
        // fmt::println("   vm_E                   = {}", vm_E);
        // fmt::println("   rev(vm_E) * I_3d       = {}", vm_dual_manual_E);
        // fmt::println("   dual(vm_E)             = {}", vm_dual_E);
        // fmt::println("");
        // fmt::println("   vm_uneven              = {}", vm_uneven);
        // fmt::println("   rev(vm_uneven) * Im_3d = {}", vm_dual_uneven_manual);
        // fmt::println("   dual(vm_uneven)        = {}", vm_dual_uneven);
        // fmt::println("");
        // fmt::println("   vm_U                   = {}", vm_U);
        // fmt::println("   rev(vm_U) * Im_3d_U    = {}", vm_dual_manual_U);
        // fmt::println("   dual(vm_U)             = {}", vm_dual_U);
        // fmt::println("");
        // fmt::println("   v                      = {}", v);
        // fmt::println("   rev(v)*I_3d            = {} - bivec ", v_dual_manual);
        // fmt::println("   dual(v)                = {} - bivec ", v_dual);
        // fmt::println("");
        // fmt::println("   B                      = {}", B);
        // fmt::println("   rev(B)*I_3d            = {} - vec", B_dual_manual);
        // fmt::println("   dual(B)                = {} - vec", B_dual);
        // fmt::println("");

        CHECK(vm_dual == vm_dual_manual);
        CHECK(vm_dual_even == vm_dual_even_manual);
        CHECK(vm_dual_uneven == vm_dual_uneven_manual);
        CHECK(vm_dual_E == vm_dual_manual_E);
        CHECK(vm_dual_U == vm_dual_manual_U);
        CHECK(dual(v) == bivec3d{1.0, 2.0, 3.0});
        CHECK(dual(B) == vec3d{10.0, 20.0, 30.0});
        CHECK(dual(scalar3d(5)) == pscalar3d(5));
        CHECK(dual(pscalar3d(6)) == scalar3d(6));

        // dual properties
        CHECK(dual(3.0 * vm) == 3.0 * dual(vm));
        CHECK(dual(vm + vm2) == dual(vm) + dual(vm2));
        CHECK(dual(dual(vm)) == vm);
        CHECK(dual(I_3d) == 1);
        CHECK(dual(v) == bivec3d(1.0, 2.0, 3.0));

        // inner and outer products are in G^n are dual (A. Macdonald, p. 111):
        //
        // dual(dot(A,B)) = wdg(A,dual(B)), dual(wdg(A,B)) = dot(A,dual(B))

        auto a = vec3d{2.0, -3.0, 2.0};
        auto b = vec3d{-5.0, 1.0, -7.0};

        // fmt::println("   a                 = {}", a);
        // fmt::println("   b                 = {}", b);
        // fmt::println("   dual(b)           = {}", dual(b));
        // fmt::println("   dot(a, b)         = {}", dot(a, b));
        // fmt::println("   wdg(a, b)         = {}", wdg(a, b));
        // fmt::println("");
        // fmt::println("   dual(dot(a, b))   = {}", dual(dot(a, b)));
        // fmt::println("   wdg(a, dual(b))   = {}", wdg(a, dual(b)));
        // fmt::println("");
        // fmt::println("   dual(wdg(a, b))   = {}", dual(wdg(a, b)));
        // fmt::println("   dot(a, dual(b))   = {}", dot(a, dual(b)));

        CHECK(dual(scalar3d{dot(a, b)}) == wdg(a, dual(b)));
        CHECK(dual(scalar3d(dot(a, b))) == wdg(a, dual(b)));
        // CHECK(dual(wdg(a, b)) == dot(a, dual(b)));
        CHECK(wdg(a, dual(b)) == dot(a, b) * I_3d);

        // just to silence unused variable warning
        CHECK(v_dual_manual == rev(v) * I_3d);
        CHECK(v_dual == dual(v));
        CHECK(B_dual_manual == rev(B) * I_3d);
        CHECK(B_dual == dual(B));
    }

    TEST_CASE("MVec3d: product tests")
    {
        fmt::println("MVec3d: product tests");
        // tests based on "The inner products of geometric algebra", Leo Dorst


        auto s1 = scalar3d{2.0};
        auto v1 = vec3d{1.0, -3.0, 2.0};
        auto b1 = bivec3d{2.0, -4.0, 1.0};
        auto ps1 = pscalar3d{-2.0};
        auto M1 = mvec3d{s1, v1, b1, ps1};

        auto s2 = scalar3d{-1.0};
        auto v2 = vec3d{2.0, 1.5, -1.0};
        auto b2 = bivec3d{-2.0, 10.5, 3.0};
        auto ps2 = pscalar3d{3.0};
        auto M2 = mvec3d{s2, v2, b2, ps2};

        auto s3 = scalar3d{-10.0};
        auto v3 = vec3d{-2.0, 4.0, -10.0};
        auto b3 = bivec3d{1.0, -2.0, 3.0};
        auto ps3 = pscalar3d{-12.0};
        auto M3 = mvec3d{s3, v3, b3, ps3};

        auto V1 = mvec3d{v1};
        auto V2 = mvec3d{v2};
        auto V3 = mvec3d{v3};
        auto B12 = mvec3d{wdg(v1, v2)};

        CHECK(nrm_sq(M1) == nrm_sq(s1) + nrm_sq(v1) + nrm_sq(b1) + nrm_sq(ps1));
        CHECK(nrm_sq(M2) == nrm_sq(s2) + nrm_sq(v2) + nrm_sq(b2) + nrm_sq(ps2));

        CHECK(dot(M1, M2) == dot(s1, s2) + dot(v1, v2) + dot(b1, b2) + dot(ps1, ps2));
        CHECK(dot(M1, M1) == nrm_sq(M1));
        CHECK(dot(M2, M2) == nrm_sq(M2));

        // wedge product
        CHECK(wdg(v1, v2) == gr2(wdg(mvec3d(v1), mvec3d(v2))));
        // check full permissible range of arguments, even the ones delivering 0)
        CHECK(wdg(b1, b2) == scalar3d(0.0));
        CHECK(wdg(v1, ps1) == scalar3d(0.0));
        CHECK(wdg(ps1, v1) == scalar3d(0.0));
        CHECK(wdg(b1, ps1) == scalar3d(0.0));
        CHECK(wdg(ps1, b1) == scalar3d(0.0));
        CHECK(wdg(ps1, ps1) == scalar3d(0.0));
        CHECK(wdg(b1, b2) == nrm_sq(b1 * b2 - cmt(b1, b2) + dot(b1, b2)));
        CHECK(wdg(v1, ps1) == nrm_sq(0.5 * (v1 * ps1 + rev(ps1) * v1)));
        CHECK(wdg(ps1, v1) == nrm_sq(0.5 * (ps1 * v1 + v1 * rev(ps1))));


        CHECK((mvec3d{s2, v2, b2, pscalar3d{0.0}} << mvec3d{s1}) ==
              mvec3d{s1 * s2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

        CHECK((mvec3d{v2} << mvec3d{s1}) ==
              mvec3d{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});

        CHECK((wdg(v2, v3) >> v1) == wdg(v1 << v2, v3) + wdg(gr_inv(v2), v1 << v3));
        CHECK((wdg(v1, v2) << b1) == (v1 << (v2 << b1)));

        // contractions - check full permissible range of arguments, even the ones
        // delivering 0 as a result
        CHECK((ps1 << s1) == 0.0);
        CHECK((ps1 << v1) == 0.0);
        CHECK((ps1 << b1) == 0.0);
        CHECK((b1 << s1) == 0.0);
        CHECK((b1 << v1) == 0.0);
        CHECK((v1 << s1) == 0.0);
        //
        CHECK((s1 >> ps1) == 0.0);
        CHECK((v1 >> ps1) == 0.0);
        CHECK((b1 >> ps1) == 0.0);
        CHECK((s1 >> b1) == 0.0);
        CHECK((v1 >> b1) == 0.0);
        CHECK((s1 >> v1) == 0.0);

        // 2.3.2
        // 1st with vectors lifted to full multivectors
        CHECK(V1 * M2 == (M2 >> V1) + wdg(V1, M2));
        CHECK(M2 * V1 == (V1 << M2) + wdg(M2, V1));
        CHECK((M2 >> V1) == -(V1 << gr_inv(M2))); // replacing the right contraction
        CHECK(V1 * M2 == -(V1 << gr_inv(M2)) + wdg(V1, M2)); // yep, replacing works

        // 2nd with vectors and multivectors directly
        CHECK(v1 * M2 == (M2 >> v1) + wdg(v1, M2));
        CHECK(M2 * v1 == (v1 << M2) + wdg(M2, v1));
        CHECK((M2 >> v1) == -(v1 << gr_inv(M2)));
        CHECK(wdg(v1, M2) == 0.5 * (v1 * M2 + gr_inv(M2) * v1));
        CHECK(wdg(M2, v1) == 0.5 * (M2 * v1 + v1 * gr_inv(M2)));

        // 3rd with vectors and bivectors directly
        CHECK(v1 * b1 == (b1 >> v1) + wdg(v1, b1));
        CHECK(b1 * v1 == (v1 << b1) + wdg(b1, v1));
        CHECK((b1 >> v1) == -(v1 << gr_inv(b1)));


        // express contraction and wedge product through the geometric product
        CHECK(wdg(V1, M2) == 0.5 * (V1 * M2 + gr_inv(M2) * V1));
        CHECK(wdg(M2, V1) == 0.5 * (M2 * V1 + V1 * gr_inv(M2)));

        // product dualities
        CHECK(dual((v1 << v2)) == wdg(v1, dual(v2)));
        CHECK(dual(wdg(v1, v2)) == (v1 << dual(v2)));
        CHECK(dual((M1 << M2)) == wdg(M1, dual(M2)));
        CHECK(dual(wdg(M1, M2)) == (M1 << dual(M2)));

        // cross-check direct implementation of rwdg by comparing with wdg
        CHECK(rwdg(M1, M2) == cmpl(wdg(cmpl(M1), cmpl(M2))));
    }

    TEST_CASE("MVec3d: simple applications, complements, contraction, expansions")
    {
        fmt::println("MVec3d: simple applications, complements, contraction, expansions");

        auto s1 = scalar3d{2.0};
        auto v1 = vec3d{1.0, -3.0, 2.0};
        auto b1 = bivec3d{2.0, -4.0, 1.0};
        auto ps1 = pscalar3d{-2.0};
        auto M1 = mvec3d{s1, v1, b1, ps1};


        // simple projections
        auto v = vec3d{4.0, 3.5, 5.0};
        // auto B = e31_3d;
        auto B = e31_3d + 0.1 * e12_3d;
        // auto B = e31_3d + 0.01 * e12_3d + 0.1 * e23_3d;


        auto v_in_B = project_onto(v, B);
        auto v_perp_B = reject_from(v, B);

        // fmt::println("   v         = {:.3f}", v);
        // fmt::println("   B         = {:.3f}", B);
        // fmt::println("");
        // fmt::println("   v_in_B    = {:.3f}", v_in_B);
        // fmt::println("   v_perp_B  = {:.3f}", v_perp_B);
        // fmt::println("   v << B    = {:.3f}", v << B);
        // fmt::println("   B >> v    = {:.3f}", B >> v);
        // fmt::println("");

        // x^b = 0 for every point x in B
        CHECK(nrm_sq(wdg(v_in_B, B)) < eps);

        // v_perp_B should be proportional to the normal vector n of B
        // n = dual(B) or n = cmpl(B)
        // thus, wdg(v_perp_B, n) == 0 is required
        CHECK(nrm_sq(wdg(v_perp_B, dual(B))) < eps);
        CHECK(nrm_sq(wdg(v_perp_B, cmpl(B))) < eps);

        // v_in_B and v_perp_B should be perpendicular to each other
        CHECK(nrm_sq(dot(v_in_B, v_perp_B)) < eps);

        // v should be the sum of v_in_B and v_perp_B
        CHECK(v == v_in_B + v_perp_B);

        // check complements
        // fmt::println("");
        // auto s = scalar3d(1);
        // auto e1 = e1_3d;
        // auto e2 = e2_3d;
        // auto e3 = e3_3d;
        // auto e23 = e23_3d;
        // auto e31 = e31_3d;
        // auto e12 = e12_3d;
        // auto ps = I_3d;
        // fmt::println("   u:  s, e1, e2, e3, e23, e31, e12, ps \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}, \n"
        //              "{:.3f}, {:.3f}, {:.3f}, {:.3f}",
        //              s, e1, e2, e3, e23, e31, e12, ps);
        // fmt::println("");
        // fmt::println("   cmpl(u):  cmpl(s), cmpl(e1), cmpl(e2), cmpl(e3) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              cmpl(s), cmpl(e1), cmpl(e2), cmpl(e3));
        // fmt::println("");
        // fmt::println("   cmpl(u):  cmpl(e23), cmpl(e31), cmpl(e12), cmpl(ps) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              cmpl(e23), cmpl(e31), cmpl(e12), cmpl(ps));
        // fmt::println("");
        // fmt::println("   rev(u)*I_3d:  rev(s)*I_3d, rev(e1)*I_3d, rev(e2)*I_3d, "
        //              "rev(e3)*I_3d \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              rev(s) * I_3d, rev(e1) * I_3d, rev(e2) * I_3d, rev(e3) * I_3d);
        // fmt::println("");
        // fmt::println("   rev(u)*I_3d:  rev(e23)*I_3d, rev(e31)*I_3d, "
        //              "rev(e12)*I_3d, rev(ps)*I_3d \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              rev(e23) * I_3d, rev(e31) * I_3d, rev(e12) * I_3d, rev(ps) *
        //              I_3d);
        // fmt::println("");
        // fmt::println("   I_3d*rev(u):  I_3d*rev(s), I_3d*rev(e1), I_3d*rev(e2), "
        //              "I_3d*rev(e3) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              I_3d * rev(s), I_3d * rev(e1), I_3d * rev(e2), I_3d * rev(e3));
        // fmt::println("");
        // fmt::println("   I_3d*rev(u):  I_3d*rev(e23), I_3d*rev(e31), "
        //              "I_3d*rev(e13), I_3d*rev(ps) \n"
        //              "   = {:.3f}: {:.3f}, {:.3f}, {:.3f}",
        //              I_3d * rev(e23), I_3d * rev(e31), I_3d * rev(e12), I_3d *
        //              rev(ps));
        // fmt::println("");

        CHECK(cmpl(scalar3d(1.0)) == I_3d);
        CHECK(cmpl(scalar3d(1.0)) == rev(scalar3d(1.0)) * I_3d);
        CHECK(cmpl(e1_3d) == e23_3d);
        CHECK(cmpl(e1_3d) == rev(e1_3d) * I_3d);
        CHECK(cmpl(e2_3d) == e31_3d);
        CHECK(cmpl(e2_3d) == rev(e2_3d) * I_3d);
        CHECK(cmpl(e3_3d) == e12_3d);
        CHECK(cmpl(e3_3d) == rev(e3_3d) * I_3d);
        CHECK(cmpl(e23_3d) == e1_3d);
        CHECK(cmpl(e23_3d) == rev(e23_3d) * I_3d);
        CHECK(cmpl(e31_3d) == e2_3d);
        CHECK(cmpl(e31_3d) == rev(e31_3d) * I_3d);
        CHECK(cmpl(e12_3d) == e3_3d);
        CHECK(cmpl(e12_3d) == rev(e12_3d) * I_3d);
        CHECK(cmpl(I_3d) == scalar2d(1.0));
        CHECK(cmpl(I_3d) == rev(I_3d) * I_3d);
        //
        CHECK(cmpl(cmpl(M1)) == M1);
        CHECK(cmpl(cmpl(mvec3d_e(s1, b1))) == mvec3d_e(s1, b1));
        CHECK(cmpl(cmpl(mvec3d_u(v1, ps1))) == mvec3d_u(v1, ps1));
        //
        CHECK(wdg(scalar3d(5), cmpl(scalar3d(5))) / nrm_sq(scalar3d(5)) == I_3d);
        CHECK(wdg(cmpl(scalar3d(5)), scalar3d(5)) / nrm_sq(scalar3d(5)) == I_3d);
        CHECK(wdg(v, cmpl(v)) / nrm_sq(v) == I_3d);
        CHECK(wdg(cmpl(v), v) / nrm_sq(v) == I_3d);
        CHECK(wdg(B, cmpl(B)) / nrm_sq(B) == I_3d);
        CHECK(wdg(cmpl(B), B) / nrm_sq(B) == I_3d);
        CHECK(wdg(pscalar3d(3), cmpl(pscalar3d(3))) / nrm_sq(pscalar3d(3)) == I_3d);
        CHECK(wdg(cmpl(pscalar3d(3)), pscalar3d(3)) / nrm_sq(pscalar3d(3)) == I_3d);

        // check contractions: <<, >> and rwdg( u, compl(v) )
        // fmt::println("   v         = {:.3f}", v);
        // fmt::println("   cmpl(v)   = {:.3f}", cmpl(v));
        // fmt::println("");
        // fmt::println("   v << B = {:.3f}", v << b);
        // fmt::println("   B >> v = {:.3f}", b >> v);
        // fmt::println("");
        // fmt::println("   rwdg(b, cmpl(v))  = {:.3f}", rwdg(b, cmpl(v)));
        // fmt::println("   rwdg(cmpl(v), b)  = {:.3f}", rwdg(cmpl(v), b));
        // fmt::println("");
        CHECK((B >> v) == rwdg(B, cmpl(v)));
        CHECK((v << B) == rwdg(cmpl(v), B));

        // check expansions: v ^ cmpl(b)  and  cmpl(b) ^ v
        // (create new bivector that contains v and is perpendicular to b)
        //
        // also checks the duality correspondence:
        //      cmpl(v >> b) == v ^ cmpl(b)
        //      cmpl(v >> b) == cmpl(b) ^ v
        // fmt::println("   v      = {:.3f}", v);
        // fmt::println("   b      = {:.3f}", b);
        // fmt::println("   v << b = {:.3f}", v << b);
        // fmt::println("   b >> v = {:.3f}", b >> v);
        // fmt::println("");
        // fmt::println("   wdg(cmpl(b), v)         = {:.3f}", wdg(cmpl(b), v));
        // fmt::println("   n=cmpl(wdg(cmpl(b), v)) = {:.3f}", cmpl(wdg(cmpl(b), v)));
        // fmt::println("");
        // fmt::println("   wdg(v, cmpl(b))         = {:.3f}", wdg(v, cmpl(b)));
        // fmt::println("   n=cmpl(wdg(v, cmpl(b))) = {:.3f}", cmpl(wdg(v, cmpl(b))));
        // vector is in plane defined by the expansion
        CHECK(wdg(v, wdg(v, cmpl(B))) == 0.0);
        CHECK(wdg(v, wdg(cmpl(B), v)) == 0.0);
        // duality of the contraction and the wedge product (based on dual)
        CHECK(dual(v << B) == wdg(v, dual(B)));
        CHECK(dual(B >> v) == wdg(dual(B), v));
        // duality of the contraction and the wedge product (based on complement)
        CHECK(cmpl(v << B) == wdg(v, cmpl(B)));
        CHECK(cmpl(B >> v) == wdg(cmpl(B), v));
    }

    TEST_CASE("G<4,0,0>: ega_4d")
    {
        fmt::println("");
        fmt::println("G<4,0,0>: ega_4d");
        // 4d euklidean geometric algebra
        const algebra<4> alg;
        CHECK(alg.p() == 4);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 0);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   ega_4d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   ega_4d: basis_name = {}", fmt::join(alg.basis_name, ", "));
        fmt::println("");
    }

    TEST_CASE("G<4,0,0>: defining basic types and ctor checks")
    {
        fmt::println("G<4,0,0>: defining basic types and ctor checks");

        auto s1 = scalar4d{5.0};
        auto v1 = vec4d{1.0, 2.0, 3.0, 4.0};
        auto b1 = bivec4d{-1.0, 2.0, -3.0, 4.0, 5.0, -6.0};
        auto t1 = trivec4d{-1.0, -2.0, -3.0, -4.0};
        auto ps1 = pscalar4d{-5.0};

        auto mv1 = mvec4d{s1};
        auto mv2 = mvec4d{v1};
        auto mv3 = mvec4d{b1};
        auto mv4 = mvec4d{t1};
        auto mv5 = mvec4d{ps1};

        auto mv6a = mvec4d{s1, b1, ps1};
        auto mv6b = mvec4d_e{s1};
        auto mv6c = mvec4d_e{b1};
        auto mv6d = mvec4d_e{ps1};
        auto mv6e = mvec4d_e{s1, b1, ps1};
        auto mv6f = mvec4d{mv6e};

        auto mv7a = mvec4d{v1, t1};
        auto mv7b = mvec4d_u{v1};
        auto mv7c = mvec4d_u{t1};
        auto mv7d = mvec4d_u{v1, t1};
        auto mv7e = mvec4d{mv7d};

        auto mv8 = mvec4d{s1, v1, b1, t1, ps1};

        // fmt::println("   mv1  = {}", mv1);
        // fmt::println("   mv2  = {}", mv2);
        // fmt::println("   mv3  = {}", mv3);
        // fmt::println("   mv4  = {}", mv4);
        // fmt::println("   mv5  = {}", mv5);
        // fmt::println("   mv6a = {}", mv6a);
        // fmt::println("   mv6b = {}", mv6b);
        // fmt::println("   mv6c = {}", mv6c);
        // fmt::println("   mv6d = {}", mv6d);
        // fmt::println("   mv6e = {}", mv6e);
        // fmt::println("   mv6f = {}", mv6f);
        // fmt::println("   mv7a = {}", mv7a);
        // fmt::println("   mv7b = {}", mv7b);
        // fmt::println("   mv7c = {}", mv7c);
        // fmt::println("   mv7d = {}", mv7d);
        // fmt::println("   mv7e = {}", mv7e);
        // fmt::println("   mv8  = {}", mv8);

        CHECK(gr0(mv1) == s1);
        CHECK(gr1(mv2) == v1);
        CHECK(gr2(mv3) == b1);
        CHECK(gr3(mv4) == t1);
        CHECK(gr4(mv5) == ps1);

        CHECK(gr0(mv6a) == s1);
        CHECK(gr1(mv6a) == vec4d{});
        CHECK(gr2(mv6a) == b1);
        CHECK(gr3(mv6a) == trivec4d{});
        CHECK(gr4(mv6a) == ps1);

        CHECK(gr0(mv6b) == s1);
        CHECK(gr2(mv6b) == bivec4d{});
        CHECK(gr4(mv6b) == pscalar4d{});

        CHECK(gr0(mv6c) == scalar4d{});
        CHECK(gr2(mv6c) == b1);
        CHECK(gr4(mv6c) == pscalar4d{});

        CHECK(gr0(mv6d) == scalar4d{});
        CHECK(gr2(mv6d) == bivec4d{});
        CHECK(gr4(mv6d) == ps1);

        CHECK(gr0(mv6e) == s1);
        CHECK(gr2(mv6e) == b1);
        CHECK(gr4(mv6e) == ps1);

        CHECK(gr0(mv6f) == s1);
        CHECK(gr1(mv6f) == vec4d{});
        CHECK(gr2(mv6f) == b1);
        CHECK(gr3(mv6f) == trivec4d{});
        CHECK(gr4(mv6f) == ps1);

        CHECK(gr0(mv7a) == scalar4d{});
        CHECK(gr1(mv7a) == v1);
        CHECK(gr2(mv7a) == bivec4d{});
        CHECK(gr3(mv7a) == t1);
        CHECK(gr4(mv7a) == pscalar4d{});

        CHECK(gr1(mv7b) == v1);
        CHECK(gr3(mv7b) == trivec4d{});

        CHECK(gr1(mv7c) == vec4d{});
        CHECK(gr3(mv7c) == t1);

        CHECK(gr1(mv7d) == v1);
        CHECK(gr3(mv7d) == t1);

        CHECK(gr0(mv7e) == scalar4d{});
        CHECK(gr1(mv7e) == v1);
        CHECK(gr2(mv7e) == bivec4d{});
        CHECK(gr3(mv7e) == t1);
        CHECK(gr4(mv7e) == pscalar4d{});

        CHECK(gr0(mv8) == s1);
        CHECK(gr1(mv8) == v1);
        CHECK(gr2(mv8) == b1);
        CHECK(gr3(mv8) == t1);
        CHECK(gr4(mv8) == ps1);
    }

} // TEST_SUITE("Euclidean Geometric Algebra (EGA)")
