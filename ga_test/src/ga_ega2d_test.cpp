// Copyright 2024-2025, Daniel Hug. All rights reserved.

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


/////////////////////////////////////////////////////////////////////////////////////////
// Euclideon geometric algebra of two-dimensional space: ega2d
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("EGA 2D Tests")
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
        CHECK(reflect_on(p, right_dual(x_axis_2d)) == vec2d{4, -1});

        // coordinate axis reflected on perpendicular axis yield their negatives
        CHECK(reflect_on(y_axis_2d, right_dual(x_axis_2d)) == -y_axis_2d);
        CHECK(reflect_on(x_axis_2d, right_dual(y_axis_2d)) == -x_axis_2d);

        // coordinate axis reflected on itself remains itself (identity)
        CHECK(reflect_on(x_axis_2d, right_dual(x_axis_2d)) == x_axis_2d);
        CHECK(reflect_on(y_axis_2d, right_dual(y_axis_2d)) == y_axis_2d);
    }

    TEST_CASE("Vec2d: operations - rotations")
    {
        fmt::println("Vec2d: operations - rotations");

        std::vector<std::tuple<double, Vec2d<double>>> v;

        // fmt::println("");
        for (int i = -12; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = vec2d(std::cos(phi), std::sin(phi));
            v.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad2deg(phi), c, angle(e1_2d, c));
            CHECK(c == rotate(e1_2d, get_rotor(I_2d, phi)));
        }
        // fmt::println("");
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

        CHECK((e12_2d >> e1_2d) == e2_2d);
        CHECK((e12_2d >> e2_2d) == -e1_2d);
        CHECK((e1_2d << e12_2d) == -e2_2d);
        CHECK((e2_2d << e12_2d) == e1_2d);
        CHECK((e12_2d >> (2.0 * e1_2d - 3.0 * e2_2d)) == 3.0 * e1_2d + 2.0 * e2_2d);
        CHECK(((2.0 * e1_2d - 3.0 * e2_2d) << e12_2d) == -3.0 * e1_2d - 2.0 * e2_2d);
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

    TEST_CASE("MVec2d: geometric product - link to inner and outer products")
    {
        fmt::println("MVec2d: geometric product - link to inner and outer products");
        vec2d u{1.0, 2.0};

        scalar2d s{3};
        vec2d v{-3.0, 2.5};
        pscalar2d ps{4.0};

        // fmt::println("");
        // fmt::println("u = {}", u);
        // fmt::println("s = {}", s);
        // fmt::println("v = {}", v);
        // fmt::println("ps = {}", ps);
        // fmt::println("");

        // fmt::println("");
        // fmt::println("scalar case: ");
        CHECK(u * s == gr1((s >> u) + wdg(u, s)));
        CHECK(u * s == gr1(rwdg(s, right_dual(u)) + wdg(u, s)));

        // fmt::println("");
        // fmt::println("u * s = {}", u * s);
        // fmt::println("");
        // fmt::println("(s >> u) = {}", (s >> u));
        // fmt::println("wdg(u, s) = {}", wdg(u, s));
        // fmt::println("(s >> u) + wdg(u, s) = {}", (s >> u) + wdg(u, s));
        // fmt::println("gr1((s >> u) + wdg(u, s)) = {}", gr1((s >> u) + wdg(u, s)));
        // fmt::println("");
        // fmt::println("right_dual(u) = {}", right_dual(u));
        // fmt::println("rwdg(s,right_dual(u))= {}", rwdg(s, right_dual(u)));
        // fmt::println("wdg(u, s) = {}", wdg(u, s));
        // fmt::println("rwdg(s,right_dual(u)) + wdg(u, s) = {}",
        //              rwdg(s, right_dual(u)) + wdg(u, s));
        // fmt::println("gr1(rwdg(s,right_dual(u)) + wdg(u, s)) = {}",
        //              gr1(rwdg(s, right_dual(u)) + wdg(u, s)));
        // fmt::println("");

        CHECK(s * u == gr1((u << s) + wdg(s, u)));
        CHECK(s * u == gr1(rwdg(left_dual(u), s) + wdg(s, u)));

        // fmt::println("");
        // fmt::println("s * u = {}", s * u);
        // fmt::println("");
        // fmt::println("(u << s) = {}", (u << s));
        // fmt::println("wdg(s, u) = {}", wdg(s, u));
        // fmt::println("(u << s) + wdg(s, u) = {}", (u << s) + wdg(s, u));
        // fmt::println("gr1((u << s) + wdg(s, u)) = {}", gr1((u << s) + wdg(s, u)));
        // fmt::println("");
        // fmt::println("left_dual(u) = {}", left_dual(u));
        // fmt::println("rwdg(left_dual(u), s)= {}", rwdg(left_dual(u), s));
        // fmt::println("wdg(s, u) = {}", wdg(s, u));
        // fmt::println("rwdg(left_dual(u), s) + wdg(s, u) = {}",
        //              rwdg(left_dual(u), s) + wdg(s, u));
        // fmt::println("gr1(rwdg(left_dual(u), s) + wdg(s, u)) = {}",
        //              gr1(rwdg(left_dual(u), s) + wdg(s, u)));
        // fmt::println("");

        // fmt::println("");
        // fmt::println("vector case: ");
        CHECK(u * v == (v >> u) + wdg(u, v));
        CHECK(u * v == rwdg(v, right_dual(u)) + wdg(u, v));

        // fmt::println("");
        // fmt::println("u * v = {}", u * v);
        // fmt::println("");
        // fmt::println("(v >> u) = {}", (v >> u));
        // fmt::println("wdg(u, v) = {}", wdg(u, v));
        // fmt::println("(v >> u) + wdg(u, v) = {}", (v >> u) + wdg(u, v));
        // fmt::println("");
        // fmt::println("right_dual(u) = {}", right_dual(u));
        // fmt::println("rwdg(v,right_dual(u))= {}", rwdg(v, right_dual(u)));
        // fmt::println("wdg(u, v) = {}", wdg(u, v));
        // fmt::println("rwdg(v,right_dual(u)) + wdg(u, v) = {}",
        //              rwdg(v, right_dual(u)) + wdg(u, v));
        // fmt::println("");

        CHECK(v * u == (u << v) + wdg(v, u));
        CHECK(v * u == rwdg(left_dual(u), v) + wdg(v, u));

        // fmt::println("");
        // fmt::println("v * u = {}", v * u);
        // fmt::println("");
        // fmt::println("(u << v) = {}", (u << v));
        // fmt::println("wdg(v, u) = {}", wdg(v, u));
        // fmt::println("(u << v) + wdg(v, u) = {}", (u << v) + wdg(v, u));
        // fmt::println("");
        // fmt::println("left_dual(u) = {}", left_dual(u));
        // fmt::println("rwdg(left_dual(u), v)= {}", rwdg(left_dual(u), v));
        // fmt::println("wdg(v, u) = {}", wdg(v, u));
        // fmt::println("rwdg(left_dual(u), v) + wdg(v, u) = {}",
        //              rwdg(left_dual(u), v) + wdg(v, u));
        // fmt::println("");

        // fmt::println("");
        // fmt::println("pscalar case: ");
        CHECK(u * ps == gr1((ps >> u) + wdg(u, ps)));
        CHECK(u * ps == gr1(rwdg(ps, right_dual(u)) + wdg(u, ps)));

        // fmt::println("");
        // fmt::println("u * ps = {}", u * ps);
        // fmt::println("");
        // fmt::println("(ps >> u) = {}", (ps >> u));
        // fmt::println("wdg(u, ps) = {}", wdg(u, ps));
        // fmt::println("(ps >> u) + wdg(u, ps) = {}", (ps >> u) + wdg(u, ps));
        // fmt::println("gr1((ps >> u) + wdg(u, ps)) = {}", gr1((ps >> u) + wdg(u, ps)));
        // fmt::println("");
        // fmt::println("right_dual(u) = {}", right_dual(u));
        // fmt::println("rwdg(ps,right_dual(u))= {}", rwdg(ps, right_dual(u)));
        // fmt::println("wdg(u, ps) = {}", wdg(u, ps));
        // fmt::println("rwdg(ps,right_dual(u)) + wdg(u, ps) = {}",
        //              rwdg(ps, right_dual(u)) + wdg(u, ps));
        // fmt::println("gr1(rwdg(ps,right_dual(u)) + wdg(u, ps)) = {}",
        //              gr1(rwdg(ps, right_dual(u)) + wdg(u, ps)));
        // fmt::println("");

        CHECK(ps * u == gr1((u << ps) + wdg(ps, u)));
        CHECK(ps * u == gr1(rwdg(left_dual(u), ps) + wdg(ps, u)));

        // fmt::println("");
        // fmt::println("ps * u = {}", ps * u);
        // fmt::println("");
        // fmt::println("(u << ps) = {}", (u << ps));
        // fmt::println("wdg(ps, u) = {}", wdg(ps, u));
        // fmt::println("(u << ps) + wdg(ps, u) = {}", (u << ps) + wdg(ps, u));
        // fmt::println("gr1((u << ps) + wdg(ps, u)) = {}", gr1((u << ps) + wdg(ps, u)));
        // fmt::println("");
        // fmt::println("left_dual(u) = {}", left_dual(u));
        // fmt::println("rwdg(left_dual(u), ps)= {}", rwdg(left_dual(u), ps));
        // fmt::println("wdg(ps, u) = {}", wdg(ps, u));
        // fmt::println("rwdg(left_dual(u), ps) + wdg(ps, u) = {}",
        //              rwdg(left_dual(u), ps) + wdg(ps, u));
        // fmt::println("gr1(rwdg(left_dual(u), ps) + wdg(ps, u)) = {}",
        //              gr1(rwdg(left_dual(u), ps) + wdg(ps, u)));
        // fmt::println("");

        // fmt::println("");
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
        CHECK(rotate(u, get_rotor(I_2d, angle_uv)) ==
              v); // 2d rotation with double product
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
        // dual(A) = cmpl(A) in spaces of odd dimension
        //         = rcmpl(A) in spaces of even dimension
        //

        auto vm_dual_manual = rev(vm) * Im_2d;
        auto vm_dual = right_dual(vm);

        auto vm_dual_even_manual = rev(vm_even) * Im_2d;
        auto vm_dual_even = right_dual(vm_even);

        auto vm_dual_manual_E = rev(vm_E) * Im_2d_E;
        auto vm_dual_E = right_dual(vm_E);

        auto v_dual_manual = rev(v) * I_2d;
        auto v_dual = right_dual(v);

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
        CHECK(right_dual(scalar2d(5)) == pscalar2d(5));
        CHECK(right_dual(pscalar2d(5)) == scalar2d(5));
        CHECK(right_dual(I_2d) == 1);

        // dual properties (A. Macdonald, p. 110):
        //
        // a) dual(aA) = a dual(A)
        // b) dual(A + B) = dual(A) + dual(B)
        // c) dual(dual(A)) = (-1)^(n*(n-1)/2) A   (with n as dimension of the (sub)space)
        // d) |dual(B)| = |B|
        // e) if B is a j-blade then dual(B) is an (n-j)-blade
        // f) if A is a j-vector then dual(A) is an (n-j)-vector
        //    (remember: a j-vector is a sum of j-blades, which are outer products)

        CHECK(right_dual(3.0 * vm) == 3.0 * right_dual(vm));
        CHECK(right_dual(vm + vm2) == right_dual(vm) + right_dual(vm2));
        CHECK(right_dual(I_2d) == 1);
        CHECK(right_dual(v) == (vec2d(-v.y, v.x))); // = rcmpl(v)

        // inner and outer products are in G^n are dual (A. Macdonald, p. 111):
        //
        // dual(dot(A,B)) = wdg(A,dual(B)), dual(wdg(A,B)) = dot(A,dual(B))

        auto a = vec2d{5.0, 1.0};
        auto b = vec2d{3.0, 3.0};

        // fmt::println("   a                 = {}", a);
        // fmt::println("   b                 = {}", b);
        // fmt::println("   dual(b)           = {}", right_dual(b));
        // fmt::println("   dot(a, b)         = {}", dot(a, b));
        // fmt::println("   wdg(a, b)         = {}", wdg(a, b));
        // fmt::println("");
        // fmt::println("   dual(dot(a, b))   = {}", right_dual(dot(a, b)));
        // fmt::println("   wdg(a, dual(b))   = {}", wdg(a, right_dual(b)));
        // fmt::println("");
        // fmt::println("   dual(wdg(a, b))   = {}", right_dual(wdg(a, b)));
        // fmt::println("   dot(a, dual(b))   = {}", dot(a, right_dual(b)));

        CHECK(right_dual(scalar2d(dot(a, b))) == wdg(a, right_dual(b)));
        CHECK(wdg(a, right_dual(b)) == dot(a, b) * I_2d);
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

        CHECK(gr1(V1) == v1);
        CHECK(gr1(V2) == v2);
        CHECK(gr1(V3) == v3);
        CHECK(gr2(B) == wdg(v1, v2));

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

        // commutator and contraction are equal for a specific case with bivectors
        CHECK(cmt(v1, ps1) == (ps1 >> v1));
        CHECK(cmt(ps1, v1) == (v1 << ps1));
        CHECK(cmt(v1, ps1) == -cmt(ps1, v1));
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
        // n = right_dual(u) or n = cmpl(u)
        // thus, wdg(v_perp_B, n) == 0 is required
        CHECK(nrm_sq(wdg(v_perp_u, right_dual(u))) < eps);
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

        // complements are defined for basis elements only
        // => magnitude has to be covered separately for non-normalized elements

        // left complements
        CHECK(wdg(lcmpl(s1), s1) / nrm_sq(s1) == I_2d);
        CHECK(wdg(lcmpl(v1), v1) / nrm_sq(v1) == I_2d);
        CHECK(wdg(lcmpl(ps1), ps1) / nrm_sq(ps1) == I_2d);

        // right complements
        CHECK(wdg(s1, rcmpl(s1)) / nrm_sq(s1) == I_2d);
        CHECK(wdg(v1, rcmpl(v1)) / nrm_sq(v1) == I_2d);
        CHECK(wdg(ps1, rcmpl(ps1)) / nrm_sq(ps1) == I_2d);

        // correspondence of complements with geometric products:
        // left and right duals correspond to left and right complements in ega2d
        // (non-degenerate metric = identity in ega2d)
        CHECK((I_2d * rev(s1)) == left_dual(s1));
        CHECK((I_2d * rev(v1)) == left_dual(v1));
        CHECK((I_2d * rev(ps1)) == left_dual(ps1));

        CHECK((rev(s1) * I_2d) == right_dual(s1));
        CHECK((rev(v1) * I_2d) == right_dual(v1));
        CHECK((rev(ps1) * I_2d) == right_dual(ps1));
    }

    TEST_CASE("EGA2D: congruence tests")
    {
        fmt::println("EGA2D: congruence tests");

        // Test scalars
        scalar2d s1{5.0};
        scalar2d s2{-3.0}; // different sign
        scalar2d s3{2.5};  // same sign as s1
        scalar2d s_zero{0.0};

        CHECK(is_congruent2d(s1, s2) == true);         // different signs are congruent
        CHECK(is_congruent2d(s1, s3) == true);         // same signs are congruent
        CHECK(is_congruent2d(s1, s_zero) == false);    // zero vs non-zero
        CHECK(is_congruent2d(s_zero, s_zero) == true); // zero vs zero

        // Test vectors
        vec2d v1{1.0, 0.0};  // x-axis
        vec2d v2{0.0, 1.0};  // y-axis (perpendicular)
        vec2d v3{2.0, 0.0};  // parallel to v1
        vec2d v4{-1.5, 0.0}; // antiparallel to v1
        vec2d v_zero{0.0, 0.0};

        CHECK(is_congruent2d(v1, v2) == false);        // perpendicular vectors
        CHECK(is_congruent2d(v1, v3) == true);         // parallel vectors
        CHECK(is_congruent2d(v1, v4) == true);         // antiparallel vectors
        CHECK(is_congruent2d(v1, v_zero) == false);    // zero vs non-zero
        CHECK(is_congruent2d(v_zero, v_zero) == true); // zero vs zero

        // Test pseudoscalars (max grade in EGA2D)
        pscalar2d p1{1.0};
        pscalar2d p2{2.0};  // different magnitude
        pscalar2d p3{-1.0}; // opposite sign
        pscalar2d p_zero{0.0};

        CHECK(is_congruent2d(p1, p2) == true); // all non-zero pseudoscalars congruent
        CHECK(is_congruent2d(p1, p3) == true); // different signs still congruent
        CHECK(is_congruent2d(p1, p_zero) == false);    // zero vs non-zero
        CHECK(is_congruent2d(p_zero, p_zero) == true); // zero vs zero

        // fmt::println("   All EGA2D congruence tests passed");
    }


    TEST_CASE("G<2,0,0>: ega2d - trial transwedge product")
    {
        fmt::println("G<2,0,0>: ega2d - trial transwedge product");
        // fmt::println("");

        scalar2d s{3.0};
        vec2d v{1.0, 2.0};
        pscalar2d ps{4.0};

        // template:
        // tr_wdg(k,a,b) = wdg(rwdg(lcmpl(c),a),rwdg(b,right_dual(c)))
        // with c from B_k = {{1},{e1,e2},{e12}} for ega2d

        // c = e12, a = e12
        {
            auto a = e12_2d;
            auto b = s;
            auto c = e12_2d;

            CHECK(wdg(rwdg(lcmpl(c), a), rwdg(b, right_dual(c))) == 0.0);
        }
        {
            auto a = e12_2d;
            auto b = v;
            auto c = e12_2d;

            CHECK(wdg(rwdg(lcmpl(c), a), rwdg(b, right_dual(c))) == 0.0);
        }
        {
            auto a = e12_2d;
            auto b = ps;
            auto c = e12_2d;

            CHECK(wdg(rwdg(lcmpl(c), a), rwdg(b, right_dual(c))) == b);
        }

        // component extraction by using the left complement (e.g. for twdg)

        // transwedge for k = 0: wdg( rwdg(lcmpl(c),a), rwdg(b,right_dual(c)) )
        // scalar c = s
        CHECK(rwdg(lcmpl(e1_2d), s) == 0.0); // zero (left contraction)
        CHECK(rwdg(lcmpl(e2_2d), s) == 0.0); // zero (left contraction)

        CHECK(rwdg(s, rcmpl(e1_2d)) == 0.0); // zero (right contraction)
        CHECK(rwdg(s, rcmpl(e2_2d)) == 0.0); // zero (right contraction)

        CHECK(rwdg(lcmpl(e1_2d), v) == 1.0); // extract v.x - remove e1 (left contr.)
        CHECK(rwdg(lcmpl(e2_2d), v) == 2.0); // extract v.y - remove e2 (left contr.)

        CHECK(rwdg(v, rcmpl(e1_2d)) == 1.0); // extract v.x - remove e1 (right contr.)
        CHECK(rwdg(v, rcmpl(e2_2d)) == 2.0); // extract v.y - remove e2 (right contr.)

        CHECK(rwdg(lcmpl(e1_2d), ps) == -4.0 * e2_2d); // remove e1 (left contr.)
        CHECK(rwdg(lcmpl(e2_2d), ps) == 4.0 * e1_2d);  // remove e2 (left contr.)

        CHECK(rwdg(ps, rcmpl(e1_2d)) == 4.0 * e2_2d);  // remove e1 (right contr.)
        CHECK(rwdg(ps, rcmpl(e2_2d)) == -4.0 * e1_2d); // remove e2 (right contr.)


        CHECK(rwdg(scalar2d(1), I_2d) == 1); // identity for pseudoscalar
        CHECK(rwdg(I_2d, scalar2d(1)) == 1); // identity for pseudoscalar

        CHECK(rwdg(scalar2d(1), 3.0 * I_2d) == 3.0); // scalar factored out
        CHECK(rwdg(3.0 * I_2d, scalar2d(1)) == 3.0); // scalar factored out

        CHECK(wdg(rwdg(scalar2d(1), 2.0 * I_2d), rwdg(3.0 * I_2d, scalar2d(1))) == 6.0);

        // fmt::println("");
    }

    TEST_CASE("G<2,0,0>: ega2d - multivector basis trial")
    {

        fmt::println("G<2,0,0>: ega2d - multivector basis trial");

        // multivector basis components
        const std::vector<std::string> mv2d_basis_str = {"1", "e1", "e2", "e12"};

        // multivector k_basis components
        // => strings representing scalar2d, vec2d, and pscalar2d
        const std::vector<std::vector<std::string>> k_basis_str = {
            {"1"}, {"e1", "e2"}, {"e12"}};


        fmt::println("mvec_coeff = {}", mv2d_basis_str);
        fmt::println("");

        fmt::println("k_basis_str = {}", k_basis_str);
        fmt::println("");
        for (size_t k = 0; k < k_basis_str.size(); ++k) {
            fmt::println("k_basis_str[{}] = {}", k, k_basis_str.at(k));
        }
        // fmt::println("");
    }

    TEST_CASE("G<2,0,0>: ega2d - multivector variant trial")
    {

        fmt::println("G<2,0,0>: ega2d - multivector variant trial");
        fmt::println("");

        // Experimental multivector2d variant type (mv2d)
        using mv2d = std::variant<scalar2d, vec2d, pscalar2d, mvec2d_e, mvec2d>;

        // Helper function to get type name for variant
        auto get_type_name = [](const mv2d& v) -> std::string {
            return std::visit(
                [](const auto& arg) -> std::string {
                    using T = std::decay_t<decltype(arg)>;
                    if constexpr (std::is_same_v<T, scalar2d>) return "scalar2d";
                    else if constexpr (std::is_same_v<T, vec2d>) return "vec2d";
                    else if constexpr (std::is_same_v<T, pscalar2d>) return "pscalar2d";
                    else if constexpr (std::is_same_v<T, mvec2d_e>) return "mvec2d_e";
                    else if constexpr (std::is_same_v<T, mvec2d>) return "mvec2d";
                    else return "unknown";
                },
                v);
        };

        // Test creation of different types
        mv2d mv_scalar = scalar2d{2.5};
        mv2d mv_vector = vec2d{1.0, 2.0};
        mv2d mv_pscalar = pscalar2d{3.14};
        mv2d mv_even = mvec2d_e{1.0, 2.0};         // scalar + bivector
        mv2d mv_full = mvec2d{1.0, 2.0, 3.0, 4.0}; // scalar + vector + bivector

        mv2d mv_even2 = scalar2d{7.0} + pscalar2d{12.0};

        fmt::println("mv2d variant type tests:");
        fmt::println("mv_scalar: {} (type: {})", mv_scalar.index(),
                     get_type_name(mv_scalar));
        fmt::println("mv_vector: {} (type: {})", mv_vector.index(),
                     get_type_name(mv_vector));
        fmt::println("mv_pscalar: {} (type: {})", mv_pscalar.index(),
                     get_type_name(mv_pscalar));
        fmt::println("mv_even: {} (type: {})", mv_even.index(), get_type_name(mv_even));
        fmt::println("mv_full: {} (type: {})", mv_full.index(), get_type_name(mv_full));
        fmt::println("");

        // Test visitor pattern for operations
        auto print_value = [](const mv2d& v) {
            std::visit([](const auto& arg) { fmt::println("Value: {}", arg); }, v);
        };

        fmt::println("Values:");
        print_value(mv_scalar);
        print_value(mv_vector);
        print_value(mv_pscalar);
        print_value(mv_even);
        print_value(mv_full);
        fmt::println("");

        print_value(mv_even2);
        // fmt::println("");

        // Test type checking with std::holds_alternative
        CHECK(std::holds_alternative<scalar2d>(mv_scalar));
        CHECK(std::holds_alternative<vec2d>(mv_vector));
        CHECK(std::holds_alternative<pscalar2d>(mv_pscalar));
        CHECK(std::holds_alternative<mvec2d_e>(mv_even));
        CHECK(std::holds_alternative<mvec2d>(mv_full));

        // Test std::get access
        CHECK(std::get<scalar2d>(mv_scalar) == doctest::Approx(2.5).epsilon(eps));
        CHECK(std::get<vec2d>(mv_vector).x == doctest::Approx(1.0).epsilon(eps));
        CHECK(std::get<vec2d>(mv_vector).y == doctest::Approx(2.0).epsilon(eps));
        CHECK(std::get<pscalar2d>(mv_pscalar) == doctest::Approx(3.14).epsilon(eps));

        // Test mvec2d_e and mvec2d component access
        CHECK(std::get<mvec2d_e>(mv_even).c0 ==
              doctest::Approx(1.0).epsilon(eps)); // scalar part
        CHECK(std::get<mvec2d_e>(mv_even).c1 ==
              doctest::Approx(2.0).epsilon(eps)); // bivector part
        CHECK(std::get<mvec2d>(mv_full).c0 ==
              doctest::Approx(1.0).epsilon(eps)); // scalar
        CHECK(std::get<mvec2d>(mv_full).c1 == doctest::Approx(2.0).epsilon(eps)); // e1
        CHECK(std::get<mvec2d>(mv_full).c2 == doctest::Approx(3.0).epsilon(eps)); // e2
        CHECK(std::get<mvec2d>(mv_full).c3 == doctest::Approx(4.0).epsilon(eps)); // e12

        fmt::println("mv2d variant type trial implementation complete");
        fmt::println("");
    }

} // EGA 2D Tests
