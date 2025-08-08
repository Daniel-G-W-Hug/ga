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
// Euclideon geometric algebra of three-dimensional space: ega3d
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("EGA 3D Tests")
{

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

        // check inverses - odd grade multivector
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
        CHECK(v7 == normalize(vec3d(0.0, -1.0, 0.0)));
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

    TEST_CASE("Vec3d: operations - rotations")
    {
        fmt::println("Vec3d: operations - rotations");

        // fmt::println("");

        double phi = deg2rad(30.0);
        auto c = vec3d(std::cos(phi), std::sin(phi), 0.0);
        auto d = rotate(e1_3d, rotor(e12_3d, phi));
        auto cm = vec3d(std::cos(-phi), std::sin(-phi), 0.0);
        auto dm = rotate(e1_3d, rotor(e12_3d, -phi));

        // fmt::println("c = {}", c);
        // fmt::println("d = {}", d);
        // fmt::println("cm = {}", cm);
        // fmt::println("dm = {}", dm);

        CHECK(c == d);
        CHECK(cm == dm);

        // fmt::println("");
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

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}, angle = {: .2f}°", v1,
        // nrm(v1),
        //              rad2deg(angle(e1_3d, v1)));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v2,
        // nrm(v2),
        //              rad2deg(angle(e1_3d, v2)));
        // fmt::println("v3  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v3,
        // nrm(v3),
        //              rad2deg(angle(e1_3d, v3)));
        // fmt::println("");
        // fmt::println("og12[0]  = {: .4f}, nrm(og12[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              og12[0], nrm(og12[0]), rad2deg(angle(e1_3d, og12[0])));
        // fmt::println("og12[1]  = {: .4f}, nrm(og12[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              og12[1], nrm(og12[1]), rad2deg(angle(e1_3d, og12[1])));
        // fmt::println("on12[0]  = {: .4f}, nrm(on12[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              on12[0], nrm(on12[0]), rad2deg(angle(e1_3d, on12[0])));
        // fmt::println("on12[1]  = {: .4f}, nrm(on12[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              on12[1], nrm(on12[1]), rad2deg(angle(e1_3d, on12[1])));
        // fmt::println("");
        // fmt::println("og13[0]  = {: .4f}, nrm(og13[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              og13[0], nrm(og13[0]), rad2deg(angle(e1_3d, og13[0])));
        // fmt::println("og13[1]  = {: .4f}, nrm(og13[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              og13[1], nrm(og13[1]), rad2deg(angle(e1_3d, og13[1])));
        // fmt::println("on13[0]  = {: .4f}, nrm(on13[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              on13[0], nrm(on13[0]), rad2deg(angle(e1_3d, on13[0])));
        // fmt::println("on13[1]  = {: .4f}, nrm(on13[1]) = {: .4f}, angle = {:
        // .2f}°",
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

        // fmt::println("v1  = {: .4f}, nrm(v1) = {: .4f}, angle = {: .2f}°", v1,
        // nrm(v1),
        //              rad2deg(angle(e1_3d, v1)));
        // fmt::println("v2  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v2,
        // nrm(v2),
        //              rad2deg(angle(e1_3d, v2)));
        // fmt::println("v3  = {: .4f}, nrm(v2) = {: .4f}, angle = {: .2f}°", v3,
        // nrm(v3),
        //              rad2deg(angle(e1_3d, v3)));
        // fmt::println("");
        // fmt::println("og123[0]  = {: .4f}, nrm(og123[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              og123[0], nrm(og123[0]), rad2deg(angle(e1_3d, og123[0])));
        // fmt::println("og123[1]  = {: .4f}, nrm(og123[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              og123[1], nrm(og123[1]), rad2deg(angle(e1_3d, og123[1])));
        // fmt::println("og123[2]  = {: .4f}, nrm(og123[2]) = {: .4f}, angle = {:
        // .2f}°",
        //              og123[2], nrm(og123[2]), rad2deg(angle(e1_3d, og123[2])));
        // fmt::println("");
        // fmt::println("on123[0]  = {: .4f}, nrm(on123[0]) = {: .4f}, angle = {:
        // .2f}°",
        //              on123[0], nrm(on123[0]), rad2deg(angle(e1_3d, on123[0])));
        // fmt::println("on123[1]  = {: .4f}, nrm(on123[1]) = {: .4f}, angle = {:
        // .2f}°",
        //              on123[1], nrm(on123[1]), rad2deg(angle(e1_3d, on123[1])));
        // fmt::println("on123[2]  = {: .4f}, nrm(on123[2]) = {: .4f}, angle = {:
        // .2f}°",
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
        bivec3d B{-3.4, 2.1, 0.7};
        pscalar3d ps{1.0};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3d{sd};
        auto t = scalar3d{st};

        CHECK(wdg(wdg(v1, v2), v3) == wdg(v1, wdg(v2, v3)));  // wdg is associative
        CHECK(wdg(v1, v1) == bivec3d{});                      // wdg=0 for collin. vectors
        CHECK(wdg(v1, v2) == -wdg(v2, v1));                   // anticommutative for vect.
        CHECK(wdg(v1, v2 + v3) == wdg(v1, v2) + wdg(v1, v3)); // wdg distributes over add.
        CHECK(wdg(v1 + v2, v3) == wdg(v1, v3) + wdg(v2, v3)); // wdg distributes over add.
        CHECK(wdg(sd * v1, v2) == wdg(v1, sd * v2)); // scalars can be factored out of wdg
        CHECK(wdg(sd * v1, v2) == sd * wdg(v1, v2)); // scalars can be factored out of wdg
        CHECK(wdg(s, t) == wdg(t, s)); // wdg between scalars equivalent to scalar mult.
        CHECK(wdg(s, t) == s * t);
        CHECK(wdg(s, v1) == wdg(v1, s));       // wdg between scalar and vector
        CHECK(wdg(s, B) == wdg(B, s));         // wdg between scalar and bivector
        CHECK(wdg(scalar3d{1.0}, I_3d) == ps); // wdg between scalar and trivector
        CHECK(wdg(I_3d, scalar3d{1.0}) == ps); // (=pseudoscalar); 1.0 is neutral element
                                               // for the wedge product
        CHECK(wdg(s, v1) == sd * v1);          // wdg between scalar and vector
    }

    TEST_CASE("MVec3d: regressive wedge product - basic properties")
    {
        fmt::println("MVec3d: regressive wedge product - basic properties");

        bivec3d B1{1.0, 2.0, 1.0};
        bivec3d B2{0.5, 3.0, 2.0};
        bivec3d B3{-2.0, 6.0, 3.0};
        vec3d v1{-3.4, 2.1, 0.7};
        vec3d v2{1, 2, 3};
        vec3d z = {-2, 3, 7};
        scalar3d sc{2.0};

        double sd = 2.3;
        double st = -5.1;
        auto s = pscalar3d{sd};
        auto t = pscalar3d{st};

        CHECK(rwdg(rwdg(B1, B2), B3) == rwdg(B1, rwdg(B2, B3))); // rwdg is associative
        CHECK(rwdg(B1, B1) == vec3d{});       // wdg=0 for coplanar bivectors
        CHECK(rwdg(B1, B2) == -rwdg(B2, B1)); // rwdg is anticommutative for bivectors
        CHECK(rwdg(B1, B2 + B3) ==
              rwdg(B1, B2) + rwdg(B1, B3)); // rwdg distributes over add.
        CHECK(rwdg(B1 + B2, B3) ==
              rwdg(B1, B3) + rwdg(B2, B3)); // rwdg distributes over add.
        CHECK(rwdg(sd * B1, B2) ==
              rwdg(B1, sd * B2)); // scalars can be factored out of wdg
        CHECK(rwdg(sd * B1, B2) ==
              sd * rwdg(B1, B2));        // scalars can be factored out of wdg
        CHECK(rwdg(s, t) == rwdg(t, s)); // rwdg between pscalars equivalent to smult.
        CHECK(rwdg(s, t) == double(s) * double(t));
        CHECK(rwdg(s, B1) == rwdg(B1, s)); // rwdg between pscalar and bivector
        CHECK(rwdg(s, v1) == rwdg(v1, s)); // rwdg between pscalar and vector
        CHECK(rwdg(sc, I_3d) == sc);       // rwdg between pscalar and scalar
        CHECK(rwdg(I_3d, sc) == sc);       // (=pseudoscalar); I_3d is neutral element
                                           // for the wedge product
        CHECK(rwdg(s, B1) == sd * B1);     // rwdg between pseudoscalar and vector
        // (I_3d is neutral element: s = sd * I_3d)

        // congruence
        CHECK(is_congruent3d(v1, 2.3 * v1));
        CHECK(is_congruent3d(v1, -12.3 * v1));
        CHECK(is_congruent3d(wdg(v1, v2), wdg(v2, v1)));
        CHECK_FALSE(is_congruent3d(v1, v2));

        // common factor axiom
        auto lhs = rwdg(wdg(v1, z), wdg(v2, z));
        auto rhs = rwdg(wdg(wdg(v1, v2), z), z);
        CHECK(lhs.x == doctest::Approx(rhs.x).epsilon(eps));
        CHECK(lhs.y == doctest::Approx(rhs.y).epsilon(eps));
        CHECK(lhs.z == doctest::Approx(rhs.z).epsilon(eps));
        CHECK(is_congruent3d(rwdg(wdg(v1, z), wdg(v2, z)), z));
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
        // assumption: the geometric product of a fully populated multivector serves
        // as
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
        // fmt::println("   gpr12_m_sym  = 0.5*(gpr12_d + gpr21_d) = {}",
        // gpr12_m_sym); fmt::println("   gpr12_m_asym = 0.5*(gpr12_m - gpr21_m) =
        // {}", gpr12_m_asym); fmt::println(""); fmt::println("   gpr12_d = b1 * b2 =
        // {} ", gpr12_d); fmt::println("   gpr21_d = b2 * b1) = {} ", gpr21_d);
        // fmt::println("   gpr12_d_sym  = 0.5*(gpr12_d + gpr21_d) = {}",
        // gpr12_d_sym); fmt::println("   gpr12_d_asym = 0.5*(gpr12_d - gpr21_d) =
        // {}", gpr12_d_asym); fmt::println("");

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

    TEST_CASE("MVec3d_E/_U: modelling even and odd parts of 3d algebra - basics")
    {
        fmt::println("MVec3d_E/_U: modelling even and odd parts of 3d algebra - basics");

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

        // full 3d multivector - odd content
        mvec3d vm_odd{0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 300.0};
        // odd grade 3d multivector
        mvec3d_u vm_U{1.0, 2.0, 3.0, 300.0};

        ////////////////////////////////////////////////////////////////////////////////
        // duality (as defined in Lengyel, "PGA illuminated"):
        ////////////////////////////////////////////////////////////////////////////////
        //
        // dual(A) = cmpl(A) in spaces of odd dimension
        //         = rcmpl(A) in spaces of even dimension
        //

        auto vm_dual_manual = rev(vm) * Im_3d;
        auto vm_dual = dual(vm);

        auto vm_dual_even_manual = rev(vm_even) * Im_3d;
        auto vm_dual_even = dual(vm_even);

        auto vm_dual_odd_manual = rev(vm_odd) * Im_3d;
        auto vm_dual_odd = dual(vm_odd);

        // result is odd, naming chosen for consistency
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
        // fmt::println("   vm_odd              = {}", vm_odd);
        // fmt::println("   rev(vm_odd) * Im_3d = {}", vm_dual_odd_manual);
        // fmt::println("   dual(vm_odd)        = {}", vm_dual_odd);
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
        CHECK(vm_dual_odd == vm_dual_odd_manual);
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

        CHECK(gr1(V1) == v1);
        CHECK(gr1(V2) == v2);
        CHECK(gr1(V3) == v3);
        CHECK(gr2(B12) == wdg(v1, v2));

        CHECK(nrm_sq(M1) == nrm_sq(s1) + nrm_sq(v1) + nrm_sq(b1) + nrm_sq(ps1));
        CHECK(nrm_sq(M2) == nrm_sq(s2) + nrm_sq(v2) + nrm_sq(b2) + nrm_sq(ps2));
        CHECK(nrm_sq(M3) == nrm_sq(s3) + nrm_sq(v3) + nrm_sq(b3) + nrm_sq(ps3));

        CHECK(dot(M1, M2) == dot(s1, s2) + dot(v1, v2) + dot(b1, b2) + dot(ps1, ps2));
        CHECK(dot(M1, M1) == nrm_sq(M1));
        CHECK(dot(M2, M2) == nrm_sq(M2));
        CHECK(dot(M3, M3) == nrm_sq(M3));

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

        // commutator and contraction are equal for a specific case with bivectors
        CHECK(cmt(v1, b1) == (b1 >> v1));
        CHECK(cmt(b1, v1) == (v1 << b1));
        CHECK(cmt(v1, b1) == -cmt(b1, v1));
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
        //              rev(s) * I_3d, rev(e1) * I_3d, rev(e2) * I_3d, rev(e3) *
        //              I_3d);
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
        //              I_3d * rev(s), I_3d * rev(e1), I_3d * rev(e2), I_3d *
        //              rev(e3));
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
        // complements are defined for basis elements only
        // => magnitude has to be covered separately for non-normalized elements

        // left complements = complements in spaces of odd dimension
        CHECK(wdg(cmpl(s1), s1) / nrm_sq(s1) == I_3d);
        CHECK(wdg(cmpl(v1), v1) / nrm_sq(v1) == I_3d);
        CHECK(wdg(cmpl(b1), b1) / nrm_sq(b1) == I_3d);
        CHECK(wdg(cmpl(ps1), ps1) / nrm_sq(ps1) == I_3d);

        // right complements = complements in spaces of odd dimension
        CHECK(wdg(s1, cmpl(s1)) / nrm_sq(s1) == I_3d);
        CHECK(wdg(v1, cmpl(v1)) / nrm_sq(v1) == I_3d);
        CHECK(wdg(b1, cmpl(b1)) / nrm_sq(b1) == I_3d);
        CHECK(wdg(ps1, cmpl(ps1)) / nrm_sq(ps1) == I_3d);

        // correspondence of complements with geometric products:
        // duals correspond to complements in ega3d (non-degenerate metric = identity)
        CHECK((I_3d * rev(s1)) == dual(s1));
        CHECK((I_3d * rev(v1)) == dual(v1));
        CHECK((I_3d * rev(b1)) == dual(b1));
        CHECK((I_3d * rev(ps1)) == dual(ps1));

        CHECK((rev(s1) * I_3d) == dual(s1));
        CHECK((rev(v1) * I_3d) == dual(v1));
        CHECK((rev(b1) * I_3d) == dual(b1));
        CHECK((rev(ps1) * I_3d) == dual(ps1));

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

    TEST_CASE("EGA3D: congruence tests")
    {
        fmt::println("EGA3D: congruence tests");

        // Test scalars
        scalar3d s1{5.0};
        scalar3d s2{-3.0}; // different sign
        scalar3d s3{2.5};  // same sign as s1
        scalar3d s_zero{0.0};

        CHECK(is_congruent3d(s1, s2) == true);         // different signs are congruent
        CHECK(is_congruent3d(s1, s3) == true);         // same signs are congruent
        CHECK(is_congruent3d(s1, s_zero) == false);    // zero vs non-zero
        CHECK(is_congruent3d(s_zero, s_zero) == true); // zero vs zero

        // Test vectors
        vec3d v1{1.0, 0.0, 0.0};  // x-axis
        vec3d v2{0.0, 1.0, 0.0};  // y-axis (perpendicular)
        vec3d v3{2.0, 0.0, 0.0};  // parallel to v1
        vec3d v4{-1.5, 0.0, 0.0}; // antiparallel to v1
        vec3d v_zero{0.0, 0.0, 0.0};

        CHECK(is_congruent3d(v1, v2) == false);        // perpendicular vectors
        CHECK(is_congruent3d(v1, v3) == true);         // parallel vectors
        CHECK(is_congruent3d(v1, v4) == true);         // antiparallel vectors
        CHECK(is_congruent3d(v1, v_zero) == false);    // zero vs non-zero
        CHECK(is_congruent3d(v_zero, v_zero) == true); // zero vs zero

        // Test bivectors
        bivec3d b1{1.0, 0.0, 0.0};  // e23 plane
        bivec3d b2{0.0, 1.0, 0.0};  // e31 plane (different)
        bivec3d b3{2.0, 0.0, 0.0};  // e23 plane (parallel to b1)
        bivec3d b4{-1.5, 0.0, 0.0}; // e23 plane (antiparallel to b1)
        bivec3d b_zero{0.0, 0.0, 0.0};

        CHECK(is_congruent3d(b1, b2) == false);        // different planes
        CHECK(is_congruent3d(b1, b3) == true);         // same plane
        CHECK(is_congruent3d(b1, b4) == true);         // same plane, opposite orientation
        CHECK(is_congruent3d(b1, b_zero) == false);    // zero vs non-zero
        CHECK(is_congruent3d(b_zero, b_zero) == true); // zero vs zero

        // Test pseudoscalars
        pscalar3d p1{1.0};
        pscalar3d p2{2.0};  // different magnitude
        pscalar3d p3{-1.0}; // opposite sign
        pscalar3d p_zero{0.0};

        CHECK(is_congruent3d(p1, p2) == true); // all non-zero pseudoscalars congruent
        CHECK(is_congruent3d(p1, p3) == true); // different signs still congruent
        CHECK(is_congruent3d(p1, p_zero) == false);    // zero vs non-zero
        CHECK(is_congruent3d(p_zero, p_zero) == true); // zero vs zero

        fmt::println("   All EGA3D congruence tests passed");
    }

} // EGA 3D Tests
