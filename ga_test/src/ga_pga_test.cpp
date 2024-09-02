// author: Daniel Hug, 2024

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_pga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)

TEST_SUITE("Projective Geometric Algebra (PGA)")
{

    TEST_CASE("algebra<2,0,1> - pga_2d")
    {
        fmt::println("");
        fmt::println("algebra<2, 0, 1> - pga_2d:");
        // 2d projective geometric algebra (models 2d Euclidean algebra)
        // embedded in 3d representational space
        const algebra<2, 0, 1> alg;
        CHECK(alg.p() == 2);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 1);
        CHECK(alg.dim_space() == 3);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 8);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 4); // == dim_space + 1
        fmt::println("   pga_2d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   pga_2d: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("algebra<3, 0, 1> pga_3d")
    {
        fmt::println("");
        fmt::println("algebra<3, 0, 1> - pga_3d:");
        // 3d projective geometric algebra  (models 3d Euclidean algebra)
        // embedded in 4d representational space
        const algebra<3, 0, 1> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 1);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   pga_3d: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   pga_3d: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("pga_2d<2,0,1> - defining basic types and ctor checks")
    {
        fmt::println("pga_2d<2,0,1>: defining basic types and ctor checks");

        auto mv1 = mvec2dp{scalar(5.0)};
        auto mv2 = mvec2dp{vec2dp{1.0, 2.0, 1.0}};
        auto mv3 = mvec2dp{bivec2dp{-1.0, 2.0, 1.0}};
        auto mv4 = mvec2dp{pscalar2dp(-5.0)};
        auto mv5a = mvec2dp_e{scalar(5.0), bivec2dp{-1.0, 2.0, 1.0}};
        auto mv5 = mvec2dp{mv5a};
        auto mv6a = mvec2dp_u{vec2dp{1.0, 2.0, 1.0}, pscalar2dp{-5.0}};
        auto mv6 = mvec2dp{mv6a};

        // fmt::println("   mv1  = {}", mv1);
        // fmt::println("   mv2  = {}", mv2);
        // fmt::println("   mv3  = {}", mv3);
        // fmt::println("   mv4  = {}", mv4);
        // fmt::println("   mv5a = {}", mv5a);
        // fmt::println("   mv5  = {}", mv5);
        // fmt::println("   mv6a = {}", mv6a);
        // fmt::println("   mv6  = {}", mv6);

        CHECK(gr0(mv1) == scalar(5.0));
        CHECK(gr1(mv2) == vec2dp{1.0, 2.0, 1.0});
        CHECK(gr2(mv3) == bivec2dp{-1.0, 2.0, 1.0});
        CHECK(gr3(mv4) == pscalar2dp{-5.0});
        CHECK(gr0(mv5a) == scalar(5.0));
        CHECK(gr2(mv5a) == bivec2dp{-1.0, 2.0, 1.0});
        CHECK(gr0(mv5) == scalar(5.0));
        CHECK(gr2(mv5) == bivec2dp{-1.0, 2.0, 1.0});
        CHECK(gr1(mv6a) == vec2dp{1.0, 2.0, 1.0});
        CHECK(gr3(mv6a) == pscalar2dp(-5.0));
        CHECK(gr1(mv6) == vec2dp{1.0, 2.0, 1.0});
        CHECK(gr3(mv6) == pscalar2dp(-5.0));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec2dp<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2dp: default init")
    {
        fmt::println("Vec2dp: default init");
        vec2dp v;
        // fmt::println("   v = {}", v);
        CHECK(v.x == 0.0);
        CHECK(v.y == 0.0);
        CHECK(v.z == 0.0);
    }
    TEST_CASE("Vec2dp: with curly braced intializer")
    {
        fmt::println("Vec2dp: with curly braced intializer");
        vec2dp v{0.0, 0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(v.x == 0.0);
        CHECK(v.y == 0.0);
        CHECK(v.z == 0.0);
    }
    TEST_CASE("Vec2dp: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("Vec2dp: cp ctor & cp assign incl. type deduction");
        vec2dp v1{1.0, 2.0, 3.0}; // init with double (type deduction)
        vec2dp v2{v1};            // cp ctor
        vec2dp v3 = v2;           // cp assign
        vec2dp v4 = -v2;          // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(v1.x == 1.0);
        CHECK(v1.y == 2.0);
        CHECK(v1.z == 3.0);
        CHECK(v2.x == 1.0);
        CHECK(v2.y == 2.0);
        CHECK(v2.z == 3.0);
        CHECK(v3.x == 1.0);
        CHECK(v3.y == 2.0);
        CHECK(v3.z == 3.0);
        CHECK(v4 == -v2);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }


    TEST_CASE("Vec2dp: fmt & cout printing")
    {
        fmt::println("Vec2dp: fmt & cout printing");

        vec2dp pf{1.0f, 2.00001f, 3.0f};
        vec2dp pd{1.0, 2.00001, 3.0};

        // std::cout << "       cout: pf = " << pf << std::endl;
        fmt::println("       fmt:  pf = {}", pf);
        fmt::println("       fmt:  pf = {:.8f}", pf);

        // std::cout << "       cout: pd = " << pd << std::endl;
        fmt::println("       fmt:  pd = {}", pd);
        fmt::println("       fmt:  pd = {:.8f}", pd);

        std::vector<Vec2dp<double>> vp1{{1.0, 1.0, 1.0}, {1.5, 2.0, 3.0}};
        fmt::println("       fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("       fmt: vp1 = {:.e}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("Vec2dp: comparison float")
    {
        fmt::println("Vec2dp: comparison float");

        Vec2dp<float> v1f{1.0f, 2.0f, 3.0f};
        Vec2dp<float> v2f{2.0f, 4.0f, 3.0f};
        Vec2dp<float> v3f{1.0f, 2.0000001f, 3.0f};
        Vec2dp<float> v4f{v1f};

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

    TEST_CASE("Vec2dp: comparison double")
    {
        fmt::println("Vec2dp: comparison double");

        Vec2dp<double> v1d{1.0, 2.0, 3.0};
        Vec2dp<double> v2d{2.0, 4.0, 3.0};
        Vec2dp<double> v3d{1.0, 2.0000000000000001, 3.0};
        Vec2dp<double> v4d{v1d};

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

    TEST_CASE("Vec2dp: vector space and linearity tests")
    {
        fmt::println("Vec2dp: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        vec2dp p0;
        vec2dp p1{1.0, 2.0, 3.0};
        vec2dp p2{2.0, 4.0, 6.0};
        vec2dp p3{3.0, 6.0, 9.0};
        vec2dp p4 = -p1; // assignment using unary minus
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

    TEST_CASE("Vec2dp: inner product properties")
    {
        fmt::println("Vec2dp: inner product properties");

        double a = 2.35;
        vec2dp u{1.0, 2.0, 1.0};
        vec2dp v{-0.5, 3.0, 0.5};
        vec2dp w{3.0, 6.0, -3.0};

        CHECK(dot(a * u, v) == a * dot(u, v));
        CHECK(dot(u + v, w) == dot(u, w) + dot(v, w));
        CHECK(dot(u, v) == dot(v, u));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec2dp<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2dp: operations - norm, inverse, dot")
    {
        fmt::println("Vec2dp: operations - norm, inverse, dot");

        vec2dp v1{2.0, 1.0, 2.0};
        vec2dp v2{normalize(v1)};

        vec2dp v3{2.0, 6.0, -4.0};
        vec2dp v4{inv(v3)};

        auto v3m = mvec2dp{v3};
        auto v4m = mvec2dp{v4};

        fmt::println("v1 = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        fmt::println("v2 = normalize(v1) = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        fmt::println("v3 = {: .4f}, nrm(v1) = {: .4f}", v3, nrm(v3));
        fmt::println(
            "v4 = inv(v3) = {: .4f}, nrm(v3) = {: .4f}, nrm(v3)*nrm(v4) = {: .4f}", v4,
            nrm(v4), nrm(v3) * nrm(v4));
        fmt::println("");
        fmt::println("v3m = {}, nrm(v3m) = {}", v3m, nrm(v3m));
        fmt::println("v4m = {}, nrm(v4m) = {}", v4m, nrm(v4m));
        fmt::println("0.5*(v3m*v4m + v4m*v3m) = {}", 0.5 * (v3m * v4m + v4m * v3m));
        fmt::println("0.5*(v3m*v4m - v4m*v3m) = {}", 0.5 * (v3m * v4m - v4m * v3m));

        CHECK(std::abs(magn_sq(v1) - 9.0) < eps);
        CHECK(std::abs(nrm_sq(v1) - 5.0) < eps);
        CHECK(std::abs(magn_sq(v2) - 1.0) < eps);
        CHECK(std::abs(nrm_sq(v3) - 40.0) < eps);
        CHECK(std::abs(dot(v4, v3) - 1.0) < eps);
    }

    TEST_CASE("Vec2dp: operations - angle I")
    {
        fmt::println("Vec2dp: operations - angle");

        vec2dp v1{1.0, 0.0, 0.0};
        vec2dp v2{normalize(vec2dp(1.0, 1.0, 0.0))};
        vec2dp v3{0.0, 1.0, 0.0};
        vec2dp v4{normalize(vec2dp(-1.0, 1.0, 0.0))};
        vec2dp v5{-1.0, 0.0, 0.0};
        vec2dp v6{normalize(vec2dp(-1.0, -1.0, 0.0))};
        vec2dp v7{0.0, -1.0, 0.0};
        vec2dp v8{normalize(vec2dp(1.0, -1.0, 0.0))};

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
        CHECK(v6 == normalize(vec2dp(-1.0, -1.0, 0.0)));
        CHECK(v8 == normalize(vec2dp(1.0, -1.0, 0.0)));
    }

    TEST_CASE("Vec2dp: operations - angle II")
    {
        fmt::println("Vec2dp: operations - angle II");

        std::vector<std::tuple<double, Vec2dp<double>>> v1;
        std::vector<std::tuple<double, Vec2dp<double>>> v2;
        std::vector<std::tuple<double, Vec2dp<double>>> v3;

        // only positive angles are easy to implement vs. the 2d case

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2dp<double>(std::cos(phi), std::sin(phi), 0.0);
            v1.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e1_2dp, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2dp<double>(std::cos(phi + pi / 2), std::sin(phi + pi / 2), 0.0);
            v2.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e2_2dp, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2dp<double>(std::cos(phi + pi / 4), std::sin(phi + pi / 4), 0.0);
            v3.push_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .3f},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e1_2dp + e2_2dp, c));
        }
        // fmt::println("");

        for (auto const& [phi, c] : v1) {
            CHECK(std::abs(phi - angle(e1_2dp, c)) < eps);
        }
        for (auto const& [phi, c] : v2) {
            CHECK(std::abs(phi - angle(e2_2dp, c)) < eps);
        }
        auto ref_vec = normalize(e1_2dp + e2_2dp);
        for (auto const& [phi, c] : v3) {
            CHECK(std::abs(phi - angle(ref_vec, c)) < eps);
        }
    }

    TEST_CASE("Vec2dp: operations - wedge")
    {
        fmt::println("Vec2dp: operations - wedge");

        vec2dp v1{1.0, 0.0, 0.0};
        vec2dp v2{normalize(vec2dp(1.0, 1.0, 0.0))};
        vec2dp v3{0.0, 1.0, 0.0};
        vec2dp v4{normalize(vec2dp(-1.0, 1.0, 0.0))};
        vec2dp v5{-1.0, 0.0, 0.0};
        vec2dp v6{normalize(vec2dp(-1.0, -1.0, 0.0))};
        vec2dp v7{0.0, -1.0, 0.0};
        vec2dp v8{normalize(vec2dp(1.0, -1.0, 0.0))};

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

    TEST_CASE("Vec2dp: operations - project / reject / reflect (vector - vector)")
    {
        fmt::println("Vec2dp: operations - project / reject / reflect (vector - vector)");

        vec2dp v1{5.0, 1.0, 1.0};
        vec2dp v2{2.0, 2.0, 1.0};

        vec2dp v2u = normalize(v2);

        vec2dp v3{project_onto(v1, v2)};
        vec2dp v4{reject_from(v1, v2)};
        vec2dp v5{v3 + v4};

        vec2dp v6{project_onto_unitized(v1, v2u)};
        vec2dp v7{reject_from_unitized(v1, v2u)};
        vec2dp v8{v6 + v7};

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
        // fmt::println("v6 = project_onto_unitized(v1, v2u) = {: .4f},"
        //              " nrm(v6) = {: .4f}",
        //              v6, nrm(v6));
        // fmt::println("v7 = reject_from_unitized(v1, v2u)  = {: .4f},"
        //              " nrm(v7) = {: .4f}",
        //              v7, nrm(v7));
        // fmt::println("v8 = v6 + v7                        = {: .4f},"
        //              " nrm(v8) = {: .4f}",
        //              v8, nrm(v8));
        // // this helps to understand, why the dot-product is sufficient
        // auto w = wdg(v1, v2);
        // auto i = inv(v2);
        // fmt::println("wdg(v1,v2)         = {: .4f}", w);
        // fmt::println("inv(v2)            = {: .4f}", i);
        // fmt::println("wdg(v1,v2)*inv(v2) = {: .4f}", w * i);
        // fmt::println("");

        vec2dp v{4.0, 1.0, 1.0};
        vec2dp b{e2_2dp};
        auto B = bivec2dp{e12_2dp};

        // auto UB = bivec2dp{e23_2dp + e12_2dp};
        // fmt::println("v   = {}", v);
        // fmt::println("b   = {}", b);
        // fmt::println("B   = {}", B);
        // fmt::println("UB  = {}", UB);
        // fmt::println("reflect_on_vec(v,b)  = {}", reflect_on_vec(v, b));
        // fmt::println("reflect_on_hyp(v,e3_2dp)  = {}", reflect_on_hyp(v, e3_2dp));
        // fmt::println("reflect_on(v,B)  = {}", reflect_on(v, B));
        // fmt::println("reflect_on(UB,B) = {}", reflect_on(UB, B));
        // fmt::println("");

        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);
        CHECK(v6 + v7 == v8);
        CHECK(v8 == v1);

        // just to suppress unused variable warnings
        CHECK(b == e2_2dp);

        CHECK(reflect_on_hyp(v, e3_2dp) == reflect_on(v, B));

        // checking time required
        //
        // auto start = std::chrono::system_clock::now();
        // for (size_t i = 0; i < 10000000; ++i) {
        //     vec2dp v = reject_from(v1, v2);
        // }
        // auto end = std::chrono::system_clock::now();
        // auto elapsed =
        // std::chrono::duration_cast<std::chrono::milliseconds>(end
        // -
        // start); fmt::println("The measurement took {}", elapsed);
    }

    TEST_CASE("Vec2dp: operations - project / reject / reflect (vector - bivector)")
    {
        fmt::println(
            "Vec2dp: operations - project / reject / reflect (vector - bivector)");

        vec2dp v1{5.0, 3.0, 1.0};
        bivec2dp v2 = wdg(vec2dp{0.0, 0.0, 2.0}, vec2dp{2.0, 0.0, 0.0});

        vec2dp v3{project_onto(v1, v2)};
        vec2dp v4{reject_from(v1, v2)};
        vec2dp v5{v3 + v4};

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

        // from Macdonald, "Linear and Geometric Algebra", Exercise 7.14, p. 129
        CHECK(reflect_on(wdg(e1_2dp + e3_2dp, e2_2dp), e12_2dp) ==
              wdg(e1_2dp - e3_2dp, e2_2dp));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // MVec2dp<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2dp: default init")
    {
        fmt::println("MVec2dp: default init");
        // default initialization
        mvec2dp v;
        // fmt::println("   v = {}", v);
        CHECK(v.c0 == 0.0);
        CHECK(v.c1 == 0.0);
        CHECK(v.c2 == 0.0);
        CHECK(v.c3 == 0.0);
        CHECK(v.c4 == 0.0);
        CHECK(v.c5 == 0.0);
        CHECK(v.c6 == 0.0);
        CHECK(v.c7 == 0.0);
    }
    TEST_CASE("MVec2dp: with curly braced intializer")
    {
        fmt::println("MVec2dp: with curly braced intializer");
        // default initialization
        mvec2dp v{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        // fmt::println("   v = {}", v);
        CHECK(v.c0 == 0.0);
        CHECK(v.c1 == 1.0);
        CHECK(v.c2 == 2.0);
        CHECK(v.c3 == 3.0);
        CHECK(v.c4 == 23.0);
        CHECK(v.c5 == 31.0);
        CHECK(v.c6 == 12.0);
        CHECK(v.c7 == 123.0);
    }

    TEST_CASE("MVec2dp: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("MVec2dp: cp ctor & cp assign incl. type deduction");
        // default initialization
        mvec2dp v1{0.0,  1.0,  2.0,  3.0,
                   23.0, 31.0, 12.0, 123.0}; // init with double (type deduction)
        mvec2dp v2{v1};                      // cp ctor
        mvec2dp v3 = v2;                     // cp assign
        mvec2dp v4 = -v3;                    // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(v2.c0 == 0.0);
        CHECK(v2.c1 == 1.0);
        CHECK(v2.c2 == 2.0);
        CHECK(v2.c3 == 3.0);
        CHECK(v2.c4 == 23.0);
        CHECK(v2.c5 == 31.0);
        CHECK(v2.c6 == 12.0);
        CHECK(v2.c7 == 123.0);

        CHECK(v3.c0 == 0.0);
        CHECK(v3.c1 == 1.0);
        CHECK(v3.c2 == 2.0);
        CHECK(v3.c3 == 3.0);
        CHECK(v3.c4 == 23.0);
        CHECK(v3.c5 == 31.0);
        CHECK(v3.c6 == 12.0);
        CHECK(v3.c7 == 123.0);

        CHECK(v4 == -v3);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }

    TEST_CASE("MVec2dp: defining basic types and ctor checks")
    {
        fmt::println("MVec2dp: defining basic types and ctor checks");

        auto mv1 = mvec2dp{scalar(5.0)};
        auto mv2 = mvec2dp{vec2dp{1.0, 2.0, 1.0}};
        auto mv3 = mvec2dp{bivec2dp{-1.0, 2.0, 1.0}};
        auto mv4 = mvec2dp{pscalar2dp(-5.0)};
        auto mv5a = mvec2dp_e{scalar(5.0), bivec2dp{-1.0, 2.0, 1.0}};
        auto mv5 = mvec2dp{mv5a};
        auto mv6a = mvec2dp_u{vec2dp{1.0, 2.0, 1.0}, pscalar2dp{-5.0}};
        auto mv6 = mvec2dp{mv6a};

        // fmt::println("   mv1  = {}", mv1);
        // fmt::println("   mv2  = {}", mv2);
        // fmt::println("   mv3  = {}", mv3);
        // fmt::println("   mv4  = {}", mv4);
        // fmt::println("   mv5a = {}", mv5a);
        // fmt::println("   mv5  = {}", mv5);
        // fmt::println("   mv6a = {}", mv6a);
        // fmt::println("   mv6  = {}", mv6);

        CHECK(gr0(mv1) == scalar(5.0));
        CHECK(gr1(mv2) == vec2dp{1.0, 2.0, 1.0});
        CHECK(gr2(mv3) == bivec2dp{-1.0, 2.0, 1.0});
        CHECK(gr3(mv4) == pscalar2dp{-5.0});
        CHECK(gr0(mv5a) == scalar(5.0));
        CHECK(gr2(mv5a) == bivec2dp{-1.0, 2.0, 1.0});
        CHECK(gr0(mv5) == scalar(5.0));
        CHECK(gr2(mv5) == bivec2dp{-1.0, 2.0, 1.0});
        CHECK(gr1(mv6a) == vec2dp{1.0, 2.0, 1.0});
        CHECK(gr3(mv6a) == pscalar2dp(-5.0));
        CHECK(gr1(mv6) == vec2dp{1.0, 2.0, 1.0});
        CHECK(gr3(mv6) == pscalar2dp(-5.0));
    }

    TEST_CASE("MVec2dp: fmt & cout printing")
    {
        fmt::println("MVec2dp: fmt & cout printing");

        mvec2dp pf{1.0f, 2.00001f, 0.0f, 3.0f, 1.0f, 2.00001f, 0.0f, 3.0f};
        mvec2dp pd{1.0, 2.00001, 0.0, 3.0, 1.0, 2.00001, 0.0, 3.0};

        std::cout << "    cout: pf = " << pf << std::endl;
        fmt::println("    fmt:  pf = {}", pf);
        fmt::println("    fmt:  pf = {:.8f}", pf);

        std::cout << "    cout: pd = " << pd << std::endl;
        fmt::println("    fmt:  pd = {}", pd);
        fmt::println("    fmt:  pd = {:.8f}", pd);

        std::vector<MVec2dp<double>> vp1{{1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0},
                                         {0.5, 1.5, 2.0, 2.5, 1.0, 1.0, 1.0, 2.0}};
        fmt::println("    fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("    fmt: vp1 = {:.e}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("MVec2dp: vector space and linearity tests")
    {
        fmt::println("MVec2dp: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        mvec2dp p0;
        mvec2dp p1{0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
        mvec2dp p2{0.0, 2.0, 4.0, 6.0, 0.0, 2.0, 4.0, 6.0};
        mvec2dp p3{0.0, 3.0, 6.0, 9.0, 0.0, 3.0, 6.0, 9.0};
        mvec2dp p4 = -p1; // assignment using unary minus
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
    // MVec2dp<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2dp: geometric product tests - vec * vec")
    {
        fmt::println("MVec2dp: geometric product tests - vec * vec");

        // ab = dot(a,b) + wdg(a,b) = gr0(ab) + gr2(ab)
        //
        // dot(a,b) = 0.5*(ab + ba)   (symmetric part)
        // wdg(a,b) = 0.5*(ab - ba)   (antisymmetric part)

        vec2dp a{1.0, 2.0, 3.0};
        vec2dp b{0.5, 3.0, -2.0};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        mvec2dp mva{a};
        mvec2dp mvb{b};
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

    TEST_CASE("MVec2dp: geometric product tests - bivec * vec")
    {
        fmt::println("MVec2dp: geometric product tests - bivec * vec");

        // Ab = dot(A,b) + wdg(A,b) = gr1(Ab) + gr3(Ab)
        //
        // dot(A,b) = 0.5*(Ab - Aa)   (antisymmetric part)
        // wdg(A,b) = 0.5*(Ab + Aa)   (symmetric part)

        bivec2dp a{1.0, 2.0, 3.0};
        vec2dp b{0.5, 3.0, -2.0};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        mvec2dp mva{a};
        mvec2dp mvb{b};
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

        CHECK(dot_ab == gr1(mvab));
        CHECK(dot_ab == gr1(mvab_asym));
        CHECK(wdg_ab == gr3(mvab));
        CHECK(wdg_ab == gr3(mvab_sym));
    }

    TEST_CASE("MVec2dp: geometric product tests - vec * bivec")
    {
        fmt::println("MVec2dp: geometric product tests - vec * bivec");

        // a*B = dot(a,B) + wdg(a,B) = gr1(aB) + gr3(aB)
        //
        // dot(a,B) = 0.5*(aB - Ba)   (antisymmetric part)
        // wdg(a,B) = 0.5*(aB + Ba)   (symmetric part)

        vec2dp a{1.0, 2.0, 3.0};
        bivec2dp b{0.5, 3.0, -2.0};
        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        mvec2dp mva{a};
        mvec2dp mvb{b};
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

        CHECK(dot_ab == gr1(mvab));
        CHECK(dot_ab == gr1(mvab_asym));
        CHECK(wdg_ab == gr3(mvab));
        CHECK(wdg_ab == gr3(mvab_sym));
    }

    TEST_CASE("MVec2dp: geometric product tests - equivalence tests")
    {
        fmt::println("MVec2dp: geometric product tests - equivalence tests");

        vec2dp a{1.0, 2.0, 3.0};
        vec2dp b{0.5, 3.0, -4.0};
        mvec2dp mva{a};
        mvec2dp mvb{b};

        bivec2dp A{1.0, 2.0, 3.0};
        bivec2dp B{0.5, 3.0, -4.0};
        mvec2dp mvA{A};
        mvec2dp mvB{B};

        auto dot_ab = dot(a, b);
        auto wdg_ab = wdg(a, b);

        auto dot_Ab = dot(A, b);
        auto wdg_Ab = wdg(A, b);

        auto dot_aB = dot(a, B);
        auto wdg_aB = wdg(a, B);

        mvec2dp_e ab = a * b;
        mvec2dp abm = mva * mvb;
        mvec2dp abd{scalar{dot_ab}, wdg_ab};

        mvec2dp_u Ab = A * b;
        mvec2dp Abm = mvA * mvb;
        mvec2dp Abd{dot_Ab, wdg_Ab};

        mvec2dp_u aB = a * B;
        mvec2dp aBm = mva * mvB;
        mvec2dp aBd{dot_aB, wdg_aB};

        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   ab  = mvec2dp_e(a * b)            = {}", ab);
        // fmt::println("   abm = mva * mvb                  = {}", abm);
        // fmt::println("   abd = mvec2dp(dot(a,b), wdg(a,b)) = {}", abd);
        // fmt::println("");
        // fmt::println("   A                                = {}", A);
        // fmt::println("   mvA                              = {}", mvA);
        // fmt::println("   b                                = {}", b);
        // fmt::println("   mvb                              = {}", mvb);
        // fmt::println("   Ab  = mvec2dp_u(A * b)            = {}", Ab);
        // fmt::println("   Abm = mvA * mvb                  = {}", Abm);
        // fmt::println("   Abd = mvec2dp(dot(A,b), wdg(A,b)) = {}", Abd);
        // fmt::println("");
        // fmt::println("   a                                = {}", a);
        // fmt::println("   mva                              = {}", mva);
        // fmt::println("   B                                = {}", B);
        // fmt::println("   mvB                              = {}", mvB);
        // fmt::println("   aB  = mvec2dp_u(a * B)            = {}", aB);
        // fmt::println("   aBm = mva * mvB                  = {}", aBm);
        // fmt::println("   aBd = mvec2dp(dot(a,B), wdg(a,B)) = {}", aBd);
        // fmt::println("");

        CHECK(gr0(ab) == gr0(abm));
        CHECK(gr1(abm) == vec2dp{});
        CHECK(gr2(ab) == gr2(abm));
        CHECK(gr3(abm) == pscalar2dp{0.0});

        CHECK(gr0(ab) == gr0(abd));
        CHECK(gr1(abd) == vec2dp{});
        CHECK(gr2(ab) == gr2(abd));
        CHECK(gr3(abd) == pscalar2dp{0.0});

        CHECK(gr0(Abm) == 0);
        CHECK(gr1(Ab) == gr1(Abm));
        CHECK(gr2(Abm) == bivec2dp{});
        CHECK(gr3(Ab) == gr3(Abm));

        CHECK(gr0(Abd) == 0);
        CHECK(gr1(Ab) == gr1(Abd));
        CHECK(gr2(Abd) == bivec2dp{});
        CHECK(gr3(Ab) == gr3(Abd));

        CHECK(gr0(aBm) == 0);
        CHECK(gr1(aB) == gr1(aBm));
        CHECK(gr2(aBm) == bivec2dp{});
        CHECK(gr3(aB) == gr3(aBm));

        CHECK(gr0(aBd) == 0);
        CHECK(gr1(aB) == gr1(aBd));
        CHECK(gr2(aBd) == bivec2dp{});
        CHECK(gr3(aB) == gr3(aBd));
    }

    TEST_CASE("MVec2dp: assignment tests")
    {
        fmt::println("MVec2dp: assignment tests");

        vec2dp v1{1.0, 2.0, 3.0};
        vec2dp v2{0.5, 1.0, 1.5};
        vec2dp v3{0.5, 1.0, -4.5};
        bivec2dp b1{1.0, 2.0, 3.0};

        mvec2dp mv1{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
        mvec2dp mv2{0.0, 0.5, 1.0, 1.5, 11.5, 15.5, 6.0, 61.5};
        mvec2dp mv3{mv1};
        mvec2dp mv4 = mv2;

        mvec2dp mv5(scalar(5.0));
        mvec2dp mv6(pscalar2dp(6.0));
        mvec2dp mv7{v1};
        mvec2dp mv8{b1};
        mvec2dp mv9{scalar{dot(v1, v3)}, wdg(v1, v3)};

        mvec2dp mv10{v1, pscalar2dp(10.0)};
        // This must not compile! Implict conversion to vec2dp possible
        // possible solution: explicitly deleted constructor for MVec2dp
        // mvec2dp mv11{b1, pscalar2dp_t(10.0)};

        // this does not compile (which is fine, a base cannot convert to derived)
        // mvec2dp mv12{scalar(10.0), v1};

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

    TEST_CASE("MVec2dp: bivector product properties")
    {
        fmt::println("MVec2dp: bivector product properties");

        bivec2dp b1{1.0, 2.0, 3.0};
        mvec2dp mb1{b1};
        bivec2dp b2{-3.0, 1.0, 2.0};
        mvec2dp mb2{b2};

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
    // MVec2dp_E<T> and MVec2dp_U<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("MVec2dp_E/_U: modelling even and uneven parts of 3d algebra - basics")
    {
        fmt::println(
            "MVec2dp_E/_U: modelling even and uneven parts of 3d algebra - basics");

        // defining a complex number in all three forms as multivector
        auto u = normalize(vec2dp{1.0, 0.0, 0.0});
        auto v =
            normalize(vec2dp(std::cos(pi / 12), std::sin(pi / 12), 0.0)); // unit vec +15%
        auto angle_uv = angle(u, v);
        auto B = wdg(u, v); // normalized bivector describing the plane spanned by u and v

        auto my_exp = exp(-B, angle_uv);
        auto my_rot = rotor(B, 2.0 * angle_uv);

        // definition of rotor used here: B = u^v
        // => B determines the meaning of the positive sign of the rotation
        //
        auto R_m =
            mvec2dp(exp(-B, angle_uv)); // Rotor formed by u and v (normalized bivector)
        auto Rr_m = mvec2dp(rev(R_m));  // and its reverse

        auto c = vec2dp{1.0, 1.0, 1.0};
        auto c_m = mvec2dp{c};

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
        // fmt::println("   angle(u,v)            = {: .3}°", rad_to_deg(angle_uv));
        // fmt::println("   sin(angle(u,v))       = {: .3}", std::sin(angle_uv));
        // fmt::println("");
        // fmt::println("   c                     = {: .3}", c);
        // fmt::println("");
        // fmt::println("Implemented as full multivector operation:");
        // fmt::println("   R_m  = mvec2dp(exp(-B,angle_uv))  = {: .3}", R_m);
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
        // fmt::println("   angle(c, c_rot_l) = {: .3}°", rad_to_deg(angle_c_c_rot));
        // fmt::println("   angle(projected)  = {: .3}°", rad_to_deg(angle_proj));
        // fmt::println("");
        // fmt::println("direct calclulation:");
        // fmt::println("   c_rot = rotate(c,R)          = {: .3}", rotate(c, R));

        // auto normalized_vec = normalize(vec2dp{1.0, 1.0, 0.0});
        // fmt::println("   normalized_vec = {}", normalized_vec);
        fmt::println("   normalized_vec = {}", normalize(vec2dp{1.0, 1.0, 0.0}));

        CHECK(nrm(rotate(c, R)) == nrm(c));
        CHECK(gr1(c_rot_m) == rotate(c, R));
        // n I_2dp approach:
        CHECK(rotate(vec2dp{1.0, 0.0, 0.0}, rotor(e3_2dp * I_2dp, -pi / 4)) ==
              normalize(vec2dp{1.0, 1.0, 0.0}));
        // using a bivector directly:
        CHECK(rotate(vec2dp{1.0, 0.0, 0.0}, rotor(e12_2dp, pi / 4)) ==
              normalize(vec2dp{1.0, 1.0, 0.0}));

        // direct rotation of a bivector
        CHECK(rotate(bivec2dp{0.0, 0.0, 1.0}, rotor(e23_2dp, pi / 2)) == -e31_2dp);

        // example see Macdonald "Linear and Geometric Algebra", Exercise 7.12, p. 127
        auto Bv = wdg(e2_2dp,
                      e1_2dp + std::sqrt(3.0) * e3_2dp); // bivector describing the plane
        CHECK(std::abs(magn(Bv) - 2.0) < eps);
        CHECK(rotate(Bv, rotor(e31_2dp, pi / 3)) == -2.0 * e12_2dp);

        // just to silence unused variable warnings
        CHECK(my_exp == exp(-B, angle_uv));
        CHECK(my_rot == rotor(B, 2.0 * angle_uv));
        CHECK(c_rot_r == gr1(c_rot_u_r));
        CHECK(angle_c_c_rot == angle(c, c_rot_l));
        CHECK(angle_proj == angle(c_proj, c_rot_proj));
    }

    //     TEST_CASE("MVec2dp: dualization")
    //     {
    //         fmt::println("MVec2dp: dualization");

    //         vec2dp v{1.0, 2.0, 3.0};                                   // 3d vector
    //         bivec2dp B{10.0, 20.0, 30.0};                              // 3d bivector
    //         mvec2dp vm{100.0, 1.0, 2.0, 3.0, 10.0, 20.0, 30.0, 300.0}; // full 3d
    //         multivector mvec2dp vm2{-20.0, 3.0,  7.0,    -4.0,
    //                     -1.0,  40.0, -330.0, -70.0}; // full 3d multivector

    //         // full 3d multivector - even content
    //         mvec2dp vm_even{100.0, 0.0, 0.0, 0.0, 10.0, 20.0, 30.0, 0.0};
    //         // even grade 3d multivector
    //         mvec2dp_e vm_E{100.0, 10.0, 20.0, 30.0};

    //         // full 3d multivector - uneven content
    //         mvec2dp vm_uneven{0.0, 1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 300.0};
    //         // uneven grade 3d multivector
    //         mvec2dp_u vm_U{1.0, 2.0, 3.0, 300.0};

    // #if defined(_HD_GA_HESTENES_DORAN_LASENBY_DUALITY)
    //         ////////////////////////////////////////////////////////////////////////////////
    //         // duality as defined by Hestenes or by Doran, Lasenby in "GA for
    //         physicists":
    //         // (same subspace as for Macdonld's definition below, but other resulting
    //         signs)
    //         // (=> will have influence on all formulae concerning duality)
    //         ////////////////////////////////////////////////////////////////////////////////
    //         //
    //         // dual(A) = I*A
    //         //

    //         auto vm_dual_manual = Im_2dp * vm;
    //         auto vm_dual = dual3d(vm);

    //         auto vm_dual_even_manual = Im_2dp * vm_even;
    //         auto vm_dual_even = dual3d(vm_even);

    //         auto vm_dual_uneven_manual = Im_2dp * vm_uneven;
    //         auto vm_dual_uneven = dual3d(vm_uneven);

    //         // result is uneven, naming chosen for consistency
    //         auto vm_dual_manual_E = I_2dp * vm_E;
    //         auto vm_dual_E = dual3d(vm_E);

    //         // result is even, naming chosen for consistency
    //         auto vm_dual_manual_U = Im_2dp_U * vm_U;
    //         auto vm_dual_U = dual3d(vm_U);

    //         auto v_dual_manual = I_2dp * v;
    //         auto v_dual = dual3d(v);

    //         auto B_dual_manual = I_2dp * B;
    //         auto B_dual = dual3d(B);

    //         // fmt::println("   I_2dp    = {}", I_2dp);
    //         // fmt::println("   Im_2dp   = {}", Im_2dp);
    //         // fmt::println("   Im_2dp_U = {}", Im_2dp_U);
    //         // fmt::println("");
    //         // fmt::println("   v             = {}", v);
    //         // fmt::println("   B             = {}", B);
    //         // fmt::println("");
    //         // fmt::println("   vm            = {}", vm);
    //         // fmt::println("   Im_2dp*vm      = {}", vm_dual_manual);
    //         // fmt::println("   dual3d(vm)      = {}", vm_dual);
    //         // fmt::println("");
    //         // fmt::println("   vm_even       = {}", vm_even);
    //         // fmt::println("   Im_2dp*vm_even = {}", vm_dual_even_manual);
    //         // fmt::println("   dual3d(vm_even) = {}", vm_dual_even);
    //         // fmt::println("");
    //         // fmt::println("   vm_E          = {}", vm_E);
    //         // fmt::println("   Im_2dp_E*vm_E  = {}", vm_dual_manual_E);
    //         // fmt::println("   dual3d(vm_E)    = {}", vm_dual_E);
    //         // fmt::println("");
    //         // fmt::println("   vm_uneven       = {}", vm_uneven);
    //         // fmt::println("   Im_2dp*vm_uneven = {}", vm_dual_uneven_manual);
    //         // fmt::println("   dual3d(vm_uneven) = {}", vm_dual_uneven);
    //         // fmt::println("");
    //         // fmt::println("   vm_U          = {}", vm_U);
    //         // fmt::println("   Im_2dp_U*vm_U  = {}", vm_dual_manual_U);
    //         // fmt::println("   dual3d(vm_U)    = {}", vm_dual_U);
    //         // fmt::println("");
    //         // fmt::println("   v               = {}", v);
    //         // fmt::println("   I_2dp * v        = {} - bivec ", v_dual_manual);
    //         // fmt::println("   dual3d(v)         = {} - bivec ", v_dual);
    //         // fmt::println("");
    //         // fmt::println("   B               = {}", B);
    //         // fmt::println("   I_2dp * B        = {} - vec", B_dual_manual);
    //         // fmt::println("   dual3d(B)         = {} - vec", B_dual);

    //         CHECK(vm_dual == vm_dual_manual);
    //         CHECK(vm_dual_even == vm_dual_even_manual);
    //         CHECK(vm_dual_uneven == vm_dual_uneven_manual);
    //         CHECK(vm_dual_E == vm_dual_manual_E);
    //         CHECK(vm_dual_U == vm_dual_manual_U);
    //         CHECK(dual3d(v) == bivec2dp{1.0, 2.0, 3.0});
    //         CHECK(dual3d(B) == -vec2dp{10.0, 20.0, 30.0});
    //         CHECK(dual3d(scalar(5)) == pscalar2dp(5));
    //         CHECK(dual3d(pscalar2dp(5)) == scalar(-5));

    //         // dual properties
    //         CHECK(dual3d(3.0 * vm) == 3.0 * dual3d(vm));
    //         CHECK(dual3d(vm + vm2) == dual3d(vm) + dual3d(vm2));
    //         CHECK(dual3d(dual3d(vm)) == -vm);
    //         CHECK(dual3d(I_2dp) == -1.0);

    //         // inner and outer products are in G^n are dual (A. Macdonald, p. 111):
    //         //
    //         // dual(dot(A,B)) = wdg(A,dual(B)), dual(wdg(A,B)) = dot(A,dual(B))
    //         //
    //         // or alternatively:
    //         //
    //         // (a.b)^* = a^(b^*), (a^b)^* = a.(b^*)

    //         auto a = vec2dp{2.0, -3.0, 2.0};
    //         auto b = vec2dp{-5.0, 1.0, -7.0};

    //         // fmt::println("   a                 = {}", a);
    //         // fmt::println("   b                 = {}", b);
    //         // fmt::println("   dual3d(b)         = {}", dual3d(b));
    //         // fmt::println("   dot(a, b)         = {}", dot(a, b));
    //         // fmt::println("   wdg(a, b)         = {}", wdg(a, b));
    //         // fmt::println("");
    //         // fmt::println("   dual3d(dot(a, b)) = {}", dual3d(dot(a, b)));
    //         // fmt::println("   wdg(a, dual3d(b)) = {}", wdg(a, dual3d(b)));
    //         // fmt::println("");
    //         // fmt::println("   dual3d(wdg(a, b)) = {}", dual3d(wdg(a, b)));
    //         // fmt::println("   dot(a, dual3d(b)) = {}", dot(a, dual3d(b)));

    //         CHECK(dot(a, b * I_2dp) == wdg(a, b) * I_2dp);

    //         // just to silence unused variable warning
    //         CHECK(v_dual_manual == I_2dp * v);
    //         CHECK(v_dual == dual3d(v));
    //         CHECK(B_dual_manual == I_2dp * B);
    //         CHECK(B_dual == dual3d(B));

    // #else
    //         ////////////////////////////////////////////////////////////////////////////////
    //         // duality (as defined in Macdonald, "Linear and geometric algebra"):
    //         ////////////////////////////////////////////////////////////////////////////////
    //         //
    //         // dual(A) = A/I = A*I^(-1) = A*rev(I)
    //         //

    //         auto vm_dual_manual = vm * rev(Im_2dp);
    //         auto vm_dual = dual3d(vm);

    //         auto vm_dual_even_manual = vm_even * rev(Im_2dp);
    //         auto vm_dual_even = dual3d(vm_even);

    //         auto vm_dual_uneven_manual = vm_uneven * rev(Im_2dp);
    //         auto vm_dual_uneven = dual3d(vm_uneven);

    //         // result is uneven, naming chosen for consistency
    //         auto vm_dual_manual_E = vm_E * rev(I_2dp);
    //         auto vm_dual_E = dual3d(vm_E);

    //         // result is even, naming chosen for consistency
    //         auto vm_dual_manual_U = vm_U * rev(Im_2dp_U);
    //         auto vm_dual_U = dual3d(vm_U);

    //         auto v_dual_manual = v * rev(I_2dp);
    //         auto v_dual = dual3d(v);

    //         auto B_dual_manual = B * rev(I_2dp);
    //         auto B_dual = dual3d(B);

    //         // fmt::println("   I_2dp                 = {}", I_2dp);
    //         // fmt::println("   Im_2dp                = {}", Im_2dp);
    //         // fmt::println("   Im_2dp_U              = {}", Im_2dp_U);
    //         // fmt::println("");
    //         // fmt::println("   v                    = {}", v);
    //         // fmt::println("   B                    = {}", B);
    //         // fmt::println("");
    //         // fmt::println("   vm                   = {}", vm);
    //         // fmt::println("   vm*rev(Im_2dp)        = {}", vm_dual_manual);
    //         // fmt::println("   dual3d(vm)           = {}", vm_dual);
    //         // fmt::println("");
    //         // fmt::println("   vm_even              = {}", vm_even);
    //         // fmt::println("   vm_even*rev(Im_2dp)   = {}", vm_dual_even_manual);
    //         // fmt::println("   dual3d(vm_even)      = {}", vm_dual_even);
    //         // fmt::println("");
    //         // fmt::println("   vm_E                 = {}", vm_E);
    //         // fmt::println("   vm_E * rev(I_2dp)     = {}", vm_dual_manual_E);
    //         // fmt::println("   dual3d(vm_E)         = {}", vm_dual_E);
    //         // fmt::println("");
    //         // fmt::println("   vm_uneven            = {}", vm_uneven);
    //         // fmt::println("   vm_uneven*rev(Im_2dp) = {}", vm_dual_uneven_manual);
    //         // fmt::println("   dual3d(vm_uneven)    = {}", vm_dual_uneven);
    //         // fmt::println("");
    //         // fmt::println("   vm_U                 = {}", vm_U);
    //         // fmt::println("   vm_U*rev(Im_2dp_U)    = {}", vm_dual_manual_U);
    //         // fmt::println("   dual3d(vm_U)         = {}", vm_dual_U);
    //         // fmt::println("");
    //         // fmt::println("   v                    = {}", v);
    //         // fmt::println("   v*rev(I_2dp)          = {} - bivec ", v_dual_manual);
    //         // fmt::println("   dual3d(v)            = {} - bivec ", v_dual);
    //         // fmt::println("");
    //         // fmt::println("   B                    = {}", B);
    //         // fmt::println("   B*rev(I_2dp)          = {} - vec", B_dual_manual);
    //         // fmt::println("   dual3d(B)            = {} - vec", B_dual);

    //         CHECK(vm_dual == vm_dual_manual);
    //         CHECK(vm_dual_even == vm_dual_even_manual);
    //         CHECK(vm_dual_uneven == vm_dual_uneven_manual);
    //         CHECK(vm_dual_E == vm_dual_manual_E);
    //         CHECK(vm_dual_U == vm_dual_manual_U);
    //         CHECK(dual3d(v) == -bivec2dp{1.0, 2.0, 3.0});
    //         CHECK(dual3d(B) == vec2dp{10.0, 20.0, 30.0});
    //         CHECK(dual3d(scalar(5)) == pscalar2dp(-5));
    //         CHECK(dual3d(pscalar2dp(6)) == scalar(6));

    //         // dual properties (A. Macdonald, p. 110):
    //         //
    //         // a) dual(aA) = a dual(A)
    //         // b) dual(A + B) = dual(A) + dual(B)
    //         // c) dual(dual(A)) = (-1)^(n*(n-1)/2) A   (with n as dimension of the
    //         (sub)space)
    //         // d) |dual(B)| = |B|
    //         // e) if B is a j-blade then dual(B) is an (n-j)-blade
    //         // f) if A is a j-vector then dual(A) is an (n-j)-vector
    //         //    (remember: a j-vector is a sum of j-blades, which are outer products)

    //         CHECK(dual3d(3.0 * vm) == 3.0 * dual3d(vm));
    //         CHECK(dual3d(vm + vm2) == dual3d(vm) + dual3d(vm2));
    //         CHECK(dual3d(dual3d(vm)) == -vm);
    //         CHECK(dual3d(I_2dp) == 1);
    //         CHECK(dual3d(v) == -bivec2dp(1.0, 2.0, 3.0));

    //         // inner and outer products are in G^n are dual (A. Macdonald, p. 111):
    //         //
    //         // dual(dot(A,B)) = wdg(A,dual(B)), dual(wdg(A,B)) = dot(A,dual(B))

    //         auto a = vec2dp{2.0, -3.0, 2.0};
    //         auto b = vec2dp{-5.0, 1.0, -7.0};

    //         // fmt::println("   a                 = {}", a);
    //         // fmt::println("   b                 = {}", b);
    //         // fmt::println("   dual3d(b)         = {}", dual3d(b));
    //         // fmt::println("   dot(a, b)         = {}", dot(a, b));
    //         // fmt::println("   wdg(a, b)         = {}", wdg(a, b));
    //         // fmt::println("");
    //         // fmt::println("   dual3d(dot(a, b)) = {}", dual3d(dot(a, b)));
    //         // fmt::println("   wdg(a, dual3d(b)) = {}", wdg(a, dual3d(b)));
    //         // fmt::println("");
    //         // fmt::println("   dual3d(wdg(a, b)) = {}", dual3d(wdg(a, b)));
    //         // fmt::println("   dot(a, dual3d(b)) = {}", dot(a, dual3d(b)));

    //         CHECK(dual3d(scalar{dot(a, b)}) == wdg(a, dual3d(b)));
    //         CHECK(dual3d(dot(a, b)) == wdg(a, dual3d(b)));
    //         CHECK(dual3d(wdg(a, b)) == dot(a, dual3d(b)));

    //         // just to silence unused variable warning
    //         CHECK(v_dual_manual == v * rev(I_2dp));
    //         CHECK(v_dual == dual3d(v));
    //         CHECK(B_dual_manual == B * rev(I_2dp));
    //         CHECK(B_dual == dual3d(B));
    // #endif
    //     }

} // TEST_SUITE("Projective Geometric Algebra (PGA)")