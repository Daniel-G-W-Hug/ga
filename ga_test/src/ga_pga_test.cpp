// author: Daniel Hug, 2024 & 2025

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

    //////////////////////////////////////////////////////////////////////////////////////
    // projective geometric algebra 2d: pga2dp (embedded in a 3d representational space)
    //////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("algebra<2,0,1> - pga2dp")
    {
        fmt::println("algebra<2, 0, 1> - pga2dp:");
        // 2d projective geometric algebra (2dp) models two-dimensional Euclidean algebra
        // embedding it in a three-dimensional representational space
        const algebra<2, 0, 1> alg;
        CHECK(alg.p() == 2);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 1);
        CHECK(alg.dim_space() == 3);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 8);            // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 4); // == dim_space + 1
        fmt::println("   pga2dp: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   pga2dp: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("G<2,0,1> - pga2dp defining basic types and ctor checks")
    {
        fmt::println("G<2,0,1>: pga2dp defining basic types and ctor checks");

        auto s = scalar2dp{5.0};
        auto v = vec2dp{1.0, 2.0, 1.0};
        auto B = bivec2dp{-1.0, 2.0, 1.0};
        auto ps = pscalar2dp{-5.0};

        auto mv0 = mvec2dp{s};
        auto mv1 = mvec2dp{v};
        auto mv2 = mvec2dp{B};
        auto mv3 = mvec2dp{ps};

        auto mv5a = mvec2dp{s, B};
        auto mv5b = mvec2dp_e{s};
        auto mv5c = mvec2dp_e{B};
        auto mv5d = mvec2dp_e{s, B};
        auto mv5e = mvec2dp{mv5d};

        auto mv6a = mvec2dp{v, ps};
        auto mv6b = mvec2dp_u{v};
        auto mv6c = mvec2dp_u{ps};
        auto mv6d = mvec2dp_u{v, ps};
        auto mv6e = mvec2dp{mv6d};

        auto mv7 = mvec2dp{s, v, B, ps};

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
        CHECK(gr1(mv5a) == vec2dp{});
        CHECK(gr2(mv5a) == B);
        CHECK(gr3(mv5a) == pscalar2dp{});

        CHECK(gr0(mv5b) == s);
        CHECK(gr2(mv5b) == bivec2dp{});

        CHECK(gr0(mv5c) == scalar2dp{});
        CHECK(gr2(mv5c) == B);

        CHECK(gr0(mv5d) == s);
        CHECK(gr2(mv5d) == B);

        CHECK(gr0(mv5e) == s);
        CHECK(gr1(mv5e) == vec2dp{});
        CHECK(gr2(mv5e) == B);
        CHECK(gr3(mv5e) == pscalar2dp{});

        CHECK(gr0(mv6a) == scalar2dp{});
        CHECK(gr1(mv6a) == v);
        CHECK(gr2(mv6a) == bivec2dp{});
        CHECK(gr3(mv6a) == ps);

        CHECK(gr1(mv6b) == v);
        CHECK(gr3(mv6b) == pscalar2dp{});

        CHECK(gr1(mv6c) == vec2dp{});
        CHECK(gr3(mv6c) == ps);

        CHECK(gr1(mv6d) == v);
        CHECK(gr3(mv6d) == ps);

        CHECK(gr0(mv6e) == scalar2dp{});
        CHECK(gr1(mv6e) == v);
        CHECK(gr2(mv6e) == bivec2dp{});
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

    ////////////////////////////////////////////////////////////////////////////////
    // Vec2dp<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec2dp: default init")
    {
        fmt::println("Vec2dp: default init");
        vec2dp v;
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.x) < eps);
        CHECK(std::abs(v.y) < eps);
        CHECK(std::abs(v.z) < eps);
    }

    TEST_CASE("Vec2dp: with curly braced intializer")
    {
        fmt::println("Vec2dp: with curly braced intializer");
        vec2dp v{0.0, 0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.x) < eps);
        CHECK(std::abs(v.y) < eps);
        CHECK(std::abs(v.z) < eps);
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
        fmt::println("       fmt: vp1 = {:e}", fmt::join(vp1, ", "));
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

        // fmt::println("v1 = {: .4f}, nrm(v1) = {: .4f}", v1, nrm(v1));
        // fmt::println("v2 = normalize(v1) = {: .4f}, nrm(v2) = {: .4f}", v2, nrm(v2));
        // fmt::println("v3 = {: .4f}, nrm(v1) = {: .4f}", v3, nrm(v3));
        // fmt::println(
        //     "v4 = inv(v3) = {: .4f}, nrm(v3) = {: .4f}, nrm(v3)*nrm(v4) = {: .4f}",
        //     v4, nrm(v4), nrm(v3) * nrm(v4));
        // fmt::println("");
        // fmt::println("v3m = {}, nrm(v3m) = {}", v3m, nrm(v3m));
        // fmt::println("v4m = {}, nrm(v4m) = {}", v4m, nrm(v4m));
        // fmt::println("0.5*(v3m*v4m + v4m*v3m) = {}", 0.5 * (v3m * v4m + v4m * v3m));
        // fmt::println("0.5*(v3m*v4m - v4m*v3m) = {}", 0.5 * (v3m * v4m - v4m * v3m));

        CHECK(std::abs(nrm_sq(v1) - 9.0) < eps);
        CHECK(std::abs(bulk_nrm_sq(v1) - 5.0) < eps);
        CHECK(std::abs(nrm_sq(v2) - 1.0) < eps);
        CHECK(std::abs(bulk_nrm_sq(v3) - 40.0) < eps);
        CHECK(std::abs(dot(v4, v3) - 1.0) < eps);
    }

    TEST_CASE("Vec2dp: operations - bulk_dual, weight_dual")
    {
        fmt::println("Vec2dp: operations - bulk_dual, weight_dual");

        auto s = scalar2dp{5.0};
        auto v = vec2dp{1.0, 2.0, 1.0};
        auto B = bivec2dp{-1.0, 2.0, 1.0};
        auto ps = pscalar2dp{-5.0};

        auto v2 = vec2dp{3.0, -2.0, 1.0};
        auto B2 = bivec2dp{-3.0, 2.0, 5.0};

        CHECK(bulk_dual(s) == pscalar2dp{s});
        CHECK(weight_dual(s) == pscalar2dp{0.0});

        CHECK(bulk_dual(v) == bivec2dp{1.0, 2.0, 0.0});
        CHECK(weight_dual(v) == bivec2dp{0.0, 0.0, 1.0});

        CHECK(bulk_dual(B) == vec2dp{0.0, 0.0, 1.0});
        CHECK(weight_dual(B) == vec2dp{-1.0, 2.0, 0.0});

        CHECK(bulk_dual(ps) == scalar2dp{0.0});
        CHECK(weight_dual(ps) == scalar2dp{ps});

        // duality of wdg and rwdg based on complements
        CHECK(bulk_dual(wdg(v, v2)) == rwdg(bulk_dual(v), bulk_dual(v2)));
        CHECK(bulk_dual(wdg(v, B)) == rwdg(bulk_dual(v), bulk_dual(B)));
        CHECK(bulk_dual(wdg(B, v)) == rwdg(bulk_dual(B), bulk_dual(v)));

        // contractions
        CHECK(bulk_contraction(v, v2) == rwdg(v, bulk_dual(v2)));
        CHECK(bulk_contraction(B, B2) == rwdg(B, bulk_dual(B2)));
        CHECK(bulk_contraction(B, v) == rwdg(B, bulk_dual(v))); // contracts v onto B

        CHECK(weight_contraction(v, v2) == rwdg(v, weight_dual(v2)));
        CHECK(weight_contraction(B, B2) == rwdg(B, weight_dual(B2)));
        CHECK(weight_contraction(B, v) == rwdg(B, weight_dual(v)));

        // expansions
        CHECK(weight_expansion(v, v2) == wdg(v, weight_dual(v2)));
        CHECK(weight_expansion(B, B2) == wdg(B, weight_dual(B2)));
        CHECK(weight_expansion(v, B) == wdg(v, weight_dual(B)));

        CHECK(bulk_expansion(v, v2) == wdg(v, bulk_dual(v2)));
        CHECK(bulk_expansion(B, B2) == wdg(B, bulk_dual(B2)));
        CHECK(bulk_expansion(v, B) == wdg(v, bulk_dual(B)));
    }

    TEST_CASE("Vec2dp: operations - angle I")
    {
        fmt::println("Vec2dp: operations - angle I");

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
            v1.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .g},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e1_2dp, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2dp<double>(std::cos(phi + pi / 2), std::sin(phi + pi / 2), 0.0);
            v2.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .g},"
            //              " angle={: .4f}",
            //              i, phi, rad_to_deg(phi), c, angle(e2_2dp, c));
        }
        // fmt::println("");

        for (int i = 0; i <= 12; ++i) {
            double phi = i * pi / 12;
            auto c = Vec2dp<double>(std::cos(phi + pi / 4), std::sin(phi + pi / 4), 0.0);
            v3.emplace_back(std::make_tuple(phi, c));
            // fmt::println("   i={: 3}: phi={: .4f}, phi={: 4.0f}°, c={: .g},"
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

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar2dp{sd};
        auto t = scalar2dp{st};

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

        CHECK(wdg(v1, v1) == bivec2dp{});                     // wdg=0 for collin. vectors
        CHECK(wdg(v1, v2) == -wdg(v2, v1));                   // anticommutative for vect.
        CHECK(wdg(wdg(v1, v2), v3) == wdg(v1, wdg(v2, v3)));  // wdg is associative
        CHECK(wdg(v1, v2 + v3) == wdg(v1, v2) + wdg(v1, v3)); // wdg distributes over add.
        CHECK(wdg(v1 + v2, v3) == wdg(v1, v3) + wdg(v2, v3)); // wdg distributes over add.
        CHECK(wdg(sd * v1, v2) == wdg(v1, sd * v2)); // scalars can be factored out of wdg
        CHECK(wdg(sd * v1, v2) == sd * wdg(v1, v2)); // scalars can be factored out of wdg
        CHECK(wdg(s, t) == wdg(t, s));   // wdg between scalars equivalent to scalar mult.
        CHECK(wdg(s, v1) == wdg(v1, s)); // wdg between scalar and vector
        CHECK(wdg(s, v1) == sd * v1);    // wdg between scalar and vector

        CHECK(std::abs(bulk_nrm(wdg(v1, v1)) - sin(angle(v1, v1))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v2)) - sin(angle(v1, v2))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v3)) - sin(angle(v1, v3))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v4)) - sin(angle(v1, v4))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v5)) - sin(angle(v1, v5))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v6)) - sin(angle(v1, v6))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v7)) - sin(angle(v1, v7))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v8)) - sin(angle(v1, v8))) < eps);
    }

    TEST_CASE("Vec2dp: operations - project / reject / reflect (vector - vector)")
    {
        fmt::println("Vec2dp: operations - project / reject / reflect (vector - vector)");

        vec2dp v1{5.0, 1.0, 1.0};
        vec2dp v2{2.0, 2.0, 1.0};

        vec2dp v3{project_onto(v1, v2)};
        vec2dp v4{reject_from(v1, v2)};
        vec2dp v5{v3 + v4};

        // fmt::println("sqrt(26)={}", std::sqrt(26.0));
        // fmt::println("v1  = {: .4f}, bulk_nrm(v1)  = {: .4f}, "
        //              "nrm(v1) = {: .4f}",
        //              v1, bulk_nrm(v1), nrm(v1));
        // fmt::println("v2  = {: .4f}, bulk_nrm(v2)  = {: .4f}, "
        //              "nrm(v2) = {: .4f}",
        //              v2, bulk_nrm(v2), nrm(v2));
        // fmt::println("v2n = {: .4f}, bulk_nrm(v2n) = {: .4f}, "
        //              "nrm(v2n) = {: .4f}",
        //              v2n, bulk_nrm(v2n), nrm(v2n));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .4f}, bulk_nrm(v3) = {: .4f}, "
        //              "nrm(v3) = {: .4f}",
        //              v3, bulk_nrm(v3), nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .4f}, bulk_nrm(v4) = {: .4f}, "
        //              "nrm(v4) = {: .4f}",
        //              v4, bulk_nrm(v4), nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .4f}, bulk_nrm(v5) = {: .4f}, "
        //              "nrm(v5) = {: .4f}",
        //              v5, bulk_nrm(v5), nrm(v5));
        // fmt::println("");
        // // this helps to understand, why the dot-product is sufficient
        // auto w = wdg(v1, v2);
        // auto i = inv(v2);
        // fmt::println("wdg(v1,v2)         = {: .4f}", w);
        // fmt::println("inv(v2)            = {: .4f}", i);
        // fmt::println("wdg(v1,v2)*inv(v2) = {: .4f}", w * i);
        // fmt::println("");

        // vec2dp v{4.0, 1.0, 0.0};
        // vec2dp b{e1_2dp};
        // auto B = bivec2dp{e23_2dp};

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

        // // just to suppress unused variable warnings
        // CHECK(b == e1_2dp);

        // CHECK(reflect_on_hyp(v, e2_2dp) == reflect_on(v, B));

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
        // fmt::println(
        //     "Vec2dp: operations - project / reject / reflect (vector - bivector)");

        // vec2dp v1{5.0, 3.0, 1.0};
        // bivec2dp v2 = wdg(vec2dp{0.0, 0.0, 2.0}, vec2dp{2.0, 0.0, 0.0});

        // vec2dp v3{project_onto(v1, v2)};
        // vec2dp v4{reject_from(v1, v2)};
        // vec2dp v5{v3 + v4};

        // fmt::println("v1  = {: .4f}, bulk_nrm(v1) = {: .4f}", v1, bulk_nrm(v1));
        // fmt::println("v2  = {: .4f}, bulk_nrm(v2) = {: .4f}", v2, bulk_nrm(v2));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .4f}, bulk_nrm(v3) = {: .4f}", v3,
        //              bulk_nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .4f}, bulk_nrm(v4) = {: .4f}", v4,
        //              bulk_nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .4f}, bulk_nrm(v5) = {: .4f}", v5,
        //              bulk_nrm(v5));

        // CHECK(v3 + v4 == v5);
        // CHECK(v5 == v1);

        // from Macdonald, "Linear and Geometric Algebra", Exercise 7.14, p. 129
        // CHECK(reflect_on(wdg(e1_2dp + e3_2dp, e2_2dp), e12_2dp) ==
        //       wdg(e1_2dp - e3_2dp, e2_2dp));
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
        CHECK(std::abs(v.c0) < eps);
        CHECK(std::abs(v.c1) < eps);
        CHECK(std::abs(v.c2) < eps);
        CHECK(std::abs(v.c3) < eps);
        CHECK(std::abs(v.c4) < eps);
        CHECK(std::abs(v.c5) < eps);
        CHECK(std::abs(v.c6) < eps);
        CHECK(std::abs(v.c7) < eps);
    }
    TEST_CASE("MVec2dp: with curly braced intializer")
    {
        fmt::println("MVec2dp: with curly braced intializer");
        // default initialization
        mvec2dp v{0.0, 1.0, 2.0, 3.0, 23.0, 31.0, 12.0, 123.0};
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

    TEST_CASE("MVec2dp: defining basic types and ctor checks")
    {
        fmt::println("MVec2dp: defining basic types and ctor checks");

        auto mv1 = mvec2dp{scalar2dp(5.0)};
        auto mv2 = mvec2dp{vec2dp{1.0, 2.0, 1.0}};
        auto mv3 = mvec2dp{bivec2dp{-1.0, 2.0, 1.0}};
        auto mv4 = mvec2dp{pscalar2dp(-5.0)};
        auto mv5a = mvec2dp_e{scalar2dp(5.0), bivec2dp{-1.0, 2.0, 1.0}};
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

        CHECK(gr0(mv1) == scalar2dp(5.0));
        CHECK(gr1(mv2) == vec2dp{1.0, 2.0, 1.0});
        CHECK(gr2(mv3) == bivec2dp{-1.0, 2.0, 1.0});
        CHECK(gr3(mv4) == pscalar2dp{-5.0});
        CHECK(gr0(mv5a) == scalar2dp(5.0));
        CHECK(gr2(mv5a) == bivec2dp{-1.0, 2.0, 1.0});
        CHECK(gr0(mv5) == scalar2dp(5.0));
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
        fmt::println("    fmt: vp1 = {:e}", fmt::join(vp1, ", "));
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

    TEST_CASE("MVec2dp: geometric product - combinatorial tests")
    {
        fmt::println("MVec2dp: geometric product - combinatorial tests");

        scalar2dp s1{3.0};
        scalar2dp s2{-1.5};
        vec2dp v1{1.0, 2.0, 3.0};
        vec2dp v2{0.5, 3.0, -2.0};
        bivec2dp b1{-1.0, -3.0, -5.0};
        bivec2dp b2{2.0, 0.5, -4.0};
        pscalar2dp ps1{-5.0};
        pscalar2dp ps2{2.0};

        // checks all ctor combinations and equivalence of simplified products with
        // geometric product of a fully populated multivector
        //
        // assumption: the geometric product of a fully populated multivector serves as
        //             reference and is correctly implemented
        //             results are relative to full geometric product

        // MVec2dp_E * MVec2dp
        CHECK(mvec2dp{s1, b1} * mvec2dp{s1, v1, b1, ps1} ==
              mvec2dp{mvec2dp_e{s1, b1} * mvec2dp{s1, v1, b1, ps1}});

        // MVec2dp * MVec2dp_E
        CHECK(mvec2dp{s1, v1, b1, ps1} * mvec2dp{s1, b1} ==
              mvec2dp{mvec2dp{s1, v1, b1, ps1} * mvec2dp_e{s1, b1}});

        // MVec2dp_E * MVec2dp_E
        CHECK(mvec2dp{s1, b1} * mvec2dp{s2, b2} ==
              mvec2dp{mvec2dp_e{s1, b1} * mvec2dp_e{s2, b2}});

        // MVec2dp_U * MVec2dp_U
        CHECK(mvec2dp{v1, ps1} * mvec2dp{v2, ps2} ==
              mvec2dp{mvec2dp_u{v1, ps1} * mvec2dp_u{v2, ps2}});

        // MVec2dp_E * MVec2dp_U
        CHECK(mvec2dp{s1, b1} * mvec2dp{v2, ps2} ==
              mvec2dp{mvec2dp_e{s1, b1} * mvec2dp_u{v2, ps2}});

        // MVec2dp_U * MVec2dp_E
        CHECK(mvec2dp{v1, ps1} * mvec2dp{s2, b2} ==
              mvec2dp{mvec2dp_u{v1, ps1} * mvec2dp_e{s2, b2}});

        // MVec2dp_E * BiVec2dp
        CHECK(mvec2dp{s1, b1} * mvec2dp{b2} == mvec2dp{mvec2dp_e{s1, b1} * mvec2dp{b2}});

        // BiVec2dp * MVec2dp_E
        CHECK(mvec2dp{b1} * mvec2dp{s2, b2} == mvec2dp{mvec2dp{b1} * mvec2dp_e{s2, b2}});

        // MVec2dp_E * Vec2dp
        CHECK(mvec2dp{s1, b1} * mvec2dp{v2} == mvec2dp{mvec2dp_e{s1, b1} * mvec2dp{v2}});

        // Vec2dp * MVec2dp_E
        CHECK(mvec2dp{v1} * mvec2dp{s2, b2} == mvec2dp{mvec2dp{v1} * mvec2dp_e{s2, b2}});

        // MVec2dp_U * BiVec2dp
        CHECK(mvec2dp{v1, ps1} * mvec2dp{b2} ==
              mvec2dp{mvec2dp_u{v1, ps1} * mvec2dp{b2}});

        // BiVec2dp * MVec2dp_U
        CHECK(mvec2dp{b1} * mvec2dp{v2, ps2} ==
              mvec2dp{mvec2dp{b1} * mvec2dp_u{v2, ps2}});

        // BiVec2dp * BiVec2dp
        CHECK(mvec2dp{b1} * mvec2dp{b2} == mvec2dp{b1 * b2});

        // BiVec2dp * Vec2dp
        CHECK(mvec2dp{b1} * mvec2dp{v1} == mvec2dp{b1 * v1});

        // Vec2dp * BiVec2dp
        CHECK(mvec2dp{v1} * mvec2dp{b1} == mvec2dp{v1 * b1});

        // Vec2dp * Vec2dp
        CHECK(mvec2dp{v1} * mvec2dp{v2} == mvec2dp{v1 * v2});


        // PScalar2dp * MVec2dp
        CHECK(mvec2dp{ps1} * mvec2dp{s2, v2, b2, ps2} ==
              mvec2dp{ps1 * mvec2dp{s2, v2, b2, ps2}});

        // MVec2dp * Pscalar2dp
        CHECK(mvec2dp{s1, v1, b1, ps1} * mvec2dp{ps2} ==
              mvec2dp{mvec2dp{s1, v1, b1, ps1} * ps2});

        // PScalar2dp * MVec2dp_E
        CHECK(mvec2dp{ps1} * mvec2dp{s2, b2} == mvec2dp{ps1 * mvec2dp_e{s2, b2}});

        // MVec2dp_E * Pscalar2dp
        CHECK(mvec2dp{s1, b1} * mvec2dp{ps2} == mvec2dp{mvec2dp_e{s1, b1} * ps2});

        // PScalar2dp * MVec2dp_U
        CHECK(mvec2dp{ps1} * mvec2dp{v2, ps2} == mvec2dp{ps1 * mvec2dp_u{v2, ps2}});

        // MVec2dp_U * Pscalar2dp
        CHECK(mvec2dp{v1, ps1} * mvec2dp{ps2} == mvec2dp{mvec2dp_u{v1, ps1} * ps2});

        // PScalar2dp * BiVec2dp
        CHECK(mvec2dp{ps1} * mvec2dp{b2} == mvec2dp{ps1 * b2});

        // BiVec2dp * Pscalar2dp
        CHECK(mvec2dp{b1} * mvec2dp{ps2} == mvec2dp{b1 * ps2});

        // PScalar2dp * Vec2dp
        CHECK(mvec2dp{ps1} * mvec2dp{v2} == mvec2dp{ps1 * v2});

        // Vec2dp * Pscalar2dp
        CHECK(mvec2dp{v1} * mvec2dp{ps2} == mvec2dp{v1 * ps2});

        // Pscalar2dp * Pscalar2dp
        CHECK(mvec2dp{ps1} * mvec2dp{ps2} == mvec2dp{ps1 * ps2});

        // Scalar * Scalar
        CHECK(mvec2dp{s1} * mvec2dp{s2} == mvec2dp{s1 * s2});
    }

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
        // auto dot_ab = dot(a, b);
        auto dot_ab = (a >> b);
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
        // auto dot_ab = dot(a, b);
        auto dot_ab = (a << b);
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

        // auto dot_Ab = dot(A, b);
        auto dot_Ab = (A >> b);
        auto wdg_Ab = wdg(A, b);

        // auto dot_aB = dot(a, B);
        auto dot_aB = (a << B);
        auto wdg_aB = wdg(a, B);

        mvec2dp_e ab = a * b;
        mvec2dp abm = mva * mvb;
        mvec2dp abd{scalar2dp{dot_ab}, wdg_ab};

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

        mvec2dp mv5(scalar2dp(5.0));
        mvec2dp mv6(pscalar2dp(6.0));
        mvec2dp mv7{v1};
        mvec2dp mv8{b1};
        mvec2dp mv9{scalar2dp{dot(v1, v3)}, wdg(v1, v3)};

        mvec2dp mv10{v1, pscalar2dp(10.0)};
        // This must not compile! Implict conversion to vec2dp possible
        // possible solution: explicitly deleted constructor for MVec2dp
        // mvec2dp mv11{b1, pscalar2dp_t(10.0)};

        // this does not compile (which is fine, a base cannot convert to derived)
        // mvec2dp mv12{scalar2dp(10.0), v1};

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
        // fmt::println("   normalized_vec = {}", normalize(vec2dp{1.0, 1.0, 0.0}));

        CHECK(std::abs(bulk_nrm(rotate(c, R)) - bulk_nrm(c)) < eps);
        CHECK(gr1(c_rot_m) == rotate(c, R));
        CHECK(rotate(c, R) == rotate_opt(c, R));
        // using a bivector directly:
        CHECK(rotate(vec2dp{1.0, 0.0, 0.0}, rotor(e12_2dp, pi / 4)) ==
              normalize(vec2dp{1.0, 1.0, 0.0}));

        // direct rotation of a bivector
        // CHECK(rotate(bivec2dp{1.0, 0.0, 0.0}, rotor(e31_2dp, pi / 2)) == -e12_2dp);

        // example see Macdonald "Linear and Geometric Algebra", Exercise 7.12, p. 127
        // auto Bv = wdg(e2_2dp,
        //               e1_2dp + std::sqrt(3.0) * e3_2dp); // bivector describing the
        //               plane
        // CHECK(std::abs(nrm(Bv) - 2.0) < eps);
        // CHECK(rotate(Bv, rotor(e31_2dp, pi / 3)) == -2.0 * e12_2dp);

        // just to silence unused variable warnings
        CHECK(my_exp == exp(-B, angle_uv));
        CHECK(my_rot == rotor(B, 2.0 * angle_uv));
        CHECK(rotate(B, my_rot) == rotate_opt(B, my_rot));
        CHECK(c_rot_r == gr1(c_rot_u_r));
        CHECK(angle_c_c_rot == angle(c, c_rot_l));
        CHECK(angle_proj == angle(c_proj, c_rot_proj));
    }

    TEST_CASE("MVec2dp: complement operation")
    {
        fmt::println("MVec2dp: complement operation");

        auto s = scalar2dp{5.0};
        auto v = vec2dp{1.0, 2.0, 3.0};
        auto B = bivec2dp{10.0, 20.0, 30.0};
        auto ps = pscalar2dp{-3.0};

        auto mv = mvec2dp{s, v, B, ps};
        auto mv_e = mvec2dp{s, B};
        auto mv_u = mvec2dp{v, ps};

        auto s2 = scalar2dp{2.5};
        auto v2 = vec2dp{0.5, 1.0, 1.5};
        auto B2 = bivec2dp{5.0, 10.0, 15.0};
        auto ps2 = pscalar2dp{-1.5};

        CHECK(wdg(s, cmpl(s)) == nrm_sq(s) * I_2dp);
        CHECK(wdg(v, cmpl(v)) == nrm_sq(v) * I_2dp);
        CHECK(wdg(B, cmpl(B)) == nrm_sq(B) * I_2dp);
        CHECK(wdg(ps, cmpl(ps)) == nrm_sq(ps) * I_2dp);

        // complement properties
        CHECK(cmpl(cmpl(s)) == s);
        CHECK(cmpl(cmpl(v)) == v);
        CHECK(cmpl(cmpl(B)) == B);
        CHECK(cmpl(cmpl(ps)) == ps);

        // complement values
        CHECK(cmpl(scalar2dp(1.0)) == I_2dp);
        CHECK(cmpl(scalar2dp(1.0)) == rev(scalar2dp(1.0)) * I_2dp);
        CHECK(cmpl(e1_2dp) == e23_2dp);
        CHECK(cmpl(e2_2dp) == e31_2dp);
        CHECK(cmpl(e3_2dp) == e12_2dp);
        CHECK(cmpl(e23_2dp) == e1_2dp);
        CHECK(cmpl(e31_2dp) == e2_2dp);
        CHECK(cmpl(e12_2dp) == e3_2dp);
        CHECK(cmpl(I_2dp) == scalar2d(1.0));
        //
        CHECK(cmpl(cmpl(mv)) == mv);
        CHECK(cmpl(cmpl(mv_e)) == mv_e);
        CHECK(cmpl(cmpl(mv_u)) == mv_u);
        //
        CHECK(wdg(s, cmpl(s)) / nrm_sq(s) == I_2dp);
        CHECK(wdg(cmpl(s), s) / nrm_sq(s) == I_2dp);
        CHECK(wdg(v, cmpl(v)) / nrm_sq(v) == I_2dp);
        CHECK(wdg(cmpl(v), v) / nrm_sq(v) == I_2dp);
        CHECK(wdg(B, cmpl(B)) / nrm_sq(B) == I_2dp);
        CHECK(wdg(cmpl(B), B) / nrm_sq(B) == I_2dp);
        CHECK(wdg(ps, cmpl(ps)) / nrm_sq(ps) == I_2dp);
        CHECK(wdg(cmpl(ps), ps) / nrm_sq(ps) == I_2dp);

        // linearity of the complement operation
        double a = 2.0;
        CHECK(cmpl(a * s) == a * cmpl(s));
        CHECK(cmpl(a * v) == a * cmpl(v));
        CHECK(cmpl(a * B) == a * cmpl(B));
        CHECK(cmpl(a * ps) == a * cmpl(ps));

        CHECK(cmpl(s + s2) == cmpl(s) + cmpl(s2));
        CHECK(cmpl(v + v2) == cmpl(v) + cmpl(v2));
        CHECK(cmpl(B + B2) == cmpl(B) + cmpl(B2));
        CHECK(cmpl(ps + ps2) == cmpl(ps) + cmpl(ps2));

        // check regressive wedge product expressed in terms of
        // wegde product and complement operation
        CHECK(rwdg(B, B2) == cmpl(wdg(cmpl(B), cmpl(B2))));
        CHECK(rwdg(B, v) == cmpl(wdg(cmpl(B), cmpl(v))));
        CHECK(rwdg(v, B) == cmpl(wdg(cmpl(v), cmpl(B))));

        CHECK(cmpl(rwdg(B, B2)) == wdg(cmpl(B), cmpl(B2)));
        CHECK(cmpl(wdg(v, v2)) == rwdg(cmpl(v), cmpl(v2)));
    }

    TEST_CASE("2dp: objects - convenience types")
    {
        fmt::println("2dp: objects - convenience types");

        auto v2d = Vec2d<double>(2.0, 1.5);
        auto p2d = Point2d<double>(2.0, 1.5);
        auto v = Vec2dp<double>(6.0, 4.5, 3.0);
        auto p = Point2dp<double>(6.0, 4.5, 3.0);
        auto p1 = Point2d<double>(1.0, 2.0);
        auto p2 = Point2d<double>(2.0, 3.0);
        auto l1 = wdg(p1, p2);

        // a Point2d is-a Vec2dp, thus all operations defined for Vec2dp should work
        // directly for Point2dp
        CHECK(dot(v, v) == dot(p, p));
        CHECK(p.unitize() == unitize(Point2dp(v)));
        CHECK(p2d == Point2d(v2d));
        CHECK(p.unitize() == Point2dp(p2d));
        CHECK(p.unitize() == Point2dp(v2d));
        CHECK(dot(BiVec2dp<double>(1, 7, 3), BiVec2dp<double>(4, 5, 6)) == -18.0);

        auto vector = Vector2d<double>(1.0, 2.0);

        // fmt::println("vector = {}", vector);

        CHECK(nrm_sq(vector) == 5.0);

        //
        // fmt::println("p1 = {}", p1);
        // fmt::println("p2 = {}", p2);
        // fmt::println("l1 = {}", l1);
        //
        // Line2d ctors
        auto l2 = Line2d(p1, p2);
        auto l3 = Line2d(BiVec2dp<double>{-1, 1, -1});
        auto l4 = Line2d(p1, Vec2d<double>{1, 1});
        auto l5 = Line2d<double>(-1, 1, -1);
        auto l6 = Line2d<double>();
        CHECK(l1 == l2);
        CHECK(l1 == l3);
        CHECK(l1 == l4);
        CHECK(l1 == l5);
        CHECK(l6 == BiVec2dp<double>(0, 0, 0));
    }

    TEST_CASE("2dp: objects - bulk and weight")
    {
        fmt::println("2dp: objects - bulk and weight");

        auto p1 = vec2dp{-2.0, 1.0, 1.0};
        auto p2 = vec2dp{2.0, 1.0, 2.0};
        auto p3 = vec2dp{2.0, 1.0, 1.0};

        auto l1 = bivec2dp{2.0, 1.0, 1.0};
        auto l2 = bivec2dp{-2.0, 1.0, 3.0};
        auto l3 = wdg(p1, p3);

        // bulk and weight
        CHECK(bulk(p1) == vec2dp{p1.x, p1.y, 0.0});
        CHECK(weight(p1) == vec2dp{0.0, 0.0, p1.z});
        CHECK(bulk(l1) == bivec2dp{0.0, 0.0, l1.z});
        CHECK(weight(l1) == bivec2dp{l1.x, l1.y, 0.0});

        // bulk_nrm_sq and weight_nrm_sq
        CHECK(bulk_nrm_sq(p1) == p1.x * p1.x + p1.y * p1.y);
        CHECK(weight_nrm_sq(p1) == p1.z * p1.z);
        CHECK(bulk_nrm_sq(l1) == l1.z * l1.z);
        CHECK(weight_nrm_sq(l1) == l1.x * l1.x + l1.y * l1.y);

        // bulk_nrm and weight_nrm
        CHECK(bulk_nrm(p1) == std::sqrt(p1.x * p1.x + p1.y * p1.y));
        CHECK(weight_nrm(p1) == std::sqrt(p1.z * p1.z));
        CHECK(bulk_nrm(l1) == std::sqrt(l1.z * l1.z));
        CHECK(weight_nrm(l1) == std::sqrt(l1.x * l1.x + l1.y * l1.y));

        // geom_nrm
        CHECK(geom_nrm(p1).c0 / geom_nrm(p1).c1 == bulk_nrm(p1) / weight_nrm(p1));
        CHECK(geom_nrm(l1).c0 / geom_nrm(l1).c1 == bulk_nrm(l1) / weight_nrm(l1));
        CHECK(3.0 * geom_nrm(l1).c0 == geom_nrm(l2).c0);
        CHECK(geom_nrm(l1).c1 == geom_nrm(l2).c1);

        // attitude
        CHECK(att(p1) == p1.z);
        CHECK(att(p1) == rwdg(p1, cmpl(e3_2dp)));
        CHECK(att(p2) == p2.z);
        CHECK(att(p2) == rwdg(p2, cmpl(e3_2dp)));
        CHECK(att(l1) == vec2dp{l1.y, -l1.x, 0.0});
        CHECK(att(l1) == rwdg(l1, cmpl(e3_2dp)));
        CHECK(att(l2) == rwdg(l2, cmpl(e3_2dp)));
        CHECK(att(l3) == rwdg(l3, cmpl(e3_2dp)));

        // intersections
        auto l4 = bivec2dp{1.0, 1.0, 3.0};
        auto l5 = bivec2dp{1.0, 1.0, 0.0};
        auto int_sec =
            rwdg(l4, l5); // should intersect at infinity, i.e. att(int_sec) == 0.0
        // fmt::println("   att(l4) = {}", att(l4));
        // fmt::println("   att(l5) = {}", att(l5));
        // fmt::println("   normalize(int_sec(l4,l5)) = {}", normalize(int_sec));
        CHECK(att(int_sec) == 0.0);

        // fmt::println("   att(l1) = {}", att(l1));
        // fmt::println("   att(l2) = {}", att(l2));
        // fmt::println("   unitize(rwdg(l1, l2)) = {}", unitize(rwdg(l1, l2)));
        CHECK(unitize(rwdg(l1, l2)) == Vec2dp{0.5, -2.0, 1.0});
    }

    TEST_CASE("2dp: objects - euclidean distance")
    {
        fmt::println("2dp: objects - euclidean distance");

        auto p1 = vec2dp{-1.0, -1.0, 1.0};
        auto p2 = vec2dp{1.0, -1.0, 1.0};
        auto p3 = vec2dp{1.0, 1.0, 1.0};

        auto l1 = wdg(p1, p2);
        auto l2 = wdg(p2, p3);
        auto l3 = wdg(p1, p3);


        auto dp1p2 = euclidean_distance2dp(p1, p2);
        auto dp2p3 = euclidean_distance2dp(p2, p3);
        auto dp1p3 = euclidean_distance2dp(p1, p3);

        auto dp3l1 = euclidean_distance2dp(p3, l1);
        auto dp1l2 = euclidean_distance2dp(p1, l2);
        auto dp2l3 = euclidean_distance2dp(p2, l3);

        // fmt::println("");
        // fmt::println("dp1p2 = {}, dist = {}", dp1p2, dp1p2.c0 / dp1p2.c1);
        // fmt::println("dp2p3 = {}, dist = {}", dp2p3, dp2p3.c0 / dp2p3.c1);
        // fmt::println("dp1p3 = {}, dist = {}", dp1p3, dp1p3.c0 / dp1p3.c1);
        // fmt::println("");
        // fmt::println("dp3l1 = {}, dist = {}", dp3l1, dp3l1.c0 / dp3l1.c1);
        // fmt::println("dp1l2 = {}, dist = {}", dp1l2, dp1l2.c0 / dp1l2.c1);
        // fmt::println("dp2l3 = {}, dist = {}", dp2l3, dp2l3.c0 / dp2l3.c1);
        // fmt::println("");

        CHECK(dp1p2 == dp2p3);
        CHECK(dp3l1 == dp1l2);
        CHECK(dp1p3.c0 / dp1p3.c1 - 2.0 * std::sqrt(2.0) < eps);
        CHECK(dp2l3.c0 / dp2l3.c1 + std::sqrt(2.0) < eps);

        // auto res = rwdg(pscalar2dp(2.5), horizon_2dp);
        // fmt::println("res = {}", res);
    }

    TEST_CASE("pga_2dp<2,0,1> - product tests")
    {
        fmt::println("pga_2dp<2,0,1> - product tests");
        // tests based on "The inner products of geometric algebra", Leo Dorst


        auto s1 = scalar2dp{2.0};
        auto v1 = vec2dp{1.0, -3.0, 0.0};
        // auto v1 = vec2dp{1.0, 0.0, 0.0};
        auto b1 = bivec2dp{2.0, -4.0, 1.0};
        auto ps1 = pscalar2dp{-2.0};

        auto s2 = scalar2dp{-1.0};
        auto v2 = vec2dp{2.0, 1.5, 0.0};
        // auto v2 = vec2dp{0.0, 1.0, 0.0};
        auto b2 = bivec2dp{-2.0, 10.5, 3.0};
        auto ps2 = pscalar2dp{3.0};

        auto s3 = scalar2dp{-10.0};
        auto v3 = vec2dp{-2.0, 4.0, -10.0};
        auto b3 = bivec2dp{1.0, -2.0, 3.0};
        auto ps3 = pscalar2dp{-12.0};

        // dot product
        CHECK(dot(scalar2dp(1.0), scalar2dp(1.0)) == scalar2dp(1.0));
        CHECK(dot(e1_2dp, e1_2dp) == scalar2dp(1.0));
        CHECK(dot(e2_2dp, e2_2dp) == scalar2dp(1.0));
        CHECK(dot(e3_2dp, e3_2dp) == scalar2dp(0.0));
        CHECK(dot(e23_2dp, e23_2dp) == scalar2dp(0.0));
        CHECK(dot(e31_2dp, e31_2dp) == scalar2dp(0.0));
        CHECK(dot(e12_2dp, e12_2dp) == scalar2dp(-1.0)); // differs from Lengyel!
        CHECK(dot(pscalar2dp(1.0), pscalar2dp(1.0)) == scalar2dp(0.0));

        // regressive dot product
        CHECK(rdot(scalar2dp(1.0), scalar2dp(1.0)) == pscalar2dp(0.0));
        CHECK(rdot(e1_2dp, e1_2dp) == pscalar2dp(0.0));
        CHECK(rdot(e2_2dp, e2_2dp) == pscalar2dp(0.0));
        CHECK(rdot(e3_2dp, e3_2dp) == pscalar2dp(-1.0));
        CHECK(rdot(e23_2dp, e23_2dp) == pscalar2dp(1.0));
        CHECK(rdot(e31_2dp, e31_2dp) == pscalar2dp(1.0));
        CHECK(rdot(e12_2dp, e12_2dp) == pscalar2dp(0.0));
        CHECK(rdot(pscalar2dp(1.0), pscalar2dp(1.0)) == pscalar2dp(1.0));


        // wedge product
        // check full permissible range of arguments, even the ones delivering 0)
        CHECK(wdg(b1, b2) == scalar2dp(0.0));
        CHECK(wdg(v1, ps1) == scalar2dp(0.0));
        CHECK(wdg(ps1, v1) == scalar2dp(0.0));
        CHECK(wdg(b1, ps1) == scalar2dp(0.0));
        CHECK(wdg(ps1, b1) == scalar2dp(0.0));
        CHECK(wdg(ps1, ps1) == scalar2dp(0.0));

        // fmt::println("   b1           = {}", b1);
        // fmt::println("   b2           = {}", b2);
        // fmt::println("   b1*b2        = {}", b1 * b2);
        // fmt::println("   cmt(b1, b2)  = {}", cmt(b1, b2));
        // fmt::println("   b1<<b2       = {}", b1 << b2);
        CHECK(wdg(b1, b2) == nrm_sq(b1 * b2 - cmt(b1, b2) - (b1 << b2)));
        CHECK(wdg(v1, ps1) == nrm_sq(0.5 * (v1 * ps1 + rev(ps1) * v1)));
        CHECK(wdg(ps1, v1) == nrm_sq(0.5 * (ps1 * v1 + v1 * rev(ps1))));


        CHECK((v1 << wdg(v2, v3)) == wdg(v1 << v2, v3) + wdg(gr_inv(v2), v1 << v3));
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

        // 3rd with vectors and bivectors directly
        CHECK(v1 * b1 == (v1 << b1) + wdg(v1, b1));
        CHECK(b1 * v1 == (b1 >> v1) + wdg(b1, v1));
        CHECK((b1 >> v1) == -(v1 << gr_inv(b1)));

        // fmt::println("   v1 << b1        = {}", v1 << b1);
        // fmt::println("   v1 * b1         = {}", v1 * b1);
        // fmt::println("   gr_inv(b1) * v1 = {}", gr_inv(b1) * v1);
        CHECK((v1 << b1) == gr1(0.5 * (v1 * b1 - gr_inv(b1) * v1)));
        CHECK((b1 >> v1) == gr1(0.5 * (b1 * v1 - v1 * gr_inv(b1))));

        // fmt::println("   wdg(v1, b1)     = {}", wdg(v1, b1));
        // fmt::println("   v1 * b1         = {}", v1 * b1);
        // fmt::println("   gr_inv(b1) * v1 = {}", gr_inv(b1) * v1);
        CHECK(wdg(v1, b1) == gr3(0.5 * (v1 * b1 + gr_inv(b1) * v1)));
        CHECK(wdg(b1, v1) == gr3(0.5 * (b1 * v1 + v1 * gr_inv(b1))));

        CHECK(wdg(v1, (b1 << b2)) == ((v1 << b1) << b2) + (gr_inv(b1) << wdg(v1, b2)));

        // cross-check direct implementation of rwdg by comparing with wdg
        CHECK(rwdg(b1, b2) == cmpl(wdg(cmpl(b1), cmpl(b2))));
    }

    TEST_CASE(
        "G<2,0,1> - pga2dp simple applications, complements, contraction, expansions")
    {
        fmt::println("G<2,0,1> - pga2dp simple applications, complements, contraction, "
                     "expansions");

        auto s1 = scalar2dp{2.0};
        auto v1 = vec2dp{1.0, -3.0, 0.0};
        auto b1 = bivec2dp{2.0, -4.0, 1.0};
        auto ps1 = pscalar2dp{-2.0};
        auto M1 = mvec2dp{s1, v1, b1, ps1};


        // simple projections
        auto v = vec2dp{4.0, 3.5, 0.0};
        auto u = vec2dp{1.0, 2.0, 0.0};
        auto B = e12_2dp;


        auto v_in_u = project_onto(v, u);
        auto v_perp_u = reject_from(v, u);

        auto v_in_B = project_onto(v, B);
        auto v_perp_B = reject_from(v, B);

        // fmt::println("   v         = {:g}", v);
        // fmt::println("   u         = {:g}", u);
        // fmt::println("   B         = {:g}", B);
        // fmt::println("   cmpl(v)   = {:g}", cmpl(v));
        // fmt::println("   cmpl(B)   = {:g}", cmpl(B));
        // fmt::println("");
        // fmt::println("   v_in_u    = {:g}", v_in_u);
        // fmt::println("   v_perp_u  = {:g}", v_perp_u);
        // fmt::println("   v << u    = {:g}", v << u);
        // fmt::println("   u >> v    = {:g}", u >> v);
        // fmt::println("");
        // fmt::println("");
        // fmt::println("   v_in_B    = {:g}", v_in_B);
        // fmt::println("   v_perp_B  = {:g}", v_perp_B);
        // fmt::println("   v << B    = {:g}", v << B);
        // fmt::println("   B >> v    = {:g}", B >> v);
        // fmt::println("");

        // x^B = 0 for every point x in B
        CHECK(nrm_sq(wdg(v_in_B, B)) < eps);

        // v_perp_B should be proportional to the normal vector n of B
        // n = cmpl(B)
        // thus, wdg(v_perp_B, n) == 0 is required
        CHECK(nrm_sq(wdg(v_perp_B, cmpl(B))) < eps);

        // v_in_B and v_perp_B should be perpendicular to each other
        CHECK(nrm_sq(dot(v_in_B, v_perp_B)) < eps);

        // v should be the sum of v_in_B and v_perp_B
        CHECK(v == v_in_B + v_perp_B);

        // v should be the sum of v_in_u and v_perp_u
        CHECK(v == v_in_u + v_perp_u);

        // check complements
        CHECK(cmpl(scalar2dp(1.0)) == I_2dp);
        CHECK(cmpl(e1_2dp) == e23_2dp);
        CHECK(cmpl(e2_2dp) == e31_2dp);
        CHECK(cmpl(e3_2dp) == e12_2dp);
        CHECK(cmpl(e23_2dp) == e1_2dp);
        CHECK(cmpl(e31_2dp) == e2_2dp);
        CHECK(cmpl(e12_2dp) == e3_2dp);
        CHECK(cmpl(I_2dp) == scalar2d(1.0));
        //
        CHECK(cmpl(cmpl(M1)) == M1);
        CHECK(cmpl(cmpl(mvec2dp_e(s1, b1))) == mvec2dp_e(s1, b1));
        CHECK(cmpl(cmpl(mvec2dp_u(v1, ps1))) == mvec2dp_u(v1, ps1));
        //
        CHECK(wdg(scalar2dp(5), cmpl(scalar2dp(5))) / nrm_sq(scalar2dp(5)) == I_2dp);
        CHECK(wdg(cmpl(scalar2dp(5)), scalar2dp(5)) / nrm_sq(scalar2dp(5)) == I_2dp);
        CHECK(wdg(v, cmpl(v)) / nrm_sq(v) == I_2dp);
        CHECK(wdg(cmpl(v), v) / nrm_sq(v) == I_2dp);
        CHECK(wdg(B, cmpl(B)) / nrm_sq(B) == I_2dp);
        CHECK(wdg(cmpl(B), B) / nrm_sq(B) == I_2dp);
        CHECK(wdg(pscalar2dp(3), cmpl(pscalar2dp(3))) / nrm_sq(pscalar2dp(3)) == I_2dp);
        CHECK(wdg(cmpl(pscalar2dp(3)), pscalar2dp(3)) / nrm_sq(pscalar2dp(3)) == I_2dp);

        // check contractions: <<, >> and rwdg( u, compl(v) )
        // fmt::println("   v         = {:g}", v);
        // fmt::println("   cmpl(v)   = {:g}", cmpl(v));
        // fmt::println("");
        // fmt::println("   v << B = {:g}", v << B);
        // fmt::println("   B >> v = {:g}", B >> v);
        // fmt::println("");
        // fmt::println("   rwdg(B, cmpl(v))  = {:g}", rwdg(B, cmpl(v)));
        // fmt::println("   rwdg(cmpl(v), B)  = {:g}", rwdg(cmpl(v), B));
        // fmt::println("");
        CHECK((v << B) == rwdg(B, cmpl(v)));
        CHECK((B >> v) == rwdg(cmpl(v), B));

        // check expansions: v ^ cmpl(B)  and  cmpl(B) ^ v
        // (create new bivector that contains v and is perpendicular to B)
        //
        // also checks the duality correspondence:
        //      cmpl(v >> B) == v ^ cmpl(B)
        //      cmpl(v >> B) == cmpl(B) ^ v
        // fmt::println("   v      = {:g}", v);
        // fmt::println("   B      = {:g}", B);
        // fmt::println("   v << B = {:g}", v << B);
        // fmt::println("   B >> v = {:g}", B >> v);
        // fmt::println("");
        // fmt::println("   wdg(cmpl(B), v)         = {:g}", wdg(cmpl(B), v));
        // fmt::println("   n=cmpl(wdg(cmpl(B), v)) = {:g}", cmpl(wdg(cmpl(B), v)));
        // fmt::println("");
        // fmt::println("   wdg(v, cmpl(B))         = {:g}", wdg(v, cmpl(B)));
        // fmt::println("   n=cmpl(wdg(v, cmpl(B))) = {:g}", cmpl(wdg(v, cmpl(B))));
        // vector is in plane defined by the expansion
        CHECK(wdg(v, wdg(v, cmpl(B))) == 0.0);
        CHECK(wdg(v, wdg(cmpl(B), v)) == 0.0);
        // duality of the contraction and the wedge product (based on complement)
        CHECK(cmpl(B >> v) == wdg(v, cmpl(B)));
        CHECK(cmpl(v << B) == wdg(cmpl(B), v));
    }

    TEST_CASE("G<2,0,1> - pga2dp join and meet (wdg, rwdg)")
    {
        fmt::println("G<2,0,1> - pga2dp join and meet (wdg, rwdg)");

        auto p1 = point2d{2, 0};
        auto p2 = point2d{4, 3};
        auto l12 = join(p1, p2); // = wdg(p1, p2)

        auto p3 = point2d{0, 6};
        auto p4 = point2d{2, 6};
        auto l34 = join(p3, p4); // = wdg(p3, p4)

        auto p1p = point2dp{p1}; // all works with projective points as well
        auto p2p = point2dp{p2};
        auto l12p = join(p1p, p2p);
        auto p3p = point2dp{p3};
        auto p4p = point2dp{p4};
        auto l34p = join(p3p, p4p);

        auto p5 = meet(l12, l34);
        auto p5p = meet(l12p, l34p);

        auto p6 = point2d{0, 3}; // the expansion: line perpendicular to l12 through p6
        auto lp6 = expand(p6, l12);

        // fmt::println("");
        // fmt::println("p1 = {}, p2 = {}", p1, p2);
        // fmt::println("l12  = join(p1, p2)   = wdg(p1, p2) = {}", l12);
        // fmt::println("att(l12) = {}", att(l12));
        // fmt::println("");
        // fmt::println("l12p = join(p1p, p2p) = wdg(p1p, p2p) = {}", l12p);
        // fmt::println("att(l12p) = {}", att(l12p));
        // fmt::println("");
        // fmt::println("p5 = meet(l12, l34) = {}", p5);
        // fmt::println("p5p = meet(l12p, l34p) = {}", p5p);
        // fmt::println("");
        // fmt::println("lp6  = {}", lp6);
        // fmt::println("att(lp6) = {}", att(lp6));
        // fmt::println("");

        auto ln1 = join(point2d{0, 3}, point2d{3, 3}).unitize();
        auto ln2 = join(point2d{0, 6}, point2d{3, 6}).unitize();
        auto pt12 = meet(ln1, ln2); // gives direction to intersection at infinity
                                    // weighted by the distance of the lines

        auto ln1p = unitize(join(point2dp{0, 3, 1}, point2dp{3, 3, 1}));
        auto ln2p = unitize(join(point2dp{0, 6, 1}, point2dp{3, 6, 1}));
        auto pt12p = meet(ln1p, ln2p);

        // fmt::println("");
        // fmt::println("ln1      = {}, ln2      = {}", ln1, ln2);
        // fmt::println("att(ln1) = {}, att(ln2) = {}", att(ln1), att(ln2));
        // fmt::println("pt12  = meet(ln1, ln2) = {}, att(pt12) = {}", pt12, att(pt12));
        // fmt::println("");
        // fmt::println("ln1p      = {}, ln2p      = {}", ln1p, ln2p);
        // fmt::println("att(ln1p) = {}, att(ln2p) = {}", att(ln1p), att(ln2p));
        // fmt::println("pt12p  = meet(ln1p, ln2p) = {}", pt12p);
        // fmt::println("");

        CHECK(join(p1, p2) == wdg(p1, p2));
        CHECK(point2dp(p5) == point2dp(p5p));
        CHECK(point2dp(unitize(meet(l12, l34))) == point2dp(unitize(rwdg(l12, l34))));
        CHECK(dot(att(l12), att(lp6)) == 0);

        CHECK(pt12 == pt12p);
    }


    //////////////////////////////////////////////////////////////////////////////////////
    // projective geometric algebra 3d: pga3dp (embedded in a 4d representational space)
    //////////////////////////////////////////////////////////////////////////////////////


    TEST_CASE("algebra<3, 0, 1> pga3dp")
    {
        fmt::println("algebra<3, 0, 1> - pga3dp:");
        // 3d projective geometric algebra (3dp) models three-dimensional Euclidean
        // algebra embedding it in a four-dimensional representational space
        const algebra<3, 0, 1> alg;
        CHECK(alg.p() == 3);
        CHECK(alg.n() == 0);
        CHECK(alg.z() == 1);
        CHECK(alg.dim_space() == 4);                 // dim_space == p+n+z
        CHECK(alg.num_components() == 16);           // num_components == 2^dim
        CHECK(alg.num_components_grade.size() == 5); // == dim_space + 1
        fmt::println("   pga3dp: dim_grade = {}",
                     fmt::join(alg.num_components_grade, ", "));
        fmt::println("   pga3dp: basis_name = {}", fmt::join(alg.basis_name, ", "));
    }

    TEST_CASE("G<3,0,1> - pga3dp defining basic types and ctor checks")
    {
        fmt::println("G<3,0,1>: defining basic types and ctor checks");

        auto s = scalar3dp{5.0};
        auto v = vec3dp{1.0, 2.0, 3.0, 1.0};
        auto B = bivec3dp{-1.0, 2.0, 1.0, -10.0, 20.0, 10.0};
        auto t = trivec3dp{3.0, 6.0, 9.0, 3.0};
        auto ps = pscalar3dp{-5.0};

        auto mv0 = mvec3dp{s};
        auto mv1 = mvec3dp{v};
        auto mv2 = mvec3dp{B};
        auto mv3 = mvec3dp{t};
        auto mv4 = mvec3dp{ps};

        auto mv5a = mvec3dp{s, B, ps};
        auto mv5b = mvec3dp_e{s};
        auto mv5c = mvec3dp_e{B};
        auto mv5d = mvec3dp_e{ps};
        auto mv5e = mvec3dp_e{s, B, ps};
        auto mv5f{mv5e};           // cp assign to even grade
        auto mv5g = mvec3dp{mv5e}; // assign to mv

        auto mv5h = mvec3dp_e{s, B};
        auto mv5i = mvec3dp_e{B, ps};
        auto mv5j = mvec3dp_e{s, ps};

        auto mv6a = mvec3dp{v, t};
        auto mv6b = mvec3dp_u{v};
        auto mv6c = mvec3dp_u{t};
        auto mv6d = mvec3dp_u{v, t};
        auto mv6e{mv6d};           // cp assign to uneven grade
        auto mv6f = mvec3dp{mv6d}; // assign to mv

        auto mv7 = mvec3dp{s, v, B, t, ps};

        // fmt::println("   mv1  = {}", mv1);
        // fmt::println("   mv2  = {}", mv2);
        // fmt::println("   mv3  = {}", mv3);
        // fmt::println("   mv4  = {}", mv4);
        // fmt::println("   mv5a = {}", mv5a);
        // fmt::println("   mv5b = {}", mv5b);
        // fmt::println("   mv5c = {}", mv5c);
        // fmt::println("   mv5d = {}", mv5d);
        // fmt::println("   mv5e = {}", mv5e);
        // fmt::println("   mv5f = {}", mv5f);
        // fmt::println("   mv5g = {}", mv5g);
        // fmt::println("   mv6a = {}", mv6a);
        // fmt::println("   mv6b = {}", mv6b);
        // fmt::println("   mv6c = {}", mv6c);
        // fmt::println("   mv6d = {}", mv6d);
        // fmt::println("   mv6e = {}", mv6e);
        // fmt::println("   mv6f = {}", mv6f);
        // fmt::println("   mv7  = {}", mv7);

        CHECK(gr0(mv0) == s);
        CHECK(gr1(mv1) == v);
        CHECK(gr2(mv2) == B);
        CHECK(gr3(mv3) == t);
        CHECK(gr4(mv4) == ps);

        CHECK(gr0(mv5a) == s);
        CHECK(gr1(mv5a) == vec3dp{});
        CHECK(gr2(mv5a) == B);
        CHECK(gr3(mv5a) == trivec3dp{});
        CHECK(gr4(mv5a) == ps);

        CHECK(gr0(mv5b) == s);
        CHECK(gr2(mv5b) == bivec3dp{});
        CHECK(gr4(mv5b) == pscalar3dp{});

        CHECK(gr0(mv5c) == scalar3dp{});
        CHECK(gr2(mv5c) == B);
        CHECK(gr4(mv5c) == pscalar3dp{});

        CHECK(gr0(mv5d) == scalar3dp{});
        CHECK(gr2(mv5d) == bivec3dp{});
        CHECK(gr4(mv5d) == ps);

        CHECK(gr0(mv5e) == s);
        CHECK(gr2(mv5e) == B);
        CHECK(gr4(mv5e) == ps);

        CHECK(mv5f == mv5e);

        CHECK(gr0(mv5h) == s);
        CHECK(gr2(mv5h) == B);
        CHECK(gr4(mv5h) == pscalar3dp{});

        CHECK(gr0(mv5i) == scalar3dp{});
        CHECK(gr2(mv5i) == B);
        CHECK(gr4(mv5i) == ps);

        CHECK(gr0(mv5j) == s);
        CHECK(gr2(mv5j) == bivec3dp{});
        CHECK(gr4(mv5j) == ps);

        CHECK(gr0(mv6a) == scalar3dp{});
        CHECK(gr1(mv6a) == v);
        CHECK(gr2(mv6a) == bivec3dp{});
        CHECK(gr3(mv6a) == t);
        CHECK(gr4(mv6a) == pscalar3dp{});

        CHECK(gr1(mv6b) == v);
        CHECK(gr3(mv6b) == trivec3dp{});

        CHECK(gr1(mv6c) == vec3dp{});
        CHECK(gr3(mv6c) == t);

        CHECK(gr1(mv6d) == v);
        CHECK(gr3(mv6d) == t);

        CHECK(mv6e == mv6d);

        CHECK(gr0(mv6f) == scalar3dp{});
        CHECK(gr1(mv6f) == v);
        CHECK(gr2(mv6f) == bivec3dp{});
        CHECK(gr3(mv6f) == t);
        CHECK(gr4(mv6f) == pscalar3dp{});

        CHECK(gr0(mv7) == s);
        CHECK(gr1(mv7) == v);
        CHECK(gr2(mv7) == B);
        CHECK(gr3(mv7) == t);
        CHECK(gr4(mv7) == ps);

        CHECK(gr(s) == 0);
        CHECK(gr(v) == 1);
        CHECK(gr(B) == 2);
        CHECK(gr(t) == 3);
        CHECK(gr(ps) == 4);
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3dp<T> basic test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3dp: default init")
    {
        fmt::println("Vec3dp: default init");
        vec3dp v;
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.x) < eps);
        CHECK(std::abs(v.y) < eps);
        CHECK(std::abs(v.z) < eps);
        CHECK(std::abs(v.w) < eps);
    }

    TEST_CASE("Vec3dp: with curly braced intializer")
    {
        fmt::println("Vec3dp: with curly braced intializer");
        vec3dp v{0.0, 0.0, 0.0, 0.0};
        // fmt::println("   v = {}", v);
        CHECK(std::abs(v.x) < eps);
        CHECK(std::abs(v.y) < eps);
        CHECK(std::abs(v.z) < eps);
        CHECK(std::abs(v.w) < eps);
    }

    TEST_CASE("Vec3dp: cp ctor & cp assign incl. type deduction")
    {
        fmt::println("Vec3dp: cp ctor & cp assign incl. type deduction");
        vec3dp v1{1.0, 2.0, 3.0, 1.0}; // init with double (type deduction)
        vec3dp v2{v1};                 // cp ctor
        vec3dp v3 = v2;                // cp assign
        vec3dp v4 = -v2;               // cp assign with unary minus

        // fmt::println("   v1 = {}", v1);
        // fmt::println("   v2 = {}", v2);
        // fmt::println("   v3 = {}", v3);
        // fmt::println("   v4 = {}", v4);

        CHECK(std::abs(v1.x - 1.0) < eps);
        CHECK(std::abs(v1.y - 2.0) < eps);
        CHECK(std::abs(v1.z - 3.0) < eps);
        CHECK(std::abs(v1.w - 1.0) < eps);
        CHECK(std::abs(v2.x - 1.0) < eps);
        CHECK(std::abs(v2.y - 2.0) < eps);
        CHECK(std::abs(v2.z - 3.0) < eps);
        CHECK(std::abs(v2.w - 1.0) < eps);
        CHECK(std::abs(v3.x - 1.0) < eps);
        CHECK(std::abs(v3.y - 2.0) < eps);
        CHECK(std::abs(v3.z - 3.0) < eps);
        CHECK(std::abs(v3.w - 1.0) < eps);
        CHECK(v4 == -v2);

        // check direct assignment operators (sequence of tests decisive!)
        CHECK((v3 += v2) == 2.0 * v1);
        CHECK((v3 -= v1) == v1);
        CHECK((v3 *= 2.0) == 2.0 * v1);
        CHECK((v3 /= 2.0) == v1);
    }

    TEST_CASE("Vec3dp: fmt & cout printing")
    {
        fmt::println("Vec3dp: fmt & cout printing");

        vec3dp pf{1.0f, 2.00001f, 3.0f, 1.0f};
        vec3dp pd{1.0, 2.00001, 3.0, 1.0};

        // std::cout << "       cout: pf = " << pf << std::endl;
        fmt::println("       fmt:  pf = {}", pf);
        fmt::println("       fmt:  pf = {:.8f}", pf);

        // std::cout << "       cout: pd = " << pd << std::endl;
        fmt::println("       fmt:  pd = {}", pd);
        fmt::println("       fmt:  pd = {:.8f}", pd);

        std::vector<Vec3dp<double>> vp1{{1.0, 1.0, 1.0, 1.0}, {1.5, 2.0, 3.0, 1.0}};
        fmt::println("       fmt: vp1 = {}", fmt::join(vp1, ", "));
        fmt::println("       fmt: vp1 = {:e}", fmt::join(vp1, ", "));
        fmt::println("");

        CHECK(nrm_sq(pf - pd) < eps);
    }

    TEST_CASE("Vec3dp: comparison float")
    {
        fmt::println("Vec3dp: comparison float");

        Vec3dp<float> v1f{1.0f, 2.0f, 3.0f, 1.0f};
        Vec3dp<float> v2f{2.0f, 4.0f, 3.0f, 1.0f};
        Vec3dp<float> v3f{1.0f, 2.0000001f, 3.0f, 1.0f};
        Vec3dp<float> v4f{v1f};

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

    TEST_CASE("Vec3dp: comparison double")
    {
        fmt::println("Vec3dp: comparison double");

        Vec3dp<double> v1d{1.0, 2.0, 3.0, 1.0};
        Vec3dp<double> v2d{2.0, 4.0, 3.0, 1.0};
        Vec3dp<double> v3d{1.0, 2.0000000000000001, 3.0, 1.0};
        Vec3dp<double> v4d{v1d};

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

    TEST_CASE("Vec3dp: vector space and linearity tests")
    {
        fmt::println("Vec3dp: vector space and linearity tests");

        // a vector space has scalar multiplication and vector addition defined
        // and is closed under these operations
        //
        // a (linear) vector space fulfills operations tested against below:

        vec3dp p0;
        vec3dp p1{1.0, 2.0, 3.0, 1.0};
        vec3dp p2{2.0, 4.0, 6.0, 2.0};
        vec3dp p3{3.0, 6.0, 9.0, 3.0};
        vec3dp p4 = -p1; // assignment using unary minus
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

    TEST_CASE("Vec3dp: inner product properties")
    {
        fmt::println("Vec3dp: inner product properties");

        double a = 2.35;
        vec3dp u{1.0, 2.0, 1.0, 1.0};
        vec3dp v{-0.5, 3.0, 0.5, 1.0};
        vec3dp w{3.0, 6.0, -3.0, 1.0};

        CHECK(dot(a * u, v) == a * dot(u, v));
        CHECK(dot(u + v, w) == dot(u, w) + dot(v, w));
        CHECK(dot(u, v) == dot(v, u));
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Vec3dp<T> operations test cases
    ////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Vec3dp: operations - norm, inverse, dot")
    {
        fmt::println("Vec3dp: operations - norm, inverse, dot");

        vec3dp v1{2.0, 1.0, 2.0, -2.0};
        vec3dp v1n{normalize(v1)};
        vec3dp v1u{unitize(v1)};

        vec3dp v3{2.0, 6.0, -4.0, 2.0};
        vec3dp v4{inv(v3)};

        auto v3m = mvec3dp{v3};
        auto v4m = mvec3dp{v4};

        // fmt::println("v1                  = {: .4f}, nrm(v1)         = {: .4f}", v1,
        //              nrm(v1));
        // fmt::println("v1                  = {: .4f}, bulk_nrm(v1)    = {: .4f}", v1,
        //              bulk_nrm(v1));
        // fmt::println("v1                  = {: .4f}, weight_nrm(v1)  = {: .4f}", v1,
        //              weight_nrm(v1));
        // fmt::println("v1n = normalize(v1) = {: .4f}, nrm(v1n)        = {: .4f}", v1n,
        //              nrm(v1n));
        // fmt::println("v1n                 = {: .4f}, bulk_nrm(v1n)   = {: .4f}", v1n,
        //              bulk_nrm(v1n));
        // fmt::println("v1n                 = {: .4f}, weight_nrm(v1n) = {: .4f}", v1n,
        //              weight_nrm(v1n));
        // fmt::println("v1u = unitize(v1)   = {: .4f}, nrm(v1u)        = {: .4f}", v1u,
        //              nrm(v1u));
        // fmt::println("v1u                 = {: .4f}, bulk_nrm(v1u)   = {: .4f}", v1u,
        //              bulk_nrm(v1u));
        // fmt::println("v1u                 = {: .4f}, weight_nrm(v1u) = {: .4f}", v1u,
        //              weight_nrm(v1u));
        // fmt::println("");
        // fmt::println("v3                  = {: .4f}, nrm(v1) = {: .4f}", v3, nrm(v3));
        // fmt::println("v4 = inv(v3)        = {: .4f}, "
        //              "nrm(v3) = {: .4f}, nrm(v3)*nrm(v4) = {: .4f}",
        //              v4, nrm(v4), nrm(v3) * nrm(v4));
        // fmt::println("");
        // fmt::println("v3m = {}, nrm(v3m) = {}", v3m, nrm(v3m));
        // fmt::println("v4m = {}, nrm(v4m) = {}", v4m, nrm(v4m));
        // fmt::println("0.5*(v3m*v4m + v4m*v3m) = {}", 0.5 * (v3m * v4m + v4m * v3m));
        // fmt::println("0.5*(v3m*v4m - v4m*v3m) = {}", 0.5 * (v3m * v4m - v4m * v3m));

        CHECK(std::abs(nrm_sq(v1) - 13) < eps);
        CHECK(std::abs(bulk_nrm_sq(v1) - 9) < eps);
        CHECK(std::abs(weight_nrm_sq(v1) - 4) < eps);

        CHECK(std::abs(nrm_sq(v1n) - 1) < eps);
        CHECK(std::abs(weight_nrm_sq(v1u) - 1) < eps);

        CHECK(std::abs(bulk_nrm_sq(v3) - 56) < eps);
        CHECK(std::abs(nrm(v4 * v3) - 1.0) < eps);
        CHECK(std::abs(dot(v4, v3) - 1.0) < eps);
        CHECK(std::abs(nrm(wdg(v4, v3))) < eps);
    }

    TEST_CASE("Vec3dp: operations - bulk_dual, weight_dual")
    {
        fmt::println("Vec3dp: operations - bulk_dual, weight_dual");

        auto s = scalar3dp{5.0};
        auto v = vec3dp{1.0, 2.0, 1.0, 1.0};
        auto B = bivec3dp{-1.0, 2.0, 1.0, -10.0, 20.0, 10.0};
        auto t = trivec3dp{3.0, 6.0, 3.0, 3.0};
        auto ps = pscalar3dp{-5.0};

        auto v2 = vec3dp{3.0, -2.0, 2.0, 5.0};
        auto B2 = bivec3dp{1.0, -2.0, 3.0, -10.0, -20.0, 30.0};

        CHECK(bulk_dual(s) == pscalar3dp{s});
        CHECK(weight_dual(s) == pscalar3dp{0.0});

        CHECK(bulk_dual(v) == trivec3dp{1.0, 2.0, 1.0, 0.0});
        CHECK(weight_dual(v) == trivec3dp{0.0, 0.0, 0.0, 1.0});

        CHECK(bulk_dual(B) == bivec3dp{10.0, -20.0, -10.0, 0.0, 0.0, 0.0});
        CHECK(weight_dual(B) == bivec3dp{0.0, 0.0, 0.0, 1.0, -2.0, -1.0});

        CHECK(bulk_dual(t) == vec3dp{0.0, 0.0, 0.0, -3.0});
        CHECK(weight_dual(t) == vec3dp{-3.0, -6.0, -3.0, 0.0});

        CHECK(bulk_dual(ps) == scalar3dp{0.0});
        CHECK(weight_dual(ps) == scalar3dp{ps});

        // duality of wdg and rwdg based on complements
        CHECK(bulk_dual(wdg(v, v2)) == rwdg(bulk_dual(v), bulk_dual(v2)));
        CHECK(bulk_dual(wdg(v, B)) == rwdg(bulk_dual(v), bulk_dual(B)));
        CHECK(bulk_dual(wdg(v, t)) == rwdg(bulk_dual(v), bulk_dual(t)));
        CHECK(bulk_dual(wdg(t, v)) == rwdg(bulk_dual(t), bulk_dual(v)));
        CHECK(bulk_dual(wdg(B, v)) == rwdg(bulk_dual(B), bulk_dual(v)));
        CHECK(bulk_dual(wdg(v, B)) == rwdg(bulk_dual(v), bulk_dual(B)));
        CHECK(bulk_dual(wdg(B, B2)) == rwdg(bulk_dual(B), bulk_dual(B2)));

        // contractions
        CHECK(bulk_contraction(v, v2) == rwdg(v, bulk_dual(v2)));
        CHECK(bulk_contraction(B, B2) == rwdg(B, bulk_dual(B2)));
        CHECK(bulk_contraction(B, v) == rwdg(B, bulk_dual(v))); // contracts v onto B

        CHECK(weight_contraction(v, v2) == rwdg(v, weight_dual(v2)));
        CHECK(weight_contraction(B, B2) == rwdg(B, weight_dual(B2)));
        CHECK(weight_contraction(B, v) == rwdg(B, weight_dual(v)));

        // expansions
        CHECK(weight_expansion(v, v2) == wdg(v, weight_dual(v2)));
        CHECK(weight_expansion(B, B2) == wdg(B, weight_dual(B2)));
        CHECK(weight_expansion(v, B) == wdg(v, weight_dual(B)));

        CHECK(bulk_expansion(v, v2) == wdg(v, bulk_dual(v2)));
        CHECK(bulk_expansion(B, B2) == wdg(B, bulk_dual(B2)));
        CHECK(bulk_expansion(v, B) == wdg(v, bulk_dual(B)));
    }

    ////
    // angle I + II
    ////

    TEST_CASE("Vec3dp: operations - wedge")
    {
        fmt::println("Vec3dp: operations - wedge");

        vec3dp v1{1.0, 0.0, 0.0, 0.0};
        vec3dp v2{normalize(vec3dp(1.0, 1.0, 0.0, 0.0))};
        vec3dp v3{0.0, 1.0, 0.0, 0.0};
        vec3dp v4{normalize(vec3dp(-1.0, 1.0, 0.0, 0.0))};
        vec3dp v5{-1.0, 0.0, 0.0, 0.0};
        vec3dp v6{normalize(vec3dp(-1.0, -1.0, 0.0, 0.0))};
        vec3dp v7{0.0, -1.0, 0.0, 0.0};
        vec3dp v8{normalize(vec3dp(1.0, -1.0, 0.0, 0.0))};

        double sd = 2.3;
        double st = -5.1;
        auto s = scalar3dp{sd};
        auto t = scalar3dp{st};

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

        CHECK(wdg(v1, v1) == bivec3dp{});                     // wdg=0 for collin. vectors
        CHECK(wdg(v1, v2) == -wdg(v2, v1));                   // anticommutative for vect.
        CHECK(wdg(wdg(v1, v2), v3) == wdg(v1, wdg(v2, v3)));  // wdg is associative
        CHECK(wdg(v1, v2 + v3) == wdg(v1, v2) + wdg(v1, v3)); // wdg distributes over add.
        CHECK(wdg(v1 + v2, v3) == wdg(v1, v3) + wdg(v2, v3)); // wdg distributes over add.
        CHECK(wdg(sd * v1, v2) == wdg(v1, sd * v2)); // scalars can be factored out of wdg
        CHECK(wdg(sd * v1, v2) == sd * wdg(v1, v2)); // scalars can be factored out of wdg
        CHECK(wdg(s, t) == wdg(t, s));   // wdg between scalars equivalent to scalar mult.
        CHECK(wdg(s, v1) == wdg(v1, s)); // wdg between scalar and vector
        CHECK(wdg(s, v1) == sd * v1);    // wdg between scalar and vector

        CHECK(std::abs(bulk_nrm(wdg(v1, v1)) - sin(angle(v1, v1))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v2)) - sin(angle(v1, v2))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v3)) - sin(angle(v1, v3))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v4)) - sin(angle(v1, v4))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v5)) - sin(angle(v1, v5))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v6)) - sin(angle(v1, v6))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v7)) - sin(angle(v1, v7))) < eps);
        CHECK(std::abs(bulk_nrm(wdg(v1, v8)) - sin(angle(v1, v8))) < eps);
    }

    TEST_CASE("Vec3dp: operations - project / reject / reflect (vector - vector)")
    {
        fmt::println("Vec3dp: operations - project / reject / reflect (vector - vector)");

        vec3dp v1{5.0, 1.0, 0.0, 1.0};
        vec3dp v2{2.0, 2.0, 0.0, 1.0};

        vec3dp v3{project_onto(v1, v2)};
        vec3dp v4{reject_from(v1, v2)};
        vec3dp v5{v3 + v4};

        // fmt::println("sqrt(26)={}", std::sqrt(26.0));
        // fmt::println("v1  = {: .4f}, bulk_nrm(v1)  = {: .4f}, "
        //              "nrm(v1) = {: .4f}",
        //              v1, bulk_nrm(v1), nrm(v1));
        // fmt::println("v2  = {: .4f}, bulk_nrm(v2)  = {: .4f}, "
        //              "nrm(v2) = {: .4f}",
        //              v2, bulk_nrm(v2), nrm(v2));
        // fmt::println("");
        // fmt::println("v3 = project_onto(v1, v2) = {: .4f}, bulk_nrm(v3) = {: .4f}, "
        //              "nrm(v3) = {: .4f}",
        //              v3, bulk_nrm(v3), nrm(v3));
        // fmt::println("v4 = reject_from(v1, v2)  = {: .4f}, bulk_nrm(v4) = {: .4f}, "
        //              "nrm(v4) = {: .4f}",
        //              v4, bulk_nrm(v4), nrm(v4));
        // fmt::println("v5 = v3 + v4              = {: .4f}, bulk_nrm(v5) = {: .4f}, "
        //              "nrm(v5) = {: .4f}",
        //              v5, bulk_nrm(v5), nrm(v5));
        // fmt::println("");

        // // this helps to understand, why the dot-product is sufficient
        // auto w = wdg(v1, v2);
        // auto i = inv(v2);
        // fmt::println("wdg(v1,v2)         = {: .4f}", w);
        // fmt::println("inv(v2)            = {: .4f}", i);
        // fmt::println("wdg(v1,v2)*inv(v2) = {: .4f}", w * i);
        // fmt::println("");

        // vec3dp v{4.0, 1.0, 0.0, 0.0};
        // vec3dp b{e1_3dp};
        // auto B = bivec3dp{e23_3dp};

        // auto UB = bivec3dp{e23_3dp + e12_3dp};
        // fmt::println("v   = {}", v);
        // fmt::println("b   = {}", b);
        // fmt::println("B   = {}", B);
        // fmt::println("UB  = {}", UB);
        // fmt::println("reflect_on_vec(v,b)  = {}", reflect_on_vec(v, b));
        // fmt::println("reflect_on_hyp(v,e3_3dp)  = {}", reflect_on_hyp(v, e3_3dp));
        // fmt::println("reflect_on(v,B)  = {}", reflect_on(v, B));
        // fmt::println("reflect_on(UB,B) = {}", reflect_on(UB, B));
        // fmt::println("");

        CHECK(v3 + v4 == v5);
        CHECK(v5 == v1);

        // // just to suppress unused variable warnings
        // CHECK(b == e1_3dp);

        // CHECK(reflect_on_hyp(v, e2_3dp) == reflect_on(v, B));

        // checking time required
        //
        // auto start = std::chrono::system_clock::now();
        // for (size_t i = 0; i < 10000000; ++i) {
        //     vec3dp v = reject_from(v1, v2);
        // }
        // auto end = std::chrono::system_clock::now();
        // auto elapsed =
        // std::chrono::duration_cast<std::chrono::milliseconds>(end
        // -
        // start); fmt::println("The measurement took {}", elapsed);
    }

    ////
    // Vec3dp: project / reject / reflect etc.
    // MVec3dp: ...
    ////

    TEST_CASE("MVec3dp: complement operation")
    {
        fmt::println("MVec3dp: complement operation");

        auto s = scalar3dp{5.0};
        auto v = vec3dp{1.0, 2.0, 3.0, 1.0};
        auto B = bivec3dp{-1.0, 2.0, 1.0, -10.0, 20.0, 10.0};
        auto t = trivec3dp{3.0, 6.0, 9.0, 3.0};
        auto ps = pscalar3dp{-5.0};

        auto mv = mvec3dp{s, v, B, t, ps};
        auto mv_e = mvec3dp_e{s, B, ps};
        auto mv_u = mvec3dp{v, t};

        auto s2 = scalar3dp{2.5};
        auto v2 = vec3dp{0.5, 1.0, 1.5, 2.0};
        auto B2 = bivec3dp{5.0, 10.0, 15.0, 50.0, 100.0, 150.0};
        auto t2 = trivec3dp{1.5, 3.0, 4.5, 6.0};
        auto ps2 = pscalar3dp{-1.5};

        CHECK(wdg(s, r_cmpl(s)) == nrm_sq(s) * I_3dp);
        CHECK(wdg(v, r_cmpl(v)) == nrm_sq(v) * I_3dp);
        CHECK(wdg(B, r_cmpl(B)) == nrm_sq(B) * I_3dp);
        CHECK(wdg(t, r_cmpl(t)) == nrm_sq(t) * I_3dp);
        CHECK(wdg(ps, r_cmpl(ps)) == nrm_sq(ps) * I_3dp);

        CHECK(wdg(l_cmpl(s), s) == nrm_sq(s) * I_3dp);
        CHECK(wdg(l_cmpl(v), v) == nrm_sq(v) * I_3dp);
        CHECK(wdg(l_cmpl(B), B) == nrm_sq(B) * I_3dp);
        CHECK(wdg(l_cmpl(t), t) == nrm_sq(t) * I_3dp);
        CHECK(wdg(l_cmpl(ps), ps) == nrm_sq(ps) * I_3dp);

        // complement properties (right and left complements are inverse to each other)
        CHECK(l_cmpl(r_cmpl(s)) == s);
        CHECK(l_cmpl(r_cmpl(v)) == v);
        CHECK(l_cmpl(r_cmpl(B)) == B);
        CHECK(l_cmpl(r_cmpl(t)) == t);
        CHECK(l_cmpl(r_cmpl(ps)) == ps);
        CHECK(l_cmpl(r_cmpl(mv)) == mv);
        CHECK(l_cmpl(r_cmpl(mv_e)) == mv_e);
        CHECK(l_cmpl(r_cmpl(mv_u)) == mv_u);

        CHECK(r_cmpl(l_cmpl(s)) == s);
        CHECK(r_cmpl(l_cmpl(v)) == v);
        CHECK(r_cmpl(l_cmpl(B)) == B);
        CHECK(r_cmpl(l_cmpl(t)) == t);
        CHECK(r_cmpl(l_cmpl(ps)) == ps);
        CHECK(r_cmpl(l_cmpl(mv)) == mv);
        CHECK(r_cmpl(l_cmpl(mv_e)) == mv_e);
        CHECK(r_cmpl(l_cmpl(mv_u)) == mv_u);

        // complement values
        CHECK(r_cmpl(scalar3dp(1.0)) == I_3dp);
        CHECK(r_cmpl(e1_3dp) == e423_3dp);
        CHECK(r_cmpl(e2_3dp) == e431_3dp);
        CHECK(r_cmpl(e3_3dp) == e412_3dp);
        CHECK(r_cmpl(e4_3dp) == e321_3dp);
        CHECK(r_cmpl(e41_3dp) == -e23_3dp);
        CHECK(r_cmpl(e42_3dp) == -e31_3dp);
        CHECK(r_cmpl(e43_3dp) == -e12_3dp);
        CHECK(r_cmpl(e23_3dp) == -e41_3dp);
        CHECK(r_cmpl(e31_3dp) == -e42_3dp);
        CHECK(r_cmpl(e12_3dp) == -e43_3dp);
        CHECK(r_cmpl(e423_3dp) == -e1_3dp);
        CHECK(r_cmpl(e431_3dp) == -e2_3dp);
        CHECK(r_cmpl(e412_3dp) == -e3_3dp);
        CHECK(r_cmpl(e321_3dp) == -e4_3dp);
        CHECK(r_cmpl(I_3dp) == scalar2d(1.0));

        CHECK(l_cmpl(scalar3dp(1.0)) == I_3dp);
        CHECK(l_cmpl(e1_3dp) == -e423_3dp);
        CHECK(l_cmpl(e2_3dp) == -e431_3dp);
        CHECK(l_cmpl(e3_3dp) == -e412_3dp);
        CHECK(l_cmpl(e4_3dp) == -e321_3dp);
        CHECK(l_cmpl(e41_3dp) == -e23_3dp);
        CHECK(l_cmpl(e42_3dp) == -e31_3dp);
        CHECK(l_cmpl(e43_3dp) == -e12_3dp);
        CHECK(l_cmpl(e23_3dp) == -e41_3dp);
        CHECK(l_cmpl(e31_3dp) == -e42_3dp);
        CHECK(l_cmpl(e12_3dp) == -e43_3dp);
        CHECK(l_cmpl(e423_3dp) == e1_3dp);
        CHECK(l_cmpl(e431_3dp) == e2_3dp);
        CHECK(l_cmpl(e412_3dp) == e3_3dp);
        CHECK(l_cmpl(e321_3dp) == e4_3dp);
        CHECK(l_cmpl(I_3dp) == scalar2d(1.0));

        // linearity of the right and left complement operation
        double a = 2.0;
        CHECK(r_cmpl(a * s) == a * r_cmpl(s));
        CHECK(r_cmpl(a * v) == a * r_cmpl(v));
        CHECK(r_cmpl(a * B) == a * r_cmpl(B));
        CHECK(r_cmpl(a * t) == a * r_cmpl(t));
        CHECK(r_cmpl(a * ps) == a * r_cmpl(ps));

        CHECK(r_cmpl(s + s2) == r_cmpl(s) + r_cmpl(s2));
        CHECK(r_cmpl(v + v2) == r_cmpl(v) + r_cmpl(v2));
        CHECK(r_cmpl(B + B2) == r_cmpl(B) + r_cmpl(B2));
        CHECK(r_cmpl(t + t2) == r_cmpl(t) + r_cmpl(t2));
        CHECK(r_cmpl(ps + ps2) == r_cmpl(ps) + r_cmpl(ps2));

        CHECK(l_cmpl(a * s) == a * l_cmpl(s));
        CHECK(l_cmpl(a * v) == a * l_cmpl(v));
        CHECK(l_cmpl(a * B) == a * l_cmpl(B));
        CHECK(l_cmpl(a * t) == a * l_cmpl(t));
        CHECK(l_cmpl(a * ps) == a * l_cmpl(ps));

        CHECK(l_cmpl(s + s2) == l_cmpl(s) + l_cmpl(s2));
        CHECK(l_cmpl(v + v2) == l_cmpl(v) + l_cmpl(v2));
        CHECK(l_cmpl(B + B2) == l_cmpl(B) + l_cmpl(B2));
        CHECK(l_cmpl(t + t2) == l_cmpl(t) + l_cmpl(t2));
        CHECK(l_cmpl(ps + ps2) == l_cmpl(ps) + l_cmpl(ps2));

        // check regressive wedge product expressed in terms of
        // wegde product and complement operation
        CHECK(rwdg(t, t2) == l_cmpl(wdg(r_cmpl(t), r_cmpl(t2))));
        CHECK(rwdg(t, B) == l_cmpl(wdg(r_cmpl(t), r_cmpl(B))));
        CHECK(rwdg(B, t) == l_cmpl(wdg(r_cmpl(B), r_cmpl(t))));
        CHECK(rwdg(t, v) == l_cmpl(wdg(r_cmpl(t), r_cmpl(v))));
        CHECK(rwdg(v, t) == l_cmpl(wdg(r_cmpl(v), r_cmpl(t))));
        CHECK(rwdg(B, B2) == l_cmpl(wdg(r_cmpl(B), r_cmpl(B2))));

        CHECK(l_cmpl(rwdg(B, B2)) == wdg(r_cmpl(B), r_cmpl(B2)));
        CHECK(l_cmpl(wdg(v, v2)) == rwdg(r_cmpl(v), r_cmpl(v2)));
    }

    TEST_CASE("3dp: objects - bulk and weight")
    {
        fmt::println("3dp: objects - bulk and weight");

        auto p1 = vec3dp{-2.0, 1.0, 1.0, 1.0};
        auto p2 = vec3dp{2.0, 1.0, 2.0, 2.0};
        auto p3 = vec3dp{2.0, 1.0, 1.0, 3.0};

        auto l1 = bivec3dp{2.0, 1.0, -1.0, 20.0, 10.0, -10.0};
        auto l2 = bivec3dp{-2.0, -1.0, 3.0, -20.0, -10.0, 30.0};
        auto l3 = wdg(p1, p3);

        auto t1 = trivec3dp{-6.0, 3.0, 3.0, 3.0};
        auto t2 = trivec3dp{6.0, 3.0, 6.0, 6.0};
        auto t3 = trivec3dp{6.0, 3.0, 3.0, 9.0};

        // bulk and weight
        CHECK(bulk(p1) == vec3dp{p1.x, p1.y, p1.z, 0.0});
        CHECK(weight(p1) == vec3dp{0.0, 0.0, 0.0, p1.w});
        CHECK(bulk(l1) == bivec3dp{0.0, 0.0, 0.0, l1.mx, l1.my, l1.mz});
        CHECK(weight(l1) == bivec3dp{l1.vx, l1.vy, l1.vz, 0.0, 0.0, 0.0});
        CHECK(bulk(t1) == trivec3dp{0.0, 0.0, 0.0, t1.w});
        CHECK(weight(t1) == trivec3dp{t1.x, t1.y, t1.z, 0.0});

        // bulk_nrm_sq and weight_nrm_sq
        CHECK(bulk_nrm_sq(p1) == p1.x * p1.x + p1.y * p1.y + p1.z * p1.z);
        CHECK(weight_nrm_sq(p1) == p1.w * p1.w);
        CHECK(bulk_nrm_sq(l1) == l1.mx * l1.mx + l1.my * l1.my + l1.mz * l1.mz);
        CHECK(weight_nrm_sq(l1) == l1.vx * l1.vx + l1.vy * l1.vy + l1.vz * l1.vz);
        CHECK(bulk_nrm_sq(t1) == t1.w * t1.w);
        CHECK(weight_nrm_sq(t1) == t1.x * t1.x + t1.y * t1.y + t1.z * t1.z);

        // bulk_nrm and weight_nrm
        CHECK(bulk_nrm(p1) == std::sqrt(p1.x * p1.x + p1.y * p1.y + p1.z * p1.z));
        CHECK(weight_nrm(p1) == std::sqrt(p1.w * p1.w));
        CHECK(bulk_nrm(l1) == std::sqrt(l1.mx * l1.mx + l1.my * l1.my + l1.mz * l1.mz));
        CHECK(weight_nrm(l1) == std::sqrt(l1.vx * l1.vx + l1.vy * l1.vy + l1.vz * l1.vz));
        CHECK(bulk_nrm(t1) == std::sqrt(t1.w * t1.w));
        CHECK(weight_nrm(t1) == std::sqrt(t1.x * t1.x + t1.y * t1.y + t1.z * t1.z));


        // geom_nrm
        CHECK(geom_nrm(p1).c0 == bulk_nrm(p1));
        CHECK(geom_nrm(p1).c1 == weight_nrm(p1));
        CHECK(geom_nrm(l1).c0 == bulk_nrm(l1));
        CHECK(geom_nrm(l1).c1 == weight_nrm(l1));
        CHECK(geom_nrm(t1).c0 == bulk_nrm(t1));
        CHECK(geom_nrm(t1).c1 == weight_nrm(t1));
        CHECK(geom_nrm(l1).c0 / geom_nrm(l1).c1 == geom_nrm(l2).c0 / geom_nrm(l2).c1);

        // attitude
        CHECK(att(p1) == p1.w);
        CHECK(att(p1) == rwdg(p1, r_cmpl(e4_3dp)));
        CHECK(att(p2) == p2.w);
        CHECK(att(p2) == rwdg(p2, r_cmpl(e4_3dp)));
        CHECK(att(l1) == vec3dp{l1.vx, l1.vy, l1.vz, 0.0});
        CHECK(att(l1) == rwdg(l1, r_cmpl(e4_3dp)));
        CHECK(att(l2) == rwdg(l2, r_cmpl(e4_3dp)));
        CHECK(att(l3) == rwdg(l3, r_cmpl(e4_3dp)));
        CHECK(att(t1) == bivec3dp{0.0, 0.0, 0.0, t1.x, t1.y, t1.z});
        CHECK(att(t1) == rwdg(t1, r_cmpl(e4_3dp)));
        CHECK(att(t2) == bivec3dp{0.0, 0.0, 0.0, t2.x, t2.y, t2.z});
        CHECK(att(t2) == rwdg(t2, r_cmpl(e4_3dp)));
    }

    TEST_CASE("3dp: objects - euclidean distance")
    {
        fmt::println("3dp: objects - euclidean distance");

        auto p1 = vec3dp{0.0, 1.0, 0.0, 1.0};
        auto p2 = vec3dp{0.0, 3.0, 0.0, -2.0};

        auto p1u = unitize(p1);
        auto p2u = unitize(p2);
        auto delta = p2u - p1u;

        // fmt::println("");
        // fmt::println("p1  = {}, p2  = {}", p1, p2);
        // fmt::println("p1u = {}, p2u = {}", p1u, p2u);
        // fmt::println("delta = {}, nrm(delta) = {}", delta, nrm(delta));
        // fmt::println("");

        auto l1 = bivec3dp{2.0, 1.0, -1.0, 20.0, 10.0, -10.0};
        auto l2 = bivec3dp{-2.0, -1.0, 3.0, -20.0, -10.0, 30.0};

        auto t1 = trivec3dp{-6.0, 3.0, 3.0, 3.0};
        auto t2 = trivec3dp{6.0, 3.0, 6.0, 6.0};
        auto t3 = trivec3dp{6.0, 3.0, 3.0, 9.0};


        auto dist_p1_p2 = euclidean_distance3dp(p1, p2);
        auto dist_p1_l1 = euclidean_distance3dp(p1, l1);
        auto dist_l1_p1 = euclidean_distance3dp(l1, p1);
        auto dist_p1_t1 = euclidean_distance3dp(p1, t1);
        auto dist_l1_l2 = euclidean_distance3dp(l1, l2);
        fmt::println("");
        fmt::println("dist_p1_p2 = {}, dist = {}", dist_p1_p2,
                     dist_p1_p2.c0 / dist_p1_p2.c1);
        fmt::println("dist_p1_l1 = {}, dist = {}", dist_p1_l1,
                     dist_p1_l1.c0 / dist_p1_l1.c1);
        fmt::println("dist_l1_p1 = {}, dist = {}", dist_l1_p1,
                     dist_l1_p1.c0 / dist_l1_p1.c1);
        fmt::println("");
        fmt::println("dist_p1_t1 = {}, dist = {}", dist_p1_t1,
                     dist_p1_t1.c0 / dist_p1_t1.c1);
        fmt::println("");
        fmt::println("dist_l1_l2 = {}, dist = {}", dist_l1_l2,
                     dist_l1_l2.c0 / dist_l1_l2.c1);
        fmt::println("");

        CHECK(nrm(delta) == dist_p1_p2.c0 / dist_p1_p2.c1);
        CHECK(dist_p1_l1 == dist_l1_p1);

        // auto res = rwdg(pscalar3dp(2.5), horizon_3dp);
        // fmt::println("res = {}", res);
    }

    TEST_CASE("pga_3dp<3,0,1> - product tests")
    {
        fmt::println("pga_3dp<3,0,1> - product tests");
        // tests based on "The inner products of geometric algebra", Leo Dorst


        // auto s1 = scalar2dp{2.0};
        // auto v1 = vec2dp{1.0, -3.0, 0.0};
        // // auto v1 = vec2dp{1.0, 0.0, 0.0};
        // auto b1 = bivec2dp{2.0, -4.0, 1.0};
        // auto ps1 = pscalar2dp{-2.0};

        // auto s2 = scalar2dp{-1.0};
        // auto v2 = vec2dp{2.0, 1.5, 0.0};
        // // auto v2 = vec2dp{0.0, 1.0, 0.0};
        // auto b2 = bivec2dp{-2.0, 10.5, 3.0};
        // auto ps2 = pscalar2dp{3.0};

        // auto s3 = scalar2dp{-10.0};
        // auto v3 = vec2dp{-2.0, 4.0, -10.0};
        // auto b3 = bivec2dp{1.0, -2.0, 3.0};
        // auto ps3 = pscalar2dp{-12.0};

        // dot product (difference due to different definition of scalar product)
        CHECK(dot(scalar2dp(1.0), scalar2dp(1.0)) == scalar2dp(1.0));
        CHECK(dot(e1_3dp, e1_3dp) == scalar3dp(1.0));
        CHECK(dot(e2_3dp, e2_3dp) == scalar3dp(1.0));
        CHECK(dot(e3_3dp, e3_3dp) == scalar3dp(1.0));
        CHECK(dot(e4_3dp, e4_3dp) == scalar3dp(0.0));

        CHECK(dot(e41_3dp, e41_3dp) == scalar3dp(0.0));
        CHECK(dot(e42_3dp, e42_3dp) == scalar3dp(0.0));
        CHECK(dot(e43_3dp, e43_3dp) == scalar3dp(0.0));
        CHECK(dot(e23_3dp, e23_3dp) == scalar3dp(-1.0)); // differs from Lengyel!
        CHECK(dot(e31_3dp, e31_3dp) == scalar3dp(-1.0)); // differs from Lengyel!
        CHECK(dot(e12_3dp, e12_3dp) == scalar3dp(-1.0)); // differs from Lengyel!

        CHECK(dot(e423_3dp, e423_3dp) == scalar3dp(0.0));
        CHECK(dot(e431_3dp, e431_3dp) == scalar3dp(0.0));
        CHECK(dot(e412_3dp, e412_3dp) == scalar3dp(0.0));
        CHECK(dot(e321_3dp, e321_3dp) == scalar3dp(-1.0)); // differs from Lengyel!
        CHECK(dot(pscalar3dp(1.0), pscalar3dp(1.0)) == scalar3dp(0.0));

        // regressive dot product
        CHECK(rdot(scalar2dp(1.0), scalar2dp(1.0)) == pscalar2dp(0.0));
        CHECK(rdot(e1_3dp, e1_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e2_3dp, e2_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e3_3dp, e3_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e4_3dp, e4_3dp) == pscalar3dp(-1.0));   // differs from Lengyel!
        CHECK(rdot(e41_3dp, e41_3dp) == pscalar3dp(-1.0)); // differs from Lengyel!
        CHECK(rdot(e42_3dp, e42_3dp) == pscalar3dp(-1.0)); // differs from Lengyel!
        CHECK(rdot(e43_3dp, e43_3dp) == pscalar3dp(-1.0)); // differs from Lengyel!
        CHECK(rdot(e23_3dp, e23_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e31_3dp, e31_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e12_3dp, e12_3dp) == pscalar3dp(0.0));
        CHECK(rdot(e423_3dp, e423_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e431_3dp, e431_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e412_3dp, e412_3dp) == pscalar3dp(1.0));
        CHECK(rdot(e321_3dp, e321_3dp) == pscalar3dp(0.0));
        CHECK(rdot(pscalar3dp(1.0), pscalar3dp(1.0)) == pscalar3dp(1.0));


        // // wedge product
        // // check full permissible range of arguments, even the ones delivering 0)
        // CHECK(wdg(b1, b2) == scalar2dp(0.0));
        // CHECK(wdg(v1, ps1) == scalar2dp(0.0));
        // CHECK(wdg(ps1, v1) == scalar2dp(0.0));
        // CHECK(wdg(b1, ps1) == scalar2dp(0.0));
        // CHECK(wdg(ps1, b1) == scalar2dp(0.0));
        // CHECK(wdg(ps1, ps1) == scalar2dp(0.0));

        // // fmt::println("   b1           = {}", b1);
        // // fmt::println("   b2           = {}", b2);
        // // fmt::println("   b1*b2        = {}", b1 * b2);
        // // fmt::println("   cmt(b1, b2)  = {}", cmt(b1, b2));
        // // fmt::println("   b1<<b2       = {}", b1 << b2);
        // CHECK(wdg(b1, b2) == nrm_sq(b1 * b2 - cmt(b1, b2) - (b1 << b2)));
        // CHECK(wdg(v1, ps1) == nrm_sq(0.5 * (v1 * ps1 + rev(ps1) * v1)));
        // CHECK(wdg(ps1, v1) == nrm_sq(0.5 * (ps1 * v1 + v1 * rev(ps1))));


        // CHECK((v1 << wdg(v2, v3)) == wdg(v1 << v2, v3) + wdg(gr_inv(v2), v1 <<
        // v3)); CHECK((wdg(v1, v2) << b1) == (v1 << (v2 << b1)));

        // // contractions - check full permissible range of arguments, even the ones
        // // delivering 0 as a result
        // CHECK((ps1 << s1) == 0.0);
        // CHECK((ps1 << v1) == 0.0);
        // CHECK((ps1 << b1) == 0.0);
        // CHECK((b1 << s1) == 0.0);
        // CHECK((b1 << v1) == 0.0);
        // CHECK((v1 << s1) == 0.0);
        // //
        // CHECK((s1 >> ps1) == 0.0);
        // CHECK((v1 >> ps1) == 0.0);
        // CHECK((b1 >> ps1) == 0.0);
        // CHECK((s1 >> b1) == 0.0);
        // CHECK((v1 >> b1) == 0.0);
        // CHECK((s1 >> v1) == 0.0);

        // // 2.3.2

        // // 3rd with vectors and bivectors directly
        // CHECK(v1 * b1 == (v1 << b1) + wdg(v1, b1));
        // CHECK(b1 * v1 == (b1 >> v1) + wdg(b1, v1));
        // CHECK((b1 >> v1) == -(v1 << gr_inv(b1)));

        // // fmt::println("   v1 << b1        = {}", v1 << b1);
        // // fmt::println("   v1 * b1         = {}", v1 * b1);
        // // fmt::println("   gr_inv(b1) * v1 = {}", gr_inv(b1) * v1);
        // CHECK((v1 << b1) == gr1(0.5 * (v1 * b1 - gr_inv(b1) * v1)));
        // CHECK((b1 >> v1) == gr1(0.5 * (b1 * v1 - v1 * gr_inv(b1))));

        // // fmt::println("   wdg(v1, b1)     = {}", wdg(v1, b1));
        // // fmt::println("   v1 * b1         = {}", v1 * b1);
        // // fmt::println("   gr_inv(b1) * v1 = {}", gr_inv(b1) * v1);
        // CHECK(wdg(v1, b1) == gr3(0.5 * (v1 * b1 + gr_inv(b1) * v1)));
        // CHECK(wdg(b1, v1) == gr3(0.5 * (b1 * v1 + v1 * gr_inv(b1))));

        // CHECK(wdg(v1, (b1 << b2)) == ((v1 << b1) << b2) + (gr_inv(b1) << wdg(v1,
        // b2)));

        // // cross-check direct implementation of rwdg by comparing with wdg
        // CHECK(rwdg(b1, b2) == cmpl(wdg(cmpl(b1), cmpl(b2))));
    }

    ////

    TEST_CASE(
        "G<3,0,1> - pga3dp simple applications, complements, contraction, expansions")
    {
        fmt::println("G<3,0,1> - pga3dp simple applications, complements, contraction,"
                     "expansions");

        auto s1 = scalar3dp{2.0};
        auto v1 = vec3dp{1.0, -3.0, 0.0, 1.0};
        auto v2 = vec3dp{-1.0, 3.0, 2.0, 1.0};
        auto b1 = bivec3dp{2.0, -4.0, 1.0, -4.0, 8.0, -2.0};
        auto t1 = trivec3dp{1.0, -5.0, 6.0, 7.0};
        auto ps1 = pscalar3dp{-2.0};
        auto M1 = mvec3dp{s1, v1, b1, t1, ps1};

        auto R = rotor(wdg(v1, e1_3dp), deg_to_rad(15));
        CHECK(rotate(v2, R) == rotate_opt(v2, R));
    }

    TEST_CASE("G<3,0,1> - pga3dp join and meet (wdg, rwdg)")
    {
        fmt::println("G<3,0,1> - pga3dp join and meet (wdg, rwdg)");

        auto p1 = point3d{2, 0, 1};
        auto p2 = point3d{4, 3, 1};
        auto l12 = join(p1, p2); // = wdg(p1, p2)

        auto p3 = point3d{0, 6, 1};
        auto pl123 = unitize(join(l12, p3)); // plane parallel to e12

        auto ln1 = join(point3d{0, 3, 0}, point3d{0, 3, 3}).unitize();
        auto pt12 = unitize(meet(ln1, pl123)); // should provide pt(0,3,1)

        auto pl123_2 = unitize(Plane3d(p1, p2, p3)); // check the ctor

        auto pl456 =
            unitize(Plane3d(point3d{2, 0, 3}, point3d{4, 3, 3}, point3d{0, 6, 3}));

        auto l_inf = meet(pl123, pl456); // gives a line at infinity
                                         // weighted by the distance of the planes

        // fmt::println("");
        // fmt::println("p1 = {}, p2 = {}, p3 = {}", p1, p2, p3);
        // fmt::println("l12      = {}, pl123      = {}", l12, pl123);
        // fmt::println("att(l12) = {}, att(pl123) = {}", att(l12), att(pl123));
        // fmt::println("pt12  = meet(l1, pl123)   = {}", pt12);
        // fmt::println("");
        // fmt::println("pl456 = {}", pl456);
        // fmt::println("l_inf = meet(pl123, pl456) = {}", l_inf);
        // fmt::println("att(l_inf)) = {}", att(l_inf));

        CHECK(join(p1, p2) == wdg(p1, p2));
        CHECK(pt12 == vec3dp{0, 3, 1, 1});
        CHECK(pl123 == pl123_2);
        CHECK(att(l_inf) == vec3dp{}); // line is at infinity
    }

} // TEST_SUITE("Projective Geometric Algebra (PGA)")