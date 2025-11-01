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
#include "ga/ga_pga.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)
using namespace hd::ga::pga; // use specific operations of PGA (Projective GA)


/////////////////////////////////////////////////////////////////////////////////////////
// Euclidean geometric algebra of two-dimensional space: ega2d, pga2dp comparison
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("PGA2DP: comparison tests")
{

    TEST_CASE("pga2dp: ega2d comparison tests")
    {

        {
            fmt::println("pga2dp: ega2d comparison tests - congruency");

            auto x = vec2d{1.0, -7.2};
            auto xv = vec2dp{1.0, -7.2, 0.0};
            auto xp = vec2dp{1.0, -7.2, 1.0};

            CHECK(is_congruent2d(x, 3.2 * x));
            CHECK(is_congruent2d(x, -1.2 * x));
            CHECK(is_congruent2dp(xv, 3.2 * xv));
            CHECK(is_congruent2dp(xv, -1.2 * xv));
            CHECK(is_congruent2dp(xp, 3.2 * xp));
            CHECK(is_congruent2dp(xp, -1.2 * xp));
            fmt::println("");
        }

        {
            fmt::println("pga2dp: ega2d comparison tests - points and lines");

            auto xp = vec2dp{1.0, -7.2, 0.0};
            auto P = vec2dp{1.0, -7.2, 1.0};

            auto xq = vec2dp{2.0, 1.4, 0.0};
            auto Q = vec2dp{2.0, 1.4, 1.0};

            auto CG = 0.5 * (P + Q); // CG = 1/n sum_{1}^{n} (P_i)

            CHECK(O_2dp + xp == P);                  // point = origion + direction
            CHECK(P - Q == xp - xq);                 // point difference = direction
            CHECK(P + Q == 2.0 * O_2dp + (xp + xq)); // weighted point
            CHECK(CG == O_2dp + 0.5 * (xp + xq));    // center of gravity
            fmt::println("");
        }

        {
            fmt::println("pga2dp: ega2d comparison tests - simple elements");

            fmt::println("");

            auto xp = vec2dp{1.0, 0.0, 0.0};
            auto P = vec2dp{1.0, 1.0, 1.0};

            auto xq = vec2dp{2.0, 1.4, 0.0};
            auto Q = vec2dp{2.0, 1.4, 1.0};


            auto M = to_val(bulk_nrm(xp) * bulk_nrm(xq) * sin(angle(xp, xq)) *
                            sign(wdg(xp, xq).z));

            fmt::println("angle(xp, xq)       = {}°", rad2deg(angle(xp, xq)));
            fmt::println("wdg(xp, xq)         = {}", wdg(xp, xq));
            fmt::println("att(wdg(xp, xq))    = {}", att(wdg(xp, xq)));
            fmt::println("bulk(wdg(xp, xq))   = {}", bulk(wdg(xp, xq)));
            fmt::println("bulk(wdg(xp, xq)).z = {}", bulk(wdg(xp, xq)).z);
            // simple bivector: vector^vector -> pure moment, w/o directional aspect
            CHECK(att(wdg(xp, xq)) == vec2dp{}); // no direction
            // pure moment M = bulk(wdg(xp, xq)).z
            CHECK(bulk(wdg(xp, xq)).z ==
                  doctest::Approx(M).epsilon(eps)); // Moment is area of parallelogram

            fmt::println("");
            // bound vector: Point ^ vector (including moment and force)
            fmt::println("wdg(P, xp)         = {}", wdg(P, xp));
            fmt::println("att(wdg(P, xp))    = {}", att(wdg(P, xp)));
            fmt::println("bulk(wdg(P, xp))   = {}", bulk(wdg(P, xp)));

            fmt::println("");
            // Point ^ Point (including moment and force)
            fmt::println("wdg(P, Q)         = {}", wdg(P, Q));
            fmt::println("att(wdg(P, Q))    = {}", att(wdg(P, Q)));
            fmt::println("bulk(wdg(P, Q))   = {}", bulk(wdg(P, Q)));
            fmt::println("");
        }

        {
            fmt::println("pga2dp: ega2d comparison tests - shift of a bound vector");

            fmt::println("");

            auto x = vec2dp{1.0, 0.0, 0.0};
            auto y = vec2dp{2.0, 1.4, 0.0};

            auto P = vec2dp{1.0, 1.0, 1.0};
            auto Q = vec2dp{3.0, 2.4, 1.0};

            auto L = wdg(P, x);
            // P ^ x + y ^ x = (P + y) ^ x = Q ^ x
            auto L_shifted = L + wdg(y, x);

            CHECK(L_shifted == wdg(Q, x));

            fmt::println("x         = {}", x);
            fmt::println("y         = {}", y);
            fmt::println("wdg(y, x) = {}", wdg(y, x)); // simple bivector as carrier
                                                       // of bound vector
            fmt::println("L         = {}", L);
            fmt::println("L_shifted = {}", L_shifted);
            fmt::println("");
        }

        {
            fmt::println("pga2dp: ega2d comparison tests");

            auto u = vec2d{2, 0.5};
            auto a = vec2d{2, 5};

            auto M = wdg(a, u);
            auto d = M * inv(u);

            auto L = u + M;

            fmt::println("u = {}, nrm(u) = {}", u, nrm(u));
            fmt::println("a = {}, nrm(a) = {}", a, nrm(a));
            fmt::println("M = wdg(a,u) = {}", M);
            fmt::println("d = M * inv(u) = {}, nrm(d) = {}", d, nrm(d));
            fmt::println("nrm(u) * nrm(d) = {}", nrm(u) * nrm(d));
            fmt::println("L = {}", L);


            fmt::println("");
            auto up = vec2dp{2, 0.5, 0.0};
            auto ap = vec2dp{2, 5, 1.0};

            auto Lp = wdg(ap, up);
            auto dp = unitize(support2dp(Lp));


            fmt::println("up = {}, bulk_nrm(up) = {}", up, bulk_nrm(up));
            fmt::println("ap = {}, bulk_nrm(up) = {}", ap, bulk_nrm(ap));
            fmt::println("dp = unitize(support2dp(Lp)) = {}, bulk_nrm(dp) = {}", dp,
                         bulk_nrm(dp));
            fmt::println("Lp = {}", Lp);
            fmt::println("att(Lp) = {}", att(Lp));   // direction vector
            fmt::println("bulk(Lp) = {}", bulk(Lp)); // resulting moment w.r.t. origion
            fmt::println("");
        }

        {
            fmt::println("pga2dp: ega2d comparison tests (contractions)");

            auto u = vec2d{1, 0};
            auto v = vec2d{0, 1};
            auto uv = wdg(u, v);

            fmt::println("u  = {}, nrm(u)  = {}", u, nrm(u));
            fmt::println("v  = {}, nrm(v)  = {}", v, nrm(v));
            fmt::println("uv = {}, nrm(uv) = {}", uv, nrm(uv));
            fmt::println("");
            fmt::println("uv >> u = {}", uv >> u);
            fmt::println("uv >> v = {}", uv >> v);
            fmt::println("u << uv = {}", u << uv);
            fmt::println("v << uv = {}", v << uv);
            fmt::println("");

            CHECK(rwdg(left_dual(u), uv) == u << uv); // left contractions
            CHECK(rwdg(left_dual(v), uv) == v << uv);
            CHECK(rwdg(uv, right_dual(u)) == uv >> u); // right contractions
            CHECK(rwdg(uv, right_dual(v)) == uv >> v);

            auto up = vec2dp{1, 0, 0};
            auto vp = vec2dp{0, 1, 0};
            auto uvp = wdg(up, vp);

            fmt::println("up  = {}, nrm(up)  = {}", up, nrm(up));
            fmt::println("vp  = {}, nrm(vp)  = {}", vp, nrm(vp));
            fmt::println("uvp = {}, nrm(uvp) = {}", uvp, nrm(uvp));
            fmt::println("");
            fmt::println("uvp >> up = {}", uvp >> up);
            fmt::println("uvp >> vp = {}", uvp >> vp);
            fmt::println("up << uvp = {}", up << uvp);
            fmt::println("vp << uvp = {}", vp << uvp);
            fmt::println("");

            CHECK(rwdg(bulk_dual(up), uvp) == up << uvp);
            CHECK(rwdg(bulk_dual(vp), uvp) == vp << uvp);
            CHECK(rwdg(bulk_dual(up), uvp) == left_bulk_contract2dp(up, uvp));
            CHECK(rwdg(bulk_dual(vp), uvp) == left_bulk_contract2dp(vp, uvp));
            CHECK(rwdg(uvp, bulk_dual(up)) == uvp >> up);
            CHECK(rwdg(uvp, bulk_dual(vp)) == uvp >> vp);
            CHECK(rwdg(uvp, bulk_dual(up)) == right_bulk_contract2dp(uvp, up));
            CHECK(rwdg(uvp, bulk_dual(vp)) == right_bulk_contract2dp(uvp, vp));
        }
    }

} // EGA2D PGA2DP Comparison Tests
