// Benchmark: sta4ds Lorentz transform of a 4-vector --- direct sandwich transform()
// vs the closed-form matrix transform_opt() (scalar and std::vector batch overloads).
//
// Standalone utility (ga + fmt, no doctest). NOT part of the test run; build and run
// it on demand via the `ga_sta_bench_transform` target. Compiled with -O3/NDEBUG
// regardless of CMAKE_BUILD_TYPE (see ga_test/CMakeLists.txt), since the relative
// cost of the two forms depends strongly on the optimization level.
//
// Three scenarios are timed:
//   BATCH  --- one fixed rotor reused for every vector (matrix can be hoisted/reused)
//   SINGLE --- a fresh rotor per call (matrix rebuilt each time; true one-off cost)
//   VECTOR --- the std::vector transform_opt overload (matrix built once, then applied)
//
// Takeaway (see ga/ga_sta4ds_ops.hpp): for one-off transforms the direct transform()
// is faster; only the batch overload (amortizing the rotor-only matrix over many
// vectors) reliably beats it.

#include "ga/ga_sta.hpp"

#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

using namespace hd::ga;
using namespace hd::ga::sta;

int main()
{
    std::mt19937 rng(12345);
    std::uniform_real_distribution<value_t> dist(-1.0, 1.0);

    // a generic unit rotor = boost(g14, phi) composed with rotor(g12, theta):
    // populates all 8 even-grade coefficients c0..c7
    auto const R = get_boost(g14_4ds, 0.5) * get_rotor(g12_4ds, 0.6);

    size_t const N = 2'000'000;
    std::vector<vec4ds> in;
    in.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        in.emplace_back(vec4ds{dist(rng), dist(rng), dist(rng), dist(rng)});
    }

    // a distinct rotor per vector -> matrix cannot be hoisted/reused (single-shot cost)
    std::vector<mvec4ds_e> rots;
    rots.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        rots.emplace_back(get_boost(g14_4ds, dist(rng)) * get_rotor(g12_4ds, dist(rng)));
    }

    int const reps = 20;

    // scenario BATCH: one fixed rotor R reused for all vectors (matrix reuse possible)
    auto bench_batch = [&](char const* name, auto&& fn) {
        auto const t0 = std::chrono::steady_clock::now();
        value_t acc = 0.0;
        for (int r = 0; r < reps; ++r) {
            for (auto const& v : in) {
                auto const w = fn(v, R);
                acc += w.x + w.y + w.z + w.w;
            }
        }
        auto const t1 = std::chrono::steady_clock::now();
        double const ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
        std::printf("  %-22s %6.3f ns/vec   (acc=%.3f)\n", name, ns / (double(N) * reps),
                    acc);
    };

    // scenario SINGLE: a different rotor for every vector (matrix recomputed each call)
    auto bench_single = [&](char const* name, auto&& fn) {
        auto const t0 = std::chrono::steady_clock::now();
        value_t acc = 0.0;
        for (int r = 0; r < reps; ++r) {
            for (size_t i = 0; i < N; ++i) {
                auto const w = fn(in[i], rots[i]);
                acc += w.x + w.y + w.z + w.w;
            }
        }
        auto const t1 = std::chrono::steady_clock::now();
        double const ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
        std::printf("  %-22s %6.3f ns/vec   (acc=%.3f)\n", name, ns / (double(N) * reps),
                    acc);
    };

    // scenario VECTOR-OVERLOAD: the std::vector transform_opt (matrix built once)
    auto bench_vecov = [&]() {
        auto const t0 = std::chrono::steady_clock::now();
        value_t acc = 0.0;
        for (int r = 0; r < reps; ++r) {
            auto const out = transform_opt(in, R);
            for (auto const& w : out) {
                acc += w.x + w.y + w.z + w.w;
            }
        }
        auto const t1 = std::chrono::steady_clock::now();
        double const ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
        std::printf("  %-22s %6.3f ns/vec   (acc=%.3f) [incl. alloc+sum]\n",
                    "transform_opt[vector]", ns / (double(N) * reps), acc);
    };

    std::printf("sta4ds transform benchmark: N=%zu vectors, reps=%d\n\n", N, reps);
    for (int k = 0; k < 2; ++k) {
        std::printf("BATCH  (one rotor, many vectors):\n");
        bench_batch("transform",
                    [](auto const& v, auto const& Rot) { return transform(v, Rot); });
        bench_batch("transform_opt",
                    [](auto const& v, auto const& Rot) { return transform_opt(v, Rot); });
        bench_vecov();
        std::printf("SINGLE (rotor varies per call):\n");
        bench_single("transform",
                     [](auto const& v, auto const& Rot) { return transform(v, Rot); });
        bench_single("transform_opt", [](auto const& v, auto const& Rot) {
            return transform_opt(v, Rot);
        });
        std::printf("--\n");
    }
    return 0;
}
