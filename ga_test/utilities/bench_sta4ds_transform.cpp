// Benchmark: sta4ds Lorentz transform of a 4-vector --- direct sandwich transform()
// vs the closed-form matrix transform_opt() (scalar and std::vector batch overloads).
//
// Standalone utility (ga + fmt, no doctest). NOT part of the test run; build and run
// it on demand via the `ga_sta_bench_transform` target. Compiled with -O3/NDEBUG
// regardless of CMAKE_BUILD_TYPE (see ga_test/CMakeLists.txt), since the relative
// cost of the two forms depends strongly on the optimization level.
//
// Two scenarios are timed, each reporting ns/vec, the speedup relative to the direct
// transform() baseline, the fastest method, and a one-line verdict:
//   BATCH  --- one fixed rotor reused for every vector (matrix can be hoisted/reused);
//              also times the std::vector transform_opt overload (matrix built once)
//   SINGLE --- a fresh rotor per call (matrix rebuilt each time; true one-off cost)
//
// Takeaway (see ga/ga_sta4ds_ops.hpp): for one-off transforms the direct transform()
// is faster; only the batch overload (amortizing the rotor-only matrix over many
// vectors) reliably beats it.

#include "ga/ga_sta.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using namespace hd::ga;
using namespace hd::ga::sta;

namespace {

struct Result {
    std::string name;
    double ns_per_vec;
};

// Print one scenario block: a header, the per-method rows (ns/vec + speedup vs the
// first/baseline row + a "fastest" marker), and a plain-language verdict line.
void report(char const* title, char const* subtitle, std::vector<Result> const& rows,
            char const* clause)
{
    double const base = rows.front().ns_per_vec; // row 0 is the transform() baseline
    double best = rows.front().ns_per_vec;
    for (auto const& r : rows) best = std::min(best, r.ns_per_vec);

    std::printf("%s - %s\n", title, subtitle);
    std::printf("  %-24s %9s  %8s   %s\n", "method", "ns/vec", "speedup", "note");
    std::string fastest;
    for (auto const& r : rows) {
        bool const is_best = (r.ns_per_vec == best);
        if (is_best) fastest = r.name;
        std::printf("  %-24s %9.3f  %7.2fx   %s\n", r.name.c_str(), r.ns_per_vec,
                    base / r.ns_per_vec, is_best ? "<- fastest" : "");
    }
    // Headline the most informative ratio: how much the closed form beats the
    // baseline when it wins, or how much slower it is when the baseline wins.
    if (best < base) {
        std::printf("  => %s is %.2fx faster than transform(); %s\n\n", fastest.c_str(),
                    base / best, clause);
    }
    else {
        double worst = base;
        for (auto const& r : rows) worst = std::max(worst, r.ns_per_vec);
        std::printf("  => transform() is fastest (transform_opt %.2fx slower); %s\n\n",
                    worst / base, clause);
    }
}

} // namespace

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
    double checksum = 0.0; // accumulated so the timed work cannot be optimized away

    // --- timing kernels: each returns ns per vector ---

    // BATCH: one fixed rotor R reused for all vectors (matrix reuse possible)
    auto time_batch = [&](auto&& fn) -> double {
        auto const t0 = std::chrono::steady_clock::now();
        value_t acc = 0.0;
        for (int r = 0; r < reps; ++r)
            for (auto const& v : in) {
                auto const w = fn(v, R);
                acc += w.x + w.y + w.z + w.w;
            }
        auto const t1 = std::chrono::steady_clock::now();
        checksum += acc;
        return std::chrono::duration<double, std::nano>(t1 - t0).count() /
               (double(N) * reps);
    };

    // SINGLE: a different rotor for every vector (matrix recomputed each call)
    auto time_single = [&](auto&& fn) -> double {
        auto const t0 = std::chrono::steady_clock::now();
        value_t acc = 0.0;
        for (int r = 0; r < reps; ++r)
            for (size_t i = 0; i < N; ++i) {
                auto const w = fn(in[i], rots[i]);
                acc += w.x + w.y + w.z + w.w;
            }
        auto const t1 = std::chrono::steady_clock::now();
        checksum += acc;
        return std::chrono::duration<double, std::nano>(t1 - t0).count() /
               (double(N) * reps);
    };

    // BATCH via the std::vector transform_opt overload (matrix built once per call)
    auto time_vecov = [&]() -> double {
        auto const t0 = std::chrono::steady_clock::now();
        value_t acc = 0.0;
        for (int r = 0; r < reps; ++r) {
            auto const out = transform_opt(in, R);
            for (auto const& w : out) acc += w.x + w.y + w.z + w.w;
        }
        auto const t1 = std::chrono::steady_clock::now();
        checksum += acc;
        return std::chrono::duration<double, std::nano>(t1 - t0).count() /
               (double(N) * reps);
    };

    auto transform_fn = [](auto const& v, auto const& Rot) { return transform(v, Rot); };
    auto opt_fn = [](auto const& v, auto const& Rot) { return transform_opt(v, Rot); };

    // warmup pass (discarded) to settle caches / CPU frequency before measuring
    (void)time_batch(transform_fn);
    (void)time_single(opt_fn);
    (void)time_vecov();

#ifdef NDEBUG
    char const* mode = "-O3 / NDEBUG (optimized)";
#else
    char const* mode = "DEBUG build -- timings NOT meaningful, rebuild optimized";
#endif

    std::printf("sta4ds transform benchmark   (N=%zu vecs x %d reps, %s, double)\n",
                N, reps, mode);
    std::printf("============================================================="
                "==========\n\n");

    report("BATCH ", "one rotor reused for many vectors",
           {{"transform (baseline)", time_batch(transform_fn)},
            {"transform_opt", time_batch(opt_fn)},
            {"transform_opt[vector]", time_vecov()}},
           "closed-form wins when the rotor matrix is built once and reused");
    std::printf("  [vector] figure includes the per-batch allocation + summation "
                "overhead\n\n");

    report("SINGLE", "a fresh rotor per call",
           {{"transform (baseline)", time_single(transform_fn)},
            {"transform_opt", time_single(opt_fn)}},
           "direct sandwich wins one-offs -- the matrix is rebuilt every call and "
           "never amortised");

    std::printf("(checksum %.3f -- ignore; prevents dead-code elimination)\n", checksum);
    return 0;
}
