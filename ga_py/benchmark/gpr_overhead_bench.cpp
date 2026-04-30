// C++ baseline for the Python-dispatch overhead benchmark.
//
// Measures the wall-clock cost of N = 50 000 geometric products
// `vec3d * vec3d -> mvec3d_e` in a tight loop, with no Python in sight.
// The matching Python script (gpr_overhead_bench.py) does the same N
// products via the ga_py bindings, in pure Python, and via numpy — the
// ratio of times answers "what does Python cost me for this calculation?"
//
// Build mode matters: Debug (`-g`) vs Release (`-O3 -DNDEBUG`) flow
// through CMAKE_BUILD_TYPE; pair Debug-vs-Debug and Release-vs-Release for
// a fair comparison with the Python side. The Python extension `_ga_py`
// is built with the same flags out of the same CMake configure.
//
// Anti-DCE: a volatile sink accumulates one component of every product so
// the compiler can't eliminate the loop body under -O3.
//
// Run: build/ga_py/benchmark/gpr_overhead_bench [N] [repeats]
//   defaults: N = 50000, repeats = 7 (best-of-N is reported)

#include "ga/ga_ega.hpp"

#include <fmt/format.h>

#include <chrono>
#include <cstdlib>
#include <vector>

using namespace hd::ga;       // type aliases (vec3d, mvec3d_e, …)
using namespace hd::ga::ega;  // EGA-specific operators (operator*, gpr, …)

namespace {

// Deterministic data, same shape as the Python side. We don't need actual
// randomness — the timing of `gpr` is data-independent (no branches), and
// reproducibility across runs matters more.
std::vector<vec3d> make_vecs(std::size_t N, double offset)
{
    std::vector<vec3d> out;
    out.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        double const x = static_cast<double>(i) * 0.001 + offset;
        out.emplace_back(x, x + 0.5, x + 1.0);
    }
    return out;
}

// One timed pass over the whole input; returns elapsed seconds.
double one_pass(std::vector<vec3d> const& as,
                std::vector<vec3d> const& bs,
                volatile double& sink)
{
    auto const t0 = std::chrono::steady_clock::now();
    double local_sink = 0.0;
    for (std::size_t i = 0; i < as.size(); ++i) {
        // operator* on (vec3d, vec3d) returns mvec3d_e: scalar + bivector.
        // Field layout: c0 = scalar; c1, c2, c3 = bivector components.
        mvec3d_e const c = as[i] * bs[i];
        local_sink += c.c0 + c.c1 + c.c2 + c.c3;
    }
    auto const t1 = std::chrono::steady_clock::now();
    sink = local_sink;  // observable side-effect → loop body cannot be DCE'd
    return std::chrono::duration<double>(t1 - t0).count();
}

}  // namespace


int main(int argc, char** argv)
{
    std::size_t N = (argc > 1) ? std::strtoull(argv[1], nullptr, 10) : 50000;
    int repeats   = (argc > 2) ? std::atoi(argv[2]) : 7;

    auto const as = make_vecs(N, 0.0);
    auto const bs = make_vecs(N, 100.0);

    volatile double sink = 0.0;
    double best = 1e9;
    for (int r = 0; r < repeats; ++r) {
        double const t = one_pass(as, bs, sink);
        if (t < best) best = t;
    }

    double const ns_per_call = best * 1e9 / static_cast<double>(N);

    fmt::print("# C++ gpr overhead — N={} elements, best of {}\n", N, repeats);
    fmt::print("# vec3d * vec3d -> mvec3d_e  (geometric product)\n");
    fmt::print("\n");
    fmt::print("  C++ tight loop:    {:7.2f} ms total   {:6.1f} ns/call\n",
               best * 1e3, ns_per_call);
    fmt::print("  (sink={:.6e} — anti-DCE, ignore the value)\n",
               static_cast<double>(sink));
    return 0;
}
