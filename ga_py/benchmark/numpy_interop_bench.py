"""Benchmark the numpy buffer-protocol surface vs manual attribute access.

Measures the per-element overhead of the Tier 1 §6.3 surface (`__array__`
plus the ndarray ctor) against the pre-Tier-1 idiom (read fields, repack).

Headline question: when is `np.array(v)` actually faster than reading the
fields by hand?

Run from the project root, with the wrapper venv active:

    PYTHONPATH=build/ga_py:ga_py/python \\
        ga_py/.venv/bin/python ga_py/benchmark/numpy_interop_bench.py

Compares "OLD" (manual unpack) against "NEW" (buffer-protocol) idioms
for three representative payload sizes (3, 8, 16 components). Output is
a flat text table — diff-friendly for regression spotting.

Build mode matters: under `-O3 -DNDEBUG` the per-call overhead halves
relative to `-g`. See ga_py/benchmark/results_*.txt for captured runs
on macOS / arm64 and ga_py/README.md §6.3 for the interpretation.
"""
from __future__ import annotations

import sys
import timeit

import numpy as np

import ga_py


N = 10_000
REPEAT = 7


def fmt(t: float) -> str:
    return f"{t * 1e3:7.2f} ms" if t > 1e-3 else f"{t * 1e6:7.2f} us"


def bench(label: str, fn) -> None:
    fn()  # warmup, also fail-fast on broken setup
    t = min(timeit.repeat(fn, number=1, repeat=REPEAT))
    print(f"  {label:42s} {fmt(t)}")


def section(title: str) -> None:
    print()
    print("=" * 64)
    print(title)
    print("=" * 64)


def main() -> None:
    print(f"# numpy buffer-protocol benchmark — N={N} elements, "
          f"best of {REPEAT}")
    print(f"# Python {sys.version.split()[0]}, NumPy {np.__version__}")

    # ----- vec3d: 3 components ------------------------------------------------
    section("vec3d (3 components) — list <-> (N,3) array")
    vs = [ga_py.ega.vec3d(float(i), float(i + 1), float(i + 2))
          for i in range(N)]
    bench("OLD list comp [v.x, v.y, v.z]",
          lambda: np.array([[v.x, v.y, v.z] for v in vs]))
    bench("NEW np.stack([np.array(v) ...])",
          lambda: np.stack([np.array(v) for v in vs]))
    arr = np.stack([np.array(v) for v in vs])
    bench("OLD cls(*row.tolist())",
          lambda: [ga_py.ega.vec3d(*row.tolist()) for row in arr])
    bench("NEW cls(row)",
          lambda: [ga_py.ega.vec3d(row) for row in arr])

    # ----- mvec3d: 8 components ----------------------------------------------
    section("mvec3d (8 components) — list <-> (N,8) array")
    ms = [ga_py.ega.mvec3d(*[float(i + j) for j in range(8)])
          for i in range(N)]
    bench("OLD list comp on .c0..c7",
          lambda: np.array([[m.c0, m.c1, m.c2, m.c3,
                             m.c4, m.c5, m.c6, m.c7]
                            for m in ms]))
    bench("NEW np.stack([np.array(m) ...])",
          lambda: np.stack([np.array(m) for m in ms]))
    arr8 = np.stack([np.array(m) for m in ms])
    bench("OLD cls(*row.tolist())",
          lambda: [ga_py.ega.mvec3d(*row.tolist()) for row in arr8])
    bench("NEW cls(row)",
          lambda: [ga_py.ega.mvec3d(row) for row in arr8])

    # ----- mvec4d: 16 components ---------------------------------------------
    section("mvec4d (16 components) — list <-> (N,16) array")
    m4s = [ga_py.ega.mvec4d(*[float(i + j) for j in range(16)])
           for i in range(N)]
    bench("OLD list comp on .c0..c15",
          lambda: np.array([[m.c0, m.c1, m.c2, m.c3, m.c4, m.c5, m.c6, m.c7,
                             m.c8, m.c9, m.c10, m.c11,
                             m.c12, m.c13, m.c14, m.c15]
                            for m in m4s]))
    bench("NEW np.stack([np.array(m) ...])",
          lambda: np.stack([np.array(m) for m in m4s]))
    arr16 = np.stack([np.array(m) for m in m4s])
    bench("OLD cls(*row.tolist())",
          lambda: [ga_py.ega.mvec4d(*row.tolist()) for row in arr16])
    bench("NEW cls(row)",
          lambda: [ga_py.ega.mvec4d(row) for row in arr16])

    # ----- single-call latency -----------------------------------------------
    section("Single-element overhead (per call) — vec3d")
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    bench("OLD np.array([v.x, v.y, v.z])",
          lambda: np.array([v.x, v.y, v.z]))
    bench("NEW np.array(v)",
          lambda: np.array(v))


if __name__ == "__main__":
    main()
