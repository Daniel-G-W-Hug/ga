"""Quantify the Python-dispatch tax for a tight loop of geometric products.

Companion to gpr_overhead_bench.cpp. Both run N = 50 000 geometric products
`vec3d * vec3d -> mvec3d_e` with deterministic data, so the C++/Python ratio
isolates the cost of "calling out to C++ N times from Python" for a calc
that's small per call (vec3d gpr is ~9 multiplies + 6 adds + 4-double result).

Three Python paths are timed:

  1. ga_py bindings        — one C++ call per element via `a * b`. The dispatch
                             tax is the fee you pay per element.
  2. pure Python (floats)  — hand-coded gpr formula on float tuples. Shows what
                             "no C++ at all" costs in CPython.
  3. numpy vectorised      — same calc, but as broadcast operations on a
                             pre-flattened (N, 3) array. The "right way" to do
                             this in Python once you commit to numpy.

Run from the project root, with the wrapper venv active and `_ga_py` built:

    PYTHONPATH=build/ga_py:ga_py/python \\
        ga_py/.venv/bin/python ga_py/benchmark/gpr_overhead_bench.py

For a fair comparison, run the C++ side in the same build mode:

    build/ga_py/benchmark/gpr_overhead_bench

Captured runs live in ga_py/benchmark/gpr_overhead_results_{debug,release}.txt.
"""
from __future__ import annotations

import sys
import timeit

import numpy as np

import ga_py


N = 50_000
REPEAT = 7


def fmt_ms(t: float) -> str:
    return f"{t * 1e3:7.2f} ms"


def fmt_ns(t: float, n: int) -> str:
    return f"{t * 1e9 / n:6.1f} ns/call"


def make_vecs(offset: float) -> list:
    """Match the C++ side's data: linearly increasing x with given offset.

    Vec3d(i*0.001 + offset, i*0.001 + offset + 0.5, i*0.001 + offset + 1.0)
    for i in [0, N).
    """
    return [ga_py.ega.vec3d(i * 0.001 + offset,
                            i * 0.001 + offset + 0.5,
                            i * 0.001 + offset + 1.0)
            for i in range(N)]


def make_float_tuples(offset: float) -> list:
    return [(i * 0.001 + offset,
             i * 0.001 + offset + 0.5,
             i * 0.001 + offset + 1.0)
            for i in range(N)]


def make_np_array(offset: float) -> np.ndarray:
    i = np.arange(N, dtype=np.float64)
    x = i * 0.001 + offset
    return np.stack([x, x + 0.5, x + 1.0], axis=1)  # (N, 3)


# --- Path 1: ga_py bindings ---------------------------------------------------

def gpr_via_bindings(as_, bs):
    """One C++ call per element via the bound `__mul__`."""
    return [a * b for a, b in zip(as_, bs)]


# --- Path 2: pure Python -----------------------------------------------------

def gpr_pure_python(as_, bs):
    """Hand-coded geometric product on float tuples, no C++ at all.

    For vec3d a, b in EGA3D, a * b = (a · b) + (a ∧ b), giving:
      scalar:        ax*bx + ay*by + az*bz
      bivec yz part: ay*bz - az*by
      bivec zx part: az*bx - ax*bz
      bivec xy part: ax*by - ay*bx
    Returned as a 4-tuple (s, byz, bzx, bxy) — same shape as mvec3d_e's
    field layout.
    """
    out = []
    append = out.append
    for (ax, ay, az), (bx, by, bz) in zip(as_, bs):
        s = ax * bx + ay * by + az * bz
        byz = ay * bz - az * by
        bzx = az * bx - ax * bz
        bxy = ax * by - ay * bx
        append((s, byz, bzx, bxy))
    return out


# --- Path 3: numpy vectorised ------------------------------------------------

def gpr_numpy(a_arr: np.ndarray, b_arr: np.ndarray) -> np.ndarray:
    """Vectorised gpr: (N, 3) × (N, 3) -> (N, 4) with [s, byz, bzx, bxy]."""
    ax, ay, az = a_arr[:, 0], a_arr[:, 1], a_arr[:, 2]
    bx, by, bz = b_arr[:, 0], b_arr[:, 1], b_arr[:, 2]
    return np.stack([
        ax * bx + ay * by + az * bz,
        ay * bz - az * by,
        az * bx - ax * bz,
        ax * by - ay * bx,
    ], axis=1)


def bench(label: str, fn) -> float:
    fn()  # warmup
    t = min(timeit.repeat(fn, number=1, repeat=REPEAT))
    print(f"  {label:24s} {fmt_ms(t)}   {fmt_ns(t, N)}")
    return t


def main() -> None:
    print(f"# Python gpr overhead — N={N} elements, best of {REPEAT}")
    print(f"# Python {sys.version.split()[0]}, NumPy {np.__version__}")
    print(f"# vec3d * vec3d -> mvec3d_e  (geometric product)")
    print()

    as_bound = make_vecs(0.0)
    bs_bound = make_vecs(100.0)
    as_tup = make_float_tuples(0.0)
    bs_tup = make_float_tuples(100.0)
    a_np = make_np_array(0.0)
    b_np = make_np_array(100.0)

    t_bind = bench("ga_py bindings",
                   lambda: gpr_via_bindings(as_bound, bs_bound))
    t_pure = bench("pure Python",
                   lambda: gpr_pure_python(as_tup, bs_tup))
    t_np = bench("numpy vectorised",
                 lambda: gpr_numpy(a_np, b_np))

    print()
    print("  ratios (lower is better):")
    print(f"    bindings / numpy        = {t_bind / t_np:6.1f}×")
    print(f"    pure python / numpy     = {t_pure / t_np:6.1f}×")
    print(f"    bindings / pure python  = {t_bind / t_pure:6.2f}×")


if __name__ == "__main__":
    main()
