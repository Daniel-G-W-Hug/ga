# ga_py micro-benchmarks

Reproducible runs that quantify two questions discussed in
[ga_py/README.md §6.3 / §5.5](../README.md):

1. **Is `np.array(v)` actually faster than `[v.x, v.y, v.z]`?**
   (numpy buffer-protocol overhead) — `numpy_interop_bench.py`
2. **What does calling C++ from Python cost vs. running the same calc in
   tight C++?** (Python dispatch tax) — `gpr_overhead_bench.{cpp,py}`

Results captured here let us argue future ideas (Tier 2 batch APIs,
zero-copy aliasing, etc.) from data, not intuition, and serve as a
regression baseline.

## Run them

From the project root, with the wrapper venv active and `_ga_py` built:

```bash
# numpy buffer-protocol benchmark (Python only)
PYTHONPATH=build/ga_py:ga_py/python \
    ga_py/.venv/bin/python ga_py/benchmark/numpy_interop_bench.py

# Python-dispatch tax benchmark (C++ + Python sides)
build/ga_py/benchmark/gpr_overhead_bench
PYTHONPATH=build/ga_py:ga_py/python \
    ga_py/.venv/bin/python ga_py/benchmark/gpr_overhead_bench.py
```

To refresh a baseline (note the suffix matches the active build mode):

```bash
PYTHONPATH=build/ga_py:ga_py/python \
    ga_py/.venv/bin/python ga_py/benchmark/numpy_interop_bench.py \
    > ga_py/benchmark/results_debug.txt   # or results_release.txt

(build/ga_py/benchmark/gpr_overhead_bench && \
 PYTHONPATH=build/ga_py:ga_py/python \
   ga_py/.venv/bin/python ga_py/benchmark/gpr_overhead_bench.py) \
   > ga_py/benchmark/gpr_overhead_results_debug.txt
```

The C++ benchmark target picks up `CMAKE_BUILD_TYPE` so a Debug-vs-Debug
or Release-vs-Release comparison is fair by construction. Switch via
`cmake -DCMAKE_BUILD_TYPE=Release build && cmake --build build --target
_ga_py gpr_overhead_bench`.

## Reference results

### macOS / arm64 / Apple M-series — Python 3.14, NumPy 2.4

Captured during the §6.3 Tier 1 work.

| file                                                                   | what it captures                  |
| ---------------------------------------------------------------------- | --------------------------------- |
| [`results_debug.txt`](results_debug.txt)                               | numpy interop, `-g`               |
| [`results_release.txt`](results_release.txt)                           | numpy interop, `-O3 -DNDEBUG`     |
| [`gpr_overhead_results_debug.txt`](gpr_overhead_results_debug.txt)     | C++ vs Python gpr, `-g`           |
| [`gpr_overhead_results_release.txt`](gpr_overhead_results_release.txt) | C++ vs Python gpr, `-O3 -DNDEBUG` |

### Windows 11 x86-64 / MSVC 19.44 (VS 2022) — Python 3.13, NumPy 2.4

Captured on the same codebase, Debug and Release MSVC builds.

| file                                                                                       | what it captures                  |
| ------------------------------------------------------------------------------------------ | --------------------------------- |
| [`results_debug_windows.txt`](results_debug_windows.txt)                                   | numpy interop, MSVC Debug         |
| [`results_release_windows.txt`](results_release_windows.txt)                               | numpy interop, MSVC Release       |
| [`gpr_overhead_results_debug_windows.txt`](gpr_overhead_results_debug_windows.txt)         | C++ vs Python gpr, MSVC Debug     |
| [`gpr_overhead_results_release_windows.txt`](gpr_overhead_results_release_windows.txt)     | C++ vs Python gpr, MSVC Release   |

## Headlines — numpy buffer protocol (`numpy_interop_bench.py`)

### numpy interop — macOS arm64 (reference)

- **Per-call overhead halves under `-O3`.** `np.array(v)` is ~1.4 µs in
  Debug, ~0.7 µs in Release; manual `[v.x, v.y, v.z]` is ~0.7 µs / ~0.4 µs.
- **NEW (`np.array(v)`) loses to OLD (manual unpack) on small types in
  both modes.** vec3d is ~2× slower per call. Crossover lives near 8
  components — at mvec3d the two paths are within noise.
- **At 16 components (mvec4d), NEW matches or slightly beats OLD**, but
  only by ~10%.
- The buffer protocol's fixed cost (ndarray construction + boundary
  crossing) doesn't go away with `-O3`; only inlining and assertion removal
  shrink. Same shape across modes; Release just scales everything by ~2.

### numpy interop — Windows x86-64 (MSVC)

- **Single-call overhead is nearly identical to macOS** (Release: OLD 0.40 µs
  / NEW 0.70 µs vs macOS OLD 0.38 µs / NEW 0.71 µs). The nanobind boundary
  crossing cost is platform-independent.
- **The list→array crossover does not happen on Windows** even at 16
  components. On macOS, NEW matches OLD around mvec4d; on Windows MSVC the
  `np.stack` path remains ~60% slower than the manual unpack path even at 16
  components (Release: 11.20 ms vs 6.87 ms). The likely cause is that
  NumPy's internal array-construction path is more heavily optimized for
  ARM64 NEON than for x86-64 with MSVC's CRT.
- **The array→object (cls(row)) path** also stays slower on Windows at all
  sizes — at 16 components NEW is 3.02 ms vs OLD 2.56 ms, whereas macOS sees
  NEW edge out OLD (3.11 ms vs 3.52 ms). Same structural reason.
- **Debug vs Release scaling** is ~1.7–1.8× on Windows, consistent with
  macOS (~2×). MSVC's Debug mode is less pessimizing than Clang's.

## Headlines — Python dispatch tax (`gpr_overhead_bench.{cpp,py}`)

50 000 geometric products `vec3d * vec3d -> mvec3d_e` (a small per-call
calc: ~9 multiplies + 6 adds + 4-double result).

### gpr dispatch — macOS arm64 (reference)

| path                     | Debug ns/call | Release ns/call | Release × vs C++ |
| ------------------------ | ------------: | --------------: | ---------------: |
| C++ tight loop           |          26.7 |         **1.0** |    1× (baseline) |
| ga_py bindings           |         731.6 |            93.6 |  **~94× slower** |
| pure Python (float math) |         202.0 |           200.4 |     ~200× slower |
| numpy vectorised         |           7.5 |             8.7 |      ~9× slower  |

### gpr dispatch — Windows x86-64 / MSVC (new)

| path                     | Debug ns/call | Release ns/call | Release × vs C++ |
| ------------------------ | ------------: | --------------: | ---------------: |
| C++ tight loop           |          13.6 |         **1.1** |    1× (baseline) |
| ga_py bindings           |         462.8 |            88.9 |  **~81× slower** |
| pure Python (float math) |         146.3 |           156.0 |     ~142× slower |
| numpy vectorised         |          20.8 |            22.9 |      ~21× slower |

- **Per-call binding cost in Release is nearly identical across platforms**
  (~89 ns Windows vs ~94 ns macOS). This confirms the nanobind ABI-crossing
  overhead is structural and platform-independent.
- **Debug cost is ~37% lower on Windows** (463 ns vs 732 ns). MSVC's Debug
  mode omits fewer inlining decisions than Clang's `-g`, so the nanobind
  templated machinery runs faster even without `-O3`.
- **NumPy is ~2.4–2.6× slower on Windows x86-64** than macOS arm64 for the
  vectorised path (23 ns vs 9 ns Release). Apple's arm64 NumPy benefits from
  NEON SIMD and Apple's tuned BLAS; the x86-64 AVX path lags behind even
  though the C++ tight loop is comparable.
- **Pure Python loops run ~25% faster on Windows** (156 ns vs 200 ns
  Release), likely a CPython 3.13 vs 3.14 difference — CPython 3.13
  introduced a new specialising adaptive interpreter that cut loop overhead.
- **The Release C++ baseline (~1 ns/call) matches across platforms** —
  both compilers vectorise the loop fully. The math cost itself is not
  a differentiator.
- **The dispatch tax is per-call overhead, not math.** For larger calls
  (mvec4d × mvec4d) the tax amortises better; for `vec3d` it's the
  dominant cost.

### Where these numbers come from

The current numbers reflect a **dispatch-order optimisation** in the
emitter (`ga_bindgen/src/emit_nanobind.py::collect_binary_operators`):
within each Python dunder bucket, same-type rhs is registered first, then
`double`, then everything else. nanobind tries overloads in registration
order, so this puts the GA hot path (`vec3d * vec3d`, etc.) at the front
of the candidate list. Earlier captures showed the pre-reorder cost as
**1395 ns/call (Debug) and 140 ns/call (Release)** — i.e. the reorder
shaved ~48% in Debug and ~28% in Release, purely by skipping 5 failed
type-checks per call on the geometric-product path.

Beyond this, no further code-gen lever was identified — the remaining
~95 ns / call in Release is structural Python-ABI overhead (entry, arg
parsing, result-object allocation, refcount), which an emitter cannot
shrink. A Tier 2 batch API (one Python call doing N elements' work) is
the only path to amortise the tax.

### What this means in practice

- **For interactive or coarse-grained use of ga_py**: the dispatch tax is
  invisible. 1000 gpr calls take 0.14 ms in Release — totally fine.
- **For tight inner loops on small types in Python**: prefer numpy
  vectorisation. Build the (N, 3) arrays once (e.g. via the §6.3 buffer
  protocol), then do the math in numpy.
- **For chains of GA operations**: the tax compounds — every operator is a
  separate crossing. Either drop into C++ for the hot path, or vectorise
  with numpy.
- **A future Tier 2 batch API** (one C++ call doing N elements' work) is
  the only way to pay the dispatch tax once instead of N times. Until then,
  numpy is the closest stand-in.

See [§6.3 in the main README](../README.md#63-numpy-buffer-protocol--tier-1)
for the user-facing prose treatment.
