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

Captured during the §6.3 Tier 1 work — macOS / arm64 / Apple M-series,
Python 3.14, NumPy 2.4.

| file                                                                   | what it captures                  |
| ---------------------------------------------------------------------- | --------------------------------- |
| [`results_debug.txt`](results_debug.txt)                               | numpy interop, `-g`               |
| [`results_release.txt`](results_release.txt)                           | numpy interop, `-O3 -DNDEBUG`     |
| [`gpr_overhead_results_debug.txt`](gpr_overhead_results_debug.txt)     | C++ vs Python gpr, `-g`           |
| [`gpr_overhead_results_release.txt`](gpr_overhead_results_release.txt) | C++ vs Python gpr, `-O3 -DNDEBUG` |

## Headlines — numpy buffer protocol (`numpy_interop_bench.py`)

- **Per-call overhead halves under `-O3`.** `np.array(v)` is ~1.2 µs in
  Debug, ~0.7 µs in Release; manual `[v.x, v.y, v.z]` is ~0.5 µs / ~0.4 µs.
- **NEW (`np.array(v)`) loses to OLD (manual unpack) on small types in
  both modes.** vec3d is ~2× slower per call. Crossover lives near 8
  components — at mvec3d the two paths are within noise.
- **At 16 components (mvec4d), NEW matches or slightly beats OLD**, but
  only by ~10%.
- The buffer protocol's fixed cost (ndarray construction + boundary
  crossing) doesn't go away with `-O3`; only inlining and assertion removal
  shrink. Same shape across modes; Release just scales everything by ~2.

## Headlines — Python dispatch tax (`gpr_overhead_bench.{cpp,py}`)

50 000 geometric products `vec3d * vec3d -> mvec3d_e` (a small per-call
calc: ~9 multiplies + 6 adds + 4-double result):

| path                     | Debug ns/call | Release ns/call | Release × vs C++ |
| ------------------------ | ------------: | --------------: | ---------------: |
| C++ tight loop           |          26.7 |         **1.0** |    1× (baseline) |
| ga_py bindings           |         731.6 |            93.6 |  **~94× slower** |
| pure Python (float math) |         202.0 |           200.4 |     ~200× slower |
| numpy vectorised         |           7.5 |             8.7 |      ~9× slower  |

- **Per-call binding cost: ~94 ns in Release, ~730 ns in Debug.** The
  Release figure is dominated by the ~30 ns Python entry, ~25 ns nanobind
  arg unpack + type check, ~30 ns Python result-object alloc, and ~10 ns
  refcount/return. Math itself is ~3 ns. Build mode scales everything
  ~7–10× because nanobind's templated machinery inlines aggressively under
  `-O3`.
- **The Release C++ baseline (1.0 ns/call) is essentially loop overhead +
  one rolled-up SIMD-fused gpr** — the compiler vectorises the loop. For
  per-element comparison purposes, treat it as the lower bound the math
  alone would cost.
- **NumPy vectorised wins decisively whenever you can batch** — within
  ~9× of tight C++ in Release, and unaffected by build mode (it's already
  optimised C). For tight loops on small types in Python, this is the
  fast path.
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
