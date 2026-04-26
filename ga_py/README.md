# ga_py — Python wrapper around the `ga/` library

`ga_py` exposes the geometric algebra library as a Python module. Types (`vec3d`, `bivec3d`, `mvec3d`, …), operators (`+`, `-`, `*`, `^`, `<<`, `>>`), and free functions (`gpr`, `wdg`, `dot`, `nrm`, …) all map 1:1 to their C++ counterparts.

```python
import ga_py

v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
B = ga_py.ega.bivec3d(0.0, 0.0, 1.0)
m = v + B                                    # cross-grade addition → mvec3d
print(m, m.gr1(), m.gr2())                   # extract grade components
```

The package is **generated** from the C++ headers by `ga_bindgen/` — see [TODO/considerations_python_wrapper.md](../TODO/considerations_python_wrapper.md) for the design and [ga_bindgen/README.md](../ga_bindgen/README.md) for the toolchain.

## Two virtual environments — don't mix them

This project keeps the wrapper and the binding generator in separate venvs. They serve different roles, have different dependencies, and must not be conflated:

| Venv path | Used for | Dependencies |
| --- | --- | --- |
| `ga_py/.venv` | **Running the wrapper and its tests.** This is the venv for everyday development of `ga_py/` itself. | `pytest`, `hypothesis`, `numpy` |
| `build/spike_libclang/.venv` | **Regenerating bindings via `ga_bindgen`.** Only contributors touching `ga/*.hpp` or the generator need this one. | `libclang`, `jinja2`, `nanobind` |

If you run `pytest ga_py/tests/` from the wrong venv, `pytest` won't be on PATH (or `hypothesis` will be missing) — `conftest.py` prints a hint pointing you here. If you run `python ga_bindgen/src/scan.py` from the wrong venv, libclang won't be importable.

**Activate the right one for the task:**

```bash
# For wrapper development / running tests:
source ga_py/.venv/bin/activate

# For regenerating bindings (separate session — don't reuse the above):
source build/spike_libclang/.venv/bin/activate
```

The `ga_py/.venv` is created once via:

```bash
python3 -m venv ga_py/.venv
ga_py/.venv/bin/pip install pytest hypothesis numpy
```

## Performance — what's compiled, what's interpreted

A common worry about Python wrappers is that the math will be re-implemented in interpreted Python and pay an interpreter tax per operation. **That is not the case here.**

When Python evaluates `v1 + v2`:

```
Python:  v1 + v2
            │
            ▼  nanobind dispatch — type check + unwrap (one C call, ~few hundred ns)
C++:     [](vec3d const& a, vec3d const& b) { return a + b; }
            │
            ▼  compiled C++ — ADL resolves to operator+ for Vec3_t<double, vec3d_tag>
C++:     hd::ga::operator+(...)
            │
            ▼  inlined and optimized exactly like a normal C++ call
            
Result wrapped back into a Python `vec3d` object.
```

The lambdas in [src/generated/](src/generated/) are real C++ code, compiled with full optimization, instantiating the GA templates for `value_t = double`. There is **no Python re-implementation of the math** anywhere. The only per-call overhead is the Python ↔ C nanobind dispatch, which is unavoidable for any Python binding regardless of how the C++ is packaged.

### "Why isn't there a separate `libga_compiled.a` that the wrapper links against?"

It would not change runtime speed. The same template instantiations would run with the same machine code whether they live inside `_ga_py.so` directly (today) or in a separate static archive (alternative). The only thing such a library can improve is **compile time**: today every `bindings_*.cpp` re-instantiates the same `Vec3_t<double, vec3d_tag>` etc., and a precompiled archive could do the work once. For the current ~50 binding files compiling in ~30s, this is not yet painful, so the precompiled-archive optimization is deferred. See [TODO/considerations_python_wrapper.md §10b](../TODO/considerations_python_wrapper.md) for the location-of-such-a-library question (which would land in `ga_py/` first if it ever does).

### Where Python *would* be slower than raw C++

Two places — both unrelated to whether the math itself is "compiled":

1. **Per-call dispatch overhead.** Each Python → C++ crossing costs a few hundred ns. For million-call inner loops on tiny types (e.g., adding millions of `vec3d`s in a Python `for` loop), this dominates. The fix is **batching** (e.g., numpy-array-of-multivectors operations), which is the Tier 2 numpy work deferred per the design doc.
2. **Object lifetime overhead.** Every result is a new Python object with a refcounted handle. Negligible for normal use; adds up in tight loops.

Neither is fixed by precompiled libraries; both are fixed (when needed) by vectorized array APIs.

## Coming from EGA, working in PGA — naming differences

The `ga_py.ega` and `ga_py.pga` submodules deliberately do **not** mirror each other one-to-one. PGA's degenerate metric (one null basis vector) splits several Euclidean concepts into two halves, and PGA models *rigid motions* rather than *rotations*. If you transfer an EGA workflow naively into `ga_py.pga`, several familiar function names will not be there. Calling them raises `TypeError`. The differences below are by C++ design, not binding gaps:

| EGA function | PGA equivalent(s) | Why different |
| --- | --- | --- |
| `nrm`, `nrm_sq`, `normalize` | `bulk_nrm`, `weight_nrm`, `geom_nrm` (and `_sq` variants); `bulk_normalize`, `unitize` | Degenerate metric splits the norm. You must pick which part you want: bulk (Euclidean part), weight (ideal part), or the geometric norm combining both. `unitize` enforces weight = 1. |
| `dual`, `l_dual`, `r_dual` | `bulk_dual`, `weight_dual`, `l_bulk_dual`, `r_bulk_dual`, `l_weight_dual`, `r_weight_dual` | Same reason: dualization splits along the bulk/weight decomposition. |
| `rotate`, `rotate_opt` | `move2dp`, `move3dp`, `move2dp_opt`, `move3dp_opt` | PGA composes rotations and translations into a single *motor*, applied via `move*`. The naming is honest — it is a rigid motion, not just a rotation. |
| `cross` | (none) | The 3D vector cross product is not generally defined on PGA elements. Use `wdg` (wedge) followed by an appropriate complement. |
| `gs_orthogonal`, `gs_orthonormal` | (none) | Gram–Schmidt requires the Euclidean inner product. |
| `get_rotor` | `get_motor`, `get_motor_from_lines`, `get_motor_from_planes` | Rotors live in EGA; motors are their PGA generalization with translation. |

PGA also adds many functions with no EGA counterpart (`join`, `meet`, `att`, `expand`, `bulk`, `weight`, `unitize`, `rgpr`, `rrev`, `rdot`, `rwdg`, `rcmt`, the `central_proj*`/`ortho_proj*`/`ortho_antiproj*` families, the `*_contract*` and `*_expand*` families, `get_motor_from_*`, ...). These follow PGA's geometric-incidence and bulk/weight conventions; see the C++ headers (`ga/ga_pga*.hpp`) for the definitions.

Functions that exist in **both** submodules with the same name (`wdg`, `dot`, `cmpl`, `rev`, `inv`, `gr_inv`, `cmt`, `conj`, `exp`, `sqrt`, `angle`, `is_congruent`, `project_onto`, `reject_from`, `reflect_on`, `rwdg`, `rtwdg1`, `twdg1`, `l_cmpl`, `r_cmpl`) dispatch on argument type — calling `ga_py.ega.wdg(vec3d, vec3d)` and `ga_py.pga.wdg(vec3dp, vec3dp)` both Just Work; their meanings are analogous but the resulting types live in different algebras.

## Gotchas worth knowing up front

These are not bugs — they are consequences of mirroring a strongly-typed C++ library 1:1 in Python. You will hit them sooner or later; better to know about them on day 1.

### 1. Typed scalars vs. plain floats

Some "scalar-valued" operations return a **typed wrapper** (`scalar3d`, `pscalar3dp`, ...) while others return a **plain Python `float`**. The split mirrors the C++ signature:

- C++ `Scalar3d<T>` / `PScalar3d<T>` return → Python typed wrapper.
- C++ `T` return → Python `float`.

| Returns plain `float` | Returns a typed wrapper |
| --- | --- |
| `ega.nrm(v)`, `ega.nrm_sq(v)` | `ega.dot(v, v)` → `scalar3d` |
| `pga.bulk_nrm_sq(p)`, `pga.weight_nrm_sq(p)` | `pga.dot(p, q)` → `scalar3dp` |
| | `pga.bulk_nrm(p)` → `scalar3dp` |
| | `pga.weight_nrm(p)` → `pscalar3dp` |
| | `pga.att(p)`, `gr0(M)`, `gr1(M)`, ... |
| | `pga.geom_nrm(p)` → `dualnum3dp` (bulk + weight together) |

Typed wrappers preserve the C++ type-safety guarantees — you cannot accidentally mix a `scalar2d` with a `pscalar3d` in arithmetic. They also do **not** participate in arithmetic with Python numbers:

```python
import ga_py
v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
ns = ga_py.ega.nrm_sq(v)       # plain float (14.0) — fine to do arithmetic
d  = ga_py.ega.dot(v, v)       # scalar3d(14) — typed; needs float() to mix
ns + 1.0                       # 15.0  ✓
d + 1.0                        # TypeError
float(d) + 1.0                 # 15.0  ✓
```

If you need a plain Python number from a typed wrapper — for `assert n < 1e-9`, for `math.sqrt(n)`, for `+`/`-`/`*` with a Python number — wrap in `float(...)`. The bound types implement `__float__`, so the conversion is cheap and explicit.

### 2. Wedge `^` has lower precedence than `+` and `*` in Python

C++ users habitually write `a + b ^ c` expecting it to mean `a + (b ^ c)`. In Python, `^` is XOR's symbol and binds *less tightly* than `+` and `*`:

```python
a + b ^ c       # Python parses as (a + b) ^ c   — different from the C++ habit
a * b ^ c       # Python parses as (a * b) ^ c
v << b ^ c      # Python parses as (v << b) ^ c
```

Defensively parenthesize wedges:

```python
a + (b ^ c)     # always safe
(v >> u) + (u ^ v)
```

This is not a bug in `ga_py` — it is Python's operator table, the same for all programs that use `^`. Surface it in your own examples.

### 3. Nanobind submodule attributes don't behave like Python submodules

`ga_py.ega` and `ga_py.pga` are exposed as attributes on the `ga_py` package (they are `nanobind` sub-modules of the C extension). They work fine for attribute access:

```python
import ga_py
v = ga_py.ega.vec3d(1, 2, 3)        # works
e = ga_py.ega                       # works as alias
```

But you cannot use the `from ... import` form, because they are not registered in `sys.modules`:

```python
import ga_py.ega                    # ModuleNotFoundError
from ga_py import ega               # works (attribute access, not submodule import)
```

If you want `ega` as a short name, use `from ga_py import ega` or `e = ga_py.ega`.

## Layout

```text
ga_py/
├── CMakeLists.txt              # nanobind + scikit-build target
├── README.md                   # this file
├── src/
│   ├── module.cpp              # hand-written entry point (NB_MODULE)
│   └── generated/              # produced by ga_bindgen — do not edit
│       ├── bindings_*.cpp      # one file per bound user type
│       ├── register_all.cpp    # dispatcher in inheritance-topological order
│       └── bindings_list.cmake # generated source list, included by CMakeLists.txt
├── python/
│   └── ga_py/__init__.py       # re-exports ega/pga submodules
└── tests/
    ├── conftest.py             # tolerance fixtures, component extraction
    ├── test_constants.py       # T1: constants verification (143 tests)
    ├── test_grade_lookup.py    # T2: gr() / rgr() coverage (30 tests)
    ├── test_identities_ega.py  # T3: EGA algebraic identities via hypothesis (15)
    ├── test_identities_pga.py  # T4: PGA algebraic identities via hypothesis (20)
    ├── test_cross_check.py     # T5: replays C++ reference cases via JSON
    └── data/
        └── ga_test_cases.json  # captured (op, args, expected) tuples from C++
```

## Cross-check tests against C++ (T5)

`tests/test_cross_check.py` loads a JSON snapshot of (operation, inputs, expected output) tuples emitted by the C++ side and verifies the Python wrapper produces bit-equivalent results within tolerance. This is the most direct equivalence test in the suite.

The C++ exporter is `ga_test/src/export_python_cases.cpp`, an executable separate from the doctest binaries. It uses the production `ga/` headers, runs each operation, and writes the result.

**Regenerating the JSON snapshot** (only needed when the case set changes — new ops, new types, new representative inputs; the snapshot is committed to git so end users don't need to regenerate it):

```bash
cmake --build build --target ga_export_python_cases
build/ga_test/ga_export_python_cases ga_py/tests/data/ga_test_cases.json
```

After regenerating, run `pytest ga_py/tests/test_cross_check.py` to verify the wrapper still matches.

## Building

### Recommended: top-level integrated build

The Python wrapper is gated by the `_GA_BUILD_PYTHON` CMake option (off by default). Enable it once and the wrapper builds alongside everything else:

```bash
# One-time wrapper venv setup (also where nanobind lives for the build):
python3 -m venv ga_py/.venv
ga_py/.venv/bin/pip install nanobind pytest hypothesis numpy

# Configure + build (from project root):
cmake -S . -B build -D_GA_BUILD_PYTHON=ON
cmake --build build
```

`ga_py/CMakeLists.txt` automatically detects `ga_py/.venv/bin/python` and uses its `nanobind` — no need to set `Python_EXECUTABLE` or `CMAKE_PREFIX_PATH` manually.

Run tests:

```bash
PYTHONPATH="build/ga_py:ga_py/python" ga_py/.venv/bin/pytest ga_py/tests/
```

### Alternative: standalone build (just the wrapper)

For quick iteration on the wrapper alone, without re-configuring the rest of the project:

```bash
cmake -S ga_py -B build/ga_py
cmake --build build/ga_py
```

Same auto-detection of `ga_py/.venv` applies.

## Regenerating bindings after C++ changes

When `ga/` gains new types, operators, or functions, regenerate from the headers:

```bash
python3 ga_bindgen/src/scan.py                  # re-scan headers → manifest.json
python3 ga_bindgen/src/emit_nanobind.py --all   # emit all types + all free functions
cmake --build build                             # recompile
```

Generated files are committed to git so PR diffs show exactly what new Python surface a C++ change exposes.

## Refreshing cross-check test data

`tests/test_cross_check.py` replays a JSON snapshot of (op, args, expected) tuples emitted by C++. The snapshot is committed; refresh it via the convenience target whenever the case set changes:

```bash
cmake --build build --target regenerate_python_test_data
```

This builds `ga_export_python_cases` (an executable that walks representative GA operations through the C++ library) and runs it, writing `ga_py/tests/data/ga_test_cases.json`. Both the executable and the convenience target only exist when `_GA_BUILD_PYTHON=ON`.
