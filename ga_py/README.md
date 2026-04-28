# ga_py — Python wrapper for the geometric algebra library

## 1. What is ga_py?

`ga_py` is a Python module that exposes the C++ geometric algebra library (`ga/`) to Python. Types (`vec3d`, `bivec3d`, `mvec3d`, …), operators (`+`, `-`, `*`, `^`, `<<`, `>>`), and free functions (`gpr`, `wdg`, `dot`, `nrm`, …) all map 1:1 to their C++ counterparts. The two algebras live in submodules:

- `ga_py.ega` — Euclidean GA (mirrors `hd::ga::ega`)
- `ga_py.pga` — Projective GA (mirrors `hd::ga::pga`)

Free functions in the `hd::ga` top-level namespace are re-exported at `ga_py.*` (so `ga_py.deg2rad(180.0)` works).

```python
import ga_py

v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
B = ga_py.ega.bivec3d(0.0, 0.0, 1.0)
m = v + B                                    # cross-grade addition → mvec3d
print(m, m.gr1(), m.gr2())                   # extract grade components
```

The value type is fixed to `double` in v1 (see [§6.4](#64-float-value_t-support)).

## 2. Relationship to the C++ library

`ga/` is the source of truth. `ga_py` is a real C++ extension module compiled from bindings that are *generated* from the `ga/` headers — there is no Python re-implementation of the math.

```text
ga/*.hpp ──libclang──► ga_bindgen/manifest.json ──emitter──► ga_py/src/generated/*.cpp ──nanobind+C++ compile──► _ga_py.so / _ga_py.pyd
```

When Python evaluates `v1 + v2`, nanobind dispatches to a generated lambda that calls the corresponding `hd::ga` operator — compiled, inlined, and optimized exactly like a normal C++ call. The only per-call cost is the Python ↔ C nanobind crossing (a few hundred ns), unavoidable for any Python binding. For details see [§5.5](#55-performance--whats-compiled-vs-interpreted).

The binding generator lives in [`ga_bindgen/`](../ga_bindgen/README.md). End users do not need it; it is run by contributors when the C++ headers change. Generated bindings are committed to git so end users build `ga_py` without libclang installed.

## 3. Quick start — running your first Python script

This section assumes you only want to *use* the library. You do **not** need libclang/LLVM, nor any of the binding-generation toolchain.

### 3.1 Prerequisites

| Platform | Required |
| --- | --- |
| **macOS** | `brew install fmt` · Python 3.10+ (`brew install python`) · CMake 3.28+ |
| **Linux** | `fmt` (apt/distro pkg) · Python 3.10+ · CMake 3.28+ |
| **Windows** | [vcpkg](https://github.com/microsoft/vcpkg) at `C:/vcpkg/` with `vcpkg install fmt` · Python 3.10+ on PATH · MSVC / Visual Studio 2022 · CMake 3.28+ |

### 3.2 One-time setup

> **Working directory for sections 3.2 – 3.5: the project root** — i.e. the directory containing `ga/`, `ga_py/`, `build/`, `CMakeLists.txt`, and this README's parent. All commands below use paths relative to that root; running them from anywhere else will fail with `file not found` or build into the wrong location.

```bash
# 1) Clone the repo and cd into the project root (this is the cwd for everything below)
git clone <repo-url> ga
cd ga                                  # project root — keep this as cwd

# 2) Create the wrapper venv and install runtime + test dependencies
#    (run from the project root)
python3 -m venv ga_py/.venv
ga_py/.venv/bin/pip install nanobind pytest hypothesis numpy

# 3) Configure and build with the Python wrapper enabled
#    (run from the project root)
cmake -S . -B build -D_GA_BUILD_PYTHON=ON
cmake --build build
```

`ga_py/CMakeLists.txt` automatically picks up `ga_py/.venv/bin/python` and uses its `nanobind` — no need to set `Python_EXECUTABLE` or `CMAKE_PREFIX_PATH`.

**Windows** (Command Prompt) — same idea, different paths. cwd is still the project root (e.g. `C:\path\to\ga`):

```bat
rem cwd: project root, e.g.  cd C:\path\to\ga
python -m venv ga_py\.venv
ga_py\.venv\Scripts\pip install nanobind pytest hypothesis numpy
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake -D_GA_BUILD_PYTHON=ON
cmake --build build --target _ga_py --config Debug
```

After the build, the compiled extension lands at (paths relative to the project root):

- macOS / Linux: `build/ga_py/_ga_py.<tag>.so`
- Windows: `build\ga_py\Debug\_ga_py.cp3xx-win_amd64.pyd`

### 3.3 Make Python find the module

There is no `pip install .` / wheel workflow yet (see [§6.1](#61-packaging-as-pip-installable-wheels)). `PYTHONPATH` must cover **both** the directory holding the compiled `_ga_py` extension and the `ga_py/python` directory that holds the pure-Python `ga_py/` package wrapper. The relative paths below assume the project root as cwd:

```bash
# macOS / Linux (bash / zsh) — cwd: project root
export PYTHONPATH="$PWD/build/ga_py:$PWD/ga_py/python"
```

```bat
rem Windows Command Prompt — cwd: project root
set PYTHONPATH=%CD%\build\ga_py\Debug;%CD%\ga_py\python
```

```powershell
# Windows PowerShell — cwd: project root
$env:PYTHONPATH = "$PWD\build\ga_py\Debug;$PWD\ga_py\python"
```

(Using `$PWD` / `%CD%` resolves the relative paths to absolute ones, so `PYTHONPATH` keeps working if you later `cd` elsewhere in the same shell.)

### 3.4 Hello-world script

A ready-to-run demo lives at [`ga_py/demo/hello_ga.py`](demo/hello_ga.py). It exercises the two core GA operations — the wedge product and the geometric product — and verifies the foundational identity `v1 * v2 = dot(v1, v2) + wdg(v1, v2)`. Run it from the project root:

```bash
# cwd: project root
PYTHONPATH="$PWD/build/ga_py:$PWD/ga_py/python" \
    ga_py/.venv/bin/python ga_py/demo/hello_ga.py
```

(If you've already exported `PYTHONPATH` per §3.3, the leading `PYTHONPATH=…` is redundant.)

Expected output:

```text
m = MVec3d(0,1,2,3,0,0,1,0)
|v| = 3.7416573867739413
v1 ^ v2 = BiVec3d(1,-1.5,2)
v1 * v2 = MVec3d_E(8.5,1,-1.5,2)
gr0(v1*v2) = 8.5   dot(v1, v2) = 8.5
gr2(v1*v2) = BiVec3d(1,-1.5,2)    wdg(v1, v2) = BiVec3d(1,-1.5,2)
line = BiVec2dp(3,4,-2)
```

If this runs without `ImportError` and the `gr0`/`gr2` values match `dot`/`wdg`, the install is good.

### 3.5 Verify with the test suite

```bash
# cwd: project root
PYTHONPATH="build/ga_py:ga_py/python" ga_py/.venv/bin/pytest ga_py/tests/
```

411 tests, ~1.5 s. They cover constants, grade lookup, EGA/PGA algebraic identities (via `hypothesis`), and a JSON cross-check against the C++ reference outputs.

### 3.6 Where to find the API reference

The Python surface mirrors the C++ headers exactly. Until `.pyi` stub generation lands ([§6.2](#62-stub-generation-pyi-and-mypy)), the authoritative reference is the C++ source:

- EGA: `ga/ga_ega2d_ops.hpp`, `ga/ga_ega3d_ops.hpp`, `ga/ga_usr_consts.hpp`
- PGA: `ga/ga_pga2dp_ops.hpp`, `ga/ga_pga3dp_ops.hpp`, `ga/ga_usr_consts.hpp`

Use `dir(ga_py.ega)` / `dir(ga_py.pga)` for a quick interactive listing.

## 4. Using the library

### 4.1 Constants and basis vectors

Every constant declared in `ga/ga_usr_consts.hpp` is exposed as a module attribute on the matching submodule:

```python
v = 2 * ga_py.ega.e1 + 3 * ga_py.ega.e2     # 2*e1 + 3*e2
I = ga_py.ega.I_3d                          # pseudoscalar
```

107 constants are bound today; new ones in `ga_usr_consts.hpp` appear automatically when bindings are regenerated.

### 4.2 Coming from EGA, working in PGA — naming differences

The `ga_py.ega` and `ga_py.pga` submodules deliberately do **not** mirror each other one-to-one. PGA's degenerate metric (one null basis vector) splits several Euclidean concepts into two halves, and PGA models *rigid motions* rather than *rotations*. If you transfer an EGA workflow naively into `ga_py.pga`, several familiar function names will not be there. Calling them raises `TypeError`. The differences below are by C++ design, not binding gaps:

| EGA function | PGA equivalent(s) | Why different |
| --- | --- | --- |
| `nrm`, `nrm_sq`, `normalize` | `bulk_nrm`, `weight_nrm`, `geom_nrm` (and `_sq` variants); `bulk_normalize`, `unitize` | Degenerate metric splits the norm. You must pick which part you want: bulk (Euclidean part), weight (ideal part), or the geometric norm combining both. `unitize` enforces weight = 1. |
| `dual`, `l_dual`, `r_dual` | `bulk_dual`, `weight_dual`, `l_bulk_dual`, `r_bulk_dual`, `l_weight_dual`, `r_weight_dual` | Same reason: dualization splits along the bulk/weight decomposition. |
| `rotate`, `rotate_opt` | `move2dp`, `move3dp`, `move2dp_opt`, `move3dp_opt` | PGA composes rotations and translations into a single *motor*, applied via `move*`. The naming is honest — it is a rigid motion, not just a rotation. |
| `cross` | (none) | The 3D vector cross product is not generally defined on PGA elements. Use `wdg` (wedge) followed by an appropriate complement. |
| `gs_orthogonal`, `gs_orthonormal` | (none) | Gram–Schmidt requires the Euclidean inner product. |
| `get_rotor` | `get_motor`, `get_motor_from_lines`, `get_motor_from_planes` | Rotors live in EGA; motors are their PGA generalization with translation. |

PGA also adds many functions with no EGA counterpart (`join`, `meet`, `att`, `expand`, `bulk`, `weight`, `unitize`, `rgpr`, `rrev`, `rdot`, `rwdg`, `rcmt`, the `central_proj*`/`ortho_proj*`/`ortho_antiproj*` families, the `*_contract*` and `*_expand*` families, `get_motor_from_*`, …). These follow PGA's geometric-incidence and bulk/weight conventions; see the C++ headers (`ga/ga_pga*.hpp`) for the definitions.

Functions that exist in **both** submodules with the same name (`wdg`, `dot`, `cmpl`, `rev`, `inv`, `gr_inv`, `cmt`, `conj`, `exp`, `sqrt`, `angle`, `is_congruent`, `project_onto`, `reject_from`, `reflect_on`, `rwdg`, `rtwdg1`, `twdg1`, `l_cmpl`, `r_cmpl`) dispatch on argument type — calling `ga_py.ega.wdg(vec3d, vec3d)` and `ga_py.pga.wdg(vec3dp, vec3dp)` both Just Work; their meanings are analogous but the resulting types live in different algebras.

### 4.3 Gotchas worth knowing up front

These are not bugs — they are consequences of mirroring a strongly-typed C++ library 1:1 in Python. You will hit them sooner or later; better to know about them on day 1.

#### Typed scalars vs. plain floats

Some "scalar-valued" operations return a **typed wrapper** (`scalar3d`, `pscalar3dp`, …) while others return a **plain Python `float`**. The split mirrors the C++ signature:

- C++ `Scalar3d<T>` / `PScalar3d<T>` return → Python typed wrapper.
- C++ `T` return → Python `float`.

| Returns plain `float` | Returns a typed wrapper |
| --- | --- |
| `ega.nrm(v)`, `ega.nrm_sq(v)` | `ega.dot(v, v)` → `scalar3d` |
| `pga.bulk_nrm_sq(p)`, `pga.weight_nrm_sq(p)` | `pga.dot(p, q)` → `scalar3dp` |
| | `pga.bulk_nrm(p)` → `scalar3dp` |
| | `pga.weight_nrm(p)` → `pscalar3dp` |
| | `pga.att(p)`, `gr0(M)`, `gr1(M)`, … |
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

#### Wedge `^` has lower precedence than `+` and `*` in Python

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

#### Nanobind submodule attributes don't behave like Python submodules

`ga_py.ega` and `ga_py.pga` are exposed as attributes on the `ga_py` package (they are `nanobind` sub-modules of the C extension). They work fine for attribute access:

```python
import ga_py
v = ga_py.ega.vec3d(1, 2, 3)        # works
e = ga_py.ega                       # works as alias
```

But you cannot use the `from … import` form, because they are not registered in `sys.modules`:

```python
import ga_py.ega                    # ModuleNotFoundError
from ga_py import ega               # works (attribute access, not submodule import)
```

If you want `ega` as a short name, use `from ga_py import ega` or `e = ga_py.ega`.

## 5. For contributors

This section covers everyday wrapper development, regenerating bindings after C++ changes, and the technical background on why the math is fast. End users do not need any of it.

### 5.1 The two virtual environments — don't mix them

Wrapper development and binding regeneration use separate venvs with different dependencies:

| Venv path | Used for | Dependencies |
| --- | --- | --- |
| `ga_py/.venv` | Running the wrapper and its tests; everyday development of `ga_py/` itself | `nanobind`, `pytest`, `hypothesis`, `numpy` |
| `build/spike_libclang/.venv` | Regenerating bindings via `ga_bindgen` (only contributors touching `ga/*.hpp` or the generator need this) | `libclang` |

If you run `pytest ga_py/tests/` from the wrong venv, `pytest` won't be on PATH (or `hypothesis` will be missing) — `conftest.py` prints a hint pointing here. If you run `python ga_bindgen/src/scan.py` from the wrong venv, `libclang` won't be importable.

```bash
# macOS / Linux
source ga_py/.venv/bin/activate           # wrapper / tests
source build/spike_libclang/.venv/bin/activate   # binding regeneration (separate session)
```

```bat
rem Windows
ga_py\.venv\Scripts\activate
build\spike_libclang\.venv\Scripts\activate
```

### 5.2 Layout

```text
ga_py/
├── CMakeLists.txt              # nanobind + scikit-build target
├── README.md                   # this file
├── src/
│   ├── module.cpp              # hand-written entry point (NB_MODULE)
│   └── generated/              # produced by ga_bindgen — do not edit
│       ├── bindings_*.cpp      # one file per bound user type (46 today)
│       ├── bindings_constants_*.cpp
│       ├── bindings_functions_*.cpp
│       ├── register_all.cpp    # dispatcher in inheritance-topological order
│       └── bindings_list.cmake # generated source list, included by CMakeLists.txt
├── python/
│   └── ga_py/__init__.py       # re-exports ega/pga submodules + Python forwarders
│                               #   for fully-generic C++ templates that the
│                               #   binding generator cannot enumerate
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

The `python/ga_py/__init__.py` shim contains pure-Python forwarders for a handful of `ga/` free functions written as fully-generic C++ templates (`template<typename A, typename B> auto f(A&&, B&&)`) — there is no overload list for the generator to enumerate. The forwarders compose already-bound C++ functions; they contain **no Python math**.

### 5.3 Build modes

#### Recommended: top-level integrated build

The Python wrapper is gated by the `_GA_BUILD_PYTHON` CMake option (off by default). Enable it once and the wrapper builds alongside everything else:

```bash
cmake -S . -B build -D_GA_BUILD_PYTHON=ON
cmake --build build
```

#### Standalone wrapper build (faster iteration)

```bash
cmake -S ga_py -B build/ga_py
cmake --build build/ga_py
```

Same auto-detection of `ga_py/.venv` applies.

### 5.4 Regenerating bindings after C++ changes

When `ga/` gains new types, operators, or functions:

```bash
source build/spike_libclang/.venv/bin/activate     # libclang venv
python3 ga_bindgen/src/scan.py                     # ga/*.hpp → manifest.json
python3 ga_bindgen/src/emit_nanobind.py --all      # manifest → bindings_*.cpp
cmake --build build                                # recompile
```

Generated files are committed so PR diffs show exactly what new Python surface a C++ change exposes. Full toolchain reference: [`ga_bindgen/README.md`](../ga_bindgen/README.md).

### 5.5 Performance — what's compiled vs interpreted

A common worry about Python wrappers is that the math will be re-implemented in interpreted Python and pay an interpreter tax per operation. **That is not the case here.**

When Python evaluates `v1 + v2`:

```text
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

The lambdas in `src/generated/` are real C++ code, compiled with full optimization, instantiating the GA templates for `value_t = double`. There is **no Python re-implementation of the math** anywhere. The only per-call overhead is the Python ↔ C nanobind dispatch, which is unavoidable for any Python binding regardless of how the C++ is packaged.

#### Why isn't there a separate `libga_compiled.a` that the wrapper links against?

It would not change runtime speed. The same template instantiations would run with the same machine code whether they live inside `_ga_py.so` directly (today) or in a separate static archive (alternative). The only thing such a library can improve is **compile time**: today every `bindings_*.cpp` re-instantiates the same `Vec3_t<double, vec3d_tag>` etc., and a precompiled archive could do the work once. For the current ~50 binding files compiling in ~30 s, this is not yet painful, so the precompiled-archive optimization is deferred.

#### Where Python *would* be slower than raw C++

Two places — both unrelated to whether the math itself is "compiled":

1. **Per-call dispatch overhead.** Each Python → C++ crossing costs a few hundred ns. For million-call inner loops on tiny types (e.g., adding millions of `vec3d`s in a Python `for` loop), this dominates. The fix is **batching**, addressed by the Tier 1 numpy buffer-protocol work in [§6.3](#63-numpy-buffer-protocol--tier-1).
2. **Object lifetime overhead.** Every result is a new Python object with a refcounted handle. Negligible for normal use; adds up in tight loops.

Neither is fixed by precompiled libraries; both are fixed (when needed) by vectorized array APIs.

### 5.6 Refreshing cross-check test data (T5)

`tests/test_cross_check.py` replays a JSON snapshot of `(op, args, expected)` tuples emitted by the C++ side and verifies the Python wrapper produces bit-equivalent results within tolerance. The C++ exporter is `ga_test/src/export_python_cases.cpp`, an executable separate from the doctest binaries that uses the production `ga/` headers.

The snapshot is committed to git. Refresh it whenever the case set changes:

```bash
cmake --build build --target regenerate_python_test_data
```

This builds `ga_export_python_cases` and runs it, writing `ga_py/tests/data/ga_test_cases.json`. Both the executable and the convenience target only exist when `_GA_BUILD_PYTHON=ON`. After regenerating, run `pytest ga_py/tests/test_cross_check.py` to verify.

## 6. Open items / future work

Items deliberately deferred or not started in v1. Listed in rough priority order.

### 6.1 Packaging as pip-installable wheels

The wrapper is currently built via direct CMake; there is no `pip install ga_py` flow. To get there:

- Add `pyproject.toml` with `scikit-build-core` as the build backend; declare `nanobind` (build-only) and `numpy` (runtime, for buffer protocol) as dependencies.
- Set up `cibuildwheel` on GitHub Actions for Python 3.10 × {macOS Intel + ARM, Linux x86_64, Windows x86_64} = 4 wheels.
- Publish a source distribution (sdist) too.
- Additional Python versions (3.11/3.12/3.13) extend the cibuildwheel matrix without design changes.

`nanobind`, `fmt`, and `ga/` are all header-only — no vcpkg/brew dependencies on the wheel-building side, so CI stays simple.

### 6.2 Stub generation `.pyi` and mypy

`ga_bindgen/src/emit_stubs.py` does not yet exist. Once added it would:

- Drive the IDE / type-check experience for users (autocomplete, hover types).
- Unblock a `mypy --strict` test pass on the example scripts (currently held until stubs exist).

### 6.3 numpy buffer protocol — Tier 1

Each bound user type should expose `__array__` / the buffer protocol so numpy arrays of GA elements can be constructed cheaply. Tier 2 (vectorized array-of-multivector ops) remains deferred indefinitely — Tier 1 is the cheap, broadly useful piece.

### 6.4 Float `value_t` support

Currently locked to `double`. Float (and possibly `long double`) support would require a build-flag-driven manifest instantiation, doubling the binding surface. Revisit if user demand justifies it.

### 6.5 `rk4_step` deferral and the scipy recipe

The C++ free function `hd::ga::rk4_step<VecType>(…)` (Runge–Kutta 4th-order time integrator, declared in `ga/ga_usr_utilities.hpp`) is **not bound** in `ga_py`. Three obstacles, in increasing order of difficulty:

1. **`std::mdspan` has no nanobind type-caster.** All three array arguments (`u`, `uh`, `rhs`) are `std::mdspan<VecType, std::dextents<size_t, N>>`. nanobind ships casters for `std::vector`, `std::array`, numpy arrays via the buffer protocol — but not for C++23 `std::mdspan` (a non-owning view that the binding layer would have to materialise over caller-supplied storage). The bindgen sees the parameter, fails type resolution, and silently skips the function.
2. **Mutates its arguments instead of returning a value.** `u` and `uh` are written through their mdspans — idiomatic for tight numerical loops but doesn't fit the Python convention of returning new state.
3. **Templated on `VecType`.** Each GA vector type (`vec2d`, `vec3d`, `vec2dp`, `vec3dp`, …) gives a different concrete instantiation that has to be enumerated and bound separately. Least hard of the three — the same machinery already handles `gpr<T>`, `wdg<T>`, etc.

**Do we even need to bind it?** Probably not. Python users who need to integrate ODEs over GA-valued state already have superior tools:

- `scipy.integrate.solve_ivp` — modern unified ODE solver with adaptive step size, dense output, event detection. Methods: `RK45` (default; Dormand–Prince adaptive 5th order with embedded 4th-order error estimator), `RK23`, `DOP853`, `Radau`, `BDF`, `LSODA`.
- `scipy.integrate.RK45` as a class for lower-level sub-step control.
- A 6-line numpy fixed-step RK4 if a user really wants non-adaptive RK4.
- `numba` / `jax` JIT-compiled integrators if scipy isn't fast enough.

The natural Python workflow:

1. Convert GA-valued state to a flat numpy array (each `vec3d` becomes 3 doubles, etc.).
2. Wrap the GA-valued `rhs(t, u)` in a closure that unpacks/repacks numpy ↔ GA objects.
3. Call `solve_ivp(rhs, t_span, u0, method="RK45")`.
4. Reconstruct GA objects from the result.

The packing/unpacking is per-component and becomes essentially zero-overhead once Tier 1 numpy buffer protocol ([§6.3](#63-numpy-buffer-protocol--tier-1)) lands.

**Recommended path: defer permanently.** Concrete cleanup:

- [ ] Add `rk4_get_time` to `SKIP_FREE_FN_NAMES` in `ga_bindgen/src/emit_nanobind.py` (it's bound today but useless without `rk4_step`).
- [ ] Regenerate bindings, verify the function is gone.
- [ ] Add a "Integrating GA-valued ODEs in Python" section to this README with a `solve_ivp` example using a `vec3d`-valued state vector.

If demand later justifies a native binding, the right approach is a hand-written `nb::ndarray<vec3d, nb::ndim<1>>` wrapper in `ga_py/src/module.cpp` that builds an `std::mdspan` over the numpy storage and calls the existing `rk4_step` — zero-copy, real performance, but the first hand-written binding in the project.

### 6.6 Demo / teaching scripts

`ga_py/demo/hello_ga.py` (the §3.4 install-check) is the first one. Topical follow-ups along the same lines are pending: `demo/ega_basics.py`, `demo/pga_basics.py`, `demo/rotors_and_motors.py`, etc. — minimal scripts demonstrating the `import ga_py` workflow on focused themes.

### 6.7 Documentation site

Sphinx + autodoc was deliberately deferred until docstrings exist (also not yet started). The natural ordering when picked up: docstrings → autodoc → Sphinx site.

### 6.8 Excluded namespace `hd::ga::detail`

The scanner has an explicit allow-list (`TARGET_NAMESPACES = {"hd::ga", "hd::ga::ega", "hd::ga::pga"}` in `ga_bindgen/src/scan.py`). Anything in `hd::ga::detail` is filtered out at scan time and never reaches Python. This is intentional — `detail` contains tolerance helpers (`safe_epsilon`), runtime guards (`check_division_by_zero`, `check_normalization`, `check_unitization`) that the public API already calls, the `extended_testing_enabled()` compile-time flag, and `fmt::formatter` specializations reached transitively through `__str__`. None are part of the public C++ API and none belong in the Python surface.

If a `detail` helper later turns out to be useful from Python, the right move is **not** to relax the scanner's allow-list — it's to **promote that specific symbol out of `detail`** in the C++ library. `hd::ga::detail` is the universally-recognised C++ idiom for "internal; subject to change without notice"; breaking that convention to suit a binding generator would muddy the public surface for C++ users too.
