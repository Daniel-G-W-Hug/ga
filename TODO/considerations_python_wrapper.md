# Python Wrapper (`ga_py/`) — Considerations & Plan

**Status (as of 2026-04-26):** v1 is implemented. All 46 user types are bound, all five testing tiers are in place, and the suite (411 tests) passes against the C++ reference library in ~1.5 s.

This document keeps the original design rationale (sections 1–10b) plus a short index of what is **still open** for future work — see [§11](#11-outstanding-work). The detailed deferred-bindings analysis is in [§12](#12-future-work--deferred-bindings) (currently just `rk4_step`).

---

## 1. Goal

Provide a Python wrapper around the `ga/` header-only C++ library so that Python users can access the full GA functionality without touching C++. The wrapper:

- is imported as `import ga_py`
- exposes submodules `ga_py.ega` and `ga_py.pga` mirroring the `hd::ga::ega` / `hd::ga::pga` C++ namespaces
- covers *only* the contents of `ga/` (not `ga_view/`, `ga_lua/`, `ga_prdxpr/`)
- is generated from the C++ headers, not hand-authored, so additions in C++ propagate to Python by rerunning a code-generation step

The C++ implementation remains the lead. The Python wrapper tracks it.

---

## 2. Primary user

A Python developer who wants to use the GA library **instead of** the C++ interface. Translating C++ documentation/examples 1:1 is the dominant workflow. Naming, types, and operator behavior should match C++ as closely as idiomatic Python allows.

Secondary concerns (scientific computing, teaching) are acknowledged but not driving decisions.

---

## 3. Decisions locked in

| Area | Decision |
| --- | --- |
| Binding technology | **nanobind** (modern, header-only, C++20/23-friendly, ~½ the binding line-count of pybind11) |
| Type model | Mirror C++ 1:1 — each concrete user type (`vec2d`, `bivec3d`, `mvec2dp_e`, …) is a distinct Python class |
| Submodule layout | `ga_py.ega` and `ga_py.pga` mirroring C++ namespaces; convenience re-exports at `ga_py` top level |
| Value type | `double` only (v1). Float support deferred — revisit if user demand justifies it. |
| numpy integration | **Tier 1 only**: buffer protocol / `__array__` on types; constructors accept `np.ndarray`. No vectorized array-of-multivector ops. |
| Operators in Python | Keep the same symbols as C++ (`*`, `^`, `<<`, `>>`, unary `-`, `~`). Document the Python precedence gotcha (see §9). |
| Constants (`ga_usr_consts.hpp`) | Exposed in the most convenient form: e.g. `ga_py.ega.e1`, `ga_py.ega.e12`, … as module attributes |
| Error handling | Custom `ga_py.GaError(RuntimeError)`; nanobind exception translator maps `hd::ga` exceptions |
| Missing-operation handling (v1) | When a Python call reaches an unbound combination: raise `GaError("Operation not yet available in Python wrapper — requires C++ extension")`. No request-tracking infrastructure. |
| Printing | `__repr__` = eval-able form (`vec3d(1.0, 2.0, 3.0)`); `__str__` forwards to existing `fmt` formatter (`Vec3d(1,2,3)`). Jupyter `_repr_latex_` is future work. |
| Python version | 3.10+ only (modern typing, covers real users; 3.9 is EOL) |
| Docstrings | None in v1 — separate effort later |
| Demo/teaching scripts | Minimal, just enough to demonstrate the `import ga_py` workflow |
| Platforms | macOS + Linux + Windows — all mandatory |
| Python-version matrix (v1) | **Build for Python 3.10 on each platform only.** Skip the 3.10/3.11/3.12/3.13 matrix for now — one wheel per platform. Can add more Python versions later without design changes. |
| Testing strategy | Cross-check equivalence against C++ reference cases + property-based tests (see §7). Do **not** reproduce `ga_test/` scope. |
| Source of truth for bindings | **C++ headers, parsed with libclang**. Not ga_lua.hpp (incomplete). Not hand-authored manifests (drift). |
| Generation model | Generated binding source committed to git; regeneration is an explicit step (CMake target) run by contributors, not an implicit build dependency |
| Tool location | **`ga_bindgen/`** as a separate top-level directory (rationale in §5) |

---

## 4. Architecture overview

Actual layout as implemented (some originally-planned items deferred — see [§11](#11-outstanding-work)):

```text
ga/                       (header-only C++23 library, source of truth)
ga_bindgen/               (scanner + binding generator, Python tool)
    src/
        scan.py           libclang → manifest.json
        emit_nanobind.py  manifest.json → ga_py/src/generated/*.cpp
        clang_setup.py    libclang probe (homebrew vs pip-bundled)
        model.py          dataclasses for the manifest schema
    manifest.json         committed snapshot of last-known API
ga_py/
    CMakeLists.txt        nanobind + ga/ include
    src/
        module.cpp        hand-written NB_MODULE entrypoint
        generated/        produced by ga_bindgen — do not edit
            bindings_<type>.cpp        one per bound user type (46 today)
            bindings_constants_*.cpp   per submodule
            bindings_functions_*.cpp   per submodule
            register_all.cpp
            bindings_list.cmake
    python/ga_py/__init__.py   re-exports submodules
    tests/                pytest, see §7
        data/ga_test_cases.json   T5 cross-check snapshot

    [DEFERRED] pyproject.toml (scikit-build-core)  — see §11
    [DEFERRED] stubs/ (.pyi files)                 — see §11
    [DEFERRED] examples/                           — see §11
```

### Dependency chain

```text
ga/*.hpp  ──(libclang parse)──►  ga_bindgen/manifest.json
                                        │
                                        ├──(emit)──►  ga_py/src/generated/*.cpp
                                        │
                                        ├──(emit)──►  ga_py/stubs/*.pyi    [DEFERRED]
                                        └──(audit)──►  ga_lua coverage      [DEFERRED]
```

---

## 5. Binding generator — `ga_bindgen/`

### 5.1 Why a separate directory (not under `ga_prdxpr/`)

`ga_prdxpr` generates GA product *expressions* (strings of C++ code) and is architecturally specialized for that domain. `ga_bindgen` parses C++ headers with libclang and emits nanobind glue. Different input, different output, different toolchain. Coupling them would muddle both.

If `ga_prdxpr` later grows into full C++ class generation (Project A from the design discussion), the two remain independent tools that can share config data if useful.

### 5.2 Role

Implemented:

1. Parse `ga/*.hpp` with libclang, enumerating class templates, free functions, operators, their signatures, template parameters, namespaces.
2. Produce an API manifest (`manifest.json`) — a canonical snapshot of every binding-relevant symbol.
3. Emit nanobind binding source (`bindings_*.cpp`) from the manifest.

Deferred (see [§11](#11-outstanding-work)):

1. Emit type stub files (`*.pyi`) from the manifest.
2. Produce a delta report against the previous manifest — "new since last run", "removed since last run", "changed signature".
3. Audit `ga_lua/src/ga_lua.hpp` against the manifest. (The 2026-04 spike confirmed gaps via crude string search; a per-overload audit is the next step.)

### 5.3 Structure

```text
ga_bindgen/
    README.md
    src/
        scan.py            libclang → manifest.json
        emit_nanobind.py   manifest.json → ga_py/src/generated/
        clang_setup.py     find_libclang() + compile flags
        model.py           dataclasses: TypeAlias, Constructor, Function, Overload, ...
    manifest.json          committed — last-known API snapshot

    [DEFERRED] emit_stubs.py
    [DEFERRED] audit_lua.py
    [DEFERRED] reports/
```

### 5.4 Dependencies

- Python 3.10+
- `libclang` Python bindings (`clang` on PyPI, plus a system/homebrew libclang the runtime probes for via `clang_setup.find_libclang()`)
- f-strings used directly (no jinja2 templates needed in practice)
- No runtime dependencies in the generated code beyond nanobind + the `ga/` headers

### 5.5 Regeneration workflow

```bash
# 1) Re-scan headers (only when ga/*.hpp changed)
python3 ga_bindgen/src/scan.py

# 2) Re-emit bindings
python3 ga_bindgen/src/emit_nanobind.py --all

# 3) Recompile the Python module
cmake --build build/ga_py
```

**Generated files are committed to git.** This means:

- PR reviewers see exactly what new Python surface a C++ change created.
- End users building `ga_py/` don't need libclang installed.
- CI doesn't need libclang on the wheel-building runners.

---

## 6. Handling ga_lua completeness (side benefit) — DEFERRED

`ga_lua/src/ga_lua.hpp` is hand-maintained and known incomplete. The plan is to invert the relationship — let `ga_bindgen` produce the authoritative manifest and audit Lua against it.

Status: a 2026-04 string-search spike confirmed gaps (`rk4_step`, `rk4_get_time`, `Hz2radps`, `radps2Hz`, `rpm2radps`, `radps2rpm`, `*_metric_view` helpers). A per-overload `audit_lua.py` driven by the manifest is the next step (see [§11](#11-outstanding-work)). The Python wrapper does not depend on it.

---

## 7. Testing strategy

### 7.1 Principles

- Re-running the full `ga_test/` suite in Python is duplication for no gain. Instead, sample representatively, randomize exhaustively (`hypothesis`), and cross-check actively against the C++ reference output.
- Internal `hd::ga::detail::*` helpers are not tested (correctly excluded per §10b).
- The binding generator itself (`ga_bindgen`) is not tested as part of this suite; it is a separate tool with its own concerns.

### 7.2 Implemented tiers — 411 tests, ~1.5 s total

| Tier | Goal | Tests | File |
| --- | --- | --- | --- |
| Smoke (T0 baseline) | Each bound type constructs / prints / compares / dispatches operators | 114 | `test_vec3d.py`, `test_coverage.py`, `test_scalars.py`, `test_multivectors.py`, `test_pga_primitives.py` |
| T1 Constants | Every bound constant has the value declared in `ga_usr_consts.hpp` (107 + alias relations + drift detection) | 144 | `test_constants.py` |
| T2 Grade lookup | `gr()` per graded type, `rgr()` per PGA graded type | 30 | `test_grade_lookup.py` |
| T3 EGA identities | Algebraic invariants randomized via `hypothesis` (dot symmetry/linearity, wedge anti-symmetry, gpr split, inv, rev involution, rotor unitarity) | 15 | `test_identities_ega.py` |
| T4 PGA identities | Bulk/weight decomposition, complement involution, unitize, motor sandwich, distance, join/meet | 20 | `test_identities_pga.py` |
| T5 Cross-check | Replays a JSON snapshot of (op, args, expected) tuples emitted by `ga_test/src/export_python_cases.cpp` and verifies bit-for-bit (within tolerance) equivalence with C++ | 83 | `test_cross_check.py`, `data/ga_test_cases.json` |

`conftest.py` carries the relative-tolerance helpers (`approx_eq`, `EPS_TIGHT`, `EPS_DEFAULT`, `EPS_LOOSE`) and the `components()` extractor used by T1 and T5.

### 7.3 Deferred

- **T6 — `mypy --strict` on the example scripts.** Held until `ga_bindgen` emits `.pyi` stubs (see [§11](#11-outstanding-work)). The test file itself (`test_stubs_mypy.py`) is one subprocess call; the precondition is the stubs.

---

## 8. Constants exposure

`ga/ga_usr_consts.hpp` holds algebra-specific constants (basis vectors, identities, etc.). Expose them as module attributes to match the obvious Python idiom:

```python
import ga_py
v = 2 * ga_py.ega.e1 + 3 * ga_py.ega.e2
```

Attributes live directly on `ga_py.ega` / `ga_py.pga` submodules. Generated from the manifest so new constants in `ga_usr_consts.hpp` appear automatically.

---

## 9. Operator mapping and precedence gotcha

| C++ | Python dunder | Python symbol |
| --- | --- | --- |
| `+` | `__add__` | `+` |
| `-` (binary) | `__sub__` | `-` |
| `-` (unary) | `__neg__` | `-` |
| `*` (gpr) | `__mul__` | `*` |
| `^` (wdg) | `__xor__` | `^` |
| `<<` (lcontract) | `__lshift__` | `<<` |
| `>>` (rcontract) | `__rshift__` | `>>` |
| `~` (TBD, if used) | `__invert__` | `~` |

**Precedence gotcha — document prominently in Python docs:**

In Python, `^` (wedge) has **lower** precedence than `+` and `*`, unlike C++ where it's already low but the parenthesization habits differ. Concretely:

```python
a + b ^ c       # Python: parses as (a + b) ^ c      — DIFFERENT FROM C++ HABIT
a * b ^ c       # Python: parses as (a * b) ^ c
a << b ^ c      # Python: parses as (a << b) ^ c
```

Users should parenthesize wedge products defensively. This is not a bug — it's Python's operator table, and it's consistent for all Python programs using `^`. Surface it in docs and examples.

---

## 10. Printing / formatting

- `__repr__` → eval-able form, e.g. `vec3d(1.0, 2.0, 3.0)` — unambiguous, round-trippable, shown in interactive REPL.
- `__str__` → forward directly to the existing `fmt` formatters in `ga/detail/ga_fmt_support.hpp`, producing output like `Vec3d(1,2,3)`. Reuse, don't re-author.
- `__format__(spec)` → support `f"{v:.3f}"` if cheap; otherwise defer.
- `_repr_latex_` (Jupyter math rendering) → deferred to v1.1.

---

## 10b. Excluded namespaces — `hd::ga::detail`

The scanner has an explicit allow-list (`TARGET_NAMESPACES = {"hd::ga", "hd::ga::ega", "hd::ga::pga"}` in `ga_bindgen/src/scan.py`). Anything in `hd::ga::detail` is filtered out at scan time and never reaches Python. This section records what's in there so we don't lose track of it as the library evolves.

### What's in `hd::ga::detail` today

**No types** (zero `class` / `struct` definitions). Only free functions and function templates serving as implementation helpers:

| Symbol | File | Purpose |
| --- | --- | --- |
| `safe_epsilon<T>()`, `safe_epsilon<T, U>()` | `ga/detail/ga_error_handling.hpp` | Picks the right floating-point tolerance for the value type(s) |
| `check_division_by_zero<T>(divisor, op_name)` | `ga/detail/ga_error_handling.hpp` | Runtime guard — throws if denominator is too small |
| `check_normalization<T>(magnitude, type_name)` | `ga/detail/ga_error_handling.hpp` | Throws if a `normalize()` would divide by ~0 |
| `check_unitization<T>(weight_norm, type_name)` | `ga/detail/ga_error_handling.hpp` | Throws if a `unitize()` would divide by ~0 |
| `extended_testing_enabled()` | `ga/detail/ga_error_handling.hpp` | Compile-time flag for `_HD_GA_EXTENDED_TEST_DIV_BY_ZERO` |
| `suppress_negative_zero(value)` | `ga/detail/fmt/` | Formatting helper — turns `-0.0` into `0.0` for display |
| (various `fmt::formatter<...>` specializations) | `ga/detail/fmt/` | Hooks GA types into the `fmt` library |

### Why these are not bound

- `safe_epsilon`, `extended_testing_enabled`, `suppress_negative_zero` — used internally to enforce invariants and produce clean output. The public API already calls them where needed. Python users would never invoke them directly.
- `check_division_by_zero` / `check_normalization` / `check_unitization` — these throw `std::invalid_argument` from inside C++. Nanobind's exception translator surfaces those as Python `RuntimeError` (or our planned `GaError`). End users see clean error messages without ever needing to call these guards from Python.
- The `fmt::formatter` specializations are reached transitively through bound `__str__` methods that call `fmt::format("{}", v)`.

### When the exclusion would be wrong

If a `detail` helper turns out to be useful from Python (e.g. exposing the configured epsilon for debugging, or a `set_epsilon()` mutator), the right move is **not** to relax the scanner's allow-list. It is to **promote that specific symbol out of `detail`** in the C++ library — `hd::ga::detail` is the universally-recognised C++ idiom for "internal; not part of the API; subject to change without notice", and breaking that convention to suit a binding generator would muddy the public surface for C++ users too.

### Audit cadence

Whenever new helpers land in `hd::ga::detail`, this table should be reviewed. The Lua-coverage audit tool (`ga_bindgen/audit_lua.py` — planned, see §5.2) will, by symmetric reasoning, also surface anything new in `detail` that callers in `ga/` start to leak through their public signatures.

---

## 11. Outstanding work

This section enumerates everything still open. Implemented work has been culled from the rest of the document; the items below are what remains.

### 11.1 Packaging (not started)

The wrapper is currently built via direct CMake (see `ga_py/README.md`). For end-user distribution as wheels:

- Add `pyproject.toml` with `scikit-build-core` as the build backend; declare `nanobind` (build-only) and `numpy` (runtime, for buffer protocol) as dependencies.
- Set up `cibuildwheel` on GitHub Actions for Python 3.10 × {macOS Intel + ARM, Linux x86_64, Windows x86_64} = 4 wheels.
- Publish source distribution (sdist) too.
- Additional Python versions (3.11/3.12/3.13) extend the cibuildwheel matrix without design changes.
- nanobind, fmt, and `ga/` are all header-only — no vcpkg/brew dependencies on the wheel-building side, so CI stays simple.

### 11.2 Stub generation `.pyi` (not started)

`ga_bindgen/src/emit_stubs.py` does not yet exist. Once added:

- Drives the IDE/type-check experience for users.
- Unblocks T6 (mypy --strict on the example scripts) — see [§7.3](#73-deferred).

### 11.3 numpy buffer protocol — Tier 1 (not started)

Per the §3 decision, each bound user type should expose `__array__` / the buffer protocol so numpy arrays of GA elements can be constructed cheaply. Tier 2 (vectorized array-of-multivector ops) remains deferred indefinitely.

### 11.4 `audit_lua.py` (spike done, full version pending)

The 2026-04 string-search spike confirmed Lua-wrapper gaps; a per-overload audit driven by the manifest is the next step. Output target: `ga_bindgen/reports/lua_coverage_*.md`.

### 11.5 `rk4_step` deferral cleanup

`rk4_get_time` is bound but useless without `rk4_step` (which cannot be auto-bound today — see [§12](#12-future-work--deferred-bindings)). The decided path is **defer permanently** and add a `scipy.integrate.solve_ivp` recipe to `ga_py/README.md`.

Action items:

- [ ] Add `rk4_get_time` to `SKIP_FREE_FN_NAMES` in `ga_bindgen/src/emit_nanobind.py` (or remove the binding via a more targeted filter).
- [ ] Regenerate bindings, verify the function is gone.
- [ ] Add a README section "Integrating GA-valued ODEs in Python" with a `solve_ivp` example using a `vec3d`-valued state vector.

### 11.6 Demo / teaching scripts (deferred per §3)

`ga_py/examples/ega_basics.py`, `ga_py/examples/pga_basics.py`, etc. — minimal scripts demonstrating the `import ga_py` workflow.

### 11.7 Documentation site

Sphinx + autodoc was deliberately deferred until docstrings exist (also deferred per §3). When picked up, the natural ordering is: docstrings → autodoc → Sphinx site.

### 11.8 Float `value_t` support

Currently locked to `double`. Float (and possibly `long double`) support would require a build-flag-driven manifest instantiation, doubling the binding surface. Revisit if user demand justifies it.

### 11.9 Manifest delta report

`ga_bindgen` does not currently produce `reports/delta_YYYY-MM-DD.md`. Useful for PR reviewers to see "new since last run", "removed since last run", "changed signature" diffs. Low priority — `git diff` on `manifest.json` mostly serves the same purpose.

---

## 12. Future work — deferred bindings

Items deliberately not bound in v1, with the reasoning preserved here so the next maintainer doesn't re-litigate the decision (or re-discover the obstacle from scratch).

### 12.1 `rk4_step` (Runge–Kutta 4th-order time integrator)

**Status:** the C++ free function `hd::ga::rk4_step<VecType>(...)` (declared in `ga/ga_usr_utilities.hpp`) is **not bound** in `ga_py`. The companion helper `rk4_get_time` *is* bound but is useless on its own (it just returns the time at sub-step 0/1/2/3 — meaningful only when paired with the integrator).

**Why it can't be auto-bound today.** Three obstacles, in increasing order of difficulty:

1. **`std::mdspan` has no nanobind type-caster.** All three array arguments (`u`, `uh`, `rhs`) are `std::mdspan<VecType, std::dextents<size_t, N>>`. nanobind ships casters for `std::vector`, `std::array`, numpy arrays via the buffer protocol, etc., but not for C++23 `std::mdspan` (it's brand new and a non-owning view that the binding layer would have to materialise over caller-supplied storage). The bindgen sees the parameter, fails `resolve_param_to_user_type`, and silently skips the function.
2. **Mutates its arguments instead of returning a value.** `u` and `uh` are written through their mdspans — that's idiomatic for tight numerical loops but doesn't fit the Python convention of returning new state from a function call.
3. **Templated on `VecType`.** Each GA vector type (`vec2d`, `vec3d`, `vec2dp`, `vec3dp`, …) gives a different concrete instantiation that has to be enumerated and bound separately. This is the *least* hard of the three (the same machinery already handles `gpr<T>`, `wdg<T>`, etc.).

**Do we even need to bind it?** Probably not. Python users who need to integrate ODEs over GA-valued state already have superior tools:

- **`scipy.integrate.solve_ivp`** — the modern unified ODE solver (scipy ≥ 1.0). Default method is `RK45` (Dormand–Prince adaptive 5th order with embedded 4th-order error estimator). Other built-in methods: `RK23`, `DOP853`, `Radau`, `BDF`, `LSODA`. Adaptive step size, dense output, event detection — all features `rk4_step` does not have.
- **`scipy.integrate.RK45`** as a class (lower-level access if you want to drive sub-stepping yourself).
- **Fixed-step RK4** — if a user really wants a non-adaptive RK4 specifically, it is a 6-line numpy function; no library needed.
- **`numba` / `jax`** — JIT-compiled integrators if performance matters more than scipy can deliver.

The natural Python workflow is:

1. Convert GA-valued state to a flat numpy array (each `vec3d` becomes 3 doubles, etc.).
2. Wrap the GA-valued `rhs(t, u)` in a closure that unpacks/repacks numpy ↔ GA objects via the bound types.
3. Call `solve_ivp(rhs, t_span, u0, method="RK45")`.
4. Reconstruct GA objects from the result.

The packing/unpacking is per-component and can be batched once we add the Tier-1 numpy buffer-protocol work mentioned in §3 (where each GA type exposes `__array__` / buffer protocol). At that point the integration loop is essentially zero-overhead from the GA wrapper's side.

**Recommendation.** Three coherent options, in order of preference:

1. **Defer permanently and document.** Remove `rk4_get_time` from the bindings (it's useless without `rk4_step`), and add a recipe to `ga_py/README.md` showing how to use `scipy.integrate.solve_ivp` on a GA-valued ODE. This is the cheapest *and* gives Python users the better tool. Cost: ~1 hour (delete + add `SKIP_FREE_FN_NAMES` entry + write README example).
2. **Provide a thin `std::vector<VecType>` overload in C++** so the bindgen picks it up automatically:

   ```cpp
   template <typename VecType>
   inline std::vector<VecType> rk4_step(
       std::vector<VecType> const& u,
       std::vector<VecType> & uh_buf,
       std::vector<VecType> const& rhs,
       value_t dt, size_t rk_step);
   ```

   `<nanobind/stl/vector.h>` is already included; bindgen would emit per-`VecType` overloads automatically. Cost: ~10 lines of C++ wrapper. Downside: per-call heap allocation, slower than scipy in hot loops, slower than option 3 in any loop.
3. **Hand-write a numpy-buffer-protocol binding** in `ga_py/src/module.cpp` (the only hand-written file). Accept `nb::ndarray<vec3d, nb::ndim<1>>`, build an `std::mdspan` over `arr.data()` and `arr.shape(0)`, call the existing `rk4_step` directly. Zero-copy, real performance. Downside: hand-written code that will diverge from C++ over time; has to be repeated per `VecType`; introduces the *first* hand-written binding in the project.

Given that scipy already provides better integrators, **option 1 is the recommended path.** Keep the C++ implementation (it's still useful for native C++ users and for the `ga_view`/`ga_lua` consumers), but stop exposing the orphaned `rk4_get_time` in Python and point users to scipy.

Action items for this option are tracked under [§11.5](#115-rk4_step-deferral-cleanup). If demand later justifies a native binding, prefer option 3 over option 2 (zero-copy beats heap-thrashing for this kind of inner loop).

---

## 13. Appendix — example of generated nanobind code (sketch)

For a single type `vec3d` with three operators, the generator would emit roughly:

```cpp
// AUTO-GENERATED BY ga_bindgen — DO NOT EDIT
// Source: ga/ga_ega3d_ops.hpp, ga/ga_usr_types.hpp
// Regenerate via: cmake --build . --target regenerate_py_bindings

#include <nanobind/nanobind.h>
#include <nanobind/operators.h>
#include "ga/ga_ega.hpp"

namespace nb = nanobind;
using namespace hd::ga;
using namespace hd::ga::ega;

void bind_vec3d(nb::module_& m) {
    nb::class_<vec3d>(m, "vec3d")
        .def(nb::init<>())
        .def(nb::init<value_t, value_t, value_t>(), "x"_a, "y"_a, "z"_a)
        .def_rw("x", &vec3d::x)
        .def_rw("y", &vec3d::y)
        .def_rw("z", &vec3d::z)
        .def("__repr__", [](vec3d const& v) {
            return fmt::format("vec3d({}, {}, {})", v.x, v.y, v.z);
        })
        .def("__str__", [](vec3d const& v) { return fmt::format("{}", v); })
        .def(-nb::self)
        .def(nb::self + nb::self)
        .def(nb::self - nb::self)
        // geometric product vec3d * vec3d (if defined):
        // .def("__mul__", [](vec3d const& a, vec3d const& b) { return gpr(a, b); })
        // wedge:
        .def("__xor__", [](vec3d const& a, vec3d const& b) { return wdg(a, b); })
        // ... all overloads enumerated by the generator from the manifest
        ;
}
```

The generator's job is to produce that code, for every type, every overload, every algebra, from the libclang-derived manifest — and to keep it in sync as `ga/` evolves.
