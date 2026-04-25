# Python Wrapper (`ga_py/`) — Considerations & Plan

**Status:** design document. Not yet implemented.
**Prerequisite:** header scanner + binding generator (see `ga_bindgen/` section).
**Author context:** result of interactive design discussion, intended for further refinement and as implementation guidance.

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

```
ga/                       (unchanged — header-only C++23 library, source of truth)
ga_bindgen/               (NEW — scanner + binding generator, Python tool)
ga_py/                    (NEW — the Python module; most .cpp files generated)
    CMakeLists.txt
    pyproject.toml        (scikit-build-core)
    src/
        module.cpp         (hand-written entrypoint: NB_MODULE + dispatch)
        generated/
            bindings_ega.cpp       (generated from ga/ga_ega*.hpp)
            bindings_pga.cpp       (generated from ga/ga_pga*.hpp)
            bindings_common.cpp    (generated: value_t, usr_consts, utilities)
    stubs/
        ga_py/__init__.pyi         (generated)
        ga_py/ega.pyi              (generated)
        ga_py/pga.pyi              (generated)
    examples/
        ega_basics.py
        pga_basics.py
    tests/
        test_equivalence.py    (cross-check vs C++ reference cases)
        test_identities.py     (hypothesis property tests)
        test_smoke.py
```

### Dependency chain

```
ga/*.hpp  ──(libclang parse)──►  ga_bindgen/manifest.json
                                        │
                                        ├──(emit)──►  ga_py/src/generated/*.cpp
                                        ├──(emit)──►  ga_py/stubs/*.pyi
                                        └──(audit)──►  reports vs ga_lua.hpp
```

---

## 5. Binding generator — `ga_bindgen/`

### 5.1 Why a separate directory (not under `ga_prdxpr/`)

`ga_prdxpr` generates GA product *expressions* (strings of C++ code) and is architecturally specialized for that domain. `ga_bindgen` parses C++ headers with libclang and emits nanobind glue. Different input, different output, different toolchain. Coupling them would muddle both.

If `ga_prdxpr` later grows into full C++ class generation (Project A from the design discussion), the two remain independent tools that can share config data if useful.

### 5.2 Role

1. Parse `ga/*.hpp` with libclang, enumerating: class templates, free functions, operators, their signatures, template parameters, namespaces.
2. Produce an **API manifest** (`manifest.json`) — a canonical snapshot of every binding-relevant symbol.
3. Emit nanobind binding source (`bindings_*.cpp`) from the manifest.
4. Emit type stub files (`*.pyi`) from the manifest.
5. Produce a **delta report** against the previous manifest — "new since last run," "removed since last run," "changed signature."
6. **Bonus:** audit `ga_lua/src/ga_lua.hpp` against the manifest and produce a coverage report — what's in `ga/` but not bound in Lua. (This is the completeness check the user asked for; same tool, double use.)

### 5.3 Structure (sketch)

```
ga_bindgen/
    README.md
    pyproject.toml         (Python tool, not a wheel)
    src/
        scan.py            (libclang → manifest.json)
        emit_nanobind.py   (manifest.json → ga_py/src/generated/)
        emit_stubs.py      (manifest.json → ga_py/stubs/)
        audit_lua.py       (manifest.json + ga_lua.hpp → coverage report)
        model.py           (dataclasses: Type, Function, Overload, etc.)
        templates/         (jinja2 templates for emitted code)
    manifest.json          (committed — last-known API snapshot)
    reports/
        delta_YYYY-MM-DD.md
        lua_coverage_YYYY-MM-DD.md
```

### 5.4 Dependencies

- Python 3.10+
- `libclang` Python bindings (`clang` on PyPI, or system libclang via `clang.cindex`)
- `jinja2` for templates (optional — could use f-strings)
- No runtime dependencies in the generated code beyond nanobind + the ga/ headers

### 5.5 Regeneration workflow

```bash
cd build
cmake --build . --target regenerate_py_bindings
# Runs: python ga_bindgen/src/scan.py  → manifest.json
#       python ga_bindgen/src/emit_nanobind.py
#       python ga_bindgen/src/emit_stubs.py
# Diff appears in: ga_py/src/generated/*.cpp, ga_py/stubs/*.pyi, ga_bindgen/manifest.json
# Contributor reviews diff, commits alongside their C++ change
```

**Generated files are committed to git.** This means:
- PR reviewers see exactly what new Python surface a C++ change created.
- End users building `ga_py/` don't need libclang installed.
- CI doesn't need libclang on the wheel-building runners.

---

## 6. Handling ga_lua completeness (side benefit)

`ga_lua/src/ga_lua.hpp` was hand-maintained and is known incomplete. Rather than mine it for the Python manifest (which would propagate its gaps), we invert the relationship:

- `ga_bindgen` parses C++ headers → authoritative manifest.
- `ga_bindgen/audit_lua.py` compares manifest against ga_lua.hpp.
- Produces `reports/lua_coverage_*.md` listing symbols present in `ga/` but missing in ga_lua.
- This becomes the to-do list for retrofitting Lua completeness (if desired).

No code changes to `ga_lua/` required as part of this project — the audit is read-only.

---

## 7. Testing strategy

### 7.1 What we don't do

- Re-run the full `ga_test/` suite in Python. That's duplication for no gain.

### 7.2 What we do

**Tier 1 — Cross-check equivalence (primary)**
- Export a curated set of (input, expected output) pairs from `ga_test/` as JSON.
- Python test harness loads JSON, performs the same operation via `ga_py`, asserts numerical equality within tolerance.
- Catches: binding dispatch errors, sign flips, wrong overload selection, type coercion bugs.

**Tier 2 — Property-based testing (`hypothesis`)**
- Algebraic identities that must hold for arbitrary inputs, e.g.:
  - `(a·b)² + (a∧b)² ≈ |a⟑b|²` for vectors
  - `a ⟑ rev(a) / nrm_sq(a) ≈ 1` for non-null multivectors
  - `l_cmpl(r_cmpl(a)) ≈ a` (complement involution in appropriate algebras)
  - Motor sandwich preserves norms in PGA
- Hypothesis generates random inputs within specified domains; failures surface minimal counterexamples.

**Tier 3 — Smoke tests**
- Every bound type: can be constructed, printed, equality-tested, combined with basic operators.
- Runs fast; first line of defense against binding regressions.

**Tier 4 — Stub validation**
- `mypy --strict` passes on the example scripts.
- Protects the IDE/type-check experience.

### 7.3 Implementation plan for testing (proposal for a new session)

Status as of writing: 114 smoke / coverage / scalar / multivector / pga-primitive tests in `ga_py/tests/` already pass. They confirm bindings load, types construct, fields access, basic operators dispatch correctly. **Tier 3 (smoke) is essentially done.** What follows is the proposed plan for Tiers 1–2 and 4, in priority order.

#### Test file layout (create / extend in `ga_py/tests/`)

```
ga_py/tests/
├── conftest.py                 # NEW: shared fixtures, EPS, common imports
├── test_vec3d.py               # existing (13)
├── test_coverage.py            # existing (35)
├── test_scalars.py             # existing (32 — incl. comparison ops)
├── test_multivectors.py        # existing (12)
├── test_pga_primitives.py      # existing (22)
│
├── test_identities_ega.py      # NEW: Tier-2 algebraic identities for EGA
├── test_identities_pga.py      # NEW: Tier-2 algebraic identities for PGA
├── test_cross_check.py         # NEW: Tier-1, loads ga_test reference cases
├── data/                       # NEW
│   └── ga_test_cases.json      # NEW: exported reference cases from ga_test
├── test_constants.py           # NEW: verify the 107 constants have known values
└── test_grade_lookup.py        # NEW: verify gr() / rgr() give correct grades
```

The first session in the new context should run the existing 114 tests as the starting baseline, then add tests in the priority order below.

#### Priority order — small, valuable steps first

**Step T1 — Constants verification (~30 min)**
File: `test_constants.py`. Cheap, high-value: verifies all 107 bound constants have the values the C++ source declares.

```python
def test_basis_2d_identity():
    assert ga_py.ega.e1_2d.x == 1.0 and ga_py.ega.e1_2d.y == 0.0
    assert (ga_py.ega.e1_2d ^ ga_py.ega.e2_2d) == ga_py.ega.I_2d

def test_basis_3d_orthonormality():
    for i, ei in enumerate([ga_py.ega.e1_3d, ga_py.ega.e2_3d, ga_py.ega.e3_3d]):
        for j, ej in enumerate([ga_py.ega.e1_3d, ga_py.ega.e2_3d, ga_py.ega.e3_3d]):
            assert float(ga_py.ega.dot(ei, ej)) == (1.0 if i == j else 0.0)
    # e1 × e2 = e3
    assert ga_py.ega.cross(ga_py.ega.e1_3d, ga_py.ega.e2_3d) == ga_py.ega.e3_3d
```

Catches: any silent corruption in the constants emission, regression in operator dispatch on basis values.

**Step T2 — Grade lookup verification (~15 min)**
File: `test_grade_lookup.py`. Tiny but completes the type-safety story.

```python
@pytest.mark.parametrize("typed_obj,expected_grade", [
    (ga_py.ega.scalar2d(1.0),    0),  (ga_py.ega.vec2d(1,2),       1),
    (ga_py.ega.pscalar2d(1.0),   2),
    (ga_py.ega.scalar3d(1.0),    0),  (ga_py.ega.vec3d(1,2,3),     1),
    (ga_py.ega.bivec3d(1,2,3),   2),  (ga_py.ega.pscalar3d(1.0),   3),
    (ga_py.pga.scalar3dp(1.0),   0),  (ga_py.pga.trivec3dp(1,2,3,4), 3),
    (ga_py.pga.pscalar3dp(1.0),  4),
])
def test_gr_returns_correct_grade(typed_obj, expected_grade):
    assert ga_py.gr(typed_obj) == expected_grade
```

**Step T3 — EGA algebraic identities via hypothesis (~half day)**
File: `test_identities_ega.py`. Add `pip install hypothesis` to the venv. Write strategy generators for `vec3d`, `bivec3d`, `mvec3d`, then test identities that must hold for any input.

Concrete identity examples:

```python
from hypothesis import given, strategies as st
finite = st.floats(min_value=-100, max_value=100, allow_nan=False, allow_infinity=False)

@st.composite
def vec3d_st(draw):
    return ga_py.ega.vec3d(draw(finite), draw(finite), draw(finite))

@given(vec3d_st(), vec3d_st())
def test_dot_wedge_split(a, b):
    """For Euclidean vectors: |dot(a,b)|² + |wdg(a,b)|² == |a|²·|b|²"""
    dot_sq = float(ga_py.ega.dot(a, b)) ** 2
    wedge_sq = float(ga_py.ega.nrm_sq(ga_py.ega.wdg(a, b)))
    expected = float(ga_py.ega.nrm_sq(a)) * float(ga_py.ega.nrm_sq(b))
    assert abs(dot_sq + wedge_sq - expected) < 1e-9 * max(abs(expected), 1.0)

@given(vec3d_st())
def test_double_complement_is_identity(a):
    """l_cmpl(r_cmpl(a)) == a in EGA — wait, EGA3D uses single cmpl. Adjust per-algebra."""
    # in 3D EGA: cmpl(cmpl(a)) is sign-dependent on grade
    pass

@given(vec3d_st(), vec3d_st(), vec3d_st())
def test_wdg_associative(a, b, c):
    """Wedge product is associative: (a^b)^c == a^(b^c)"""
    lhs = ga_py.ega.wdg(ga_py.ega.wdg(a, b), c)
    # mvec3d_e wedge with vec3d returns mvec3d_? — check shape; this may need API tweak
    pass
```

Identity inventory to start with (each one becomes a test):

- `(a·b)² + (a∧b)² == |a|²·|b|²` for vec×vec
- `a · b == b · a` (dot symmetry)
- `a ∧ b == −(b ∧ a)` (wedge anti-symmetry)
- `gpr(a, a) == |a|²` for vectors (`a²` = scalar = squared magnitude)
- `inv(a) ⟑ a == 1` for non-null
- `rev(rev(a)) == a` for any grade
- `rotate(a, R) preserves nrm_sq(a)` when R is a unit rotor
- `(R ⟑ R^~) == 1` for unit rotors

**Step T4 — PGA algebraic identities (~half day)**
File: `test_identities_pga.py`. Same pattern for PGA-specific operations:

- `unitize(unitize(a)) == unitize(a)`  (idempotent)
- `bulk(a) + weight(a) == a`  (bulk + weight decomposition)
- `geom_nrm(unitize(a)).c1 == 1.0`  (unitized objects have weight 1)
- Motor identity: `M ⟇ M~ == 1` for unit motors
- Sandwich preserves grade: `move(point, M)` returns a `point*dp`, not a different grade
- `meet(join(p1, p2), p3) == p3`-on-line check (incidence)
- `dist(p, p) ≈ 0` (same point has zero distance)
- `cmpl(cmpl(a)) == ±a` (with sign depending on algebra dimension)

**Step T5 — Cross-check harness (~1–2 days)**
File: `test_cross_check.py`. Most expensive step, highest confidence.

Approach:

1. Add a small **harness in `ga_test/`** that emits a JSON file of selected (operation, inputs, expected output) tuples. Use existing test cases verbatim. CMake target `make export_python_test_cases`.
2. JSON schema:

   ```json
   {
     "version": 1,
     "cases": [
       { "id": "ega3d_dot_basic_001",
         "op": "dot",
         "submodule": "ega",
         "args": [
           {"type": "vec3d", "fields": [1.0, 2.0, 3.0]},
           {"type": "vec3d", "fields": [4.0, 5.0, 6.0]}
         ],
         "expected": {"type": "scalar3d", "fields": [32.0]},
         "tolerance": 1e-12
       }
     ]
   }
   ```

3. `test_cross_check.py` loads JSON, reconstructs Python objects via the type registry, calls the bound function, compares with `is_congruent`-style relative tolerance.

Picking cases: target one case per (algebra, op, type-pair) combination — that's a few hundred entries, not exhaustive but representative. Coverage of every overload is unrealistic; coverage of every dispatch path is.

**Step T6 — Stub validation (deferred until `.pyi` files exist)**
File: `tests/test_stubs_mypy.py`. After `ga_bindgen` gains stub emission (item 6 in §12 of this doc):

```python
def test_mypy_strict_passes_on_examples():
    result = subprocess.run(
        ["mypy", "--strict", "ga_py/examples/"],
        capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout
```

#### Success criterion

When all of T1–T5 pass, the Python wrapper has:

- 114 smoke tests covering construction / printing / dispatch (existing)
- ~50 constant verification tests (T1)
- ~20 grade lookup tests (T2)
- ~30 EGA algebraic identities × hypothesis (T3) — each with arbitrary input shrinkage
- ~30 PGA algebraic identities × hypothesis (T4)
- ~200 cross-check cases against C++ reference output (T5)

Total: roughly 400 tests + property-based fuzzing across both algebras. That's the threshold for confident "the bindings produce identical results to the C++ library."

#### What we deliberately don't do

- Re-run the entire `ga_test/` suite in Python. Each `ga_test` case is a sample; we want **representative** samples plus **property-based** invariants, which is much more efficient than 1:1 reproduction.
- Test internal `hd::ga::detail::*` helpers (correctly excluded per §10b).
- Test the binding generator itself (`ga_bindgen` is a separate tool with its own concerns).

#### Tools to install in the venv before starting

```bash
pip install pytest hypothesis numpy
# (numpy will become useful when the Tier 1 buffer-protocol work
# from §3 is implemented; the tests can already use it for fixture
# data even if the wrapper doesn't expose buffers yet.)
```

#### Open questions (for the new session)

1. **Tolerance policy.** Single relative tolerance everywhere, or per-operation? The C++ side uses `is_congruent` with ~1e-9 to 1e-12 depending on op. Python tests should mirror that. Likely answer: introduce a `pytest` fixture that returns the right epsilon for the algebra/op being tested.

2. **Whether to also bind / test the `move*_opt` family** beyond `move*dp`. They're already bound but I haven't smoke-checked them; should the cross-check cases include the optimized variants?

3. **Hypothesis strategy for non-degenerate inputs.** For `inv(a)`, we need `nrm_sq(a)` not too close to zero. For motors, we need the rotation axis non-ideal. Need bespoke `@st.composite` strategies that filter via `hypothesis.assume(...)`.

4. **Where the JSON export from `ga_test/` lives.** Inside `ga_test/` (so it's regenerated alongside the C++ tests), or in `ga_py/tests/data/` (committed snapshot)? Recommend: emit during `ga_test` build, commit the snapshot in `ga_py/tests/data/`.

5. **Should the cross-check harness also run on the actual `ga_test/` doctest cases by parsing them?** Probably not — explicit JSON is cleaner and decouples the two test harnesses.

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

## 11. Packaging

### 11.1 Build system

- `scikit-build-core` as the Python build backend.
- `pyproject.toml` declares dependencies: `nanobind` (build-only), `numpy` (runtime, for buffer protocol interop).
- CMake integration reuses the existing CMake infrastructure; `ga_py/CMakeLists.txt` references `ga/` as include directory.
- nanobind is header-only and vendored via pip install; fmt is header-only; ga/ is header-only. **No vcpkg/brew dependencies on the wheel-building side** → simple CI.

### 11.2 CI / wheel matrix (v1)

Python 3.10 × {macOS (Intel + ARM), Linux (x86_64), Windows (x86_64)} = 4 wheels. `cibuildwheel` on GitHub Actions is standard for this. Source distribution (sdist) also published for users who want to build themselves.

Additional Python versions (3.11/3.12/3.13) added post-v1 by extending the cibuildwheel version list — no design changes needed.

---

## 12. Open questions

Items that are **not yet decided** and need further thought or a spike to resolve:

1. **Manifest JSON schema.** Custom schema (simplest) vs. reuse of an existing tool's format (e.g., cppyy, castxml). Recommendation: custom, kept small. Decide during Spike.
2. **Template-instantiation freeze.** The C++ types are templates on `value_t`. The manifest must record that bindings only instantiate `T = double`. Implementation detail — handle in `emit_nanobind.py` by emitting only the `double` instantiation.
3. **One bindings.cpp per algebra vs. all-in-one.** Compile-time cost vs. link complexity. Start per-algebra; revisit if it matters.
4. **Documentation site.** Sphinx + autodoc — now or later? Lean *later*; docstrings are deferred, so autodoc would have nothing to render.
5. **How to verify numerical equivalence in cross-check tests.** Bit-exact vs. tolerance-based. Expect small differences from FP ordering in rare cases; use tolerance like the existing `is_congruent` tests.
6. **What belongs in `module.cpp` (hand-written) vs. generated.** Entrypoint and exception translator are hand-written. Beyond that, everything generated. Revisit if generation runs into awkward corners.
7. **Regression protection for manifest drift.** If someone edits generated files by hand, regeneration overwrites them. Mitigation: generated files carry a clear "DO NOT EDIT — regenerate via ga_bindgen" header.

---

## 13. Effort estimate (rough)

| Phase | Effort |
| --- | --- |
| Spike: libclang parse of `ga/ga_ega.hpp` → dump functions/types → sanity check | 1–2 days |
| `ga_bindgen` full scanner + manifest | 3–5 days |
| `ga_bindgen` emit nanobind bindings (for a single algebra first) | 3–5 days |
| `ga_bindgen` emit stub files | 1–2 days |
| `ga_bindgen` Lua coverage audit report | 1 day |
| `ga_py/` build scaffolding (CMake + pyproject + scikit-build-core) | 2–3 days |
| CI setup: cibuildwheel × 3 platforms | 2–3 days |
| Test suite (cross-check + hypothesis + smoke) | 3–5 days |
| Example scripts + minimal docs | 1–2 days |
| Iteration / bug fixing / platform-specific fixes | buffer: 3–5 days |

**Total first delivery:** ~4–6 focused weeks, depending on libclang-parsing surprises and platform-specific wheel issues.

---

## 14. Plan — concrete next steps

1. **Document this design** (this file) — **done.**
2. **Spike — libclang viability.** — **done, 2026-04-24.** Results in `build/spike_libclang/SPIKE_RESULTS.md`.
   - Clean parse of `ga/ga_ega.hpp` and `ga/ga_pga.hpp` (zero errors).
   - Enumerated 21 class templates, 39 user typedefs, 154 unique non-operator function names (1137 overloads), 21 unique operators (774 overloads).
   - All user types (`vec3d`, `bivec2dp`, `mvec3dp_u`, …) resolved to their `*_t<value_t>` instantiations.
   - Artifact: `build/spike_libclang/manifest_spike.json` — raw dump, reusable as seed for the real manifest schema.
   - **Key pitfall logged:** pip-bundled `libclang` package has no resource dir — must point `clang.cindex` at homebrew/system libclang and pass `-isystem` for libc++ and SDK paths explicitly. Real generator needs a `find_libclang()` probe step.
3. **Spike — audit ga_lua.hpp** — **done, crude version.** Confirmed genuine gaps in Lua wrapper: `rk4_step`, `rk4_get_time`, `Hz2radps`, `radps2Hz`, `rpm2radps`, `radps2rpm`, `*_metric_view` helpers. Crude string-search under-reports presence for ops bound via `sol::meta_function::*`; a per-overload audit is a small extension.
4. **Next — promote spike to `ga_bindgen/`.** Move `spike.py` to `ga_bindgen/src/scan.py`, define manifest schema formally, add `emit_nanobind.py` that emits bindings for one type end-to-end (e.g., `vec3d`). Build a tiny `ga_py/` stub compiling that single binding and importable in Python. This confirms the full pipeline end-to-end before scaling to all 39 types.
5. **Then — scale.** Generate for all types in EGA → then PGA → stubs → tests → packaging.

Immediate next step (4) is ~1–2 days. Total project remains ~4–6 weeks per §13.

---

## 15. Appendix — example of generated nanobind code (sketch)

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
