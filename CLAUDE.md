# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this
repository.

## Commit policy (MANDATORY)

**Never commit on your own initiative.** After completing and verifying a change, present
the result and wait for the user's explicit request to commit. This includes follow-up
work: an instruction like "apply the fixes" or an earlier "commit everything" authorizes
that step only -- it does NOT extend to committing later changes in the same session.
When in doubt, leave the working tree for review and say it is ready.

## Project Paths and Build System

**CRITICAL PATH INFORMATION:**

- **Project Root**: `/Users/hud3bh/prg/cpp/pj/ga/` (absolute path)
- **Build Directory**: `/Users/hud3bh/prg/cpp/pj/ga/build/` (absolute path)
- **Working Directory**: Always work from the build directory when running executables or
  when creating temporary files
- **Source Files**: Always in `/Users/hud3bh/prg/cpp/pj/ga/[module]/src/` or
  `/Users/hud3bh/prg/cpp/pj/ga/ga_prdxpr/src_prdxpr/`
- **File Organization**: Keep all temporary/debug files in the build directory (`build/`)
  rather than the source tree to maintain clean project organization.

**Build Directory Structure:**

```text
/Users/hud3bh/prg/cpp/pj/ga/build/                  # Main build directory (working directory)
├── ga_lua/ga_lua                                   # Lua interface executable
├── ga_test/                                        # doctest-based test suite
│   ├── ga_ega_test                                 #   Euclidean GA tests (2D/3D)
│   ├── ga_pga_test                                 #   Projective GA tests (2dp/3dp)
│   ├── ga_sta_test                                 #   Space-Time Algebra (STA4D) tests
│   ├── ga_appl2dp_test                             #   PGA2D applications (kinematics/frame trees)
│   ├── ga_appl3dp_test                             #   PGA3D applications (generic geometry/kinematics/mechanics)
│   ├── ga_integrator_test                          #   ODE integrators: RK4 vs ABM2 (+ timing)
│   ├── ga_stencil_test                             #   FD stencil generator (detail/ga_stencil.hpp)
│   ├── ga_export_python_cases                      #   emits cross-check JSON for ga_py tests (python_utilities/)
│   └── ga_sta_bench_transform                      #   STA transform_opt micro-benchmark
├── ga_prdxpr/                                      # Code generator + its own tests/tools
│   ├── ga_prdxpr                                   #   main product-expression generator
│   ├── ga_prdxpr_parser_test                       #   sandwich expression parser tests
│   ├── ga_prdxpr_rule_generator_test               #   rule/metric/dual/complement table dump
│   └── ga_prdxpr_viscmp                            #   visual comparison helper
└── ga_py/_ga_py.cpython-<ver>-<plat>.so            # nanobind Python extension (opt-in, see below)
```

**Notes:**

- The `ga_py` extension is built only when CMake is configured with
  `-D_GA_BUILD_PYTHON=ON` (see the "Python wrapper" section). Without that flag the
  `build/ga_py/` directory is absent.
- `ga_lua` requires Lua + sol2 — the target is skipped when the optional dependency is
  not found (see "Dependencies").
- Source modules without a runtime binary: `ga/` (header-only library), `ga_bindgen/`
  (libclang scanner, run via its own venv — see "Python wrapper").
- **This is the PUBLIC repo.** The visualization viewer (`ga_view`), specific application
  test bundles, and the `ga_docu` LaTeX sources live in a
  private superset repo that embeds this one as a submodule; the public repo ships only the
  generic library, generators, bindings, generic tests, and a generated `ga_docu/ga_docu.pdf`.

**File Path Rules for Code:**

- When reading source files from executables: `../[module]/src/filename` (relative to
  build directory)
- NEVER guess paths - use these established patterns
- All executables run from `/Users/hud3bh/prg/cpp/pj/ga/build/` directory

This project uses CMake (minimum version 3.28) with C++23 standard.

**IMPORTANT: Always use out-of-source builds exclusively in the `build/` directory.**

```bash
# Configure (from project root)
mkdir -p build
cd build
cmake ..

# Build (from build directory)  
cmake --build .

# Clean build (if needed)
cd .. && rm -rf build && mkdir build && cd build && cmake ..
```

**Critical Build Requirements:**

- ALL build artifacts must be generated in `build/` directory only
- NEVER run cmake directly in source directories (ga_lua/, ga_test/, etc.)
- NEVER run make/ninja commands outside `build/` directory
- Build directory structure: `build/ga_lua/`, `build/ga_test/`, `build/ga_prdxpr/`, etc.
- Always run builds from the project root's `build/` directory

**Build Workflow:**

1. Always start from project root directory (`/Users/hud3bh/prg/cpp/pj/ga/`)
2. Ensure `build/` directory exists: `mkdir -p build`
3. Navigate to build directory: `cd build`
4. Configure: `cmake ..`
5. Build: `cmake --build .`
6. Run executables: `cd ga_lua && ./ga_lua && cd ..`, `cd ga_test && ./ga_ega_test &&
   cd..`, etc.
7. For refactoring with reference output data always match reference output 100%
   character-identical

## Running Tests

Execute test binaries from the build directory (ensure you're in `build/` directory):

```bash
cd ga_test && ./ga_ega_test && cd ..    # Tests for Euclidean GA
cd ga_test && ./ga_pga_test && cd ..    # Tests for Projective GA
```

**Runtime note (Debug builds):** `ga_appl3dp_test` (generic geometry/mechanics, ~7670
assertions) is the long pole in a `Debug` build (the default `CMAKE_BUILD_TYPE`) — it can
take a couple of minutes, so a 2-minute command timeout may kill it mid-run. Give it a
generous timeout (>=300 s) and run it on its own rather than chained behind the other
suites; the remaining suites finish in a few seconds each. For a fast full-suite pass,
configure a separate `Release` build (`cmake -DCMAKE_BUILD_TYPE=Release ..`), where it
drops to a few seconds — at the cost of a longer one-off compile.

## Running Applications

Execute applications from the build directory:

```bash
cd ga_lua && ./ga_lua && cd ..                    # Interactive Lua shell
cd ga_lua && ./ga_lua script.lua && cd ..         # Run Lua script
cd ga_prdxpr && ./ga_prdxpr --help && cd ..       # Code generator - see README.md for details
```

**Note:** For complete `ga_prdxpr` usage documentation, see the
[README.md](README.md#ga-product-expression-generator-ga_prdxpr) file.

## Project Architecture

This is a header-only geometric algebra library (`ga/`) with accompanying modules:

- **ga/**: Core header-only GA library supporting both Euclidean (EGA) and Projective
  (PGA) geometric algebra
- **ga_test/**: Test suite using doctest framework
- **ga_lua/**: Lua scripting interface for GA operations
- **ga_prdxpr/**: Code generator for GA product expressions (EGA2D, EGA3D, PGA2DP, PGA3DP)

### Notation Conventions

When writing documentation, comments, or mathematical expressions for this codebase, use
the following notation:

| Symbol       | Unicode | LaTeX                         | Meaning                                     |
| ------------ | ------- | ----------------------------- | ------------------------------------------- |
| ⟑            | U+27D1  | `\wedgedot` (MnSymbol)        | Geometric product (EGA/STA)                 |
| ⟇            | U+27C7  | `\veedot` (MnSymbol)          | Regressive geometric product (PGA)          |
| $\tilde{M}$  | —       | `\tilde{}`                    | Reverse of M (`rev(M)` in code)             |
| $\utilde{M}$ | —       | `\utilde{}` / `\undertilde{}` | Regressive reverse of M (`rrev(M)` in code) |

**Notes:**

- **Markdown files** (`.md`): use `\utilde{}` — works natively in LaTeX (Markdown
  previewers) without a package import
- **LaTeX documentation** (`ga_docu`): use `\undertilde{}` from the `accents` package
  (`\usepackage{accents}`)
- **LaTeX math (`ga_docu` `.tex`)**: write the geometric product as `\wedgedot` (⟑) and the
  regressive product as `\veedot` (⟇), both from `MnSymbol` — not as bare juxtaposition.
  `\wedgedot` is the EGA/STA product, `\veedot` the PGA product; e.g. rotor sandwich
  `R \wedgedot X \wedgedot \rev{R}`, motor sandwich `M \veedot X \veedot \rrev{M}`. Prose
  and `.md` may use the Unicode ⟑/⟇ directly.
- For unit motors, `rev(M)` and `rrev(M)` coincide as group inverses, but the symbols and
  function names are kept distinct throughout this codebase
- Motor composition (PGA) uses ⟇: `M₁ ⟇ M₂` corresponds to `rgpr(M1, M2)` in code
- The rigid-motion sandwich differs by algebra — the product and the reverse must match:
  - **EGA / STA** (rotor): $R ⟑ X ⟑ \tilde{R}$ — geometric product `⟑` with reverse `rev`;
    corresponds to `rotate(X, R)` (EGA) / `transform(X, R)` (STA) in code
  - **PGA** (motor): $M ⟇ X ⟇ \utilde{M}$ — regressive product `⟇` with regressive reverse
    `rrev`; corresponds to `move2dp(X, M)` / `move3dp(X, M)` in code

### Markdown Formatting

When writing or editing `.md` files in this repository, wrap prose at **column 90** to
match the editor's auto-wrap setting (VS Code line break at column 90). This applies to
paragraphs, bullet items, blockquotes, and lead-in sentences before code blocks.

**Exceptions — keep on a single line even if wider than 90 columns:**

- **Tables** — Markdown table rows must stay on one physical line; let columns be as wide
  as needed. Do not soft-wrap table cells.
- **Code fences** — preserve the original code's line lengths verbatim; do not reflow code
  to fit 90 columns.
- **URLs and reference-style link targets** — let them overrun rather than break the link.
- **Headings** — keep on one line even if a long heading exceeds 90.

### C++ Formatting (clang-format — MANDATORY before writing/committing)

**Run clang-format on every C++ file you generate or edit, before considering the change
done.** The repo style is a **committed repo-local `.clang-format` at the repo root**
(added so co-workers get identical formatting regardless of their global config; it is a
verbatim copy of the maintainer's former global `~/.clang-format`). clang-format finds it
by searching upward from the file's directory, so the file must live **inside the repo
tree** when formatted (a temp file under `/tmp/` silently falls back to LLVM defaults — see
the "clang-format gotcha" note later in this file).

**Version baseline:** the whole tree was reformatted once with **clang-format 22.x** (the
config alone does not pin output — different clang-format versions reorder includes and
break lines differently). Use a **matching major version (22.x)** so a format-on-save does
not re-churn the repo. If you must use a newer major, expect a one-time repo-wide reflow and
raise it before committing.

```bash
clang-format -i path/to/file.cpp   # format in place (picks up the repo-root .clang-format)
```

Why this matters: the editor reformats on save with the same config, so unformatted output
gets reflowed later — which **churns line numbers and re-sorts `#include` groups**,
producing noisy, hard-to-review diffs and stale `file:line` references. Formatting up
front avoids that.

Key settings to match even when hand-writing (full config in the repo-root `.clang-format`):

- **`ColumnLimit: 90`** — the hard wrap limit (this is the minimum to respect).
- `IndentWidth: 4`, `UseTab: Never`; `BreakBeforeBraces: Stroustrup`; `IndentCaseLabels:
  true`.
- `PointerAlignment: Left` / `ReferenceAlignment: Left` → `Type* p`, `Type& r` (compatible
  with the east-const convention below: `Type const& name`).
- Includes are **sorted within each group** — keep includes grouped and let clang-format
  settle the order rather than fighting it.
- `EmptyLineBeforeAccessModifier: Always`, `EmptyLineAfterAccessModifier: Always`;
  `MaxEmptyLinesToKeep: 2`.

#### Comment formatting — protect display formulas from reflow

clang-format reflows each run of consecutive `//` lines as one paragraph (ReflowComments),
which **collapses an indented formula or matrix into the surrounding prose** when nothing
separates them. To keep a display block readable:

- **Bound every display block with blank `//` lines.** An empty comment line is a
  paragraph boundary clang-format will not cross, so the block survives intact. Without
  it, a `//   X = ...` line wedged between two prose lines gets reflowed into the text.
- **Indent display lines** (4–5 spaces after `//`) so equations read as set-apart blocks.
- **Keep multi-line blocks contiguous** — matrix rows / multi-line equations take blank
  `//` lines only at the block's *outer* boundaries, never between the rows.
- **A definition list is the opposite of a matrix** — its items are *separate* blocks, so
  put a blank `//` between every item and keep each item to one line; otherwise reflow
  merges adjacent items and orphans the wrapped words. Aligning the `:` reads well (see
  the `joint{2,3}dp` enums).

```cpp
// prose introducing the step ... :
//
//     P_new = P (x) exp(0.5 * B_rel * dt)        [P = body->parent motor]
//
// prose continuing after the formula ...
```

The mechanics headers (`ga/ga_pga{2,3}dp_ops_mechanics.hpp`) follow this throughout.

### Library Usage Patterns

The library provides two main entry points:

- `#include "ga/ga_ega.hpp"` for Euclidean GA (2D, 3D, 4D)
- `#include "ga/ga_pga.hpp"` for Projective GA (2dp, 3dp)

Critical usage requirements:

1. All GA code must be wrapped in try-catch blocks to handle exceptions
2. Use namespace `hd::ga` and either `hd::ga::ega` or `hd::ga::pga`
3. **Always enclose left/right contractions (`<<` and `>>` operators) in parentheses due
   to operator precedence**
4. The base scalar type (float/double) is configured globally in `ga/ga_value_t.hpp`

### Header layering — where to find a `ga/` operation

Each algebra's operations are split across a fixed set of `ga/ga_<alg>_ops*.hpp` headers
(`<alg>` ∈ `ega2d, ega3d, pga2dp, pga3dp, sta4ds`), layered by what they build on. **Every
file lists the functions it provides in a `// provides ... operations:` comment right after
the `namespace hd::ga{,::pga}` declaration** — read that block first when hunting for a
function; it is the authoritative per-file index. The split:

| File | Provides (operations defined there) |
| ---- | ----------------------------------- |
| `ga_<alg>_ops_basics.hpp` | involutions (`gr_inv`, `rev`, `rrev`, `conj`); complements (`l_cmpl`/`r_cmpl`, `cmpl`); duals (`*_bulk_dual`, `*_weight_dual`); `bulk`/`weight`; norms (`bulk_nrm{,_sq}`, `weight_nrm{,_sq}`, `geom_nrm{,_sq}`); `bulk_normalize`, `unitize` |
| `ga_<alg>_ops_products.hpp` | `dot`/`rdot`; `wdg`/`join`, `rwdg`/`meet`; contractions (`<<`, `>>`, `*_bulk/weight_contract`); expansions (`*_bulk/weight_expand`); `cmt`/`rcmt`; `operator*`(=`gpr`)/`rgpr`; `inv`/`rinv` |
| `ga_<alg>_ops.hpp` | higher-level ops built on basics+products: `angle`; **`exp`/`log`/`sqrt`** (w.r.t. `gpr` for EGA/STA, `rgpr` for PGA); `get_motor*`; `move{2,3}dp`/`rotate`; projections/rejections, `reflect_on`/`invert_on`, `expand`, `att`, `dist*`, `is_congruent` |
| `ga_<alg>_ops_mechanics.hpp` (PGA2DP/3DP only) | rigid-body dynamics: `Inertia{2,3}dp`, `pose`/`motor` converters (`motor_from_pose3dp`, `pose3dp_from_motor`), moving-frame kinematics, `static_/kinematic_/dynamic_system{2,3}dp`, force elements (`grounded_spring`), and the `extra_wrenches()` subclass hook for application-specific wrenches |
| `ga_<alg>_ops_constraints.hpp` (PGA2DP/3DP only) | the opt-in `closed_loop_system{2,3}dp` KKT layer |

So e.g. `exp`/`log` live in `ga_<alg>_ops.hpp` (not basics/products); the pose↔motor and
`rcmt`-velocity-field helpers used together with them are in `ops_mechanics.hpp` (PGA) and
`ops_products.hpp` (`rcmt`). EGA/STA have no `mechanics`/`constraints` headers.

#### Constants, user types, utilities, the solver (the non-`ops` headers)

The per-algebra `ops*.hpp` above are the operations. The cross-cutting building blocks
people reach for repeatedly live in a fixed set of headers — check here before grepping
the whole tree or hand-rolling a helper:

| File | Provides |
| ---- | -------- |
| `ga/ga_value_t.hpp` | the scalar type `value_t` (float/double switch); `eps` (equality) and `eps_congruent`; `is_congruent` + `detail::coeffs_congruent` |
| `ga/ga_usr_consts.hpp` | named constants per algebra: basis blades (`e1_3dp`, `e23_3dp`, …), **origins** (`O_2dp`, `O_3dp`), projection/attitude blades (`e423_3dp`, …), and `pi` |
| `ga/ga_usr_types.hpp` | the user value-type aliases (`vec3dp`, `bivec3dp`, `mvec3dp{,_e,_u}`, `scalar2d`, … — the `value_t` instantiations of the templates) |
| `ga/ga_usr_types_mechanics.hpp` | physics aliases (`Inertia{2,3}dp`, `pose{2,3}dp`, kinematic frame/system types); **included after** the physics `ops` headers (it aliases templates they define) |
| `ga/ga_usr_utilities.hpp` | `deg2rad`/`rad2deg`/`rpm2radps`; **`rk4_step`** (an `mdspan` form and a `std::vector` form); step/easing helpers `linear_step`/`smooth_step`/`smoother_step` |
| `ga/detail/ga_solver.hpp` | the small dense linear solver — `lu_solve`, `lstsq_solve` (least-squares / minimum-norm via normal equations), `kkt_solve` (constrained KKT); used by the physics assembly and the closed-loop layer |
| `ga/detail/ga_stencil.hpp` | finite-difference stencil generator — `stencil_t` (Fornberg-style weights, order, truncation error via `lu_decomp`/`lu_backsubs`; explicit + compact schemes), `factorial`; infrastructure for discretized field derivatives (numerical nabla, planned STA electrodynamics) |

`ga/detail/` also holds the foundation the user types are built on (`ga_core_types`,
`ga_ega_types`, `ga_pga_types`, `ga_sta_types`, `ga_error_handling`, `ga_fmt_support`).
The two umbrella headers `ga/ga_ega.hpp` and `ga/ga_pga.hpp` include the whole stack in
the right order, so user code includes only those (see "Library Usage Patterns"). When a
numeric helper seems missing, the solver and `rk4_step` above are the usual answer —
don't re-derive them.

### Type System Architecture

- Template types: `Scalar2d<T>`, `Vec2d<T>`, `BiVec2d<T>`, `MVec2d_E<T>`, etc.
- User convenience types: `scalar2d`, `vec2d`, `bivec2d`, `mvec2d`, etc. (based on
  `value_t`)
- Multivector types differentiate between even (`_E`), odd (`_U`), and full (`MVec`)
  multivectors
- Function call arguments for scalar and pscalar use `Scalar2d<T> a`, etc. directly
  instead of `Scalar2d<T> const& a` (cheap to copy), while for vector and bivector, etc.
  use `Vec2d<T> const&` for unmodifiable arguments
- **CRITICAL**: ALWAYS use "east const" convention (see dedicated section below)

### East const convention (MANDATORY)

**This codebase uses "east const" (const on the right) convention throughout.**

When writing or modifying ANY C++ code in this repository, you MUST follow these rules:

**✅ CORRECT (East const):**

```cpp
// Function parameters - references
void foo(std::string const& str);
void bar(AlgebraConfig const& config);
void baz(prd_rules const& rules);
Vec3d<T> const& getVector();

// Function parameters - pointers
void process(int const* data);
char const* getMessage();

// Local variables and member variables
std::vector<int> const values = {1, 2, 3};
auto const& item = container[0];
size_t const n = basis.size();
mvec_coeff const& basis = config.multivector_basis;

// Return types
std::string const& getName() const;
std::vector<int> const getValues();
```

**❌ INCORRECT (West const - DO NOT USE):**

```cpp
// DO NOT write any of these patterns:
const std::string& str           // Wrong! Use: std::string const& str
const AlgebraConfig& config      // Wrong! Use: AlgebraConfig const& config
const auto& item                 // Wrong! Use: auto const& item
const size_t n                   // Wrong! Use: size_t const n
```

**EXCEPTIONS (Keep "const" on the left in these cases ONLY):**

1. **Static and extern declarations with initialization:**

   ```cpp
   static const std::string one_str{"1"};     // OK - static initialization
   extern const prd_rules gpr_rules;          // OK - extern declaration
   inline const std::string& one_str() { ... } // OK - inline function return
   ```

2. **Member function const qualifiers:**

   ```cpp
   void getValue() const;              // OK - member function qualifier
   int operator[](size_t i) const;     // OK - const member function
   ```

3. **Const iterator types (standard library convention):**

   ```cpp
   std::vector<int>::const_iterator it;  // OK - STL convention
   ```

**WHY East const?**

- **Consistency**: Type qualifiers read left-to-right: `int const*` = "pointer to const
  int"
- **Clarity**: `T const&` clearly shows "const reference to T"
- **No Ambiguity**: With templates, `T const` works uniformly regardless of what T is
- **This Project's Standard**: All existing code in ga/ and ga_prdxpr/ follows this
  convention

**IMPORTANT FOR CODE GENERATION:**

When you generate new C++ code (functions, classes, templates), you MUST:

- Use east const for ALL function parameters taking references: `Type const&`
- Use east const for ALL function parameters taking pointers: `Type const*`
- Use east const for ALL local const variables: `auto const`, `size_t const`, etc.
- Use east const for ALL return types: `std::string const&`, `Type const`, etc.

**This is NOT optional** - all code must follow this convention for consistency with the
existing codebase.

### Function-body style for ops_products

In `ga/*_ops_products.hpp`, primitive product functions follow a uniform
canonical form (enforced project-wide; produced by
`ga_prdxpr --output=code`):

- **Flat-constructor returns** use temp-vars:

  ```cpp
  ctype const c0 = ...;
  ctype const c1 = ...;
  return Type<ctype>(c0, c1, ...);
  ```

  Applies to `MVec`, `MVec_E`, `MVec_U`, `Vec`, `BiVec`, `TriVec`, `Plane`,
  `Line`, etc. Single-value (`Scalar`, `PScalar`) returns stay inline.
  Composite returns (`MVec_E(Scalar(...), BiVec(...))`) keep nested
  sub-constructors.

- **Zero-result functions** use the `[[maybe_unused]]` form with
  anonymous parameters and `return Type<ctype>(0.0);` (paren-init, not
  `{0.0}` brace-init).

- **Strong-typed `Scalar`/`PScalar` operands** must be wrapped with
  `ctype(...)` when used in coefficient arithmetic, e.g.
  `ctype const c0 = ctype(s) * M.c0;`. Without the wrap the strong type
  doesn't implicitly convert to `ctype`.

- **Semantic delegations are intentional**: ~32 functions across the
  four algebras delegate to other primitives (e.g.
  `gpr(v,v) = MVec_E(dot(v,v), wdg(v,v))`). These encode GA identities
  and are kept by design — codegen emits the expanded form there, so a
  re-validation will report ~32 mismatches as expected, not as bugs.

### Dependencies

Required: fmt library (header-only)
Optional: doctest (testing), Lua + sol2 (scripting), readline
(for enhanced interactive scripting experience)

Install on macOS: `brew install fmt doctest lua readline`

### Windows Dependency Management (vcpkg + System Lua)

**Critical Windows Build Requirements:**

The Windows build uses a hybrid approach: vcpkg for most dependencies + direct system Lua
static linking. This avoids ABI compatibility issues between MSVC and different Lua
implementations.

**Build Command Sequence** (PowerShell — the default Windows shell; in Git Bash use the
POSIX equivalents `rm -rf build && mkdir build && cd build`):

```powershell
cd <your checkout root>
Remove-Item -Recurse -Force build -ErrorAction SilentlyContinue
mkdir build; cd build
cmake .. -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build . --target ga_lua --config Debug
```

**Expected Success Indicators:**

- ✅ `"Using system-installed Lua with static linking (no vcpkg)"`
- ✅ `"✓ Found static system lua: headers at C:/Users/[USERNAME]/AppData/Local/Programs/Lua/include, lib at C:/Users/[USERNAME]/AppData/Local/Programs/Lua/lib/lua54.lib"`
- ✅ `"ga_lua: Using sol2 headers from ...vcpkg_installed/x64-windows/include"`
- ✅ `"ga_lua: Include order: Lua (...) then sol2 (...)"`
- ✅ Final: `"ga_lua.vcxproj -> ...ga_lua.exe"` (NO DLL copying message)

**Troubleshooting Windows Build Issues:**

1. **CMake Variable Issues**: Check `ga_lua: DEBUG - LUA_INCLUDE_DIR=` and
   `LUA_LIBRARIES=` are populated
2. **Include Order Problems**: Lua headers MUST come before sol2 in
   `target_include_directories`
3. **Static Linking Verification**: Should link to `lua54.lib`, NOT `lua54.dll`
4. **vcpkg Integration**: Ensure `CMAKE_TOOLCHAIN_FILE` points to vcpkg toolchain

**Key Architecture Files:**

- `cmake/vcpkg_dependencies.cmake` - System Lua detection + vcpkg integration
- `ga_lua/CMakeLists.txt` - Include ordering + static linking (no DLL copying)
- `vcpkg.json` - Dependencies without Lua (fmt, doctest, sol2 only)

**Prerequisites:**

- vcpkg installed at `C:/vcpkg/`
- System Lua 5.4+ installed at `C:/Users/[USERNAME]/AppData/Local/Programs/Lua/`
- MSVC compiler (Visual Studio 2022 recommended)

### Development Notes

- Compiler definitions: `-D_HD_GA_EXTENDED_TEST_DIV_BY_ZERO` (extended testing),
  `-D_HD_GA_HAVE_READLINE`(readline support)
- The library supports switching between debug/release builds via CMAKE_BUILD_TYPE
- MSVC uses `/bigobj` flag due to template instantiation complexity
- Additional compiler definitions start with `_HD_GA_` in order to be consistent
  with already existing definitions
- Use the existing cmake infrastructure as much as possible to test changes

## GA Product Expression Generator (ga_prdxpr/)

The `ga_prdxpr/` directory contains a sophisticated **code generator** that produces
optimized C++ expressions for geometric algebra operations. It supports four complete
algebras: EGA2D, EGA3D, PGA2DP, and PGA3DP, plus STA4D (G(1,3,0), algebra name
`sta4ds`). STA4D's `.cases` are populated for `dot, gpr, wdg, l/r contraction,
rwdg, cmt` (354 signatures, all pasted into the library); `rdot, rcmt, rgpr`
remain empty by design (the regressive products are fleshed out for PGA).

### C++ Code Generation (`--output=code`)

`ga_prdxpr --output=code` emits ready-to-paste C++ implementations for primitive
products in the canonical form documented under "Function-body style for
ops_products" above. Source files `ga/*_ops_products.hpp` are ~95% char-identical
with this output post-`clang-format`; the remaining ~5% are intentional semantic
delegations.

- Module: `ga_prdxpr/src_prdxpr/codegen/` (`ga_codegen_types.{hpp,cpp}`,
  `ga_codegen_emitter.{hpp,cpp}`).
- Validation pipeline: `ga_prdxpr/src_prdxpr/codegen/tools/`
  (`diff_codegen.py`, `inline_to_tempvars.py`, `expand_delegations.py`).
- Default output (`ga_prdxpr` with no args) is a stable reference (currently 25081
  lines, see the README's "Verification" section) — `code` is opt-in only.
- See `ga_prdxpr/README.md` for usage, the invocation matrix, and the
  "Open Codegen Work" section (sandwich product codegen, STA4D rollout).

### Codegen gotchas: non-primitive products and hand-coded duals

- **Expansions are non-primitive**: `l_expand`/`r_expand` have no entry in the
  emitter's `product_to_cpp_function` map, so `--output=code` skips them (logs
  `// SKIP ... non-primitive product 'l_expand'`) even though their `.cases`
  exist for `coeffs`/`tables`. Implement them in the library as thin wrappers
  over dual+wedge: `l_expand(a,b) = wdg(l_dual(a), b)`,
  `r_expand(a,b) = wdg(a, r_dual(b))`.
- **Duals are complement∘metric** (`dual(a) = cmpl(G·a)`; degenerate PGA splits into
  bulk/weight variants). Hand-transcribing the complement/dual tables from
  `ga_prdxpr_rule_generator_test` output is error-prone — a flipped-sign bivector dual in
  sta4ds made the contraction identity `a << b == rwdg(l_dual(a), b)` fail *only* at grade
  2 (the rule generator was correct; the hand copy wasn't). Guard every algebra with a
  **transcription gate** test: `dual(e) == nrm_sq(e) * cmpl(e)` per unit basis blade
  (bulk/weight and l/r variants for PGA). The contraction is the metric interior product
  (`= dot` on equal grades); the `rwdg(dual)` form is exact at every grade only when the
  dual is correct. **`--output=code` now emits the complements and duals**
  (`l_cmpl`/`r_cmpl`/`cmpl`, `l_dual`/`r_dual`/`dual`, PGA `bulk`/`weight`; graded inputs
  as flat forms, aggregates `mv`/`mv_e`/`mv_u` as grade-wise delegations) — so they can be
  regenerated and spliced into `*_ops_basics.hpp` (see
  `utilities/splice_generated_code.py` below) instead of hand-transcribed. The
  transcription gate stays valuable as a guard.

### Key Components

**Configuration System:**

- `ga_prdxpr_*.hpp`: User defined information on algebra, rules and coefficients
- `ga_prdxpr_*_config.cpp`: Algebra-specific configurations defining basis, coefficients,
  and product rules to be used for product generation
- **AlgebraConfig**: Complete algebra specification (basis elements, complement rules,
  etc.)
- **ProductDefinition**: Configures product types (geometric, wedge, contraction,
  sandwich, etc.)
- **ProductCase**: Individual expression cases with specific coefficient combinations

**Generation Engine:**

- `generator/ga_prdxpr_generator.cpp`: Main generation logic with algebra-specific
  handlers
- **Dimensional dispatch**: 2D, 3D, 4D specialized implementations
- **Product-specific handlers**: Each ProductType gets specialized mathematical treatment

### Critical Mathematical Patterns

**Complement Transformations for regressive products:**

- for odd dimensinal algebras ega3d, pga2dp: `cmpl(operation(cmpl(A), cmpl(B)))`
- for even dimensinal algebras ega2d, pga3dp, sta4ds: `l_cmpl(operation(r_cmpl(A), r_cmpl(B)))`

**Coefficient Usage Patterns:**

```cpp
A/B                   // Pure multivector × multivector operations
A_even/B_even         // Even part operations (rotors)
R_even/R_rev_even     // Rotor operations (EGA)
R_odd/R_rrev_odd      // Motor operations (PGA)
svps/svps             // Symmetric scalar/vector/pseudoscalar
svps1/svps2           // Asymmetric patterns (e.g., v1.x*v2.y)
```

**Sandwich Products (Two-Step Pattern):**

1. `rotor/motor * object → intermediate_result` (with parentheses via
   `brace_switch::use_braces`)
2. `intermediate_result * rev/rrev(rotor/motor) → final_result`

**Section ordering inside `.cases` blocks:**

Single-grade entries are grouped by section comments in this order, top-down:

- 4D algebras (pga3dp, sta4d): `// mv` → `// mv_e` → `// mv_u` → `// ps`
  → `// trivec` → `// bivec` → `// vec` → `// s`
- 3D algebras (ega3d, pga2dp): same order, no `// trivec` section
- 2D algebra (ega2d): no `// mv_u` and no `// trivec` (those filters
  don't exist in `filter_2d`)

Used uniformly across `gpr`, `wdg`, `rgpr`, `rwdg`, `cmt`, `dot`, etc.

### Algebra-Specific Behaviors

**EGA2D/EGA3D (Euclidean):**

- Standard complement sequences
- Rotor-based sandwich products with `rev()` operations based on geometric product `gpr()`
- Coefficient pattern: R_even/R_rev_even

**PGA2DP/PGA3DP (Projective):**

- Regressive complement transformations
- Motor-based sandwich products with `rrev()` operations based on regressive geometric
  product `rgpr()`

**Success Criteria:**

- Character-perfect match with reference implementation
- All product types working (geometric, wedge, contraction, sandwich, etc.)
- Mathematical accuracy verified against known GA identities

### Common Pitfalls

1. **Coefficient Mismatches**: Each ProductType needs algebra-specific coefficient
   analysis
2. **Operator Precedence**: Always parenthesize left/right contractions (`<<`, `>>`) in GA
   expressions
3. **Output Format**: Case descriptions must match mathematical reality exactly

### Working with ga_prdxpr configs

A small set of conventions and tools makes editing `ga_prdxpr_<algebra>_config.cpp`
predictable.

#### `case_name` format

Both forms are accepted by the parser, the codegen, and the validator:

- Infix:      `"X * Y -> Z"`, `"X ^ Y -> Z"`, `"X << Y -> Z"`, `"X >> Y -> Z"`
- Functional: `"gpr(X,Y) -> Z"`, `"wdg(X,Y) -> Z"`, `"dot(X,Y) -> Z"`, …

The functional form is the standardised style across the configs. The codegen reads
only the result token after `->`; the actual C++ function name comes from
`ProductConfig.product_name`, never from `case_name`. So `case_name` is purely a
descriptive label — it can be reformatted freely.

#### Typed-zero convention (`-> 0` vs `-> 0 ps`)

When a case is identically zero, only two result-side forms are allowed:

- **Bare `-> 0`** — scalar zero (codegen emits `Scalar*<ctype>(0.0)`). This is the
  default for any zero result.
- **`-> 0 ps`** — pseudoscalar zero (codegen emits `PScalar*<ctype>(0.0)`). Use this
  only when the natural return type would be pseudoscalar.

Other typed forms (`-> 0 s`, `-> 0 vec`, `-> 0 bivec`, …) are NOT used:

- `-> 0 s` is redundant with bare `-> 0`.
- `-> 0 vec`/`-> 0 bivec`/etc. break the convention.

Validator check G enforces this (see below).

#### Validator: checks A–G

Every `ga_prdxpr` run validates each declared `OutputCase` and emits warnings to
stderr (with an end-of-run summary). Implementation lives in
[ga_prdxpr_generator.cpp](ga_prdxpr/src_prdxpr/generator/ga_prdxpr_generator.cpp)
inside `validate_case`.

| Check | Catches                                                                                                                                        |
| ----- | ---------------------------------------------------------------------------------------------------------------------------------------------- |
| **A** | `case_name` LHS/RHS token does not match `left_filter_name`/`right_filter_name`                                                                |
| **B** | Result token after `->` is not a known filter for this algebra                                                                                 |
| **C** | Computed non-zero components fall outside the declared result type's basis support (suggests the minimal sufficient type)                      |
| **D** | Functional `case_name` function name does not match `ProductConfig.product_name` (e.g. a `gpr(...)` entry sitting inside a `wdg` config)       |
| **E** | Declared a non-zero result type but the actual computation is identically zero (suggests `-> 0` or `-> 0 ps`)                                  |
| **F** | Coefficient with `_even`/`_odd` suffix paired with a filter other than `mv_e`/`mv_u` (parity mismatch — typically silently produces all zeros) |
| **G** | Typed-zero uses anything other than `ps` (e.g. `-> 0 s` or `-> 0 vec`)                                                                         |

Warnings are non-fatal — generation continues so a single run surfaces every issue.

#### Workflow for adding new product entries

The proven recipe (used to expand wdg, rwdg, dot, rdot, cmt across algebras):

1. **Find the matching `gpr` block** for the same algebra. `gpr` is the canonical
   complete enumeration of `(LHS, RHS)` pair combinations — it covers every grade
   pair that has a meaning in the algebra.
2. **Copy** the relevant section (typically the single-grade lines under
   `// ps` / `// trivec` / `// bivec` / `// vec` / `// s`) into the target
   product's block.
3. **Rename** `gpr(` to the target product (`wdg(`, `cmt(`, `rwdg(`, …).
   Initial result types come from `gpr` and will mostly be wrong — the validator
   will fix them.
4. **Run `ga_prdxpr` for the algebra**. Validator check E lists every entry whose
   computation is identically zero, and check C lists every entry whose declared
   type is too narrow.
5. **Apply suggested fixes**. Zero results → `-> 0` (or `-> 0 ps` if the suggested
   type is pseudoscalar). Tighten over-wide types (e.g. `mv_e` → `bivec`) when
   only one grade survives — `cmt` typically picks the antisymmetric grade out of
   `gpr`'s wider result.
6. **Re-run** until clean. Most products converge in 1–2 iterations.

Single-grade × single-grade pair completeness is the standard scope; only
`gpr`/`wdg` carry the `// mv_e` and `// mv_u` cross-grade sections.

#### Discovering non-zero entries (the "tight scalar" scan)

Use this technique when you don't know which `(LHS, RHS)` pairs have non-zero
results — typical for `l_*_expand` / `r_*_expand` (PGA expansions) and similar
products where most grade combinations yield zero.

Trick: declare every candidate pair with the **narrowest non-zero type**
(`-> s`). Then for each entry the validator emits exactly one of two
warnings, which together give you the full picture:

- **Check E** fires when the computation is identically zero — you can drop
  the entry (or keep as `-> 0` / `-> 0 ps` if you want completeness).
- **Check C** fires when the result has components outside `s`'s mask, and
  its `suggested minimal result type: 'X'` line tells you the right type
  to declare.

Why `-> s` and not `-> mv`: with `-> mv` everything fits, so check C never
fires — you'd see only the zeros, not the suggested types. The narrowest
non-zero type forces check C to suggest the minimal sufficient type for
every non-zero entry.

Recipe:

1. **Fill the empty (or partial) block** with all single-grade pairs declared
   as `-> s`, using the standard coefficient pattern (`svBps1`/`svBps2` for
   same-grade pairs, `svBps`/`svBps` for cross-grade — substitute
   `svBtps`/`svBtps1`/`svBtps2` for 4D algebras).
2. **Run `ga_prdxpr --algebra=<a> --product=<p>`** and collect:
   - `grep "identically zero"` → the zero entries (drop or mark `-> 0`).
   - `grep "suggested minimal"` → the non-zero entries with their proper
     result types.
3. **Filter** the non-zero set by the result-type criterion you care about
   (e.g. for expansions: keep only `vec`/`bivec`/`trivec` results — drop
   entries that produce `s` or `ps`, which carry less geometric meaning).
4. **Replace the scan block** with the filtered list, with proper section
   labels (`// ps`, `// trivec`, `// bivec`, `// vec`, `// s`).

This was the workflow used to populate all PGA bulk/weight expansions and
ega3d / sta4d `l_expand` / `r_expand`. Typically converges in one
scan + one finalisation edit per block.

#### Algebra-specific completeness scope

Not every product is meaningful in every algebra. Convention:

- **`gpr`, `wdg`, `dot`, `cmt`** — defined for every algebra; complete coverage
  expected (`mv,mv` plus all relevant pair combinations).
- **`rgpr`, `rwdg`, `rdot`, `rcmt`** — regressive products are most useful in PGA
  (where the dual structure carries the projective geometry). They appear in
  other algebras too (e.g. `rwdg` is used in EGA applications), but the more
  complete coverage of `rgpr`, `rdot`, and `rcmt` is added to PGA only.
- **`cmt` scope**: just `cmt(mv, mv)` plus all single-grade × single-grade pairs
  (no `mv_e`/`mv_u` cross-grade entries). Deliberate scope decision.

#### Validation utilities

Audit scripts complementing the in-process validator live in
[ga_prdxpr/src_prdxpr/algebras/validation_utilities/](ga_prdxpr/src_prdxpr/algebras/validation_utilities/):

- `pair_coverage.sh` — pair-count matrix per (algebra, product).
- `compare_products.sh <p1> <p2> <alg>` — diff two products' pair sets within one
  algebra.
- `compare_algebras.sh <product>` — diff one product's pair sets across same-dim
  algebras.
- `library_coverage.py` — compare what `ga_prdxpr --output=code` emits against
  what's in the library headers: `ga/ga_<algebra>_ops_products.hpp` **plus**
  `ga_<algebra>_ops_basics.hpp` (where the generated complements/duals are
  spliced). Reports which generated functions
  still need to be copy-pasted into the library; supports `--show-code` for
  paste-ready blocks and `--diff` for body comparison after clang-format
  normalisation. Uses `--algebra=ALGEBRAS` (comma-separated, mirrors
  `ga_prdxpr`'s CLI).

The validator finds **wrong** declarations; these scripts find **missing** ones.

A companion tool that *applies* generator output (rather than only reporting mismatches)
lives in [ga_prdxpr/src_prdxpr/utilities/](ga_prdxpr/src_prdxpr/utilities/):

- `splice_generated_code.py` — regenerate selected products via
  `ga_prdxpr --output=code`, clang-format them, and splice the matching function blocks
  into a library header **by signature** (comments, namespaces, hand-written delegations
  and unrelated products are left intact). Used to keep `dot` / contractions
  (`*_ops_products.hpp`) and the complements / duals (`*_ops_basics.hpp`) in sync.

**clang-format gotcha (bit us once):** clang-format finds the project style (the repo-root
`.clang-format`) by searching **upward from the input file's directory**. A temp file
under `/tmp/` is not under the repo, so clang-format silently falls back to LLVM defaults
(2-space indent, `Type const &name`) and the output will not match the library. Always
clang-format files **inside the repo tree** — the splicer writes its temp file into the
repo for exactly this reason.

### Key Files for Modifications in ga_prdxpr subfolder

Source tree is split by role under `ga_prdxpr/src_prdxpr/`:

- `generator/ga_prdxpr_generator.cpp`: core generation logic and product-specific handlers
- `algebras/ga_prdxpr_<alg>_config.cpp`: per-algebra `ProductConfig` builders — the
  primary user-editable surface for adding products or tweaking coefficient patterns
- `algebras/ga_prdxpr_<alg>.{hpp,cpp}`: per-algebra basis, coefficient strings, and
  rule-set instantiation at static-init time
- `rules/ga_prdxpr_rule_generator.cpp`: rule-generation engine (basis multiplication,
  metric / dual / complement rule generation) — shared across algebras
- `core/ga_prdxpr_common.{hpp,cpp}`: foundational types and table operations
- `codegen/ga_codegen_*.{hpp,cpp}`: `--output=code` C++ emitter

**CMakeLists ordering constraint**: in `SOURCES`, every
`algebras/ga_prdxpr_<alg>_config.cpp` must be listed *before* the matching
`algebras/ga_prdxpr_<alg>.cpp`. The `*.cpp` TU's static initializer calls
`get_<alg>_algebra_config()` (defined in the `_config.cpp` TU), which reads
namespace-scope `const` globals that have internal linkage — one copy per TU. Listing the
config TU first preserves the dynamic-init order the layout relies on.

**Important**: The ga_prdxpr system is a **complete, production-ready geometric algebra
code generator** that produces mathematically accurate, optimized C++ expressions for all
supported algebras.

### Static Initialization Order Safety in ga_prdxpr

Generating complement/dual rules at static-init time caused a static-init-order
fiasco (segfaults). Root causes: inter-dependent static string constants, and
`const` rule tables initialized from `extern` rules not yet constructed.

Safety patterns (keep these when adding rules):

1. **Constants as inline functions**, not static variables:

   ```cpp
   // UNSAFE: static const std::string one_str{"1"s};
   inline const std::string& one_str() { static const std::string s{"1"s}; return s; }
   ```

2. **Staged init**: generate all rules in one pass (`generate_algebra_rules(config)`),
   then assign to the individual `const` tables — no cross-references during init.
3. **No extern references in `const` initializers**: give derived (dual) rule tables
   explicit literal values, not `= some_extern_rules`.

## Python wrapper (ga_py / ga_bindgen)

`ga_bindgen/` scans `ga/*.hpp` with libclang into `manifest.json`; `ga_py/` turns that
into a nanobind extension (`_ga_py`) exposing submodules `ega` / `pga` / `sta` plus the
top-level `hd::ga` free functions.

**Two venvs — never conflate:** `ga_bindgen/.venv` (scan/emit, needs `libclang`) vs
`ga_py/.venv` (build + tests, needs `pytest hypothesis numpy nanobind`).

**Regeneration chain** (after changing the C++ API; steps 2 and 5 need CMake configured
with `-D_GA_BUILD_PYTHON=ON`):

```bash
ga_bindgen/.venv/bin/python ga_bindgen/src/scan.py            # -> manifest.json
ga_bindgen/.venv/bin/python ga_bindgen/src/emit_nanobind.py   # -> ga_py/src/generated/*.cpp
cmake --build build --target _ga_py                           # rebuild the extension
PYTHONPATH=build/ga_py:ga_py/python ga_py/.venv/bin/python ga_bindgen/src/emit_stubs.py  # -> *.pyi
cmake --build build --target regenerate_python_test_data      # -> cross-check JSON
PYTHONPATH=build/ga_py:ga_py/python ga_py/.venv/bin/python -m pytest ga_py/tests/
```

**Hand-sync surface.** The generated layer (`register_all.cpp`, `bindings_*.cpp`) is
automatic, but adding a new algebra / type / function requires manual edits to these
NON-generated files — each one surfaced as a test failure when missed:

- `ga_py/src/module.cpp` — create the submodule and pass it to `register_all`
- `ga_py/python/ga_py/__init__.py` — import/forward the submodule and update `__all__`
- `ga_bindgen/src/emit_stubs.py` — add `<submod>.pyi` to `STUB_FILES`
- `ga_test/python_utilities/export_python_cases.cpp` — representative cross-check cases
- `ga_py/tests/conftest.py` — `FIELD_ORDER` for new struct-shaped types
- `ga_py/tests/test_cross_check.py` — `_TYPE_LOCATIONS` submodule scan loop
- `ga_py/tests/test_constants.py` — `EXPECTED` table **and** the total-count assertion
- `ga_py/tests/test_{coverage,scalars,grade_lookup}.py` — per-type parametrize lists
- `ga_py/tests/test_stubs.py` — `STUBS` dict and the `__format__` parametrize

Helpers used only by the C++ ops (not public API) belong in `hd::ga::detail`, which
`scan.py` excludes from the binding — see the STA rotor section below.

**Binding scope.** The generator binds user type-aliases (the GA value types), free
functions, operators, constants, and **pure-data structs** — concrete (non-template)
aggregates whose members are *all* public fields (the physics PODs `pose2dp`/`pose3dp`,
`kin_state2dp`/`kin_state3dp`). Stateful classes (`static_/kinematic_/dynamic_system*`)
are NOT bound by design — reproduce their behaviour in Python from the bound primitives
(see `ga_py/tests/test_merry_go_round.py`, which reconstructs `kinematic_system2dp`). A
new pure-data struct is auto-bound by a plain regeneration and does NOT trip the hand-sync
failures above — the coverage/cross-check lists don't assert completeness, and their
float-arity ctor pattern doesn't fit a field-typed ctor (`pose2dp(vec2dp, double)`), so
cover a new struct with a dedicated test instead. Give the C++ type an fmt formatter (in
`ga/detail/fmt/`) so the generated `__str__`/`__format__` (which call `fmt::format("{}",
v)`) work.

**Adding a FIELD to an already-bound pure-data struct** (e.g. new `joint_state{2,3}dp`
members) is a hand-sync trap even though it needs no manifest/scope edits: the generated
binding emits a single *all-fields, positional, no-defaults* constructor, so the arity grows
(`joint_state2dp` went 5 -> 8 args) and **every positional call site in `ga_py/tests/`
breaks at runtime** (`incompatible function arguments`). After regenerating, grep the tests
for the struct ctor and append the new trailing args (and assert the new fields round-trip).
If you only want C++ state that should NOT widen the Python ctor, store it OUTSIDE the bound
struct — e.g. `dynamic_system`'s applied-wrench / driven-joint specs live in side
`std::unordered_map`s on the (unbound) system, not as `joint_state` fields, precisely to
avoid this churn.

**libclang version resilience.** The pip `libclang` package is pinned ≤18.x while the
system LLVM (Homebrew) may be much newer. `scan.py` does NOT require a matching dylib — it
tolerates unknown libclang `CursorKind` ids (e.g. the C++20 parenthesized aggregate-init
expression, which appears in member initialisers like `pose(origin, angle)`), skipping
them since they only occur inside bodies the scanner never collects. A newer system LLVM
is fine; do NOT hardcode a versioned `libclang` path in `clang_setup.py` (keep it
system-agnostic across macOS/Windows).

## Congruence Testing (`is_congruent`)

Two GA elements are **congruent** if they span the same subspace up to a non-zero
scalar — scalar multiples, **regardless of sign or magnitude**.

**Algorithm (unified `A = k*B`)** — all graded overloads delegate to the shared
`hd::ga::detail::coeffs_congruent<N>(a, b, rel_tol)` in `ga/ga_value_t.hpp`:

1. Component-wise zero detection (robust for degenerate PGA metrics), using an absolute
   `eps` floor independent of `rel_tol`.
2. Fix the scale factor `k` from the **largest-magnitude** component of `b` (NOT the first
   non-zero): dividing by the dominant component keeps `k` well-conditioned, so a small
   pivot's rounding is never amplified onto the dominant entries (amplification factor
   `max_component / pivot_component` is unbounded for a first-non-zero pivot).
3. Verify `a = k*b` for all components within a **relative** tolerance.
4. Scalars and top-grade pseudoscalars: any two non-zero ones are congruent (own overloads,
   not via the helper).

**Two tolerance lessons (both essential):**

- **Relative**, not absolute — scale the per-component threshold by the operands
  (`rel_tol * max(|components|, 1)`), since absolute `eps` is meaningless across magnitudes.
- **The default is `eps_congruent` (1e-12), not machine `eps`.** Congruence is almost always
  asked of two *independently computed* quantities (results of wedge/dual/motor chains)
  whose difference has accumulated error far beyond a single rounding; machine `eps`
  (~1.11e-15) rejects geometrically-identical operands. `1e-12` is still ~12 orders below any
  real geometric signal. Pass a tighter/looser `tolerance` only for a specific reason.

**Notes:**

- Use type-specific overloads with `requires`, not generic templates — same-grade
  high-order wedges are undefined (e.g. `wdg(BiVec3d, BiVec3d)` → grade 4 in 3D).
- Provided per graded type in every algebra (EGA2D/3D, PGA2DP/3DP), as the grade
  structure allows.

## Inverses: `inv()` (geometric) vs `rinv()` (regressive)

Two distinct multiplicative inverses live in `*_ops_products.hpp`:

- **`inv(u)`** — inverse w.r.t. the geometric product: `u ⟑ inv(u) = inv(u) ⟑ u = 1`
  (scalar identity). Implemented (Hitzer–Sangwine 2016 formula for the multivector cases)
  for **every graded + multivector type in all five algebras**, with two metric caveats:
  the **PGA pseudoscalars are NOT invertible** (degenerate metric — intentionally
  omitted), but the **STA pseudoscalar IS** (`I_4ds² = -1`, so `inv(I_4ds) = -I_4ds`).
- **`rinv(u)`** — inverse w.r.t. the regressive geometric product `rgpr` (`⟇`):
  `u ⟇ rinv(u) = rinv(u) ⟇ u = Iₙ` (the pseudoscalar is the rgpr identity). Defined as the
  geometric inverse carried through the complement (dual) map:
  - odd-dim (pga2dp): `rinv(u) = cmpl(inv(cmpl(u)))`
  - even-dim (pga3dp): `rinv(u) = l_cmpl(inv(r_cmpl(u))) = r_cmpl(inv(l_cmpl(u)))`
    (both forms coincide)

  **`rinv()` exists for PGA2DP and PGA3DP ONLY** — `rgpr`/`⟇` is implemented only in PGA
  (regressive products are PGA-scoped by design; EGA/STA have no `rgpr`, so `rinv` there
  is neither defined nor testable). Its degeneracy is the **mirror image** of `inv()`'s:
  the **scalar has no `rinv`** (dual to the pseudoscalar having no `inv`), while the
  **pseudoscalar gains an `rinv`** (dual to the scalar's `inv`). `rinv()` is bound into
  `ga_py` (regenerate via the scan chain if its signature set changes).

## STA4D rotor operations (`ga_sta4ds_ops.hpp`)

STA transforms via the geometric-product rotor sandwich `X' = R ⟑ X ⟑ rev(R)` (like
ega3d, not pga3dp's regressive motors). Gotchas worth remembering:

- **Causal character & rotor branch read the geometric square, not `nrm_sq`.**
  `is_timelike`/`is_spacelike`/`is_lightlike` and the `exp` rotation-vs-boost split test
  the sign of `B² = gr0(X ⟑ X)` (the blade times itself), computed by
  `detail::sta4ds_geom_sq`. Under the exomorphism metric `nrm_sq` is the *reverse-norm*
  `⟨rev(X) X⟩₀` = the extended metric `P`, which is sign-flipped from `B²` at grades 2,3 —
  so `nrm_sq` must NOT drive causal character. (`detail::sta4ds_geom_sq` is a hand-written
  transcription of `gr0(X⟑X)`; a test pins it to the geometric product.) The `P`
  (metric) vs `B²` (geometric square) distinction is the key subtlety here.

- **Versor norm vs grade-wise `nrm_sq`.** `sqrt(rotor)` / normalizing a rotor
  (`MVec4ds_E`) must use the **versor norm** `sqrt(gr0(rev(X) ⟑ X))`, *not*
  `nrm_sq(MVec4ds_E)`. The latter is the grade-wise sum with the signed Lorentzian
  bivector norm, which for a unit rotor evaluates to `cos(2a)` (rotation) or
  `cosh(2φ)` (boost) — **not 1**. `sqrt(R) = X / sqrt(gr0(rev(X)·X))` with `X = 1 + R`
  is exact for simple rotors of either kind; guard the `R = -1` (2π) degenerate case.

- **Closed-form `transform_opt` helpers live in `hd::ga::detail`.** The per-grade
  sandwich matrices (4×4 for vec **and** trivec — identical in STA; 6×6 for bivec) are
  built once from the rotor coefficients and shared by the scalar and `std::vector`
  batch overloads. Putting the builder helpers in `hd::ga::detail` keeps them out of the
  bindgen scan **and** avoids creating a `hd::ga::sta::detail` namespace, which would
  shadow the `detail::safe_epsilon` etc. used unqualified inside `hd::ga::sta`. Scalar
  `transform_opt` is slower than `transform()` for one-offs (matrix build not amortised);
  only the batch overload wins (~3× at `-O3`). Bench: `ga_test/utilities/`.

- **Planned physics tier**: an STA differential-operator/electrodynamics layer (vector
  derivative ∇, Maxwell `∇F = J`, Lorentz force + proper-time rotor EOM) is planned;
  `ga/detail/ga_stencil.hpp` (FD stencil weights) is its numeric groundwork.

- **`get_rotor` and `get_boost` use OPPOSITE half-angle signs by design.** `get_rotor`
  negates (`-theta/2`, matching ega3d); `get_boost` does not (`+phi/2`). Each sign is the
  one that makes a positive parameter a positive-sense active transform (a `+theta`
  right-handed rotation; a boost with `gamma = cosh(phi)`). This is NOT an inconsistency —
  "unifying" the two silently reverses one transform. Called out at both function comments
  and in `ga_docu` ("Motion in STA"). (`exp`/`log` are sign-consistent — the `+/-` split is
  only in these two convenience builders.)

## Rigid-body dynamics (`dynamic_system{2,3}dp`, `ga_pga{2,3}dp_ops_mechanics.hpp`)

**Three kinds of force element** feed the joint-space generalised force `tau` (all additive;
with none attached the gravity/bias path is byte-unchanged):

- `set_joint_spring_damper(idx, k, c, q0)` — a linear spring + damper on a joint
  *coordinate*: `tau_j += -k(q-q0) - c q̇` (torsional for revolute, linear for prismatic).
- `set_applied_wrench(idx, fn)` — an external wrench given as a function of *time*
  `fn(t) -> bivecNdp` (world frame), projected onto the supporting joint screws.
- `add_grounded_spring(idx, anchor_b, [p0_world,] k, c)` — a *configuration-dependent*
  spring/damper tying a body-fixed point `anchor_b` to a ground anchor `p0_world` with
  anisotropic world-axis stiffness `k` + isotropic damping `c`. Its restoring wrench
  `wdg(P, F)` is recomputed from the live pose/velocity each RK4 sub-step (NOT a function of
  time). Because it acts at a physical point, ONE spring yields both a translational
  stiffness AND, via its lever arm, a tilt stiffness `k·l²` (e.g. two radial springs at
  axial `±l` give the emergent tilt stiffness). Contributes its potential to
  `potential_energy()`.

**Application-specific force elements via subclassing.** `dynamic_system{2,3}dp` exposes a
`protected virtual extra_wrenches()` (returns `{}` in the base); `assemble_mass_bias` folds
each returned `(frame, world-wrench)` pair onto the frame's supporting joints exactly like
the grounded-spring path. A subclass overrides it to inject configuration-dependent wrenches
(e.g. a contact/penalty model) without the base knowing anything about them — the generic
base stays free of any application vocabulary.

**`assemble_mass_bias` moving-base invariant (do NOT regress).** The assembler sums each
body's inertia over **dof joints ∪ kinematically-driven joints**, projecting every
inertia-bearing body onto its *ancestor* dof joints. A driven joint is a moving base: its
body inertia loads the joints above it and its prescribed velocity feeds their Newton-Euler
bias (so a driven spinning rotor produces the clamped gyroscopic/centrifugal dynamics).
Refactoring the body list back to dof-joints-only **silently drops driven-joint inertia** →
a singular mass matrix when a driven joint carries the only inertia below a dof joint (e.g. a
motor-clamped rotor spin). Only indirectly tested (via a driven spinning rotor), so guard
it deliberately.

## Geometric Algebra Mathematical Foundations

### Core GA Conventions

- **G(p,n,z)**: p positive, n negative, z null (degenerate) basis vectors.
- **Grades**: scalar 0, vector 1, bivector 2, trivector 3, … pseudoscalar n.
- **Products**: geometric (full, carries metric), wedge `^` (antisymmetric, 0 on
  repeats), dot (symmetric, extended-metric contraction).
- Basis is 1-based (`e1,e2,…` / `g1,…`); the rule generator handles 0-based too.

### PGA motion: velocity twists are dimension-dependent (do NOT assume bivector)

The motor rate / velocity "twist" `Omega = 2 * Mdot ⟇ rrev(M)` (the `exp` generator of a
motor) has a grade that **depends on the dimension** — do not assume it is always a
bivector (a common, wrong assumption carried over from 3D):

- **pga2dp**: motors are odd-grade (`MVec2dp_U`), so `Mdot ⟇ rrev(M)` is odd → the twist
  `Omega` is a **vector** (`Vec2dp`); `exp(Vec2dp)` builds the motor. (The kinematics
  `twist2dp` alias is a `vec2dp` for exactly this reason.)
- **pga3dp**: motors are even-grade (`MVec3dp_E`), so the twist `Omega` is a **bivector**
  (`BiVec3dp`); `exp(BiVec3dp)` builds the motor.

(This mirrors the odd/even-dimensional split already noted for complements.)

The PGA rate-of-change of a point is `Xdot = rcmt(Omega, X)` — **argument order matters**:
`rcmt(Omega, X) = -rcmt(X, Omega)` (twist first); the sign/order differs from EGA's
`cmt(r, Omega_E)`. Moving-frame kinematic fields (see the `ga_docu` PDF, "Moving coordinate
systems"): velocity `rcmt(Omega, X)`, centripetal `rcmt(Omega, rcmt(Omega, r))`, Coriolis
`2*rcmt(Omega, v_rel)`, frame/Euler `rcmt(Omega_dot, r)`. Full derivations in the `ga_docu`
PDF (modelling-motion chapter).

### Supported Algebra Types

**EGA2D - G(2,0,0)**: Euclidean 2D

- Basis: `{1, e1, e2, e12}`
- Metric: `{+1, +1}`
- Extended Metric: `{1, 1, 1, 1}`

**EGA3D - G(3,0,0)**: Euclidean 3D

- Basis: `{1, e1, e2, e3, e23, e31, e12, e123}`
- Metric: `{+1, +1, +1}`
- Extended Metric: `{1, 1, 1, 1, 1, 1, 1, 1}`

**PGA2DP - G(2,0,1)**: Projective 2D

- Basis: `{1, e1, e2, e3, e31, e32, e12, e321}`
- Metric: `{+1, +1, 0}`
- Extended Metric: `{1, 1, 1, 0, 0, 0, 1, 0}`

**PGA3DP - G(3,0,1)**: Projective 3D

- Basis: `{1, e1, e2, e3, e4, e41, e42, e43, e23, e31, e12, e423, e431, e412, e321, e1234}`
- Metric: `{+1, +1, +1, 0}`
- Extended Metric: `{1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0}`

**STA4D - G(1,3,0)**: Space-Time Algebra (library type `sta4ds`, `*_4ds` types)

- Basis: `{1, g1, g2, g3, g4, g14, g24, g34, g23, g31, g12, g234, g314, g124, g123, g1234}`
- Metric: `{-1, -1, -1, +1}` (space-like g1,g2,g3; time-like g4 with g4² = +1) —
  i.e. G(1,3,0): 1 positive (g4), 3 negative (g1,g2,g3)
- Extended Metric: `{1, -1, -1, -1, +1, -1, -1, -1, +1, +1, +1, +1, +1, +1, -1, -1}`
  (the wedge exomorphism `P = ∏ gᵢ`, same construction as EGA/PGA — see below)

### Extended Metric Calculation Rules

**Standard Algebras (EGA, PGA):**

- **Vectors**: Use base metric signature directly
- **Higher grades**: Product of constituent vector metrics
- **Pseudoscalar**: Determinant of metric tensor
- **Mixed signatures**: Apply conforming property G(a ∧ b) = G(a) ∧ G(b)

**Space-Time Algebra (STA4D):** NO special rules — STA follows the **Standard Algebras**
rule above. The extended metric is the wedge exomorphism `P(e_S) = ∏_{i∈S} g_i` (product
of the constituent vector metrics), seeded by the signature `{-1,-1,-1,+1}`:

- **Vectors**: g1, g2, g3 → -1; g4 → +1
- **Bivectors**: with g4 (g14, g24, g34) → -1; without g4 (g23, g31, g12) → +1
- **Trivectors**: with g4 (g234, g314, g124) → +1; without g4 (g123) → -1
- **Pseudoscalar**: determinant = -1

CRITICAL: these are the *metric exomorphism* values `P`, NOT the blade squares
`e_S ⟑ e_S` (which carry the reversion sign and are sign-flipped at grades 2,3). An
earlier `is_minkowski` override stored the blade squares in the metric slot — that was a
bug, now removed. Do NOT reintroduce "special" STA metric rules. The blade square (the
**geometric square** `B² = gr0(X⟑X)`) is a separate quantity used only for causal
character / rotors — see "STA4D rotor operations" above.

### Automatic Rule Generation System

The rule generator builds, per algebra from its `AlgebraConfig` (basis, metric
signature, canonical ordering):

- geometric-product rules via `multiply_basis_elements()` (sign = #swaps to
  canonical order),
- the extended (dot) metric for every basis blade,
- complement and dual rule tables.

Output is validated character-identical against the reference for all five
algebras. `ga_prdxpr_rule_generator_test` prints these tables — the source of
truth for the hand-coded complement/dual functions (guard them with the
transcription gate above).
