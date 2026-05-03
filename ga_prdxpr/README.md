# Copyright

Copyright 2024-2026, Daniel Hug. All rights reserved.

# Configurable GA Product Expression Generator

A standalone, configurable geometric algebra product expression generator that produces
character-identical output to the reference implementation while providing a clean,
maintainable configuration-driven architecture.

## Features

- **Complete Coverage**: Five geometric algebras — EGA2D, EGA3D, PGA2DP, PGA3DP (full),
  STA4D (skeleton, see *Open Codegen Work*)
- **All Product Families**: geometric, wedge, contraction, expansion, regressive variants,
  sandwich
- **Character-Identical Reference Output**: the bare invocation produces a stable
  11559-line searchable reference
- **C++ Code Generation**: `--output=code` emits ready-to-paste primitive product
  implementations matching `ga/*_ops_products.hpp` byte-for-byte (~95% of cases; the rest
  are intentional semantic delegations)
- **Configuration-Driven**: per-algebra configs in
  [src_prdxpr/algebras/](src_prdxpr/algebras/) are the user-editable surface; the
  generator engine is reused unchanged across algebras
- **Standalone**: only depends on fmt

## Supported Algebras

### EGA2D (Euclidean 2D)

- 10 products: gpr, cmt, wdg, dot, rwdg, l_contract, r_contract, l_expand, r_expand,
  rotate

### EGA3D (Euclidean 3D)

- 10 products: gpr, cmt, wdg, dot, rwdg, l_contract, r_contract, l_expand, r_expand,
  rotate

### PGA2DP (Projective 2D+)

- 16 products: gpr, cmt, wdg, dot, rwdg, rdot, rgpr, bulk/weight contractions/expansions,
  move2dp

### PGA3DP (Projective 3D+)

- 16 products: gpr, cmt, wdg, dot, rwdg, rdot, rgpr, bulk/weight contractions/expansions,
  move3dp

## Building

### Prerequisites

- CMake 3.28 or higher
- C++23 compatible compiler
- fmt library

### Build Instructions

```bash
# Install fmt library (macOS)
brew install fmt

# Build (from project root)
mkdir -p build && cd build
cmake ..
cmake --build . --target ga_prdxpr ga_prdxpr_rule_generator_test

# Run (executable lives in build/ga_prdxpr/)
cd ga_prdxpr
./ga_prdxpr > output.txt
```

### Verification

The bare invocation must produce a stable 11559-line reference:

```bash
./ga_prdxpr | wc -l   # 11559
```

A round-trip byte-identity check against the C++ source (`ga/*_ops_products.hpp`) is
documented under [Validating against existing
source](#validating-against-existing-source).

## Project Layout

The generator is split into clearly-scoped subfolders. The only folder a user
typically edits to add a new product or tweak an existing one is
[src_prdxpr/algebras/](src_prdxpr/algebras/).

```text
src_prdxpr/
├── ga_prdxpr_main.cpp          # entry point — composes per-algebra ProductConfig
│                                 lists and calls the generator
├── cli/                         # argument parser (--algebra, --product, --output)
├── core/                        # foundational types (mvec_coeff, prd_rules,
│                                 ProductConfig, OutputCase, AlgebraData) +
│                                 table operations and printing helpers
├── rules/                       # rule-generation engine + metric/dual math
│                                 (basis multiplication, Gram matrices, complement
│                                 and dual rules, metric export to ga_usr_consts.hpp)
├── algebras/                    # ★ user-editable: per-algebra configuration
│   ├── ga_prdxpr_<alg>.{hpp,cpp}        # basis vectors, coefficient strings,
│   │                                       rule-set instantiation at startup
│   └── ga_prdxpr_<alg>_config.{hpp,cpp} # one get_<alg>_<prd>_config() per product;
│                                          composes OutputCases that the generator
│                                          turns into output
├── sandwich/                    # sandwich-product transformation pipeline
│                                 (parser → simplifier → n-ary AST → transformer)
├── generator/                   # ConfigurableGenerator orchestrator: dispatches a
│                                 ProductConfig to coeffs/tables/code output
├── codegen/                     # --output=code emitter (TypeRegistry + emit_function)
│   └── tools/                   # Python validation pipeline (diff_codegen.py, ...)
├── tests/                       # standalone ga_prdxpr_rule_generator_test executable
└── reference_output/            # golden .txt fixtures used during development
```

### Where to make changes

| You want to...                                          | Edit in...                                                                                                                               |
| ------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| Add a new product to an existing algebra                | `algebras/ga_prdxpr_<alg>_config.cpp` (new `get_*_config()`) + register in `ga_prdxpr_main.cpp`                                           |
| Add a new coefficient pattern (e.g. `R_even`, `svBps`)  | `algebras/ga_prdxpr_<alg>.hpp`                                                                                                           |
| Add a brand-new algebra                                 | mirror the pair `algebras/ga_prdxpr_sta4d{,_config}.{hpp,cpp}`, then register it in `ga_prdxpr_main.cpp` and `CMakeLists.txt`            |
| Tweak how products are *rendered* (formatting, headers) | `generator/ga_prdxpr_generator.cpp`                                                                                                      |
| Adjust the C++ codegen output style                     | `codegen/ga_codegen_emitter.cpp` (+ `codegen/ga_codegen_types.cpp` for new types)                                                        |
| Change the rule-generation math                         | `rules/ga_prdxpr_rule_generator.cpp`                                                                                                     |

### Configuration system

Each algebra in [src_prdxpr/algebras/](src_prdxpr/algebras/) defines:

- Multivector basis (e.g. `mv2d_basis = {"1","e1","e2","e12"}`) and metric signature
- Coefficient string maps (`mv2d_coeff_A`, `mv3d_coeff_R_even`, `mv3dp_coeff_svBps`, ...)
- A list of `ProductConfig`s (one per product variant) built from `OutputCase`s that
  pair a coefficient pattern with a result-type filter (`mv`, `vec`, `bivec`, ...)

`ProductConfig` and `OutputCase` are defined in
[src_prdxpr/core/ga_prdxpr_config_types.hpp](src_prdxpr/core/ga_prdxpr_config_types.hpp);
`AlgebraConfig` (used by the rule-generation engine) is defined in
[src_prdxpr/rules/ga_prdxpr_rule_generator.hpp](src_prdxpr/rules/ga_prdxpr_rule_generator.hpp).

### Key behaviors

1. **Sandwich Products**: two-step `M ⟑ X ⟑ rev(M)` with parenthesization preserved
   between steps via `brace_switch::use_braces`
2. **Brace Switch**: automatic parentheses on intermediate expressions so that the second
   step parses correctly
3. **Basis Tables**: generated from `AlgebraConfig` via `generate_algebra_rules()` — no
   hand-maintained tables
4. **Separator Logic**: section dividers are emitted only when human-readable output is
   requested (suppressed in `--output=code` to keep the output clang-format-clean)

## Mathematical Accuracy

The generator preserves all critical mathematical insights from the reference
implementation:

- **PGA Regressive Products**: Uses correct complement transformation sequences
- **Bulk/Weight Operations**: Implements proper dual transformation order
- **Sandwich Products**: Maintains two-step pattern with correct coefficient usage
- **Mathematical Identities**: All GA identities preserved through exact table generation

## Usage Example

How `ga_prdxpr_main.cpp` composes a per-algebra product list and hands it to the
generator (excerpt for PGA2DP):

```cpp
auto pga2dp_algebra = create_pga2dp_algebra_data();
std::vector<ProductConfig> pga2dp_configs = {
    get_pga2dp_gpr_config(),
    get_pga2dp_cmt_config(),
    // ... more configs ...
    get_pga2dp_sandwich_rgpr_config(),
};
generate_algebra_products(generator, pga2dp_configs, pga2dp_algebra,
                          /*is_last_algebra=*/false, options);
```

The `get_*_config()` functions live in
[src_prdxpr/algebras/ga_prdxpr_pga2dp_config.cpp](src_prdxpr/algebras/ga_prdxpr_pga2dp_config.cpp).

## C++ Code Generation (`--output=code`)

Beyond coefficient expressions and basis tables, ga_prdxpr can emit
ready-to-paste C++ function implementations for primitive products.
Opt-in only — never appears in the default output:

```bash
ga_prdxpr --algebra=ega2d --output=code                  # one algebra
ga_prdxpr --algebra=pga3dp --products=gpr --output=code  # one product
ga_prdxpr --output=all                                   # everything
```

### Invocation matrix

| Invocation             | coeffs | tables | metrics | code |
| ---------------------- | ------ | ------ | ------- | ---- |
| `ga_prdxpr` (no args)  | ✓      | ✓      | ✓       | —    |
| `--output=all`         | ✓      | ✓      | ✓       | ✓    |
| `--output=code`        | —      | —      | —       | ✓    |
| `--output=coeffs,code` | ✓      | —      | —       | ✓    |

The bare invocation produces the canonical 11559-line searchable
reference (coeffs + tables + metrics). `code` is opt-in only so the
default output stays stable; ask for it explicitly when you want C++
implementations.

### What is emitted

Flat-constructor returns with temp-var aliasing (`ctype const cN = ...; return
Type<ctype>(c0, c1, ...)`), east-const, and the `[[maybe_unused]]` form for zero-result
cases.

**Currently covered:** gpr, wdg, dot, cmt, rgpr, rwdg, rdot, rcmt, twdg1, rtwdg1,
left/right contraction (`<<` / `>>`).

**Skipped:**

- Sandwich products — separate two-step structure, hand-rolled per-algebra in
  `generate_sandwich_case`. Codegen-side handling deferred (see *Open Codegen Work*
  below).
- Wrapper-only product families (PGA `l_bulk_contract` and friends, EGA
  `l_expand`/`r_expand`) — they have no per-case C++ overloads.
- Display-only "alternative" product variants.

### Codegen architecture

- `src_prdxpr/codegen/ga_codegen_types.{hpp,cpp}`: `TypeRegistry` mapping filter names
  (`vec`, `bivec`, `mv_e`, ...) to C++ types, storage kinds (Indexed / Named / Composite /
  SingleValue), and basis-index layouts. One per algebra (ega2d, ega3d, pga2dp, pga3dp,
  sta4d).
- `src_prdxpr/codegen/ga_codegen_emitter.{hpp,cpp}`: `emit_function()` builds a complete
  C++ overload from one `OutputCase`. Detects zero results (explicit `-> 0` marker, `-> 0
  X` typed override, or all-zero computed coefficients), wraps bare Scalar/PScalar
  operands with `ctype(...)`, and dispatches to body builders by storage kind.
- `src_prdxpr/codegen/tools/`: validation pipeline (`diff_codegen.py`) and one-shot source
  transformers (`inline_to_tempvars.py`, `expand_delegations.py`). See `tools/README.md`
  for usage details.

### Validating against existing source

`ga/*_ops_products.hpp` are 95% char-identical (598/630) with codegen output
post-`clang-format`. The remaining ~5% are intentional semantic delegations (e.g.
`gpr(v,v) = MVec_E(dot(v,v), wdg(v,v))`) which encode GA identities and are kept by
design.

```bash
# round-trip check for one algebra (run from build/)
./ga_prdxpr/ga_prdxpr --algebra=ega3d --output=code | \
    clang-format -assume-filename=x.cpp > /tmp/codegen.cpp
clang-format ../ga/ga_ega3d_ops_products.hpp > /tmp/source.cpp
python3 ../ga_prdxpr/src_prdxpr/codegen/tools/diff_codegen.py \
    /tmp/source.cpp /tmp/codegen.cpp
```

## Open Codegen Work

### STA4D rollout

The codegen plumbing for STA4D is in place — there is a `TypeRegistry` entry using
`Scalar4d` / `Vec4d` / `BiVec4d` / `TriVec4d` / `PScalar4d` / `MVec4d_E` / `MVec4d_U` /
`MVec4d` (see
[src_prdxpr/codegen/ga_codegen_types.cpp](src_prdxpr/codegen/ga_codegen_types.cpp)
`build_sta4d`). What's still missing is the per-product *content*:

1. Populate the `.cases = {}` arrays in
   [src_prdxpr/algebras/ga_prdxpr_sta4d_config.cpp](src_prdxpr/algebras/ga_prdxpr_sta4d_config.cpp).
   A reasonable starting point is to mirror the PGA3DP structure (same 4D dimensionality,
   similar product family) and drop the bulk_*/weight_* contractions which are
   PGA-specific.
2. (Optional) If the future `ga/ga_sta4d_ops_products.hpp` uses a different type-name
   convention (e.g. `ScalarSta4d`), edit the `cpp_type` strings in `build_sta4d()` —
   that's the only place that pins names.
3. Run `ga_prdxpr --algebra=sta4d --output=code` and pipe the output to a fresh
   `ga/ga_sta4d_ops_products.hpp` (and the basics/only headers as needed by the existing
   per-algebra pattern).
4. Wire up the test build under `ga_test/`. Note that `ga_test/src/ga_sta4d_test.hpp`
   already exists but `#include`s `ga/ga_sta.hpp` which doesn't exist yet — that include
   is the placeholder for the new STA4D library headers.

### Sandwich product codegen

`generate_sandwich_case` in
[src_prdxpr/generator/ga_prdxpr_generator.cpp](src_prdxpr/generator/ga_prdxpr_generator.cpp)
is hand-rolled per algebra (vector + bivector for ega3d/pga2dp; +trivector for
pga3dp/sta4d). Code generation for sandwich products needs its own emitter that handles
the two-step `M ⟑ X ⟑ rev(M)` structure and the substantial simplification that source
applies to the second step. For the simplified target form see
`ga_pga2dp_ops_physics.hpp`'s `move2dp(X, M)` and `ga_pga3dp_ops_physics.hpp`'s
`move3dp(X, M)`.

## Extending the System

### Adding a new product type

1. In [src_prdxpr/algebras/](src_prdxpr/algebras/), add a new `get_<alg>_<prd>_config()`
   function in `ga_prdxpr_<alg>_config.cpp` (and declare it in the matching
   `_config.hpp`). It returns a `ProductConfig` with one `OutputCase` per result-type
   filter you want to emit.
2. Append the new config to the product list for that algebra in
   [src_prdxpr/ga_prdxpr_main.cpp](src_prdxpr/ga_prdxpr_main.cpp).
3. If the product needs a different basis table than gpr/wdg/dot, teach
   `generator/ga_prdxpr_generator.cpp::get_basis_table_for_product()` how to build it.
4. Build, run, and compare output against expected mathematical
   identities.

### Adding a new algebra

1. Mirror the file pair `algebras/ga_prdxpr_sta4d.{hpp,cpp}` and
   `algebras/ga_prdxpr_sta4d_config.{hpp,cpp}` — the header declares basis vectors /
   coefficients / extern rules, the implementation triggers `generate_algebra_rules()` at
   static-init time, and the `_config.cpp` exposes `get_<alg>_algebra_config()` plus the
   per-product `ProductConfig` builders.
2. Add a `TypeRegistry` entry for the new algebra in
   [src_prdxpr/codegen/ga_codegen_types.cpp](src_prdxpr/codegen/ga_codegen_types.cpp) if
   you want `--output=code` support.
3. Register the algebra block in `ga_prdxpr_main.cpp` (mirroring the existing
   EGA2D/EGA3D/... blocks).
4. Add the new files to `CMakeLists.txt` — note that `*_config.cpp` **must** be listed
   before `*.cpp` in `SOURCES` so that the internal-linkage basis globals defined per-TU
   are dynamically initialized before the rule-instantiation TU's static init runs.

## Verification and Testing

The bare invocation must remain 11559-line stable; deviations indicate a regression. The
C++ code generator is independently validated by piping `--output=code` through
`clang-format` and diffing against the matching `ga/*_ops_products.hpp` (see [Validating
against existing source](#validating-against-existing-source)).
