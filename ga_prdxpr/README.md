# Configurable GA Product Expression Generator

A standalone, configurable geometric algebra product expression generator that produces character-identical output to the reference implementation while providing a clean, maintainable configuration-driven architecture.

## Features

- **Complete Coverage**: Supports all 4 geometric algebras (EGA2D, EGA3D, PGA2DP, PGA3DP)
- **52 Total Products**: All product types including geometric, wedge, contractions, expansions, and sandwich products
- **Character-Identical Output**: 100% match with reference implementation (6928 lines)
- **Standalone**: No external dependencies except fmt library
- **Configuration-Driven**: Clean separation between mathematical operations and output formatting
- **Maintainable**: Easy to add new products or modify existing ones

## Supported Algebras

### EGA2D (Euclidean 2D)

- 10 products: gpr, cmt, wdg, dot, rwdg, lcontract, rcontract, lexpand, rexpand, rotor_gpr

### EGA3D (Euclidean 3D)

- 10 products: gpr, cmt, wdg, dot, rwdg, lcontract, rcontract, lexpand, rexpand, rotor_gpr

### PGA2DP (Projective 2D+)

- 16 products: gpr, cmt, wdg, dot, rwdg, rdot, rgpr, bulk/weight contractions/expansions, regressive_sandwich

### PGA3DP (Projective 3D+)

- 16 products: gpr, cmt, wdg, dot, rwdg, rdot, rgpr, bulk/weight contractions/expansions, regressive_sandwich

## Building

### Prerequisites

- CMake 3.28 or higher
- C++23 compatible compiler
- fmt library

### Build Instructions

```bash
# Install fmt library (macOS)
brew install fmt

# Build
mkdir build
cd build
cmake ..
cmake --build .

# Run
./ga_prdxpr_configurable > output.txt
```

### Verification

The output should be exactly 6928 lines and match the reference implementation character-for-character:

```bash
wc -l output.txt  # Should show: 6928 output.txt
```

## Architecture

### Core Components

- **`ga_prdxpr_configurable_main.cpp`**: Clean main function with helper functions
- **`ga_prdxpr_configurable_generator.cpp`**: Core generation logic
- **Configuration Files**: Separate config files for each algebra (`*_config.cpp/hpp`)
- **Mathematical Libraries**: Reference mathematical implementations (`ga_prdxpr_*.cpp/hpp`)

### Configuration System

Each algebra has its own configuration defining:

- Product types and their mathematical basis tables
- Coefficient mappings (e.g., `mv2d_coeff_R_even`, `mv3dp_coeff_svBps`)
- Filter specifications for different multivector types
- Sandwich product settings with brace handling

### Key Features

1. **Sandwich Products**: Two-step mathematical operations with proper parenthesization
2. **Brace Switch**: Automatic parentheses handling for complex intermediate expressions
3. **Basis Tables**: Mathematical product tables generated using proven transformation rules
4. **Separator Logic**: Exact output formatting matching reference implementation

## Mathematical Accuracy

The generator preserves all critical mathematical insights from the reference implementation:

- **PGA Regressive Products**: Uses correct complement transformation sequences
- **Bulk/Weight Operations**: Implements proper dual transformation order
- **Sandwich Products**: Maintains two-step pattern with correct coefficient usage
- **Mathematical Identities**: All GA identities preserved through exact table generation

## Usage Example

```cpp
// Generate all products for PGA2DP
auto pga2dp_algebra = create_pga2dp_algebra_data();
std::vector<ProductConfig> pga2dp_configs = {
    get_pga2dp_gpr_config(),
    get_pga2dp_cmt_config(),
    // ... more configs
    get_pga2dp_regressive_sandwich_config()
};
generate_algebra_products(generator, pga2dp_configs, pga2dp_algebra, false);
```

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

Flat-constructor returns with temp-var aliasing
(`ctype const cN = ...; return Type<ctype>(c0, c1, ...)`), east-const,
and the `[[maybe_unused]]` form for zero-result cases.

**Currently covered:**
gpr, wdg, dot, cmt, rgpr, rwdg, rdot, rcmt, twdg1, rtwdg1,
left/right contraction (`<<` / `>>`).

**Skipped:**

- Sandwich products — separate two-step structure, hand-rolled per-algebra
  in `generate_sandwich_case`. Codegen-side handling deferred (see
  *Open Codegen Work* below).
- Wrapper-only product families (PGA `l_bulk_contract` and friends, EGA
  `l_expand`/`r_expand`) — they have no per-case C++ overloads.
- Display-only "alternative" product variants.

### Codegen architecture

- `src_prdxpr/codegen/ga_codegen_types.{hpp,cpp}`: `TypeRegistry` mapping
  filter names (`vec`, `bivec`, `mv_e`, ...) to C++ types, storage kinds
  (Indexed / Named / Composite / SingleValue), and basis-index layouts.
  One per algebra; STA4D entry not yet added (see *Open Codegen Work*).
- `src_prdxpr/codegen/ga_codegen_emitter.{hpp,cpp}`: `emit_function()`
  builds a complete C++ overload from one `OutputCase`. Detects zero
  results (explicit `-> 0` marker, `-> 0 X` typed override, or all-zero
  computed coefficients), wraps bare Scalar/PScalar operands with
  `ctype(...)`, and dispatches to body builders by storage kind.
- `src_prdxpr/codegen/tools/`: validation pipeline (`diff_codegen.py`)
  and one-shot source transformers (`inline_to_tempvars.py`,
  `expand_delegations.py`). See `tools/README.md` for usage details.

### Validating against existing source

`ga/*_ops_products.hpp` are 95% char-identical (598/630) with codegen
output post-`clang-format`. The remaining ~5% are intentional semantic
delegations (e.g. `gpr(v,v) = MVec_E(dot(v,v), wdg(v,v))`) which encode
GA identities and are kept by design.

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

`get_sta4d_*_config()` functions exist in `ga_prdxpr_sta4d_config.cpp`
but their `.cases = {}` arrays are empty. To generate STA4D code:

1. Populate the `.cases` arrays. A reasonable starting point is to
   mirror the PGA3DP structure (same 4D dimensionality, similar
   product family) and drop the bulk_*/weight_* contractions which are
   PGA-specific.
2. Add an STA4D entry to the `TypeRegistry` in
   `src_prdxpr/codegen/ga_codegen_types.cpp` — decide on the type-name
   convention (`Scalar4d` vs `ScalarSta4d`, etc.).
3. Run `ga_prdxpr --algebra=sta4d --output=code` and pipe to a fresh
   `ga/ga_sta4d_ops_products.hpp` (and the basics/only headers as
   needed by the existing per-algebra pattern).
4. Wire up the test build under `ga_test/`. Note that
   `ga_test/src/ga_sta4d_test.hpp` already exists but `#include`s
   `ga/ga_sta.hpp` which doesn't exist yet — that include is the
   placeholder for the new STA4D library headers.

### Sandwich product codegen

`generate_sandwich_case` in `ga_prdxpr_generator.cpp` is hand-rolled
per algebra (vector + bivector for ega3d/pga2dp; +trivector for
pga3dp/sta4d). Code generation for sandwich products needs its own
emitter that handles the two-step `M ⟑ X ⟑ rev(M)` structure and the
substantial simplification that source applies to the second step.
For the simplified target form see `ga_pga2dp_ops_physics.hpp`'s
`move2dp(X, M)` and `ga_pga3dp_ops_physics.hpp`'s `move3dp(X, M)`.

## File Structure

```text
ga_prdxpr_configurable/
├── CMakeLists.txt                        # Standalone build configuration
├── README.md                             # This file
├── ga_prdxpr_configurable_main.cpp       # Main application
├── ga_prdxpr_configurable_generator.*    # Core generator
├── ga_prdxpr_config_types.hpp            # Configuration data structures
├── ga_prdxpr_*_config.*                  # Algebra-specific configurations
└── ga_prdxpr_*.*                         # Mathematical reference implementations
```

## Extending the System

### Adding a New Product Type

1. Add configuration function in appropriate `*_config.cpp`
2. Add to the product list in `main.cpp`
3. Implement basis table generation if needed
4. Test output matches expected mathematical results

### Adding a New Algebra

1. Create new `ga_prdxpr_newAlgebra_config.*` files
2. Implement `create_newAlgebra_algebra_data()` function
3. Add mathematical reference implementation files
4. Add to main function with product configurations

## Verification and Testing

The system maintains 100% character-identical output with the reference implementation, ensuring mathematical correctness. Every generated expression has been verified against known geometric algebra identities.

## Copyright

Copyright 2024-2025, Daniel Hug. All rights reserved.
