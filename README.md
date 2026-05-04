# Overview on modules

- **ga**: header-only library for geometric algebra in 2d and 3d
- **ga_test**: tests for the ga header-only library
- **ga_lua**: lua scripting interface based on the ga library
- **ga_view**: qt-based visualization of 2d ga objects and operations
- **ga_prdxpr**: code generator for GA product expressions
- **ga_py**: Python wrapper around the `ga/` library — see
  [`ga_py/README.md`](ga_py/README.md) for installation and usage
- **ga_bindgen**: contributor tool that scans `ga/` headers with libclang and emits the
  `ga_py` bindings — see [`ga_bindgen/README.md`](ga_bindgen/README.md)

## Dependencies

- ga: fmt (must be included by every user — header-only lib!)
- ga_test: fmt, doctest
- ga_lua: fmt, lua, sol2
- ga_view: fmt, qt6
- ga_prdxpr: fmt
- ga_py: fmt, Python 3.10+, Python packages `nanobind` (build) and
  `pytest`/`hypothesis`/`numpy` (tests). End users do **not** need libclang.
- ga_bindgen: Python 3.10+, libclang (provided by a full LLVM install; the pip `libclang`
  package alone is insufficient on macOS because its resource directory is missing). Only
  contributors regenerating bindings need this.

## Installing dependencies on macOS

```bash
brew install fmt         # always required
brew install doctest     # only required for ga_test
brew install qt6         # only required for ga_view
brew install lua         # lua and sol2 only required for ga_lua
# sol2 to be installed from: https://github.com/ThePhD/sol2
brew install python      # only required for ga_py (>= 3.10)
brew install llvm        # only required for ga_bindgen (provides usable libclang)
```

For `ga_py` (Python wrapper) install and venv setup, see [`ga_py/README.md
§3`](ga_py/README.md#3-quick-start--running-your-first-python-script).

## Installing dependencies on Windows

```text
fmt:      https://github.com/fmtlib/fmt
doctest:  https://github.com/doctest/doctest
lua:      winget install --id DEVCOM.Lua   (need an x64 lua binary or linking will fail)
qt6:      official installer
sol2:     https://github.com/ThePhD/sol2
python:   winget install --id Python.Python.3.12   (only for ga_py; any 3.10+ works)
LLVM:     https://releases.llvm.org/   (only for ga_bindgen)
```

For `ga_py` (Python wrapper) install and venv setup, see [`ga_py/README.md
§3`](ga_py/README.md#3-quick-start--running-your-first-python-script).

## Usage of the ga library

0) Enclose the headers

    ```cpp
    "#include ga/ega.hpp"
    ```

    or

    ```cpp
    "#include ga/pga.hpp"
    ```

    as needed.

1) Make all available that is needed for GA, e.g. by

    ```cpp
    using namespace hd::ga;
    using namespace hd::ga::ega;
    ```

    or

    ```cpp
    using namespace hd::ga;
    using namespace hd::ga::pga;
    ```

2) Create a try catch block to catch any exceptions thrown by GA code

    ```cpp
    try {

        // call GA stuff here

    } catch (std::exception const& e) {

        // whichever alternative suits your needs
        qDebug() << "Something went wrong: " << e.what();
        std::cout << "Something went wrong: " << e.what();
        fmt::print("Something went wrong: {}", e.what());

        // ... or similar

    }
    ```

3) Due to low operator priority it is important to always enclose the left and right
   contraction performed by operator<<() and operator>>() in braces!

4) To compile use

    ```bash
    cd build
    cmake ..
    cmake --build .
    ```

## GA Product Expression Generator (ga_prdxpr)

The `ga_prdxpr` tool generates optimized C++ expressions for geometric algebra operations
across multiple algebras. It supports selective output generation through command-line
arguments.

### Quick Reference - Common Workflows

```bash
# View only product tables for one algebra (no coefficient expressions)
cd build/ga_prdxpr && ./ga_prdxpr --algebra=ega2d --output=tables && cd ../..

# View only coefficient expressions for one product (no tables)
cd build/ga_prdxpr && ./ga_prdxpr --algebra=ega3d --products=gpr --output=coeffs && cd ../..

# Compare one product across multiple algebras
cd build/ga_prdxpr && ./ga_prdxpr --algebra=ega2d,ega3d --products=gpr --output=coeffs && cd ../..
```

### All Command-Line Options

```bash
cd build/ga_prdxpr && ./ga_prdxpr --help && cd ../..         # Show help with examples
cd build/ga_prdxpr && ./ga_prdxpr --list && cd ../..         # List available algebras and products
cd build/ga_prdxpr && ./ga_prdxpr && cd ../..                # Generate everything (searchable reference)
```

**Filtering Options:**

- `--algebra=NAMES` or `--a=NAMES`: Generate only specified algebras (comma-separated)
  - Example: `--algebra=ega2d,ega3d`
  - Available: `ega2d`, `ega3d`, `pga2dp`, `pga3dp`, `sta4d`
  - Special keywords: `all` (all algebras), `none` (no algebras)

- `--products=NAMES` or `--p=NAMES`: Generate only specified products (comma-separated)
  - Example: `--products=gpr,wdg,dot`
  - Common products: `gpr`, `wdg`, `dot`, `cmt`, `l_contract`, `r_contract`, `sandwich_gpr`
  - Special keywords: `all` (all products), `none` (no products)

- `--output=TYPES` or `--o=TYPES`: Show only specified output types (comma-separated)
  - Example: `--output=coeffs,tables`
  - Types: `coeffs` (coefficient expressions), `tables` (basis multiplication tables),
    `metrics` (C++ array exports)
  - Special keywords: `all` (all types), `none` (no types)

- `--show-sym-asym=BOOL`: Show symmetric/asymmetric parts of basis tables
  - Values: `true`, `false`, `yes`, `no`, `1`, `0`
  - Default: `false`
  - Example: `--show-sym-asym=true`

### Use Case Examples

```bash
# Generate everything for EGA2D only
cd build/ga_prdxpr && ./ga_prdxpr --algebra=ega2d && cd ../..

# Show only basis tables for all algebras
cd build/ga_prdxpr && ./ga_prdxpr --output=tables && cd ../..

# Show tables with symmetric/asymmetric decomposition
cd build/ga_prdxpr && ./ga_prdxpr --output=tables --show-sym-asym=true && cd ../..

# Compare geometric product across all algebras (expressions only)
cd build/ga_prdxpr && ./ga_prdxpr --products=gpr --output=coeffs && cd ../..

# Generate metric exports for PGA algebras only
cd build/ga_prdxpr && ./ga_prdxpr --algebra=pga2dp,pga3dp --output=metrics && cd ../..

# View wedge and dot products for EGA3D (tables and expressions)
cd build/ga_prdxpr && ./ga_prdxpr --algebra=ega3d --products=wdg,dot && cd ../..
```

### C++ Code Generation (`--output=code`)

Beyond coefficient expressions and basis tables, ga_prdxpr can emit
ready-to-paste C++ function implementations for diverse products.
Opt-in only — never appears in the default output:

```bash
ga_prdxpr --algebra=ega2d --output=code                  # one algebra
ga_prdxpr --algebra=pga3dp --products=gpr --output=code  # one product
ga_prdxpr --output=all                                   # everything
```

### Available Algebras

- **ega2d**: Euclidean Geometric Algebra 2D
- **ega3d**: Euclidean Geometric Algebra 3D
- **pga2dp**: Projective Geometric Algebra 2D
- **pga3dp**: Projective Geometric Algebra 3D
<!-- - **sta4d**: Space-Time Algebra 4D -->

### Common Product Types

- Geometric product (available as `operator*()` in the library)
- `wdg` - Wedge (outer) product
- `dot` - Inner (dot) product
- `cmt` - Commutator product
- Left (bulk) contraction (implemented as `operator<<()` in the library)
- Right (bulk) contraction (implemented as `operator>>()` in the library)
- Left expansion (EGA `l_expand`; PGA uses `l_bulk_expand2dp()`/`l_weight_expand2dp()`)
- Right expansion (EGA `r_expand`; PGA uses `r_bulk_expand2dp()`/`r_weight_expand2dp()`)
- `rgpr` - Regressive geometric product
- `rwdg` - Regressive wedge product
- `rdot` - Regressive inner (dot) product
- `rcmt` - Regressive commutator product
- Geometric sandwich product R ⟑ X ⟑ rev(R) → implemented as `rotate()` in EGA
- Regressive sandwich product M ⟇ X ⟇ rrev(M) → implemented as `move2dp()`/`move3dp()` in PGA

**Note**: PGA algebras have additional products for bulk/weight operations
(`l_bulk_contract2dp()`, `r_weight_expand3dp()`, etc.). Use `--list` or
`--algebra=pga2dp` / `--algebra=pga3dp` to see all available products.

### Default Behavior

With no arguments, `ga_prdxpr` generates complete output for all algebras, products, and
output types. This serves as a comprehensive searchable reference.
