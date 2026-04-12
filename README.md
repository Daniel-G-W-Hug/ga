# Overview on modules

- ga: header-only library for geometric algebra in 2d and 3d
- ga_test: tests for ga header only library
- ga_lua: lua scripting based on the ga library
- ga_view: qt-based visualization of 2d ga objects and operations
- ga_prdxpr: generate product expressions to be used in implementation

## Dependencies

- ga: fmt (must be included by every user - header only lib!)
- ga_test: fmt, doctest
- ga_lua: fmt, lua, sol2
- ga_view: fmt, qt6
- ga_prdxpr: fmt

## Installing dependencies on macOS

```bash
    brew install fmt         // always required
    brew istall doctest      // only required for ga_test
    brew install qt6         // only required for ga_view
    brew install lua         // lua and sol2 only required for ga_lua
    sol2 to be installed from: <https://github.com/ThePhD/sol2>
```

## Installing dependencies on Windows

```bash
    fmt to be installed from: https://github.com/fmtlib/fmt
    doctest to be installed from: https://github.com/doctest/doctest
    winget install --id DEVCOM.Lua
    (you will need a x64 lua binary or linking will not work)
    qt6 to be installed from official installer
    sol2 to be installed from: <https://github.com/ThePhD/sol2>
```

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
  - Common products: `gpr`, `wdg`, `dot`, `cmt`, `lcontract`, `rcontract`, `sandwich_gpr`
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

### Available Algebras

- **ega2d**: Euclidean Geometric Algebra 2D
- **ega3d**: Euclidean Geometric Algebra 3D
- **pga2dp**: Projective Geometric Algebra 2D
- **pga3dp**: Projective Geometric Algebra 3D
<!-- - **sta4d**: Space-Time Algebra 4D -->

### Common Product Types

- `gpr` - Geometric product
- `wdg` - Wedge (outer) product
- `dot` - Inner (dot) product
- `cmt` - Commutator product
- `lcontract` - Left contraction
- `rcontract` - Right contraction
- `lexpand` - Left expansion
- `rexpand` - Right expansion
- `rgpr` - Regressive geometric product
- `rwdg` - Regressive wedge product
- `rdot` - Regressive inner (dot) product
- `rcmt` - Regressive commutator product
- `sandwich_gpr` - Sandwich product (geometric)
- `sandwich_rgpr` - Sandwich product (regressive)

**Note**: PGA algebras have additional products for bulk/weight operations
(`left_bulk_contract`, `right_weight_expand`, etc.). Use `--algebra=pga2dp` or
`--algebra=pga3dp` to see all products for those algebras.

### Default Behavior

With no arguments, `ga_prdxpr` generates complete output for all algebras, products, and
output types. This serves as a comprehensive searchable reference.
