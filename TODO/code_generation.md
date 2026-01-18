# Code Generation Extension for ga_prdxpr

## Overview

This document outlines a plan to extend ga_prdxpr with a `--output=code` option that
generates copy-paste-ready C++ function implementations, eliminating manual transcription
errors.

## Architectural Principle: Single Source of Truth

**Critical Design Constraint**: The existing `ProductConfig` with its `.cases` array
(defined in files like `ga_prdxpr_pga2dp_config.cpp`) is the **single specification
point** for all product definitions. The code generation feature must:

1. **Consume existing configuration** - Use `ProductConfig` and `AlgebraData` as input
2. **Be a pure output transformation** - No parallel specification system
3. **Integrate with existing generator flow** - Same data path as `--output=coeffs`

The user defines products once in the config files, then selects output format via command
line:

- `--output=coeffs` → Human-readable coefficient expressions (current)
- `--output=code` → Copy-paste-ready C++ functions (new)
- `--output=tables` → Basis multiplication tables (current)

## Problem Statement

### Current Workflow Issues

1. **Manual Translation Required**: The current `--output=coeffs` produces human-readable
   coefficient expressions that must be manually translated into C++ code.

2. **Non-Uniform Implementation Patterns**: The existing implementations use different
   patterns:

   | Input Types      | Return Type | Constructor Pattern                                    |
   |------------------|-------------|--------------------------------------------------------|
   | `MVec, MVec`     | `MVec`      | `MVec<ctype>(c0, c1, c2, c3, c4, c5, c6, c7)`          |
   | `MVec_E, MVec_E` | `MVec_E`    | `MVec_E<ctype>(Scalar<ctype>(...), BiVec<ctype>(...))` |
   | `MVec_U, MVec_U` | `MVec_E`    | `MVec_E<ctype>(Scalar<ctype>(...), BiVec<ctype>(...))` |
   | `BiVec, BiVec`   | `MVec_E`    | Delegates to `dot()` and `cmt()`                       |
   | `Vec, Vec`       | `MVec_E`    | Delegates to `dot()` and `wdg()`                       |
   | `BiVec, Vec`     | `MVec_U`    | Delegates to `(v << B)` and `wdg()`                    |
   | `PScalar, BiVec` | `Vec/BiVec` | Scalar multiplication delegation                       |

3. **Variable Naming Inconsistency**: Reference output uses:
   - Full multivectors: `A.c0`, `B.c0`
   - Bivectors: `B.x`, `B.y`, `B.z` or `B1.x`, `B2.x`
   - Vectors: `v.x`, `v.y`, `v.z` or `v1.x`, `v2.x`
   - Scalars: `s`, `s1`, `s2`
   - Pseudoscalars: `ps`, `ps1`, `ps2`

4. **Coefficient-to-Position Mapping**: Output shows basis labels (`| e1`, `| e31`) but
   implementation uses positional indices (`c0`, `c1`, etc.) or named components (`x`,
   `y`, `z`).

## Proposed Solution

### New Output Type: `--output=code`

Generate complete, compilable C++ function bodies that match the existing code style.

### Architecture: Integration with Existing Generator

The key insight is that the existing generator already computes coefficient expressions
for each product case. The code generation module hooks into this existing flow:

```text
┌─────────────────────────────────────────────────────────────────────────────┐
│                         EXISTING FLOW (unchanged)                           │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ProductConfig          AlgebraData           Generator                     │
│  (from *_config.cpp)    (coefficients,        (computes expressions)        │
│                          filters, basis)                                    │
│         │                     │                      │                      │
│         └─────────────────────┴──────────────────────┘                      │
│                               │                                             │
│                               ▼                                             │
│                    ┌──────────────────────┐                                 │
│                    │  Coefficient Results │                                 │
│                    │  (per basis element) │                                 │
│                    └──────────────────────┘                                 │
│                               │                                             │
└───────────────────────────────┼─────────────────────────────────────────────┘
                                │
                ┌───────────────┼───────────────┐
                │               │               │
                ▼               ▼               ▼
         ┌──────────┐    ┌──────────┐    ┌──────────┐
         │ coeffs   │    │  tables  │    │   code   │  ← NEW
         │ output   │    │  output  │    │  output  │
         └──────────┘    └──────────┘    └──────────┘
```

### File Structure

```text
ga_prdxpr/
├── src_prdxpr/
│   ├── ga_prdxpr_options.hpp/cpp      # Add OutputType::Code
│   ├── ga_prdxpr_generator.hpp/cpp    # Add emit_code() alongside emit_coeffs()
│   ├── ga_prdxpr_main.cpp             # Wire up should_show_code()
│   │
│   ├── codegen/                       # NEW: Code generation module
│   │   ├── ga_codegen.hpp             # Public interface (single include)
│   │   ├── ga_codegen_types.hpp       # Type metadata structures
│   │   ├── ga_codegen_types.cpp       # Type registry implementation
│   │   ├── ga_codegen_emitter.hpp     # C++ code emitter
│   │   └── ga_codegen_emitter.cpp     # Emitter implementation
│   │
│   ├── ga_prdxpr_pga2dp_config.cpp    # UNCHANGED - still the source of truth
│   ├── ga_prdxpr_ega3d_config.cpp     # UNCHANGED
│   └── [other algebra configs...]     # UNCHANGED
```

## Detailed Design

### Integration Point: Leveraging Existing ProductCase Data

The existing `ProductCase` structure (from `ga_prdxpr_config.hpp`) already contains:

```cpp
struct ProductCase {
    std::string case_name;       // e.g., "mv_e * mv_e -> mv_e"
    std::string left_coeff;      // e.g., "A_even"
    std::string right_coeff;     // e.g., "B_even"
    std::string left_filter;     // e.g., "mv_e"
    std::string right_filter;    // e.g., "mv_e"
    // ... sandwich product fields
};
```

The `case_name` field already encodes:

- Left operand type (before operator)
- Operator symbol (`*`, `^`, or function name like `dot`)
- Right operand type (after operator)
- Result type (after `->`)

**Key Insight**: We can parse `case_name` to extract type information, then use the
already-computed coefficient expressions to generate the C++ code.

### Phase 1: Type Metadata System

Create a type registry that maps filter names to C++ type information.

```cpp
// ga_codegen_types.hpp

namespace codegen {

// Describes how a type stores its coefficients
enum class StorageKind {
    Indexed,      // c0, c1, c2, ... (MVec, MVec_E, MVec_U)
    Named,        // x, y, z (Vec, BiVec, TriVec)
    SingleValue   // Scalar, PScalar (implicit conversion to value_t)
};

// Maps filter names (from ProductConfig) to C++ type names
struct TypeInfo {
    std::string filter_name;         // "mv_e", "vec", "bivec", "s", "ps"
    std::string cpp_type_name;       // "MVec2dp_E", "Vec2dp", "BiVec2dp", etc.
    StorageKind storage;
    bool pass_by_value;              // true for Scalar, PScalar

    // For composite return types (MVec_E = Scalar + BiVec)
    bool uses_nested_constructors;
    std::vector<std::string> nested_types;  // {"Scalar2dp", "BiVec2dp"}

    // Coefficient component mapping
    // Maps basis element index to component accessor
    // e.g., for Vec2dp: {0 -> "x", 1 -> "y", 2 -> "z"}
    // e.g., for MVec2dp_E: {0 -> "c0", 1 -> "c1", 2 -> "c2", 3 -> "c3"}
    std::map<int, std::string> coefficient_map;
};

// Per-algebra type registry
class TypeRegistry {
public:
    explicit TypeRegistry(std::string const& algebra_name);

    // Lookup type info by filter name
    TypeInfo const& get_type(std::string const& filter_name) const;

    // Get the result type for a given case_name
    TypeInfo const& get_result_type(std::string const& case_name) const;

    // Map basis element to coefficient access expression
    // e.g., ("vec", "e1", "v") -> "v.x"
    // e.g., ("mv_e", "e31", "A") -> "A.c1"
    std::string coefficient_access(std::string const& filter_name,
                                   std::string const& basis_element,
                                   std::string const& var_name) const;
};

// Factory function - creates registry for a specific algebra
TypeRegistry create_type_registry(std::string const& algebra_name,
                                  AlgebraData const& data);

} // namespace codegen
```

### Phase 2: Case Name Parser and Signature Builder

Parse the existing `case_name` strings to build function signatures.

```cpp
// ga_codegen_emitter.hpp

namespace codegen {

// Parsed information from case_name like "mv_e * mv_e -> mv_e"
struct ParsedCase {
    std::string left_filter;      // "mv_e"
    std::string right_filter;     // "mv_e"
    std::string result_filter;    // "mv_e"
    std::string op_symbol;        // "*", "^", or function name like "dot"
    bool is_infix;                // true for *, ^; false for dot(), wdg(), etc.
};

// Parse case_name from ProductCase
ParsedCase parse_case_name(std::string const& case_name);

// Determine parameter names based on filter types
// Uses existing conventions: A/B for mv, v/v1/v2 for vec, B/B1/B2 for bivec, etc.
struct ParamNames {
    std::string left_param;       // "A", "v", "B", "s", "ps", "M"
    std::string right_param;      // "B", "v", "B2", "s", "ps", "M"
};

ParamNames determine_param_names(std::string const& left_filter,
                                  std::string const& right_filter,
                                  std::string const& left_coeff,
                                  std::string const& right_coeff);

// Describes a complete function signature
struct FunctionSignature {
    std::string function_name;       // "operator*", "wdg", "rgpr", etc.
    std::string lhs_cpp_type;        // "MVec2dp_E<T>", "BiVec2dp<T>", etc.
    std::string rhs_cpp_type;
    std::string return_cpp_type;

    ParamNames params;

    // Whether parameters are const ref or by value
    bool lhs_by_value;               // true for Scalar, PScalar
    bool rhs_by_value;
};

// Build signature from ProductConfig case and TypeRegistry
FunctionSignature build_signature(ProductCase const& case_def,
                                  ProductConfig const& product,
                                  TypeRegistry const& types);

} // namespace codegen
```

### Phase 3: Code Emitter

Generate properly formatted C++ code using the coefficient expressions already computed by
the existing generator.

```cpp
// ga_codegen_emitter.hpp (continued)

namespace codegen {

// Coefficient data from existing generator output
// This is what the current emit_coeffs() produces internally
struct CoefficientResult {
    std::string basis_element;    // "1", "e1", "e31", "e321", etc.
    std::string expression;       // "A.c0 * B.c0 - A.c3 * B.c3"
    bool is_zero;                 // true if expression simplifies to 0
};

class CodeEmitter {
public:
    struct Options {
        bool include_comments = true;
        bool include_template_header = true;
        int indent_width = 4;
        int max_line_length = 88;  // Match project style
    };

    explicit CodeEmitter(Options const& opts = {});

    // Generate a complete function from:
    // - ProductCase (from config)
    // - ProductConfig (for product name/operator)
    // - TypeRegistry (for type mappings)
    // - Coefficient results (from existing generator computation)
    std::string emit_function(
        ProductCase const& case_def,
        ProductConfig const& product,
        TypeRegistry const& types,
        std::vector<CoefficientResult> const& coefficients) const;

private:
    Options opts_;

    // Helper methods
    std::string emit_template_header() const;
    std::string emit_signature(FunctionSignature const& sig) const;
    std::string emit_body(FunctionSignature const& sig,
                          TypeRegistry const& types,
                          std::vector<CoefficientResult> const& coefficients) const;

    // Format long expressions with line continuation
    std::string format_expression(std::string const& expr,
                                  std::string const& indent) const;

    // Determine if result can use nested constructors
    // e.g., MVec_E(Scalar(...), BiVec(...)) vs MVec_E(c0, c1, c2, c3)
    bool should_use_nested_constructors(TypeInfo const& result_type,
                                        std::vector<CoefficientResult> const& coeffs) const;
};

} // namespace codegen
```

### Phase 4: Integration with Existing Generator

The key integration point is in the existing generator loop. Currently, for each product
case, the generator computes coefficient expressions and calls `emit_coeffs()`. We add a
parallel path for `emit_code()` that uses the same computed coefficients.

```cpp
// In ga_prdxpr_generator.cpp - modified generation loop

void generate_product_cases(AlgebraData const& algebra,
                            ProductConfig const& product,
                            GeneratorOptions const& options)
{
    // Create type registry once per algebra (lazy initialization)
    static std::map<std::string, codegen::TypeRegistry> registries;
    if (registries.find(algebra.name) == registries.end()) {
        registries[algebra.name] = codegen::create_type_registry(algebra.name, algebra);
    }
    auto const& types = registries[algebra.name];

    codegen::CodeEmitter emitter;

    for (auto const& case_def : product.cases) {
        // EXISTING: Compute coefficient expressions (unchanged)
        auto coefficients = compute_coefficients(algebra, product, case_def);

        // EXISTING: Output coeffs format
        if (options.should_show_coeffs()) {
            emit_coeffs(algebra, product, case_def, coefficients);
        }

        // NEW: Output code format
        if (options.should_show_code()) {
            std::string code = emitter.emit_function(case_def, product, types, coefficients);
            fmt::println("{}", code);
        }
    }
}
```

**Critical**: The `compute_coefficients()` function already exists and produces the
coefficient expressions. The code emitter simply transforms those expressions into C++
syntax.

### Phase 5: Output Options Integration

Add the new output type alongside existing options.

```cpp
// In ga_prdxpr_options.hpp

enum class OutputType {
    Coeffs,   // Human-readable coefficient expressions
    Tables,   // Basis multiplication tables
    Metrics,  // Extended metric exports
    Code      // NEW: C++ function implementations
};

// In GeneratorOptions
bool should_show_code() const
{
    return output_types.empty() || output_types.count(OutputType::Code) > 0;
}
```

```cpp
// In ga_prdxpr_options.cpp - parse_output_list()

else if (type == "code") {
    options_.output_types.insert(OutputType::Code);
}
```

## Coefficient-to-Code Transformation Details

### Current Coefficient Output Format

The existing `--output=coeffs` produces output like this:

```text
pga2dp: mv_e * mv_e -> mv_e

        | 1           A.c0*B.c0 - A.c3*B.c3
        | e31         A.c0*B.c1 + A.c1*B.c0 - A.c2*B.c3 + A.c3*B.c2
        | e32         A.c0*B.c2 + A.c1*B.c3 + A.c2*B.c0 - A.c3*B.c1
        | e12         A.c0*B.c3 + A.c3*B.c0
```

### Transformation Process

The code emitter performs these transformations:

1. **Parse case_name**: `"mv_e * mv_e -> mv_e"` → types and operator
2. **Map filter to C++ type**: `"mv_e"` → `MVec2dp_E<T>`
3. **Map basis to position**: `"e31"` → position 1 in BiVec2dp (`.x` component)
4. **Group by grade**: scalar coefficients, bivector coefficients
5. **Build constructor**: `MVec2dp_E(Scalar2dp(...), BiVec2dp(...))`

### Basis-to-Component Mapping (PGA2DP Example)

| Basis Element | Grade | Type       | Component |
|---------------|-------|------------|-----------|
| `1`           | 0     | Scalar2dp  | (value)   |
| `e1`          | 1     | Vec2dp     | `.x`      |
| `e2`          | 1     | Vec2dp     | `.y`      |
| `e3`          | 1     | Vec2dp     | `.z`      |
| `e31`         | 2     | BiVec2dp   | `.x`      |
| `e32`         | 2     | BiVec2dp   | `.y`      |
| `e12`         | 2     | BiVec2dp   | `.z`      |
| `e321`        | 3     | PScalar2dp | (value)   |

### Coefficient Access Mapping

For input operands, the TypeRegistry provides mappings:

| Filter  | Variable | Coefficient | Access Expression                  |
|---------|----------|-------------|------------------------------------|
| `mv_e`  | `A`      | `c0`        | `A.c0` (scalar)                    |
| `mv_e`  | `A`      | `c1-c3`     | `A.c1`, `A.c2`, `A.c3` (bivector)  |
| `vec`   | `v`      | `x,y,z`     | `v.x`, `v.y`, `v.z`                |
| `bivec` | `B`      | `x,y,z`     | `B.x`, `B.y`, `B.z`                |
| `s`     | `s`      | (value)     | `value_t{s}`                       |
| `ps`    | `ps`     | (value)     | `value_t{ps}`                      |

## Implementation Plan

### Step 1: Add OutputType::Code

**Files to modify:**

- `ga_prdxpr_options.hpp`: Add `Code` to `OutputType` enum
- `ga_prdxpr_options.cpp`: Add parsing for `--output=code`, update
  `is_valid_output_type()`
- `ga_prdxpr_options.hpp`: Add `should_show_code()` method to `GeneratorOptions`

**Changes:**

1. Add `Code` to `OutputType` enum in `ga_prdxpr_options.hpp`
2. Add `should_show_code()` method to `GeneratorOptions`
3. Update argument parser to accept `--output=code`
4. Update help text and documentation

### Step 2: Create Type Metadata System

**New files to create:**

- `codegen/ga_codegen_types.hpp`: Type metadata structures
- `codegen/ga_codegen_types.cpp`: TypeRegistry implementation

**Implementation tasks:**

1. Define `StorageKind` enum (Indexed, Named, SingleValue)
2. Define `TypeInfo` struct with coefficient mapping
3. Implement `TypeRegistry` class with lookup methods
4. Create factory function `create_type_registry()` for each algebra

**Type metadata for each algebra:**

- EGA2D: Scalar2d, Vec2d, BiVec2d, PScalar2d, MVec2d, MVec2d_E
- EGA3D: Add TriVec3d, MVec3d_U
- PGA2DP: Vec2dp, BiVec2dp, etc. with degenerate metric awareness
- PGA3DP: Full 4D projective types
- STA4D: Spacetime algebra types

### Step 3: Implement Case Parser and Signature Builder

**New files to create:**

- `codegen/ga_codegen_emitter.hpp`: ParsedCase, ParamNames, FunctionSignature structs

**Implementation tasks:**

1. Implement `parse_case_name()` to extract types and operator from case_name string
2. Implement `determine_param_names()` using existing conventions (A/B for mv, v/v1/v2 for
   vec)
3. Implement `build_signature()` to create complete function signature from ProductCase

**Case name parsing examples:**

- `"mv_e * mv_e -> mv_e"` → left=mv_e, op=*, right=mv_e, result=mv_e
- `"vec ^ bivec -> ps"` → left=vec, op=^, right=bivec, result=ps
- `"dot(bivec, bivec) -> s"` → left=bivec, op=dot, right=bivec, result=s

### Step 4: Implement Code Emitter

**New files to create:**

- `codegen/ga_codegen_emitter.cpp`: CodeEmitter implementation
- `codegen/ga_codegen.hpp`: Public interface (single include)

**Implementation tasks:**

1. Implement `emit_function()` main entry point
2. Implement helper methods:
   - `emit_template_header()`: Template declaration with requires clause
   - `emit_signature()`: Function signature with east-const style
   - `emit_body()`: Return statement with constructors
   - `format_expression()`: Line wrapping for long expressions
3. Handle all constructor patterns:
   - Direct indexed: `MVec(c0, c1, c2, ...)`
   - Nested typed: `MVec_E(Scalar(...), BiVec(...))`
   - Delegation: `MVec_E(dot(a,b), wdg(a,b))`

### Step 5: Generator Integration

**Files to modify:**

- `ga_prdxpr_generator.hpp`: Add code generation interface
- `ga_prdxpr_generator.cpp`: Integrate CodeEmitter in generation loop
- `CMakeLists.txt`: Add codegen/ source files

**Implementation tasks:**

1. Add TypeRegistry caching per algebra in generator
2. Create CodeEmitter instance in generation loop
3. Call `emitter.emit_function()` when `should_show_code()` is true
4. Ensure coefficient computation is shared between coeffs and code output

### Step 6: Testing and Validation

**Testing approach:**

1. Generate code for all products in one algebra (e.g., pga2dp)
2. Compare generated code with existing implementations character-by-character
3. Verify generated code compiles by creating a test harness
4. Run existing test suite with generated code substituted

**Regression tests:**

- Add test cases for each algebra
- Verify output stability across generator changes

### Step 7: Documentation

**Files to update:**

- `README.md`: Add `--output=code` to usage section
- `ga_prdxpr_options.cpp`: Update `print_help()` and `print_list()`

**Documentation tasks:**

1. Add usage examples for `--output=code`
2. Document the codegen module architecture
3. Add workflow examples combining filters with code output

## Handling Non-Uniformity

The key insight is that non-uniformity exists for good reasons:

1. **Performance**: Delegation avoids code duplication and enables optimizations
2. **Readability**: `dot(a, b) + wdg(a, b)` is clearer than expanded coefficients
3. **Maintainability**: Changes to `dot()` automatically propagate

### Strategy: Preserve Intentional Delegation

The code generator should:

1. **Detect delegation opportunities** based on the product definition
2. **Generate delegation code** when mathematically equivalent
3. **Fall back to direct coefficients** only when necessary

Example detection rules:

```cpp
// In pattern builder
CodePattern determine_pattern(ProductDefinition const& product,
                              std::string const& lhs_type,
                              std::string const& rhs_type)
{
    // Rule: vec * vec = dot + wdg (even grade result)
    if (is_vector(lhs_type) && is_vector(rhs_type) &&
        product.type == ProductType::Geometric) {
        return CodePattern::DelegateToProducts;
    }

    // Rule: bivec * bivec = -dot + cmt (in 3D)
    if (is_bivector(lhs_type) && is_bivector(rhs_type) &&
        product.type == ProductType::Geometric) {
        return CodePattern::DelegateToProducts;
    }

    // Rule: ps * ps = 0 in degenerate algebras
    if (is_pseudoscalar(lhs_type) && is_pseudoscalar(rhs_type) &&
        algebra.is_degenerate) {
        return CodePattern::ZeroResult;
    }

    // Default: direct coefficients
    return CodePattern::DirectCoefficients;
}
```

## Example Output

For `--algebra=pga2dp --products=gpr --output=code`:

```cpp
// pga2dp geometric product mv_e * mv_e -> mv_e
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(MVec2dp_E<T> const& A,
                                                        MVec2dp_E<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(
        Scalar2dp<ctype>(A.c0 * B.c0 - A.c3 * B.c3),
        BiVec2dp<ctype>(A.c0 * B.c1 + A.c1 * B.c0 - A.c2 * B.c3 + A.c3 * B.c2,
                        A.c0 * B.c2 + A.c1 * B.c3 + A.c2 * B.c0 - A.c3 * B.c1,
                        A.c0 * B.c3 + A.c3 * B.c0));
}

// pga2dp geometric product vec * vec -> mv_e
// Uses delegation pattern: gpr(v1, v2) = dot(v1, v2) + wdg(v1, v2)
template <typename T, typename U>
    requires(numeric_type<T> && numeric_type<U>)
constexpr MVec2dp_E<std::common_type_t<T, U>> operator*(Vec2dp<T> const& a,
                                                        Vec2dp<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2dp_E<ctype>(dot(a, b), wdg(a, b));
}
```

## Future Extensions

### Validation Mode

Add `--output=validate --file=path` to compare generated code against existing
implementations:

```bash
# Validate pga2dp products against existing implementation
cd build/ga_prdxpr && ./ga_prdxpr --algebra=pga2dp --output=validate \
    --file=../ga/ga_pga2dp_ops_products.hpp && cd ../..
```

Output would show:

- Matching functions (green checkmarks)
- Mismatched functions with diff highlighting
- Missing functions in implementation
- Extra functions not in generator

### Diff Mode

Show side-by-side comparison between generated and existing code:

```bash
cd build/ga_prdxpr && ./ga_prdxpr --algebra=pga2dp --products=gpr \
    --output=diff --file=../ga/ga_pga2dp_ops_products.hpp && cd ../..
```

### Update Mode

Automatically update existing files with generated code (with backup):

```bash
cd build/ga_prdxpr && ./ga_prdxpr --algebra=pga2dp --output=update \
    --file=../ga/ga_pga2dp_ops_products.hpp && cd ../..
```

This would:

1. Create backup of original file
2. Parse existing file to find function boundaries
3. Replace matching functions with generated versions
4. Preserve non-product code (includes, namespace, etc.)

### Test Generation

Generate unit tests for each product function:

```bash
cd build/ga_prdxpr && ./ga_prdxpr --algebra=pga2dp --output=tests && cd ../..
```

Would generate doctest test cases verifying:

- Known mathematical identities (e.g., `gpr(v, rev(v)) = dot(v, v)`)
- Zero results for degenerate cases
- Associativity and distributivity properties

## Risk Assessment

| Risk                                   | Mitigation                                         |
|----------------------------------------|----------------------------------------------------|
| Generated code doesn't compile         | Extensive testing against existing implementations |
| Performance regression from naive code | Preserve delegation patterns                       |
| Maintenance burden of codegen module   | Clean modular design with clear interfaces         |
| Type system complexity                 | Comprehensive TypeMetadata with validation         |

## Success Criteria

1. Generated code compiles without modification
2. Generated code matches existing implementations exactly (where not delegating)
3. Delegation patterns are correctly detected and preserved
4. All five algebras (EGA2D, EGA3D, PGA2DP, PGA3DP, STA4D) supported
5. Clear module boundaries that don't increase coupling with existing code

## Summary

This extension adds a `--output=code` option to ga_prdxpr that generates copy-paste-ready
C++ function implementations. The key design principles are:

1. **Single Source of Truth**: ProductConfig remains the only place to define products.
   Code generation is a pure output transformation that consumes existing configuration.

2. **Integration with Existing Flow**: The code emitter hooks into the existing generator
   loop, using already-computed coefficient expressions. No parallel computation path.

3. **Modular Architecture**: The codegen/ subdirectory encapsulates all code generation
   logic with clear interfaces, minimizing coupling with existing generator code.

4. **Preserve Intentional Patterns**: The system detects and preserves delegation patterns
   (e.g., `dot + wdg` for `vec * vec`) that exist in the current implementations for
   performance and maintainability reasons.

**Usage Example:**

```bash
# Generate C++ code for all pga2dp products
cd build/ga_prdxpr && ./ga_prdxpr --algebra=pga2dp --output=code && cd ../..

# Generate code for specific product only
cd build/ga_prdxpr && ./ga_prdxpr --algebra=pga2dp --products=gpr --output=code && cd ../..

# Combine with other output types
cd build/ga_prdxpr && ./ga_prdxpr --algebra=pga2dp --output=coeffs,code && cd ../..
```
