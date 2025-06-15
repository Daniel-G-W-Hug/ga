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

## File Structure

```
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
