# Vector2D with nested_formatter Example

This example demonstrates how to implement a custom formatter using `fmt::nested_formatter` for a simple 2D vector type.

## Vector2D Type Definition

```cpp
#include <fmt/format.h>
#include <cmath>

// Simple 2D vector template with public components
template<typename T>
struct Vector2D {
    T x;
    T y;
    
    // Constructors
    constexpr Vector2D() : x(T{}), y(T{}) {}
    constexpr Vector2D(T x_val, T y_val) : x(x_val), y(y_val) {}
    
    // Basic operations
    constexpr Vector2D operator+(const Vector2D& other) const {
        return Vector2D(x + other.x, y + other.y);
    }
    
    T magnitude() const {
        return std::sqrt(x * x + y * y);
    }
};

// Type aliases
using Vector2f = Vector2D<float>;
using Vector2d = Vector2D<double>;
using Vector2i = Vector2D<int>;
```

## Formatter Implementation with nested_formatter

```cpp
template<typename T>
struct fmt::formatter<Vector2D<T>> : fmt::nested_formatter<T> {
    // Parse method - delegates to nested_formatter for component formatting
    constexpr auto parse(format_parse_context& ctx) {
        return fmt::nested_formatter<T>::parse(ctx);
    }
    
    // Format method - handles vector structure, delegates component formatting
    template<typename FormatContext>
    auto format(const Vector2D<T>& vec, FormatContext& ctx) const {
        // Format as (x, y)
        auto out = fmt::format_to(ctx.out(), "(");
        
        // Use nested formatter for x component
        out = fmt::nested_formatter<T>::format(vec.x, ctx);
        
        out = fmt::format_to(ctx.out(), ", ");
        
        // Update context and use nested formatter for y component
        ctx.advance_to(out);
        out = fmt::nested_formatter<T>::format(vec.y, ctx);
        
        return fmt::format_to(out, ")");
    }
};
```

## Key Points

1. **Inheritance**: The formatter inherits from `fmt::nested_formatter<T>` where `T` is the component type
2. **Parse delegation**: `parse()` method delegates to the nested formatter to handle format specifications
3. **Format delegation**: Component formatting is handled by calling `fmt::nested_formatter<T>::format()`
4. **Context management**: Use `ctx.advance_to(out)` when formatting multiple components
5. **Automatic support**: All format specifications (`:08`, `.3f`, `:x`, etc.) work automatically

## Usage Examples

```cpp
int main() {
    // Create vectors
    Vector2f float_vec(3.14159f, 2.71828f);
    Vector2d double_vec(1.41421, 1.73205);
    Vector2i int_vec(42, -17);
    
    // Basic output using fmt::print
    fmt::print("Float vector: {}\n", float_vec);           // (3.14159, 2.71828)
    fmt::print("Double vector: {}\n", double_vec);         // (1.41421, 1.73205)
    fmt::print("Integer vector: {}\n", int_vec);           // (42, -17)
    
    // Format specifications work automatically
    fmt::print("2 decimals: {:.2f}\n", float_vec);         // (3.14, 2.72)
    fmt::print("Zero-padded: {:08}\n", int_vec);           // (00000042, -0000017)
    fmt::print("Hex: {:x}\n", int_vec);                    // (2a, -11)
    fmt::print("Scientific: {:e}\n", double_vec);          // (1.414210e+00, 1.732050e+00)
    
    // Creating formatted strings using fmt::format
    std::string vec_str = fmt::format("{:.3f}", float_vec);
    fmt::print("Formatted string: '{}'\n", vec_str);       // '(3.142, 2.718)'
    
    // Use in logging or configuration
    std::string log_entry = fmt::format("Position: {:.2f}", float_vec);
    fmt::print("Log: {}\n", log_entry);                    // Log: Position: (3.14, 2.72)
    
    return 0;
}
```

## Benefits of nested_formatter

- **No manual format parsing**: The nested formatter handles all format specifications
- **Type consistency**: Component formatting matches standard library behavior
- **Automatic precision**: Format specs like `.3f` apply to both components
- **Width and alignment**: Format specs like `:08` work correctly
- **Base conversions**: Format specs like `:x` for hex work automatically
- **Extensibility**: Easy to modify structure without touching component formatting logic

This approach separates concerns: your formatter handles the vector's structure (parentheses, comma), while nested_formatter handles the component formatting details.