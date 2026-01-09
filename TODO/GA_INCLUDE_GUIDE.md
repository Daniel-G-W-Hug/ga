# GA Library Include Guide

This guide helps you choose the optimal include strategy for your geometric algebra applications.

## Quick Reference

| Use Case | Include Header | Compilation Speed | Features |
|----------|----------------|-------------------|----------|
| **Type definitions only** | `#include "ga/ga_minimal.hpp"` | 🚀 Fastest (60-70% faster) | All types, no operations |
| **2D Euclidean GA only** | `#include "ga/ga_ega2d_only.hpp"` | 🏃 Fast (30-40% faster) | Complete 2D EGA |
| **3D Euclidean GA only** | `#include "ga/ga_ega3d_only.hpp"` | 🏃 Fast (25-35% faster) | Complete 3D EGA |
| **2DP Projective GA only** | `#include "ga/ga_pga2dp_only.hpp"` | 🏃 Fast (35-45% faster) | Complete 2DP PGA |
| **3DP Projective GA only** | `#include "ga/ga_pga3dp_only.hpp"` | 🏃 Fast (30-40% faster) | Complete 3DP PGA |
| **All Euclidean GA** | `#include "ga/ga_ega.hpp"` | 🐌 Standard | Complete EGA (2D, 3D, 4D) |
| **All Projective GA** | `#include "ga/ga_pga.hpp"` | 🐌 Standard | Complete PGA (2DP, 3DP) |

## Detailed Usage Patterns

### 1. Minimal Include (Data Storage & Interfaces)

```cpp
// users header.hpp - Interface definitions only
#include "ga/ga_minimal.hpp"

class GeometryProcessor {
    hd::ga::vec3d position;
    hd::ga::bivec3d rotation;
    hd::ga::point3d target;
    
public:
    void setPosition(const hd::ga::vec3d& pos);
    hd::ga::vec3d getDirection() const;
};
```

```cpp
// users implementation.cpp - Include operations when needed
#include "ga/ga_ega3d_only.hpp"  // Now operations are available

using namespace hd::ga;
using namespace hd::ga::ega;

hd::ga::vec3d GeometryProcessor::getDirection() const {
    return normalize(target - position);  // Operations now available
}
```

**Benefits:**

- Fastest compilation for headers
- Clean separation of interface/implementation
- Perfect for template libraries and header-only interfaces

### 2. Dimension-Specific Includes

#### 2D Computer Graphics

```cpp
#include "ga/ga_ega2d_only.hpp"

using namespace hd::ga;
using namespace hd::ga::ega;

// 2D rotations using complex numbers (represented as even multivectors)
vec2d rotatePoint(const vec2d& point, double angle) {
    auto rotor = exp(pscalar2d{angle * 0.5});  // Half-angle rotation
    return (rotor * point * ~rotor).vector_part();
}
```

#### 3D Computer Graphics & Robotics

```cpp
#include "ga/ga_ega3d_only.hpp"

using namespace hd::ga;
using namespace hd::ga::ega;

// 3D rotations using quaternions (represented as even multivectors)
vec3d rotateVector(const vec3d& vector, const bivec3d& axis, double angle) {
    auto rotor = exp(axis * (angle * 0.5));
    return (rotor * vector * ~rotor).vector_part();
}
```

#### 2D Projective Geometry

```cpp
#include "ga/ga_pga2dp_only.hpp"

using namespace hd::ga;
using namespace hd::ga::pga;

// Point-line incidence and intersection
line2d createLine(const point2d& p1, const point2d& p2) {
    return wdg(p1,p2);  // Wedge product creates line through points
}

point2d intersectLines(const line2d& l1, const line2d& l2) {
    return rwdg(l1,l2);  // Regressive product finds intersection
}
```

#### 3D Scene Management & CAD

```cpp
#include "ga/ga_pga3dp_only.hpp"

using namespace hd::ga;
using namespace hd::ga::pga;

// Plane-line-point relationships
plane3d createPlane(const point3d& p1, const point3d& p2, const point3d& p3) {
    return wdg(wdg(p1,p2),p3);  // Three points define a plane
}

point3d intersectLineWithPlane(const line3d& line, const plane3d& plane) {
    return rwdg(line,plane);  // Regressive product finds intersection
}
```

### 3. Mixed Algebra Applications

When you need both EGA and PGA:

```cpp
// Option A: Use full includes
#include "ga/ga_ega.hpp"
#include "ga/ga_pga.hpp"

// Option B: Use selective includes
#include "ga/ga_ega3d_only.hpp"
#include "ga/ga_pga3dp_only.hpp"
```

### CMake Integration

```cmake

# For development with full debugging
target_compile_definitions(your_target PRIVATE 
    _HD_GA_EXTENDED_TEST_DIV_BY_ZERO
)
```

## Performance Recommendations

### Development Phase

- Use full includes (`ga_ega.hpp`, `ga_pga.hpp`) for maximum functionality
- Enable verbose error checking and type info for printing

### Production Phase

- Switch to selective includes for your specific use case

### Library Development

- Use `ga_minimal.hpp` in public headers
- Include operations only in implementation files
- Provide multiple include options for library users

## Common Include Patterns

### Game Engine (3D Graphics)

```cpp
// Core engine headers
#include "ga/ga_minimal.hpp"

// Rendering system
#include "ga/ga_ega3d_only.hpp"

// Physics/collision system  
#include "ga/ga_pga3dp_only.hpp"
```

### Computer Vision

```cpp
// 2D image processing
#include "ga/ga_ega2d_only.hpp"
#include "ga/ga_pga2dp_only.hpp"

// 3D reconstruction
#include "ga/ga_pga3dp_only.hpp"
```

### Mathematical Research

```cpp
// Full functionality for experimentation
#include "ga/ga_ega.hpp" 
#include "ga/ga_pga.hpp"
```

## Migration Guide

### From Full Includes

Replace:

```cpp
#include "ga/ga_ega.hpp"
```

With dimension-specific:

```cpp
#include "ga/ga_ega3d_only.hpp"  // If you only use 3D
```

### From Individual Headers

Replace multiple includes:

```cpp
#include "ga/detail/ga_mvec3d.hpp"
#include "ga/ga_ega3d_ops.hpp"
// ... many individual includes
```

With single selective include:

```cpp
#include "ga/ga_ega3d_only.hpp"
```

## Troubleshooting

### Compilation Errors

- **Missing operations**: Switch from `ga_minimal.hpp` to operation-enabled headers
- **Missing types**: Ensure you're including the right dimension (2D vs 3D vs 4D)
- **Missing geometric objects**: PGA types need `ga_pga_*_only.hpp` includes

### Performance Issues

- **Slow compilation**: Switch to selective includes
- **Template instantiation**: Use `ga_minimal.hpp` in headers

### Linking Issues

- **Multiple definitions**: Ensure consistent include patterns across translation units
- **Missing symbols**: Check that all required headers are included

---

*For more information, see the main library documentation and examples in the `ga_test/` directory.*
