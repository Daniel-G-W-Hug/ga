#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace codegen {

// How a result type's coefficients are laid out in its constructor.
enum class StorageKind {
    SingleValue, // Scalar2d, PScalar3dp, ... -- one value
    Named,       // Vec2d(x,y), BiVec3dp(vx,vy,vz,mx,my,mz), ... -- named components
    Indexed,     // MVec2d(c0,c1,c2,c3), MVec2d_E(c0,c1) -- positional c0..cN
    Composite,   // MVec3d_E(Scalar3d(...), BiVec3d(...)) -- nested sub-types
};

struct TypeInfo {
    // C++ type name without template parameters: "Vec3dp", "MVec2dp_E", "Scalar3d"
    std::string cpp_type;

    StorageKind storage;

    // Basis indices (positions in algebra.basis) that this type spans, in
    // constructor-argument order. For Composite, this is the flat union of all
    // sub-types' basis indices, also in constructor order.
    std::vector<std::size_t> basis_indices;

    // For Composite only: ordered list of sub-type filter names.
    // e.g. mv_e (PGA2DP) -> {"s", "bivec"}
    // e.g. mv_e (PGA3DP) -> {"s", "bivec", "ps"}
    std::vector<std::string> sub_filter_names;

    // True for Scalar/PScalar (passed by value), false for everything else (const&).
    bool pass_by_value = false;
};

class TypeRegistry {
  public:

    explicit TypeRegistry(std::string const& algebra_name);

    // Algebra name this registry was built for ("ega2d", "ega3d", ...).
    std::string const& algebra() const { return algebra_; }

    // True if this filter has a registered TypeInfo.
    bool has(std::string const& filter_name) const;

    // Lookup TypeInfo by filter name. Throws std::out_of_range if missing.
    TypeInfo const& get(std::string const& filter_name) const;

  private:

    std::string algebra_;
    std::map<std::string, TypeInfo> types_;
};

} // namespace codegen
