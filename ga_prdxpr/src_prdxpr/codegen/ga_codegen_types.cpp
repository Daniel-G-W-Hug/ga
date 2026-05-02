// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "codegen/ga_codegen_types.hpp"

#include <stdexcept>

namespace codegen {

namespace {

// Build the mv_e/mv_u basis_indices for a Composite type by concatenating the
// basis_indices of its sub-types. Sub-types must already be inserted into `m`.
TypeInfo make_composite(std::string const& cpp_type,
                        std::vector<std::string> const& sub_filters,
                        std::map<std::string, TypeInfo> const& m)
{
    TypeInfo info;
    info.cpp_type = cpp_type;
    info.storage = StorageKind::Composite;
    info.sub_filter_names = sub_filters;
    for (auto const& sub : sub_filters) {
        auto const& sub_info = m.at(sub);
        info.basis_indices.insert(info.basis_indices.end(),
                                  sub_info.basis_indices.begin(),
                                  sub_info.basis_indices.end());
    }
    return info;
}

TypeInfo make_single(std::string const& cpp_type, std::size_t basis_index)
{
    return TypeInfo{.cpp_type = cpp_type,
                    .storage = StorageKind::SingleValue,
                    .basis_indices = {basis_index},
                    .sub_filter_names = {},
                    .pass_by_value = true};
}

TypeInfo make_named(std::string const& cpp_type,
                    std::vector<std::size_t> const& basis_indices)
{
    return TypeInfo{.cpp_type = cpp_type,
                    .storage = StorageKind::Named,
                    .basis_indices = basis_indices,
                    .sub_filter_names = {},
                    .pass_by_value = false};
}

TypeInfo make_indexed(std::string const& cpp_type,
                      std::vector<std::size_t> const& basis_indices)
{
    return TypeInfo{.cpp_type = cpp_type,
                    .storage = StorageKind::Indexed,
                    .basis_indices = basis_indices,
                    .sub_filter_names = {},
                    .pass_by_value = false};
}

// EGA2D: basis = {1, e1, e2, e12} (4 elements)
// MVec2d_E uses FLAT (Indexed) construction MVec2d_E<ctype>(c0, c1) -- diverges from
// the other algebras which use nested Composite construction for their _E/_U types.
std::map<std::string, TypeInfo> build_ega2d()
{
    std::map<std::string, TypeInfo> m;
    m["s"] = make_single("Scalar2d", 0);
    m["vec"] = make_named("Vec2d", {1, 2});
    m["ps"] = make_single("PScalar2d", 3);
    m["mv_e"] = make_indexed("MVec2d_E", {0, 3});
    m["mv"] = make_indexed("MVec2d", {0, 1, 2, 3});
    return m;
}

// EGA3D: basis = {1, e1, e2, e3, e23, e31, e12, e123} (8 elements)
std::map<std::string, TypeInfo> build_ega3d()
{
    std::map<std::string, TypeInfo> m;
    m["s"] = make_single("Scalar3d", 0);
    m["vec"] = make_named("Vec3d", {1, 2, 3});
    m["bivec"] = make_named("BiVec3d", {4, 5, 6});
    m["ps"] = make_single("PScalar3d", 7);
    m["mv_e"] = make_composite("MVec3d_E", {"s", "bivec"}, m);
    m["mv_u"] = make_composite("MVec3d_U", {"vec", "ps"}, m);
    m["mv"] = make_indexed("MVec3d", {0, 1, 2, 3, 4, 5, 6, 7});
    return m;
}

// PGA2DP: basis = {1, e1, e2, e3, e23, e31, e12, e321} (8 elements)
std::map<std::string, TypeInfo> build_pga2dp()
{
    std::map<std::string, TypeInfo> m;
    m["s"] = make_single("Scalar2dp", 0);
    m["vec"] = make_named("Vec2dp", {1, 2, 3});
    m["bivec"] = make_named("BiVec2dp", {4, 5, 6});
    m["ps"] = make_single("PScalar2dp", 7);
    m["mv_e"] = make_composite("MVec2dp_E", {"s", "bivec"}, m);
    m["mv_u"] = make_composite("MVec2dp_U", {"vec", "ps"}, m);
    m["mv"] = make_indexed("MVec2dp", {0, 1, 2, 3, 4, 5, 6, 7});
    return m;
}

// PGA3DP: basis = {1, e1, e2, e3, e4, e41, e42, e43, e23, e31, e12, e423, e431, e412,
//                  e321, e1234} (16 elements)
// BiVec3dp components: (vx, vy, vz, mx, my, mz) -- non-uniform naming.
std::map<std::string, TypeInfo> build_pga3dp()
{
    std::map<std::string, TypeInfo> m;
    m["s"] = make_single("Scalar3dp", 0);
    m["vec"] = make_named("Vec3dp", {1, 2, 3, 4});
    m["bivec"] = make_named("BiVec3dp", {5, 6, 7, 8, 9, 10});
    m["trivec"] = make_named("TriVec3dp", {11, 12, 13, 14});
    m["ps"] = make_single("PScalar3dp", 15);
    m["mv_e"] = make_composite("MVec3dp_E", {"s", "bivec", "ps"}, m);
    m["mv_u"] = make_composite("MVec3dp_U", {"vec", "trivec"}, m);
    m["mv"] = make_indexed("MVec3dp",
                           {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15});
    return m;
}

} // namespace

TypeRegistry::TypeRegistry(std::string const& algebra_name) : algebra_(algebra_name)
{
    if (algebra_name == "ega2d") {
        types_ = build_ega2d();
    }
    else if (algebra_name == "ega3d") {
        types_ = build_ega3d();
    }
    else if (algebra_name == "pga2dp") {
        types_ = build_pga2dp();
    }
    else if (algebra_name == "pga3dp") {
        types_ = build_pga3dp();
    }
    else {
        // sta4d intentionally omitted: code generation for it lands once the user
        // supplies the configuration. Falls through to empty registry.
        throw std::invalid_argument(
            "TypeRegistry: no registry available for algebra '" + algebra_name + "'");
    }
}

bool TypeRegistry::has(std::string const& filter_name) const
{
    return types_.find(filter_name) != types_.end();
}

TypeInfo const& TypeRegistry::get(std::string const& filter_name) const
{
    auto it = types_.find(filter_name);
    if (it == types_.end()) {
        throw std::out_of_range("TypeRegistry(" + algebra_ +
                                "): no type for filter '" + filter_name + "'");
    }
    return it->second;
}

} // namespace codegen
