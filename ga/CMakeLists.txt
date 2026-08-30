# incrementally add test files needed

set(LIB_NAME ga)

set(HEADERS
    #
    detail/fmt/ga_fmt_core.hpp
    detail/fmt/ga_fmt_ega.hpp
    detail/fmt/ga_fmt_pga.hpp
    detail/fmt/ga_fmt_sta.hpp
    detail/ga_fmt_support.hpp
    #
    detail/ga_error_handling.hpp
    detail/ga_core_types.hpp
    detail/ga_foundation.hpp
    detail/ga_ega_types.hpp
    detail/ga_pga_types.hpp
    detail/ga_sta_types.hpp
    detail/ga_solver.hpp
    detail/ga_stencil.hpp
    #
    detail/type_t/ga_scalar_t.hpp
    detail/type_t/ga_vec2_t.hpp
    detail/type_t/ga_vec3_t.hpp
    detail/type_t/ga_vec4_t.hpp
    detail/type_t/ga_bvec6_t.hpp
    detail/type_t/ga_mvec2_t.hpp
    detail/type_t/ga_mvec4_t.hpp
    detail/type_t/ga_mvec8_t.hpp
    detail/type_t/ga_mvec16_t.hpp
    detail/type_t/ga_type_tags.hpp
    detail/type_t/ga_type2d.hpp
    detail/type_t/ga_type2dp.hpp
    detail/type_t/ga_type3d.hpp
    detail/type_t/ga_type3dp.hpp
    detail/type_t/ga_type4ds.hpp
    #
    detail/type_t/ga_mvec2d.hpp
    detail/type_t/ga_mvec2d_e.hpp
    #
    detail/type_t/ga_mvec3d.hpp
    detail/type_t/ga_mvec3d_e.hpp
    detail/type_t/ga_mvec3d_u.hpp
    #
    detail/type_t/ga_mvec2dp.hpp
    detail/type_t/ga_mvec2dp_e.hpp
    detail/type_t/ga_mvec2dp_u.hpp
    detail/type_t/ga_dualnum2dp.hpp
    #
    detail/type_t/ga_mvec3dp.hpp
    detail/type_t/ga_mvec3dp_e.hpp
    detail/type_t/ga_mvec3dp_u.hpp
    detail/type_t/ga_dualnum3dp.hpp
    #
    detail/type_t/ga_mvec4ds.hpp
    detail/type_t/ga_mvec4ds_e.hpp
    detail/type_t/ga_mvec4ds_u.hpp
    #
    # user facing includes
    #
    ga_ega.hpp
    ga_pga.hpp
    ga_sta.hpp
    #
    ga_usr_consts.hpp
    ga_usr_types.hpp
    ga_usr_types_mechanics.hpp
    ga_usr_utilities.hpp
    ga_algebra.hpp
    ga_value_t.hpp
    #
    ga_ega2d_ops_basics.hpp
    ga_ega2d_ops_products.hpp
    ga_ega2d_ops.hpp
    ga_ega3d_ops_basics.hpp
    ga_ega3d_ops_products.hpp
    ga_ega3d_ops.hpp
    #
    ga_pga2dp_ops_basics.hpp
    ga_pga2dp_ops_products.hpp
    ga_pga2dp_ops.hpp
    ga_pga2dp_ops_mechanics.hpp
    ga_pga2dp_ops_constraints.hpp
    #
    ga_pga3dp_ops_basics.hpp
    ga_pga3dp_ops_products.hpp
    ga_pga3dp_ops.hpp
    ga_pga3dp_ops_mechanics.hpp
    ga_pga3dp_ops_constraints.hpp
    #
    ga_sta4ds_ops_basics.hpp
    ga_sta4ds_ops_products.hpp
    ga_sta4ds_ops.hpp
    )

add_library(${LIB_NAME} INTERFACE ${HEADERS})   #dep: fmt (must be provided by user!)

# Root-relative includes ("ga/ga_ega.hpp") resolve against the PARENT of this dir (the
# repo root). Export it as an INTERFACE include so every target that links `ga` gets it
# automatically -- crucially cross-scope, e.g. a target added by an enclosing build, where
# a GA_ROOT set inside this repository is not visible.
target_include_directories(${LIB_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)
