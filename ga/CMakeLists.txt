# incrementally add test files needed

set(LIB_NAME ga)

set(HEADERS
    #
    detail/ga_fmt_support.hpp
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
    detail/type_t/ga_type_2d.hpp
    detail/type_t/ga_type_2dp.hpp
    detail/type_t/ga_type_3d.hpp
    detail/type_t/ga_type_3dp.hpp
    detail/type_t/ga_type_4d.hpp
    #
    detail/ga_mvec2d.hpp
    detail/ga_mvec2d_e.hpp
    #
    detail/ga_mvec3d.hpp
    detail/ga_mvec3d_e.hpp
    detail/ga_mvec3d_u.hpp
    #
    detail/ga_mvec4d.hpp
    detail/ga_mvec4d_e.hpp
    detail/ga_mvec4d_u.hpp
    #
    detail/ga_mvec2dp.hpp
    detail/ga_mvec2dp_e.hpp
    detail/ga_mvec2dp_u.hpp
    #
    detail/ga_mvec3dp.hpp
    detail/ga_mvec3dp_e.hpp
    detail/ga_mvec3dp_u.hpp
    #
    # user facing includes
    #
    ga_ega.hpp
    ga_pga.hpp
    #
    ga_usr_consts.hpp
    ga_usr_utilities.hpp
    ga_algebra.hpp
    ga_value_t.hpp
    #
    ga_ega_2d_ops.hpp
    ga_ega_3d_ops.hpp
    #
    detail/ga_pga_2dp_objects.hpp
    ga_pga_2dp_ops.hpp
    #
    detail/ga_pga_3dp_objects.hpp
    ga_pga_3dp_ops.hpp
    #
    ga_pga_common_ops.hpp
    )

add_library(${LIB_NAME} INTERFACE ${HEADERS})   #dep: fmt (must be provided by user!)
