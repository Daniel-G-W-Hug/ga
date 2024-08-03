# incrementally add test files needed

set(LIB_NAME ga)

set(HEADERS ga.hpp
    ga_cfg_algebra.hpp
    ga_cfg_value_t.hpp ga_cfg_strong_t.hpp
    ga_cfg_scalar.hpp
    ga_cfg_pscalar2d.hpp ga_cfg_vec2d.hpp ga_cfg_mvec2d.hpp ga_cfg_mvec2d_e.hpp
    ga_cfg_2d_ops.hpp
    ga_cfg_pscalar3d.hpp ga_cfg_vec3d.hpp ga_cfg_mvec3d.hpp ga_cfg_mvec3d_e.hpp
    ga_cfg_mvec3d_u.hpp ga_cfg_bivec3d.hpp
    ga_cfg_3d_ops.hpp
    ga_usr_consts.hpp ga_usr_types.hpp ga_usr_utilities.hpp
    )

add_lnibrary(${LIB_NAME} INTERFACE ${HEADERS})   #dep: fmt (must be provided by user)
