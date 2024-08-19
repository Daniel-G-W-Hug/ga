# incrementally add test files needed

set(LIB_NAME ga)

set(HEADERS ga_ega.hpp
    detail/ga_algebra.hpp
    detail/ga_value_t.hpp detail/ga_Scalar_t.hpp
    detail/a_scalar.hpp
    detail/ga_pscalar2d.hpp detail/ga_vec2d.hpp detail/ga_mvec2d.hpp
    detail/ga_mvec2d_e.hpp detail/ga_ega_2d_ops.hpp
    detail/ga_pscalar3d.hpp detail/ga_vec3d.hpp detail/ga_mvec3d.hpp
    detail/ga_mvec3d_e.hpp detail/ga_mvec3d_u.hpp detail/ga_bivec3d.hpp
    detail/ga_ega_3d_ops.hpp
    ga_usr_consts.hpp ga_usr_types.hpp ga_usr_utilities.hpp
    )

add_library(${LIB_NAME} INTERFACE ${HEADERS})   #dep: fmt (must be provided by user)
