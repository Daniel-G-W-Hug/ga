# incrementally add test files needed

set(LIB_NAME ga)

set(HEADERS
    #
    detail/fmt/ga_fmt_core.hpp
    detail/fmt/ga_fmt_cga.hpp
    detail/fmt/ga_fmt_ega.hpp
    detail/fmt/ga_fmt_pga.hpp
    detail/fmt/ga_fmt_sta.hpp
    detail/fmt/ga_fmt_geodesics.hpp
    detail/fmt/ga_fmt_mechanics.hpp
    detail/ga_fmt_support.hpp
    #
    detail/ga_error_handling.hpp
    detail/ga_core_types.hpp
    detail/ga_foundation.hpp
    detail/ga_ega_types.hpp
    detail/ga_pga_types.hpp
    detail/ga_sta_types.hpp
    detail/ga_cga_types.hpp
    detail/ga_solver.hpp
    detail/ga_stencil.hpp
    #
    detail/type_t/ga_scalar_t.hpp
    detail/type_t/ga_vec2_t.hpp
    detail/type_t/ga_vec3_t.hpp
    detail/type_t/ga_vec4_t.hpp
    detail/type_t/ga_vec5_t.hpp
    detail/type_t/ga_bvec6_t.hpp
    detail/type_t/ga_bvec10_t.hpp
    detail/type_t/ga_mvec2_t.hpp
    detail/type_t/ga_mvec4_t.hpp
    detail/type_t/ga_mvec8_t.hpp
    detail/type_t/ga_mvec16_t.hpp
    detail/type_t/ga_mvec32_t.hpp
    detail/type_t/ga_type_tags.hpp
    detail/type_t/ga_type2d.hpp
    detail/type_t/ga_type2dc.hpp
    detail/type_t/ga_type2dp.hpp
    detail/type_t/ga_type3d.hpp
    detail/type_t/ga_type3dc.hpp
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
    detail/type_t/ga_mvec2dc.hpp
    detail/type_t/ga_mvec2dc_e.hpp
    detail/type_t/ga_mvec2dc_u.hpp
    #
    detail/type_t/ga_mvec3dc.hpp
    detail/type_t/ga_mvec3dc_e.hpp
    detail/type_t/ga_mvec3dc_u.hpp
    #
    # user facing includes
    #
    ga_ega.hpp
    ga_pga.hpp
    ga_sta.hpp
    ga_cga.hpp
    #
    ga_usr_consts.hpp
    ga_usr_fd.hpp
    ga_usr_geodesics.hpp
    ga_usr_types.hpp
    ga_usr_types_mechanics.hpp
    ga_usr_utilities.hpp
    ga_algebra.hpp
    ga_value_t.hpp
    #
    ga_ega2d_ops_basics.hpp
    ga_ega2d_ops_products.hpp
    ga_ega2d_ops.hpp
    ga_ega2d_ops_calculus.hpp
    ga_ega3d_ops_basics.hpp
    ga_ega3d_ops_products.hpp
    ga_ega3d_ops.hpp
    ga_ega3d_ops_calculus.hpp
    #
    ga_pga2dp_ops_basics.hpp
    ga_pga2dp_ops_products.hpp
    ga_pga2dp_ops.hpp
    ga_pga2dp_ops_mechanics.hpp
    ga_pga2dp_ops_constraints.hpp
    ga_pga2dp_ops_contact.hpp
    #
    ga_pga3dp_ops_basics.hpp
    ga_pga3dp_ops_products.hpp
    ga_pga3dp_ops.hpp
    ga_pga3dp_ops_mechanics.hpp
    ga_pga3dp_ops_constraints.hpp
    ga_pga3dp_ops_contact.hpp
    #
    ga_sta4ds_ops_basics.hpp
    ga_sta4ds_ops_products.hpp
    ga_sta4ds_ops.hpp
    ga_sta4ds_ops_calculus.hpp
    #
    ga_cga2dc_ops_basics.hpp
    ga_cga2dc_ops_products.hpp
    ga_cga2dc_ops.hpp
    #
    ga_cga3dc_ops_basics.hpp
    ga_cga3dc_ops_products.hpp
    ga_cga3dc_ops.hpp
    )

add_library(${LIB_NAME} INTERFACE ${HEADERS})   #dep: fmt (must be provided by user!)

# Root-relative includes ("ga/ga_ega.hpp") resolve against the PARENT of this dir (the
# repo root). Export it as an INTERFACE include so every target that links `ga` gets it
# automatically -- crucially cross-scope, e.g. a target added by an enclosing build, where
# a GA_ROOT set inside this repository is not visible.
target_include_directories(${LIB_NAME} INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/..)

# THE DIVISION GUARD: ON BY DEFAULT, OFF PER TARGET WHERE SPEED IS THE GOAL.
#
# With _HD_GA_EXTENDED_TEST_DIV_BY_ZERO defined, every library division, normalization
# and unitization checks its divisor against a safe epsilon and THROWS on a near-zero
# (ga/detail/ga_error_handling.hpp), so a degenerate element fails where it happens
# instead of surfacing as a NaN three functions later. That is a development safety
# gate, and it is the default for every consumer: the definition is an INTERFACE usage
# requirement of `ga`, so linking the library is enough to get it -- in this tree and in
# any enclosing build that adds it via add_subdirectory.
#
# It is not free. The check is a compare-and-branch, but a throw inside an innermost
# loop can stop the compiler inlining and vectorizing the loop around it. With the
# throw INLINE that was measured at 1.63x on a pair-summation control loop (byte-
# identical output); since the throw has been outlined into a cold, non-returning helper
# the remaining cost is what the two benches below report -- build and run them:
#
#     cmake --build . --target ga_bench_division_guard ga_bench_division_guard_off
#     ./ga_test/ga_bench_division_guard && ./ga_test/ga_bench_division_guard_off
#
# The same kernel (a pairwise inverse-cube vector sum, the shape of any N-body-like
# inner loop), built once with the guard and once without through the function below.
# The ratio between the two ns/op figures is the cost on THIS compiler; the number
# recorded in BUILD_USAGE.md ("Build Options") is the last measurement.
#
# TO GO FOR SPEED: call ga_division_guard(<your target> OFF) on the targets that need
# it -- and only on those. The choice is per TRANSLATION UNIT (the library is header-
# only, so the check is instantiated in the consumer), hence per target, never per
# call site; a target that keeps the guard while linking the same library is exactly
# what "guarded application over a fast library" looks like. Two binaries are needed to
# compare, which is why the bench is a pair: the same inline templates with two
# different bodies in ONE binary would be an ODR violation.
#
# Mechanism: the generator expression reads the CONSUMING target's own
# GA_DIVISION_GUARD_OFF property (a TARGET_PROPERTY expression without a target name
# evaluates on the target being compiled), so the default needs nothing from the
# consumer and the opt-out is a property set on the consumer alone.
define_property(TARGET PROPERTY GA_DIVISION_GUARD_OFF
    BRIEF_DOCS "TRUE removes the ga division guard from this target's translation units"
    FULL_DOCS "Set through ga_division_guard(<target> OFF); see ga/CMakeLists.txt.")
target_compile_definitions(${LIB_NAME} INTERFACE
    $<$<NOT:$<BOOL:$<TARGET_PROPERTY:GA_DIVISION_GUARD_OFF>>>:_HD_GA_EXTENDED_TEST_DIV_BY_ZERO>)

# ga_division_guard(<target> ON|OFF): keep (default) or drop the guard for one target.
function(ga_division_guard target state)
    if(state)
        set_property(TARGET ${target} PROPERTY GA_DIVISION_GUARD_OFF FALSE)
    else()
        set_property(TARGET ${target} PROPERTY GA_DIVISION_GUARD_OFF TRUE)
    endif()
endfunction()
