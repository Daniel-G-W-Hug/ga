#pragma once

// author: Daniel Hug, 2024 & 2025

#include <algorithm> // std::clamp

// PGA-specific operations are in namespace hd::ga::pga
#include "ga_pga_2dp_ops.hpp" // pga 2d operations
#include "ga_pga_3dp_ops.hpp" // pga 3d operations

////////////////////////////////////////////////////////////////////////////////
// commonly used operations across pga2dp and pga3dp
////////////////////////////////////////////////////////////////////////////////

namespace hd::ga::pga {

////////////////////////////////////////////////////////////////////////////////
// Projective contractions (same for 2dp and 3dp):
//
// bulk_contraction(a,b) = rwdg(a, bulk_dual(b))
//
// weight_contraction(a,b) = rwdg(a, weight_dual(b))
//
// The contraction subracts the grades of the operands.
//
// When the metric is the identity, these two contractions are identical and are
// called the contraction (they are identical with the right contraction).
// When the metric is degenerate they produce different results.
//
// In general a contraciton throws away parts that are perpendicular to each other.
// The result of bulk_contraction(B,v) lies in B and is perpendicular to v.
////////////////////////////////////////////////////////////////////////////////

template <typename arg1, typename arg2>
decltype(auto) bulk_contraction(arg1&& a, arg2&& b)
{
    // return rwdg(a, bulk_dual(b));
    return rwdg(std::forward<arg1>(a), bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) weight_contraction(arg1&& a, arg2&& b)
{
    // return rwdg(a, weight_dual(b));
    return rwdg(std::forward<arg1>(a), weight_dual(std::forward<arg2>(b)));
}


////////////////////////////////////////////////////////////////////////////////
// Projective expansions (same for 2dp and 3dp):
//
// bulk_expansion(a,b) = wdg(a, bulk_dual(b))       (dual to weight_contraction)
//
// weight_expansion(a,b) = wdg(a, weight_dual(b))   (dual to bulk_contraction)
//
// The expansion subtracts the antigrades of the objects.
//
// When the metric is the identity, these two expansions are identical and are
// called the expansion (they are identical with the left complement of the left
// contraction). When the metric is degenerate they produce different results.
//
// The expansion takes the parts of one object that are parallel to the other
// object and combines them with the space that is perpendicular to that other
// dualized object.
////////////////////////////////////////////////////////////////////////////////

template <typename arg1, typename arg2> decltype(auto) bulk_expansion(arg1&& a, arg2&& b)
{
    // return wdg(a, bulk_dual(b));
    return wdg(std::forward<arg1>(a), bulk_dual(std::forward<arg2>(b)));
}

template <typename arg1, typename arg2>
decltype(auto) weight_expansion(arg1&& a, arg2&& b)
{
    // return wdg(a, weight_dual(b));
    return wdg(std::forward<arg1>(a), weight_dual(std::forward<arg2>(b)));
}

////////////////////////////////////////////////////////////////////////////////
// Projections (same for 2dp and 3dp):
//
// ortho_proj(a, b)     = rwdg(b, weight_expansion(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' contained in b)
// REQUIRES: gr(a) < gr(b)
//
//
// central_proj(a, b)        = rwdg(b, bulk_expansion(a, b) )
// (a projected centrally towards origin onto b, effectively creating a new a'
// contained in b)
// REQUIRES: gr(a) < gr(b)
//
// ortho_antiproj(a, b) = wdg(b, weight_contraction(a, b) )
// (a projected orthogonally onto b, effectively creating a new a' containing b)
// REQUIRES: gr(a) > gr(b)
//
////////////////////////////////////////////////////////////////////////////////

template <typename arg1, typename arg2> decltype(auto) ortho_proj(arg1&& a, arg2&& b)
{

    // REQUIRES: gr(a) < gr(b), or does not compile!
    // project the smaller grade object onto to larger grade object
    return rwdg(std::forward<arg2>(b),
                weight_expansion(std::forward<arg1>(a), std::forward<arg2>(b)));
}

template <typename arg1, typename arg2> decltype(auto) central_proj(arg1&& a, arg2&& b)
{
    // REQUIRES: gr(a) < gr(b), or does not compile!
    // project the smaller grade object onto to larger grade object
    return rwdg(std::forward<arg2>(b),
                bulk_expansion(std::forward<arg1>(a), std::forward<arg2>(b)));
}

template <typename arg1, typename arg2> decltype(auto) ortho_antiproj(arg1&& a, arg2&& b)
{
    return wdg(std::forward<arg2>(b),
               weight_contraction(std::forward<arg1>(a), std::forward<arg2>(b)));
}


// return the point nearest to the origin by projecting the origin onto line
template <typename arg1> decltype(auto) support2dp(arg1&& a)
{
    // REQUIRES: a line (BiVec2dp) as argument
    return ortho_proj(origin_2dp, std::forward<arg1>(a));
}

// return the point nearest to the origin by projecting the origin onto plane or line
template <typename arg1> decltype(auto) support3dp(arg1&& a)
{
    // REQUIRES: a line (BiVec3dp) or a plane (TriVec3dp) as argument
    return ortho_proj(origin_3dp, std::forward<arg1>(a));
}

////////////////////////////////////////////////////////////////////////////////
// angle operations 2dp
////////////////////////////////////////////////////////////////////////////////

// return the angle between of two vectors, i.e. directions to points at infinity
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(Vec2dp<T> const& v1, Vec2dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    if ((weight_nrm_sq(v1) != 0.0) || (weight_nrm_sq(v2) != 0.0)) {
        // the angle between points not at infinity or points not at infinity and a
        // direction towards infinity is defined as zero
        return 0.0;
    }

    // angle is defined only between directions towards points at infinity

    ctype nrm_prod = bulk_nrm(v1) * bulk_nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle " +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between two bivectors, i.e. between two lines
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(BiVec2dp<T> const& B1,
                                                BiVec2dp<U> const& B2)
{
    //// TODO: angle not consistently working in all cases yet -> define tests

    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(weight_contraction(B1, B2));
    ctype nrm_prod = ctype(weight_nrm(B1) * weight_nrm(B2));
    // fmt::println("contr: {}, nrm_prod = {}", contr, nrm_prod);
    if (nrm_prod != 0.0) {
        // fmt::println("angle - using nrm_prod");
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        // fmt::println("angle - not using nrm_prod");
        // fmt::println("after clamp: {}", std::clamp(contr, ctype(-1.0), ctype(1.0)));
        // fmt::println("acos: {}", std::acos(std::clamp(contr, ctype(-1.0),
        // ctype(1.0))));
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}

////////////////////////////////////////////////////////////////////////////////
// angle operations 3dp
////////////////////////////////////////////////////////////////////////////////

// return the angle between of two vectors, i.e. directions to points at infinity
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(Vec3dp<T> const& v1, Vec3dp<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;

    if ((weight_nrm_sq(v1) != 0.0) || (weight_nrm_sq(v2) != 0.0)) {
        // the angle between points not at infinity or points not at infinity and a
        // direction towards infinity is defined as zero
        return 0.0;
    }

    // angle is defined only between directions towards points a infinity

    ctype nrm_prod = bulk_nrm(v1) * bulk_nrm(v2);
#if defined(_HD_GA_EXTENDED_TEST_DIV_BY_ZERO)
    if (nrm_prod < std::numeric_limits<ctype>::epsilon()) {
        throw std::runtime_error(
            "vector norm product too small for calculation of angle " +
            std::to_string(nrm_prod) + "\n");
    }
#endif
    // std::clamp must be used to take care of numerical inaccuracies
    return std::acos(std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0)));
}

// return the angle between two bivectors, i.e. lines
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(BiVec3dp<T> const& B1,
                                                BiVec3dp<U> const& B2)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(weight_contraction(B1, B2));
    ctype nrm_prod = ctype(weight_nrm(B1) * weight_nrm(B2));
    // fmt::println("contr: {}, nrm_prod = {}", contr, nrm_prod);
    if (nrm_prod != 0.0) {
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}

// return the angle between a trivector and a bivector, i.e. a plane and a line
// range of angle: 0 <= angle <= pi/2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(TriVec3dp<T> const& t,
                                                BiVec3dp<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(bulk_nrm(weight_contraction(t, B)));
    ctype nrm_prod = ctype(weight_nrm(t) * weight_nrm(B));
    // fmt::println("contr: {}, nrm_prod = {}", contr, nrm_prod);
    if (nrm_prod != 0.0) {
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}

// return the angle between a trivector and a bivector, i.e. a plane and a line
// range of angle: 0 <= angle <= pi/2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(BiVec3dp<T> const& B,
                                                TriVec3dp<U> const& t)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(bulk_nrm(weight_contraction(t, B)));
    ctype nrm_prod = ctype(weight_nrm(t) * weight_nrm(B));
    // fmt::println("contr: {}, nrm_prod = {}", contr, nrm_prod);
    if (nrm_prod != 0.0) {
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}

// return the angle between two trivectors, i.e. two planes
// range of angle: 0 <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr std::common_type_t<T, U> angle(TriVec3dp<T> const& t1,
                                                TriVec3dp<U> const& t2)
{
    using ctype = std::common_type_t<T, U>;
    ctype contr = ctype(weight_contraction(t1, t2));
    ctype nrm_prod = ctype(weight_nrm(t1) * weight_nrm(t2));
    // fmt::println("contr: {}, nrm_prod = {}", contr, nrm_prod);
    if (nrm_prod != 0.0) {
        return std::acos(std::clamp(contr / nrm_prod, ctype(-1.0), ctype(1.0)));
    }
    else {
        return std::acos(std::clamp(contr, ctype(-1.0), ctype(1.0)));
    }
}

} // namespace hd::ga::pga