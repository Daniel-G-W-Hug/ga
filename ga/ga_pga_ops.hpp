#pragma once

// author: Daniel Hug, 2024 & 2025

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

} // namespace hd::ga::pga