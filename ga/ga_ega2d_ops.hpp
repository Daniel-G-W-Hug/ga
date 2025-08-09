#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_ega2d_ops_basics.hpp"   // ega2d ops basics
#include "ga_ega2d_ops_products.hpp" // ega2d ops products


namespace hd::ga::ega {

/////////////////////////////////////////////////////////////////////////////////////////
// provides ega2d functionality that is based on ega2d ops basics and products:
//
// - angle(), angle_to_re()         -> angle operations
// - exp()                          -> exponential function
// - get_rotor()                    -> provide a rotor
// - rotate()                       -> rotate object with rotor
// - project_onto(), reject_from()  -> projection and rejection
// - reflect_on(), reflect_on_vec() -> reflections
// - gs_orthogonal()                -> Gram-Schmidt-orthogonalization
// - is_congruent2d()               -> Same up to a scalar factor (is same subspace)
//
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// angle operations
////////////////////////////////////////////////////////////////////////////////

// return the angle between two vectors
// range of angle: -pi <= angle <= pi
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
inline std::common_type_t<T, U> angle(Vec2d<T> const& v1, Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    using std::numbers::pi;

    ctype nrm_prod = nrm(v1) * nrm(v2);
    hd::ga::detail::check_division_by_zero<T, U>(nrm_prod, "vector division");

    // std::clamp must be used to take care of numerical inaccuracies
    auto cos_angle = std::clamp(ctype(dot(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0));
    auto sin_angle = std::clamp(ctype(wdg(v1, v2)) / nrm_prod, ctype(-1.0), ctype(1.0));
    // wdg() in 2d contains magnitude and orientation, but works this easy only in 2d,
    // because it is already a scalar value
    // (for 3d to be as effective, the 3d vectors would need to be transformed
    //  to a plane, the angle measured w.r.t. to the pseudoscalar of the plane)

    // fmt::println("   c = {: .4f}, s = {: .4f}, wdg = {: .4f}, nrm_wdg = {: .4f}",
    //              cos_angle, sin_angle, wdg(v1, v2), nrm(wdg(v1, v2)));

    if (cos_angle >= 0.0) {
        // quadrant I or IV
        return std::asin(sin_angle);
    }
    else if (cos_angle < 0.0 && sin_angle >= 0.0) {
        // quadrant II
        return pi - std::asin(sin_angle);
    }
    else {
        // cos_angle < 0.0 && sin_angle < 0.0)
        // quadrant III
        return -pi - std::asin(sin_angle);
    }
}

// MVec2d_E<T> is used to model complex numbers using multivectors from the
//             even subalgebra
// returns the angle of the complex number w.r.t. the real axis
// range of angle: -pi <= angle <= pi
template <typename T>
    requires(std::floating_point<T>)
constexpr T angle_to_re(MVec2d_E<T> const& v)
{

    // TODO: replace by atan2 funcion

    using std::numbers::pi;
    if (v.c0 > 0.0) {
        // quadrant I & IV
        return std::atan(v.c1 / v.c0);
    }
    if (v.c0 < 0.0 && v.c1 >= 0.0) {
        // quadrant II
        return std::atan(v.c1 / v.c0) + pi;
    }
    if (v.c0 < 0.0 && v.c1 < 0.0) {
        // quadrant III
        return std::atan(v.c1 / v.c0) - pi;
    }
    if (v.c0 == 0.0) {
        // on y-axis
        if (v.c1 > 0.0) {
            return pi / 2.0;
        }
        if (v.c1 < 0.0) {
            return -pi / 2.0;
        }
    }
    return 0.0; // zero as input => define 0 as corresponding angle
}

////////////////////////////////////////////////////////////////////////////////
// 2d rotation operations
////////////////////////////////////////////////////////////////////////////////

// exponential function for setup of complex numbers and rotations
// as geometric multivectors with a scalar and a bivector part
//
// r = 1 is the vector length of the complex number in polar form
// theta is the bivector angle (i.e. a multiple of the bivector I_2d)
// such that uv = r exp(I_2d, theta) = a + I_2d b
// with r = |u| |v| = sqrt(a^2 + b^2) = 1
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> exp([[maybe_unused]] PScalar2d<T> I, U theta)
{
    // PScalar2d<T> is just used here for a unique overload of exp()
    // and to make the function signature similar to the 3d case
    using ctype = std::common_type_t<T, U>;
    return MVec2d_E<ctype>(Scalar2d<ctype>(std::cos(theta)),
                           PScalar2d<ctype>(std::sin(theta)));
}

//////////////////////////////////////////////////////////////////////////////////////////
// Inputs:
//         - a 2d pseudoscalar representing the plane of 2d space
//         - a rotation angle in the plane of 2d space
// Output:
//         - a rotor representing the requested rotation,
//           when applying the sandwich product with the rotor as in rotate(v,rotor)
//
//////////////////////////////////////////////////////////////////////////////////////////
//
// implemented here to make it formally the same with the 3d case (and potentially higher
// dimensional applications). In 2d the rotation can be directly expressed with less
// effort as
//
// exp(I_2d, -theta) * v = v * exp(I_2d, theta)  to express a 2d rotation of vector v by
//                                               the angle theta
//
//////////////////////////////////////////////////////////////////////////////////////////
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d_E<std::common_type_t<T, U>> get_rotor([[maybe_unused]] PScalar2d<T> I,
                                                       U theta)
{
    // PScalar2d<T> is just used here to be able to overload exp
    // and to make the function similar to the 3d case
    using ctype = std::common_type_t<T, U>;
    ctype half_angle = -0.5 * theta;
    return MVec2d_E<ctype>(Scalar2d<ctype>(std::cos(half_angle)),
                           PScalar2d<ctype>(std::sin(half_angle)));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> rotate(Vec2d<T> const& v,
                                                 MVec2d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(rotor * v * rev(rotor));
}

template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr MVec2d<std::common_type_t<T, U>> rotate(MVec2d<T> const& M,
                                                  MVec2d_E<U> const& rotor)
{
    using ctype = std::common_type_t<T, U>;
    return MVec2d<ctype>(rotor * M * rev(rotor));
}


////////////////////////////////////////////////////////////////////////////////
// Vec2d<T> projections, rejections and reflections
////////////////////////////////////////////////////////////////////////////////

// projection of v1 onto v2
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> project_onto(Vec2d<T> const& v1,
                                                       Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return ctype(dot(v1, v2)) * inv(v2);
}

// rejection of v1 from v2
// v_perp = gr1(wdg(v1,v2) * inv(v2))
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> reject_from(Vec2d<T> const& v1,
                                                      Vec2d<U> const& v2)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(v1 - project_onto(v1, v2));

    // works, but is more effort compared to solution via projection and vector difference
    // return Vec2d<ctype>(gr1(wdg(v1, v2) * inv(v2)));
}

////////////////////////////////////////////////////////////////////////////////
// reflections
////////////////////////////////////////////////////////////////////////////////
// reflect a vector u on a hyperplane B which is orthogonal (or dual) to vector nB
// choose: B = dual(nB) with nB as a normal to hyperplane B
//
// hyperplane: a n-1 dimensional subspace in a space of dimension n (e.g. a line through
// the origin in 2d space)
// orthogonal to vector b means the hyperplane is dual to b
// (call dual(nb) to create b from the normal vector nb)
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> reflect_on(Vec2d<T> const& u, Vec2d<U> const& B)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(-B * u * inv(B));
}

// reflect a vector u another vector
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr Vec2d<std::common_type_t<T, U>> reflect_on_vec(Vec2d<T> const& u,
                                                         Vec2d<U> const& b)
{
    using ctype = std::common_type_t<T, U>;
    return Vec2d<ctype>(b * u * inv(b));
}

////////////////////////////////////////////////////////////////////////////////
// Gram-Schmidt-Orthogonalization:
////////////////////////////////////////////////////////////////////////////////
// input:  two linear independent vectors u and v in 2d
// output: two orthogonal vectors with the first one being u and the second one a vector
// perpendicular to u in the orientation of v, both forming an orthogonal system
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::vector<Vec2d<std::common_type_t<T, U>>> gs_orthogonal(Vec2d<T> const& u,
                                                                     Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec2d<ctype>> basis;
    basis.push_back(u);
    basis.emplace_back(reject_from(v, u));
    return basis;
}

// input:  two linear independent vectors u and v in 2d
// output: two orthonormal vectors with the first one being normalize(u) and
// the second one a normalized vector perpendicular to u in the orientation of v,
// both forming an orthogonal system
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
constexpr std::vector<Vec2d<std::common_type_t<T, U>>> gs_orthonormal(Vec2d<T> const& u,
                                                                      Vec2d<U> const& v)
{
    using ctype = std::common_type_t<T, U>;
    std::vector<Vec2d<ctype>> basis;
    Vec2d<ctype> u_unitized{normalize(u)};
    basis.push_back(u_unitized);
    basis.emplace_back(normalize(reject_from(v, u_unitized)));
    return basis;
}

////////////////////////////////////////////////////////////////////////////////
// test congruence (same up to a scalar factor, i.e. representing same subspace)
////////////////////////////////////////////////////////////////////////////////

// For scalars: all non-zero scalars represent the same 0-dimensional subspace
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
bool is_congruent2d(Scalar2d<T> a, Scalar2d<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero scalars are congruent (represent the same 0-dimensional subspace)
    return true;
}

// For vectors: use unified A = k*B component-wise approach
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
bool is_congruent2d(Vec2d<T> const& a, Vec2d<U> const& b, value_t tolerance = eps)
{
    using ctype = std::common_type_t<T, U>;

    // Handle zero cases using component-wise check
    bool a_is_zero = (std::abs(a.x) < tolerance) && (std::abs(a.y) < tolerance);
    bool b_is_zero = (std::abs(b.x) < tolerance) && (std::abs(b.y) < tolerance);

    if (a_is_zero && b_is_zero) {
        return true; // Both are effectively zero
    }
    if (a_is_zero || b_is_zero) {
        return false; // Only one is zero
    }

    // Find scale factor k where a = k*b, checking all components
    ctype k = 0.0;
    bool k_found = false;

    // Find first non-zero component pair to establish k
    if (std::abs(b.x) > tolerance) {
        k = a.x / b.x;
        k_found = true;
    }
    else if (std::abs(b.y) > tolerance) {
        k = a.y / b.y;
        k_found = true;
    }

    if (!k_found) return false; // All components of b are zero, but a is not

    // Check if a = k*b for all components using relative tolerance
    value_t rel_tol = tolerance * std::max({std::abs(a.x), std::abs(a.y), std::abs(b.x),
                                            std::abs(b.y), value_t(1.0)});
    return (std::abs(a.x - k * b.x) < rel_tol) && (std::abs(a.y - k * b.y) < rel_tol);
}

// For pseudoscalars: all non-zero pseudoscalars in 2D represent the same subspace
template <typename T, typename U>
    requires(std::floating_point<T> && std::floating_point<U>)
bool is_congruent2d(PScalar2d<T> a, PScalar2d<U> b, value_t tolerance = eps)
{
    // Handle zero cases
    if (std::abs(T(a)) < tolerance && std::abs(U(b)) < tolerance) {
        return true; // Both are effectively zero
    }
    if (std::abs(T(a)) < tolerance || std::abs(U(b)) < tolerance) {
        return false; // Only one is zero
    }

    // All non-zero pseudoscalars in 2D are congruent (represent the full 2D space)
    return true;
}

} // namespace hd::ga::ega