#pragma once

// author: Daniel Hug, 2024

// strong type inspired by Jonathan Boccara's blog:
// https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
//
// and by the blog of foonathan
// https://www.foonathan.net/2016/10/strong-typedefs/


#include <cmath>     // std::abs
#include <concepts>  // std::floating_point
#include <iostream>  // std::ostream
#include <limits>    // std::numeric_limits
#include <stdexcept> // std::runtime_error
#include <string>    // std::string, std::to_string

// provide common type for Scalar<T>, PScalar2d<T>, PScalar3d<T>

namespace hd::ga {

template <typename T, typename Tag>
    requires(std::floating_point<T>)
class Scalar_t {

  public:

    Scalar_t() = default; // default initialization
    explicit Scalar_t(T const& val) : value(val) {}
    explicit Scalar_t(T&& val) noexcept(std::is_nothrow_move_constructible<T>::value) :
        value(std::move(val))
    {
    }

    // these operators bring implicit conversion to the underlying type T
    // (hint: can be suppressed with marking the operator explicit)
    operator T&() noexcept { return value; }
    operator T const&() const noexcept { return value; }

    friend void swap(Scalar_t& lhs, Scalar_t& rhs) noexcept
    {
        using std::swap;
        swap(static_cast<T&>(lhs), static_cast<T&>(rhs));
    }

  private:

    T value{};
};

////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag> core operations
////////////////////////////////////////////////////////////////////////////////

// unary minus
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr Scalar_t<T, Tag> operator-(Scalar_t<T, Tag> v)
{
    return Scalar_t<T, Tag>(-T(v));
}

// adding pseudoscalars
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator+(Scalar_t<T, Tag> v1,
                                                                   Scalar_t<U, Tag> v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(v1) + ctype(v2));
}

// substracting pseudoscalars
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator-(Scalar_t<T, Tag> v1,
                                                                   Scalar_t<U, Tag> v2)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(v1) - ctype(v2));
}

// multiply a pseudoscalar with a scalar (in both constellations)
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator*(Scalar_t<T, Tag> v,
                                                                   U s)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(ctype(v) * s);
}

template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator*(T s,
                                                                   Scalar_t<U, Tag> v)
{
    using ctype = std::common_type_t<T, U>;
    return Scalar_t<ctype, Tag>(s * ctype(v));
}

// devide a pseudoscalar by a scalar
template <typename T, typename U, typename Tag>
    requires(std::floating_point<T> && std::floating_point<U>)
inline constexpr Scalar_t<std::common_type_t<T, U>, Tag> operator/(Scalar_t<T, Tag> v,
                                                                   U s)
{
    using ctype = std::common_type_t<T, U>;
    if (std::abs(s) < 5.0 * std::numeric_limits<T>::epsilon()) {
        throw std::runtime_error("scalar too small, division by zero" +
                                 std::to_string(s) + "\n");
    }
    ctype inv = ctype(1.0) / s; // for multiplicaton with inverse value
    return Scalar_t<ctype, Tag>(ctype(v) * inv);
}

// return squared magnitude
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T sq_nrm(Scalar_t<T, Tag> s)
{
    return T(s) * T(s);
}

// return absolute magnitude
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T nrm(Scalar_t<T, Tag> s)
{
    return std::abs(T(s));
}

// return the scalar value as value_t (for use in scripting)
template <typename T, typename Tag>
    requires(std::floating_point<T>)
inline constexpr T to_val(Scalar_t<T, Tag> s)
{
    return T(s);
}

////////////////////////////////////////////////////////////////////////////////
// Scalar_t<T, Tag> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T, typename Tag>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Scalar_t<T, Tag> v)
{
    os << "(" << T(v) << ")";
    return os;
}

} // namespace hd::ga