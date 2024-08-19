#pragma once

// author: Daniel Hug, 2024

#include <concepts> // std::floating_point
#include <iostream>

#include "ga_value_t.hpp"


namespace hd::ga {

////////////////////////////////////////////////////////////////////////////////
// Scalar<T> basic operations
////////////////////////////////////////////////////////////////////////////////

// return the value of the scalar as value_t (for use in scripting)
template <typename T>
    requires(std::floating_point<T>)
inline constexpr value_t to_val(Scalar<T> s)
{
    return value_t(s);
}

////////////////////////////////////////////////////////////////////////////////
// Scalar<T> printing support via iostream
////////////////////////////////////////////////////////////////////////////////
template <typename T>
    requires(std::floating_point<T>)
std::ostream& operator<<(std::ostream& os, Scalar<T> v)
{
    os << "(" << T(v) << ")";
    return os;
}
} // namespace hd::ga