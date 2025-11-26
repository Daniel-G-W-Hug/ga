#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <cmath>   // std::cos, std::sin
#include <numbers> // math constants like pi

#include "detail/type_t/ga_scalar_t.hpp"
#include "ga_value_t.hpp"

namespace hd::ga {

using std::numbers::pi; // make pi available for users

// conversion between degrees and radians [360° = 2*pi rad]
constexpr value_t deg2rad(value_t value) { return value / 180.0 * pi; }
constexpr value_t rad2deg(value_t value) { return value * 180.0 / pi; }

// conversion between revolutions per minute [rpm] and radians per second [rad/s = radps]
constexpr value_t rpm2radps(value_t value) { return value * 2.0 * pi / 60.0; };
constexpr value_t radps2rpm(value_t value) { return value * 60.0 / (2.0 * pi); };

// conversion between revolutions per second (rps) = Hertz [1 rps = 1 Hz = 1 1/s^(-1)] and
// radians per second [rad/s = radps]
constexpr value_t Hz2radps(value_t value) { return value * 2.0 * pi; };
constexpr value_t radps2Hz(value_t value) { return value / (2.0 * pi); };

// sign function for floating point types
// returns +1.0 for value >= 0.0 and -1.0 for value < 0.0
// HINT: this is intentionally NOT the typical signum function
//       which returns 0.0 for value == 0.0
template <typename T>
    requires(std::floating_point<T>)
constexpr T sign(T value)
{
    if (value >= 0.0) {
        return 1.0;
    }
    return -1.0;
}

// sign function overload for Scalar_t types
template <typename T, typename Tag>
    requires(std::floating_point<T>)
constexpr T sign(Scalar_t<T, Tag> s)
{
    return sign(T(s));
}

// Templates for is_even and is_odd work with any integer type
template <typename T>
    requires(std::is_integral_v<T>) // is_even: argument must be an integral type
bool is_even(T arg)
{
    return arg % 2 == 0;
}

template <typename T>
    requires(std::is_integral_v<T>) // is_odd: argument must be an integral type
bool is_odd(T arg)
{
    return arg % 2 != 0;
}

template <typename T>
    requires(std::floating_point<T>)
class CosSin {

  public:

    // ctor
    CosSin() = delete; // always require an argument for initialization
    CosSin(T const& x) : m_c(std::cos(x)), m_s(std::sin(x)) {}

    // simple gettern (no setter, this must be done via the ctor)
    T const& c() const { return m_c; }
    T const& s() const { return m_s; }

  private:

    T m_c;
    T m_s;
};

} // namespace hd::ga