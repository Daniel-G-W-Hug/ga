#pragma once

// author: Daniel Hug, 2024

// strong type inspired by Jonathan Boccara's blog:
// https://www.fluentcpp.com/2016/12/08/strong-types-for-strong-interfaces/
//
// and by the blog of foonathan
// https://www.foonathan.net/2016/10/strong-typedefs/


#include <concepts> // std::floating_point

namespace hd::ga {

template <typename T, typename Tag>
    requires(std::floating_point<T>)
class Strong_t {

  public:

    Strong_t() = default; // default initialization
    explicit Strong_t(T const& val) : value(val) {}
    explicit Strong_t(T&& val) noexcept(std::is_nothrow_move_constructible<T>::value) :
        value(std::move(val))
    {
    }

    // these operators bring implicit conversion to the underlying type T
    // (hint: can be suppressed with marking the operator explicit)
    operator T&() noexcept { return value; }
    operator T const&() const noexcept { return value; }

    friend void swap(Strong_t& lhs, Strong_t& rhs) noexcept
    {
        using std::swap;
        swap(static_cast<T&>(lhs), static_cast<T&>(rhs));
    }

  private:

    T value{};
};

} // namespace hd::ga