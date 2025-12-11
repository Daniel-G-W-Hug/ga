#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include <concepts>
#include <type_traits>

/////////////////////////////////////////////////////////////////////////////////////////
// type tags for consistent type definitions
/////////////////////////////////////////////////////////////////////////////////////////

namespace hd::ga {

struct default_tag {}; // used as reference for inheriting ctors from base types

struct ega2d_tag {};
struct scalar2d_tag : public ega2d_tag {};
struct vec2d_tag : public ega2d_tag {};
struct pscalar2d_tag : public ega2d_tag {};
struct mvec2d_e_tag : public ega2d_tag {};
struct mvec2d_tag : public ega2d_tag {};

struct ega3d_tag {};
struct scalar3d_tag : public ega3d_tag {};
struct vec3d_tag : public ega3d_tag {};
struct bivec3d_tag : public ega3d_tag {};
struct pscalar3d_tag : public ega3d_tag {};
struct mvec3d_e_tag : public ega3d_tag {};
struct mvec3d_u_tag : public ega3d_tag {};
struct mvec3d_tag : public ega3d_tag {};

struct ega4d_tag {};
struct scalar4d_tag : public ega4d_tag {};
struct vec4d_tag : public ega4d_tag {};
struct bivec4d_tag : public ega4d_tag {};
struct trivec4d_tag : public ega4d_tag {};
struct pscalar4d_tag : public ega4d_tag {};
struct mvec4d_e_tag : public ega4d_tag {};
struct mvec4d_u_tag : public ega4d_tag {};
struct mvec4d_tag : public ega4d_tag {};

struct pga2dp_tag {};
struct scalar2dp_tag : public pga2dp_tag {};
struct vec2dp_tag : public pga2dp_tag {};
struct bivec2dp_tag : public pga2dp_tag {};
struct pscalar2dp_tag : public pga2dp_tag {};
struct mvec2dp_e_tag : public pga2dp_tag {};
struct mvec2dp_u_tag : public pga2dp_tag {};
struct mvec2dp_tag : public pga2dp_tag {};
struct dual_number2dp_tag : public pga2dp_tag {};

struct pga3dp_tag {};
struct scalar3dp_tag : public pga3dp_tag {};
struct vec3dp_tag : public pga3dp_tag {};
struct bivec3dp_tag : public pga3dp_tag {};
struct trivec3dp_tag : public pga3dp_tag {};
struct pscalar3dp_tag : public pga3dp_tag {};
struct mvec3dp_e_tag : public pga3dp_tag {};
struct mvec3dp_u_tag : public pga3dp_tag {};
struct mvec3dp_tag : public pga3dp_tag {};
struct dual_number3dp_tag : public pga3dp_tag {};


/////////////////////////////////////////////////////////////////////////////////////////
// concept to check for consistent type definitions checks via template specialization
//
// USAGE: to check whether two tags come from the same base class, use the concept as:
//
// template<typename T, typename U>
//     requires(same_base_class<T, U>)
// ... your code added here
//
// HINT: add any new types here or the code will not compile, if same_base_class<T,U>
//       is used anywhere for unknown added type tags
/////////////////////////////////////////////////////////////////////////////////////////

// Primary template - will trigger static_assert for unknown types
template <typename T> struct base_class {
    static_assert(std::is_same_v<T, void>,
                  "Type doesn't belong to any known base class family");
};

// Specialization for each base class family
template <typename T>
    requires std::is_base_of_v<ega2d_tag, T>
struct base_class<T> {
    using type = ega2d_tag;
};

template <typename T>
    requires std::is_base_of_v<ega3d_tag, T>
struct base_class<T> {
    using type = ega3d_tag;
};

template <typename T>
    requires std::is_base_of_v<ega4d_tag, T>
struct base_class<T> {
    using type = ega4d_tag;
};

template <typename T>
    requires std::is_base_of_v<pga2dp_tag, T>
struct base_class<T> {
    using type = pga2dp_tag;
};

template <typename T>
    requires std::is_base_of_v<pga3dp_tag, T>
struct base_class<T> {
    using type = pga3dp_tag;
};

// hint: just add newly defined tag classes here

// Helper alias template
template <typename T> using base_class_t = typename base_class<T>::type;

// Concept to check if two types have the same base class
template <typename T, typename U>
concept same_base_class = std::same_as<base_class_t<T>, base_class_t<U>>;


/////////////////////////////////////////////////////////////////////////////////////////
// concept to check for numerically permissible types
// (with more precision/clarity for template type checking)
//
// USAGE:
//
// template<typename T>
//     requires(numeric_type<T>)
// ... your code added here
//
// Accepts:
//        Floating-point: float, double, long double
//        Integer types: short, int, long, long long, unsigned short, unsigned int,
//        unsigned long, unsigned long long
//
// Rejects:
//        bool, char, signed char, unsigned char, char8_t, char16_t, char32_t, wchar_t
//
/////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
concept numeric_type = std::floating_point<T> ||
                       (std::integral<T> && !std::same_as<std::remove_cv_t<T>, bool> &&
                        !std::same_as<std::remove_cv_t<T>, char> &&
                        !std::same_as<std::remove_cv_t<T>, signed char> &&
                        !std::same_as<std::remove_cv_t<T>, unsigned char> &&
                        !std::same_as<std::remove_cv_t<T>, char8_t> &&
                        !std::same_as<std::remove_cv_t<T>, char16_t> &&
                        !std::same_as<std::remove_cv_t<T>, char32_t> &&
                        !std::same_as<std::remove_cv_t<T>, wchar_t>);

} // namespace hd::ga