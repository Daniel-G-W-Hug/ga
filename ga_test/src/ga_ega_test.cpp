// Copyright 2024-2025, Daniel Hug. All rights reserved.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"

#include <chrono>
#include <iostream>
#include <tuple>
#include <variant>
#include <vector>

#include "fmt/chrono.h"  // chrono support
#include "fmt/format.h"  // formatting
#include "fmt/ostream.h" // ostream support
#include "fmt/ranges.h"  // support printing of (nested) containers & tuples

// include functions to be tested
#include "ga/ga_ega.hpp"

using namespace hd::ga;      // use ga types, constants, etc.
using namespace hd::ga::ega; // use specific operations of EGA (Euclidean GA)

// Include dimension-specific test files
#include "ga_ega2d_test.cpp"
#include "ga_ega3d_test.cpp"
#include "ga_ega4d_test.cpp"


/////////////////////////////////////////////////////////////////////////////////////////
// Generic Euclideon geometric algebra (valid for ega independent of dimension)
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("EGA")
{

    TEST_CASE("G<2,0,0>: ega2d")
    {

        fmt::println("");

        // multivector basis components
        const std::vector<std::string> mv2d_basis_str = {"1", "e1", "e2", "e12"};

        // multivector k_basis components
        // => strings representing scalar2d, vec2d, and pscalar2d
        const std::vector<std::vector<std::string>> k_basis_str = {
            {"1"}, {"e1", "e2"}, {"e12"}};


        fmt::println("mvec_coeff = {}", mv2d_basis_str);
        fmt::println("");

        fmt::println("k_basis_str = {}", k_basis_str);
        fmt::println("");
        for (size_t k = 0; k < k_basis_str.size(); ++k) {
            fmt::println("k_basis_str[{}] = {}", k, k_basis_str.at(k));
        }
        fmt::println("");

        // Experimental multivector2d variant type (mv2d)
        using mv2d = std::variant<scalar2d, vec2d, pscalar2d, mvec2d_e, mvec2d>;

        // Helper function to get type name for variant
        auto get_type_name = [](const mv2d& v) -> std::string {
            return std::visit(
                [](const auto& arg) -> std::string {
                    using T = std::decay_t<decltype(arg)>;
                    if constexpr (std::is_same_v<T, scalar2d>) return "scalar2d";
                    else if constexpr (std::is_same_v<T, vec2d>) return "vec2d";
                    else if constexpr (std::is_same_v<T, pscalar2d>) return "pscalar2d";
                    else if constexpr (std::is_same_v<T, mvec2d_e>) return "mvec2d_e";
                    else if constexpr (std::is_same_v<T, mvec2d>) return "mvec2d";
                    else return "unknown";
                },
                v);
        };

        // Test creation of different types
        mv2d mv_scalar = scalar2d{2.5};
        mv2d mv_vector = vec2d{1.0, 2.0};
        mv2d mv_pscalar = pscalar2d{3.14};
        mv2d mv_even = mvec2d_e{1.0, 2.0};         // scalar + bivector
        mv2d mv_full = mvec2d{1.0, 2.0, 3.0, 4.0}; // scalar + vector + bivector

        mv2d mv_even2 = scalar2d{7.0} + pscalar2d{12.0};

        fmt::println("mv2d variant type tests:");
        fmt::println("mv_scalar: {} (type: {})", mv_scalar.index(),
                     get_type_name(mv_scalar));
        fmt::println("mv_vector: {} (type: {})", mv_vector.index(),
                     get_type_name(mv_vector));
        fmt::println("mv_pscalar: {} (type: {})", mv_pscalar.index(),
                     get_type_name(mv_pscalar));
        fmt::println("mv_even: {} (type: {})", mv_even.index(), get_type_name(mv_even));
        fmt::println("mv_full: {} (type: {})", mv_full.index(), get_type_name(mv_full));
        fmt::println("");

        // Test visitor pattern for operations
        auto print_value = [](const mv2d& v) {
            std::visit([](const auto& arg) { fmt::println("Value: {}", arg); }, v);
        };

        fmt::println("Values:");
        print_value(mv_scalar);
        print_value(mv_vector);
        print_value(mv_pscalar);
        print_value(mv_even);
        print_value(mv_full);
        fmt::println("");

        print_value(mv_even2);
        fmt::println("");

        // Test type checking with std::holds_alternative
        CHECK(std::holds_alternative<scalar2d>(mv_scalar));
        CHECK(std::holds_alternative<vec2d>(mv_vector));
        CHECK(std::holds_alternative<pscalar2d>(mv_pscalar));
        CHECK(std::holds_alternative<mvec2d_e>(mv_even));
        CHECK(std::holds_alternative<mvec2d>(mv_full));

        // Test std::get access
        CHECK(std::get<scalar2d>(mv_scalar) == doctest::Approx(2.5).epsilon(eps));
        CHECK(std::get<vec2d>(mv_vector).x == doctest::Approx(1.0).epsilon(eps));
        CHECK(std::get<vec2d>(mv_vector).y == doctest::Approx(2.0).epsilon(eps));
        CHECK(std::get<pscalar2d>(mv_pscalar) == doctest::Approx(3.14).epsilon(eps));

        // Test mvec2d_e and mvec2d component access
        CHECK(std::get<mvec2d_e>(mv_even).c0 ==
              doctest::Approx(1.0).epsilon(eps)); // scalar part
        CHECK(std::get<mvec2d_e>(mv_even).c1 ==
              doctest::Approx(2.0).epsilon(eps)); // bivector part
        CHECK(std::get<mvec2d>(mv_full).c0 ==
              doctest::Approx(1.0).epsilon(eps)); // scalar
        CHECK(std::get<mvec2d>(mv_full).c1 == doctest::Approx(2.0).epsilon(eps)); // e1
        CHECK(std::get<mvec2d>(mv_full).c2 == doctest::Approx(3.0).epsilon(eps)); // e2
        CHECK(std::get<mvec2d>(mv_full).c3 == doctest::Approx(4.0).epsilon(eps)); // e12

        fmt::println("mv2d variant type implementation complete");
        fmt::println("");
    }
}
