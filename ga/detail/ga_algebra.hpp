#pragma once

// author: Daniel Hug, 2024

#include <array>
#include <cstdint> // uint8_t
#include <string>

namespace hd::ga {

// create an algebra G(P,N,Z) with:
//
// P generators for numbers that square to +1
// N generators for numbers that square to -1
// Z generators for numbers that square to  0
//
// used to provide frequently used values of the algebra
// should be assigned to a const variable
template <uint8_t P, uint8_t N = 0, uint8_t Z = 0>
    requires(P + N + Z >= 2 && P + N + Z <= 4 && P >= 2 && P <= 4 && N == 0 &&
             Z <= 1) // no implementation for other algebras yet
struct algebra {
    constexpr static uint8_t p() { return P; }; // number of +1 generators
    constexpr static uint8_t n() { return N; }; // number of -1 generators
    constexpr static uint8_t z() { return Z; }; // number of  0 generators

    constexpr static uint8_t dim_space() { return P + N + Z; }; // dimension of the space

    constexpr static uint8_t num_components()
    {
        return 1 << (dim_space()); // the number of basis components == 2^dim_space
    };

    constexpr static std::array<const uint8_t, dim_space() + 1> num_components_grade =
        []() -> std::array<const uint8_t, dim_space() + 1> {
        if constexpr (dim_space() == 2) return {1, 2, 1};
        if constexpr (dim_space() == 3) return {1, 3, 3, 1};
        if constexpr (dim_space() == 4) return {1, 4, 6, 4, 1};
    }();

    constexpr static std::array<const char[6], num_components()> basis_name =
        []() -> std::array<const char[6], num_components()> {
        //
        // this defines the basis components of the multivectors:
        //
        // ega:
        if constexpr (dim_space() == 2 && (p() == 2 && n() == 0 && z() == 0))
            return {"    1", "   e1", "   e2", "  e12"};
        if constexpr (dim_space() == 3 && (p() == 3 && n() == 0 && z() == 0))
            return {"    1", "   e1", "   e2", "   e3",
                    "  e23", "  e31", "  e12", " e123"};
        //
        // pga:
        if constexpr (dim_space() == 3 && (p() == 2 && n() == 0 && z() == 1))
            return {"    1", "   e1", "   e2", "   e3",
                    "  e23", "  e31", "  e12", " e321"};
        if constexpr (dim_space() == 4 && (p() == 3 && n() == 0 && z() == 1))
            return {"    1", "   e1", "   e2", "   e3", "   e4", "  e41",
                    "  e42", "  e43", "  e23", "  e31", "  e12", " e423",
                    " e431", " e412", " e321", "e1234"};
    }();
};

} // namespace hd::ga