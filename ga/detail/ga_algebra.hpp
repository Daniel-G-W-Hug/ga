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
             Z == 0) // no implementation for other algebras yet
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
        if constexpr (dim_space() == 2) return {"    1", "   e1", "   e2", "  e12"};
        if constexpr (dim_space() == 3)
            return {"    1", "   e1", "   e2", "   e3",
                    "  e23", "  e31", "  e12", " e123"};
        if constexpr (dim_space() == 4)
            return {"    1", "   e0", "   e1", "   e2", "   e3", "  e01",
                    "  e02", "  e03", "  e23", "  e31", "  e12", " e023",
                    " e031", " e012", " e123", "e0123"};
    }();
};

} // namespace hd::ga