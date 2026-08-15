#pragma once

// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include <array>   // std::array
#include <cstdint> // uint8_t

namespace hd::ga {

/////////////////////////////////////////////////////////////////////////////////////////
// create an algebra agebra<P,N,Z> = G(P,N,Z) with:
//      P generators for numbers that square to +1
//      N generators for numbers that square to -1
//      Z generators for numbers that square to  0
//
// provides infomation on the dimension and basis vectors of the used algebra
/////////////////////////////////////////////////////////////////////////////////////////

template <uint8_t P, uint8_t N = 0, uint8_t Z = 0>
    requires(P + N + Z >= 2) && (P + N + Z <= 5) && (P >= 1) && (P <= 4) && (N <= 3) &&
            (Z <= 1) // no implementation for other algebras yet
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
        if constexpr (dim_space() == 2) {
            return {1, 2, 1};
        }
        if constexpr (dim_space() == 3) {
            return {1, 3, 3, 1};
        }
        if constexpr (dim_space() == 4) {
            return {1, 4, 6, 4, 1};
        }
        if constexpr (dim_space() == 5) {
            return {1, 5, 10, 10, 5, 1};
        }
    }();

    constexpr static std::array<const char[7], num_components()> basis_name =
        []() -> std::array<const char[7], num_components()> {
        //
        // this defines the basis components of the multivectors:
        //
        // ega: algebra<2,0,0>
        if constexpr (dim_space() == 2 && (p() == 2 && n() == 0 && z() == 0)) {
            return {"    1", "   e1", "   e2", "  e12"};
        }
        //
        // ega: algebra<3,0,0>
        if constexpr (dim_space() == 3 && (p() == 3 && n() == 0 && z() == 0)) {
            return {"    1", "   e1", "   e2", "   e3",
                    "  e23", "  e31", "  e12", " e123"};
        }
        //
        // pga2dp:  algebra<2,0,1>
        if constexpr (dim_space() == 3 && (p() == 2 && n() == 0 && z() == 1)) {
            return {"    1", "   e1", "   e2", "   e3",
                    "  e31", "  e32", "  e12", " e321"};
        }
        //
        // pga3dp:  algebra<3,0,1>
        if constexpr (dim_space() == 4 && (p() == 3 && n() == 0 && z() == 1)) {
            return {"    1", "   e1", "   e2", "   e3", "   e4", "  e41",
                    "  e42", "  e43", "  e23", "  e31", "  e12", " e423",
                    " e431", " e412", " e321", "e1234"};
        }
        //
        // sta4ds:  algebra<1,3,0> "mostly negative"-convention
        if constexpr (dim_space() == 4 && (p() == 1 && n() == 3 && z() == 0)) {
            return {"   1",  "   g1", "   g2", "   g3", "   g4", "  g14",
                    "  g24", "  g34", "  g23", "  g31", "  g12", " g234",
                    " g314", " g124", " g123", "g1234"};
        }
        //
        // cga2dc:  algebra<3,1,0> (conformal 2d) — STORED in the null basis
        // {e1, e2, e3, e4} with e3*e3 = e4*e4 = 0 and e3.e4 = -1; the signature
        // G(3,1,0) refers to the diagonalized form of that metric
        if constexpr (dim_space() == 4 && (p() == 3 && n() == 1 && z() == 0)) {
            return {"    1", "   e1", "   e2", "   e3", "   e4", "  e31",
                    "  e32", "  e12", "  e14", "  e24", "  e34", " e314",
                    " e324", " e124", " e321", "e1234"};
        }
        //
        // cga3dc:  algebra<4,1,0> (conformal 3d) — STORED in the null basis
        // {e1, e2, e3, e4, e5} with e4*e4 = e5*e5 = 0 and e4.e5 = -1; the
        // signature G(4,1,0) refers to the diagonalized form of that metric
        if constexpr (dim_space() == 5 && (p() == 4 && n() == 1 && z() == 0)) {
            return {"     1", "    e1", "    e2", "    e3", "    e4", "    e5", "   e41",
                    "   e42", "   e43", "   e23", "   e31", "   e12", "   e15", "   e25",
                    "   e35", "   e45", "  e415", "  e425", "  e435", "  e235", "  e315",
                    "  e125", "  e423", "  e431", "  e412", "  e321", " e4235", " e4315",
                    " e4125", " e3215", " e1234", "e12345"};
        }
    }();
};

} // namespace hd::ga