// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "algebras/ga_prdxpr_cga3dc_config.hpp"
#include "algebras/ga_prdxpr_cga3dc.hpp"

////////////////////////////////////////////////////////////////////////////////
// cga3dc algebra configuration
//
// Basis vectors: e1, e2, e3 Euclidean; e4 projective (w); e5 round (u).
// The null pair satisfies e4*e4 = e5*e5 = 0 and e4.e5 = -1, i.e. the metric
// has an off-diagonal coupling — carried by metric_matrix (metric_signature
// holds its diagonal, as required by the validation).
//
// Diagonalization (exact, dyadic coefficients): internal orthogonal basis
// (e1, e2, e3, e-, e+) with signature (+1, +1, +1, -1, +1) and
//
//     e4 = 0.5*(e- - e+)          e- =  e4 + 0.5*e5
//     e5 = e- + e+                e+ = -e4 + 0.5*e5
//
// The normalized point embedding q = x*e1 + y*e2 + z*e3 + e4
// + 0.5*(x²+y²+z²)*e5 is null: q ⟑ q = 0.
////////////////////////////////////////////////////////////////////////////////

AlgebraConfig get_cga3dc_algebra_config()
{
    // Extract basis prefix from vector basis and validate consistency
    std::string const prefix = extract_basis_prefix(mv3dc_basis_kvec[1]);
    validate_basis_consistency(mv3dc_basis, mv3dc_basis_kvec, prefix, one_str());

    AlgebraConfig config;
    config.basis_vectors = {"e1", "e2", "e3", "e4", "e5"};
    config.metric_signature = {1, 1, 1, 0, 0}; // diagonal of metric_matrix
    config.multivector_basis = mv3dc_basis;
    config.scalar_name = one_str();
    config.basis_prefix = prefix;

    config.metric_matrix = {{1, 0, 0, 0, 0},
                            {0, 1, 0, 0, 0},
                            {0, 0, 1, 0, 0},
                            {0, 0, 0, 0, -1},
                            {0, 0, 0, -1, 0}};
    config.basis_change = {{1.0, 0.0, 0.0, 0.0, 0.0},
                           {0.0, 1.0, 0.0, 0.0, 0.0},
                           {0.0, 0.0, 1.0, 0.0, 0.0},
                           {0.0, 0.0, 0.0, 0.5, -0.5},
                           {0.0, 0.0, 0.0, 1.0, 1.0}};
    config.basis_change_inv = {{1.0, 0.0, 0.0, 0.0, 0.0},
                               {0.0, 1.0, 0.0, 0.0, 0.0},
                               {0.0, 0.0, 1.0, 0.0, 0.0},
                               {0.0, 0.0, 0.0, 1.0, 0.5},
                               {0.0, 0.0, 0.0, -1.0, 0.5}};
    config.diag_signature = {1, 1, 1, -1, 1};

    return config;
}
