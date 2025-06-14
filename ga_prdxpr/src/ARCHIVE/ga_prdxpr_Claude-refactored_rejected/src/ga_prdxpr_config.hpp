#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_common.hpp"
#include <map>
#include <vector>
#include <string>

////////////////////////////////////////////////////////////////////////////////
// Configuration structures for algebra-specific settings
////////////////////////////////////////////////////////////////////////////////

struct ProductCase {
    std::string description;        // "mv * mv -> mv"
    std::string left_coeff_key;     // "A" -> maps to coefficients["A"]  
    std::string right_coeff_key;    // "B" -> maps to coefficients["B"]
    std::string left_filter;        // "mv" -> filter_2d::mv
    std::string right_filter;       // "mv" -> filter_2d::mv  
    bool enabled = true;
    std::string comment = "";       // for future use or documentation
};

struct AlgebraConfig {
    std::string name;
    std::string description;
    mvec_coeff basis;
    std::map<std::string, prd_rules> product_rules;  // "gpr", "wdg", "dot", etc.
    std::map<std::string, mvec_coeff> coefficients;  // "A", "B", "M", "svps", etc.
    std::map<std::string, prd_rules> complement_rules; // "lcmpl", "rcmpl"
    std::vector<std::string> filter_names;           // available filter types
    int dimension;  // 2d, 3d, 4d
};

enum class ProductType {
    // Specific product types that follow original implementation exactly
    geometric_product,   // gpr: geometric product 
    commutator_product,  // cmt: commutator product (asymmetric part of gpr)
    wedge_product,       // wdg: wedge product (outer product)
    inner_product,       // dot: inner product
    
    // Advanced products
    regressive_wedge,    // rwdg(A,B) = lcmpl(wdg(rcmpl(A), rcmpl(B)))
    left_contraction,    // A << B = rwdg(lcmpl(A), B)
    right_contraction,   // A >> B = rwdg(A, rcmpl(B))
    left_expansion,      // lexpand(A,B) = wdg(lcmpl(A), B)
    right_expansion,     // rexpand(A,B) = wdg(A, rcmpl(B))
    regressive_inner,    // rdot: regressive inner product
    regressive_geometric,// rgpr: regressive geometric product
    right_bulk_contraction,   // PGA bulk contractions
    left_bulk_contraction,
    right_weight_contraction, // PGA weight contractions  
    left_weight_contraction,
    right_bulk_expansion,     // PGA bulk expansions
    left_bulk_expansion,
    right_weight_expansion,   // PGA weight expansions
    left_weight_expansion,
    sandwich,            // rotor * object * rev(rotor)
    regressive_sandwich  // motor * object * rrev(motor)
};

struct ProductDefinition {
    std::string key;              // "gpr", "wdg", "dot"
    std::string display_name;     // "geometric product"
    std::string rule_key;         // maps to AlgebraConfig.product_rules (for direct products)
    std::string operator_str;     // "*", "^", "<<", ">>"
    ProductType type = ProductType::geometric_product;
    bool needs_asymmetric = false;  // for commutator products
    bool enabled = true;
    std::vector<ProductCase> cases;  // product-specific cases
};

////////////////////////////////////////////////////////////////////////////////
// Configuration namespace - contains all algebra configurations
////////////////////////////////////////////////////////////////////////////////

namespace configs {
    extern const AlgebraConfig ega2d_config;
    extern const AlgebraConfig ega3d_config;
    extern const AlgebraConfig pga2dp_config;
    extern const AlgebraConfig pga3dp_config;
    
    // Product definitions for each algebra type (now include their own cases)
    extern const std::vector<ProductDefinition> ega2d_products;
    extern const std::vector<ProductDefinition> ega3d_products; 
    extern const std::vector<ProductDefinition> pga2dp_products;
    extern const std::vector<ProductDefinition> pga3dp_products;

    ////////////////////////////////////////////////////////////////////////////////
    // Helper functions for configuration management
    ////////////////////////////////////////////////////////////////////////////////

    // Get filter enum value from string (algebra-specific)
    filter_2d get_filter_2d_from_string(const std::string& filter_name);
    filter_3d get_filter_3d_from_string(const std::string& filter_name);
    filter_4d get_filter_4d_from_string(const std::string& filter_name);

    // Validate configuration completeness
    bool validate_algebra_config(const AlgebraConfig& config);
}