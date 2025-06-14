#pragma once

// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_config.hpp"
#include "ga_prdxpr_common.hpp"

////////////////////////////////////////////////////////////////////////////////
// ProductExpressionGenerator - unified generator for all algebras
////////////////////////////////////////////////////////////////////////////////

class ProductExpressionGenerator {
public:
    // Generate all products for all algebras
    void generate_all_algebras();
    
    // Generate all products for a specific algebra
    void generate_algebra(const AlgebraConfig& config, 
                         const std::vector<ProductDefinition>& products);
    
    // Generate a specific product type for an algebra
    void generate_product_type(const AlgebraConfig& config,
                              const ProductDefinition& product);

private:
    // Generate a single product case (mirrors one current function call)
    void generate_product_case(const AlgebraConfig& config,
                              const ProductDefinition& product,
                              const ProductCase& case_def,
                              const prd_table& basis_tab);
    
    // Helper functions for different dimensions
    void generate_case_2d(const AlgebraConfig& config,
                         const ProductDefinition& product,
                         const ProductCase& case_def,
                         const prd_table& basis_tab);
                         
    void generate_case_3d(const AlgebraConfig& config,
                         const ProductDefinition& product,
                         const ProductCase& case_def,
                         const prd_table& basis_tab);
                         
    void generate_case_4d(const AlgebraConfig& config,
                         const ProductDefinition& product,
                         const ProductCase& case_def,
                         const prd_table& basis_tab);
    
    // Get basis product table for a specific product type
    prd_table get_basis_product_table(const AlgebraConfig& config,
                                     const ProductDefinition& product);
    
    // Generate sandwich product calculations (two-step process)
    void generate_sandwich_products(const AlgebraConfig& config,
                                   const ProductDefinition& product,
                                   const prd_table& basis_tab);
};

////////////////////////////////////////////////////////////////////////////////
// Convenience functions for generating specific algebras
////////////////////////////////////////////////////////////////////////////////

void generate_ega2d_expressions();
void generate_ega3d_expressions();
void generate_pga2dp_expressions();
void generate_pga3dp_expressions();