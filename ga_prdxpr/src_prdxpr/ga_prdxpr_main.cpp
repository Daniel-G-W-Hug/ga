// Copyright 2024-2025, Daniel Hug. All rights reserved.
// Complete configurable GA product expression generator for all four algebras

#include "ga_prdxpr_ega2d_config.hpp"
#include "ga_prdxpr_ega3d_config.hpp"
#include "ga_prdxpr_generator.hpp"
#include "ga_prdxpr_pga2dp_config.hpp"
#include "ga_prdxpr_pga3dp_config.hpp"
#include <fmt/core.h>
#include <vector>

using namespace configurable;

// Helper function to generate all products for an algebra with proper separators
void generate_algebra_products(ConfigurableGenerator& generator,
                               const std::vector<ProductConfig>& configs,
                               const AlgebraData& algebra_data, bool is_last_algebra)
{
    for (size_t i = 0; i < configs.size(); ++i) {
        generator.generate_product_expressions(algebra_data, configs[i]);

        // Add separator after each product (except the last product of the last algebra)
        if (!(is_last_algebra && i == configs.size() - 1)) {
            fmt::println(
                "-------------------------------------------------------------------");
            fmt::println("");
        }
    }
}

int main()
{
    try {
        ConfigurableGenerator generator;

        // EGA2D - All 10 products
        {
            auto ega2d_algebra = create_ega2d_algebra_data();
            std::vector<ProductConfig> ega2d_configs = {
                get_ega2d_gpr_config(),         get_ega2d_cmt_config(),
                get_ega2d_wdg_config(),         get_ega2d_dot_config(),
                get_ega2d_dot_alt_config(),     get_ega2d_rwdg_config(),
                get_ega2d_lcontract_config(),   get_ega2d_rcontract_config(),
                get_ega2d_lexpand_config(),     get_ega2d_rexpand_config(),
                get_ega2d_sandwich_gpr_config()};
            generate_algebra_products(generator, ega2d_configs, ega2d_algebra, false);
        }

        // EGA3D - All 10 products
        {
            auto ega3d_algebra = create_ega3d_algebra_data();
            std::vector<ProductConfig> ega3d_configs = {
                get_ega3d_gpr_config(),       get_ega3d_cmt_config(),
                get_ega3d_wdg_config(),       get_ega3d_dot_config(),
                get_ega3d_rwdg_config(),      get_ega3d_lcontract_config(),
                get_ega3d_rcontract_config(), get_ega3d_lexpand_config(),
                get_ega3d_rexpand_config(),   get_ega3d_sandwich_gpr_config()};
            generate_algebra_products(generator, ega3d_configs, ega3d_algebra, false);
        }

        // PGA2DP - All 16 products
        {
            auto pga2dp_algebra = create_pga2dp_algebra_data();
            std::vector<ProductConfig> pga2dp_configs = {
                get_pga2dp_gpr_config(),
                get_pga2dp_cmt_config(),
                get_pga2dp_wdg_config(),
                get_pga2dp_dot_config(),
                get_pga2dp_rwdg_config(),
                get_pga2dp_rdot_config(),
                get_pga2dp_rgpr_config(),
                get_pga2dp_right_bulk_contract_config(),
                get_pga2dp_right_weight_contract_config(),
                get_pga2dp_left_bulk_contract_config(),
                get_pga2dp_left_weight_contract_config(),
                get_pga2dp_right_bulk_expand_config(),
                get_pga2dp_right_weight_expand_config(),
                get_pga2dp_left_bulk_expand_config(),
                get_pga2dp_left_weight_expand_config(),
                get_pga2dp_sandwich_rgpr_config()};
            generate_algebra_products(generator, pga2dp_configs, pga2dp_algebra, false);
        }

        // PGA3DP - All 16 products
        {
            auto pga3dp_algebra = create_pga3dp_algebra_data();
            std::vector<ProductConfig> pga3dp_configs = {
                get_pga3dp_gpr_config(),
                get_pga3dp_cmt_config(),
                get_pga3dp_wdg_config(),
                get_pga3dp_dot_config(),
                get_pga3dp_rwdg_config(),
                get_pga3dp_rdot_config(),
                get_pga3dp_rgpr_config(),
                get_pga3dp_right_bulk_contract_config(),
                get_pga3dp_right_weight_contract_config(),
                get_pga3dp_left_bulk_contract_config(),
                get_pga3dp_left_weight_contract_config(),
                get_pga3dp_right_bulk_expand_config(),
                get_pga3dp_right_weight_expand_config(),
                get_pga3dp_left_bulk_expand_config(),
                get_pga3dp_left_weight_expand_config(),
                get_pga3dp_sandwich_rgpr_config()};
            generate_algebra_products(generator, pga3dp_configs, pga3dp_algebra, true);
        }

        fmt::println(
            "-------------------------------------------------------------------");
        fmt::println("");

        return 0;
    }
    catch (const std::exception& e) {
        fmt::println("Error: {}", e.what());
        return 1;
    }
}