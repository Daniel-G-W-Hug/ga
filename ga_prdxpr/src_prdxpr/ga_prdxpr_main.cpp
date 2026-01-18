// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.
// Complete configurable GA product expression generator for all four algebras

#include "ga_prdxpr_ega2d_config.hpp"
#include "ga_prdxpr_ega3d_config.hpp"
#include "ga_prdxpr_generator.hpp"
#include "ga_prdxpr_metric_export.hpp"
#include "ga_prdxpr_options.hpp"
#include "ga_prdxpr_pga2dp_config.hpp"
#include "ga_prdxpr_pga3dp_config.hpp"
#include "ga_prdxpr_rule_generator.hpp"
#include "ga_prdxpr_sta4d_config.hpp"
#include <fmt/core.h>
#include <vector>

using namespace configurable;

// Helper function to generate all products for an algebra with proper separators
void generate_algebra_products(ConfigurableGenerator& generator,
                               std::vector<ProductConfig> const& configs,
                               AlgebraData const& algebra_data, bool is_last_algebra,
                               GeneratorOptions const& options)
{
    bool first_product = true;

    for (size_t i = 0; i < configs.size(); ++i) {
        ProductConfig const& config = configs[i];

        // Filter products if requested
        if (!options.should_generate_product(config.product_name)) {
            continue;
        }

        // Add separator before product (except first)
        if (!first_product) {
            fmt::println(
                "-------------------------------------------------------------------");
            fmt::println("");
        }
        first_product = false;

        generator.generate_product_expressions(algebra_data, config, options);
    }

    // Add separator after last product if not the last algebra
    if (!first_product && !is_last_algebra) {
        fmt::println(
            "-------------------------------------------------------------------");
        fmt::println("");
    }
}

int main(int argc, char const* argv[])
{
    try {
        // Parse command-line arguments
        ArgumentParser parser(argc, argv);

        // Handle special flags
        if (parser.should_show_help()) {
            parser.print_help();
            return 0;
        }

        if (parser.should_show_list()) {
            parser.print_list();
            return 0;
        }

        if (parser.has_error()) {
            fmt::println("Error: {}", parser.get_error_message());
            fmt::println("");
            parser.print_help();
            return 1;
        }

        GeneratorOptions const& options = parser.get_options();
        ConfigurableGenerator generator;

        // EGA2D
        if (options.should_generate_algebra("ega2d") &&
            (options.should_show_coeffs() || options.should_show_tables())) {
            auto ega2d_algebra = create_ega2d_algebra_data();
            std::vector<ProductConfig> ega2d_configs = {
                get_ega2d_gpr_config(),         get_ega2d_gpr_alt_config(),
                get_ega2d_twdg1_config(),       get_ega2d_cmt_config(),
                get_ega2d_wdg_config(),         get_ega2d_dot_config(),
                get_ega2d_dot_alt_config(),

                get_ega2d_lcontract_config(),   get_ega2d_rcontract_config(),

                get_ega2d_lexpand_config(),     get_ega2d_rexpand_config(),

                get_ega2d_rgpr_config(),        get_ega2d_rgpr_alt_config(),
                get_ega2d_rtwdg1_config(),      get_ega2d_rcmt_config(),
                get_ega2d_rwdg_config(),        get_ega2d_rdot_config(),

                get_ega2d_sandwich_gpr_config()};

            generate_algebra_products(generator, ega2d_configs, ega2d_algebra, false,
                                      options);
        }

        // EGA3D
        if (options.should_generate_algebra("ega3d") &&
            (options.should_show_coeffs() || options.should_show_tables())) {
            auto ega3d_algebra = create_ega3d_algebra_data();
            std::vector<ProductConfig> ega3d_configs = {
                get_ega3d_gpr_config(),         get_ega3d_gpr_alt_config(),
                get_ega3d_twdg1_config(),       get_ega3d_cmt_config(),
                get_ega3d_wdg_config(),         get_ega3d_dot_config(),

                get_ega3d_lcontract_config(),   get_ega3d_rcontract_config(),

                get_ega3d_lexpand_config(),     get_ega3d_rexpand_config(),

                get_ega3d_rgpr_config(),        get_ega3d_rgpr_alt_config(),
                get_ega3d_rtwdg1_config(),      get_ega3d_rcmt_config(),
                get_ega3d_rwdg_config(),        get_ega3d_rdot_config(),

                get_ega3d_sandwich_gpr_config()};

            generate_algebra_products(generator, ega3d_configs, ega3d_algebra, false,
                                      options);
        }

        // PGA2DP
        if (options.should_generate_algebra("pga2dp") &&
            (options.should_show_coeffs() || options.should_show_tables())) {
            auto pga2dp_algebra = create_pga2dp_algebra_data();
            std::vector<ProductConfig> pga2dp_configs = {
                get_pga2dp_gpr_config(),
                get_pga2dp_gpr_alt_config(),
                get_pga2dp_twdg1_config(),
                get_pga2dp_cmt_config(),
                get_pga2dp_wdg_config(),
                get_pga2dp_dot_config(),

                get_pga2dp_left_bulk_contract_config(),
                get_pga2dp_right_bulk_contract_config(),

                get_pga2dp_left_weight_contract_config(),
                get_pga2dp_right_weight_contract_config(),

                get_pga2dp_left_bulk_expand_config(),
                get_pga2dp_right_bulk_expand_config(),

                get_pga2dp_left_weight_expand_config(),
                get_pga2dp_right_weight_expand_config(),

                get_pga2dp_rgpr_config(),
                get_pga2dp_rgpr_alt_config(),
                get_pga2dp_rtwdg1_config(),
                get_pga2dp_rcmt_config(),
                get_pga2dp_rwdg_config(),
                get_pga2dp_rdot_config(),

                get_pga2dp_sandwich_rgpr_config()};

            generate_algebra_products(generator, pga2dp_configs, pga2dp_algebra, false,
                                      options);
        }

        // PGA3DP
        if (options.should_generate_algebra("pga3dp") &&
            (options.should_show_coeffs() || options.should_show_tables())) {
            auto pga3dp_algebra = create_pga3dp_algebra_data();
            std::vector<ProductConfig> pga3dp_configs = {
                get_pga3dp_gpr_config(),
                get_pga3dp_gpr_alt_config(),
                get_pga3dp_twdg1_config(),
                get_pga3dp_cmt_config(),
                get_pga3dp_wdg_config(),
                get_pga3dp_dot_config(),

                get_pga3dp_left_bulk_contract_config(),
                get_pga3dp_right_bulk_contract_config(),

                get_pga3dp_left_weight_contract_config(),
                get_pga3dp_right_weight_contract_config(),

                get_pga3dp_left_bulk_expand_config(),
                get_pga3dp_right_bulk_expand_config(),

                get_pga3dp_left_weight_expand_config(),
                get_pga3dp_right_weight_expand_config(),

                get_pga3dp_rgpr_config(),
                get_pga3dp_rgpr_alt_config(),
                get_pga3dp_rtwdg1_config(),
                get_pga3dp_rcmt_config(),
                get_pga3dp_rwdg_config(),
                get_pga3dp_rdot_config(),

                get_pga3dp_sandwich_rgpr_config()};

            generate_algebra_products(generator, pga3dp_configs, pga3dp_algebra, true,
                                      options);
        }

        // STA4D
        if (options.should_generate_algebra("sta4d") &&
            (options.should_show_coeffs() || options.should_show_tables())) {
            auto sta4d_algebra = create_sta4d_algebra_data();
            std::vector<ProductConfig> sta4d_configs = {
                get_sta4d_gpr_config(),           get_sta4d_cmt_config(),
                get_sta4d_wdg_config(),           get_sta4d_dot_config(),

                get_sta4d_left_contract_config(), get_sta4d_right_contract_config(),

                get_sta4d_left_expand_config(),   get_sta4d_right_expand_config(),

                get_sta4d_rgpr_config(),          get_sta4d_rcmt_config(),
                get_sta4d_rwdg_config(),          get_sta4d_rdot_config(),

                get_sta4d_sandwich_gpr_config()};

            generate_algebra_products(generator, sta4d_configs, sta4d_algebra, true,
                                      options);
        }

        // Generate extended metric exports for all algebras (if requested)
        if (options.should_show_metrics()) {

            // Export metrics for each algebra (filtered by algebra selection)
            if (options.should_generate_algebra("ega2d")) {
                AlgebraConfig ega2d_config = get_ega2d_algebra_config();
                print_metrics_for_algebra(ega2d_config);
            }

            if (options.should_generate_algebra("ega3d")) {
                AlgebraConfig ega3d_config = get_ega3d_algebra_config();
                print_metrics_for_algebra(ega3d_config);
            }

            if (options.should_generate_algebra("pga2dp")) {
                AlgebraConfig pga2dp_config = get_pga2dp_algebra_config();
                print_metrics_for_algebra(pga2dp_config);
            }

            if (options.should_generate_algebra("pga3dp")) {
                AlgebraConfig pga3dp_config = get_pga3dp_algebra_config();
                print_metrics_for_algebra(pga3dp_config);
            }

            if (options.should_generate_algebra("sta4d")) {
                AlgebraConfig sta4d_config = get_sta4d_algebra_config();
                print_metrics_for_algebra(sta4d_config);
            }
        }

        return 0;
    }
    catch (const std::exception& e) {
        fmt::println("Error: {}", e.what());
        return 1;
    }
}