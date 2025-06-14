// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prdxpr_generator.hpp"
#include "ga_prdxpr_config.hpp"
#include "ga_prdxpr_ega2d.hpp"
#include "ga_prdxpr_ega3d.hpp"
#include "ga_prdxpr_pga2dp.hpp"
#include "ga_prdxpr_pga3dp.hpp"

// New simplified main function using the configuration-driven system
// This replaces the repetitive function calls in the original main

int main()
{
    try {
        // Generate specific algebras individually - now including sandwich products
        generate_ega2d_expressions();  // EGA2D complete ✅ (will include sandwich at end)
        generate_ega3d_expressions();  // EGA3D complete ✅ (will include sandwich at end)
        generate_pga2dp_expressions();  // PGA2DP complete ✅ (will include sandwich at end) 
        generate_pga3dp_expressions();  // PGA3DP complete ✅ (will include sandwich at end)
        
    } catch (std::exception const& e) {
        fmt::println("Exception: {}", e.what());
        return 1;
    }

    return 0;
}