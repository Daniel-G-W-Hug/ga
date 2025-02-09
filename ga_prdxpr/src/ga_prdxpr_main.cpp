// author: Daniel Hug, 2025

#include "ga_prd_ega2d.hpp"
#include "ga_prd_ega3d.hpp"

#include "ga_prd_ega2dp.hpp"

#include "ga_prdxpr_transformer.hpp" // just for ParseError


// multivector basis, product rules, simplification rules and coefficients need to be
// provided in header files included above by the user

int main()
{

    try {

        // regular products with two operands

        // ega2d
        generate_and_print_ega2d_gpr();    // geometric product
        generate_and_print_ega2d_wdg();    // wegde product
        generate_and_print_ega2d_lcontr(); // left contraction
        generate_and_print_ega2d_rcontr(); // right contraction
        generate_and_print_ega2d_dot();    // scalar product

        generate_and_print_ega2d_rwdg(); // regressive wedge product


        // ega3d
        generate_and_print_ega3d_gpr();    // geometric product
        generate_and_print_ega3d_wdg();    // wegde product
        generate_and_print_ega3d_lcontr(); // left contraction
        generate_and_print_ega3d_rcontr(); // right contraction
        generate_and_print_ega3d_dot();    // scalar product

        generate_and_print_ega3d_rwdg(); // regressive wedge product


        // ega2dp
        generate_and_print_ega2dp_gpr();    // geometric product
        generate_and_print_ega2dp_wdg();    // wegde product
        generate_and_print_ega2dp_lcontr(); // left contraction
        generate_and_print_ega2dp_rcontr(); // right contraction
        generate_and_print_ega2dp_dot();    // scalar product

        generate_and_print_ega2dp_rdot(); // regressive scalar product
        generate_and_print_ega2dp_rwdg(); // regressive wedge product

        // sandwich products of type rotor * object * rev(rotor)
        generate_and_print_ega2d_rotor();
        generate_and_print_ega3d_rotor();
    }

    catch (ParseError const& e) {
        fmt::println("Parse error: {}", e.what());
    }
    catch (std::exception const& e) {
        fmt::println("Exception: {}", e.what());
    }

    return 0;
}
