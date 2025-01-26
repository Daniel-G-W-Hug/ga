// author: Daniel Hug, 2024

#include "ga_prd_ega2d.hpp"
#include "ga_prd_ega3d.hpp"

int main()
{
    // ega2d
    generate_and_print_ega2d_gpr();    // geometric product
    generate_and_print_ega2d_wdg();    // wegde product
    generate_and_print_ega2d_lcontr(); // left contraction
    generate_and_print_ega2d_rcontr(); // right contraction
    generate_and_print_ega2d_dot();    // scalar product

    // ega3d
    // generate_and_print_ega3d_gpr();    // geometric product
    // generate_and_print_ega3d_wdg();    // wegde product
    // generate_and_print_ega3d_lcontr(); // left contraction
    // generate_and_print_ega3d_rcontr(); // right contraction
    // generate_and_print_ega3d_dot();    // scalar product

    return 0;
}
