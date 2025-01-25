// author: Daniel Hug, 2024

#include "ga_prd_ega2d.hpp"

int main()
{
    generate_and_print_product("ega2d geometric product", gpr_ega2d_basis_tab,
                               mv_ega2d_basis, mv_ega2d_lcoeff, mv_ega2d_rcoeff);
    generate_and_print_product("ega2d wegde product", wdg_ega2d_basis_tab, mv_ega2d_basis,
                               mv_ega2d_lcoeff, mv_ega2d_rcoeff);
    generate_and_print_product("ega2d dot product", dot_ega2d_basis_tab, mv_ega2d_basis,
                               mv_ega2d_lcoeff, mv_ega2d_rcoeff);
    generate_and_print_product("ega2d left contraction", lcontr_ega2d_basis_tab,
                               mv_ega2d_basis, mv_ega2d_lcoeff, mv_ega2d_rcoeff);
    generate_and_print_product("ega2d right contraction", rcontr_ega2d_basis_tab,
                               mv_ega2d_basis, mv_ega2d_lcoeff, mv_ega2d_rcoeff);

    return 0;
}
