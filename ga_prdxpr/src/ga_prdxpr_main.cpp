// author: Daniel Hug, 2024

#include "ga_prd_ega2d.hpp"
#include "ga_prd_ega3d.hpp"

int main()
{
    // ega2d
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

    // ega3d
    generate_and_print_product("ega3d geometric product", gpr_ega3d_basis_tab,
                               mv_ega3d_basis, mv_ega3d_lcoeff, mv_ega3d_rcoeff);
    generate_and_print_product("ega3d wedge product", wdg_ega3d_basis_tab, mv_ega3d_basis,
                               mv_ega3d_lcoeff, mv_ega3d_rcoeff);
    generate_and_print_product("ega3d dot product", dot_ega3d_basis_tab, mv_ega3d_basis,
                               mv_ega3d_lcoeff, mv_ega3d_rcoeff);
    generate_and_print_product("ega3d left contraction", lcontr_ega3d_basis_tab,
                               mv_ega3d_basis, mv_ega3d_lcoeff, mv_ega3d_rcoeff);
    generate_and_print_product("ega3d right contraction", rcontr_ega3d_basis_tab,
                               mv_ega3d_basis, mv_ega3d_lcoeff, mv_ega3d_rcoeff);

    return 0;
}
