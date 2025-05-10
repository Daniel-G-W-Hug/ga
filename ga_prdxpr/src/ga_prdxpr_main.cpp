// Copyright 2024-2025, Daniel Hug. All rights reserved.

#include "ga_prd_ega2d.hpp"
#include "ga_prd_ega3d.hpp"

#include "ga_prd_pga2dp.hpp"
#include "ga_prd_pga3dp.hpp"

// multivector basis, product rules, simplification rules and coefficients need to be
// provided in header files included above by the user

int main()
{

    try {

        // regular products with two operands

        ///////////////////////////////////////////////////////////////////////
        // ega2d
        ///////////////////////////////////////////////////////////////////////
        generate_and_print_ega2d_gpr(); // geometric product
        generate_and_print_ega2d_cmt(); // commutator product (asymmetric part of gpr)
        generate_and_print_ega2d_wdg(); // wegde product (=outer product)
        generate_and_print_ega2d_dot(); // dot product (=inner product)

        // ega2d - regressive products
        generate_and_print_ega2d_rwdg(); // regressive wedge product

        // ega2d - interior products (contractions, expansions)
        generate_and_print_ega2d_lcontract(); // left contraction
        generate_and_print_ega2d_rcontract(); // right contraction

        generate_and_print_ega2d_lexpand(); // left expansion
        generate_and_print_ega2d_rexpand(); // right expansion
        // => expansions don't yield anything interesting in ega2d

        ///////////////////////////////////////////////////////////////////////
        // ega3d
        ///////////////////////////////////////////////////////////////////////
        generate_and_print_ega3d_gpr(); // geometric product
        generate_and_print_ega3d_cmt(); // commutator product (asymmetric part of gpr)
        generate_and_print_ega3d_wdg(); // wegde product (=outer product)
        generate_and_print_ega3d_dot(); // dot product (=inner product)

        // ega3d - regressive products
        generate_and_print_ega3d_rwdg(); // regressive wedge product

        // ega3d - interior products (contractions, expansions)
        generate_and_print_ega3d_lcontract(); // left contraction
        generate_and_print_ega3d_rcontract(); // right contraction

        generate_and_print_ega3d_lexpand(); // left expansion
        generate_and_print_ega3d_rexpand(); // right expansion
        // => expansions only have one intresting output in ega3d:
        //    expand(vec,bivec) => bivec
        //    which yields a plane that contains the vec and is perpendicular to
        //    the argument bivec, i.e. the reference plane


        ///////////////////////////////////////////////////////////////////////
        // pga2dp
        ///////////////////////////////////////////////////////////////////////
        generate_and_print_pga2dp_gpr(); // geometric product
        generate_and_print_ega2dp_cmt(); // commutator product (asymmetric part of gpr)
        generate_and_print_pga2dp_wdg(); // wegde product (=outer product)
        generate_and_print_pga2dp_dot(); // dot product (=inner product)

        // pga2dp - regressive products
        generate_and_print_pga2dp_rwdg(); // regressive wedge product
        generate_and_print_pga2dp_rdot(); // regressive dot product
        generate_and_print_pga2dp_rgpr(); // regressive geometric product

        // pga2dp - interior products (contractions, expansions)
        generate_and_print_pga2dp_rbulk_contract();   // right bulk contraction
        generate_and_print_pga2dp_rweight_contract(); // right weight contraction
        generate_and_print_pga2dp_lbulk_contract();   // left bulk contraction
        generate_and_print_pga2dp_lweight_contract(); // left weight contraction

        generate_and_print_pga2dp_rbulk_expand();   // right bulk expansion
        generate_and_print_pga2dp_rweight_expand(); // right weight expansion
        generate_and_print_pga2dp_lbulk_expand();   // left bulk expansion
        generate_and_print_pga2dp_lweight_expand(); // left weight expansion


        ///////////////////////////////////////////////////////////////////////
        // pga3dp
        ///////////////////////////////////////////////////////////////////////
        generate_and_print_pga3dp_gpr(); // geometric product
        generate_and_print_pga3dp_cmt(); // commutator product (asymmetric part of gpr)
        generate_and_print_pga3dp_wdg(); // wegde product (=outer product)
        generate_and_print_pga3dp_dot(); // dot product (=inner product)

        generate_and_print_pga3dp_rwdg(); // regressive wedge product
        generate_and_print_pga3dp_rdot(); // regressive dot product
        generate_and_print_pga3dp_rgpr(); // regressive geometric product

        // pga3dp - interior products (contractions, expansions)
        generate_and_print_pga3dp_rbulk_contract();   // right bulk contraction
        generate_and_print_pga3dp_rweight_contract(); // right weight contraction
        generate_and_print_pga3dp_lbulk_contract();   // left bulk contraction
        generate_and_print_pga3dp_lweight_contract(); // left weight contraction

        generate_and_print_pga3dp_rbulk_expand();   // right bulk expansion
        generate_and_print_pga3dp_rweight_expand(); // right weight expansion
        generate_and_print_pga3dp_lbulk_expand();   // left bulk expansion
        generate_and_print_pga3dp_lweight_expand(); // left weight expansion


        ///////////////////////////////////////////////////////////////////////
        // sandwich products
        ///////////////////////////////////////////////////////////////////////

        // ega - sandwich products of type rotor * object * rev(rotor)
        //       based on the geometric product (keeps the origin stable)
        generate_and_print_ega2d_rotor_gpr();
        generate_and_print_ega3d_rotor_gpr();

        // pga - sandwich products of type motor * object * rev(motor)
        //       based on the regressive geometric product (keeps the horizon stable)
        generate_and_print_pga2dp_motor_rgpr();
        generate_and_print_pga3dp_motor_rgpr();
    }

    catch (std::exception const& e) {
        fmt::println("Exception: {}", e.what());
    }

    return 0;
}
