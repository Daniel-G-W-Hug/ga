// Copyright 2024-2026, Daniel Hug. All rights reserved.
// Licensed under the terms specified in LICENSE.txt file.

#include "ga_test_arith_check.hpp"

#include <cmath>

/////////////////////////////////////////////////////////////////////////////////////////
// CGA contracts: the arithmetic of every value type, the multivector constructors from
// parts, the four-part round/flat split, and the operations of ga_cga{2,3}dc_ops*.hpp
// that the topical tests did not reach -- each pinned against a documented identity or
// an independently known answer
/////////////////////////////////////////////////////////////////////////////////////////

TEST_SUITE("CGA contracts")
{
    // the four parts of a CGA object partition its coefficients: they sum to the object,
    // their norms are the plain norms of the parts, and the conformal conjugate is the
    // documented signed sum  cconj(u) = rb + rw - fb - fw
    auto const check_split = [](auto const& u) {
        auto const rb = round_bulk(u);
        auto const rw = round_weight(u);
        auto const fb = flat_bulk(u);
        auto const fw = flat_weight(u);
        CHECK(rb + rw + fb + fw == u);
        CHECK(cconj(u) == rb + rw - fb - fw);
        CHECK(cconj(cconj(u)) == u);
        using hd::ga::detail::coeff_sq;
        CHECK(round_bulk_nrm_sq(u) == doctest::Approx(coeff_sq(rb)));
        CHECK(round_weight_nrm_sq(u) == doctest::Approx(coeff_sq(rw)));
        CHECK(flat_bulk_nrm_sq(u) == doctest::Approx(coeff_sq(fb)));
        CHECK(flat_weight_nrm_sq(u) == doctest::Approx(coeff_sq(fw)));
        CHECK(coeff_sq(u) ==
              doctest::Approx(coeff_sq(rb) + coeff_sq(rw) + coeff_sq(fb) + coeff_sq(fw)));
        CHECK(round_bulk_nrm(u) == doctest::Approx(std::sqrt(round_bulk_nrm_sq(u))));
        CHECK(round_weight_nrm(u) == doctest::Approx(std::sqrt(round_weight_nrm_sq(u))));
        CHECK(flat_bulk_nrm(u) == doctest::Approx(std::sqrt(flat_bulk_nrm_sq(u))));
        CHECK(flat_weight_nrm(u) == doctest::Approx(std::sqrt(flat_weight_nrm_sq(u))));
        CHECK(center_nrm(u) == doctest::Approx(std::sqrt(center_nrm_sq(u))));
    };

    TEST_CASE("cga2dc / cga3dc: arithmetic contract of every value type")
    {
        arith_check::check<scalar2dc>("scalar2dc");
        arith_check::check<vec2dc>("vec2dc");
        arith_check::check<bivec2dc>("bivec2dc");
        arith_check::check<trivec2dc>("trivec2dc");
        arith_check::check<pscalar2dc>("pscalar2dc");
        arith_check::check<mvec2dc_e>("mvec2dc_e");
        arith_check::check<mvec2dc_u>("mvec2dc_u");
        arith_check::check<mvec2dc>("mvec2dc");

        arith_check::check<scalar3dc>("scalar3dc");
        arith_check::check<vec3dc>("vec3dc");
        arith_check::check<bivec3dc>("bivec3dc");
        arith_check::check<trivec3dc>("trivec3dc");
        arith_check::check<quadvec3dc>("quadvec3dc");
        arith_check::check<pscalar3dc>("pscalar3dc");
        arith_check::check<mvec3dc_e>("mvec3dc_e");
        arith_check::check<mvec3dc_u>("mvec3dc_u");
        arith_check::check<mvec3dc>("mvec3dc");
    }

    TEST_CASE("cga2dc / cga3dc: a multivector built from parts is the sum of the parts")
    {
        vec2dc const v2(1.0, -2.0, 0.5, 3.0);
        bivec2dc const b2(1.0, 2.0, 3.0, 4.0, 5.0, 6.0);
        trivec2dc const t2(2.0, -1.0, 0.5, 4.0);
        CHECK(mvec2dc_u(v2) + mvec2dc_u(t2) == mvec2dc_u(v2, t2));
        CHECK(mvec2dc(scalar2dc(2.0), b2, pscalar2dc(-1.0)) ==
              mvec2dc(scalar2dc(2.0)) + mvec2dc(b2) + mvec2dc(pscalar2dc(-1.0)));

        vec3dc const v(1.0, -2.0, 0.5, 3.0, -1.5);
        bivec3dc const b(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
        trivec3dc const t(-1, 2, -3, 4, -5, 6, -7, 8, -9, 10);
        quadvec3dc const q(0.5, 1.5, -2.5, 3.5, 4.5);
        scalar3dc const s(2.0);
        pscalar3dc const ps(-0.75);

        CHECK(mvec3dc_u(v) + mvec3dc_u(t) == mvec3dc_u(v, t));
        CHECK(mvec3dc_u(t) + mvec3dc_u(ps) == mvec3dc_u(t, ps));
        CHECK(mvec3dc_u(v) + mvec3dc_u(ps) == mvec3dc_u(v, ps));
        CHECK(mvec3dc_u(v, t) + mvec3dc_u(ps) == mvec3dc_u(v, ps) + mvec3dc_u(t));
        CHECK(mvec3dc_e(b) + mvec3dc_e(q) == mvec3dc_e(b, q));
        CHECK(mvec3dc_e(s) + mvec3dc_e(q) == mvec3dc_e(s, q));
        CHECK(mvec3dc(s, b, q) == mvec3dc(s) + mvec3dc(b) + mvec3dc(q));
        CHECK(mvec3dc(v, t, ps) == mvec3dc(v) + mvec3dc(t) + mvec3dc(ps));
    }

    TEST_CASE("cga2dc / cga3dc: the round/flat split of every graded type")
    {
        check_split(vec2dc(1.0, -2.0, 0.5, 3.0));
        check_split(bivec2dc(1.0, 2.0, 3.0, 4.0, 5.0, 6.0));
        check_split(trivec2dc(2.0, -1.0, 0.5, 4.0));

        check_split(vec3dc(1.0, -2.0, 0.5, 3.0, -1.5));
        check_split(bivec3dc(1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
        check_split(trivec3dc(-1, 2, -3, 4, -5, 6, -7, 8, -9, 10));
        check_split(quadvec3dc(0.5, 1.5, -2.5, 3.5, 4.5));
    }

    TEST_CASE("cga2dc / cga3dc: the center norm is the distance of the centre")
    {
        // documented for a unitized round object u: dot(u, cconj(u)) and center_nrm_sq(u)
        // are the squared distance origin <-> centre
        auto const p3 = unitize(round_point3dc(1.0, 2.0, 2.0, 0.5)); // |c| = 3
        CHECK(center_nrm_sq(p3) == doctest::Approx(9.0));
        CHECK(center_nrm(p3) == doctest::Approx(3.0));
        CHECK(value_t(dot(p3, cconj(p3))) == doctest::Approx(9.0));

        auto const p2 = unitize(round_point2dc(3.0, 4.0, 0.25)); // |c| = 5
        CHECK(center_nrm_sq(p2) == doctest::Approx(25.0));
        CHECK(center_nrm(p2) == doctest::Approx(5.0));
        CHECK(value_t(dot(p2, cconj(p2))) == doctest::Approx(25.0));
    }

    TEST_CASE("cga2dc / cga3dc: is_degenerate compares the round weight, scale-free")
    {
        // a round point has a centre; the same coefficients without the round weight do
        // not -- whatever the object's scale
        auto const p3 = round_point3dc(1.0, 2.0, 2.0, 0.5);
        CHECK_FALSE(is_degenerate(p3));
        vec3dc const no_centre(p3.x, p3.y, p3.z, 0.0, p3.u);
        CHECK(is_degenerate(no_centre));
        CHECK(is_degenerate(1.0e6 * no_centre));
        CHECK_FALSE(is_degenerate(1.0e-6 * p3));

        auto const d = dipole3dc(1.0, 0.0, 0.0, 0.5, 0.0, 0.0, 1.0);
        CHECK_FALSE(is_degenerate(d));
        CHECK(is_degenerate(d - round_weight(d)));

        auto const p2 = round_point2dc(3.0, 4.0, 0.25);
        CHECK_FALSE(is_degenerate(p2));
        CHECK(is_degenerate(vec2dc(p2.x, p2.y, 0.0, p2.w)));
    }

    TEST_CASE("cga2dc / cga3dc: inversion is an outermorphism")
    {
        // inverting a point pair, a circle or a sphere through given points gives the one
        // through the inverted points
        auto const S = sphere3dc(0.0, 0.0, 0.0, 1.0);
        auto const p1 = round_point3dc(2.0, 0.0, 0.0, 0.0);
        auto const p2 = round_point3dc(0.0, 3.0, 0.0, 0.0);
        auto const p3 = round_point3dc(0.0, 0.0, -1.5, 0.0);
        auto const p4 = round_point3dc(1.0, 1.0, 1.0, 0.0);
        auto const i1 = invert_on(p1, S);
        auto const i2 = invert_on(p2, S);
        auto const i3 = invert_on(p3, S);
        auto const i4 = invert_on(p4, S);
        CHECK(is_congruent(invert_on(wdg(p1, p2), S), wdg(i1, i2)));
        CHECK(is_congruent(invert_on(wdg(wdg(p1, p2), p3), S), wdg(wdg(i1, i2), i3)));
        CHECK(is_congruent(invert_on(wdg(wdg(wdg(p1, p2), p3), p4), S),
                           wdg(wdg(wdg(i1, i2), i3), i4)));
        // and it is an involution
        auto const B = wdg(p1, p2);
        CHECK(is_congruent(invert_on(invert_on(B, S), S), B));

        auto const uc = circle2dc(0.0, 0.0, 1.0);
        auto const q1 = round_point2dc(2.0, 0.0, 0.0);
        auto const q2 = round_point2dc(0.0, -3.0, 0.0);
        CHECK(is_congruent(invert_on(wdg(q1, q2), uc),
                           wdg(invert_on(q1, uc), invert_on(q2, uc))));
    }

    TEST_CASE("cga2dc / cga3dc: unitize divides a single-coefficient round weight out "
              "signed")
    {
        // a round point, a sphere and a 2dc circle carry their round weight in ONE
        // coefficient: whatever its sign, the unitized object has round weight +1
        auto const rw = [](auto const& u) { return value_t(round_weight_nrm_sq(u)); };
        auto const p2 = -2.0 * round_point2dc(1.0, 2.0, 0.5);
        CHECK(unitize(p2) == unitize(round_point2dc(1.0, 2.0, 0.5)));
        CHECK(round_weight(unitize(p2)) == vec2dc(0.0, 0.0, 1.0, 0.0));
        auto const c2 = -3.0 * circle2dc(1.0, 2.0, 0.5);
        CHECK(unitize(c2) == unitize(circle2dc(1.0, 2.0, 0.5)));
        CHECK(round_weight(unitize(c2)) == trivec2dc(0.0, 0.0, 0.0, 1.0));
        auto const p3 = -2.0 * round_point3dc(1.0, 2.0, 3.0, 0.5);
        CHECK(unitize(p3) == unitize(round_point3dc(1.0, 2.0, 3.0, 0.5)));
        CHECK(round_weight(unitize(p3)) == vec3dc(0.0, 0.0, 0.0, 1.0, 0.0));
        auto const s3 = -4.0 * sphere3dc(1.0, 2.0, 3.0, 0.5);
        CHECK(unitize(s3) == unitize(sphere3dc(1.0, 2.0, 3.0, 0.5)));
        CHECK(round_weight(unitize(s3)) == quadvec3dc(0.0, 0.0, 0.0, 0.0, 1.0));

        // a dipole's round weight has several coefficients: its norm becomes one and
        // the orientation is kept
        auto const d = dipole3dc(1.0, 0.0, 0.0, 0.5, 0.0, 0.6, 0.8);
        CHECK(rw(unitize(-d)) == doctest::Approx(1.0));
        CHECK(unitize(-d) == -unitize(d));
    }

    TEST_CASE("cga3dc: par of a dipole is its imaginary partner")
    {
        // same centre and attitude, radius squared negated; applied twice, the identity
        auto const d = dipole3dc(1.0, 2.0, 3.0, 0.5, 0.0, 0.6, 0.8);
        CHECK(par(par(d)) == d);
        CHECK(radius_sq(par(d)) == doctest::Approx(-radius_sq(d)));
        CHECK(radius_sq(d) == doctest::Approx(0.25));
    }

    TEST_CASE("cga2dc / cga3dc: is_close on the remaining graded types")
    {
        CHECK(is_close(scalar3dc(3.0), scalar3dc(3.0 + 3.0e-14)));
        CHECK_FALSE(is_close(scalar3dc(3.0), scalar3dc(3.1)));
        bivec3dc const b(1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
        CHECK(is_close(b, b * (1.0 + 1.0e-14)));
        CHECK_FALSE(is_close(b, -b));
        trivec3dc const t(-1, 2, -3, 4, -5, 6, -7, 8, -9, 10);
        CHECK(is_close(t, t * (1.0 + 1.0e-14)));
        CHECK_FALSE(is_close(t, 2.0 * t));
        quadvec3dc const q(0.5, 1.5, -2.5, 3.5, 4.5);
        CHECK(is_close(q, q * (1.0 + 1.0e-14)));
        CHECK_FALSE(is_close(q, quadvec3dc(0.5, 1.5, -2.5, 3.5, 4.6)));

        mvec2dc_u const u(vec2dc(1.0, -2.0, 0.5, 3.0), trivec2dc(2.0, -1.0, 0.5, 4.0));
        CHECK(is_close(u, u * (1.0 + 1.0e-14)));
        CHECK_FALSE(is_close(u, -u));
    }
}
