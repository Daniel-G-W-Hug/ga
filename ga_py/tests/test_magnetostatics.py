"""The magnetostatics of ga_usr_magnetostatics.hpp, through the bindings.

The C++ test suites gate the mathematics (ga_appl2dp_test, "magnetostatics 2d", and
ga_appl3dp_test, "magnetostatics 3d"); this file checks that the Python surface reaches
the data types and the functions and returns the same physics. References are closed
forms written out classically here, not values copied from the C++ output. The span
overloads (b_line_currents, b_segments) and the templated stress_wrench and box_wrench
are not bound; their sums are one Python loop over the bound kernels. mu0 is bound
(since 2026-11-28); MU0 below is the independent value it is checked against.
"""

import math

import pytest

import ga_py
from ga_py import ega

MU0 = 4.0e-7 * math.pi  # the reference value; ga_py.mu0 is checked against it


def norm(v):
    return math.sqrt(sum(c * c for c in v))


def test_line_current_field_and_repr():
    c = ega.line_current()
    c.p = ega.vec2d(0.3, -0.2)
    c.I = ega.pscalar2d(250.0)
    assert "line_current(" in repr(c)
    r = ega.vec2d(1.0, 0.0)
    b = ga_py.b_line_current(r, c)
    dx, dy = r.x - 0.3, r.y + 0.2
    k = MU0 * 250.0 / (2.0 * math.pi * (dx * dx + dy * dy))
    assert norm((b.x + k * dy, b.y - k * dx)) < 1e-14 * k * math.hypot(dx, dy)


def test_traction_is_maxwells_tensor():
    for b, n in [((0.9, -0.4, 0.3), (0.0, 0.6, 0.8)), ((0.1, 1.2, -0.7), (1.0, 0.0, 0.0))]:
        t = ga_py.traction(ega.vec3d(*b), ega.vec3d(*n))
        bn = sum(x * y for x, y in zip(b, n))
        b2 = sum(x * x for x in b)
        ref = [(bn * bi - 0.5 * b2 * ni) / MU0 for bi, ni in zip(b, n)]
        assert norm((t.x - ref[0], t.y - ref[1], t.z - ref[2])) < 1e-14 * norm(ref)


def test_long_segment_is_the_line_current():
    s = ega.current_segment()
    s.a = ega.vec3d(0.0, 0.0, -1.0e4)
    s.b = ega.vec3d(0.0, 0.0, 1.0e4)
    s.I = 100.0
    assert "current_segment(" in repr(s)
    b = ga_py.b_segment(ega.vec3d(0.5, 0.0, 0.0), s)
    # mu0 I / (2 pi d) along +y (right-handed about +z)
    assert abs(b.y - MU0 * 100.0 / (2.0 * math.pi * 0.5)) < 1e-8 * abs(b.y)
    assert abs(b.x) < 1e-20 and abs(b.z) < 1e-20


def test_dipole_both_conventions():
    mu = ega.vec3d(0.2, -0.5, 1.0)
    r = ega.vec3d(0.3, 0.1, -0.4)
    b = ga_py.b_dipole(r, mu)
    rn = norm((r.x, r.y, r.z))
    h = (r.x / rn, r.y / rn, r.z / rn)
    md = mu.x * h[0] + mu.y * h[1] + mu.z * h[2]
    c = MU0 / (4.0 * math.pi * rn**3)
    ref = [c * (3.0 * md * hi - mi) for hi, mi in zip(h, (mu.x, mu.y, mu.z))]
    assert norm((b.x - ref[0], b.y - ref[1], b.z - ref[2])) < 1e-14 * norm(ref)
    # the bivector convention: the same expression, dual(b) = B
    B = ga_py.b_dipole(r, ega.dual(mu))
    assert norm((B.x - b.x, B.y - b.y, B.z - b.z)) < 1e-14 * norm(ref)


def test_quadrature_rules():
    contour = ga_py.circle_contour(ega.vec2d(0.0, 0.0), 2.0, 64)
    assert abs(sum(q.ds for q in contour) - 4.0 * math.pi) < 1e-12
    assert "contour_point(" in repr(contour[0])
    surface = ga_py.sphere_surface(ega.vec3d(0.0, 0.0, 0.0), 2.0, 2, 16)
    assert len(surface) == 2 * 12 * 16
    assert abs(sum(q.dA for q in surface) - 16.0 * math.pi) < 1e-12
    assert "surface_point(" in repr(surface[0])


def box(c, h, m, R=None):
    return ega.magnetised_box(ega.vec3d(*c), ega.vec3d(*h),
                              R if R is not None else ega.mvec3d_e(1.0, 0.0, 0.0, 0.0),
                              ega.vec3d(*m))


def test_magnetised_box_fields_and_repr():
    b = ega.magnetised_box()
    assert "magnetised_box(" in repr(b)
    assert b.R.c0 == 1.0  # the default orientation is the identity rotor, as in C++
    b = box((0.1, 0.2, 0.3), (0.01, 0.02, 0.03), (4.0, 5.0, 6.0))
    assert (b.c.x, b.h.y, b.m.z) == (0.1, 0.02, 6.0)


def test_magnetised_box_field_inside_and_outside():
    m = 1.0e6
    cube = box((0.0, 0.0, 0.0), (0.005, 0.005, 0.005), (0.0, 0.0, m))
    # inside a cube the demagnetising field is -m/3: b = (2/3) mu0 m at its centre
    bc = ga_py.b_box(ega.vec3d(0.0, 0.0, 0.0), cube)
    assert abs(bc.z - 2.0 / 3.0 * MU0 * m) < 1e-14 * MU0 * m
    # outside it approaches its dipole mu = 8 h_x h_y h_z m at (a/r)^4: 16x per doubling
    mu = ega.vec3d(0.0, 0.0, 8 * 0.005**3 * m)

    def dev(z):
        r = ega.vec3d(0.0, 0.0, z)
        return abs(ga_py.b_box(r, cube).z / ga_py.b_dipole(r, mu).z - 1.0)

    assert 15.5 < dev(0.1) / dev(0.2) < 16.5


def test_magnetised_box_orientation_reaches_the_field():
    # a flat box turned 90 degrees about z is the box with its x and y extents swapped
    # (the magnetisation is in world coordinates), and not the unturned box
    Rz = ega.get_rotor(ega.bivec3d(0.0, 0.0, 1.0), math.pi / 2)
    m = (0.3e6, 0.8e6, 0.5e6)
    turned = box((0.01, 0.0, 0.0), (0.01, 0.004, 0.002), m, Rz)
    swapped = box((0.01, 0.0, 0.0), (0.004, 0.01, 0.002), m)
    plain = box((0.01, 0.0, 0.0), (0.01, 0.004, 0.002), m)
    r = ega.vec3d(0.02, 0.013, 0.007)
    a, b, c = (ga_py.b_box(r, x) for x in (turned, swapped, plain))
    ref = norm((b.x, b.y, b.z))
    assert norm((a.x - b.x, a.y - b.y, a.z - b.z)) < 1e-12 * ref
    assert norm((a.x - c.x, a.y - c.y, a.z - c.z)) > 0.1 * ref


def test_mu0_and_the_tolerances_are_bound():
    assert isinstance(ga_py.mu0, float) and ga_py.mu0 == MU0
    # eps and eps_congruent: 5 and 4500 machine epsilons of value_t (double)
    assert ga_py.eps == 5.0 * 2.0**-52
    assert ga_py.eps_congruent == 4.5e3 * 2.0**-52


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def test_dipole_force_and_wrench():
    # the force against the classical closed form written out here; the wrench's weight,
    # its moment about the second dipole (the couple of the torque mu2 x b1), and action
    # = reaction as one bivector
    from ga_py import pga

    p1, p2 = (0.01, -0.02, 0.0), (0.07, 0.015, 0.03)
    mu1, mu2 = (0.3, -0.2, 1.1), (-0.6, 0.9, 0.4)
    r = [b - a for a, b in zip(p1, p2)]
    rn = norm(r)
    h = [x / rn for x in r]

    def dot(a, b):
        return sum(x * y for x, y in zip(a, b))

    ci, cj = dot(h, mu1), dot(h, mu2)
    c = 3.0 * MU0 / (4.0 * math.pi * rn**4)
    ref = [c * (dot(mu1, mu2) * hk + cj * a + ci * b - 5.0 * ci * cj * hk)
           for hk, a, b in zip(h, mu1, mu2)]
    f = ga_py.dipole_force(ega.vec3d(*r), ega.vec3d(*mu1), ega.vec3d(*mu2))
    assert norm((f.x - ref[0], f.y - ref[1], f.z - ref[2])) < 1e-14 * norm(ref)
    fb = ga_py.dipole_force(ega.vec3d(*r), ega.dual(ega.vec3d(*mu1)),
                            ega.dual(ega.vec3d(*mu2)))
    assert (fb.x, fb.y, fb.z) == (f.x, f.y, f.z)

    W12 = ga_py.dipole_wrench(ega.vec3d(*p1), ega.vec3d(*mu1), ega.vec3d(*p2),
                              ega.vec3d(*mu2))
    W21 = ga_py.dipole_wrench(ega.vec3d(*p2), ega.vec3d(*mu2), ega.vec3d(*p1),
                              ega.vec3d(*mu1))
    assert norm((W12.vx - f.x, W12.vy - f.y, W12.vz - f.z)) < 1e-15 * norm(ref)
    b1 = ga_py.b_dipole(ega.vec3d(*r), ega.vec3d(*mu1))
    T2 = cross(mu2, (b1.x, b1.y, b1.z))
    M2 = pga.moment_about(pga.vec3dp(*p2, 1.0), W12)
    assert norm((M2.mx - T2[0], M2.my - T2[1], M2.mz - T2[2])) < 1e-13 * norm(T2)
    s = [getattr(W12, k) + getattr(W21, k) for k in ("vx", "vy", "vz", "mx", "my", "mz")]
    assert max(abs(x) for x in s) < 1e-13 * norm(ref)


def test_couple_is_a_torque_of_zero_weight():
    from ga_py import pga

    C = pga.couple(ega.bivec3d(0.4, -1.3, 0.25))
    assert (C.vx, C.vy, C.vz, C.mx, C.my, C.mz) == (0.0, 0.0, 0.0, 0.4, -1.3, 0.25)
    R = pga.vec3dp(-2.0, 0.8, 0.5, 1.0)
    M = pga.moment_about(R, C)
    assert (M.mx, M.my, M.mz) == pytest.approx((0.4, -1.3, 0.25))
    C2 = pga.couple(ega.pscalar2d(0.75))
    assert (C2.x, C2.y, C2.z) == (0.0, 0.0, 0.75)

