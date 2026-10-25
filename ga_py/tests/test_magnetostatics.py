"""The magnetostatics of ga_usr_magnetostatics.hpp, through the bindings.

The C++ test suites gate the mathematics (ga_appl2dp_test, "magnetostatics 2d", and
ga_appl3dp_test, "magnetostatics 3d"); this file checks that the Python surface reaches
the data types and the functions and returns the same physics. References are closed
forms written out classically here, not values copied from the C++ output. The span
overloads (b_line_currents, b_segments) and the templated stress_wrench are not bound;
their sums are one Python loop over the bound kernels.
"""

import math

import ga_py
from ga_py import ega

MU0 = 4.0e-7 * math.pi  # mu0 is a primitive constant and, like pi, not bound


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
