"""The geodesy surface of ga_usr_geodesics.hpp, through the bindings.

The C++ test suite gates the mathematics (ga_appl3dp_test, "PGA3DP: coordinate
transformation" and "PGA3DP: projections onto the tangent plane"); this file checks that
the Python surface reaches all of it and returns the same numbers: the ellipsoid and its
constants, geodetic <-> ECEF <-> ENU, the projections onto the tangent plane (including
`unproject`, hand-bound to return None beyond the limb), and the transverse Mercator.
References are closed forms and PROJ values, not values copied from the C++ output.
"""

import math

import pytest

import ga_py
from ga_py import ega, pga

DEG = math.pi / 180.0
BERLIN = ega.geo_pos(52.523765 * DEG, 13.399603 * DEG, 0.0)
KINDS = [pga.azimuthal.orthographic, pga.azimuthal.gnomonic,
         pga.azimuthal.stereographic, pga.azimuthal.perspective]
OPTIONS = pga.azimuthal_options(35786.0e3, 0.0)  # the C++ defaults, passed explicitly


def test_ellipsoid_and_constants():
    w = ga_py.wgs84
    assert w.r_equator == 6378137.0
    assert abs(1.0 / w.flattening() - 298.257223563) < 1e-6
    assert ga_py.grs80 != w and ga_py.grs80 == ega.ellipsoid(w.r_equator,
                                                             ga_py.grs80.r_pole)
    sphere = ega.ellipsoid(6371000.0, 6371000.0)
    assert sphere.e_sq() == 0.0
    assert sphere.radius(0.7) == pytest.approx(6371000.0, rel=1e-12)  # an ulp off on MSVC
    assert "ellipsoid(" in repr(w)


def test_geodetic_ecef_round_trip():
    P = pga.geo_to_ecef(BERLIN, ga_py.wgs84)
    q = pga.ecef_to_geo(P, ga_py.wgs84)
    assert abs(q.lat - BERLIN.lat) < 1e-12 and abs(q.lon - BERLIN.lon) < 1e-12
    assert abs(q.height) < 1e-6


def test_enu_motor_carries_the_origin_to_the_station():
    M = pga.enu_to_ecef_motor(BERLIN, ga_py.wgs84)
    S = pga.geo_to_ecef(BERLIN, ga_py.wgs84)
    X = pga.unitize(pga.move3dp(pga.O_3dp, M))
    assert max(abs(X.x - S.x), abs(X.y - S.y), abs(X.z - S.z)) < 1e-6
    f = pga.enu_at(BERLIN, ga_py.wgs84)
    assert abs(float(pga.dot(f.east, f.north))) < 1e-15


@pytest.mark.parametrize("kind", KINDS)
def test_station_is_the_map_origin_and_unproject_inverts(kind):
    m = pga.make_azimuthal_map(BERLIN, kind, ga_py.wgs84, OPTIONS)
    o = pga.project(m, BERLIN)
    assert math.hypot(o.x, o.y) < 1e-6
    p = ega.geo_pos(48.0 * DEG, 11.0 * DEG, 0.0)
    q = pga.unproject(m, pga.project(m, p))
    assert q is not None
    assert abs(q.lat - p.lat) < 1e-11 and abs(q.lon - p.lon) < 1e-11


def test_orthographic_on_a_sphere_is_snyder():
    R = 6371000.0
    sphere = ega.ellipsoid(R, R)
    st = ega.geo_pos(52.5 * DEG, 13.4 * DEG, 0.0)
    p = ega.geo_pos(45.0 * DEG, 20.0 * DEG, 0.0)
    m = pga.make_azimuthal_map(st, pga.azimuthal.orthographic, sphere, OPTIONS)
    x = pga.project(m, p)
    dl = p.lon - st.lon
    xs = R * math.cos(p.lat) * math.sin(dl)
    ys = R * (math.cos(st.lat) * math.sin(p.lat) -
              math.sin(st.lat) * math.cos(p.lat) * math.cos(dl))
    assert math.hypot(x.x - xs, x.y - ys) < 1e-6


def test_beyond_the_limb_there_is_no_position():
    m = pga.make_azimuthal_map(BERLIN, pga.azimuthal.orthographic, ga_py.wgs84, OPTIONS)
    assert pga.unproject(m, ega.vec2d(7.0e6, 0.0)) is None
    far = pga.to_local(m, pga.geo_to_ecef(ega.geo_pos(-40.0 * DEG, -160.0 * DEG, 0.0),
                                          ga_py.wgs84))
    assert not pga.is_drawable(m, far)
    pts = pga.limit_outline(m, 36)
    assert len(pts) == 36  # the edge of the drawable region, on the ellipsoid


def test_densify_and_clip():
    a = ega.geo_pos(50.0 * DEG, 10.0 * DEG, 0.0)
    b = ega.geo_pos(50.0 * DEG, 20.0 * DEG, 0.0)
    assert len(ga_py.densify(a, b, 1.0 * DEG)) == 9
    m = pga.make_azimuthal_map(BERLIN, pga.azimuthal.orthographic, ga_py.wgs84, OPTIONS)
    A = pga.to_local(m, pga.geo_to_ecef(ega.geo_pos(10.0 * DEG, 60.0 * DEG, 0.0),
                                        ga_py.wgs84))
    B = pga.to_local(m, pga.geo_to_ecef(ega.geo_pos(-10.0 * DEG, 120.0 * DEG, 0.0),
                                        ga_py.wgs84))
    Q = pga.clip_at(A, B, m.limit)
    assert abs(float(pga.wdg(Q, m.limit))) < 1e-6 * abs(float(pga.wdg(A, m.limit)))


def test_transverse_mercator_against_proj():
    # PROJ 9.9: cs2cs +proj=longlat +ellps=GRS80 +to +proj=utm +zone=32 +ellps=GRS80
    utm = ga_py.utm_zone(32, True)
    for lon, lat, e, n in [(9.177734, 48.775817, 513058.154680, 5402550.141164),
                           (6.0, 48.0, 276224.084601, 5320655.789070),
                           (13.0, -33.0, 873787.712841, -3658400.226820)]:
        x = ga_py.tm_forward(ega.geo_pos(lat * DEG, lon * DEG, 0.0), utm, ga_py.grs80)
        assert math.hypot(x.x - e, x.y - n) < 1e-5
        q = ga_py.tm_inverse(x, utm, ga_py.grs80)
        assert abs(q.lat - lat * DEG) < 1e-14 and abs(q.lon - lon * DEG) < 1e-14
    assert utm.k0 == 0.9996 and utm.false_easting == 500000.0
    with pytest.raises(Exception):
        ga_py.utm_zone(0, True)


def test_lambert_azimuthal_equal_area_against_proj():
    # PROJ 9.9: cs2cs +proj=longlat +ellps=GRS80 +to +proj=laea +lat_0=52 +lon_0=10
    #           +x_0=4321000 +y_0=3210000 +ellps=GRS80  (EPSG:3035, easting first)
    eu = ga_py.etrs89_laea()
    assert eu.false_easting == 4321000.0 and abs(eu.lat0 - 52.0 * DEG) < 1e-15
    for lon, lat, e, n in [(13.4, 52.52, 4551697.475620, 3273252.407167),
                           (-9.14, 38.72, 2665402.840046, 1946531.124408),
                           (33.38, 35.17, 6436679.629092, 1670516.596678)]:
        p = ega.geo_pos(lat * DEG, lon * DEG, 0.0)
        x = ga_py.laea_forward(p, eu, ga_py.grs80)
        assert math.hypot(x.x - e, x.y - n) < 1e-5
        q = ga_py.laea_inverse(x, eu, ga_py.grs80)
        assert abs(q.lat - p.lat) < 1e-14 and abs(q.lon - p.lon) < 1e-14
    c = ga_py.laea_at(BERLIN)
    o = ga_py.laea_forward(BERLIN, c, ga_py.wgs84)
    assert math.hypot(o.x, o.y) < 1e-6
    assert "lambert_azimuthal(" in repr(eu)
    with pytest.raises(Exception):  # the centre's antipode has no image
        ga_py.laea_forward(ega.geo_pos(-52.0 * DEG, -170.0 * DEG, 0.0), eu, ga_py.grs80)


def test_helmert_both_conventions_against_proj():
    # PROJ 9.9: cct +proj=helmert +x=598.1 +y=73.7 +z=418.2 +rx=0.202 +ry=0.045
    #           +rz=-2.455 +s=6.7 +convention=<c> +exact
    h = pga.make_helmert7(598.1, 73.7, 418.2, 0.202, 0.045, -2.455, 6.7)
    assert abs(h.ds - 6.7e-6) < 1e-18 and abs(h.rz - (-2.455 / 3600.0) * DEG) < 1e-18
    X = pga.vec3dp(4157222.543, 664789.307, 4774952.099, 1.0)
    for f, ref in [(pga.helmert_position_vector,
                    (4157857.450329, 664813.304462, 4775402.035203)),
                   (pga.helmert_coordinate_frame,
                    (4157839.541809, 664921.617609, 4775402.547100))]:
        Y = f(X, h)
        assert math.hypot(Y.x - ref[0], Y.y - ref[1], Y.z - ref[2]) < 1e-4
    pv, cf = pga.helmert_position_vector(X, h), pga.helmert_coordinate_frame(X, h)
    assert math.hypot(pv.x - cf.x, pv.y - cf.y, pv.z - cf.z) > 100.0  # not interchangeable
    I = pga.helmert_position_vector(X, pga.helmert7())
    assert math.hypot(I.x - X.x, I.y - X.y, I.z - X.z) < 1e-8
    assert "helmert7(" in repr(h)
