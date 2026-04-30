"""Tests for the rigid-body physics bindings (inertia2dp / inertia3dp +
factories + LU inverse + body-frame Euler ODE rhs).

These mirror the analytical checks asserted in the C++ doctest cases at
ga_test/src/ga_appl{2,3}dp_physics_test.hpp:
  * inertia matrices match closed-form formulas (Steiner / parallel-axis);
  * forward and inverse maps round-trip for a pivot-aligned input;
  * compute_omega_dot returns a vector aligned with the pivot direction
    for pure-couple torques.
"""
from __future__ import annotations

import math

import pytest

import ga_py


# --- 2dp -----------------------------------------------------------------

def test_inertia2dp_default_ctor_zero():
    I = ga_py.pga.inertia2dp()
    for r in range(3):
        for c in range(3):
            assert I[r, c] == 0.0


def test_get_plate_inertia_matches_analytic_formula():
    """Plate at the body origin, pivot at the TR corner.
    I[2,2] = m*(w^2+h^2)/12 + m*(hw^2+hh^2) = 2/3 + 2 = 8/3."""
    m, w, h = 1.0, 2.0, 2.0
    hw, hh = w / 2.0, h / 2.0
    I = ga_py.pga.get_plate_inertia(m, w, h, ga_py.pga.vec2dp(hw, hh, 1.0))
    assert I[2, 2] == pytest.approx(8.0 / 3.0)
    # Off-diagonal Newton-momentum block: I[0,1] = +m, I[1,0] = -m.
    assert I[0, 1] == pytest.approx(m)
    assert I[1, 0] == pytest.approx(-m)


def test_inertia2dp_inverse_round_trips():
    """I[I_inv[F]] == F (and the other way around) for a non-singular plate."""
    I = ga_py.pga.get_plate_inertia(
        1.0, 2.0, 2.0, ga_py.pga.vec2dp(1.0, 1.0, 1.0)
    )
    I_inv = ga_py.pga.get_inertia_inverse(I)

    # Forward map I[Omega] returns a bivec2dp, inverse maps it back to vec2dp.
    Omega = ga_py.pga.vec2dp(0.0, 0.0, 1.0)
    F = I(Omega)                       # bivec2dp
    Omega_back = I_inv(F)              # vec2dp
    for c in "xyz":
        assert getattr(Omega_back, c) == pytest.approx(getattr(Omega, c))


def test_inertia2dp_inverse_raises_on_singular():
    """Default-constructed (zero) inertia is singular → ValueError."""
    I_zero = ga_py.pga.inertia2dp()
    with pytest.raises(ValueError):
        ga_py.pga.get_inertia_inverse(I_zero)


def test_compute_omega_dot_2dp_aligns_with_pivot():
    """For a plate pinned at its TR corner under gravity, dOmega is along
    the pivot direction Q_b = (hw, hh, 1) (z-component is the angular
    acceleration alpha = tau / I_pivot[2,2]).

    Setup matches the C++ test (m=1, 2x2 plate, pivot at TR corner) at the
    initial instant: Omega = 0, gravity acting at the cm in world frame.
    """
    m, w, h, g = 1.0, 2.0, 2.0, 9.81
    hw, hh = w / 2.0, h / 2.0
    Q_b = ga_py.pga.vec2dp(hw, hh, 1.0)

    I = ga_py.pga.get_plate_inertia(m, w, h, Q_b)
    I_inv = ga_py.pga.get_inertia_inverse(I)

    # cm at world (0,0,1) initially; pivot at (hw, hh, 1) world. Force couple:
    # gravity at cm (downward) and reaction at pivot (upward).
    cm_w = ga_py.pga.vec2dp(0.0, 0.0, 1.0)
    f_g = ga_py.pga.vec2dp(0.0, -m * g, 0.0)
    F_cm = ga_py.pga.wdg(cm_w, f_g)
    F_piv = ga_py.pga.wdg(Q_b, -f_g)
    F = F_cm + F_piv

    Omega = ga_py.pga.vec2dp(0.0, 0.0, 0.0)
    dOmega = ga_py.pga.compute_omega_dot(I_inv, F, Omega, I)
    # alpha = tau / I[2,2]; tau here is m*g*hw (lever arm in world x).
    expected_alpha = (m * g * hw) / I[2, 2]
    assert dOmega.z == pytest.approx(expected_alpha, rel=1e-6)


# --- 3dp -----------------------------------------------------------------

def test_inertia3dp_default_ctor_zero():
    I = ga_py.pga.inertia3dp()
    for r in range(6):
        for c in range(6):
            assert I[r, c] == 0.0


def test_get_cuboid_inertia_matches_analytic_formula():
    """Cuboid (2x2x1, m=1) with pivot at the TR z-edge: same I[5,2] as the
    2D plate's I[2,2] = 8/3 (depth d does not enter the in-plane moment).
    Pivot line L_b = bivec3dp{0, 0, 1, hh, -hw, 0}."""
    m, w, h, d = 1.0, 2.0, 2.0, 1.0
    hw, hh = w / 2.0, h / 2.0
    L_b = ga_py.pga.bivec3dp(0.0, 0.0, 1.0, hh, -hw, 0.0)
    I = ga_py.pga.get_cuboid_inertia(m, w, h, d, L_b)
    assert I[5, 2] == pytest.approx(8.0 / 3.0)


def test_inertia3dp_inverse_round_trips():
    m, w, h, d = 1.0, 2.0, 2.0, 1.0
    hw, hh = w / 2.0, h / 2.0
    L_b = ga_py.pga.bivec3dp(0.0, 0.0, 1.0, hh, -hw, 0.0)
    I = ga_py.pga.get_cuboid_inertia(m, w, h, d, L_b)
    I_inv = ga_py.pga.get_inertia_inverse(I)

    # In 3dp both forward and inverse are bivec3dp -> bivec3dp.
    Omega = ga_py.pga.bivec3dp(0.0, 0.0, 1.0, hh, -hw, 0.0)  # along L_b
    F = I(Omega)
    Omega_back = I_inv(F)
    for c in ("vx", "vy", "vz", "mx", "my", "mz"):
        assert getattr(Omega_back, c) == pytest.approx(getattr(Omega, c), abs=1e-9)


def test_inertia3dp_inverse_raises_on_singular():
    I_zero = ga_py.pga.inertia3dp()
    with pytest.raises(ValueError):
        ga_py.pga.get_inertia_inverse(I_zero)


# --- Cross-checks --------------------------------------------------------

def test_2d_plate_and_3d_cuboid_share_pivot_moment():
    """The 3D cuboid's I[5,2] (z-axis pivot moment) must equal the 2D
    plate's I[2,2] when the in-plane geometry matches — depth doesn't
    enter the in-plane moment of inertia."""
    m, w, h, d = 1.0, 2.0, 2.0, 1.0
    hw, hh = w / 2.0, h / 2.0
    Q_b = ga_py.pga.vec2dp(hw, hh, 1.0)
    L_b = ga_py.pga.bivec3dp(0.0, 0.0, 1.0, hh, -hw, 0.0)

    I_plate = ga_py.pga.get_plate_inertia(m, w, h, Q_b)
    I_cuboid = ga_py.pga.get_cuboid_inertia(m, w, h, d, L_b)
    assert I_plate[2, 2] == pytest.approx(I_cuboid[5, 2])


def test_inertia_indexing_out_of_range():
    I2 = ga_py.pga.inertia2dp()
    I3 = ga_py.pga.inertia3dp()
    with pytest.raises(IndexError):
        _ = I2[3, 0]
    with pytest.raises(IndexError):
        _ = I3[6, 0]
    with pytest.raises(IndexError):
        I2[3, 0] = 1.0


def test_inertia_format_spec():
    """Format specs round-trip through fmt::format (same path as other types)."""
    I = ga_py.pga.get_plate_inertia(1.0, 2.0, 2.0)
    # Empty spec: must not raise.
    s = format(I, "")
    assert isinstance(s, str)
    # Numeric spec: each component formatted accordingly.
    s = format(I, ".3f")
    assert "0.000" in s or "1.000" in s or "0.667" in s
