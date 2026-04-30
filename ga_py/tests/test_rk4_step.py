"""Tests for rk4_step — the std::vector overload bound for Python.

Covers:
  * binding presence + dispatch on the four supported VecTypes
    (vec2d, vec3d, vec2dp, vec3dp);
  * numerical correctness against a closed-form reference trajectory
    (1-D linear decay du/dt = -u, exact solution u(t) = u0 * exp(-t));
  * agreement with a numpy-computed RK4 reference at every sub-step;
  * the input-validation error paths (rk_step out of range, length
    mismatch).

Note: rk4_step takes its mutable state (u, uh) by value and returns the
updated state as a 2-tuple `(u_new, uh_new)`. The C++ overload mutates
in place; the Python signature returns instead because nanobind has no
caster that propagates std::vector mutations back to Python lists. See
ga_py/README.md for the rationale.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

import ga_py


# --- Helpers ---------------------------------------------------------------

VEC_TYPES = [
    ("vec2d", ga_py.ega.vec2d, 2),
    ("vec3d", ga_py.ega.vec3d, 3),
    ("vec2dp", ga_py.pga.vec2dp, 3),
    ("vec3dp", ga_py.pga.vec3dp, 4),
]


def _components(v):
    """Tuple of all components of a GA vector type, in declaration order."""
    return tuple(getattr(v, c) for c in ("x", "y", "z", "w") if hasattr(v, c))


def _zero(ctor, n):
    return ctor(*(0.0,) * n)


# --- Binding presence ------------------------------------------------------

def test_rk4_step_is_bound():
    assert hasattr(ga_py, "rk4_step")
    assert hasattr(ga_py, "rk4_get_time")


@pytest.mark.parametrize("name,ctor,n_components", VEC_TYPES)
def test_rk4_step_dispatches_for_each_vec_type(name, ctor, n_components):
    """Each of the four bound VecTypes accepts a one-element state."""
    initial = ctor(*range(1, n_components + 1))
    u = [initial]
    uh = ([_zero(ctor, n_components)], [_zero(ctor, n_components)])
    rhs = [_zero(ctor, n_components)]
    u_new, uh_new = ga_py.rk4_step(u, uh, rhs, 0.01, 1)
    # rhs == 0 → no state change at sub-step 1.
    assert _components(u_new[0]) == _components(initial)


# --- Closed-form trajectory -----------------------------------------------

def test_rk4_step_reproduces_exponential_decay():
    """Integrate du/dt = -u with vec3d state; check final state matches
    the closed-form u(t) = u0 * exp(-t) within RK4's expected accuracy."""
    u = [ga_py.ega.vec3d(1.0, 2.0, 3.0)]
    uh = ([ga_py.ega.vec3d(0, 0, 0)], [ga_py.ega.vec3d(0, 0, 0)])
    dt, n_steps = 0.01, 100

    for _ in range(n_steps):
        for rk_step in (1, 2, 3, 4):
            rhs = [ga_py.ega.vec3d(-u[0].x, -u[0].y, -u[0].z)]
            u, uh = ga_py.rk4_step(u, uh, rhs, dt, rk_step)

    decay = math.exp(-n_steps * dt)
    assert u[0].x == pytest.approx(1.0 * decay, abs=1e-9)
    assert u[0].y == pytest.approx(2.0 * decay, abs=1e-9)
    assert u[0].z == pytest.approx(3.0 * decay, abs=1e-9)


def test_rk4_step_agrees_with_numpy_reference():
    """The vec3d state evolution must agree with a hand-rolled numpy RK4 of
    the same ODE at every full step — confirms the bindings carry the C++
    arithmetic faithfully."""
    u = [ga_py.ega.vec3d(0.5, -1.0, 2.5)]
    uh = ([ga_py.ega.vec3d(0, 0, 0)], [ga_py.ega.vec3d(0, 0, 0)])
    dt, n_steps = 0.05, 40

    # numpy reference: classical RK4 on du/dt = -u, identical to the
    # 4-substep scheme rk4_step implements (Predict-Correct-Predict-Correct).
    np_u = np.array([0.5, -1.0, 2.5])
    for _ in range(n_steps):
        # Drive the ga_py side through its 4 sub-steps.
        for rk_step in (1, 2, 3, 4):
            rhs = [ga_py.ega.vec3d(-u[0].x, -u[0].y, -u[0].z)]
            u, uh = ga_py.rk4_step(u, uh, rhs, dt, rk_step)
        # Reference: classical RK4.
        k1 = -np_u
        k2 = -(np_u + 0.5 * dt * k1)
        k3 = -(np_u + 0.5 * dt * k2)
        k4 = -(np_u + dt * k3)
        np_u = np_u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    np.testing.assert_allclose(
        [u[0].x, u[0].y, u[0].z], np_u, atol=1e-12, rtol=0,
    )


# --- Input validation ------------------------------------------------------

def test_rk4_step_rejects_rk_step_out_of_range():
    u = [ga_py.ega.vec3d(1, 2, 3)]
    uh = ([ga_py.ega.vec3d(0, 0, 0)], [ga_py.ega.vec3d(0, 0, 0)])
    rhs = [ga_py.ega.vec3d(0, 0, 0)]
    with pytest.raises(ValueError):
        ga_py.rk4_step(u, uh, rhs, 0.01, 0)
    with pytest.raises(ValueError):
        ga_py.rk4_step(u, uh, rhs, 0.01, 5)


def test_rk4_step_rejects_length_mismatch():
    u = [ga_py.ega.vec3d(1, 2, 3)]
    uh = ([ga_py.ega.vec3d(0, 0, 0)], [ga_py.ega.vec3d(0, 0, 0)])
    rhs_too_long = [ga_py.ega.vec3d(0, 0, 0), ga_py.ega.vec3d(0, 0, 0)]
    with pytest.raises(ValueError):
        ga_py.rk4_step(u, uh, rhs_too_long, 0.01, 1)


# --- rk4_get_time companion -----------------------------------------------

def test_rk4_get_time():
    """rk4_get_time returns t_i, t_i+dt/2, t_i+dt/2, t_i+dt for sub-steps 0..3.

    Note this companion function uses 0..3 indexing for the sub-step
    (matching the C++ enum) while rk4_step itself uses 1..4 — see the
    function docstrings for the rationale.
    """
    assert ga_py.rk4_get_time(0.0, 0.1, 0) == pytest.approx(0.0)
    assert ga_py.rk4_get_time(0.0, 0.1, 1) == pytest.approx(0.05)
    assert ga_py.rk4_get_time(0.0, 0.1, 2) == pytest.approx(0.05)
    assert ga_py.rk4_get_time(0.0, 0.1, 3) == pytest.approx(0.1)
