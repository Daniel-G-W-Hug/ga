"""Tests for the numpy buffer-protocol surface (Tier 1 of §6.3).

Each bound vec-shaped and inherited-shaped type exposes:
  * `__array__(self, dtype=None, copy=None) -> ndarray[float64]`  — copy-out
  * `__init__(self, arr: ndarray[float64, shape=(N,)])`            — copy-in

The contract is: `T(np.array(T(*xs))) == T(*xs)` component-wise, for every
bound type. We verify that explicitly per type, plus a few cross-cutting
properties (dtype, shape, batch construction, edge cases).

The values come back as the underlying C++ field declaration order — which
is the order the Python `__init__(*args)` overload also accepts, so the
round-trip is closed by construction.
"""
from __future__ import annotations

import numpy as np
import pytest

import ga_py


# Each entry: (constructor, N, sample-args).
# The N component values are passed positionally to the bound ctor; the
# resulting object's __array__ must yield exactly those values in order.
_VEC_TYPES = [
    (ga_py.ega.vec2d, 2, (1.5, 2.5)),
    (ga_py.ega.vec3d, 3, (1.0, 2.0, 3.0)),
    (ga_py.ega.vec4d, 4, (1.0, 2.0, 3.0, 4.0)),
    (ga_py.ega.bivec3d, 3, (0.1, 0.2, 0.3)),
    (ga_py.ega.bivec4d, 6, (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)),
    (ga_py.ega.trivec4d, 4, (7.0, 8.0, 9.0, 10.0)),
    (ga_py.ega.mvec2d, 4, (1.0, 2.0, 3.0, 4.0)),
    (ga_py.ega.mvec2d_e, 2, (5.0, 6.0)),
    (ga_py.ega.mvec3d, 8, tuple(float(i + 1) for i in range(8))),
    (ga_py.ega.mvec3d_e, 4, (1.0, 2.0, 3.0, 4.0)),
    (ga_py.ega.mvec3d_u, 4, (1.0, 2.0, 3.0, 4.0)),
    (ga_py.ega.mvec4d, 16, tuple(float(i + 1) for i in range(16))),
    (ga_py.ega.mvec4d_e, 8, tuple(float(i + 1) for i in range(8))),
    (ga_py.ega.mvec4d_u, 8, tuple(float(i + 1) for i in range(8))),
    (ga_py.pga.vec2dp, 3, (1.0, 2.0, 3.0)),
    (ga_py.pga.vec3dp, 4, (1.0, 2.0, 3.0, 4.0)),
    (ga_py.pga.bivec2dp, 3, (0.5, 1.5, 2.5)),
    (ga_py.pga.bivec3dp, 6, (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)),
    (ga_py.pga.trivec3dp, 4, (1.0, 2.0, 3.0, 4.0)),
    (ga_py.pga.mvec2dp, 8, tuple(float(i + 1) for i in range(8))),
    (ga_py.pga.mvec3dp, 16, tuple(float(i + 1) for i in range(16))),
    (ga_py.pga.dualnum2dp, 2, (3.0, 4.0)),
    (ga_py.pga.dualnum3dp, 2, (5.0, 6.0)),
    # Inherited-shape PGA primitives — fields come from the base type.
    (ga_py.pga.point2d, 2, (7.0, 8.0)),
    (ga_py.pga.point3d, 3, (1.0, 2.0, 3.0)),
    (ga_py.pga.point2dp, 3, (1.0, 2.0, 3.0)),
    (ga_py.pga.point3dp, 4, (1.0, 2.0, 3.0, 4.0)),
    (ga_py.pga.line2d, 3, (0.1, 0.2, 0.3)),
    (ga_py.pga.line3d, 6, (1.0, 2.0, 3.0, 4.0, 5.0, 6.0)),
    (ga_py.pga.plane3d, 4, (1.0, 2.0, 3.0, 4.0)),
]


@pytest.mark.parametrize("cls,n,args", _VEC_TYPES,
                         ids=[c.__name__ for c, _, _ in _VEC_TYPES])
def test_array_shape_and_dtype(cls, n, args):
    obj = cls(*args)
    arr = np.array(obj)
    assert arr.shape == (n,)
    assert arr.dtype == np.float64


@pytest.mark.parametrize("cls,n,args", _VEC_TYPES,
                         ids=[c.__name__ for c, _, _ in _VEC_TYPES])
def test_array_values_match_ctor_args(cls, n, args):
    obj = cls(*args)
    arr = np.array(obj)
    np.testing.assert_array_equal(arr, np.asarray(args, dtype=np.float64))


@pytest.mark.parametrize("cls,n,args", _VEC_TYPES,
                         ids=[c.__name__ for c, _, _ in _VEC_TYPES])
def test_ndarray_ctor_round_trip(cls, n, args):
    obj1 = cls(*args)
    arr1 = np.array(obj1)
    obj2 = cls(arr1)
    arr2 = np.array(obj2)
    np.testing.assert_array_equal(arr1, arr2)


@pytest.mark.parametrize("cls,n,args", _VEC_TYPES,
                         ids=[c.__name__ for c, _, _ in _VEC_TYPES])
def test_array_returns_independent_copy(cls, n, args):
    """Mutating the returned ndarray must not affect the C++ object."""
    obj = cls(*args)
    arr = np.array(obj)
    arr[0] = -999.0
    arr2 = np.array(obj)
    assert arr2[0] == args[0]


def test_array_accepts_dtype_and_copy_kwargs():
    """NumPy 2.x calls __array__(dtype=..., copy=...). Must not raise."""
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)
    # Both kwargs explicitly None.
    a = v.__array__(dtype=None, copy=None)
    np.testing.assert_array_equal(a, [1.0, 2.0, 3.0])
    # NumPy 2.x dtype-cast path: numpy itself converts after we return float64.
    a32 = np.asarray(v, dtype=np.float32)
    assert a32.dtype == np.float32
    np.testing.assert_array_equal(a32, np.asarray([1.0, 2.0, 3.0],
                                                  dtype=np.float32))


def test_ndarray_ctor_rejects_wrong_shape():
    with pytest.raises(TypeError):
        ga_py.ega.vec3d(np.array([1.0, 2.0]))      # too short
    with pytest.raises(TypeError):
        ga_py.ega.vec3d(np.array([1.0, 2.0, 3.0, 4.0]))  # too long


def test_ndarray_ctor_accepts_contiguous_row():
    """A row of a 2D C-contiguous array is itself contiguous."""
    big = np.arange(12.0).reshape(4, 3)
    v = ga_py.ega.vec3d(big[1])
    assert (v.x, v.y, v.z) == (3.0, 4.0, 5.0)


def test_batch_via_list_comprehension():
    """The Tier-1-cheap idiom: build a 2D ndarray from a list of GA elements."""
    vs = [ga_py.ega.vec3d(i, i + 1, i + 2) for i in range(4)]
    batch = np.stack([np.array(v) for v in vs])
    assert batch.shape == (4, 3)
    np.testing.assert_array_equal(
        batch,
        np.array([[i, i + 1, i + 2] for i in range(4)], dtype=np.float64),
    )
