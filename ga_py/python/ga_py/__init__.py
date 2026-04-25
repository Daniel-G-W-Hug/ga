"""ga_py — Python bindings for the hd::ga geometric algebra library.

The compiled extension is exposed as `_ga_py`; this package re-exports its
submodules under more idiomatic names. Import as:

    import ga_py
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)

A small number of free functions in `ga/` are written as fully-generic C++
templates (e.g. `template<typename A, typename B> auto f(A&& a, B&& b)`),
which the binding generator cannot bind by enumerating overloads — there's
no overload list, just a body that compiles for whatever types it's called
with. We expose those functions here as Python wrappers that call the
underlying composition. They contain NO Python math — every call dispatches
into the bound C++ functions.
"""
import _ga_py
from _ga_py import ega, pga  # noqa: F401

# Forward the top-level (`hd::ga::*`) functions from the compiled module to
# this package, so callers can write `ga_py.deg2rad(180.0)` directly. This
# mirrors the C++ idiom of bringing `using namespace hd::ga;` into scope.
for _attr in dir(_ga_py):
    if _attr.startswith("_") or _attr in ("ega", "pga"):
        continue
    globals()[_attr] = getattr(_ga_py, _attr)
del _attr

__all__ = ["ega", "pga"] + [n for n in dir(_ga_py)
                            if not n.startswith("_") and n not in ("ega", "pga")]


# ---------------------------------------------------------------------------
# PGA 2dp — forwarder functions (generic C++ templates)
# Source: ga/ga_pga2dp_ops_products.hpp
# ---------------------------------------------------------------------------

def _pga2dp_l_bulk_contract2dp(a, b):
    """Left bulk contraction (PGA 2dp). Equals `rwdg(bulk_dual(a), b)`."""
    return pga.rwdg(pga.bulk_dual(a), b)


def _pga2dp_l_weight_contract2dp(a, b):
    """Left weight contraction (PGA 2dp). Equals `rwdg(weight_dual(a), b)`."""
    return pga.rwdg(pga.weight_dual(a), b)


def _pga2dp_r_bulk_contract2dp(a, b):
    """Right bulk contraction (PGA 2dp). Equals `rwdg(a, bulk_dual(b))`."""
    return pga.rwdg(a, pga.bulk_dual(b))


def _pga2dp_r_weight_contract2dp(a, b):
    """Right weight contraction (PGA 2dp). Equals `rwdg(a, weight_dual(b))`."""
    return pga.rwdg(a, pga.weight_dual(b))


def _pga2dp_l_bulk_expand2dp(a, b):
    """Left bulk expansion (PGA 2dp). Equals `wdg(bulk_dual(a), b)`."""
    return pga.wdg(pga.bulk_dual(a), b)


def _pga2dp_l_weight_expand2dp(a, b):
    """Left weight expansion (PGA 2dp). Equals `wdg(weight_dual(a), b)`."""
    return pga.wdg(pga.weight_dual(a), b)


def _pga2dp_r_bulk_expand2dp(a, b):
    """Right bulk expansion (PGA 2dp). Equals `wdg(a, bulk_dual(b))`."""
    return pga.wdg(a, pga.bulk_dual(b))


def _pga2dp_r_weight_expand2dp(a, b):
    """Right weight expansion (PGA 2dp). Equals `wdg(a, weight_dual(b))`."""
    return pga.wdg(a, pga.weight_dual(b))


# Inject into the pga submodule so users access them at the same path as the
# generated C++ bindings: `ga_py.pga.l_bulk_contract2dp(a, b)`.
pga.l_bulk_contract2dp = _pga2dp_l_bulk_contract2dp
pga.l_weight_contract2dp = _pga2dp_l_weight_contract2dp
pga.r_bulk_contract2dp = _pga2dp_r_bulk_contract2dp
pga.r_weight_contract2dp = _pga2dp_r_weight_contract2dp
pga.l_bulk_expand2dp = _pga2dp_l_bulk_expand2dp
pga.l_weight_expand2dp = _pga2dp_l_weight_expand2dp
pga.r_bulk_expand2dp = _pga2dp_r_bulk_expand2dp
pga.r_weight_expand2dp = _pga2dp_r_weight_expand2dp


# ---------------------------------------------------------------------------
# PGA 2dp — projection / distance forwarders (also generic C++ templates)
# Source: ga/ga_pga2dp_ops.hpp
# ---------------------------------------------------------------------------

# Tolerance for "is the weight squared norm essentially 1?" — matches the
# rough scale of `eps` used in the C++ unitize-or-not heuristic. Cheap
# guard; the C++ source uses a comparable epsilon.
_EPS = 1e-9


def _pga2dp_ortho_proj2dp(a, b):
    """Orthogonal projection of `a` onto the larger-grade `b` (PGA 2dp).

    Equals `unitize(rwdg(b, r_weight_expand2dp(a, b)))` when not at horizon.
    REQUIRES: gr(a) < gr(b).
    """
    p = pga.rwdg(b, pga.r_weight_expand2dp(a, b))
    nrm_sq = float(pga.weight_nrm_sq(p))
    if nrm_sq > _EPS and nrm_sq != 1.0:
        p = pga.unitize(p)
    return p


def _pga2dp_central_proj2dp(a, b):
    """Central projection of `a` onto the larger-grade `b` (PGA 2dp).

    Equals `unitize(rwdg(b, r_bulk_expand2dp(a, b)))` when not at horizon.
    REQUIRES: gr(a) < gr(b).
    """
    p = pga.rwdg(b, pga.r_bulk_expand2dp(a, b))
    nrm_sq = float(pga.weight_nrm_sq(p))
    if nrm_sq > _EPS and nrm_sq != 1.0:
        p = pga.unitize(p)
    return p


def _pga2dp_ortho_antiproj2dp(a, b):
    """Orthogonal anti-projection (PGA 2dp).

    Equals `wdg(b, r_weight_contract2dp(a, b))`.
    """
    return pga.wdg(b, pga.r_weight_contract2dp(a, b))


# Per-class grade lookup for the 2dp algebra. Mirrors the `gr()` overloads in
# ga/detail/type_t/ga_type2dp.hpp. Used by dist2dp to dispatch the C++
# `if constexpr` on grade sum at runtime.
_GRADE_2DP = {
    "scalar2dp": 0,
    "vec2dp": 1, "point2dp": 1, "point2d": 1, "vector2d": 1,
    "bivec2dp": 2, "line2d": 2,
    "pscalar2dp": 3,
}


def _grade_2dp(x):
    name = type(x).__name__
    g = _GRADE_2DP.get(name)
    if g is None:
        raise TypeError(
            f"dist2dp: cannot determine 2dp grade of {name!r}; "
            "expected scalar2dp / vec2dp / bivec2dp / pscalar2dp or a derived "
            "geometric primitive (point*, vector*, line*)"
        )
    return g


def _pga2dp_dist2dp(a, b):
    """Euclidean distance between two PGA 2dp objects, returned as a
    DualNum2dp `(homogeneous_magnitude, weight)`. Source: ga_pga2dp_ops.hpp.

    Mirrors the C++ `if constexpr (gr(a)+gr(b) == 3)` dispatch at runtime via
    a per-class grade lookup.
    """
    if _grade_2dp(a) + _grade_2dp(b) == 3:
        c0 = float(pga.rwdg(a, b))
        c1 = float(pga.weight_nrm(pga.wdg(a, pga.att(b))))
    else:
        c0 = float(pga.bulk_nrm(pga.att(pga.wdg(a, b))))
        c1 = float(pga.weight_nrm(pga.wdg(a, pga.att(b))))
    return pga.dualnum2dp(c0, c1)


pga.ortho_proj2dp = _pga2dp_ortho_proj2dp
pga.central_proj2dp = _pga2dp_central_proj2dp
pga.ortho_antiproj2dp = _pga2dp_ortho_antiproj2dp
pga.dist2dp = _pga2dp_dist2dp


# ---------------------------------------------------------------------------
# PGA 3dp — forwarder functions (generic C++ templates)
# Source: ga/ga_pga3dp_ops_products.hpp
#
# 3dp is even-dimensional, so the bulk/weight duals come in left and right
# variants (l_bulk_dual / r_bulk_dual etc.). The contractions and expansions
# use the matching side.
# ---------------------------------------------------------------------------

def _pga3dp_l_bulk_contract3dp(a, b):
    """Left bulk contraction (PGA 3dp). = `rwdg(l_bulk_dual(a), b)`."""
    return pga.rwdg(pga.l_bulk_dual(a), b)


def _pga3dp_l_weight_contract3dp(a, b):
    """Left weight contraction (PGA 3dp). = `rwdg(l_weight_dual(a), b)`."""
    return pga.rwdg(pga.l_weight_dual(a), b)


def _pga3dp_r_bulk_contract3dp(a, b):
    """Right bulk contraction (PGA 3dp). = `rwdg(a, r_bulk_dual(b))`."""
    return pga.rwdg(a, pga.r_bulk_dual(b))


def _pga3dp_r_weight_contract3dp(a, b):
    """Right weight contraction (PGA 3dp). = `rwdg(a, r_weight_dual(b))`."""
    return pga.rwdg(a, pga.r_weight_dual(b))


def _pga3dp_l_bulk_expand3dp(a, b):
    """Left bulk expansion (PGA 3dp). = `wdg(l_bulk_dual(a), b)`."""
    return pga.wdg(pga.l_bulk_dual(a), b)


def _pga3dp_l_weight_expand3dp(a, b):
    """Left weight expansion (PGA 3dp). = `wdg(l_weight_dual(a), b)`."""
    return pga.wdg(pga.l_weight_dual(a), b)


def _pga3dp_r_bulk_expand3dp(a, b):
    """Right bulk expansion (PGA 3dp). = `wdg(a, r_bulk_dual(b))`."""
    return pga.wdg(a, pga.r_bulk_dual(b))


def _pga3dp_r_weight_expand3dp(a, b):
    """Right weight expansion (PGA 3dp). = `wdg(a, r_weight_dual(b))`."""
    return pga.wdg(a, pga.r_weight_dual(b))


pga.l_bulk_contract3dp = _pga3dp_l_bulk_contract3dp
pga.l_weight_contract3dp = _pga3dp_l_weight_contract3dp
pga.r_bulk_contract3dp = _pga3dp_r_bulk_contract3dp
pga.r_weight_contract3dp = _pga3dp_r_weight_contract3dp
pga.l_bulk_expand3dp = _pga3dp_l_bulk_expand3dp
pga.l_weight_expand3dp = _pga3dp_l_weight_expand3dp
pga.r_bulk_expand3dp = _pga3dp_r_bulk_expand3dp
pga.r_weight_expand3dp = _pga3dp_r_weight_expand3dp


# ---------------------------------------------------------------------------
# PGA 3dp — projection / distance forwarders (also generic C++ templates)
# Source: ga/ga_pga3dp_ops.hpp
# ---------------------------------------------------------------------------

def _pga3dp_ortho_proj3dp(a, b):
    """Orthogonal projection of `a` onto larger-grade `b` (PGA 3dp).

    Equals `unitize(rwdg(b, r_weight_expand3dp(a, b)))` when not at horizon.
    REQUIRES: gr(a) < gr(b).
    """
    p = pga.rwdg(b, pga.r_weight_expand3dp(a, b))
    nrm_sq = float(pga.weight_nrm_sq(p))
    if nrm_sq > _EPS and nrm_sq != 1.0:
        p = pga.unitize(p)
    return p


def _pga3dp_central_proj3dp(a, b):
    """Central projection of `a` onto larger-grade `b` (PGA 3dp).

    Equals `unitize(rwdg(b, r_bulk_expand3dp(a, b)))` when not at horizon.
    REQUIRES: gr(a) < gr(b).
    """
    p = pga.rwdg(b, pga.r_bulk_expand3dp(a, b))
    nrm_sq = float(pga.weight_nrm_sq(p))
    if nrm_sq > _EPS and nrm_sq != 1.0:
        p = pga.unitize(p)
    return p


def _pga3dp_ortho_antiproj3dp(a, b):
    """Orthogonal anti-projection (PGA 3dp).

    Equals `wdg(b, r_weight_contract3dp(a, b))`.
    """
    return pga.wdg(b, pga.r_weight_contract3dp(a, b))


# Per-class grade lookup for the 3dp algebra. Mirrors `gr()` overloads in
# ga/detail/type_t/ga_type3dp.hpp. Used by dist3dp to dispatch the C++
# `if constexpr (gr(a)+gr(b) == 4)` at runtime.
_GRADE_3DP = {
    "scalar3dp": 0,
    "vec3dp": 1, "point3dp": 1, "point3d": 1, "vector3d": 1,
    "bivec3dp": 2, "line3d": 2,
    "trivec3dp": 3, "plane3d": 3,
    "pscalar3dp": 4,
}


def _grade_3dp(x):
    name = type(x).__name__
    g = _GRADE_3DP.get(name)
    if g is None:
        raise TypeError(
            f"dist3dp: cannot determine 3dp grade of {name!r}; "
            "expected scalar3dp / vec3dp / bivec3dp / trivec3dp / pscalar3dp "
            "or a derived geometric primitive (point*, vector*, line*, plane3d)"
        )
    return g


def _pga3dp_dist3dp(a, b):
    """Euclidean distance between two PGA 3dp objects, returned as a
    DualNum3dp `(homogeneous_magnitude, weight)`. Source: ga_pga3dp_ops.hpp.

    Mirrors the C++ `if constexpr (gr(a)+gr(b) == 4)` dispatch at runtime via
    a per-class grade lookup. (In 3dp, the meaningful incidence threshold is
    grade-sum 4 because the pseudoscalar is grade 4, not 3 as in 2dp.)
    """
    if _grade_3dp(a) + _grade_3dp(b) == 4:
        c0 = float(pga.rwdg(a, b))
        c1 = float(pga.weight_nrm(pga.wdg(a, pga.att(b))))
    else:
        c0 = float(pga.bulk_nrm(pga.att(pga.wdg(a, b))))
        c1 = float(pga.weight_nrm(pga.wdg(a, pga.att(b))))
    return pga.dualnum3dp(c0, c1)


pga.ortho_proj3dp = _pga3dp_ortho_proj3dp
pga.central_proj3dp = _pga3dp_central_proj3dp
pga.ortho_antiproj3dp = _pga3dp_ortho_antiproj3dp
pga.dist3dp = _pga3dp_dist3dp
