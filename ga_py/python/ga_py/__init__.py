"""ga_py — Python bindings for the hd::ga geometric algebra library.

The compiled extension is exposed as `_ga_py`; this package re-exports its
submodules under more idiomatic names. Import as:

    import ga_py
    v = ga_py.ega.vec3d(1.0, 2.0, 3.0)

A small number of free functions in `ga/` are written as fully-generic C++
templates (e.g. `template<typename A, typename B> auto f(A&& a, B&& b)`),
which the binding generator cannot bind by enumerating overloads — there's
no overload list, just a body that compiles for whatever types it's called
with. The ones below are exposed here as Python wrappers, and each is a
LITERAL one-line composition of already-bound primitives, identical to the
C++ body: `l_bulk_contract2dp(a, b)` is `rwdg(bulk_dual(a), b)` in both. They
contain NO Python math, and because the composition is what dispatches, they
cover every grade pair the bound `rwdg` / `wdg` / `dual` accept.

The generic templates that are NOT one-line compositions — the projection and
antiprojection families, `dist*dp`, `try_unitize`, `is_simple`'s defaulted
tolerance — are bound in C++ instead (`ga_py/src/bindings_projections.cpp`).
They carry thresholds, a blade-target guard and an `if constexpr` dispatch,
which a Python wrapper can only reproduce by copying constants; when it did,
it drifted. Do not move them back here, and do not move these ones there:
enumerating type pairs in C++ would narrow the coverage they have now.
"""
import _ga_py
from _ga_py import ega, pga, cga, sta  # noqa: F401

# Forward the top-level (`hd::ga::*`) functions from the compiled module to
# this package, so callers can write `ga_py.deg2rad(180.0)` directly. This
# mirrors the C++ idiom of bringing `using namespace hd::ga;` into scope.
for _attr in dir(_ga_py):
    if _attr.startswith("_") or _attr in ("ega", "pga", "cga", "sta"):
        continue
    globals()[_attr] = getattr(_ga_py, _attr)
del _attr

__all__ = ["ega", "pga", "cga", "sta"] + [n for n in dir(_ga_py)
                                          if not n.startswith("_")
                                          and n not in ("ega", "pga", "cga", "sta")]


# ---------------------------------------------------------------------------
# PGA 2dp — forwarder functions (generic C++ templates)
# Source: ga/pga/ga_pga2dp_ops_products.hpp
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
# PGA 3dp — forwarder functions (generic C++ templates)
# Source: ga/pga/ga_pga3dp_ops_products.hpp
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
# STA4D — expansion forwarders (generic C++ templates)
# Source: ga/sta/ga_sta4ds_ops_products.hpp
# ---------------------------------------------------------------------------

def _sta_l_expand4ds(a, b):
    """Left expansion (STA4D). Equals `wdg(l_dual(a), b)`."""
    return sta.wdg(sta.l_dual(a), b)


def _sta_r_expand4ds(a, b):
    """Right expansion (STA4D). Equals `wdg(a, r_dual(b))`."""
    return sta.wdg(a, sta.r_dual(b))


sta.l_expand4ds = _sta_l_expand4ds
sta.r_expand4ds = _sta_r_expand4ds


# ---------------------------------------------------------------------------
# Pure-Python reconstruction of the stateful multibody system classes.
# The C++ static_/kinematic_/dynamic_system3dp classes are not bound (by
# design); `ga_py.systems` rebuilds them from the primitives above. Imported
# last, since it depends on the `pga` submodule being fully set up.
# ---------------------------------------------------------------------------
from . import systems  # noqa: E402,F401

__all__.append("systems")
