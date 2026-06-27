#!/usr/bin/env python3
"""
lua_coverage.py

Compare the GA surface exposed by the Python wrapper (`ga_py`) against what the
Lua wrapper (`ga_lua`) binds, and report the per-algebra gaps. The Python binding
is the parity target: it is the manifest-driven, POD-filtered realisation of the
public API (stateful systems and internal template helpers are already excluded
there), so "comparable in completeness to ga_py" is exactly "no items reported
missing here".

What it compares (per algebra ega / pga / sta, plus top-level free functions):

  - types       : ga_py classes            vs  ga_lua `new_usertype<NAME>`
  - functions   : ga_py module callables   vs  ga_lua `set_function("NAME", ...)`
  - constants   : ga_py module data attrs  vs  ga_lua `lua["NAME"] = ...`

Source of truth:
  - ga_py  : introspected from the *built* extension (so the build must be current;
             this auto-excludes templates / stateful systems that ga_py never binds).
  - ga_lua : parsed from `ga_lua/src/ga_lua.hpp` by simple, robust regexes.

Granularity & caveats:
  - Function and constant matching is by NAME (ga_lua binds one global `nrm` with
    per-type `sol::overload`s; ga_py exposes `ega.nrm`, `pga.nrm`, ...). So a name
    reported as "present" is bound *somewhere* in ga_lua -- this tool does NOT verify
    that every per-type overload is present. The dominant gaps (a whole algebra, e.g.
    STA, or an unbound function) show up clearly via the type/constant lists and the
    missing-name lists. For per-overload auditing, inspect the relevant
    `register_*` block by hand.
  - Operators are intentionally out of scope: ga_py exposes them as dunder methods on
    classes, ga_lua as `sol::meta_function` entries inside each `new_usertype`. They
    are bound per type, not as free functions, so they are tracked when a type is
    bound, not listed here.

Usage:
    ./lua_coverage.py                       # full per-algebra report
    ./lua_coverage.py --algebra=sta         # one algebra (ega|pga|sta|top)
    ./lua_coverage.py --algebra=ega,pga     # comma-separated
    ./lua_coverage.py --kind=types          # only types (types|functions|constants)
    ./lua_coverage.py --summary             # totals only, no item lists
    ./lua_coverage.py --present             # also list what IS already bound

Requires the ga_py extension to be built (`-D_GA_BUILD_PYTHON=ON`, then
`cmake --build build --target _ga_py`). Run with the ga_py venv so `ga_py` imports:

    ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py
"""

import argparse
import inspect
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
GA_LUA_HPP = REPO_ROOT / "ga_lua" / "src" / "ga_lua.hpp"
ALGEBRAS = ["ega", "pga", "sta", "top"]


def _import_ga_py():
    """Import ga_py from the in-tree build, returning the module or exiting."""
    sys.path.insert(0, str(REPO_ROOT / "build" / "ga_py"))
    sys.path.insert(0, str(REPO_ROOT / "ga_py" / "python"))
    try:
        import ga_py  # noqa: E402
    except ImportError as exc:
        sys.exit(
            f"cannot import ga_py ({exc}).\n"
            "Build it first:\n"
            "  cmake -S . -B build -D_GA_BUILD_PYTHON=ON\n"
            "  cmake --build build --target _ga_py\n"
            "and run this script with the ga_py venv "
            "(ga_py/.venv/bin/python).")
    return ga_py


def _classify_module(mod):
    """Split a module's public attributes into (types, functions, constants)."""
    import types as _t

    out = {"types": set(), "functions": set(), "constants": set()}
    for name in dir(mod):
        if name.startswith("_"):
            continue
        obj = getattr(mod, name)
        if isinstance(obj, _t.ModuleType):
            continue
        if inspect.isclass(obj):
            out["types"].add(name)
        elif callable(obj):
            out["functions"].add(name)
        else:
            out["constants"].add(name)
    return out


def collect_ga_py(ga_py):
    """Return {algebra: {types, functions, constants}} for ega/pga/sta/top."""
    result = {}
    for alg in ("ega", "pga", "sta"):
        sub = getattr(ga_py, alg, None)
        if sub is None:
            result[alg] = {"types": set(), "functions": set(), "constants": set()}
            continue
        result[alg] = _classify_module(sub)

    # top-level: free functions / constants that live directly on ga_py
    # (exclude the ega/pga/sta submodules themselves, handled above).
    top = _classify_module(ga_py)
    result["top"] = top
    return result


def collect_ga_lua():
    """Parse ga_lua.hpp into flat sets of bound types / functions / constants."""
    src = GA_LUA_HPP.read_text()
    return {
        "types": set(re.findall(r"new_usertype<([A-Za-z0-9_]+)>", src)),
        "functions": set(re.findall(r'set_function\("([A-Za-z0-9_]+)"', src)),
        "constants": set(re.findall(r'lua\["([A-Za-z0-9_]+)"\]\s*=', src)),
    }


def main():
    ap = argparse.ArgumentParser(
        description="Report ga_lua coverage gaps relative to the ga_py binding.")
    ap.add_argument("--algebra", default="all",
                    help="comma list of ega,pga,sta,top (default: all)")
    ap.add_argument("--kind", default="all",
                    help="comma list of types,functions,constants (default: all)")
    ap.add_argument("--summary", action="store_true",
                    help="print totals only, no per-item lists")
    ap.add_argument("--present", action="store_true",
                    help="also list items already bound in ga_lua")
    args = ap.parse_args()

    algebras = ALGEBRAS if args.algebra == "all" else args.algebra.split(",")
    kinds = (["types", "functions", "constants"]
             if args.kind == "all" else args.kind.split(","))

    ga_py = _import_ga_py()
    py = collect_ga_py(ga_py)
    lua = collect_ga_lua()

    grand = {"missing": 0, "present": 0, "target": 0}
    for alg in algebras:
        if alg not in py:
            print(f"unknown algebra '{alg}' (choose from {ALGEBRAS})")
            continue
        print(f"\n=== {alg.upper()} "
              "==============================================================")
        for kind in kinds:
            target = py[alg][kind]
            bound = lua[kind]  # ga_lua is flat (global) -> name-level membership
            present = sorted(n for n in target if n in bound)
            missing = sorted(n for n in target if n not in bound)
            grand["target"] += len(target)
            grand["present"] += len(present)
            grand["missing"] += len(missing)
            print(f"  {kind:<10} target {len(target):>3} | "
                  f"bound {len(present):>3} | MISSING {len(missing):>3}")
            if not args.summary and missing:
                for n in missing:
                    print(f"       - {n}")
            if args.present and present:
                for n in present:
                    print(f"       + {n}")

    print("\n=== TOTAL "
          "===================================================================")
    print(f"  target {grand['target']} | bound {grand['present']} | "
          f"MISSING {grand['missing']}")


if __name__ == "__main__":
    main()
