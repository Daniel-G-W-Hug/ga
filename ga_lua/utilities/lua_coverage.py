#!/usr/bin/env python3
"""
lua_coverage.py

Compare the GA surface exposed by the Python wrapper (`ga_py`) against what the
Lua wrapper (`ga_lua`) binds, and report the per-algebra gaps. The Python binding
is the parity target: it is the manifest-driven, POD-filtered realisation of the
public API (stateful systems and internal template helpers are already excluded
there), so "comparable in completeness to ga_py" is exactly "no items reported
missing here".

What it compares (per algebra ega / pga / cga / sta, plus top-level free functions):

  - types       : ga_py classes            vs  ga_lua `new_usertype<NAME>`
  - functions   : ga_py module callables   vs  ga_lua `set_function("NAME", ...)`
  - constants   : ga_py module data attrs  vs  ga_lua `lua["NAME"] = ...`
  - operators   : ga_py `a + b` / `a - b`  vs  ga_lua `sol::resolve<R(A, B)>(operator+)`
                  per ORDERED type pair (see the dispatch note below)

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
  - Operators are checked PER ORDERED PAIR, against Lua's actual dispatch rule,
    which has two steps and a base-class wrinkle. For `a + b` Lua takes `__add`
    from a's metatable; ONLY IF a has none does it fall back to b's. So:

      * if the left type binds any operator+, the pair must match one of ITS
        resolves -- a `sol::resolve<R(A, B)>(operator+)` sitting in B's block is
        then dead code, since Lua never gets there. Those are reported as
        UNREACHABLE, because they compile and read correctly.
      * if the left type binds none, the right type's set is what decides.

    `sol::base_classes` widens this: a resolve taking a base accepts a derived
    argument, so point3dp + vec3dp is reachable through vec3dp's own set even
    though point3dp binds no operator at all. The reverse of that coin is that
    two DERIVED operands with no operator between them fail outright
    ("attempt to perform arithmetic"), which no amount of base_classes fixes.

    Checking types alone missed a 340-of-428 gap in the operator surface until
    2026-08-30, which is why this kind exists.

Usage:
    ./lua_coverage.py                       # full per-algebra report
    ./lua_coverage.py --algebra=sta         # one algebra (ega|pga|cga|sta|top)
    ./lua_coverage.py --algebra=ega,pga     # comma-separated
    ./lua_coverage.py --kind=types          # types|functions|constants|operators
    ./lua_coverage.py --summary             # totals only, no item lists
    ./lua_coverage.py --present             # also list what IS already bound

Requires the ga_py extension to be built (`-D_GA_BUILD_PYTHON=ON`, then
`cmake --build build --target _ga_py`). Run with the ga_py venv so `ga_py` imports:

    ga_py/.venv/bin/python ga_lua/utilities/lua_coverage.py
"""

import argparse
import inspect
import operator as _op
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
GA_LUA_HPP = REPO_ROOT / "ga_lua" / "src" / "ga_lua.hpp"
ALGEBRAS = ["ega", "pga", "cga", "sta", "top"]
KINDS = ["types", "functions", "constants", "operators"]


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
    """Return {algebra: {types, functions, constants}} for ega/pga/cga/sta/top."""
    result = {}
    for alg in ("ega", "pga", "cga", "sta"):
        sub = getattr(ga_py, alg, None)
        if sub is None:
            result[alg] = {"types": set(), "functions": set(), "constants": set()}
            continue
        result[alg] = _classify_module(sub)

    # top-level: free functions / constants that live directly on ga_py
    # (exclude the ega/pga/cga/sta submodules themselves, handled above).
    top = _classify_module(ga_py)
    result["top"] = top
    return result


def collect_ga_lua():
    """Parse ga_lua.hpp into bound types / functions / constants.

    Functions are collected PER ALGEBRA TABLE, because that is how they are bound:
    `ega.set_function("wdg", ...)`, `pga.rgpr`, `cga.dot`, `sta.wdg`. Matching the
    name alone would report a function bound in ONE table as present for every
    algebra --- which was harmless while a single global name served all of them,
    and became a blind spot the moment the tables were introduced. Types and
    constants stay global, so those remain flat sets.
    """
    src = GA_LUA_HPP.read_text()
    # functions reach Lua two ways: C++ `<table>.set_function("name", ...)` (the
    # name may sit on the next line) and Lua-level `function <table>.name(...)`
    # definitions inside the embedded prelude (register_forwarders -> lua.script).
    funcs = {a: set() for a in ("ega", "pga", "cga", "sta", "top")}
    for tbl, name in re.findall(
            r'\b(ega|pga|cga|sta|lua)\.set_function\(\s*"([A-Za-z0-9_]+)"', src):
        funcs["top" if tbl == "lua" else tbl].add(name)
    for tbl, name in re.findall(
            r"^\s*function (?:(ega|pga|cga|sta)\.)?([A-Za-z0-9_]+)\(", src, re.MULTILINE):
        funcs[tbl or "top"].add(name)
    # ga_py exposes both value/POD types and scoped enums as classes; ga_lua binds
    # the former with new_usertype<> and the latter with new_enum<> -> count both.
    types = set(re.findall(r"new_usertype<\s*([A-Za-z0-9_]+)\s*>", src))
    types |= set(re.findall(r'new_enum<\s*([A-Za-z0-9_]+)\s*>', src))
    return {
        "types": types,
        "functions": funcs,
        "constants": set(re.findall(r'lua\[\s*"([A-Za-z0-9_]+)"\s*\]\s*=', src)),
    }


def collect_ga_py_operators(ga_py):
    """{algebra: {"A + B", "A - B", ...}} -- the pairs ga_py actually accepts.

    Determined by TRYING each ordered pair rather than from a table of type
    names, so cross-family pairs (a 2d type plus a 3d one, which share the `ega`
    submodule but no algebra) drop out on their own, and a newly bound type is
    picked up without editing this script.
    """
    out = {a: set() for a in ALGEBRAS}
    for alg in ("ega", "pga", "cga", "sta"):
        sub_mod = getattr(ga_py, alg, None)
        if sub_mod is None:
            continue
        inst = {}
        for name in dir(sub_mod):
            if name.startswith("_"):
                continue
            obj = getattr(sub_mod, name)
            if not inspect.isclass(obj):
                continue
            try:                     # value types are default-constructible;
                inst[name] = obj()   # enums and PODs that are not simply drop out
            except Exception:
                pass
        for a, av in inst.items():
            for b, bv in inst.items():
                for sym, fn in (("+", _op.add), ("-", _op.sub)):
                    try:
                        fn(av, bv)
                    except TypeError:
                        continue
                    out[alg].add(f"{a} {sym} {b}")
    return out


def _usertype_blocks(src):
    """Yield (type_name, block_text) for each `lua.new_usertype<NAME>(...)`.

    Paren matching ONLY: `<` and `>` cannot be treated as brackets here because
    the blocks bind `operator<<` / `operator>>`, whose angle brackets never close.
    Counting them silently swallows the rest of the file.
    """
    for m in re.finditer(r"lua\.new_usertype<\s*([A-Za-z0-9_]+)\s*>\(", src):
        i = m.end() - 1
        depth = 0
        j = i
        while j < len(src):
            c = src[j]
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        yield m.group(1), src[i:j + 1]


def _parse_lua_operator_tables(src):
    """-> (resolves, bases) where

    resolves[owner][sym] = {(lhs_type, rhs_type), ...} declared in owner's block
    bases[t]             = the type t derives from via sol::base_classes, or None
    """
    resolves, bases = {}, {}
    rx = re.compile(
        r"sol::resolve<\s*[A-Za-z0-9_]+\s*\(([^)]*)\)\s*>\(operator([+-])\)")
    rb = re.compile(r"sol::bases<\s*([A-Za-z0-9_]+)\s*>")
    for owner, blk in _usertype_blocks(src):
        tbl = resolves.setdefault(owner, {"+": set(), "-": set()})
        for args, sym in rx.findall(blk):
            parts = [re.sub(r"\bconst\b|&|\s", "", a) for a in args.split(",")]
            if len(parts) == 2:
                tbl[sym].add((parts[0], parts[1]))
        m = rb.search(blk)
        bases[owner] = m.group(1) if m else None
    return resolves, bases


def collect_ga_lua_operators(src):
    """(bound, unreachable) sets of "A + B" strings, per Lua's dispatch rule.

    A pair is BOUND when the metatable Lua actually consults carries a resolve
    that accepts both operands (a base-class parameter accepts a derived
    argument). A resolve is UNREACHABLE when it names a left operand whose own
    type binds that operator elsewhere -- Lua stops at the left metatable, so the
    misplaced copy can never fire.
    """
    resolves, bases = _parse_lua_operator_tables(src)

    def anc(t):
        out, cur = [], t
        while cur is not None:
            out.append(cur)
            cur = bases.get(cur)
        return out

    def matches(owner, a, b, sym):
        la, lb = anc(a), anc(b)
        return any(p in la and q in lb for p, q in resolves.get(owner, {}).get(sym, ()))

    bound, unreachable = set(), set()
    for owner, tbl in resolves.items():
        for sym, pairs in tbl.items():
            for lhs, rhs in pairs:
                item = f"{lhs} {sym} {rhs}"
                if lhs in anc(owner):
                    # owner's own type, or a base of it: sol converts the derived
                    # argument, so Lua reaches this through owner's metatable
                    bound.add(item)
                elif resolves.get(lhs, {}).get(sym):
                    # lhs is an unrelated type that binds this operator itself, so
                    # Lua stops at ITS metatable and never reaches this copy
                    unreachable.add(item)
                else:
                    # lhs binds nothing -> Lua falls back to this (the rhs) table
                    if matches(owner, lhs, rhs, sym):
                        bound.add(item)

    # a derived type that binds no operator of its own still reaches its base's
    # set through the right-operand fallback -- enumerate those too
    for a in resolves:
        for b in resolves:
            for sym in ("+", "-"):
                item = f"{a} {sym} {b}"
                if item in bound:
                    continue
                owner = a if resolves.get(a, {}).get(sym) else b
                if matches(owner, a, b, sym):
                    bound.add(item)
    return bound, unreachable


def main():
    ap = argparse.ArgumentParser(
        description="Report ga_lua coverage gaps relative to the ga_py binding.")
    ap.add_argument("--algebra", default="all",
                    help="comma list of ega,pga,cga,sta,top (default: all)")
    ap.add_argument("--kind", default="all",
                    help="comma list of types,functions,constants,operators "
                         "(default: all)")
    ap.add_argument("--summary", action="store_true",
                    help="print totals only, no per-item lists")
    ap.add_argument("--present", action="store_true",
                    help="also list items already bound in ga_lua")
    args = ap.parse_args()

    algebras = ALGEBRAS if args.algebra == "all" else args.algebra.split(",")
    kinds = KINDS if args.kind == "all" else args.kind.split(",")
    for k in kinds:
        if k not in KINDS:
            sys.exit(f"unknown kind '{k}' (choose from {KINDS})")

    ga_py = _import_ga_py()
    py = collect_ga_py(ga_py)
    lua = collect_ga_lua()

    if "operators" in kinds:
        py_ops = collect_ga_py_operators(ga_py)
        lua_ops, lua_unreachable = collect_ga_lua_operators(GA_LUA_HPP.read_text())
        for alg in py:
            py[alg]["operators"] = py_ops.get(alg, set())

    grand = {"missing": 0, "present": 0, "target": 0}
    for alg in algebras:
        if alg not in py:
            print(f"unknown algebra '{alg}' (choose from {ALGEBRAS})")
            continue
        print(f"\n=== {alg.upper()} "
              "==============================================================")
        for kind in kinds:
            target = py[alg][kind]
            # functions are per-algebra tables; types/constants stay global. A
            # top-level helper (deg2rad, sign, ...) counts for any algebra asking.
            if kind == "operators":
                bound = lua_ops
            elif kind == "functions" and alg == "top":
                # ga_py keeps a few algebra-agnostic helpers (gr, rgr) at the top
                # level; ga_lua binds them per algebra, where they dispatch on the
                # argument type just the same. Reachable, so not a gap --- count a
                # top-level target found in any table.
                bound = set().union(*lua["functions"].values())
            elif kind == "functions":
                bound = lua["functions"][alg] | lua["functions"]["top"]
            else:
                bound = lua[kind]
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

    if "operators" in kinds and lua_unreachable:
        print("\n=== UNREACHABLE "
              "=============================================================")
        print("  resolves bound in the WRONG usertype block: Lua takes __add/__sub")
        print("  from the LEFT operand's metatable, so these can never fire.")
        for item in sorted(lua_unreachable):
            print(f"       ! {item}")

    print("\n=== TOTAL "
          "===================================================================")
    print(f"  target {grand['target']} | bound {grand['present']} | "
          f"MISSING {grand['missing']}")


if __name__ == "__main__":
    main()
