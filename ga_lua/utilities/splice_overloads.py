#!/usr/bin/env python3
"""
splice_overloads.py

*Apply* counterpart to gen_lua_overloads.py: insert one algebra family's overloads
into that family's EXISTING `<table>.set_function("name", sol::overload(...))` block
in `ga_lua/src/ga_lua.hpp`, in place. Use it when a name is already bound for a family
and newly added library overloads have to join the same block.

A name is bound PER ALGEBRA TABLE -- `ega.wdg`, `pga.wdg`, `cga.wdg` and `sta.wdg` are
four separate bindings of one name, and 53 names are registered on more than one table.
So the table is part of the lookup key (FAMILY_TABLE below), and `--algebra` selects
both the manifest namespace to read AND the table to edit. Matching on the name alone
finds whichever block comes first in the file -- ega's -- and would compare, and then
splice, one family's overloads against another's.

(gen_lua_overloads.py drafts whole new blocks; this one extends existing ones. Same
find-vs-apply split as ga_prdxpr's validation_utilities vs utilities.)

For each requested function it:
  - reads the family's overloads from `ga_bindgen/manifest.json` (the namespace
    `hd::ga::<family>`), concretised to Lua aliases (shared with gen_lua_overloads.py);
  - finds the function's `set_function(...)` block, paren-matches the enclosing
    `sol::overload(`, and inserts the new `sol::resolve<...>(name)` lines before its
    close, under a `// <family>` comment;
  - is idempotent: a signature already present in the block is skipped (so re-running
    adds nothing), and a block with no `sol::overload` wrapper (a single resolve) is
    wrapped first. The presence test ignores whitespace, so it still recognises a
    signature clang-format has WRAPPED over several lines -- the long ones (the
    products, the `std::vector` batch overloads, any 7-argument constructor) always
    are, and a line-based test silently appends a duplicate for each of them.

It does NOT create new blocks (use gen_lua_overloads.py + paste for a brand-new name),
and it skips any requested name that is not already bound.

Usage:
    # auto: every already-bound name that has overloads in the family
    ./splice_overloads.py --algebra=sta
    # explicit subset
    ./splice_overloads.py --algebra=sta nrm,dot,inv,exp,log,sqrt
    ./splice_overloads.py --algebra=sta --dry-run        # report, write nothing

After splicing: clang-format ga_lua.hpp, build ga_lua, run lua_coverage.py and a smoke
script. Review the diff -- decide per overload whether it belongs in the REPL.
"""

import argparse
import importlib.util
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
HPP = REPO_ROOT / "ga_lua" / "src" / "ga_lua.hpp"
MANIFEST = REPO_ROOT / "ga_bindgen" / "manifest.json"

# reuse the concretiser / namespace map from the generator
_spec = importlib.util.spec_from_file_location("gen", HERE / "gen_lua_overloads.py")
gen = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gen)


# A concretised type the generator could not resolve leaves a bare `std` behind
# (`concretize("std::pair<Vec3dc<T>, Vec3dc<T>>")` -> `"std"`; only std::vector is
# handled). Such a signature can never match the file and would not compile if
# pasted, so it must be dropped rather than reported as a missing overload.
# std::pair returns are bound by HAND in ga_lua anyway -- sol2 unpacks a pair
# into two Lua return values only through a lambda, which is how dipole_points
# and the two_link_ik* helpers are written.
_UNRESOLVED = re.compile(r"\bstd\b(?!::)")


def family_resolves(manifest, family, name):
    """Deduped `sol::resolve<Ret(Args)>(name)` strings for one name in one family.

    Returns (resolves, n_skipped): signatures the concretiser could not express
    are skipped and counted, never emitted.
    """
    ns = gen.ALGEBRA_NS[family]
    seen, out, skipped = set(), [], 0
    for f in manifest["functions"]:
        if f["name"] != name or f["namespace"] != ns:
            continue
        for ov in f["overloads"]:
            ret = gen.concretize(ov["return_type"])
            sig = f'{ret}({", ".join(gen.param(p) for p in ov["param_types"])})'
            if sig in seen:
                continue
            seen.add(sig)
            if _UNRESOLVED.search(sig):
                skipped += 1
                continue
            out.append((sig, f"sol::resolve<{sig}>({name})"))
    return out, skipped


def _squeeze(s):
    """Drop every whitespace character, for wrap-insensitive comparison.

    The presence test below looks for a generated one-line signature inside a
    block that clang-format has already reflowed, so it MUST NOT depend on
    where the line breaks fell. A long signature is wrapped at whatever column
    the formatter chose:

        sol::resolve<trivec3dc(value_t, value_t, value_t, value_t, value_t,
                               value_t, value_t)>(circle3dc)

    which contains none of the one-line form as a substring. Collapsing runs of
    whitespace to a single space is not enough either -- a break placed right
    after '(' would leave 'trivec3dc( value_t' against 'trivec3dc(value_t'. So
    remove whitespace entirely, on both sides. Nothing in these signatures
    needs a space to stay well-formed once both sides are treated the same
    (e.g. 'vec3d const&' and 'vec3d const &' both become 'vec3dconst&').

    This mattered: without it the tool silently appended a SECOND, identical
    sol::resolve for circle3dc and dipole3dc, whose 7-argument signatures the
    formatter had split over two lines. It compiles -- sol::overload accepts
    duplicates -- so nothing downstream reports it. The products of each
    algebra carry the longest signatures in the file and wrap the same way.
    """
    return "".join(s.split())


def _overload_close(text, start):
    """Index of the ')' closing the first 'sol::overload(' at/after start, or None."""
    key = "sol::overload("
    pos = text.find(key, start)
    if pos == -1:
        return None
    i = pos + len(key)
    depth = 1
    while depth and i < len(text):
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return i - 1


# Which Lua table each family is registered on. ga_lua binds a name PER ALGEBRA
# TABLE -- ega.wdg, pga.wdg, cga.wdg and sta.wdg are four separate bindings of
# the same name -- so a block must be located by table as well as by name.
FAMILY_TABLE = {"ega": "ega", "pga": "pga", "cga": "cga", "sta": "sta", "top": "lua"}


def _block_span(text, table, name):
    """(start, end) of `<table>.set_function("name", ...);`, or None.

    The TABLE is part of the search key, not decoration. 53 names in ga_lua.hpp
    are registered on more than one table, and matching on the name alone finds
    whichever comes first in the file -- for `wdg` that is ega's block. The tool
    would then compare one family's overloads against another family's block,
    report every one of them as missing (145 of 145 for cga `wdg`), and, if run,
    splice cga overloads into the EGA table. That is wrong output, not merely a
    duplicate: the docstring's old premise that "a name binds once in ga_lua"
    stopped holding when the per-algebra tables were introduced.
    """
    m = re.search(re.escape(table) + r'\.set_function\(\s*"' + re.escape(name) + r'"',
                  text)
    if not m:
        return None
    i = text.index("set_function(", m.start()) + len("set_function(")
    depth = 1
    while depth and i < len(text):
        depth += {"(": 1, ")": -1}.get(text[i], 0)
        i += 1
    return (m.start(), i)  # end just past the matching ')'


def main():
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n", 2)[2],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("functions", nargs="?", default="",
                    help="comma list of names (default: all bound names with "
                    "family overloads)")
    ap.add_argument("--algebra", required=True, choices=list(gen.ALGEBRA_NS),
                    help="family whose overloads to splice in (ega/pga/cga/sta/top)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report changes, write nothing")
    args = ap.parse_args()

    src = HPP.read_text()
    manifest = json.loads(MANIFEST.read_text())
    # bound ON THIS FAMILY'S TABLE -- a name bound only on ega is not bound for
    # cga, and treating it as bound sends the splice into the wrong block
    table = FAMILY_TABLE[args.algebra]
    bound = set(re.findall(re.escape(table) + r'\.set_function\(\s*"([A-Za-z0-9_]+)"',
                           src))

    if args.functions:
        names = [n for n in args.functions.split(",")]
    else:
        fam_ns = gen.ALGEBRA_NS[args.algebra]
        fam_names = {f["name"] for f in manifest["functions"]
                     if f["namespace"] == fam_ns}
        names = sorted(bound & fam_names)

    changed = 0
    for name in names:
        if name not in bound:
            print(f"  {name}: not bound on {table}. "
                  f"(use gen_lua_overloads.py to add a block)")
            continue
        resolves, skipped = family_resolves(manifest, args.algebra, name)
        if skipped:
            print(f"  {name}: {skipped} overload(s) not expressible as sol::resolve "
                  f"(std::pair etc.) -- bind by hand with a lambda, skipped")
        if not resolves:
            continue
        span = _block_span(src, table, name)
        block = src[span[0]:span[1]]
        # idempotency: keep only signatures not already present in the block.
        # Compared with whitespace removed on both sides, so a signature the
        # formatter has wrapped over several lines is still recognised -- see
        # _squeeze(); a plain substring test duplicates those instead.
        flat_block = _squeeze(block)
        new = [r for sig, r in resolves if _squeeze(f"<{sig}>") not in flat_block]
        if not new:
            print(f"  {name}: already up to date")
            continue
        close = _overload_close(src, span[0])
        if close is None or close > span[1]:
            print(f"  {name}: no sol::overload wrapper -- wrap by hand, skipped")
            continue
        ins = (f",\n        // {args.algebra}\n        "
               + ",\n        ".join(new))
        if not args.dry_run:
            src = src[:close] + ins + src[close:]
        changed += 1
        print(f"  {name}: +{len(new)} {args.algebra} overload(s)")

    if not args.dry_run:
        HPP.write_text(src)
    print(f"{'(dry-run) ' if args.dry_run else ''}spliced into {changed} function(s)")


if __name__ == "__main__":
    main()
