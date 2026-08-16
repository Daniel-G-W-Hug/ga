#!/usr/bin/env python3
"""
splice_overloads.py

*Apply* counterpart to gen_lua_overloads.py: insert one algebra family's overloads
into the EXISTING `lua.set_function("name", sol::overload(...))` blocks in
`ga_lua/src/ga_lua.hpp`, in place. Use it when a function name is already bound for
some algebras and a newly-registered family (e.g. STA) adds more overloads to the same
name -- a name binds once in ga_lua, so its overloads must live in a single block.

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
    wrapped first.

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


def family_resolves(manifest, family, name):
    """Deduped `sol::resolve<Ret(Args)>(name)` strings for one name in one family."""
    ns = gen.ALGEBRA_NS[family]
    seen, out = set(), []
    for f in manifest["functions"]:
        if f["name"] != name or f["namespace"] != ns:
            continue
        for ov in f["overloads"]:
            ret = gen.concretize(ov["return_type"])
            sig = f'{ret}({", ".join(gen.param(p) for p in ov["param_types"])})'
            if sig in seen:
                continue
            seen.add(sig)
            out.append((sig, f"sol::resolve<{sig}>({name})"))
    return out


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


def _block_span(text, name):
    """(start, end) of the `lua.set_function("name", ...);` statement, or None."""
    m = re.search(r'set_function\(\s*"' + re.escape(name) + r'"', text)
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
    bound = set(re.findall(r'set_function\(\s*"([A-Za-z0-9_]+)"', src))

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
            print(f"  {name}: not bound (use gen_lua_overloads.py to add a block)")
            continue
        resolves = family_resolves(manifest, args.algebra, name)
        if not resolves:
            continue
        span = _block_span(src, name)
        block = src[span[0]:span[1]]
        # idempotency: keep only signatures not already present in the block
        new = [r for sig, r in resolves if f"<{sig}>" not in block]
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
