#!/usr/bin/env python3
"""
gen_lua_overloads.py

Draft `sol::resolve` overload blocks for `ga_lua`'s `register_functions`, generated
from `ga_bindgen/manifest.json`. Complements `lua_coverage.py`: that tool *finds* the
gaps, this one *drafts* the binding code to close them (the find-vs-apply split mirrors
ga_prdxpr's validation_utilities vs utilities).

For each requested function name it emits one `lua.set_function("name", sol::overload(
... ))` block listing every manifest overload, filtered to one algebra family, in the
exact `sol::resolve<Ret(Args)>(name)` style used throughout `ga_lua.hpp`.

How signatures are concretised (manifest stores template forms like
`Vec2d<std::common_type_t<T, U>>`):

  - the outer type name is lowercased to the Lua alias: `Vec2d`->`vec2d`,
    `MVec2d_E`->`mvec2d_e`, `PScalar3d`->`pscalar3d`, `BiVec3d`->`bivec3d`;
  - `value_t` / `double` -> `value_t`, `bool` -> `bool`, `std::vector<X>` recurses;
  - reference parameters (`const T &`) keep `const&`; value parameters (scalars,
    pseudoscalars, `value_t`) are passed by value -- matching the library's calling
    convention and the existing `sol::resolve` blocks.

Family filter (`--algebra`): selected by the manifest `namespace` of each overload
(`hd::ga::ega` / `hd::ga::pga` / `hd::ga::sta`). So a multi-algebra function (e.g.
`rwdg`, `angle`, whose name appears once per namespace) yields only the slice for the
family you ask for. The generic top-level `hd::ga` overloads (e.g. `sqrt`, `gr`, `sign`)
are NOT emitted by a family run -- handle those separately.

Usage:
    ./gen_lua_overloads.py rotate,project_onto,angle            # default --algebra=ega
    ./gen_lua_overloads.py rwdg,rgpr,rinv --algebra=pga
    ./gen_lua_overloads.py transform,time_split --algebra=sta
    ./gen_lua_overloads.py sqrt,log --algebra=ega,pga           # merge families

A name that spans several families (e.g. `sqrt`, `log`, `twdg1`) must bind ONCE in
`ga_lua` (a second `set_function("sqrt", ...)` replaces the first). Pass all the
families it covers in one `--algebra=` list so the emitted block carries every
overload; identical signatures (the manifest lists a function once per header it
appears in) are de-duplicated.

This is a DRAFTING aid: review the emitted block (return-type / by-value-vs-ref
correctness, and whether an overload is wanted in the REPL) before pasting it into the
relevant `register_*` section, then build `ga_lua` and re-run `lua_coverage.py`.
"""

import argparse
import json
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST = REPO_ROOT / "ga_bindgen" / "manifest.json"

ALGEBRA_NS = {
    "ega": "hd::ga::ega",
    "pga": "hd::ga::pga",
    "sta": "hd::ga::sta",
    "top": "hd::ga",  # top-level free functions (gr, sign, unit conversions, ...)
}


def _base(t):
    return t.replace("const", "").replace("&", "").strip()


def concretize(t):
    """Template type string -> concrete Lua alias (no reference)."""
    t = _base(t)
    if t in ("value_t", "double"):
        return "value_t"
    if t == "bool":
        return "bool"
    if t.startswith("std::common_type_t"):   # bare scalar result -> value_t
        return "value_t"
    if re.fullmatch(r"[A-Z]", t):            # bare template param (scalar arg) -> value_t
        return "value_t"
    mv = re.match(r"std::vector<(.+)>$", t)
    if mv:
        return f"std::vector<{concretize(mv.group(1))}>"
    name = re.match(r"([A-Za-z_][A-Za-z0-9_]*)", t).group(1)
    return name.lower()


def param(t):
    alias = concretize(t)
    return f"{alias} const&" if "&" in t else alias


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 2)[2],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("functions", help="comma-separated function names")
    ap.add_argument("--algebra", default="ega",
                    help="comma-separated families to merge (ega,pga,sta; default: ega). "
                    "A name with overloads in several gets ONE set_function carrying "
                    "all listed families' overloads (each name binds once in ga_lua).")
    args = ap.parse_args()

    families = args.algebra.split(",")
    bad = [a for a in families if a not in ALGEBRA_NS]
    if bad:
        ap.error(f"unknown algebra(s) {bad}; choose from {list(ALGEBRA_NS)}")
    target_ns = {ALGEBRA_NS[a] for a in families}
    # a function name appears once per namespace; merge the requested families' overloads
    fns = {}
    for f in json.loads(MANIFEST.read_text())["functions"]:
        if f["namespace"] in target_ns:
            fns.setdefault(f["name"], []).extend(f["overloads"])

    for fname in args.functions.split(","):
        overloads = fns.get(fname)
        if overloads is None:
            print(f"// {fname}: no {args.algebra} overloads (not in {target_ns})\n")
            continue
        lines = []
        seen = set()
        for ov in overloads:
            ret = concretize(ov["return_type"])
            sig = f'{ret}({", ".join(param(p) for p in ov["param_types"])})'
            if sig in seen:  # same signature can appear in several headers
                continue
            seen.add(sig)
            lines.append(f"sol::resolve<{sig}>({fname})")
        print(f'lua.set_function("{fname}",')
        print("                 sol::overload(")
        for i, line in enumerate(lines):
            tail = "));" if i == len(lines) - 1 else ","
            print(f"                     {line}{tail}")
        print()


if __name__ == "__main__":
    main()
