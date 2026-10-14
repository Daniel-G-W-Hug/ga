#!/usr/bin/env python3
"""Generate the cross-check tests: every typed product overload and every typed unary
basic operation against the same operation on the full multivector of the algebra.

For each algebra the overloads are read from the LIBRARY itself (libclang over
ga/<family>/ga_<alg>_ops_products.hpp and ..._ops_basics.hpp), so the generated tests
always cover exactly the overloads that exist -- add or remove one and a regeneration follows it; the
--check mode reports a generated file that no longer matches the headers.

The references (see ga_test/src/ga_test_crosscheck.hpp for the contract):

  wdg rwdg dot rdot cmt rcmt operator* rgpr operator<< operator>>
                     the same product on the operands lifted into the full multivector
  join / meet        wdg / rwdg on the lifted operands
  twdg1(a, b)        the grade r+s-2 part of the geometric product (Lengyel's k = 1
                     summand of  a (x) b = sum_k (-1)^(k(k-1)/2) (a x_k b) )
  rtwdg1(a, b)       the same through the complements that define every regressive
                     product: cmpl(op(cmpl a, cmpl b)) in an odd dimension,
                     l_cmpl(op(r_cmpl a, r_cmpl b)) in an even one
  inv / rinv         u (x) inv(u) == 1 == inv(u) (x) u, and u (x) rinv(u) == I with the
                     regressive product
  unary basics       the same operation on the lifted operand, for every unary function
                     of ga_<alg>_ops_basics.hpp that has a full-multivector overload
                     (involutions, complements, duals, bulk/weight, norms, normalizations)

Single fully generic templates (the PGA bulk/weight contractions and expansions, the
generic l_/r_expand) and ega3d's cross() are compositions of the products above and
are tested by hand elsewhere, as are the unary basics WITHOUT a full-multivector form
(the CGA round/flat split, STA's normalize) and those taking a dual number.

Usage (from the repo root, with the ga_bindgen venv, which carries libclang):

    ga_bindgen/.venv/bin/python ga_test/utilities/gen_crosscheck_tests.py
    ga_bindgen/.venv/bin/python ga_test/utilities/gen_crosscheck_tests.py --check

Exit: 0, or 1 when --check finds a generated file out of date.
"""

import argparse
import re
import subprocess
import sys
import zlib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "ga_test/src"
sys.path.insert(0, str(ROOT / "ga_bindgen/src"))

# algebra -> (dimension, full multivector, scalar, pseudoscalar)
ALGEBRAS = {
    "ega2d": (2, "MVec2d", "Scalar2d", "PScalar2d"),
    "ega3d": (3, "MVec3d", "Scalar3d", "PScalar3d"),
    "pga2dp": (3, "MVec2dp", "Scalar2dp", "PScalar2dp"),
    "pga3dp": (4, "MVec3dp", "Scalar3dp", "PScalar3dp"),
    "sta4ds": (4, "MVec4ds", "Scalar4ds", "PScalar4ds"),
    "cga2dc": (4, "MVec2dc", "Scalar2dc", "PScalar2dc"),
    "cga3dc": (5, "MVec3dc", "Scalar3dc", "PScalar3dc"),
}

# product -> how the call is spelled on two operands a, b
CALL = {
    "wdg": "wdg(a, b)", "rwdg": "rwdg(a, b)", "dot": "dot(a, b)", "rdot": "rdot(a, b)",
    "cmt": "cmt(a, b)", "rcmt": "rcmt(a, b)", "operator*": "a * b",
    "rgpr": "rgpr(a, b)", "operator<<": "(a << b)", "operator>>": "(a >> b)",
    "join": "join(a, b)", "meet": "meet(a, b)", "twdg1": "twdg1(a, b)",
    "rtwdg1": "rtwdg1(a, b)",
}
SAME_REF = {"wdg", "rwdg", "dot", "rdot", "cmt", "rcmt", "operator*", "rgpr",
            "operator<<", "operator>>"}
ALIAS_REF = {"join": "wdg(a, b)", "meet": "rwdg(a, b)"}
LABEL = {"operator*": "operator* (gpr)", "operator<<": "operator<< (left contraction)",
         "operator>>": "operator>> (right contraction)"}
# typed unary overloads that DIFFER from the full-multivector form on purpose, with why.
# unitize() divides a SINGLE-coefficient weight out signed (a point or a pseudoscalar
# comes back with weight +1), while unitize() of a multivector divides by the unsigned
# weight norm -- its weight has several coefficients in general, so there is no sign to
# divide by, and an object lifted into it keeps a negative weight's sign. Both are
# pinned by hand in ga_pga_contract_test.hpp.
EXEMPT = {
    ("pga2dp", "unitize", "Vec2dp"), ("pga2dp", "unitize", "Point2dp"),
    ("pga2dp", "unitize", "PScalar2dp"),
    ("pga3dp", "unitize", "Vec3dp"), ("pga3dp", "unitize", "Point3dp"),
    ("pga3dp", "unitize", "PScalar3dp"),
}

ORDER = ["wdg", "join", "rwdg", "meet", "dot", "rdot", "operator<<", "operator>>",
         "cmt", "rcmt", "operator*", "rgpr", "twdg1", "rtwdg1", "inv", "rinv"]


def grade(t, n):
    for prefix, g in (("PScalar", n), ("Scalar", 0), ("QuadVec", 4), ("TriVec", 3),
                      ("BiVec", 2), ("Vec", 1)):
        if t.startswith(prefix):
            return g
    raise ValueError(f"no grade for {t}")


def scan():
    """{algebra: {product: [(param types...)]}} from the library headers."""
    from clang_setup import compile_args  # noqa: E402  (ga_bindgen venv)
    import clang.cindex as cx  # noqa: E402

    args = compile_args()  # sets the libclang library: must precede Index.create
    src = '#include "ga/ga_pga.hpp"\n#include "ga/ga_sta.hpp"\n#include "ga/ga_cga.hpp"\n'
    tu = cx.Index.create().parse("scan.cpp", args=args, unsaved_files=[("scan.cpp", src)])
    errs = [d for d in tu.diagnostics if d.severity >= cx.Diagnostic.Error]
    if errs:
        sys.exit("parse errors:\n" + "\n".join(map(str, errs[:10])))

    def kind(c):  # libclang 18 vs a newer LLVM: unknown cursor kinds raise
        try:
            return c.kind
        except ValueError:
            return None

    found = defaultdict(lambda: defaultdict(list))
    clean = lambda t: re.sub(r"^const |\s*&$|<.*>", "", t.strip()).replace("hd::ga::", "")

    def walk(c):
        for ch in c.get_children():
            k = kind(ch)
            if k == cx.CursorKind.NAMESPACE:
                walk(ch)
                continue
            f = ch.location.file.name if ch.location.file else ""
            m = re.search(r"ga_(\w+?)_ops_(products|basics)\.hpp$", f)
            if not m or k not in (cx.CursorKind.FUNCTION_TEMPLATE,
                                  cx.CursorKind.FUNCTION_DECL) or not ch.is_definition():
                continue
            params = tuple(clean(a.type.spelling) for a in ch.get_children()
                           if kind(a) == cx.CursorKind.PARM_DECL)
            key = ch.spelling if m.group(2) == "products" else "basics:" + ch.spelling
            if params not in found[m.group(1)][key]:
                found[m.group(1)][key].append(params)

    walk(tu.cursor)
    return found


def emit(alg, prods):
    n, full, scalar, pscalar = ALGEBRAS[alg]
    odd = n % 2 == 1
    lines = []
    w = lines.append
    special = {"Point2d": "Vec2dp<value_t>(x.x, x.y, 1)",
               "Point3d": "Vec3dp<value_t>(x.x, x.y, x.z, 1)"}
    w("// Copyright 2024-2026, Daniel Hug. All rights reserved.")
    w("// Licensed under the terms specified in LICENSE.txt file.")
    w("")
    w("// GENERATED by ga_test/utilities/gen_crosscheck_tests.py -- do not edit by hand.")
    w(f"// Every typed product overload of ga_{alg}_ops_products.hpp and every typed unary")
    w(f"// operation of ga_{alg}_ops_basics.hpp against the same operation on the full")
    w("// multivector (see ga_test_crosscheck.hpp). Regenerate after adding, removing or")
    w("// renaming an overload:")
    w("//")
    w("//     ga_bindgen/.venv/bin/python ga_test/utilities/gen_crosscheck_tests.py")
    w("")
    w('#include "ga_test_crosscheck.hpp"')
    w("")
    uses = {t for sigs in prods.values() for s in sigs for t in s if len(s) <= 2}
    if "Plane3d" in uses:
        w("// a plane is built from its trivector, not from four loose coefficients")
        w("template <>")
        w("inline Plane3d<value_t> xcheck::make<Plane3d<value_t>>(xcheck::rng& r)")
        w("{")
        w("    return Plane3d<value_t>(xcheck::make<TriVec3dp<value_t>>(r));")
        w("}")
        w("")
    w(f'TEST_SUITE("{alg}: typed overloads vs the full multivector")')
    w("{")
    w(f"    // any {alg} object as the full multivector")
    w(f"    struct lift_{alg} {{")
    w(f"        template <typename X> {full}<value_t> operator()(X const& x) const")
    w("        {")
    branch = "if"
    for t, expr in special.items():
        if t in uses:
            w(f"            {branch} constexpr (std::is_same_v<X, {t}<value_t>>)")
            w(f"                return {full}<value_t>({expr});")
            branch = "else if"
    w(f"            {'else ' if branch != 'if' else ''}return {full}<value_t>(x);")
    w("        }")
    w("    };")
    total = 0
    for prod in [p for p in ORDER if p in prods and not p.startswith("basics:")]:
        sigs = [s for s in prods[prod] if all(not re.fullmatch(r"A|B|arg\d &", t) for t in s)]
        if not sigs:
            continue
        seed = zlib.crc32(f"{alg}:{prod}".encode())
        w("")
        w(f'    TEST_CASE("{alg}: {LABEL.get(prod, prod)} -- {len(sigs)} overloads")')
        w("    {")
        w(f"        xcheck::rng r({seed}u);")
        w(f"        lift_{alg} const lift;")
        if prod in ("inv", "rinv"):
            one = f"{scalar}<value_t>(1)" if prod == "inv" else f"{pscalar}<value_t>(1)"
            prd = "a * b" if prod == "inv" else "rgpr(a, b)"
            w(f"        auto const inverse_of = [](auto const& x) {{ return {prod}(x); }};")
            w(f"        auto const prd = [](auto const& a, auto const& b) {{ return {prd}; }};")
            for s in sigs:
                (t,) = s
                w(f"        xcheck::inverse<{t}<value_t>>(r, inverse_of, prd, lift, "
                  f"{one}, 1.0e6 * xcheck::tol, \"{prod}({t})\");")
                total += 1
        else:
            w(f"        auto const f = [](auto const& a, auto const& b) {{ return {CALL[prod]}; }};")
            if prod in SAME_REF:
                ref = "f"
            elif prod in ALIAS_REF:
                w(f"        auto const ref = [](auto const& a, auto const& b) {{ return {ALIAS_REF[prod]}; }};")
                ref = "ref"
            for s in sigs:
                ta, tb = s
                if prod in ("twdg1", "rtwdg1"):
                    ra, rb = grade(ta, n), grade(tb, n)
                    if prod == "twdg1":
                        k = ra + rb - 2
                        body = f"gr{k}(a * b)"
                    else:
                        k = (n - ra) + (n - rb) - 2
                        body = (f"cmpl(gr{k}(cmpl(a) * cmpl(b)))" if odd else
                                f"l_cmpl(gr{k}(r_cmpl(a) * r_cmpl(b)))")
                    if not 0 <= k <= n:
                        body = f"{scalar}<value_t>(0)"
                        w(f"        xcheck::binary<{ta}<value_t>, {tb}<value_t>>(r, f, "
                          f"[](auto const&, auto const&) {{ return {body}; }}, lift, "
                          f"\"{prod}({ta}, {tb})\");")
                    else:
                        w(f"        xcheck::binary<{ta}<value_t>, {tb}<value_t>>(r, f, "
                          f"[](auto const& a, auto const& b) {{ return {body}; }}, lift, "
                          f"\"{prod}({ta}, {tb})\");")
                else:
                    w(f"        xcheck::binary<{ta}<value_t>, {tb}<value_t>>(r, f, {ref}, "
                      f"lift, \"{prod}({ta}, {tb})\");")
                total += 1
        w("    }")

    # unary basics with a full-multivector form: the same operation on the lifted operand
    unary = defaultdict(list)
    for key, sigs in prods.items():
        if not key.startswith("basics:"):
            continue
        name = key.split(":", 1)[1]
        # a dual number has no full-multivector form to lift into: tested by hand
        one = [s for s in sigs if len(s) == 1 and not s[0].startswith("DualNum")
               and (alg, name, s[0]) not in EXEMPT]
        if any(s == (full,) for s in one):
            unary[name] = one
    for name in sorted(unary):
        sigs = unary[name]
        seed = zlib.crc32(f"{alg}:basics:{name}".encode())
        w("")
        w(f'    TEST_CASE("{alg}: {name} -- {len(sigs)} overloads")')
        w("    {")
        w(f"        xcheck::rng r({seed}u);")
        w(f"        lift_{alg} const lift;")
        w(f"        auto const f = [](auto const& x) {{ return {name}(x); }};")
        for (t,) in sigs:
            w(f"        xcheck::unary<{t}<value_t>>(r, f, lift, \"{name}({t})\");")
            total += 1
        w("    }")
    w("}")
    return "\n".join(lines) + "\n", total


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="report generated files that differ from what would be written")
    a = ap.parse_args()
    found = scan()
    stale = []
    for alg in ALGEBRAS:
        text, total = emit(alg, found[alg])
        path = OUT_DIR / f"ga_{alg}_crosscheck_test.hpp"
        tmp = path.with_suffix(".hpp.tmp")
        tmp.write_text(text)
        # format inside the repo tree so the repo's .clang-format applies
        subprocess.run(["clang-format", "-i", str(tmp)], check=True)
        text = tmp.read_text()
        tmp.unlink()
        if a.check:
            if not path.exists() or path.read_text() != text:
                stale.append(path.name)
        else:
            path.write_text(text)
        print(f"{alg:7} {total:4d} overloads -> {path.relative_to(ROOT)}")
    if stale:
        print("OUT OF DATE: " + ", ".join(stale))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
