#!/usr/bin/env python3
"""Splice `ga_prdxpr --output=code` function bodies into a library ops_products header.

Purpose
-------
When the extended metric (or any other generator input) changes, the metric-derived
products emitted by `ga_prdxpr --output=code` change too. This tool regenerates the
selected products for one algebra, clang-formats them with the project style, and
replaces *only the matching function blocks* in
`ga/ga_<algebra>_ops_products.hpp` — every comment, namespace, hand-written
delegation and unrelated product in the library is left byte-for-byte intact.

It matches functions by their normalized declaration (return type + name + argument
types), so an overload is replaced iff the generator emits one with the same
signature. Functions present in the library but not in the generated set (e.g. the
hand-written `l_expand4ds` / `r_expand4ds` delegations) are never touched.

This automates what was previously a hand-transcription step (a known source of
sign errors when copying rule-generator output by hand).

Workflow it replaces / complements
-----------------------------------
- `library_coverage.py --algebra=A --diff` tells you *which* functions differ.
- This tool *applies* the generated versions of those functions.
Run the coverage diff again afterwards; it should report `0 differ`.

Prerequisites
-------------
- `ga_prdxpr` built (default: build/ga_prdxpr/ga_prdxpr under the repo root).
- `clang-format` on PATH (or pass --clang-format). The project style is the global
  ~/.clang-format; clang-format finds it by searching up from the formatted file,
  so the temp file is written *inside the repo* on purpose.

Examples
--------
    # Dry run: show what would change for sta4ds dot + contractions
    python3 splice_generated_code.py --algebra=sta4ds \
        --products=dot,l_contract,r_contract --dry-run

    # Apply
    python3 splice_generated_code.py --algebra=sta4ds \
        --products=dot,l_contract,r_contract

Notes / gotchas
---------------
- Only *primitive* products are emitted by `--output=code`; expansions and other
  non-primitive products are skipped by the generator (and stay hand-written in the
  library — they usually delegate to dual+wdg and update automatically).
- Hand-coded duals (`l_dual`/`r_dual`) and norms (`nrm_sq`) live in
  `ga_<algebra>_ops_basics.hpp`, NOT here — update those by hand from the
  rule-generator dual tables (`ga_prdxpr_rule_generator_test`).
- This rewrites function bodies verbatim from the generator; it does not validate
  math. Always rebuild + run the algebra's test suite afterwards.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)


def split_functions(lines):
    """Yield (start, end_exclusive, block_lines) for each top-level template function.

    A function starts at a `template <` line and ends at the first line that is
    exactly `}` (closing brace at column 0)."""
    i, n = 0, len(lines)
    while i < n:
        if lines[i].startswith("template <"):
            j = i
            while j < n and lines[j].rstrip("\n") != "}":
                j += 1
            yield (i, j + 1, lines[i:j + 1])
            i = j + 1
        else:
            i += 1


def decl_key(block):
    """Normalized declaration key: the `constexpr`/`inline` line(s) up to `{`,
    whitespace-collapsed. Captures return type + name + argument types."""
    text = "".join(block)
    m = re.search(r"\n(constexpr|inline)\b(.*?)\{", text, re.DOTALL)
    if not m:
        return None
    return re.sub(r"\s+", " ", m.group(1) + m.group(2)).strip()


def main():
    ap = argparse.ArgumentParser(
        description="Splice ga_prdxpr --output=code functions into a library header.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--algebra", required=True,
                    help="algebra name (e.g. sta4ds, ega3d)")
    ap.add_argument("--products", required=True,
                    help="comma-separated products to regenerate (e.g. dot,l_contract)")
    ap.add_argument("--ga-prdxpr", default=None,
                    help="path to ga_prdxpr binary (default: build/ga_prdxpr/ga_prdxpr)")
    ap.add_argument("--lib", default=None,
                    help="library header to patch (default: ga/ga_<algebra>_ops_products.hpp)")
    ap.add_argument("--clang-format", default="clang-format",
                    help="clang-format executable (default: clang-format on PATH)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report changes without writing the library file")
    args = ap.parse_args()

    ga_prdxpr = args.ga_prdxpr or os.path.join(REPO_ROOT, "build", "ga_prdxpr", "ga_prdxpr")
    lib = args.lib or os.path.join(REPO_ROOT, "ga", f"ga_{args.algebra}_ops_products.hpp")

    if not os.path.isfile(ga_prdxpr):
        sys.exit(f"error: ga_prdxpr not found at {ga_prdxpr} (build it first)")
    if not os.path.isfile(lib):
        sys.exit(f"error: library header not found at {lib}")

    # 1) generate raw single-line code
    gen_raw = subprocess.run(
        [ga_prdxpr, f"--algebra={args.algebra}",
         f"--products={args.products}", "--output=code"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=True,
    ).stdout

    # 2) clang-format (temp file under the repo so the global ~/.clang-format is found)
    with tempfile.NamedTemporaryFile("w", suffix=".cpp", dir=REPO_ROOT, delete=False) as tf:
        tf.write(gen_raw)
        tmp_path = tf.name
    try:
        subprocess.run([args.clang_format, "-i", "--style=file", tmp_path], check=True)
        with open(tmp_path) as f:
            gen_lines = f.readlines()
    finally:
        os.unlink(tmp_path)

    gen_map = {}
    for _, _, block in split_functions(gen_lines):
        k = decl_key(block)
        if k:
            gen_map[k] = block
    print(f"generated functions: {len(gen_map)}")

    with open(lib) as f:
        lib_lines = f.readlines()

    spans = list(split_functions(lib_lines))
    changed_keys, matched_keys = [], set()
    for (start, end, block) in reversed(spans):
        k = decl_key(block)
        if k in gen_map:
            matched_keys.add(k)
            if "".join(gen_map[k]) != "".join(block):
                changed_keys.append(k)
            lib_lines[start:end] = gen_map[k]

    print(f"library functions matched: {len(matched_keys)} / {len(gen_map)}")
    print(f"library functions changed: {len(changed_keys)}")
    for k in changed_keys:
        print("   changed:", k[:90])

    missing = set(gen_map) - matched_keys
    if missing:
        print(f"WARNING: {len(missing)} generated signatures not found in library:")
        for k in sorted(missing):
            print("   missing:", k[:90])

    if args.dry_run:
        print("(dry run: library not modified)")
    elif changed_keys:
        with open(lib, "w") as f:
            f.writelines(lib_lines)
        print(f"wrote {lib}")
    else:
        print("(nothing changed: library not modified)")


if __name__ == "__main__":
    main()
