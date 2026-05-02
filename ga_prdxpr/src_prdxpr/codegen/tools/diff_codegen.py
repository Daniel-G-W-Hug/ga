#!/usr/bin/env python3
"""Pair functions by signature between source .hpp and codegen output, diff each pair.

Outputs a summary: total functions on each side, matched count, and per-function diff
chunks for mismatches. Also lists functions present on only one side.
"""
import re
import sys
import difflib

def extract_funcs(path, source_mode=False):
    """Return list of (signature_key, full_text) tuples.

    signature_key normalizes line wraps in the signature so it matches across files
    (whitespace collapsed inside the parenthesized parameter list).
    """
    text = open(path).read()
    chunks = re.split(r'(?=^template <typename)', text, flags=re.MULTILINE)
    funcs = []
    for c in chunks:
        if not c.lstrip().startswith('template <typename'):
            continue
        first_brace = c.find('{')
        if first_brace < 0:
            continue
        depth = 0
        i = first_brace
        end = -1
        while i < len(c):
            if c[i] == '{':
                depth += 1
            elif c[i] == '}':
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
            i += 1
        if end < 0:
            continue
        body = c[:end]
        # Strip comment-only lines (full-line // ...) and trailing spaces
        cleaned = []
        for line in body.splitlines():
            stripped = line.rstrip()
            if stripped.lstrip().startswith('//'):
                continue
            cleaned.append(stripped)
        body = '\n'.join(cleaned)

        # Build signature key: everything up to (and including) the opening '{' header
        sig_end = body.find('\n{', body.find('constexpr'))
        if sig_end < 0:
            continue
        sig_raw = body[:sig_end]
        # Collapse multi-line signature into one line; squash repeated spaces
        sig_key = re.sub(r'\s+', ' ', sig_raw).strip()
        funcs.append((sig_key, body))
    return funcs


# Type-only signature: drops parameter names so functions that differ only in
# parameter naming hash to the same key. Used to detect Cat 5 cases.
#
#   "...constexpr X foo(MVec3d<T> const& A, PScalar3d<U> ps)"
#   "...constexpr X foo(MVec3d<T> const& M, PScalar3d<U> B)"
#   both -> "...constexpr X foo(MVec3d<T> const&, PScalar3d<U>)"
PARAM_NAME_RE = re.compile(
    r"((?:const&|&)?)\s+[A-Za-z_][A-Za-z_0-9]*(\s*[,)])"
)

def type_only_sig(sig_key):
    # Strip trailing identifier in each parameter (after type/qualifiers, before
    # the next ',' or ')'). Matches an identifier preceded by space (possibly with
    # '&' or 'const&' immediately before the space).
    return PARAM_NAME_RE.sub(r"\1\2", sig_key)

def main():
    src_path, gen_path = sys.argv[1], sys.argv[2]
    src = extract_funcs(src_path, source_mode=True)
    gen = extract_funcs(gen_path)
    src_map = {sig: body for sig, body in src}
    gen_map = {sig: body for sig, body in gen}
    common = set(src_map) & set(gen_map)
    only_src = set(src_map) - set(gen_map)
    only_gen = set(gen_map) - set(src_map)

    print(f"=== summary ===")
    print(f"source {src_path}: {len(src)} functions")
    print(f"codegen {gen_path}: {len(gen)} functions")
    print(f"matched signatures: {len(common)}")
    print(f"only in source: {len(only_src)}")
    print(f"only in codegen: {len(only_gen)}")
    print()

    identical = 0
    mismatches = []
    for sig in sorted(common):
        if src_map[sig] == gen_map[sig]:
            identical += 1
        else:
            mismatches.append(sig)
    print(f"=== char-identical: {identical}/{len(common)} ===")
    print(f"=== mismatches: {len(mismatches)} ===")
    print()

    # Show first N mismatches in diff form
    LIMIT = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    for sig in mismatches[:LIMIT]:
        print(f"---- MISMATCH: {sig[:120]} ----")
        diff = difflib.unified_diff(
            src_map[sig].splitlines(),
            gen_map[sig].splitlines(),
            fromfile='source', tofile='codegen', lineterm='', n=2)
        for line in diff:
            print(line)
        print()

    if len(mismatches) > LIMIT:
        print(f"... and {len(mismatches) - LIMIT} more mismatches.")
    print()

    SHOW_FULL = '--full' in sys.argv
    width = 10000 if SHOW_FULL else 140

    if only_src:
        print(f"=== only in source (first 20) ===")
        for sig in sorted(only_src)[:20]:
            print(f"  {sig[:width]}")
        if len(only_src) > 20:
            print(f"  ... and {len(only_src) - 20} more.")
    if only_gen:
        print(f"=== only in codegen (first 20) ===")
        for sig in sorted(only_gen)[:20]:
            print(f"  {sig[:width]}")
        if len(only_gen) > 20:
            print(f"  ... and {len(only_gen) - 20} more.")

    # Cat 5 detection: pair only-in-source with only-in-codegen by type-only sig.
    # Such pairs have matching types & function names but differ in parameter names.
    src_types = {type_only_sig(s): s for s in only_src}
    gen_types = {type_only_sig(s): s for s in only_gen}
    cat5_pairs = []
    for tkey, src_full in src_types.items():
        if tkey in gen_types:
            cat5_pairs.append((src_full, gen_types[tkey]))
    if cat5_pairs:
        print()
        print(f"=== Cat 5 (parameter-naming divergence): {len(cat5_pairs)} pairs ===")
        for src_full, gen_full in cat5_pairs:
            print()
            print(f"  source:  {src_full[:width]}")
            print(f"  codegen: {gen_full[:width]}")

if __name__ == '__main__':
    main()
