#!/usr/bin/env python3
"""Expand Cat-4 scalar-multiplication delegations in ops_products.hpp files.

Pairs each function in a source .hpp with its codegen counterpart. For any pair
that:
  (a) shares the same (cleaned) signature,
  (b) differs in body, AND
  (c) the source body contains NO semantic-delegation function calls,
the source body is replaced with the codegen body. The source's leading
comments and signature are preserved.

Semantic-delegation calls are lowercase-name function invocations OTHER than
type-cast helpers (`ctype`, `T`, `U`). Examples kept:

    return MVec3d_E<ctype>(dot(a, b), wdg(a, b));    // dot+wdg decomposition
    return wdg(v1, v2);                              // shortcut
    return MVec3d_U<ctype>((v << B), wdg(B, v));     // operator + wdg

Examples expanded (no semantic call -- pure scalar-mult delegation):

    return X * ctype(s);
    return Type<ctype>(reordered...) * ctype(s);
    return ctype(s) * Type<ctype>(...);

Usage:
    python3 expand_delegations.py --source SRC.hpp --codegen CODEGEN.cpp [--dry-run]
"""
import argparse
import re
import sys
from pathlib import Path

# Reuse the cleaned-extractor for signature matching only.
sys.path.insert(0, str(Path(__file__).parent))
from diff_codegen import extract_funcs  # noqa: E402

ALLOWED_CASTS = {"ctype", "T", "U"}
LOWER_CALL_RE = re.compile(r"\b([a-z_][a-zA-Z_0-9]*)\s*\(")
OPERATOR_DELEG_RE = re.compile(r"\(\s*[a-zA-Z_][a-zA-Z_0-9.]*\s*(?:<<|>>)")


def has_semantic_delegation(body):
    # Only inspect the function body proper (between '{' and '}'), not the
    # template/requires preamble (which contains `requires(...)`, `numeric_type<T>`,
    # etc. that would otherwise look like calls).
    open_brace = body.find("\n{")
    close_brace = body.rfind("}")
    if open_brace < 0 or close_brace < 0:
        return False
    body_only = body[open_brace + 2 : close_brace]
    for m in LOWER_CALL_RE.finditer(body_only):
        if m.group(1) not in ALLOWED_CASTS:
            return True
    if OPERATOR_DELEG_RE.search(body_only):
        return True
    return False


def find_function_ranges(text):
    """Yield (sig_key, body_start, body_end) for each function in `text`,
    where body_start..body_end is the substring [from "{" to "}"] inclusive.

    sig_key matches what diff_codegen.extract_funcs produces (after stripping
    full-line // comments and collapsing whitespace inside the signature line).
    """
    pat = re.compile(r"^template <typename", re.MULTILINE)
    for m in pat.finditer(text):
        start = m.start()
        # Find the first '{' after the signature.
        first_brace = text.find("{", start)
        if first_brace < 0:
            continue
        depth = 0
        i = first_brace
        end = -1
        while i < len(text):
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
            i += 1
        if end < 0:
            continue
        # Build cleaned sig_key (matching diff_codegen.extract_funcs)
        chunk = text[start:end]
        cleaned_lines = []
        for line in chunk.splitlines():
            stripped = line.rstrip()
            if stripped.lstrip().startswith("//"):
                continue
            cleaned_lines.append(stripped)
        cleaned = "\n".join(cleaned_lines)
        sig_end = cleaned.find("\n{", cleaned.find("constexpr"))
        if sig_end < 0:
            continue
        sig_raw = cleaned[:sig_end]
        sig_key = re.sub(r"\s+", " ", sig_raw).strip()
        yield sig_key, first_brace, end


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--codegen", required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    src_text = open(args.source).read()
    gen_funcs = dict(extract_funcs(args.codegen))
    src_funcs_cleaned = dict(extract_funcs(args.source, source_mode=True))

    # Walk source bodies in REVERSE order (so earlier offsets stay valid as we splice).
    ranges = list(find_function_ranges(src_text))

    skipped_unmatched = 0
    skipped_identical = 0
    skipped_semantic = 0
    expansions = []  # (sig, body_start, body_end, new_body)

    for sig, b_start, b_end in ranges:
        if sig not in gen_funcs:
            skipped_unmatched += 1
            continue
        if sig not in src_funcs_cleaned:
            skipped_unmatched += 1
            continue
        if src_funcs_cleaned[sig] == gen_funcs[sig]:
            skipped_identical += 1
            continue
        src_body_clean = src_funcs_cleaned[sig]
        if has_semantic_delegation(src_body_clean):
            skipped_semantic += 1
            continue
        # Extract the codegen body's "{...}" chunk.
        gb = gen_funcs[sig]
        gen_open = gb.find("\n{")
        gen_close = gb.rfind("}")
        if gen_open < 0 or gen_close < 0:
            skipped_unmatched += 1
            continue
        new_body_chunk = gb[gen_open + 1 : gen_close + 1]  # starts at '{'
        expansions.append((sig, b_start, b_end, new_body_chunk))

    print(f"functions in source:       {len(ranges)}", file=sys.stderr)
    print(f"unmatched / no codegen:    {skipped_unmatched}", file=sys.stderr)
    print(f"already char-identical:    {skipped_identical}", file=sys.stderr)
    print(f"skipped (semantic):        {skipped_semantic}", file=sys.stderr)
    print(f"to expand:                 {len(expansions)}", file=sys.stderr)

    if args.dry_run:
        print(f"(dry-run -- no file changes)", file=sys.stderr)
        return
    if not expansions:
        return

    # Apply replacements in reverse order to keep earlier offsets valid.
    new_text = src_text
    for sig, b_start, b_end, new_chunk in reversed(expansions):
        new_text = new_text[:b_start] + new_chunk + new_text[b_end:]
    open(args.source, "w").write(new_text)
    print(f"applied:                   {len(expansions)}", file=sys.stderr)


if __name__ == "__main__":
    main()
