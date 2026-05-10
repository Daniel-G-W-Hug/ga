#!/usr/bin/env python3
"""
library_coverage.py

For each algebra, compare the set of product functions emitted by
`ga_prdxpr --output=code` against the set already present in the
hand-maintained `ga/ga_<algebra>_ops_products.hpp`. Reports which
generated functions still need to be copy-pasted into the ga library.

Signatures are normalized to `funcname(LeftType, RightType)` after
stripping template parameters (`<T>`, `<std::common_type_t<T, U>>`),
qualifiers (`const&`), `[[maybe_unused]]` markers, and parameter
names. Two functions match iff this canonical form matches.

Usage:
    ./library_coverage.py                              # all four default algebras
    ./library_coverage.py --algebra=ega3d              # just one
    ./library_coverage.py --algebra=ega2d,ega3d        # comma-separated list
    ./library_coverage.py --algebra=all                # all five (incl. sta4d)
    ./library_coverage.py --algebra=ega3d --show-code  # also print full bodies
                                                       # of every missing function
                                                       # ready to paste

The --algebra= flag mirrors ga_prdxpr's interface so usage stays consistent
across tools.
"""

import argparse
import difflib
import os
import re
import shutil
import subprocess
import sys


def _canonical_signature_from_line(line):
    """Normalize a single `constexpr ...` line into `name(LType, RType)`.

    Returns None if the line isn't a function declaration we recognise.
    """
    if not line.startswith("constexpr") or "(" not in line:
        return None

    # Strip [[maybe_unused]]
    line = re.sub(r"\[\[maybe_unused\]\]\s*", "", line)

    # Strip template args iteratively from the innermost outward.
    # `<[^<>]*>` only matches a `<...>` containing no further angle
    # brackets, so on the second pass we collapse the now-inner
    # `<std::common_type_t>` (etc.) too.
    while True:
        new_line = re.sub(r"<[^<>]*>", "", line)
        if new_line == line:
            break
        line = new_line

    # Strip qualifiers
    line = re.sub(r"\s*const\s*&\s*", " ", line)
    line = re.sub(r"\s*&\s*", " ", line)

    # Now: "constexpr ReturnType funcname(LType paramL, RType paramR)"
    # function name must allow operator* / operator<< / operator>> too.
    m = re.match(r"^constexpr\s+\S+\s+(\S+?)\s*\(([^)]*)\)", line)
    if not m:
        return None
    name = m.group(1)
    args_str = m.group(2).strip()

    # Each arg's TYPE is the first whitespace-delimited token in the
    # comma-separated chunk; the param name (if any) follows.
    args = []
    if args_str:
        for chunk in args_str.split(","):
            tok = chunk.strip().split()
            if tok:
                args.append(tok[0])

    return f"{name}({', '.join(args)})"


def _join_signature_lines(lines, start):
    """Given a constexpr line at lines[start], join continuation lines until
    the parameter-list parens balance. Returns (joined_text, end_index) where
    end_index is the last line consumed.

    Library headers wrap declarations in two distinct shapes:

      (a) Long signature wrapped after an unclosed '(':
          constexpr MVec3d_E<...> operator*(BiVec3d<T> const& B1,
                                            BiVec3d<U> const& B2)

      (b) Return-type-only first line, function name on the next:
          constexpr Scalar3d<std::common_type_t<T, U>>
          operator<<([[maybe_unused]] BiVec3d<T> const&, ...)

    We handle both: first glue continuation lines until we see a '(' (covers
    case b), then continue gluing until parens balance (covers case a).
    """
    text = lines[start]
    end = start
    # Case (b): no '(' yet — keep gluing until one shows up. Bail out if we
    # hit a ';' (statement terminator → not a function decl) or look too far.
    max_lookahead = 5
    while ("(" not in text) and (";" not in text) \
            and end + 1 < len(lines) and (end - start) < max_lookahead:
        end += 1
        text += " " + lines[end].strip()
    # Case (a): glue until parens balance.
    opens = text.count("(")
    closes = text.count(")")
    while opens > closes and end + 1 < len(lines):
        end += 1
        text += " " + lines[end].strip()
        opens += lines[end].count("(")
        closes += lines[end].count(")")
    return text, end


def extract_signatures(text):
    """Return a set of canonical signatures parsed from a chunk of C++.

    Looks at every line whose first non-whitespace token is `constexpr`,
    joining continuation lines so multi-line signatures match too.
    """
    sigs = set()
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        stripped = lines[i].lstrip()
        if stripped.startswith("constexpr"):
            joined, end = _join_signature_lines(lines, i)
            sig = _canonical_signature_from_line(joined.strip())
            if sig is not None:
                sigs.add(sig)
            i = end + 1
        else:
            i += 1
    return sigs


def extract_blocks(text):
    """Return a dict signature -> full block text from generated codegen output.

    Codegen blocks are delimited by `// <algebra> <product> :: <case_name>`
    comment lines that introduce each function. We accumulate lines into a
    block until we see the next such marker, then map the block by its
    canonical signature.
    """
    sig_to_block = {}
    current = []
    marker_re = re.compile(r"^//\s+\w+\s+\S+\s+::\s+")

    def flush(buf):
        if not buf:
            return
        block = "\n".join(buf).rstrip()
        # Find the constexpr line (joining continuation lines for safety even
        # though codegen currently emits single-line signatures) and
        # canonicalize it.
        for j, line in enumerate(buf):
            stripped = line.lstrip()
            if not stripped.startswith("constexpr"):
                continue
            joined, _ = _join_signature_lines(buf, j)
            sig = _canonical_signature_from_line(joined.strip())
            if sig is not None:
                sig_to_block[sig] = block
                return

    for line in text.splitlines():
        if marker_re.match(line):
            flush(current)
            current = [line]
        else:
            current.append(line)
    flush(current)
    return sig_to_block


def _run_ga_prdxpr(binary_path, algebra):
    proc = subprocess.run(
        [binary_path, f"--algebra={algebra}", "--output=code"],
        capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        sys.exit(f"ga_prdxpr exited {proc.returncode} for {algebra}: "
                 f"{proc.stderr.strip()}")
    return proc.stdout


def get_generated_signatures(binary_path, algebra):
    return extract_signatures(_run_ga_prdxpr(binary_path, algebra))


def get_generated_blocks(binary_path, algebra):
    return extract_blocks(_run_ga_prdxpr(binary_path, algebra))


def get_library_signatures(library_path):
    with open(library_path) as f:
        return extract_signatures(f.read())


def extract_lib_blocks(text):
    """Return dict signature -> (block_text, line_number) for a hand-
    maintained .hpp file.

    Library files don't have `// algebra :: case_name` markers; instead we
    locate each `constexpr ... name(...)` line, walk back through any
    preceding `template`/`requires` lines, and walk forward by counting
    braces until the function body closes. The returned line_number is the
    1-based line of the `constexpr` declaration (typical IDE jump target).
    """
    lines = text.splitlines()
    sig_to_block = {}
    n = len(lines)
    i = 0
    while i < n:
        stripped = lines[i].lstrip()
        if not stripped.startswith("constexpr"):
            i += 1
            continue

        constexpr_line = i + 1   # 1-based for IDE navigation

        # Join multi-line signatures into one logical line for canonicalisation.
        joined_sig, sig_end_idx = _join_signature_lines(lines, i)
        sig = _canonical_signature_from_line(joined_sig.strip())
        if sig is None:
            i += 1
            continue

        # Walk back through template/requires lines that introduce this fn.
        start = i
        while start > 0:
            prev = lines[start - 1].lstrip()
            if prev.startswith(("template ", "template<", "requires")):
                start -= 1
            else:
                break

        # Walk forward from the end of the signature, brace-counting, until
        # the body closes.
        end = sig_end_idx
        depth = 0
        seen_open = False
        while end < n:
            for ch in lines[end]:
                if ch == "{":
                    depth += 1
                    seen_open = True
                elif ch == "}":
                    depth -= 1
                    if seen_open and depth == 0:
                        break
            if seen_open and depth == 0:
                break
            end += 1

        sig_to_block[sig] = ("\n".join(lines[start:end + 1]), constexpr_line)
        i = end + 1
    return sig_to_block


_COMMENT_RE_BLOCK = re.compile(r"/\*.*?\*/", re.DOTALL)
_COMMENT_RE_LINE = re.compile(r"//[^\n]*")


def _strip_comments(code):
    code = _COMMENT_RE_BLOCK.sub("", code)
    code = _COMMENT_RE_LINE.sub("", code)
    return code


def _clang_format(code, formatter_path):
    """Run clang-format on `code`. Returns the formatted text, or `code`
    unchanged if clang-format isn't available or errors out."""
    if formatter_path is None:
        return code
    try:
        proc = subprocess.run(
            [formatter_path], input=code, capture_output=True, text=True,
            check=False)
        if proc.returncode == 0:
            return proc.stdout
    except OSError:
        pass
    return code


def _normalize_for_diff(block, formatter_path):
    # Strip comments, collapse blank lines that result, then clang-format.
    # The generated block always starts with a `// algebra ... :: ...` comment
    # marker that becomes a blank line after stripping; the library block has
    # no such marker. Without this normalisation every comparison would diff
    # by exactly one leading blank line.
    code = _strip_comments(block)
    # Collapse runs of blank lines and strip leading/trailing whitespace.
    code = re.sub(r"\n[ \t]*\n+", "\n", code).strip()
    return _clang_format(code, formatter_path).rstrip()


def _extract_signature_text(block):
    """Return the signature lines (constexpr line through closing paren of
    the parameter list) from a function block, preserving original line
    breaks. Used as a header above each diff so the reader sees the full
    declaration even when the diff itself starts mid-body.
    """
    lines = block.splitlines()
    for i, line in enumerate(lines):
        if line.lstrip().startswith("constexpr") and "(" in line:
            _, end = _join_signature_lines(lines, i)
            return "\n".join(lines[i:end + 1])
    return ""


def report_algebra(binary_path, lib_dir, algebra, show_code, diff_mode,
                   formatter_path):
    lib_file = os.path.join(lib_dir, f"ga_{algebra}_ops_products.hpp")
    print(f"=== {algebra} ===")
    if not os.path.isfile(lib_file):
        print(f"  library file not found: {lib_file}")
        print(f"  (codegen-only — every generated function would be new)")
        print()
        return

    # Always fetch generated blocks when we need full text (show_code or diff).
    need_blocks = show_code or diff_mode
    if need_blocks:
        sig_to_block = get_generated_blocks(binary_path, algebra)
        gen_sigs = set(sig_to_block.keys())
    else:
        gen_sigs = get_generated_signatures(binary_path, algebra)
        sig_to_block = None

    if diff_mode:
        with open(lib_file) as f:
            lib_blocks = extract_lib_blocks(f.read())
        lib_sigs = set(lib_blocks.keys())
        # Workspace-relative path for the location hint (one level above
        # `ga/` is the project root, e.g. `ga/ga_ega3d_ops_products.hpp`).
        proj_root = os.path.dirname(lib_dir)
        lib_rel = os.path.relpath(lib_file, proj_root)
    else:
        lib_blocks = None
        lib_sigs = get_library_signatures(lib_file)
        lib_rel = None

    missing = sorted(gen_sigs - lib_sigs)
    extra = sorted(lib_sigs - gen_sigs)

    print(f"  generated: {len(gen_sigs)} signatures")
    print(f"  library:   {len(lib_sigs)} signatures")
    print(f"  in both:   {len(gen_sigs & lib_sigs)}")

    if missing:
        print(f"\n  TO ADD to library ({len(missing)} signatures present in "
              f"generated code but not in {os.path.basename(lib_file)}):")
        for s in missing:
            print(f"    {s}")

        if show_code:
            print()
            print(f"{'=' * 78}")
            print(f"  Full source for {len(missing)} missing functions "
                  f"(copy-paste into {os.path.basename(lib_file)}):")
            print(f"{'=' * 78}")
            for s in missing:
                block = sig_to_block.get(s)
                if block is None:
                    print(f"\n// (no block found for {s} — internal error)")
                    continue
                print()
                print(block)
            print()

    if extra:
        print(f"\n  in library but NOT in generated ({len(extra)} — likely "
              f"hand-written delegations or out-of-scope helpers):")
        for s in extra:
            print(f"    {s}")
    if not missing and not extra:
        print("\n  signature sets match — library is fully covered")

    if diff_mode:
        in_both = sorted(gen_sigs & lib_sigs)
        differing = []
        for s in in_both:
            lib_text, lib_line = lib_blocks[s]
            gen_norm = _normalize_for_diff(sig_to_block[s], formatter_path)
            lib_norm = _normalize_for_diff(lib_text, formatter_path)
            if gen_norm != lib_norm:
                differing.append((s, gen_norm, lib_norm, lib_line))

        fmt_note = "" if formatter_path else "  (clang-format not found — "\
                                             "comparing comment-stripped only)"
        print(f"\n  diff: {len(in_both)} signatures in both; "
              f"{len(differing)} differ after normalisation{fmt_note}")
        if differing:
            print(f"{'=' * 78}")
            for s, gen_norm, lib_norm, lib_line in differing:
                print(f"\n--- {s} ---")
                # Print the full signature (from generated, since codegen always
                # produces a single-line form that's easy to read) so the
                # reader sees the complete declaration even when the diff
                # context starts deep in the body.
                sig_text = _extract_signature_text(sig_to_block[s])
                if sig_text:
                    print(f"  signature: {sig_text}")
                # Workspace-relative file:line is recognised as a clickable
                # link by most editors / IDEs.
                print(f"  location:  {lib_rel}:{lib_line}")
                diff = difflib.unified_diff(
                    lib_norm.splitlines(keepends=True),
                    gen_norm.splitlines(keepends=True),
                    fromfile=f"library:{os.path.basename(lib_file)}",
                    tofile=f"generated:{algebra}",
                    n=3)
                sys.stdout.writelines(diff)
            print(f"{'=' * 78}")
    print()


def find_binary(proj_root):
    """Locate the ga_prdxpr binary.

    Unix/macOS single-config builds put it at build/ga_prdxpr/ga_prdxpr;
    Windows MSVC multi-config builds put it under build/ga_prdxpr/<Config>/
    with a .exe suffix. We probe both shapes.
    """
    base = os.path.join(proj_root, "build", "ga_prdxpr")
    candidates = [
        os.path.join(base, "ga_prdxpr"),                  # Unix / Ninja / make
        os.path.join(base, "ga_prdxpr.exe"),              # MinGW / Ninja on Win
        os.path.join(base, "Debug", "ga_prdxpr.exe"),     # MSVC Debug
        os.path.join(base, "Release", "ga_prdxpr.exe"),   # MSVC Release
        os.path.join(base, "RelWithDebInfo", "ga_prdxpr.exe"),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    sys.exit(
        "ga_prdxpr binary not found. Looked in:\n  "
        + "\n  ".join(candidates)
        + f"\nBuild it first, e.g.:\n  cd {os.path.join(proj_root, 'build')}"
        f" && cmake --build . --target ga_prdxpr")


def main():
    known_algebras = ["ega2d", "ega3d", "pga2dp", "pga3dp", "sta4d"]

    # Resolve project paths up-front so the default algebra set can be derived
    # from which `ga/ga_<algebra>_ops_products.hpp` files actually exist —
    # adding a new library header automatically expands the default with no
    # script edit needed.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    proj_root = os.path.abspath(
        os.path.join(script_dir, "..", "..", "..", ".."))
    lib_dir = os.path.join(proj_root, "ga")
    default_algebras = [
        a for a in known_algebras
        if os.path.isfile(os.path.join(lib_dir, f"ga_{a}_ops_products.hpp"))
    ]

    parser = argparse.ArgumentParser(
        description="Compare generated product code against the ga library.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Available algebras:\n"
            "  ega2d   Euclidean 2D, G(2,0,0)\n"
            "  ega3d   Euclidean 3D, G(3,0,0)\n"
            "  pga2dp  Projective 2D, G(2,0,1)\n"
            "  pga3dp  Projective 3D, G(3,0,1)\n"
            "  sta4d   Space-time, G(1,3,0)\n"
            "\n"
            "Default selection is derived from existing\n"
            "ga/ga_<algebra>_ops_products.hpp files; algebras whose library\n"
            "file is absent must be requested explicitly with --algebra=NAME.\n"))
    parser.add_argument(
        "--algebra", metavar="ALGEBRAS", default=None,
        help=("comma-separated algebra list (mirrors ga_prdxpr's interface). "
              "Choose any of: " + ",".join(known_algebras) + ". "
              "Special value: 'all' selects every known algebra. "
              "Default: every algebra whose ga/ga_<name>_ops_products.hpp "
              "file exists (currently: "
              + (",".join(default_algebras) if default_algebras else "none")
              + "). "
              "Examples: --algebra=ega3d, --algebra=ega2d,ega3d, --algebra=all"))
    parser.add_argument("--show-code", action="store_true",
                        help="for each missing-from-library signature, also "
                             "print the full generated function block "
                             "(comment + template + body) so it can be "
                             "copy-pasted directly into the .hpp")
    parser.add_argument("--diff", action="store_true",
                        help="for each signature present in BOTH generated "
                             "and library, strip comments and run clang-format "
                             "on each side, then print a unified diff if the "
                             "normalised bodies differ. Picks up the system "
                             "clang-format if present; otherwise falls back to "
                             "comment-stripped textual comparison.")
    args = parser.parse_args()

    # Resolve --algebra= into a concrete list, validating each name.
    if args.algebra is None:
        algebras = list(default_algebras)
    elif args.algebra.strip() == "all":
        algebras = list(known_algebras)
    else:
        algebras = [a.strip() for a in args.algebra.split(",") if a.strip()]
        unknown = [a for a in algebras if a not in known_algebras]
        if unknown:
            sys.exit(f"unknown algebra(s): {', '.join(unknown)}\n"
                     f"choose from: {', '.join(known_algebras)} (or 'all')")

    formatter_path = shutil.which("clang-format")
    binary = find_binary(proj_root)

    for alg in algebras:
        report_algebra(binary, lib_dir, alg, args.show_code, args.diff,
                       formatter_path)


if __name__ == "__main__":
    main()
