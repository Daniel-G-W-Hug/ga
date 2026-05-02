#!/usr/bin/env python3
"""Convert inline-arg flat-constructor return statements to temp-var form.

Rewrites function bodies of the form:

    {
        using ctype = std::common_type_t<T, U>;
        return TypeName<ctype>(expr1, expr2, ..., exprN);
    }

into:

    {
        using ctype = std::common_type_t<T, U>;
        ctype const c0 = expr1;
        ctype const c1 = expr2;
        ...
        ctype const cN = exprN;
        return TypeName<ctype>(c0, c1, ..., cN);
    }

Conversion only applies when:
  * the body contains exactly the `using ctype = ...;` line and a single `return`
    statement (no extra logic like comments inside the body, [[maybe_unused]] params, etc.);
  * the return arguments are "pure expressions" -- no '<' character at the top level
    of the argument list, which would indicate nested constructor calls (Composite types
    are deliberately skipped).

Use:
    python3 inline_to_tempvars.py path/to/file.hpp           # writes back in place
    python3 inline_to_tempvars.py --dry-run path/to/file.hpp # prints diff to stdout
"""
import argparse
import re
import sys
import difflib

# Match the body shape we want to convert. Captures:
#   1: leading whitespace before `using`
#   2: whole `using ctype = ...;` line (verbatim, preserved)
#   3: leading whitespace before `return`
#   4: result type expression (e.g. "MVec2d<ctype>")
#   5: arg list (everything between the outer parens of the constructor call)
#   6: trailing whitespace + closing brace
BODY_RE = re.compile(
    r"""(\n[ \t]*)(using\ ctype\ =\ std::common_type_t<T,\ U>;)"""
    r"""(\n[ \t]*)return\ ([A-Za-z_][A-Za-z_0-9]*(?:<[^>]*>)?)\("""
    r"""((?:[^()]|\([^()]*\))*)\);"""
    r"""(\n[ \t]*\})""",
    re.MULTILINE,
)


def split_top_level_commas(args):
    """Split a comma-separated arg list at depth 0 only."""
    out = []
    cur = ""
    depth = 0
    for ch in args:
        if ch == "(" or ch == "<":
            depth += 1
        elif ch == ")" or ch == ">":
            depth -= 1
        if ch == "," and depth == 0:
            out.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur.strip())
    return out


def has_nested_type(args):
    """True if any '<' appears at depth>=0 (excluding inside parens already), meaning
    one of the args contains a templated constructor call -- a Composite or
    type-expression. Such cases are skipped."""
    return "<" in args


# Match an identifier directly followed by '(' (no whitespace). This signals a function
# call (delegation) like dot(v1, v2) or wdg(a, b), as opposed to pure arithmetic.
FUNC_CALL_RE = re.compile(r"\b[A-Za-z_][A-Za-z_0-9]*\(")


def has_function_call(args):
    """True if any of the arg expressions contain a function call. Used to preserve
    intentional delegations like `MVec2d_E<ctype>(dot(v1, v2), wdg(v1, v2))`."""
    return FUNC_CALL_RE.search(args) is not None


def transform_match(m):
    using_indent = m.group(1)
    using_line = m.group(2)
    return_indent = m.group(3)
    type_expr = m.group(4)
    args_text = m.group(5)
    closing = m.group(6)

    # Skip if Composite (nested constructors). Composite-style returns wrap
    # sub-type constructors (e.g. MVec_E(Scalar(...), BiVec(...))); their inner
    # structure is preserved in source and codegen alike.
    if has_nested_type(args_text):
        return m.group(0)

    # Skip intentional delegations: args contain function calls like dot(...), wdg(...).
    # Codegen handles those as a separate (Phase 5) case.
    if has_function_call(args_text):
        return m.group(0)

    args = split_top_level_commas(args_text)
    if len(args) < 2:
        # 0 or 1 arg -- not worth converting (would be ambiguous with delegations).
        return m.group(0)

    # Detect the return-statement indent (column of "return" keyword) so the temp-var
    # lines align with it.
    body_indent = return_indent.lstrip("\n")

    # Build replacement.
    parts = [using_indent, using_line]
    for i, expr in enumerate(args):
        parts.append("\n" + body_indent + f"ctype const c{i} = {expr};")
    parts.append("\n" + body_indent + f"return {type_expr}(")
    parts.append(", ".join(f"c{i}" for i in range(len(args))))
    parts.append(");")
    parts.append(closing)
    return "".join(parts)


def transform_text(text):
    return BODY_RE.sub(transform_match, text)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("file")
    p.add_argument("--dry-run", action="store_true",
                   help="print unified diff to stdout, do not modify file")
    args = p.parse_args()

    with open(args.file) as f:
        original = f.read()
    transformed = transform_text(original)

    if args.dry_run:
        diff = difflib.unified_diff(
            original.splitlines(),
            transformed.splitlines(),
            fromfile=args.file,
            tofile=args.file + " (transformed)",
            lineterm="",
        )
        for line in diff:
            print(line)
        return

    if original == transformed:
        sys.stderr.write(f"{args.file}: no changes\n")
        return
    with open(args.file, "w") as f:
        f.write(transformed)

    # Count number of conversions performed.
    n = sum(1 for _ in BODY_RE.finditer(original)) - sum(1 for _ in BODY_RE.finditer(transformed))
    sys.stderr.write(f"{args.file}: rewritten in place\n")


if __name__ == "__main__":
    main()
