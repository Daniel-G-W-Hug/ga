# Codegen validation tools

Helper scripts for comparing `--output=code` results against the existing C++
implementations under `ga/`. Used during Phase 4 char-identity validation; kept around
because they will be useful again when (re-)generating new algebras (e.g. STA4D) or when
auditing future codegen changes.

## diff_codegen.py

Pair functions by signature between two .cpp/.hpp files and report char-identity per
pair. Strips full-line `//` comments before comparing (per project policy: comments are
out of the char-identity goal). Both files should be passed through `clang-format` first
so that line-wrapping doesn't pollute the diff.

```bash
clang-format ga/ga_ega2d_ops_products.hpp                 > /tmp/source.cpp
ga_prdxpr --algebra=ega2d --output=code | clang-format    > /tmp/codegen.cpp
python3 diff_codegen.py /tmp/source.cpp /tmp/codegen.cpp [N] [--full]
```

- `N` (optional, default 6) — number of mismatch diffs to display.
- `--full` — print full signatures in the only-in-source / only-in-codegen lists
  instead of truncating at 140 chars.

Output sections:

- **summary**: function counts on each side, matched signatures, only-source, only-codegen.
- **char-identical**: how many matched pairs are byte-for-byte equal post-format.
- **mismatches**: pairs that share a signature but diverge in body. First N shown as
  unified diff.
- **only in source / only in codegen**: signatures present on only one side. Causes:
  hand-written aliases (e.g. PGA `join`/`meet`), Phase 6 zero-result `[[maybe_unused]]`
  functions, parameter-naming divergence (same operator, different param names → counted
  as separate signatures).

## extract_funcs.py

Standalone function-extractor used by `diff_codegen.py`. Useful in isolation when you
just want a list of function bodies from a header without the surrounding scaffolding.

```bash
python3 extract_funcs.py ga/ga_pga2dp_ops_products.hpp
```

Outputs each function as: signature line(s) + body, separated by blank lines. Strips
trailing whitespace per line; preserves comments (use `diff_codegen.py` for
comment-stripped comparison).

## Limitations

- Both tools key on the full signature *string* after whitespace collapsing. They do
  not normalize parameter names, template-arg ordering, or `const`/`&` placement, so
  cosmetic divergence shows up as "only in one side" rather than as a mismatch.
- Sandwich products are not in scope; codegen skips them entirely.
