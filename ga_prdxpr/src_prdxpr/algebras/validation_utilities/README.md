# validation_utilities

Ad-hoc audit scripts for the per-algebra `ProductConfig` cases. These complement
the per-case validator built into `ga_prdxpr` (checks A–F): the validator catches
misconfigurations in **declared** cases, while these scripts find **missing**
cases by comparing pair-set coverage across products and across algebras.

All scripts read the `ga_prdxpr_<algebra>_config.cpp` files in the parent
directory and produce text reports — they do not modify any files.

## Scripts

| Script                                | Question it answers                                                                                                                                                                                                                                                                    |
| ------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `pair_coverage.sh`                    | How many unique (LHS, RHS) operand pairs does each product declare in each algebra? At-a-glance completeness matrix.                                                                                                                                                                   |
| `compare_products.sh <p1> <p2> <alg>` | Within one algebra, which pairs does `p1` declare that `p2` doesn't, and vice versa? Used when a regressive product should mirror its non-regressive counterpart.                                                                                                                      |
| `compare_algebras.sh <product>`       | Across same-dimension algebras, which pairs does each algebra declare for `product`? Used to spot cross-algebra coverage gaps.                                                                                                                                                         |
| `library_coverage.py`                 | Which functions does `ga_prdxpr --output=code` emit that are **not yet** present in the hand-maintained `ga/ga_<algebra>_ops_products.hpp`? Tells you what still needs to be copy-pasted into the ga library after expanding the configs. Requires the `ga_prdxpr` binary to be built. |

## Examples

The scripts resolve their own paths from `BASH_SOURCE` / `__file__`, so they
can be invoked from any working directory. The `./script` form below assumes
you've first `cd`'d into this directory; otherwise pass the full path:

```bash
cd ga_prdxpr/src_prdxpr/algebras/validation_utilities

# Default product set (gpr wdg rgpr rwdg dot rdot cmt rcmt) across all algebras
./pair_coverage.sh

# Custom subset
./pair_coverage.sh gpr wdg rwdg

# Find pairs in pga3dp's gpr that are missing from its wdg
./compare_products.sh gpr wdg pga3dp

# Compare dot coverage across same-dim algebras
./compare_algebras.sh dot

# Show what's emitted by codegen but not yet in the ga library
./library_coverage.py

# Single algebra (--algebra= flag mirrors ga_prdxpr's CLI convention)
./library_coverage.py --algebra=ega3d

# Multiple algebras as a comma-separated list, or 'all' for the full set
./library_coverage.py --algebra=ega2d,ega3d
./library_coverage.py --algebra=all

# Also print the full body of every missing-from-library function
# (ready for direct copy-paste into the .hpp). Best scoped to one algebra.
./library_coverage.py --algebra=ega3d --show-code > /tmp/ega3d_to_add.cpp

# For each function present in BOTH generated code and library, strip
# comments, run clang-format on each side, and diff. Surfaces semantic
# delegations (where the library uses an identity like `wdg(v1, v2)` and
# codegen emits the expanded sum) and any actual implementation drift.
./library_coverage.py --algebra=pga3dp --diff
```

Equivalent invocation from the project root, no `cd` required:

```bash
ga_prdxpr/src_prdxpr/algebras/validation_utilities/pair_coverage.sh
```

## Workflow notes

- Counts come from regex parsing of the `.cases = {...}` arrays. Both the
  infix form (`X * Y -> Z`, `X ^ Y -> Z`) and the functional form
  (`gpr(X,Y) -> Z`, `wdg(X,Y) -> Z`) are recognised; the comparison is on
  the function-name + (LHS, RHS) tuple, so re-formatting an entry from
  infix to functional doesn't change pair counts.
- Whitespace inside parens is normalised: `wdg(mv, mv_e)` and `wdg(mv,mv_e)`
  count as the same pair.
- Sandwich products use a `dummy` case_name and don't contribute to pair
  counts — they're filtered out naturally by the regex.

## Platform notes

- **macOS / Linux**: all scripts work as-is. The bash scripts target bash 3.2
  (macOS's stock `/bin/bash`), so they don't use associative arrays,
  `pipefail`, or `set -e` — anything that would trip on `grep`'s legitimate
  "no matches" exit code.
- **Windows**: the bash scripts run in **Git Bash** or **WSL** (they use
  POSIX utilities that Git Bash bundles); they do **not** run in `cmd.exe`
  or PowerShell. The Python script runs natively on Windows — its binary
  lookup probes both Unix-style (`build/ga_prdxpr/ga_prdxpr`) and the MSVC
  multi-config layout (`build/ga_prdxpr/<Config>/ga_prdxpr.exe`), so it
  works without modification after a Visual Studio build.

## When NOT to use these

- Validating individual case correctness (wrong result type, parity mismatch,
  identically-zero result) — the in-process validator (checks A–F) handles
  that automatically on every `ga_prdxpr` run, with stderr warnings and an
  end-of-run summary.
- Performance benchmarks, codegen diff'ing, or anything beyond pair-set
  comparison — these scripts only see operand pairs, not result types,
  coefficients, or generated output.
