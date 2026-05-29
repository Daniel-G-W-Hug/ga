# ga_prdxpr code-update utilities

Helpers for **applying** `ga_prdxpr --output=code` output back into the library
headers. Complements `algebras/validation_utilities/` (which *find* mismatches);
these *apply* the fix.

## `splice_generated_code.py`

Regenerates the selected products for one algebra, clang-formats them with the
project style (the global `~/.clang-format`, found by searching up from a temp file
written inside the repo), and replaces **only the matching function blocks** in
`ga/ga_<algebra>_ops_products.hpp`. Matching is by normalized declaration (return
type + name + argument types), so comments, namespaces, hand-written delegations
(`l_expand4ds`, `r_expand4ds`, …) and unrelated products are left untouched.

Typical use — after a metric / rule change that alters metric-derived products:

```bash
# 1. see what differs
python3 ../algebras/validation_utilities/library_coverage.py --algebra=sta4ds --diff

# 2. preview the splice
python3 splice_generated_code.py --algebra=sta4ds \
    --products=dot,l_contract,r_contract --dry-run

# 3. apply, then rebuild + run the algebra's test suite
python3 splice_generated_code.py --algebra=sta4ds \
    --products=dot,l_contract,r_contract
```

Re-run the coverage diff afterwards: it should report `0 differ`.

### Complements and duals (now code-generated, in `*_ops_basics.hpp`)

`--output=code` now also emits the unary complements/duals
(`l_cmpl`/`r_cmpl`/`cmpl`, `l_dual`/`r_dual`/`dual`, PGA `bulk`/`weight`) — both the
**graded** types (s/vec/bivec/trivec/ps, flat closed forms) and the **aggregates**
(mv/mv_e/mv_u, grade-wise delegations `func(gr_{n-g}(M))`). Point this tool at the
basics header to splice them in:

```bash
python3 splice_generated_code.py --algebra=sta4ds \
    --products=l_cmpl,r_cmpl,l_dual,r_dual \
    --lib ../../../ga/ga_sta4ds_ops_basics.hpp --dry-run
```

Notes: the generated forms are flat closed forms with no inner derivation comments
(so splicing complements that currently carry comments will strip them). For EGA the
hand-coded duals delegate to the complement (`r_dual = r_cmpl`), so the generated flat
EGA duals differ in style; STA duals match byte-for-byte.

### Out of scope (still hand-coded)

- **Norms** (`nrm_sq`/`nrm`) and `normalize` stay hand-written in
  `ga_<algebra>_ops_basics.hpp`.
- Non-primitive products (expansions) are skipped by `--output=code`; they stay
  hand-written and usually update automatically (they delegate to dual+wdg).
- The tool rewrites bodies verbatim — it does not check the math. Always rebuild
  and run the test suite after applying.
