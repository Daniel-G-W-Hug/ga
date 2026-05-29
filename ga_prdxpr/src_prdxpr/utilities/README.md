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

### Out of scope (still hand-coded)

- **Duals** (`l_dual`/`r_dual`) and **norms** (`nrm_sq`) live in
  `ga_<algebra>_ops_basics.hpp`, not in `*_ops_products.hpp`. Update them by hand
  from the `ga_prdxpr_rule_generator_test` dual tables.
- Non-primitive products (expansions) are skipped by `--output=code`; they stay
  hand-written and usually update automatically (they delegate to dual+wdg).
- The tool rewrites bodies verbatim — it does not check the math. Always rebuild
  and run the test suite after applying.
