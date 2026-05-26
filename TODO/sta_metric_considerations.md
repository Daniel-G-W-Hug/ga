# STA4ds metric, dualization, and the metric / antimetric duality

Working notes on how dualization behaves in `sta4ds` = G(1,3,0), why the naive
metric-dual round trip does not return to the identity, and the principled,
metric-independent round trip that does. Also records the bridge to a future
projective spacetime algebra.

All statements below are verified in code in
[ga_test/src/ga_sta4ds_test.hpp](../ga_test/src/ga_sta4ds_test.hpp) (test cases
"MVec4ds: dualization - complement vs. pseudoscalar-multiplication",
"MVec4ds: metric / antimetric exomorphisms (G, rG) and metric-indep. dual",
"MVec4ds: left/right antidual; the wedge product fixes space vs dual", and
"MVec4ds: explicit space -> dual -> space round trip"). The ega2d analogue is in
`ga_ega2d_test.hpp` ("MVec2d: dualization - complement vs. pseudoscalar-multiplication");
the ega3d original in `ga_ega3d_test.hpp`.

> **RESUME HERE (paused 2026-05-25).** Open thread: §6 — STA's stored extended metric is
> the **blade square** `Q`, EGA's is the **pure product** `P`; they differ by the
> reversion sign `sigma(k)`, and STA's is **not a wedge exomorphism**. Undecided whether
> this is an intended convention or a **sign bug** in the `is_minkowski` override of
> `calculate_extended_metric` (ga_prdxpr). **Agreed next step: build the exomorphism
> comparison table (§13)** — one section per algebra, reading the stored arrays, tabulating
> `P`, `Q`, the stored value, and the multiplicativity verdict `g_{S∪T} == g_S·g_T`. No
> library behaviour is to change before that table + the deciding experiment in §13.

## 1. Two schools of dualization

There are two common ways to dualize a blade / multivector `A`:

1. **Complement after metric** (this library, after Lengyel): `dual(A) = cmpl(G·A)`,
   where `G` is the extended metric and `cmpl` is the purely combinatorial complement
   (defined by `A ^ cmpl(A) = I`). Split into a left and a right variant in even
   dimension (see §4).
2. **Pseudoscalar multiplication** (the "other school"): `to_dual(A) = A ⟑ I`,
   `from_dual(A) = A ⟑ inv(I)`.

In `sta4ds` the metric is non-degenerate, so both schools exist and agree on the
forward map (see §5): `A ⟑ I == r_dual(A)` exactly, at every grade.

## 2. Why pseudoscalar multiplication works here but dies in PGA

The pseudoscalar-multiplication scheme needs `I` to be invertible. In `sta4ds`
`I ⟑ I == -1`, so `inv(I) == -I` and `(A ⟑ I) ⟑ inv(I) == A`.

In a **projective / degenerate** algebra (PGA, and a future projective STA) one basis
vector squares to zero. That null vector sits inside `I`, so `I ⟑ I == 0` and `I` has
**no inverse** — the pseudoscalar-multiplication scheme is structurally impossible
there. This is not a defect: it shows that pseudoscalar multiplication was only ever
the non-degenerate shortcut for `cmpl(G·A)`. The complement-based dual survives because
`I` never appears in it.

## 3. The dual is `complement ∘ metric`

`dual(A) = cmpl(G·A)`. The complement alone is metric-blind (it maps a blade to its
complementary blade with a sign so that `A ^ cmpl(A) = I`, the same for any signature).
The metric re-enters only through the `G·A` factor. So "just use the complement" is not
enough if the dual space is meant to carry metric structure — the metric must appear.

## 4. Handedness: `l_dual` vs `r_dual` (even dimension)

`sta4ds` is 4-dimensional (even), so left and right complements / duals differ for the
odd-grade elements (vector, trivector) and agree for the even grades (scalar, bivector,
pseudoscalar):

- `l_dual(A) = l_cmpl(G·A)`,  `r_dual(A) = r_cmpl(G·A)`
- `l_dual(A) == r_dual(A)` for grades 0, 2, 4;  `l_dual(A) == -r_dual(A)` for grades 1, 3.

The non-metric complements are exact mutual inverses for any metric:
`l_cmpl(r_cmpl(A)) == r_cmpl(l_cmpl(A)) == A`.

## 5. The metric `G` and antimetric `rG` as exomorphisms

Following Lengyel, the metric carries two exomorphisms:

- metric `G`: a **wedge** exomorphism, `G(wdg(a,b)) == wdg(G(a), G(b))`
- antimetric `rG`: an **antiwedge** exomorphism, `rG(rwdg(a,b)) == rwdg(rG(a), rG(b))`

related by ordinary matrix product

```text
G · rG = det(metric) · I        (I = identity matrix)
```

For `sta4ds`: `det(metric) = (-1)(-1)(-1)(+1) = -1`, and every basis blade squares to
±1 so `G² = I`. The relation therefore collapses to

```text
rG = -G                          (the antimetric is the negated metric)
```

`det(metric)` is signature-parity invariant: `det = (-1)^q` with `q` the number of
negative directions, so both physically-equivalent signatures `(-,-,-,+)` (q=3) and
`(+,+,+,-)` (q=1) give `det = -1`. The whole construction is signature-independent.

## 6. The extended metric: pure product (EGA) vs blade square (STA) — OPEN / possible bug

This is the crux, and it is **not settled theory** — it is an open question whose
resolution is the comparison table in §13. The extended metric stored in
`ga_usr_consts.hpp` uses a **different convention in STA than in EGA**, and STA's choice
is not a wedge exomorphism. Whether that is an intended convention or a sign bug is still
to be decided.

### Two scalars per basis blade (precise definitions)

Write each basis vector's square as `g_i^2 = m_i` (the signature value, in {-1,0,+1}).
For a grade-`k` basis blade `e_S` (product of `k` distinct vectors in canonical order):

- **Pure product**:  `P(e_S) = prod_{i in S} m_i`            (just multiply the squares)
- **Blade square**:  `Q(e_S) = e_S . e_S = sigma(k) * P(e_S)` (the blade times itself)

with reversion sign `sigma(k) = (-1)^(k(k-1)/2)`. `P` and `Q` are **equal** where
`sigma(k) = +1` (grades k = 0,1 mod 4) and **opposite** where `sigma(k) = -1`
(grades 2,3 mod 4). A third related quantity: the reverse scalar product
`<rev(e_S).e_S>_0 = P` — so "pure product" is also the with-reverse norm.

Worked examples in STA `(m1=m2=m3=-1, m4=+1)`:

```text
  g14:  P = m1*m4 = (-1)(+1) = -1      Q = sigma(2)*P = (-1)(-1) = +1   ( = g14.g14 )
  g23:  P = m2*m3 = (-1)(-1) = +1      Q = sigma(2)*P = (-1)(+1) = -1   ( = g23.g23 )
```

Direct check: `g14.g14 = g1g4g1g4 = -g1^2 g4^2 = -(-1)(+1) = +1 = Q(g14)`; and
`g23.g23 = -g2^2 g3^2 = -(-1)(-1) = -1 = Q(g23)`.

### What the stored arrays actually contain

- `ega3d_metric` (and ega2d / pga): the **pure product** `P`. E.g. `e12 -> +1 = P(e12)`,
  although the blade square `Q(e12) = e12^2 = -1`. For Euclidean this is just the
  identity (all `+1`).
- `sta4ds_metric`: the **blade square** `Q`. All 16 entries match `Q` — e.g. `g14 -> +1`,
  `g23 -> -1`, `I -> -1`; the pure product would instead give `g14 -> -1`, `g23 -> +1`.

So **EGA stores `P`, STA stores `Q = sigma(k)*P`** — the single difference is `sigma(k)`.

### Why it matters: STA's stored metric is NOT a wedge exomorphism

A diagonal metric `g_S` is a genuine wedge exomorphism (Lengyel's
`G(a^b) = G(a)^G(b)`, the "conforming property" CLAUDE.md says the extended metric must
satisfy) **iff** it is multiplicative: `g_{S∪T} = g_S * g_T` for disjoint S,T.

- `P` is multiplicative by construction → EGA's stored metric **is** a wedge exomorphism.
- `Q = sigma(k)*P` is **not** multiplicative (`sigma(k+l) = (-1)^(kl) sigma(k) sigma(l)`)
  → STA's stored metric is **not** a wedge exomorphism.

Confirmed explicitly in code (`ga_sta4ds_test.hpp`, test "metric / antimetric
exomorphisms"): the **pure-product reconstruction** `G(X) = rev(Gx(X))` satisfies
`G(wdg(a,b)) == wdg(G(a),G(b))` (four blade pairs incl. mixed grade), while the library's
own metric `Gx(X) = l_cmpl(r_dual(X))` — which reproduces `sta4ds_metric` entry-for-entry
(`Gx = l_cmpl o r_dual` because `r_dual = r_cmpl o G` and `l_cmpl o r_cmpl = id`) —
**fails** it: `Gx(wdg(g1,g4)) != wdg(Gx(g1),Gx(g4))`. So when earlier notes said "sta4ds
holds a wedge exomorphism", that was the *pure-product reconstruction* `G`, NOT the stored
metric; the stored metric is explicitly confirmed to fail. (Tested for a few pairs only —
the exhaustive check is the §13 table.)

### Why no existing test caught it

The extended metric is **architecturally decoupled** from the rest: the wedge product is
metric-free, the geometric product uses the vector signature `m_i` directly, and the
complements are purely combinatorial. The stored array feeds **only** `dot` / `nrm_sq`
and the metric dual `cmpl(G.A)`. The transcription gate `dual(e) == nrm_sq(e)*cmpl(e)`
and the contraction identity `a<<b == rwdg(l_dual(a),b)` hold **by construction** for
whatever the metric is (they check internal consistency, not correctness). So the 508
passing STA assertions do **not** validate the extended metric's exomorphism property —
a non-exomorphism metric sails through untouched. The `==`/`!=` checks above are the first
ones that even look.

### Consequence and locus

Everything that looked strange in STA — the `gr_inv` in `l_dual(r_dual(u))` (§7), the
antidual / `r_undual` (§8), the `A.I == r_dual` coincidence — is a direct consequence of
this one `sigma(k)`. With the pure product `P` (as EGA uses), STA's dual would be the
with-reverse Hodge dual `rev(A).I` and `l_dual(r_dual(u)) == u` would hold exactly as in
EGA; none of that apparatus would arise.

If the stored STA metric is wrong, the locus is narrow: the `is_minkowski` override in
[calculate_extended_metric](../ga_prdxpr/src_prdxpr/rules/ga_prdxpr_rule_generator.cpp)
(lines ~277–291). It overrides the product-of-vector-metrics with the blade-square values
for grades 2,3; its comment claims to derive them "from `<A,~A>`" but `<g14,~g14>_0 = -1`
while the code sets `+1` — i.e. the emitted values are `<A.A>_0` (blade square), the
*opposite* sign. The complement / product / dual **generation** is independent and
well-tested, so those are the less likely suspects.

**Verdict pending:** is the blade-square STA metric an intentional convention (e.g. a
Dorst-style `<A.B>_0` scalar product) or a sign bug? The §13 comparison table is the test.

## 7. The naive round trip, and why it fails

**What "naive" means, explicitly.** In an even-dimensional algebra the obvious way to
go to the dual space and back is: dualize with one handedness, undualize with the
other. That is the *naive metric-dual round trip*

```text
forward  (space -> dual):   D  = r_dual(A)
backward (dual  -> space):   A' = l_dual(D) = l_dual(r_dual(A))
expectation:                 A' == A
```

In the Euclidean even algebra `ega2d` this expectation is **true** (`det = +1`):

```text
ega2d:    l_dual(r_dual(A)) == A
```

In `sta4ds` the **same** composition does **not** return `A`. Because both `l_dual` and
`r_dual` apply the metric `Gx` on the *primal* side, the round trip applies the metric
twice on the same side and picks up `det · gr_inv`:

```text
sta4ds:   l_dual(r_dual(A)) == r_dual(l_dual(A))
                            == det(metric) · gr_inv(A)
                            == -gr_inv(A)            (det = -1; gr_inv(A) = (-1)^k A)
```

so by grade: `s -> -s`, `v -> +v`, `B -> -B`, `t -> +t`, `ps -> -ps`. The two
metric duals agree with each other, but neither is the inverse of the other.

It cannot be repaired by the determinant alone:

```text
det · l_dual(r_dual(A)) == gr_inv(A)      (still grade-flipped, not A)
```

Healing needs the determinant **and** the grade involution:

```text
A == det(metric) · gr_inv( l_dual(r_dual(A)) )
```

Equivalently `det · gr_inv` is an involution, so applying the naive pair twice returns
`A`:  `l_dual(r_dual(l_dual(r_dual(A)))) == A`.

## 8. The metric-independent round trip (space → dual → space)

The principled fix is to apply the metric on the **dual** side coming back, not the
primal side again. Since `Gx² = I` (so `Gx⁻¹ = Gx`) the inverse of `r_dual = r_cmpl ∘ Gx`
is `Gx ∘ l_cmpl`, the **antidual**:

```text
forward  (space -> dual):   D = r_dual(A)                  = r_cmpl(G · A)
backward (dual  -> space):  A = antidual(D)                = l_cmpl(r_dual(l_cmpl(D)))
                                                            = Gx⁻¹ · l_cmpl(D)
                                                            = (1/det) · rG · l_cmpl(D)
```

`antidual(r_dual(A)) == A` and `r_dual(antidual(A)) == A` at every grade, for any
non-degenerate signature (it depends on the metric only through `Gx² = I`). This is the
round trip that uses **both** `G` (going out) and `rG` (coming back), as desired.

Contrast with §7: the naive backward leg is `l_dual = l_cmpl ∘ Gx` (metric on the primal
side); the correct backward leg is `Gx ∘ l_cmpl` (metric on the dual side). Same two
ingredients, opposite order.

**Left and right antidual.** Because the algebra is even-dimensional there are two duals
(`l_dual`, `r_dual`), and correspondingly two antiduals — one inverse per handedness:

```text
right antidual (inverts r_dual):   r_undual(D) = l_cmpl(r_dual(l_cmpl(D)))   = Gx o l_cmpl
left  antidual (inverts l_dual):   l_undual(D) = r_cmpl(l_dual(r_cmpl(D)))   = Gx o r_cmpl
```

The metric `Gx` sits on the dual side in both; the complement that opens and closes the
expression matches the dual's handedness. They coincide on even-grade inputs (where
`l_dual == r_dual`) and differ by a sign on odd grades, exactly as `l_dual`/`r_dual` do.
Each inverts only its own handedness: the cross pairing flips the odd grades
(`l_undual(r_dual(v)) == -v`, `r_undual(l_dual(v)) == -v`).

Two other exact, metric-independent round trips already exist and are simpler:

- complement pair (metric-free): `l_cmpl(r_cmpl(A)) == A`
- pseudoscalar group inverse: `(A ⟑ I) ⟑ inv(I) == A`, with `inv(I) == -I`

## 9. Signature duality `(-,-,-,+)` <-> `(+,+,+,-)`, and which is "the space"

Because `rG = -G`, the antimetric negates every vector square: `G(g1) = -g1` (signature
`(-,-,-,+)`) while `rG(g1) = +g1` (signature `(+,+,+,-)`). The two signatures are the
wedge-metric and antiwedge-antimetric of one another — i.e. **dual spaces**.

**Can we say which one is "the space"?** Yes, relative to the choice of wedge product —
because each signature is an exomorphism of *exactly one* product, never both:

```text
G  = (-,-,-,+)   is the WEDGE     exomorphism   (G(wdg(a,b))  == wdg(G a, G b)),  NOT antiwedge
rG = (+,+,+,-)   is the ANTIWEDGE exomorphism   (rG(rwdg(a,b)) == rwdg(rG a,rG b)), NOT wedge
```

So given the library's wedge convention, the **space** is uniquely the signature whose
metric respects the wedge product, `(-,-,-,+)`, and its **dual space** is `(+,+,+,-)`.
The assignment is not arbitrary: the wedge product singles it out.

But the whole structure is symmetric under `(wdg <-> rwdg, G <-> rG, space <-> dual)`. If
one declares `rwdg` to be the primary wedge (i.e. works in the dual algebra), then
`(+,+,+,-)` becomes *the* wedge-exomorphism = "the space" and `(-,-,-,+)` its dual. So
there is no absolute, convention-free distinction — only one relative to which product
is called the wedge. This is precisely why both signatures are physically equivalent:
each is "the space" for one of the two products, the other being its dual.

## 10. A unified round trip across the algebras

The handedness names hide a simpler truth: **`l_dual` was never the inverse of `r_dual` —
it only coincides with the inverse when the metric is the identity.** The genuine inverse
is the un-dual, with the metric applied on the dual side:

```text
r_undual(D) = G⁻¹ · l_cmpl(D) = (1/det) · rG · l_cmpl(D) = l_cmpl(r_dual(l_cmpl(D)))
```

Applying the *same* sta4ds formula `l_cmpl(r_dual(l_cmpl(D)))` verbatim in `ega2d` gives
`r_undual == l_dual` exactly (verified: `r_undual(v) == l_dual(v) == Vec2d(2,-1)`), and
`r_undual(r_dual(u)) == u` at every grade. So the familiar Euclidean identity
`l_dual(r_dual(u)) == u` is just the `G = I` special case of the universal law

```text
r_undual(r_dual(u)) == u            (l_undual(l_dual(u)) == u for the other handedness)
```

**One formula, all non-degenerate algebras.** The implementable form
`l_cmpl(r_dual(l_cmpl(·)))` needs no per-algebra branching — the metric is already inside
`r_dual`. The determinant is the only thing that varies:

| algebra            | det  | `rG`     | `r_undual` reduces to |
| ------------------ | ---- | -------- | --------------------- |
| ega2d / ega3d (E)  | `+1` | `+G = I` | `l_dual`              |
| sta4ds (Minkowski) | `-1` | `-G`     | the antidual          |

(In Euclidean signatures `r_undual` *is* `l_dual`, so `l_dual(r_dual(u)) == u` holds; in
`sta4ds` `r_undual` is the antidual and the `l_dual` composition fails.)

So the user-facing unification is conceptual: say **"`r_undual` undoes `r_dual`"**, not
"`l_dual` undoes `r_dual`". Nothing changes for Euclidean users (there `r_undual` *is*
`l_dual`); in `sta4ds` it becomes the antidual.

**PGA is the exception.** With `det = 0`, `G⁻¹` / `(1/det)` does not exist and `r_dual`
(the bulk dual) annihilates the weight, so it has no inverse — PGA must split into bulk +
weight (§11). There is no single *metric-dual* round trip there.

**The only round trip that is universal across all four** (ega, pga, sta) is the
metric-free complement pair, which round-trips through the complement rather than the
metric dual (so it carries no metric structure):

```text
l_cmpl(r_cmpl(u)) == u            (and r_cmpl(l_cmpl(u)) == u)
```

Summary: there is one metric-aware round trip for the non-degenerate algebras
(`r_undual ∘ r_dual = id`, differing only by `det`/`rG`), and one metric-free round trip
for everything including PGA (`l_cmpl ∘ r_cmpl = id`).

## 11. Bridge to projective STA (degenerate metric)

`G · rG = det · I` is the universal statement; the determinant decides the regime:

- **`sta4ds`** (non-degenerate): `det = -1 ≠ 0` => `rG = -G`, a global rescaling, no split.
- **projective STA** (degenerate, e.g. G(1,3,1)): adding one null direction makes
  `det = 0` => `G · rG = 0`. `G` and `rG` then have **disjoint support** — which is
  exactly the PGA **bulk / weight** split (`bulk_dual = cmpl(G·A)`,
  `weight_dual = cmpl(rG·A)`). The pseudoscalar-multiplication dual dies there (no
  `inv(I)`), but `complement ∘ metric` survives and simply splits in two.

So the `sta4ds` understanding directly informs the degenerate construction: the dual is
always `complement ∘ metric`; the only question is whether `G` and `rG` are proportional
(non-degenerate) or supported on complementary subspaces (degenerate).

## 12. Status in code / tests

Implemented and verified in [ga_test/src/ga_sta4ds_test.hpp](../ga_test/src/ga_sta4ds_test.hpp):

- forward dual `A ⟑ I == r_dual(A)`; `l_dual` vs `r_dual` even/odd split
- the naive `l_dual(r_dual(A)) == det · gr_inv(A)` obstruction and its `det · gr_inv` healing
- the metric / antimetric exomorphism relations (`G·rG = det·I`, `rG = -G`, wedge /
  antiwedge exomorphism laws, faithful `G = rev ∘ Gx`)
- left and right antidual (`r_undual` inverts `r_dual`, `l_undual` inverts `l_dual`)
- the wedge product uniquely selects `(-,-,-,+)` as the space, `(+,+,+,-)` as its dual
- the metric-independent `antidual` round trip, both directions, all grades

## 13. Open questions / next steps

### NEXT STEP (agreed) — the exomorphism comparison table

Build a comparison test, **one section per algebra** (ega2d, ega3d, pga2dp, pga3dp,
sta4ds), reading the stored extended-metric arrays from `ga_usr_consts.hpp` directly (not
a transcription), and tabulating for every basis blade `e_S`:

- `P(e_S) = prod_{i in S} m_i`     (pure product, computed from the signature)
- `Q(e_S) = sigma(k) * P`          (blade square)
- the stored extended-metric value `g_S`
- **multiplicativity verdict**: does `g_{S∪T} == g_S * g_T` for disjoint S,T?  (= is the
  stored metric a wedge exomorphism)
- the antimetric `rG_S = prod_{i not in S} m_i` (complement product), `g_S * rG_S == det`,
  and the antiwedge-multiplicativity.

Expected (to be confirmed): EGA/PGA stored arrays `== P`, multiplicative; STA stored array
`== Q`, multiplicativity **fails at grades 2,3**. This converts "bug vs convention" into a
mechanical, exhaustive, signature-only check, and doubles as the full exomorphism/
antiexomorphism tables for user documentation. Put it in each algebra's test file.

### The decision the table feeds: is `sta4ds_metric` a convention or a bug?

If a bug, the fix is to drop the `is_minkowski` override in `calculate_extended_metric` so
STA uses the pure product `P` like EGA/PGA. **Deciding experiment** (throwaway build, do
NOT change library behaviour without explicit sign-off): make that change, regenerate,
re-run the STA suite + transcription gate + contraction identities.

- core GA identities still hold → the special case was unnecessary (likely a sign bug);
  pure product is the fix and the whole `gr_inv` / antidual story disappears.
- core identities break → the blade square is load-bearing; learn *why* before changing.

### Later

- **Projective STA** (G(1,3,1) or G(3,1,1)): generate the rule tables via `ga_prdxpr`
  (likely needs a 5D dispatch arm in the product-case handlers; the rule generator
  itself looks dimension-general). Validate the complement / extended-metric /
  antimetric / bulk-dual / weight-dual tables with the transcription gate
  `dual(e) == nrm_sq(e) · cmpl(e)`.
- Decide whether to expose named `antidual` / `r_undual` / `l_undual` (and/or the faithful
  pure-product metric) in the library `sta4ds` ops, or keep them as documented test-only
  constructions.
- Confirm the signature-independence empirically once a `(+,+,+,-)` variant exists (only
  G(1,3,0) is instantiated today; the construction uses no signature-specific value
  beyond `Gx² = I`).
