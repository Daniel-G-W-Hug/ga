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

> **DECISION (2026-05-28), IMPLEMENTED (2026-05-29): P-unify.** The investigation (below)
> established that STA's stored extended metric *was* the **blade square** `Q`, injected by
> the `is_minkowski` override, whereas the canonical Lengyel **metric exomorphism** is the
> **pure product** `P = ∏_{i∈S} gᵢ` (seeded by the vector metric, propagated by
> `G(a∧b)=G(a)∧G(b)`; unique). The `Q` choice made the extended metric a *non*-exomorphism
> (comparison table `NO`/`NO`) and conflated two distinct objects. **Resolution (now in
> code): the extended metric is `P` for every algebra (`is_minkowski` removed), and
> `dot`/`nrm_sq` are defined from it (`P`, the reverse-norm), unifying STA with EGA/PGA.**
> The physics is preserved by re-sourcing the *causal character* and *rotor* decisions from
> the **geometric square** `B² = ⟨B B⟩₀ = gr0(X⟑X)` (which comes from the vector signature
> and never changed), not from `nrm_sq`.
>
> Under P-unify the dual is the clean exomorphism dual, so the **grade-dependent `gr_inv`**
> of §7 is gone — but it is replaced by a single **global determinant sign**, not by the
> identity:
>
> ```text
> l_dual(r_dual(A)) == r_dual(l_dual(A)) == det · A == -A     (det = -1 for STA)
> ```
>
> This equals `A` only for a Euclidean `det = +1` algebra (ega2d); for STA it is `-A`. The
> antidual still inverts `r_dual` exactly but is no longer *needed*. Likewise the old
> coincidence "`A⟑I == r_dual(A)` exactly" is gone: now `r_dual(A) == rev(A)⟑I` (differ by
> the reversion sign, exactly as in EGA — equal at grades 0,1,4, sign-flipped at 2,3).
> Corrected model, rationale and the (now-completed) **stepwise implementation plan** are in
> §13 and §14. *(This supersedes the interim 2026-05-27 "keep `Q`" reading: `Q` is the
> physical norm, but it is the **geometric square**, not the metric — so it belongs in the
> causal/rotor path, while the metric itself is the exomorphism `P`.)*

## 1. Two schools of dualization

There are two common ways to dualize a blade / multivector `A`:

1. **Complement after metric** (this library, after Lengyel): `dual(A) = cmpl(G·A)`,
   where `G` is the extended metric and `cmpl` is the purely combinatorial complement
   (defined by `A ^ cmpl(A) = I`). Split into a left and a right variant in even
   dimension (see §4).
2. **Pseudoscalar multiplication** (the "other school"): `to_dual(A) = A ⟑ I`,
   `from_dual(A) = A ⟑ inv(I)`.

In `sta4ds` the metric is non-degenerate, so both schools exist. Under P-unify they relate
by the reversion, exactly as in EGA: `r_dual(A) == rev(A) ⟑ I` at every grade. So the
plain `A ⟑ I` equals `r_dual(A)` only where the reversion sign `sigma(k) = +1` (grades
0,1,4) and equals `-r_dual(A)` at grades 2,3. *(Before P-unify, with the blade-square
metric `Q`, `A ⟑ I == r_dual(A)` held exactly at every grade — that coincidence was a
`Q`-artifact; see §6/§13.)*

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

## 6. The extended metric: pure product `P` (metric) vs blade square `Q` (geometric square) — RESOLVED

This was the crux of the investigation; it is now **settled and implemented (P-unify)**.
The two scalars `P` and `Q` below are genuinely **different objects** that coincide for
EGA/PGA but split at grades 2,3 in STA. The resolution: the **extended metric is `P`** (the
unique wedge exomorphism), for every algebra; the blade square `Q` is the **geometric
square** `gr0(X⟑X)` and belongs to the causal/rotor path, not the metric. STA used to store
`Q` (via the now-removed `is_minkowski` override); that was the bug. The section is kept as
the conceptual reference (P vs Q) that feeds `ga_docu`.

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

### What the stored arrays contain (after P-unify)

- `ega3d_metric` (and ega2d / pga): the **pure product** `P`. E.g. `e12 -> +1 = P(e12)`,
  although the blade square `Q(e12) = e12^2 = -1`. For Euclidean this is just the
  identity (all `+1`).
- `sta4ds_metric`: now also the **pure product** `P` — e.g. `g14 -> -1`, `g23 -> +1`,
  `I -> -1`. *(Before P-unify it held the blade square `Q`: `g14 -> +1`, `g23 -> -1`. The
  single difference was the reversion factor `sigma(k)`, nonzero only at grades 2,3.)*

So **every algebra now stores `P`**; `Q = sigma(k)*P` is no longer stored anywhere — it is
recomputed on demand as the geometric square `gr0(X⟑X)` where the physics needs it.

### Why `Q` cannot be the metric (it is not a wedge exomorphism)

A diagonal metric `g_S` is a genuine wedge exomorphism (Lengyel's
`G(a^b) = G(a)^G(b)`, the "conforming property" CLAUDE.md says the extended metric must
satisfy) **iff** it is multiplicative: `g_{S∪T} = g_S * g_T` for disjoint S,T.

- `P` is multiplicative by construction → **`P` is a wedge exomorphism** (so it is the
  metric).
- `Q = sigma(k)*P` is **not** multiplicative (`sigma(k+l) = (-1)^(kl) sigma(k) sigma(l)`)
  → `Q` is **not** a wedge exomorphism, so it cannot be the metric.

Confirmed in code (`ga_sta4ds_test.hpp`, test "metric / antimetric exomorphisms ... and
metric-indep. dual"): now that the library stores `P`, the library's own metric
`G(X) = l_cmpl(r_dual(X))` (which reproduces `sta4ds_metric` entry-for-entry, because
`r_dual = r_cmpl o G` and `l_cmpl o r_cmpl = id`) satisfies
`G(wdg(a,b)) == wdg(G(a),G(b))` at every grade incl. mixed. The exhaustive check is the
comparison-table test (§13), which now reads `YES`/`YES`. *(Before P-unify the library
metric was `Q` and this law failed: `Gx(wdg(g1,g4)) != wdg(Gx(g1),Gx(g4))`; only the
reversion-corrected `rev(Gx)` recovered `P`. That `rev` is exactly the `sigma(k)` that
distinguished `Q` from `P`.)*

### Why no existing test caught it (historical)

The extended metric is **architecturally decoupled** from the rest: the wedge product is
metric-free, the geometric product uses the vector signature `m_i` directly, and the
complements are purely combinatorial. The stored array fed **only** `dot` / `nrm_sq` and
the metric dual `cmpl(G·A)`. The transcription gate `dual(e) == nrm_sq(e)*cmpl(e)` and the
contraction identity `a<<b == rwdg(l_dual(a),b)` hold **by construction** for whatever the
metric is (they check internal consistency, not correctness), so a non-exomorphism metric
sailed through untouched. The gap is now closed by the **comparison-table test** (checks
multiplicativity directly, §13) and the **`detail::sta4ds_geom_sq` guard** (pins the
geometric square that feeds causal character to `gr0(X⟑X)`).

### Consequence and locus

Everything that looked strange in STA — the `gr_inv` in `l_dual(r_dual(u))` (§7), the need
for the antidual / `r_undual` (§8), the `A⟑I == r_dual` coincidence — was a direct
consequence of this one `sigma(k)`. With the pure product `P` (as EGA uses) those all
collapse: `r_dual(A) = rev(A)⟑I`, `l_dual(r_dual(A)) = det·A` (a clean global sign, no
`gr_inv`), and the antidual is no longer needed. The locus was narrow: the `is_minkowski`
override in
[calculate_extended_metric](../ga_prdxpr/src_prdxpr/rules/ga_prdxpr_rule_generator.cpp),
which overrode the product-of-vector-metrics with the blade-square values at grades 2,3
(its comment claimed to derive them "from `<A,~A>`", but it emitted `<A·A>_0`, the
geometric square — the opposite sign). The complement / product / dual **generation** was
independent and correct.

**Verdict — RESOLVED & IMPLEMENTED (2026-05-29).** The `is_minkowski` override is removed;
the extended metric is now the exomorphism `P` for every algebra, and the geometric square
`Q = gr0(X⟑X)` is recomputed where the physics needs it (causal character, rotors). The
investigation had confirmed two things: (a) `Q` is the correct *geometric square* / physical
norm, and (b) `Q` is **not** the metric exomorphism — the metric is uniquely `P`. The
override had put the geometric square into the *metric* slot (wrong object: non-exomorphism
metric, corrupted dual). See §13 for the corrected model, §14 for the (completed) plan. The
grade-3 signature-independence caveat below is now moot — the generic path computes
`P = ∏ gᵢ` directly from the signature, so it tracks any signature automatically.

## 7. The naive round trip, and why it fails

> **P-unify note (read first).** §7–§10 below were written for the **old `Q`-metric** and
> document its `det·gr_inv` obstruction and the antidual machinery that repaired it. Under
> P-unify (now implemented) the obstruction collapses to a **clean global sign**:
> `l_dual(r_dual(A)) == det · A == -A` at *every* grade — no `gr_inv`. The healing is just
> `det` (`A == det · l_dual(r_dual(A))`); the antidual still inverts `r_dual` exactly but is
> no longer required. The `gr_inv`/`det·gr_inv` formulas in §7 hold only for the superseded
> `Q`-metric; they are retained as the derivation that motivated P-unify. §8–§10 (antidual,
> signature duality, the unified `r_undual ∘ r_dual = id`) remain valid as written.

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

**P-unify implemented and verified (2026-05-29)** — ega 2583 / pga 2748 / sta 570 (28
cases) / ga_py 688 assertions pass; ega/pga rule-generator output byte-identical. In
[ga_test/src/ga_sta4ds_test.hpp](../ga_test/src/ga_sta4ds_test.hpp):

- `r_dual(A) == rev(A) ⟑ I` (the dual is pseudoscalar-mult of the reverse); `l_dual` vs
  `r_dual` even/odd split; `A ⟑ I == r_dual(A)` at grades 0,1,4 and `== -r_dual(A)` at 2,3
- the clean round trip `l_dual(r_dual(A)) == r_dual(l_dual(A)) == det · A == -A` (no
  `gr_inv`); healing by `det` alone (`A == det · l_dual(r_dual(A))`)
- the metric / antimetric exomorphism relations (`G·rG = det·I = -I`, `rG = -G`, wedge /
  antiwedge exomorphism laws — comparison table now `YES`/`YES`); the library metric
  `G = l_cmpl ∘ r_dual` is itself the exomorphism `P` (no `rev` correction needed)
- left and right antidual (`r_undual` inverts `r_dual`, `l_undual` inverts `l_dual`) — still
  exact, though no longer required for the round trip
- the wedge product uniquely selects `(-,-,-,+)` as the space, `(+,+,+,-)` as its dual
- causal character (`is_timelike`/`is_spacelike`/`is_lightlike`) and rotors sourced from the
  geometric square `gr0(X⟑X)` via `detail::sta4ds_geom_sq`, guarded by a dedicated test
- the transcription gate `dual(e) == nrm_sq(e) · cmpl(e)` holds with both sides `P`

## 13. Open questions / next steps

### DONE — the exomorphism comparison table

Implemented as a `TEST_CASE("... metric / antimetric exomorphisms (G, rG) - comparison
table")` in each algebra's test file (ega2d, ega3d, pga2dp, pga3dp, sta4ds). It reads the
stored extended-metric arrays from `ga_usr_consts.hpp` directly (not a transcription), and
tabulates for every basis blade `e_S`:

- `P(e_S) = prod_{i in S} m_i`     (pure product, computed from the signature)
- `Q(e_S) = sigma(k) * P`          (blade square)
- the stored extended-metric value `g_S`
- **multiplicativity verdict**: does `g_{S∪T} == g_S * g_T` for disjoint S,T?  (= is the
  stored metric a wedge exomorphism)
- the antimetric `rG_S = prod_{i not in S} m_i` (complement product), `g_S * rG_S == det`,
  and the antiwedge-multiplicativity.

Confirmed by the implemented tables: EGA/PGA stored arrays `== P`, multiplicative, both
verdicts `YES`; STA stored array `== Q`, multiplicativity **fails at grades 2,3**, both
verdicts `NO`. Each table also prints the antimetric `rG_S = g_{S^c}` and `g_S * rG_S`
(`== det` for EGA/PGA; for STA it alternates `sigma(k)sigma(n-k)*det` by grade — a second,
independent fingerprint of the non-multiplicativity). The tables double as the full
exomorphism / antiexomorphism reference for user documentation.

### Investigation (2a/2b): root cause and the load-bearing test

Two follow-up experiments localised the cause and characterised `Q`.

**2a — root-cause localisation (read-only, rule-generator dump).**

- The generator emits the same `Q` diagonal as `ga_usr_consts.hpp` (faithful copy), and
  the extended-metric matrix is purely diagonal (no off-diagonal coupling).
- The wedge product table is **metric-blind**: it carries only `0`/`1` literals, every
  sign living inside the blade name (orientation from reordering). It is also the same
  dimension-general 4D code as pga3dp (which passes). Combined with the existing test that
  shows `P` *does* satisfy the exomorphism with the current wedge, **wedge generation is
  ruled out** as the cause — for both `wdg` and `rwdg`.
- `calculate_extended_metric_recursive` is **not** an independent `P` source: it extracts
  `dot(e,e)`, which equals the stored metric (circular). It reproduces `Q`, and is dormant
  (never called).
- Locus confirmed: the `is_minkowski` override (lines ~277–291) overwrites the
  already-computed pure product `P` with the blade square `Q` at grades 2,3.

**2b — "is `Q` load-bearing?" (throwaway: flipped `nrm_sq` grades 2,3 from `Q` to `P`,
rebuilt, ran the STA suite, then reverted).** 54 assertions / 4 cases failed, splitting
cleanly into:

- **genuine physics** — `is_timelike`/`is_spacelike` of *every* bivector and trivector
  invert; the rotor `exp()` `cos`/`sin`/`cosh`/`sinh` structure breaks; normalization
  signs flip. The other 24 cases (wedge, gpr, complements, contraction, transform
  sandwich, spacetime split, projections, reflections) were untouched.
- **coupling artifacts** (expected, since only `nrm_sq` was flipped) — `dot(u,u)==nrm_sq`
  and the transcription gate `dual==nrm_sq·cmpl` fail, confirming `nrm_sq`/`dot`/`dual`
  are welded to one metric.

**What 2b actually shows.** The `Q` values are the **geometric square** `⟨B B⟩₀ = B²`
(the blade times itself, via the geometric product / vector signature). That quantity *is*
the physically meaningful thing for causal character (boost bivectors `B²>0`, rotation
bivectors `B²<0`) and for the rotor `cos`/`cosh` split. 2b broke because the library
computes those *from `nrm_sq`* and `nrm_sq` was defined as `Q`; flipping `nrm_sq` alone (to
`P`) without re-sourcing the causal/rotor code is what inverted the labels. **`Q` is real
and load-bearing — but it is the geometric square, not the metric.** The metric exomorphism
is uniquely `P` (§ refs to Lengyel below). The `is_minkowski` override conflated the two by
writing the geometric square into the *metric* slot.

### Decision: P-unify (2026-05-28)

Corrected model (confirmed against Lengyel's framework — the one this codebase follows):

- **The extended metric is the metric exomorphism** `G`, defined by extending the vector
  metric `𝔤`: `G(eᵢ)=gᵢ` (grade-1 seed), propagated by `G(a∧b)=G(a)∧G(b)`, giving
  `G(e_S)=∏_{i∈S} gᵢ = P`. It is **unique** ("there is exactly one valid inner product
  on the full exterior algebra", Lengyel). It is *not* metric-ignorant — the signature
  flows in at grade 1; EGA only looks metric-free because every `gᵢ=+1`.
- The antimetric `rG` (`rG(e_S)=∏_{i∉S} gᵢ`) and `G·rG = det(𝔤)·I` hold for **every**
  algebra; `det(𝔤)` is the **vector**-metric determinant (`−1` for STA).
- The metric is **independent of blade squares**. PGA3DP's extended metric is `1`/`0`
  (exomorphism / bulk-weight split) even though its bivectors square to `−1`. Likewise STA's
  metric is `P`, independent of the geometric squares `Q`.
- The exomorphism metric **defines `dot` and hence `nrm_sq`** (`s·s==1`; the rest by
  the exomorphism). This is `P` (the reverse-norm `⟨A Ã⟩₀`), consistent across all
  algebras — EGA already does this (`nrm_sq(BiVec3d)=+1=P`, not the geometric square `−1`).

**Chosen resolution:**

1. Remove the `is_minkowski` override → the extended metric is `P` for every algebra
   (the generic product-of-vector-metrics path already computes it).
2. Regenerate the metric-derived sta4ds products against `P`: `dot`, `nrm_sq`, the
   **dual** (`cmpl(P·A)` → clean exomorphism dual, `l_dual(r_dual) == det·A == -A`, *not*
   `== A` — that holds only for a Euclidean `det = +1` algebra), and the **contractions**
   `l_contract`/`r_contract` (they are metric-derived via the dual: `a<<b = rwdg(l_dual(a),
   b)`). Expansions (`l_expand`/`r_expand`) follow automatically (thin `wdg(l_dual, ·)`
   wrappers).
3. **Preserve the physics** by re-sourcing the *causal character* (`is_timelike` /
   `is_spacelike` / `is_lightlike`) and the *rotor* `exp`/`sqrt`/`transform` decisions from
   the **geometric square** `B² = gr0(B*B)` (vector-signature based, unchanged), not from
   `nrm_sq`. Net effect: causal character and rotor behaviour are **unchanged**; only the
   *value* of `nrm_sq` changes (`Q → P`, the reverse-norm).
4. After: the exomorphism comparison table reads `YES`/`YES`; `G·rG=det·I` holds; the
   transcription gate `dual==nrm_sq·cmpl` holds with both sides `P`; the `gr_inv`/antidual
   apparatus of §7–8 is no longer required (it was a `Q`-artifact).

**NOT in scope / explicitly unchanged:** the **geometric product** `gpr`, `wdg`, `cmpl`,
`rwdg`, `cmt` — these read the **vector signature** directly (`multiply_basis_elements →
vector_metric_value`), never the extended metric, so they are already correct and must NOT
be regenerated (the output would be byte-identical). The earlier worry "all geometric
products are wrong" does not hold.

The full stepwise plan is in §14. The earlier "keep `Q` / two-metric split (option C)" idea
is rejected in favour of full P-unify (one consistent `P` metric everywhere, with the
geometric square `Q` used only where it physically belongs).

### Grade-3 signature-independence caveat (moot under P-unify; kept for the record)

> **Superseded by P-unify:** removing the override (§13 Decision, §14 P1.1) makes the
> metric `P = ∏ gᵢ` straight from the signature, so there is no hard-coded grade-3 value
> left to be signature-specific. The analysis below explains *why* the old override was
> additionally fragile, and is retained as rationale only.

Separate from the convention question, the `is_minkowski` override's comment (lines
278–281) claims its hard-coded grade-2/3 values are *"independent of which signature
convention"*. That holds at **grade 2 only**: a convention switch `(−,−,−,+) ↔ (+,+,+,−)`
is a global `m_i → −m_i`, which scales a grade-`k` blade square by `(−1)^k` — invariant
for even `k`, sign-flipping for odd `k`. So:

- grade 2 (even): `Q(bivec)` is genuinely convention-independent → the hard-coded
  `has_time ? +1 : -1` is correct in both signatures. ✓
- grade 3 (odd): `Q(trivec)` **flips** with the convention, yet the override hard-codes a
  fixed `has_time ? -1 : +1`. It is therefore correct only for the instantiated `G(1,3,0)`
  `(−,−,−,+)`; under `G(3,1,0)` `(+,+,+,−)` the grade-3 metric would come out wrong.
- grade 1 is safe because it is *not* overridden — it falls through to the
  product-of-vector-metrics, which tracks the convention automatically.

This also corrects the optimistic note below that the construction "uses no signature-
specific value beyond `Gx²=I`": the override *does* bake in a signature-specific grade-3
value. Fixing it (e.g. computing `Q` from `sigma(k)·P` rather than hard-coding) is a
prerequisite for the projective-STA / `(+,+,+,-)` work.

### Later

- **Code-generate the complements and the duals in `ga_prdxpr` (`--output=code`).** Today
  `l_cmpl`/`r_cmpl` and `l_dual`/`r_dual` are **hand-transcribed** in `*_ops_basics.hpp`
  from the `ga_prdxpr_rule_generator_test` dual/complement tables — a repeated source of
  sign bugs (the bivector-dual slip that motivated this whole investigation; a fresh
  hand-edit of the 4 STA duals was needed again in P2.1). The rules and the extended metric
  already exist in the generator, so emitting C++ should be straightforward:
  - The combinatorial complement is `l_cmpl`/`r_cmpl` (basis order only); the dual is
    `cmpl(G·A)`. Both are linear per-grade maps → a flat-constructor emitter like the other
    primitives. Add `l_cmpl`/`r_cmpl`/`l_dual`/`r_dual` (and PGA `bulk_dual`/`weight_dual`)
    to the emitter's `product_to_cpp_function` map and wire the per-grade component layout.
  - Then `splice_generated_code.py` can keep `*_ops_basics.hpp` in sync exactly as it now
    does for `dot`/contractions, and the hand-transcription step (and its bug class) goes
    away. Guard with the existing transcription gate + the comparison-table + `geom_sq`
    tests.
  - Bonus: this also removes the only remaining hand-coded metric-derived ops, so a metric
    or signature change becomes a pure "regenerate + splice" operation.
- **Projective STA** (G(1,3,1) or G(3,1,1)): generate the rule tables via `ga_prdxpr`
  (likely needs a 5D dispatch arm in the product-case handlers; the rule generator
  itself looks dimension-general). Validate the complement / extended-metric /
  antimetric / bulk-dual / weight-dual tables with the transcription gate
  `dual(e) == nrm_sq(e) · cmpl(e)`.
- **DONE (2026-05-29):** `r_undual` / `l_undual` are now first-class library ops in
  [ga_sta4ds_ops_basics.hpp](../ga/ga_sta4ds_ops_basics.hpp) (8 overloads each, composed as
  `r_undual = l_cmpl∘r_dual∘l_cmpl`, `l_undual = r_cmpl∘l_dual∘r_cmpl`), auto-bound into the
  `ga_py` `sta` module. The test file uses them directly (no more local lambdas). *(Still
  open: whether to also expose the faithful pure-product metric as a named op — currently
  it is just `l_cmpl∘r_dual`, used inline in tests.)*
- Confirm signature-independence empirically once a `(+,+,+,-)` variant exists (only
  G(1,3,0) instantiated today). After P-unify the metric is `P = ∏ gᵢ` straight from the
  signature, so it tracks any signature automatically — the grade-3 caveat is moot.
- Feed the corrected §13 model into the `ga_docu` STA chapter (the conceptual content,
  not the §14 implementation steps).

## 14. Implementation plan (P-unify rework) — DONE (2026-05-29)

Two parts, each with validation gates. Part 1 reworks the generator (`ga_prdxpr`) and
proves no regression for EGA/PGA. Part 2 rebuilds `ga_sta4ds_ops*.hpp` stepwise. Nothing
in Part 1 changes EGA/PGA output; the only algebra whose generated tables change is sta4ds.

> **Executed as written, with two refinements found during implementation:** (i) the dual
> round trip is `l_dual(r_dual) == det·id == -id`, not `== id` (the `gr_inv` is replaced by
> a single global `det` sign — see §7 note / §13.2); (ii) the metric-derived products that
> changed are `dot` **and both contractions** (`l_contract`/`r_contract`) — 61 functions
> total — regenerated via `ga_prdxpr --output=code` and spliced in by the new tool
> [splice_generated_code.py](../ga_prdxpr/src_prdxpr/utilities/splice_generated_code.py)
> (clang-format-matched). The duals and `nrm_sq` are hand-coded in `*_ops_basics.hpp` and
> were edited by hand from the rule-generator dual tables.

### Part 1 — ga_prdxpr generator rework + validation

**P1.1 — Remove the `is_minkowski` override.**
In `calculate_extended_metric`
([rules/ga_prdxpr_rule_generator.cpp](../ga_prdxpr/src_prdxpr/rules/ga_prdxpr_rule_generator.cpp),
~lines 202–291) delete the Minkowski detection (`is_minkowski`, `time_slot`) and the
grade-2/3 override block. The generic product-of-vector-metrics path
(`metric_value *= metric[slot]`) then yields `P = ∏ gᵢ` for every algebra, pseudoscalar =
`det`, scalar = `1`. Net: `calculate_extended_metric` becomes a pure exomorphism builder.

**P1.2 — Regression gate (EGA/PGA must be byte-identical).**
Build + run `ga_prdxpr_rule_generator_test`. Diff the extended metric (and, for PGA, the
regressive extended metric) for ega2d, ega3d, pga2dp, pga3dp against pre-change output.
They MUST be unchanged (those never used the override). This proves removal is safe.
*Gate: zero diff for the four non-Minkowski algebras.*

**P1.3 — New sta4ds metric + antimetric.**
Confirm the sta4ds extended-metric diagonal is now
`P = {1, -1,-1,-1,+1, -1,-1,-1,+1,+1,+1, +1,+1,+1,-1, -1}` (basis order
`1,g1,g2,g3,g4,g14,g24,g34,g23,g31,g12,g234,g314,g124,g123,g1234`) and the antimetric is
`rG = -P`. Verify `G·rG = det·I = -I` and `G(a∧b)=G(a)∧G(b)` for sample blade pairs.
*Gate: comparison-table verdicts (computed from the array) flip to `YES`/`YES`.*

**P1.4 — Enumerate and regenerate only the metric-derived sta4ds products.**
Metric-INDEPENDENT (do NOT touch): `gpr`, `wdg`, `cmpl`/`l_cmpl`/`r_cmpl`, `rwdg`, `cmt`
— they read `config.metric_signature` or are combinatorial.
Metric-DERIVED (regenerate against `P`): `dot` (`multiply_basis_elements_dot` reads
`extended_metric[index]`), and the duals (`generate_l_dual_rules`/`generate_r_dual_rules`
→ `calculate_dual_rules(config, G_data, …)`). Before regenerating, grep the generator for
every consumer of `calculate_extended_metric` / `G_data` to confirm the full list
(currently: dot, l_dual, r_dual, the metric-matrix display, `ga_prdxpr_metric_calc.cpp`).
*Gate: regenerated sta4ds `dot(e_S,e_S)` equals `P_S`; regenerated `r_dual`/`l_dual`
satisfy `l_dual(r_dual(u))==u` on basis blades.*

**P1.5 — Update the reference array.**
Set `sta4ds_metric` in
[ga/ga_usr_consts.hpp](../ga/ga_usr_consts.hpp) to the `P` diagonal from P1.3 (the array is
reference-only; it feeds the comparison-table test, not the library ops at runtime).
*Gate: the sta4ds comparison-table test reads `P`, no `[g_S != P!]` flags.*

### Part 2 — ga_sta4ds_ops*.hpp library rebuild (stepwise)

The library is hand-pasted generator output + hand-written norm/rotor code. Do the steps in
order; rebuild + run `ga_sta_test` after each so breakage is attributable. Expect specific
test groups to need updating — that is the point, not a regression.

**P2.0 — Decouple causal character + rotor decisions onto `B²` (do this FIRST).**
*Before* changing any metric values, re-source the sign-of-square decisions from the
geometric square so they survive the `nrm_sq` change:
- `is_timelike`/`is_spacelike`/`is_lightlike` (Vec/BiVec/TriVec): base the sign test on
  `B² = gr0(X*X)` (geometric product), not on `nrm_sq`. For a vector this equals today's
  value; for bi/trivectors it is the geometric square (current behaviour). Result must be
  unchanged at this step (still `Q`-equivalent), so the suite stays green.
- `exp`/`sqrt`/`transform` rotor helpers: where they branch rotation-vs-boost or compute a
  half-angle from `nrm_sq`, switch to `B²` (or `gr0(rev(R)*R)` for the versor norm, per the
  STA rotor notes in CLAUDE.md). *Gate: full STA suite still green with `nrm_sq` untouched.*

**P2.1 — Swap in the `P`-based duals.**
Replace `l_dual`/`r_dual` for `BiVec4ds` and `TriVec4ds` (grades 2,3 flip sign) with the
P1.4 output; cascade through the `MVec`/`MVec_E`/`MVec_U` aggregates. *Gate:
`l_dual(r_dual(u)) == det·u == -u` and `r_dual(l_dual(u)) == -u` at every grade (a clean
global sign — EGA-like up to `det`).*

**P2.2 — Swap in the `P`-based `dot` and contractions.**
Update `dot` (mv, mv_e, mv_u, bivec, trivec) and the contractions `operator<<`/`operator>>`
to the regenerated values (61 functions total: 5 `dot` + 28 `<<` + 28 `>>`). These are
standalone codegen output in `*_ops_products.hpp`, so use
`ga_prdxpr --output=code --products=dot,l_contract,r_contract` +
`splice_generated_code.py`. *Gate: `dot(u,u)` equals `P_S` on basis blades;
`library_coverage.py --algebra=sta4ds --diff` reports `0 differ`.*

**P2.3 — Redefine `nrm_sq` to `P` (the reverse-norm).**
Update `nrm_sq(BiVec4ds)`/`nrm_sq(TriVec4ds)` to match `dot(u,u)=P`; the aggregate
`nrm_sq(MVec*)` already delegates. *Gate: `value_t(dot(u,u))==nrm_sq(u)` holds again; the
transcription gate `dual(e)==nrm_sq(e)·cmpl(e)` holds with both sides `P`.*

**P2.4 — Update the tests that encode the old `Q` numeric values.**
- "metric signature (dot/gpr/nrm_sq)": `nrm_sq`/`dot` of bi/trivectors now `P` (signs flip
  vs old asserts) — update expected values.
- "causal character / nrm / normalization": causal-character *labels* should be
  UNCHANGED (P2.0 routed them through `B²`); only `nrm_sq`/`normalize` magnitudes' signs
  change — update those numeric asserts.
- comparison-table test: flip `CHECK(!wedge_exo)`/`CHECK(!antiwedge_exo)` →
  `CHECK(wedge_exo)`/`CHECK(antiwedge_exo)`, drop the two `!=` witnesses, and reword the
  comment (now `YES`/`YES`, consistent with EGA/PGA).
- §7–8 dualization tests: the `det·gr_inv` obstruction is replaced by the clean global sign
  `l_dual(r_dual(A)) == det·A == -A`; the antidual/`r_undual` constructions still hold but
  are no longer needed for the round trip. Replace the `gr_inv`-correction cases with the
  clean `det·A` assertions (and the EGA-like `A⟑I == r_dual(rev(A))`, i.e. `A⟑I` flips
  vs `r_dual` at grades 2,3).

**P2.5 — Sweep the remaining metric-dependent ops.**
Audit `transform` sandwich, spacetime split, projections/rejections, reflections,
angle/rapidity for any direct use of `nrm_sq` sign that should now be `B²`; verify against
the suite. *Gate: full `ga_sta_test` green.*

**P2.6 — Python bindings cross-check.**
If the sta module is built, regenerate the cross-check JSON and run `ga_py` sta tests
(`regenerate_python_test_data`; the `sta` submodule), since `dot`/`nrm_sq`/`dual` values
changed. *Gate: `ga_py` sta tests green.*

### Final acceptance

- ega2d/ega3d/pga2dp/pga3dp generator output byte-identical (P1.2). ✓
- sta4ds exomorphism comparison table `YES`/`YES`; `G·rG=-I`. ✓
- clean exomorphism dual `l_dual(r_dual) == det·id == -id` (no `gr_inv`); transcription gate
  holds (`P`). ✓
- causal character + rotor behaviour physically unchanged (sourced from `B²`). ✓
- full `ga_ega_test`, `ga_pga_test`, `ga_sta_test` (+ `ga_py`) green. ✓

### Risk notes

- The only *behavioural* change visible to a library user is the **value** of
  `nrm_sq`/`dot` on bi/trivectors (now the reverse-norm `P`) and the **dual** (now the
  clean exomorphism dual). Causal-character predicates and rotors are unchanged by design.
- Keep `gpr`/`wdg`/`cmpl`/`rwdg`/`cmt` untouched; regenerating them is unnecessary and
  risks spurious diffs.
- Do P2.0 before P2.3, or the suite will show false rotor/causal breakage mid-rework.

## 15. Dependency graph: what depends on what

The single root is the **vector signature** `𝔤 = (g1²,g2²,g3²,g4²) = (-1,-1,-1,+1)` (plus
the basis and its canonical order). From it two branches grow, and **keeping them separate
is the whole point of P-unify** — they were conflated by the old `Q`-metric:

- the **extended-metric branch** `G` — the wedge exomorphism seeded by `𝔤`; feeds
  `dot`/`nrm_sq`, the duals, and (via the duals) the contractions/expansions;
- the **geometric-product branch** — `gpr` feeds the geometric square `B² = gr0(X⟑X)`,
  which feeds causal character and the rotors.

```text
vector signature 𝔤  +  basis & canonical order
        │
        ├─ combinatorial / signature-direct (NOT via the extended metric)
        │     wedge ∧          e_S ∧ e_T = ±e_{S∪T} (disjoint), 0 (repeat)      [metric-blind]
        │     complement cmpl  A ∧ r_cmpl(A) = I ; l_cmpl∘r_cmpl = id            [metric-blind]
        │     rwdg             rwdg(a,b) = l_cmpl(wdg(r_cmpl a, r_cmpl b))        [metric-blind]
        │     geometric ⟑      gᵢ⟑gᵢ = gᵢ ; reorder + repeated-index metric      [signature]
        │
        ├─ EXTENDED-METRIC branch   G = wedge-exomorphism of ∧, seeded by 𝔤
        │     extended metric G   G(eᵢ)=gᵢ, G(a∧b)=G(a)∧G(b)  ⟹  G(e_S)=∏_{i∈S} gᵢ = P
        │     antimetric rG       rG(e_S)=∏_{i∉S} gᵢ ; G·rG = det·I ; (STA: rG = -G)
        │        ├─ dot / nrm_sq  dot(e_S,e_S) = G(e_S) = P ; nrm_sq(X) = ⟨dot(X,X)⟩₀
        │        └─ l/r_dual      dual(A) = cmpl(G·A)   [= rev(A)⟑I for non-degenerate]
        │              ├─ contraction <<,>>   a<<b = rwdg(l_dual a, b) ; a>>b = rwdg(a, r_dual b)
        │              ├─ expansion l/r_expand wdg(l_dual a, b) ; wdg(a, r_dual b)
        │              └─ round trip   l_dual(r_dual(A)) = det·A ; antidual = exact inverse
        │
        └─ GEOMETRIC-PRODUCT branch   (causal / rotor physics — NOT the metric)
              geometric square  B² = gr0(X⟑X) = σ(k)·P        [detail::sta4ds_geom_sq]
                 ├─ causal char   is_timelike/spacelike/lightlike = sign(B²)
                 ├─ exp(B)        branch on causal char ; |a| = nrm(B) = √|nrm_sq(B)|
                 ├─ sqrt(R)       versor norm √(gr0(rev(R)⟑R))
                 ├─ transform     R⟑X⟑rev(R)
                 └─ inv, proj/refl  gpr  (+ nrm_sq only for scalar/vector inv)
```

### Per-operation table

| operation | built from | metric-derived? | lives in code as |
| --- | --- | --- | --- |
| `gpr` ⟑, `wdg` ∧, `cmt` | signature + basis order | no (signature-direct / combinatorial) | codegen (`--output=code`) |
| `l_cmpl`/`r_cmpl`, `rwdg` | basis order | no (metric-blind) | hand-coded (← codegen target, §13 Later) |
| extended metric `G`, antimetric `rG` | signature, via the ∧-exomorphism | — (it *is* the metric) | `calculate_extended_metric` (generator) + `sta4ds_metric` (reference array) |
| `dot`, `nrm_sq` | `G` | **yes** | `dot`: codegen; `nrm_sq`: hand (must equal `dot` diagonal) |
| `l_dual`/`r_dual` | `G` + `cmpl` | **yes** | hand-coded from rule-gen dual tables (← codegen target, §13 Later) |
| contractions `<<`/`>>` | dual + `rwdg` | **yes** (via dual) | codegen |
| expansions `l_expand`/`r_expand` | dual + `wdg` | **yes** (via dual) | hand wrappers (auto-update from dual) |
| geometric square `B²` | `gpr` (signature) | no (not the metric) | `detail::sta4ds_geom_sq` (hand; `== gr0(X⟑X)`, guarded by a test) |
| causal char `is_*like` | `B²` | no | hand |
| `exp`/`sqrt`/`get_rotor`/`get_boost` | `gpr` + causal char + `nrm` | no | hand |
| `transform` sandwich | `gpr` | no | hand + codegen sandwich matrices |
| `inv` | `gpr` (+`nrm_sq` for scalar/vector) | partly | hand (Hitzer-Sangwine for grades ≥ 2) |
| `nrm` | `sqrt(abs(nrm_sq))` (magnitude only) | no (sign-agnostic) | hand |

### When you change X, re-check Y

- **vector signature** (e.g. `(−,−,−,+)` ↔ `(+,+,+,−)`): the root — *both* branches move.
  Regenerate the extended metric, `dot`, duals, contractions; recompute the geometric square
  / causal character / rotors; `nrm_sq` flips on bi/trivectors. (Causal-character *labels*
  stay if sourced from `B²`, since `det` and `B²` are signature-parity invariant.)
- **extended metric only** (e.g. add a null direction → projective STA, `det = 0`):
  `dot`, the dual (which then *splits* into `bulk_dual`/`weight_dual`, §11), and the
  contractions change; the geometric-product branch (causal character, rotors) is untouched.
- **geometric product / reordering convention**: `gpr`, the geometric square, rotors,
  `transform`, `inv` — **and** the extended metric (it is seeded by the signature through
  the same reordering).
- **complement tables**: `cmpl`, the dual (`= cmpl∘G`), `rwdg`, and the contractions.
