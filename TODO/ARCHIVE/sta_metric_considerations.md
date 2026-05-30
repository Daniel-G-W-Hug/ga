# STA4ds metric, dualization, and the metric / antimetric duality

Reference notes on the metric structure of `sta4ds` = G(1,3,0): how the extended metric
is built, why it is a different object from the geometric (blade) square, how
dualization works, the metric / antimetric duality `G · rG = det·I`, and the bridge to a
future projective spacetime algebra. These notes are the conceptual core intended to feed
the `ga_docu` STA chapter. Every statement is exercised in the sta4ds test suite
([ga_test/src/ga_sta4ds_test.hpp](../ga_test/src/ga_sta4ds_test.hpp)), with the EGA / PGA
analogues in the other algebras' test files.

Notation: `g_i = e_i²` is the **vector metric** (the signature). `G` is the **extended
metric**, `rG` the **antimetric**, `cmpl` the (metric-free) **complement**, `⟑` the
geometric product, `∧` the wedge, `rwdg` the regressive wedge (antiwedge), `I` the
pseudoscalar.

## 1. The extended metric is a wedge exomorphism

The only metric input to an algebra is the **vector metric** `𝔤` — the signature
`(g1², g2², g3², g4²) = (−1, −1, −1, +1)` for sta4ds. The inner product of two vectors
comes from it directly.

The **extended metric** `G` lifts `𝔤` to every grade as the unique **wedge exomorphism**:

```text
G(e_i)    = g_i              (grade-1 seed = the signature)
G(a ∧ b)  = G(a) ∧ G(b)      (propagation through the wedge)
  ⟹   G(e_S) = ∏_{i∈S} g_i   (product of the constituent vector metrics)
```

Write this product `P(e_S) = ∏_{i∈S} g_i`. The scalar maps to `1` (`s · s = 1`), the
pseudoscalar to `det(𝔤)`, and there is exactly one such inner product on the whole
exterior algebra.

**The extended metric is independent of how a blade squares.** This is the central point
and the easiest to get wrong. The extended metric is a property of the *wedge* structure
seeded by the signature — not of the geometric product. The clearest evidence is PGA:
pga3dp's extended-metric diagonal is `1` / `0` (the exomorphism value, with `0` marking
the degenerate bulk/weight split) even though its bivectors square to `−1` under the
geometric product. The metric value and the blade-times-itself are simply two different
objects.

A common pitfall is to push physical *characteristics* (spacelike / timelike) into the
extended metric. That is wrong: causal character lives in the geometric square (§3), not
the metric. The extended metric encodes only the wedge mechanics seeded by the signature.

Finally, the extended metric **defines** `dot` and hence `nrm_sq` at every grade:
`dot(e_S, e_S) = G(e_S) = P`. This scalar `P` is the blade's **reverse-norm** — the scalar
part of the blade times its *reverse*, `⟨rev(A) ⟑ A⟩₀`. §3 explains the reverse-norm and
why it (not the geometric square) is the quantity a metric is built from.

## 2. The metric `G` and antimetric `rG`

Following Lengyel, the metric carries two exomorphisms:

- the metric `G`: a **wedge** exomorphism, `G(a ∧ b) = G(a) ∧ G(b)`;
- the antimetric `rG`: an **antiwedge** exomorphism, `rG(rwdg(a, b)) = rwdg(rG(a), rG(b))`,
  with `rG(e_S) = ∏_{i∉S} g_i` (the complementary product).

They are related, for **every** algebra, by the ordinary matrix product

```text
G · rG = det(𝔤) · I        (I = identity matrix)
```

For sta4ds `det(𝔤) = (−1)(−1)(−1)(+1) = −1`, and every basis blade squares to ±1 so
`G² = I`; the relation collapses to `rG = −G` (the antimetric is the negated metric).
`det(𝔤)` is signature-parity invariant (`det = (−1)^q`, `q` = number of negative
directions), so the physically-equivalent signatures `(−,−,−,+)` and `(+,+,+,−)` both give
`−1`.

## 3. The geometric square `B²` — causal character, distinct from the metric

The **geometric square** of a blade is `B² = ⟨X ⟑ X⟩₀ = gr0(X ⟑ X)` — the blade times
itself under the geometric product. It relates to the metric `P` only by the reversion
sign:

```text
B²(e_S) = σ(k) · P(e_S),    σ(k) = (−1)^{k(k−1)/2}   (k = grade)
```

So `B²` equals the metric `P` at grades `k ≡ 0, 1 (mod 4)` and is its negative at grades
2, 3. Worked examples in sta4ds (`g1² = g2² = g3² = −1`, `g4² = +1`):

```text
g14:  P = g1·g4 = −1      B² = σ(2)·P = +1   ( = g14 ⟑ g14 )
g23:  P = g2·g3 = +1      B² = σ(2)·P = −1   ( = g23 ⟑ g23 )
```

**The metric uses the *reverse*-norm, which removes that `σ(k)`.** The reverse `Ã = rev(A)`
reverses the order of the vector factors in each blade; on a grade-`k` blade that is just
the sign `rev(e_S) = σ(k)·e_S`. The metric value `P` is the scalar part of a blade times
its reverse:

```text
nrm_sq(e_S) = ⟨rev(e_S) ⟑ e_S⟩₀ = σ(k)·⟨e_S ⟑ e_S⟩₀ = σ(k)·B² = σ(k)·σ(k)·P = P
```

because `σ(k)² = 1`. So a blade has two natural "squares" that differ *exactly* by the
reversion sign: the geometric square `B² = σ(k)·P` (blade ⟑ blade) and the reverse-norm
`P` (rev(blade) ⟑ blade). The reverse-norm is the one that behaves like a norm — e.g. in
**EGA3D** (Euclidean) the bivector `e23` has geometric square `e23 ⟑ e23 = −1` but
reverse-norm `rev(e23) ⟑ e23 = (−e23) ⟑ e23 = +1`, so `nrm_sq(e23) = +1` (a unit bivector
with positive norm). The metric / `dot` / `nrm_sq` use the reverse-norm `P`; causal
character uses the geometric square `B²` (below).

`B²` is the physically meaningful quantity: a bivector is a **boost** plane when `B² > 0`
and a **rotation** plane when `B² < 0`, and the rotor exponential splits on its sign
(`cosh / sinh` vs `cos / sin`). Causal character — `is_timelike`, `is_spacelike`,
`is_lightlike` — therefore reads `B²`, **not** `nrm_sq`. The latter is the metric
reverse-norm `P`, which carries the opposite sign at grades 2, 3. Keeping `B²` and the
metric `P` apart is exactly what avoids the §1 pitfall of encoding "timelike" into the
metric.

## 4. Two schools of dualization

There are two common ways to dualize a blade / multivector `A`:

1. **complement after metric** (this library, after Lengyel): `dual(A) = cmpl(G · A)` —
   multiply by the extended metric, then take the purely combinatorial complement (defined
   by `A ∧ cmpl(A) = I`). Splits into a left and a right variant in even dimension (§5).
2. **pseudoscalar multiplication**: `to_dual(A) = A ⟑ I`, `from_dual(A) = A ⟑ inv(I)`.

In a non-degenerate algebra the two relate by reversion: `r_dual(A) = rev(A) ⟑ I`. So
plain `A ⟑ I` equals `r_dual(A)` where `σ(k) = +1` (grades 0, 1, 4) and `−r_dual(A)` at
grades 2, 3 — exactly the EGA behaviour.

Pseudoscalar multiplication needs `I` to be invertible (`sta4ds`: `I ⟑ I = −1`, so
`inv(I) = −I`). In a **degenerate** algebra (PGA, and a future projective STA) a null
vector sits inside `I`, so `I ⟑ I = 0` and `I` has no inverse — the scheme is structurally
impossible there. `complement ∘ metric` survives because `I` never appears in it; it is
the robust definition.

## 5. Handedness: `l_dual` vs `r_dual` (even dimension)

`sta4ds` is 4-dimensional (even), so the left and right complements / duals differ on the
odd grades and agree on the even grades:

- `l_dual(A) = l_cmpl(G · A)`, `r_dual(A) = r_cmpl(G · A)`;
- `l_dual(A) = r_dual(A)` for grades 0, 2, 4; `l_dual(A) = −r_dual(A)` for grades 1, 3.

The non-metric complements are exact mutual inverses for any signature:
`l_cmpl(r_cmpl(A)) = r_cmpl(l_cmpl(A)) = A`.

## 6. The dual round trip

Composing the two metric duals yields a **clean global sign**, the determinant:

```text
l_dual(r_dual(A)) = r_dual(l_dual(A)) = det(𝔤) · A
```

For a Euclidean algebra (`det = +1`) this is the identity — `l_dual` undoes `r_dual`. For
sta4ds (`det = −1`) it is `−A`: `l_dual` is the inverse only up to the sign `det`.

The genuine inverse applies the metric on the **dual** side instead of the primal side.
Since `G² = I`, the inverse of `r_dual = r_cmpl ∘ G` is `G ∘ l_cmpl`, the **un-dual**
(Lengyel: antidual) — one per handedness:

```text
r_undual(D) = l_cmpl(r_dual(l_cmpl(D)))     inverts r_dual:  r_undual(r_dual(A)) = A
l_undual(D) = r_cmpl(l_dual(r_cmpl(D)))     inverts l_dual:  l_undual(l_dual(A)) = A
```

They coincide on the even grades and differ by a sign on the odd grades. Equivalently
`r_undual = (1/det)·rG·l_cmpl = det · l_dual` (so `= −l_dual` for sta4ds).

| algebra              | det  | `r_undual` is                          |
| -------------------- | ---- | -------------------------------------- |
| ega2d / ega3d (Eucl) | `+1` | `l_dual`                               |
| sta4ds (Minkowski)   | `−1` | the antidual (`−l_dual`)               |
| pga (degenerate)     | `0`  | — no inverse; bulk/weight split (§8)   |

So the convention-clean statement is "**`r_undual` undoes `r_dual`**": for Euclidean
signatures `r_undual` *is* `l_dual`; in sta4ds it is the antidual. There is also one round
trip valid for **every** algebra (including PGA): the metric-free complement pair
`l_cmpl(r_cmpl(A)) = A`, which round-trips through the complement and so carries no metric
structure.

## 7. Signature duality `(−,−,−,+)` ↔ `(+,+,+,−)`

Because `rG = −G`, the antimetric negates every vector square: `G` realises the signature
`(−,−,−,+)` while `rG` realises `(+,+,+,−)`. Each signature is an exomorphism of **exactly
one** product:

```text
G  = (−,−,−,+)   is the WEDGE     exomorphism   (NOT the antiwedge)
rG = (+,+,+,−)   is the ANTIWEDGE exomorphism   (NOT the wedge)
```

So, relative to the library's wedge convention, the **space** is the signature whose
metric respects the wedge product, `(−,−,−,+)`, and its **dual space** is `(+,+,+,−)`. The
whole structure is symmetric under `(∧ ↔ rwdg, G ↔ rG, space ↔ dual)`: if one takes the
regressive wedge as the primary product (working in the dual algebra), the roles swap.
This is precisely why the two signatures are physically equivalent — each is "the space"
for one of the two products, the other being its dual.

## 8. Bridge to projective STA (degenerate metric)

`G · rG = det · I` is the universal statement; the determinant decides the regime:

- **`sta4ds`** (non-degenerate, `det = −1 ≠ 0`): `rG = −G` — a global rescaling, no split.
- **projective STA** (degenerate, e.g. G(1,3,1), `det = 0`): `G · rG = 0`, so `G` and `rG`
  have **disjoint support** — exactly the PGA **bulk / weight** split
  (`bulk_dual = cmpl(G·A)`, `weight_dual = cmpl(rG·A)`). Pseudoscalar multiplication dies
  (no `inv(I)`), but `complement ∘ metric` survives and simply splits in two.

The dual is always `complement ∘ metric`; the only question is whether `G` and `rG` are
proportional (non-degenerate) or supported on complementary subspaces (degenerate).

## 9. Dependency map: what is built from what

The single root is the **vector signature** `𝔤` (plus the basis and its canonical order).
Two independent branches grow from it — and keeping them apart is the whole point of
§1 / §3:

```text
vector signature 𝔤  +  basis & canonical order
        │
        ├─ combinatorial / signature-direct (NOT via the extended metric)
        │     wedge ∧          metric-blind (orientation only)
        │     complement cmpl  metric-blind, A ∧ r_cmpl(A) = I, l_cmpl∘r_cmpl = id
        │     regressive rwdg  metric-blind, rwdg(a,b) = l_cmpl(wdg(r_cmpl a, r_cmpl b))
        │     geometric ⟑      uses the signature directly (g_i ⟑ g_i = g_i)
        │
        ├─ EXTENDED-METRIC branch   G = wedge exomorphism of ∧, seeded by 𝔤
        │     extended metric G   G(e_S) = ∏_{i∈S} g_i ;  antimetric rG ;  G·rG = det·I
        │        ├─ dot / nrm_sq  dot(e_S, e_S) = G(e_S)        (the reverse-norm)
        │        └─ l/r_dual      dual(A) = cmpl(G·A)           [= rev(A) ⟑ I]
        │              ├─ contraction <<, >>   a << b = rwdg(l_dual a, b)
        │              ├─ expansion            wdg(l_dual a, b)
        │              └─ round trip   l_dual(r_dual(A)) = det·A ;  un-dual = exact inverse
        │
        └─ GEOMETRIC-SQUARE branch   (causal / rotor physics — NOT the metric)
              geometric square  B² = gr0(X ⟑ X) = σ(k)·P
                 ├─ causal character   is_timelike / is_spacelike / is_lightlike = sign(B²)
                 ├─ rotors exp / sqrt  branch on B² ;  versor norm gr0(rev(R) ⟑ R)
                 └─ transform          R ⟑ X ⟑ rev(R)
```

What moves when something changes:

- **vector signature** (`(−,−,−,+) ↔ (+,+,+,−)`): both branches move — the extended
  metric, `dot`, the duals and contractions, and the geometric square / causal character /
  rotors. (Causal-character labels are signature-parity invariant, since both `det` and
  `B²` are.)
- **degeneracy** (add a null direction → projective STA, `det = 0`): the extended-metric
  branch changes (the dual splits into bulk/weight, §8); the geometric-square branch is
  untouched.
- **complement convention**: the complements, the dual (`cmpl ∘ G`), `rwdg`, and the
  contractions.
