# Choosing the Lift Anchor `z₀` (equivalently, the integer `k`)

**Scope.** This note gives a clean, defensible *argument* for how to pick the
auxiliary-variable anchor `z₀` in the lifting module (Manual §1b, §"Delta
Elimination", §"Worked Example"). It is a procedure + rationale only — **no
fitting code yet**, by request. It is meant to replace the current placeholder
heuristic in `DetectExtremeCoefficients`, which hardcodes `SuggestedK -> 1`
(`tropical_eval.wl:603`).

---

## 1. The setup, and why `z₀` is not actually a free knob

Lifting replaces an extreme coefficient `C` in a monomial `C·xᵅ` by `zᵏ·xᵅ` for
a new auxiliary variable `z`, and inserts `∫dz δ(z − z₀)`. For the lift to be an
*exact identity* (the round-trip check in `LiftCoefficients`, Manual §1b), the
anchor must satisfy

```
z₀ᵏ = |C|        ⇔        z₀ = |C|^(1/k).
```

So `z₀` and `k` are **locked together**. We do not get to choose `z₀` freely and
then pick `k`; choosing the integer `k` *determines* `z₀`, and conversely the
only way to move `z₀` while staying exact is to change `k`. Since `k` must be a
positive integer (it is the exponent of `z` in the lifted monomial and a vertex
coordinate of the lifted Newton polytope), the achievable anchors are the
discrete ladder

```
z₀ ∈ { |C|, |C|^(1/2), |C|^(1/3), … }   →  1   as  k → ∞.
```

**The real question is therefore: which integer `k` puts `z₀` in the best
place?**

---

## 2. What `k` controls — two distinct mechanisms

`k` enters the construction in two independent ways, and they pull in opposite
directions. Keeping them separate is the whole key to the argument.

### (a) Geometry mechanism — set by `k` as an *integer*, independent of `|C|`

The lifted monomial `zᵏ·xᵅ` contributes a Newton-polytope vertex at
auxiliary-exponent `k`. The integer `k` therefore fixes:

- the shape of the lifted Newton polytope → its normal fan → the simplicial
  cones, ray matrices `M`, and determinants `|det M|`;
- whether an admissible **pivot** exists, and whether the re-cleared sector
  polynomial keeps a constant term (`HasConstantTerm`, Manual §"Delta
  Elimination" step 5);
- the integrability margins `Re ãⱼ > 0` of each surviving sector.

None of this depends on the *magnitude* of `C` — only on the integer `k` and the
exponent vectors. This is the **gate**: it decides whether lifting produces
finite-variance, integrable sectors at all.

### (b) Anchor-magnitude mechanism — set by `z₀ = |C|^(1/k)`

After the pivot substitution, each monomial `{c, e}` in a surviving sector maps
to (`tropical_eval.wl:826`):

```
{ c · z₀^(e_p/m_p) ,  e_j − e_p·m_j/m_p } .
```

The deterministic prefactor absorbs `z₀^(a_p/m_p − 1)` (`tropical_eval.wl:995`)
— that carries the bulk of the extreme scale *out of the random integrand* and
into an exact constant, which is the whole point of lifting and is why the
dominant cone in the worked example becomes a flat `10⁻⁶·(1+y′)⁻²`.

But the **residual coefficient spread inside the sampled polynomial** is

```
spread  ≈  z₀^( Δe_p / m_p ),     Δe_p = (max − min) pivot-coordinate exponent,
```

which is `1` (perfectly flat) only when the pivot equalizes the `z₀`-powers.
In secondary cones it does not: e.g. cone 2 of the toy keeps
`Q̃ = z₀⁻¹ + y′ = 10⁻⁶ + y′`, a residual scale of `z₀⁻¹`, and that cone stays
heavy-tailed (`σ/μ ≈ 410`). The closer `z₀` is to `1`, the closer **all** such
residual powers are to `O(1)`, and the flatter every surviving sector.

> **Summary of the tension.**
> Mechanism (b) wants `z₀ → 1`, i.e. `k` **large**.
> Mechanism (a) wants `k` **small**: large `k` inflates `|det M|`, multiplies the
> cone count and `EmptyDomain` drops, shrinks the margins `Re ãⱼ` toward 0
> (heavy tails / `liftnopivot`), and can tip the polytope into degeneracy
> (`liftdegenerate`).

---

## 3. The clean target for `z₀`

We want `z₀` close to `1` for residual flatness, but only as far as the geometry
gate (a) stays healthy. The principled stopping point is **self-consistency with
the detector itself**:

> The package already declares a coefficient "extreme," and worth lifting, when
> it leaves the band `[1/τ, τ]` (`DetectExtremeCoefficients`, threshold `τ`,
> default `1000`, `tropical_eval.wl:597`). By the same standard, an anchor `z₀`
> is "tame" — no longer the kind of magnitude we set out to remove — exactly
> when `z₀ ∈ [1/τ, τ]`.

So: **pick the smallest integer `k` that brings the anchor back inside the
non-extreme band.** Going beyond that buys only sub-threshold residual
shrinkage (the detector would no longer even flag `z₀`) while paying rising
geometry cost — so the *smallest* qualifying `k` is the sweet spot, not the
largest feasible one.

```
        ⌈  |log|C||  ⌉            ⌈ log₁₀|C| ⌉
k*  =   |  ────────  |     =      | ──────── |        (with τ = 10³)
        |   log τ    |            |    3     |
```

clamped to `k* ≥ 1` (if `|C|` is already inside the band, no lift). The same
formula serves large *and* small coefficients, since it uses `|log|C||`.

### This reproduces every hand-picked `k` in the manual

| Case (Manual / benchmark)        |    `|C|`  | `k* = ⌈log₁₀|C| / 3⌉` | `z₀ = |C|^(1/k*)` | Hand-picked `k` | Outcome              |
|----------------------------------|-----------|-----------------------|-------------------|-----------------|----------------------|
| Case A (`1 + 10⁶x₁² + …`)         |   `10⁶`   |        `⌈2⌉ = 2`      |      `10³`        |       `2`       | 17.95× PASS          |
| Case B primary (`… + 10⁴x₁² …`)   |   `10⁴`   |       `⌈1.33⌉ = 2`    |      `10²`        |       `2`       | 47× reliable (6.5%)  |
| Example 20 (`10⁻⁴x₁² → z²`)       |   `10⁻⁴`  |       `⌈1.33⌉ = 2`    |      `10⁻²`       |       `2`       | matches manual       |

The agreement across large- and small-coefficient cases — and the fact that the
benchmarks *abandoned* `k=1` (the current hardcoded `SuggestedK`) in favor of
exactly these `k=2` choices — is the empirical support for the rule.

---

## 4. The procedure (decision order)

For a single dominant extreme coefficient `C` (the common case):

1. **Geometry gate first (mechanism a).** Form the lifted Newton polytope for the
   candidate `k` and check that at least one surviving sector has
   `HasConstantTerm -> True` with all margins `Re ãⱼ > 0`. If *no* `k` achieves
   this (degenerate polytope, all pivots lose the constant term), lifting cannot
   reduce variance regardless of `z₀` — **stop and report**, do not lift. (This
   is Case C: `k=2` and `k=4` gave *identical* failing results because the
   `HasConstantTerm=False` heavy tail dominates and is blind to the finite
   `z₀`-rescaling.)

2. **Anchor target (mechanism b).** Among `k` that pass the gate, take the
   smallest one with `z₀ = |C|^(1/k) ∈ [1/τ, τ]`, i.e.
   `k* = max(1, ⌈|log|C|| / log τ⌉)`.

3. **Tie-break toward simpler geometry.** Prefer the `k` that admits a `|m_p| = 1`
   pivot and the largest worst-case margin `min_j Re ãⱼ` (the same ranking
   `ProcessSectorLifted` already uses, Manual §"Delta Elimination" step 5). If
   `k*` from step 2 fails the gate but `k*+1` passes, step up; never step *down*
   below the band unless forced.

4. **Sanity floor/ceiling.** `k ≥ 1` always; in practice `k` rarely needs to
   exceed `⌈log₁₀|C|/3⌉`, so a small hard ceiling (e.g. `k ≤ 6`) is a reasonable
   guard against runaway geometry.

---

## 5. Multi-rule lifts (several extreme coefficients, one shared anchor)

When several monomials are lifted against the *same* `z₀`, each rule `i` carries
a residual `cᵢ = Cᵢ / z₀^{kᵢ}` (Manual §1b; `tropical_eval.wl:679`). Exactness
holds for any `kᵢ`, but the residuals are `O(1)` — i.e. the *cleared* polynomial
is well-balanced — only when

```
z₀^{kᵢ} ≈ |Cᵢ|     ⇔     z₀ ≈ |Cᵢ|^(1/kᵢ)   for every lifted i simultaneously.
```

So the single-rule rule generalizes: fix the shared `z₀` from the *primary*
(most extreme) coefficient via §3, then choose each secondary `kᵢ` so that
`|Cᵢ|^(1/kᵢ)` lands as close as possible to that shared `z₀`. If two extreme
coefficients cannot be matched to a common anchor within the band, prefer a
**single-rule lift on the dominant one** — the benchmark Case B confirms this:
the two-rule `k=1` lift reported a spectacular 2450× *sample* variance reduction
but a 64% true error (one sector `HasConstantTerm=False`, optimistic σ), while
the single-rule `k=2` lift gave an honest 47× with 6.5% deviation.

---

## 6. Caveats / when the argument doesn't apply

- **Degenerate lifted polytope.** If all pivots yield `HasConstantTerm=False`,
  the surviving sectors have (possibly infinite) true variance and the
  sample-based σ is optimistic; no choice of `z₀` fixes this (Case C). The gate
  in §4.1 must override the anchor target in §4.2.
- **Complex polynomial exponents `Bⱼ`.** Lifting requires real `Bⱼ` (`liftcomplex`);
  the procedure is moot otherwise.
- **Symbolic / kinematic coefficients.** `|C|` is unknown at lift time, so they
  are skipped by the detector and have no `z₀` to choose.
- **The prefactor never costs variance.** Note that the large
  `z₀^(a_p/m_p − 1)` prefactor is *deterministic* — it does not enter the MC
  variance. So "is `z₀` huge?" is not by itself bad; what matters is the
  *residual spread* `z₀^(Δe_p/m_p)` left inside the sampled polynomial, which is
  what §3 controls.

---

## 7. One-line statement

> `z₀` is forced to `|C|^(1/k)`, so choosing the anchor *is* choosing the integer
> `k`. Pick the **smallest `k` that pulls `z₀` back inside the detector's own
> non-extreme band `[1/τ, τ]`** — `k* = ⌈|log|C|| / log τ⌉` — provided that `k`
> passes the geometry gate (`HasConstantTerm`, positive margins). This minimizes
> the residual coefficient spread in the sampled integrand while keeping the
> lifted geometry as well-conditioned as possible, and it reproduces every
> hand-tuned `k` in the v2 benchmarks.
