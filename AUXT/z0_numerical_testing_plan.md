# Plan: Numerical Testing for the Optimal Lift Anchor `z₀`

**Companion to** [`anchor_selection_procedure.md`](anchor_selection_procedure.md).
That note derives the analytic rule; this note plans the **empirical test** that
either confirms it or maps where it breaks. Still **no fitting code** — this is
the experiment design. Code lands only after the harness and battery below are
agreed.

---

## 1. Hypothesis under test

From the companion note, with detector threshold `τ` (default `1000`):

> **H1 (anchor rule).** Subject to the geometry gate, the variance-minimizing
> integer is the smallest one that pulls the anchor into the non-extreme band:
> `k* = max(1, ⌈|log|C|| / log τ⌉)`, giving `z₀* = |C|^(1/k*)`.

> **H2 (gate dominance).** When *no* `k` yields a surviving sector with
> `HasConstantTerm = True`, lifting cannot reduce true variance for any `z₀`
> (the geometry gate overrides the anchor target).

> **H3 (monotone-then-flat).** Holding the gate fixed, true `σ/μ` of the
> dominant surviving cone decreases as `z₀ → 1` (increasing `k`) and then
> *plateaus* once `z₀` enters `[1/τ, τ]` — so larger `k` past `k*` buys
> negligible accuracy at rising geometry cost.

The experiment must **distinguish argmin over `k`** (the empirical optimum) from
**`k*`** (the rule's prediction) and report the gap.

---

## 2. Metrics (per `(case, k)` cell)

| Metric | Source | Why it matters |
|--------|--------|----------------|
| `σ_sample = stderr · √N` | `EvaluateTropicalMC[Lifted]` (cf. `bench_lift_variance.wl:73,135`) | the headline per-sample spread; **but optimistic when `HasConstantTerm=False`** |
| `σ_true` (per sector) | `trueSigma[sd, kin, pg]` (`sandbox_lift_common.wl`, driven as in `sandbox_true_sigma.wl`) | exact NIntegrate variance — the **decisive** figure; catches heavy tails the sampler misses |
| `VarRed = (σ_unlifted / σ_lifted)²` | ratio of the above | the actual goal of lifting |
| `relErr = |MC − truth| / |truth|` | vs. truth oracle (§4) | accuracy guard — a big `VarRed` with bad `relErr` is a false win (Case B 2-rule trap) |
| gate: `HasConstantTerm`, `min_j Re ãⱼ`, `|m_p|`, `#EmptyDomain` | `ProcessSectorLifted` output | the geometry gate of H2; recorded per surviving cone |
| `z₀`, `z₀^(Δe_p/m_p)` residual spread | `LiftData` + sector polys | links measured σ back to the anchor mechanism |

**Headline reduction must be computed from `σ_true`, not `σ_sample`.** This is
non-negotiable: the benchmark already shows `σ_sample` reporting a fake 2450×
(Case B) and a fake-flat 0.0824 (Case C) where `σ_true` is infinite.

---

## 3. The `k`-sweep harness (reuse, don't rebuild)

Generalize the existing `runLifted` helper (`bench_lift_variance.wl:104–210`)
into a loop over `k`:

```
sweepK[spec, primaryRule, kList, fanFor, N, seeds] :=
  for k in kList:
    rule_k   = primaryRule with "k" -> k
    lifted   = LiftCoefficients[spec, {rule_k}]          (* exact round-trip checked here *)
    fan      = fanFor[k]   (* Automatic, or explicit when liftdegenerate fires *)
    sectors  = ProcessSectorLifted[...] over all cones   (* gate metrics *)
    for seed in seeds:                                   (* seed-robustness, §6 *)
       mc   = EvaluateTropicalMCLifted[spec, {{}}, "LiftRules"->{rule_k},
                                       "FanData"->fan, "NumSamples"->N, "Seed"->seed]
    tsig   = trueSigma per surviving sector              (* σ_true *)
    record (k, z0, gate, σ_sample(seed-stats), σ_true, relErr)
```

Notes:
- `fanFor[k]` must handle `liftdegenerate`: when the lifted Newton polytope is
  lower-dimensional, supply the explicit fan exactly as Case C does
  (`benchmark_results.md`, the `{±k,0,∓3}` lineality pattern). The harness should
  *detect* the degeneracy and either load a precomputed fan or skip with a logged
  reason — never silently fall back to the unlifted pipeline.
- One **unlifted baseline** per case (the `k`-independent reference, via the
  unit-coefficient proxy fan, `bench_lift_variance.wl:32–101`) so every `VarRed`
  shares a denominator.
- Emit one row per `(case, k, seed)` to a CSV/`Dataset`, plus a per-`(case,k)`
  aggregate. Mirror the table format of `benchmark_results.md` for continuity.

---

## 4. Truth oracle (and the NIntegrate pitfall)

Accuracy needs a trusted value of the integral. **Known hazard (recorded in
project memory): direct `NIntegrate` is itself unreliable at extreme
coefficients** — the integrand's mass lives in a `1/|C|`-wide spike that adaptive
quadrature can miss just as the unlifted MC does. Mitigations, in priority order:

1. **Closed form where it exists.** Toy 0 (`∫(1+10⁶x)⁻² = 10⁻⁶`) and other 1-D /
   separable cases have exact answers — use them as anchor truths with zero
   oracle risk.
2. **Hardened `NIntegrate`.** High `WorkingPrecision`, raised
   `MaxRecursion`/`AdaptiveSampling`, and an explicit singular-point hint at the
   monomial-balance locus `xᵅ = 1/|C|`. Cross-check two methods (e.g.
   `"GlobalAdaptive"` vs `"MultiPeriodic"`); accept the oracle only if they agree.
3. **`ValidateLiftedDecomposition` as a consistency cross-check**, not as truth —
   it compares the lifted sector sum to direct `NIntegrate` of the original
   (`DirectResult` vs `SectorSum`), so it inherits oracle (2)'s risk but is useful
   for catching dropped-sector bookkeeping errors.
4. **Disagreement protocol.** If oracle methods disagree by `> 1%`, flag the case
   `oracle-untrusted` and rank `k` by `σ_true` *only* (which needs no global
   truth), reporting accuracy as indeterminate rather than guessing.

---

## 5. Test battery

A matrix that exercises both mechanisms and both rule edges (`⌈⌉` boundaries).

### Tier 1 — single extreme coefficient, non-degenerate polytope (validates H1, H3)
Vary magnitude across the `k*` step boundaries so we *see* the ceiling kick in:

| ID | polynomial | `|C|` | predicted `k*` (τ=10³) | sweep `k` |
|----|-----------|-------|------------------------|-----------|
| T1a | `1 + 10³x₁² + x₂² + x₁x₂²` | `10³` | 1 | 1–4 |
| T1b | `1 + 10⁴x₁² + x₂² + x₁x₂²` | `10⁴` | 2 | 1–4 |
| T1c | `1 + 10⁶x₁² + x₂² + x₁x₂²` (= Case A) | `10⁶` | 2 | 1–5 |
| T1d | `1 + 10⁷x₁² + x₂² + x₁x₂²` | `10⁷` | 3 | 1–5 |
| T1e | `1 + 10⁹x₁² + x₂² + x₁x₂²` | `10⁹` | 3 | 1–6 |

Expectation: empirical argmin-`k` (by `σ_true`, gate-passing) equals the `k*`
column, and `σ_true` plateaus for `k ≥ k*` (H3).

### Tier 2 — small coefficient (symmetry check)
`10⁻⁴ x₁² → z² …` (Example 20), plus `10⁻⁶`, `10⁻⁷`. Confirms the `|log|C||`
form handles `|C| < 1` identically (predicted `k* = 2, 2, 3`).

### Tier 3 — gate-dominated / degenerate (validates H2)
`1 + 10⁸ x₁³x₂ + x₂³` (= Case C). Sweep `k = 2,3,4`; expectation: **all** fail the
gate (`HasConstantTerm=False`), `σ_true` infinite/large for every `k`, so the
harness reports "do not lift" — and crucially the *finite* `σ_sample` numbers are
`k`-independent (Case C showed k=2 ≡ k=4), demonstrating the trap.

### Tier 4 — multi-rule / shared anchor (validates §5 of companion note)
Case B (`10⁻⁴ + 10⁴x₁² + 10⁻⁴x₂² + 10⁴x₁x₂² + x₁²x₂`). Compare single-rule `k*=2`
vs two-rule `k=1`, scoring on `σ_true` **and** `relErr`. Expectation: single-rule
`k=2` wins on the honest (true-σ, accuracy-gated) score even though two-rule `k=1`
wins on the misleading `σ_sample`.

### Tier 5 — dimension scaling (cost side of the tension)
Re-run T1c at `n = 3, 4` to measure how cone count, `#EmptyDomain`, and `|det M|`
grow with `k` — quantifying the geometry cost that justifies "smallest `k`, not
largest."

---

## 6. Controls and threats to validity

- **Seed robustness.** Run each MC at ≥ 5 seeds; report mean ± spread of
  `σ_sample`. A single seed is how the unlifted Toy 0 got a deceptively narrow
  bar (Manual §"The comparison"). `σ_true` is seed-free and is the tiebreaker.
- **`HasConstantTerm=False` ⇒ report `σ_true`, mark `σ_sample` optimistic.** Never
  rank on `σ_sample` for such sectors.
- **Fixed `N`** across all cells (e.g. `N = 10⁶`, matching `benchmark_results.md`)
  so `VarRed` is comparable; record wall-clock separately (it conflates
  compile + NIntegrate, per the benchmark caveat).
- **Exactness precondition.** Assert the `LiftCoefficients` round-trip
  (`liftidentity`) passes before trusting any row — a failed identity invalidates
  the whole cell.
- **One mechanism at a time.** Tier 1 isolates the anchor (fixed polytope family,
  only `|C|` and `k` move); Tier 5 isolates geometry cost. Don't co-vary.

---

## 7. Deliverables and decision criteria

1. **`AUXT/z0_sweep_results.md`** — benchmark-style tables, one block per case,
   columns `(k, z₀, gate, σ_sample, σ_true, VarRed, relErr, argmin?)`.
2. **Verdict on H1/H2/H3** — for each tier, does empirical argmin-`k` match `k*`?
   Tabulate `Δk = argmin_empirical − k*`; the rule **passes** if `Δk = 0` for all
   gate-passing Tier-1/2 cases and the Tier-3 cases correctly report "no lift."
3. **A falsification log** — every case where the rule misses, with the observed
   optimum and a hypothesis (e.g. a `|m_p|>1` pivot forcing fractional `z₀`
   powers, or a margin `Re ãⱼ→0⁺` heavy tail not captured by the band argument).
4. **Go/no-go for the fitting step.** Only if H1 holds (or holds with a small,
   characterized correction) do we wire `k*` into `DetectExtremeCoefficients`,
   replacing the hardcoded `SuggestedK -> 1` (`tropical_eval.wl:603`). If H1 needs
   a correction term, the falsification log defines what to fit — *then* the
   fitting procedure (deferred until now) gets written against this data.

---

## 8. Phasing

| Phase | Work | Exit condition |
|-------|------|----------------|
| P0 | Generalize `runLifted` → `sweepK`; wire `trueSigma`; CSV/`Dataset` output | sweep runs end-to-end on Case A, reproduces its 17.95× at `k=2` |
| P1 | Tier 1 (anchor rule, non-degenerate) | `Δk` table for T1a–T1e |
| P2 | Tier 2 + Tier 3 (small-coeff symmetry, gate dominance) | H2 confirmed on Case C; small-coeff symmetric to large |
| P3 | Tier 4 + Tier 5 (multi-rule, cost scaling) | honest single-vs-multi verdict; cost-vs-`k` curve |
| P4 | Synthesis: `z0_sweep_results.md`, H1/H2/H3 verdict, go/no-go | rule validated or correction characterized |

---

## 9. One-line statement

> Build a `k`-sweep over a magnitude-laddered battery, score every `(case,k)` on
> **true** (not sample) variance and on accuracy against a trusted oracle, and
> check whether the empirical argmin-`k` equals the analytic
> `k* = ⌈|log|C||/log τ⌉` under the geometry gate. The result is the go/no-go —
> and, if needed, the dataset — for replacing `SuggestedK -> 1`.
