# Implementation Plan: `z₀` Test Harness

**What this is.** The build guide for the experiment specified in
[`z0_numerical_testing_plan.md`](z0_numerical_testing_plan.md) (the *what/why*)
under the rule argued in [`anchor_selection_procedure.md`](anchor_selection_procedure.md)
(the *theory*). This doc is the *how*: concrete files, function signatures, build
order, and the code-level gotchas to avoid. Still design-level — it names every
piece to write, but the `.wl` itself is written against this checklist, phase by
phase.

---

## 1. File layout

All harness code lives in `SANDBOX/` (next to the infrastructure it reuses);
results and reports land in `AUXT/`.

```
SANDBOX/
  z0_sweep_common.wl     NEW — sweepK, gate scoring, trueSigma aggregation, row schema
  z0_sweep_battery.wl    NEW — the Tier 1–5 case specs + fan providers (driver script)
  z0_sweep_report.wl     NEW — turns the results Dataset into AUXT/z0_sweep_results.md
  bench_lift_variance.wl REUSE — source of runUnlifted/runLifted patterns (copy, don't import)
  sandbox_lift_common.wl REUSE — provides trueSigma[sd, kinRules, pg]
AUXT/
  z0_sweep_results.csv   OUTPUT — one row per (case, k, seed)
  z0_sweep_results.md    OUTPUT — benchmark-style tables + H1/H2/H3 verdict
```

Package-load idiom (copy verbatim from `sandbox_true_sigma.wl:14–18`):

```wl
Module[{pkgRoot},
  pkgRoot = FileNameJoin[{DirectoryName[$InputFileName], ".."}];
  SetDirectory[pkgRoot];
  Get[FileNameJoin[{pkgRoot, "tropical_eval.wl"}]];
  Get[FileNameJoin[{pkgRoot, "SANDBOX", "sandbox_lift_common.wl"}]];   (* trueSigma *)
];
```

---

## 2. Inventory: what already exists (reuse, exact signatures)

| Piece | Signature / call | Returns | Source |
|-------|------------------|---------|--------|
| Unlifted baseline | `runUnlifted[spec, proxyPoly, caseLabel]` | assoc w/ `value, err, sigma, direct, relErr, fan, cfm*` | `bench_lift_variance.wl:38` |
| Lifted run | `runLifted[spec, liftRules, fanDataOpt, caseLabel, liftLabel, direct]` | assoc w/ `value, err, sigma, nDropped, hasConstList, …` | `bench_lift_variance.wl:108` |
| Driver | `EvaluateTropicalMCLifted[spec, {{}}, "LiftRules"->…, "FanData"->…, "NSamples"->N, "RunChecks"->False]` | assoc w/ `Results[[1]]["Re"]`, `["ReErr"]` | `tropical_eval.wl:2268` |
| Lift (exact) | `LiftCoefficients[spec, liftRules]` | `<|"LiftedSpec", "LiftData"|>` (round-trip checked) | `tropical_eval.wl` |
| Sector + gate | `ProcessSectorLifted[liftedSpec, fan[[1]], fan[[2,s]], s, liftData]` | sector assoc w/ `EmptyDomain, HasConstantTerm, FlattenedPolys, Prefactor, Dimension, DomainConstraint, PolynomialExponents` | `tropical_eval.wl:774` |
| **True σ** | `trueSigma[sectorData, kinRules_List, pg_Integer]` | `<|"I1", "I2", "Sigma", "I2Converged"|>` | `sandbox_lift_common.wl:554` |
| Fan build | `ComputeDecomposition[PolytopeVertices[poly^(-1), vars], "ShowProgress"->False]` | `{rays, simplices}` | per `bench:44,144` |
| Var-reduction | `vrFactor[sigU, sigL] = (sigU/sigL)^2` | number / `"N/A"` | `bench:209` |

`σ` convention everywhere: `sigma = stderr * Sqrt[NSamples]` (`bench:73,135`).

---

## 3. New components to build

### 3.1 `sweepK` — the core loop (in `z0_sweep_common.wl`)

```
sweepK[ spec_, primaryRule_, kList_, fanFor_, opts ]  ->  list of row-assocs
```

For each `k ∈ kList`:
1. `rule = primaryRule` with `"k" -> k`.
2. `lifted = LiftCoefficients[spec, {rule}]` — **abort the row if the
   `liftidentity` round-trip fails** (mark `exact -> False`); never score a
   non-exact lift.
3. `z0 = lifted["LiftData"]["z0"]` (exact); compute `z0Num = N[z0]`.
4. `fan = fanFor[k]` (§3.3 — handles `liftdegenerate`).
5. **Gate pass** — `ProcessSectorLifted` over every cone; collect per-cone
   `EmptyDomain`, `HasConstantTerm`, `min Re ãⱼ` (`NewExponents`), `|m_p|`
   (`PivotIndex`/`ZRow`). Derive `gateOK = ∃ surviving cone with HasConstantTerm
   ∧ all Re ãⱼ > 0`.
6. **MC** — for each `seed`: run the driver, record `value, err, sigma_sample,
   relErr`. (Seed handling: §6.)
7. **True σ** — `trueSigma` on each non-empty sector; aggregate (§3.2).
8. Emit one row per seed (raw) + one aggregate row per `k`.

`kList` per case comes from the battery, centered on the predicted
`kStar = Max[1, Ceiling[Abs[Log[10, magnitude]] / Log[10, threshold]]]` and
spanning `1 .. kStar+2` so the plateau (H3) is visible.

### 3.2 True-σ aggregation (decisive metric)

```
aggTrueSigma[sectorTrueSigmas] ->
  if Any[!#["I2Converged"]] : <|"sigmaTrue" -> Infinity, "heavyTail" -> True|>
  else                      : <|"sigmaTrue" -> Sqrt[Total[#["Sigma"]^2]], "heavyTail" -> False|>
```

Var(total) = Σ Var(sector) (sectors independent — Manual §Welford). A single
non-convergent `I2` ⇒ infinite true variance ⇒ the case is heavy-tailed and
**must not be ranked on `sigma_sample`** (this is the Case B/C trap).

### 3.3 Fan provider `fanFor[k]`

```
fanFor[k_] := Module[{auto},
  auto = tryAutoFan[liftedSpec];                 (* PolytopeVertices + ComputeDecomposition *)
  If[auto === $Failed || degenerateQ[liftedSpec],
     Lookup[explicitFans, k, $Failed],            (* precomputed, e.g. Case C {±k,0,∓3} lineality *)
     auto]
];
```

Detect degeneracy *before* trusting the auto fan (the lifted Newton polytope is
lower-dimensional when monomials are few). When `explicitFans` lacks an entry for
that `k`, **skip the cell with a logged reason** — do not silently fall back to
the unlifted pipeline (which `EvaluateTropicalMCLifted` does on empty rules).

### 3.4 Row schema (the unit of output)

```
<| "caseId", "n", "magnitude", "alpha", "k", "kStar", "z0", "z0Num",
   "seed", "value", "err", "sigmaSample", "relErr",
   "sigmaTrue", "heavyTail", "gateOK",
   "hasConstAny", "minMargin", "absMp", "nDropped", "nSectors",
   "varRedTrue"(=(sigmaTrue_unlifted/sigmaTrue)^2), "wallClock", "exact", "skipReason" |>
```

---

## 4. Build order (milestones, each with an acceptance check)

| Phase | Concrete tasks | Acceptance check |
|-------|----------------|------------------|
| **P0** | Scaffold `z0_sweep_common.wl`: package load, copy `runUnlifted`/`runLifted` bodies, write `sweepK`, `aggTrueSigma`, `fanFor`, row schema, CSV writer. Resolve **seed control** (§6). | `sweepK[specA, {2,0}-rule, {2}, Automatic, …]` **reproduces Case A's 17.95×** from `benchmark_results.md`. Gate metrics match (5/6 cones `HasConstantTerm=True`). |
| **P1** | `z0_sweep_battery.wl` Tier 1 (T1a–T1e); run sweep `k=1..kStar+2`; write CSV. | `Δk = argmin_true(k) − kStar` computed per case; table emitted. |
| **P2** | Tier 2 (small-coeff) + Tier 3 (Case C degenerate, explicit fans `{±k,0,∓3}`). | Tier 3: every `k` reports `heavyTail=True` ⇒ "no-lift"; small-coeff `kStar` mirrors large-coeff. |
| **P3** | Tier 4 (multi-rule shared anchor — residuals `cᵢ = Cᵢ/z0^{kᵢ}`) + Tier 5 (n=3,4 cost curve). | Single-rule `k=2` beats two-rule `k=1` on `sigmaTrue`∧`relErr` (not `sigmaSample`); cost vs `k` tabulated. |
| **P4** | `z0_sweep_report.wl`: render `AUXT/z0_sweep_results.md` (per-case tables + `Δk` summary + H1/H2/H3 verdict + falsification log). | Report regenerates from CSV with no manual edits; go/no-go stated. |

---

## 5. Truth oracle implementation (accuracy column)

Per `z0_numerical_testing_plan.md §4`, and accounting for the **known
NIntegrate-at-extreme-coefficients hazard**:

1. **Closed form first.** Hard-code exact truths for cases that have them
   (Toy-0 family). Store as `caseTruth[caseId]`; when present, `relErr` uses it.
2. **Hardened NIntegrate** otherwise: raise `WorkingPrecision`, `MaxRecursion`;
   pass the balance locus `xᵅ = 1/|C|` as a singular-point/`Exclusions` hint;
   run two methods and accept only on `<1%` agreement.
3. On disagreement → set `relErr -> Indeterminate`, `oracleTrusted -> False`, and
   rank that case by `sigmaTrue` alone (needs no global truth).

`ValidateLiftedDecomposition[…]` is used only as a *bookkeeping* cross-check
(dropped-sector accounting), never as the truth source.

---

## 6. Code-level gotchas (grounded in the source)

- **Option name is `"NSamples"`, not `NumSamples`** (`bench:60,122`). Easy to get
  wrong; the driver silently ignores unknown options.
- **Seed control is unverified.** `bench_lift_variance.wl` passes *no* seed to
  `EvaluateTropicalMC`. P0 must determine whether the driver/C++ MC exposes a seed
  option; if not, either (a) wrap each call in `SeedRandom[seed]` if that reaches
  the sampler, or (b) extend the driver. Seed-robustness (≥5 seeds) is meaningless
  until this is real — block P1 on it.
- **`trueSigma` needs `sandbox_lift_common.wl` loaded** *and* a sector assoc with
  `FlattenedPolys`/`Prefactor`/`Dimension`/`DomainConstraint` — i.e. the output of
  `ProcessSectorLifted`, not of `EvaluateTropicalMCLifted`. Re-process sectors
  exactly as `runLifted` does (`bench:138–168`).
- **`relErr` denominator** uses the *oracle*, not the lifted value
  (`bench:88` does `(value-direct)/direct`). Keep that, but swap `direct` for the
  closed form when available.
- **Degenerate auto-fan returns garbage, not an error**, in some cases — gate on
  `degenerateQ` explicitly; don't rely on `liftdegenerate` always firing.
- **Rank on `sigmaTrue`, gate on `gateOK`.** `sigmaSample` is reporting-only.
  Reaffirm in code: the "best `k`" selector reads `sigmaTrue` + `gateOK`, never
  `sigmaSample`.
- **`NSamples` fixed across all cells** (e.g. `10^6`) so `varRed` is comparable;
  wall-clock is logged but not a ranking input (it conflates compile + NIntegrate).

---

## 7. Definition of done

- [ ] P0 harness reproduces Case A 17.95× and gate metrics from `benchmark_results.md`.
- [ ] Seed control is real and verified (≥5 distinct seeds give distinct streams).
- [ ] CSV emitted for all Tier 1–5 cells; every row has `exact=True` or a `skipReason`.
- [ ] `AUXT/z0_sweep_results.md` regenerates from CSV and states, per tier, whether
      `argmin_true(k) == kStar` (H1), Case C reports no-lift (H2), and `sigmaTrue`
      plateaus for `k ≥ kStar` (H3).
- [ ] Falsification log lists every `Δk ≠ 0` with a diagnosed cause.
- [ ] Go/no-go recorded: if H1 holds, the follow-up is wiring `kStar` into
      `DetectExtremeCoefficients` (replacing `SuggestedK -> 1`, `tropical_eval.wl:603`);
      if not, this CSV is the dataset the deferred fitting step trains against.
```
