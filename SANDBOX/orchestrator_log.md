# Orchestrator Log — TROPICAL_MONTE_CARLOv2 Migration

Maintained by the Lead Systems Architect / Orchestrator. One entry per gate or significant event.
Plan: `PLAN_SUM/plan.md`. Phases: 0 → 1 → 2A → 2B → 3, strictly gated.

---

## Phase 0 — Setup and Baseline

**Date:** 2026-06-12

### Step 0.1 — Copy codebase ✅
- `ditto TROPICAL_MONTE_CARLO → TROPICAL_MONTE_CARLOv2` completed.
- `INTERFILES/` verified empty (was already empty; the `rm -f` glob had nothing to remove).
- `SANDBOX/` created.
- Original `TROPICAL_MONTE_CARLO/` untouched (and will remain so).
- Snapshot of the protected `ProcessSector` divergence-check + `IsDivergent` early-return branch
  (pre-removal lines 267–348) saved to `SANDBOX/protected_processsector_snapshot.txt` for
  post-Phase-2A diff verification (plan §1.3 rule 4 / §6 step A6).

### Step 0.2 — Smoke test ✅
- `wolframscript -file EXAMPLES/test_quick.wl` in the v2 copy: **PASS**.
- 3 rays / 3 sectors, all convergent; direct NIntegrate 0.3926989769 vs sector sum 0.3926986968,
  rel. err 7.13e-7. Polymake plumbing works in the new location.

### Step 0.3 — Baseline regression suite (in progress)
- **Incident:** first run of the plan's verbatim command
  (`wolframscript -code '...; RunAllTests[]'`) executed all 14 worked examples successfully but
  `RunAllTests[]` itself did **not** run: parsing the `-code` string pre-created `Global`RunAllTests`,
  which shadowed `TropicalEval`RunAllTests` (the `shdw` warning confirmed it); the final expression
  echoed unevaluated as `Global`RunAllTests[]`.
- **Fix:** invoke via runtime resolution —
  `wolframscript -code 'SetDirectory["EXAMPLES"]; Get["tropical_eval_examples.wl"]; ToExpression["RunAllTests[]"]'`.
  (`ToExpression` parses after `TropicalEval`` is first on `$ContextPath`, so the symbol resolves to
  the package definition. No package file was modified.)
- Re-run completed; full transcript (744 lines) saved to `SANDBOX/baseline_tests.txt`.
- **Result: 18 PASSED, 0 FAILED** (Tests 1–18; suite wires in RunTest1–18 only — RunTest19–22
  exist in the file but are not called by `RunAllTests`, consistent with the plan's removal map).
- All 14 worked examples passed, including both C++ codegen examples (6, 7) — g++ retried
  without `-fopenmp` as expected; NIntegrate cross-checks at 1e-4..7e-7, MC at ~1e-3.
- Retained-subset baseline (Tests 1, 2, 5, 6, 7 blocks, 74 lines) extracted to
  `SANDBOX/baseline_retained_subset.txt` — this is the regression target for Gates 2A and 3.
  Key reference values:
  - Test 1: A=2 sector sum 0.3610108194006753; A=3 sector sum 0.20308396820714944 — PASS
  - Test 2: sector sum 0.3112722115271244 − 0.12261702711022615 I — PASS
  - Test 5: 100-pt kinematic scan, max rel err 4.52e-5 — PASS
  - Test 6: Cases A/B/C sector sums 0.0007849539982702149 / 144.04704179722273 /
    0.00168677137414064 — PASS (the motivating large-coefficient cases)
  - Test 7: 3D 0.2648333702850404, 4D 0.06094511956155343 — PASS

### Gate 0 — ✅ PASSED (2026-06-12)
- [x] test_quick passes
- [x] Baseline transcript saved with all tests passing (18/18)

---

## Phase 1 — Sandbox

**Date:** 2026-06-12

### Worker dispatches
1. **Dispatch 1 (Sonnet):** created the four §5 sandbox files (`sandbox_lift_common.wl`,
   `sandbox_toy0_1d.wl`, `sandbox_toy1_eps.wl`, `sandbox_toy2_largecoeff.wl`). All exactness
   checks passed on first delivery. Worker reported variance INCREASE from lifting — flagged for
   orchestrator review rather than accepted.
2. **Orchestrator review found two defects in the measurement, not the math:**
   - `sectorMC` sigma aggregation (`StdErr/feasFrac*Sqrt[N]`) overstated lifted sigma by ~1/feasFrac
     in domain-constrained sectors (toy2 sector 3: 435x inflation).
   - toy1's `buildExtendedFan` appended z-ray (0,0,-1), which is not in the lineality space
     ν = (0,k,-1) of the degenerate lifted polytope's normal fan → cones crossed normal-fan walls
     → dominance broken → artificial constant-term loss. (Plan gap: §5.3 step 4 assumed
     Polymake handles the lifted polytope; it is degenerate for toy1 — 3 support points in R³.)
3. **Dispatch 2 (Sonnet):** fixed `sectorMC` statistics (true estimator sigma with infeasible
   draws contributing 0), replaced toy1's fan with the principled lineality-completion refinement
   (rays {(1,0,0),(0,1,0),(-1,-1,0),(0,k,-1),(0,-k,1)}, 6 cones), tightened the pivot
   admissibility test to require real atilde. All exactness checks still pass; k=1 vs k=2 now
   genuinely differ.
4. **Dispatch 3 (Sonnet):** built `sandbox_true_sigma.wl` — exact per-sample sigma
   (σ² = ∫g² − (∫g)², g = flattened integrand × indicator) per sector via NIntegrate, plus
   per-pivot diagnostic tables. Sample sigma had proven unreliable (unlifted toy1 MC misses the
   spike mass entirely: sample mean ~927 vs true 5e5).

### Verified results (orchestrator's own runs: `verify_toy*.txt`, `verify_true_sigma.txt`)
- Exactness: all five PASS lines (toy0 1e-13/1e-16; toy1 7.9e-5 both k; toy2 1.2e-5).
  §3.6 hand-check reproduces; EmptyDomain classification works.
- True sigma: TOY1 unlifted INF, lifted INF (both k); TOY2 unlifted 0.0166, lifted INF
  (one sector, cone 8, carrying 0.2% of the integral; excluding it lifted = 0.0039, 4.2x better).
- Root cause of every divergent σ²: §3.3 constant-term loss (HasConstantTerm → False).
  Structural (no admissible constant-preserving pivot) in all toy1 sectors and toy2 cone 8;
  plan's §3.5 ranking provably suboptimal in toy2 cone 7 (constant-preserving pivot p=2 exists
  but |m_p|=1 ranks above it).
- Full numbers: `SANDBOX/phase1_results.txt`.

### Gate 1 — ❌ FAILED on the variance criterion (2026-06-12) → HARD STOP per plan §5.5
- [x] All three toys PASS exactness (rel. err < 1e-3)
- [x] Variance comparison recorded in `SANDBOX/phase1_results.txt`
- [ ] Lifted sigma materially below unlifted on Toys 1 and 2 — **NO**
- **Decision:** per plan §5.5, Phase 2A is NOT started. Numbers and mitigation options
  reported to the project owner. Phases 2A/2B/3 remain blocked pending owner decision.

### Owner decision (2026-06-12, after Gate 1 report)
Owner reviewed the report and directed continuation: v2 restricts to convergent (non-divergent)
integrals, for which the pipeline's correctness is unaffected; the code must run exactly the
same for convergent examples. **Gate 1 STOP lifted by owner — proceeding to Phase 2A.**
Orchestrator amendments carried into Phase 2B, justified by Gate 1 evidence:
1. §3.5 pivot ranking reordered: constant-term preservation ranks ABOVE |m_p| = 1
   (toy2 cone 7 proved the plan's order discards strictly better pivots).
2. The lifted-fan construction must handle degenerate lifted Newton polytopes
   (lineality-completion refinement, validated in the sandbox).
The HasConstantTerm=False / infinite-σ² risk remains open and will be re-measured in the
Phase 3 benchmarks (§8.3) as the plan specifies.

---

## Phase 2A — Remove divergent machinery

**Date:** 2026-06-12. Worker dispatch 4 (Sonnet) executed A1–A10 per plan §6, bottom-up with
reload checks. `tropical_eval.wl` 3834 → 1389 lines; `tropical_eval_examples.wl` 3129 → 1494;
`test_divergent.wl` deleted; SUMMARY.txt updated; new messages `divergentinput` / `noregulator`
/ `nodivergent` defined. Regulator guards use the Missing-tolerant pattern
`!MatchQ[spec["RegulatorSymbol"], None | _Missing]` (examples omit the key).

### Gate 2A — ✅ PASSED (orchestrator-verified, not worker-claimed)
1. [x] `Get` loads clean, LOAD OK, no warnings.
2. [x] Retained regression: Tests 1, 2, 5, 6, 7 output **identical** to
   `SANDBOX/baseline_retained_subset.txt` (diff clean modulo one extraction blank line).
3. [x] `RunTest3v2` passes: Part A ε-regulated spec → `$Failed` + `noregulator`;
   Part B genuinely divergent spec → `$Failed` + `divergentinput`. Suite: 6 PASSED, 0 FAILED.
4. [x] `test_cpp.wl` passes (MC rel err 0.0015%); `test_quick.wl` passes (7.13e-7).
5. [x] Protected ProcessSector region diff vs snapshot: ONLY the planned line-273
   simplification (`a0vals = effectiveAVals;`). IsDivergent branch fully intact.
6. [x] Greps: zero references to all 15 removed symbols/patterns.

---

## Phase 2B — Add lifting capability

Amendments in force (owner-approved 2026-06-12): pivot ranking prefers constant-term
preservation ABOVE |m_p| = 1; automatic lifted-fan construction must fail cleanly on
degenerate lifted polytopes (explicit "FanData" remains the workaround, per sandbox).

**B1 (FlattenSector refactor)** — dispatch 5. Gate verified by orchestrator: suite 6 PASSED,
retained tests identical to baseline; both FlattenSector branches probed directly (divergent
branch: IsDivergent True, last-wins divVar, raw |detM| prefactor — WL `Return` propagates
correctly out of the nested Module to the function boundary, probe-confirmed).

**B2+B3 (Module 1b + ProcessSectorLifted)** — dispatch 6. Round-trips pass incl. two-rule
residuals {1, 300} and negative-coefficient residual {-1}. Toys rewired (useSandboxImpl flag);
orchestrator re-ran all three toys under the package implementation: all PASS; toy2 cone 7 now
picks the constant-preserving pivot p=2 (amended ranking), sector integral unchanged (4.3e-8);
only structural cone 8 remains HasConstantTerm=False.

**B4+B5+B6** — dispatch 7. ValidateLiftedDecomposition: toy1 rel err 1.4e-4 (3 dropped),
toy2 4.4e-5 (2 dropped). B5 keyed indicator insertion: unlifted codegen byte-identical
(orchestrator-verified diff at this gate); toy2 lifted C++ binary (orchestrator re-ran):
7.8158e-4 ± 4.0e-6 vs 7.8501e-4 → 0.44% < 1%. 4 indicator blocks emitted (cones 2/3/7/8).
B6 FeasibleFraction extension verified on all three cases.

**B7 (driver wiring + EvaluateTropicalMCLifted)** — dispatch 8, plus a fix dispatch 9.
- Dispatch 8 delivered the wiring; unlifted driver end-to-end 0.069% vs Pi/8; lifted k=2
  explicit-rule end-to-end 0.44%; kinematic lifted scan (lam = 0.5/1/2) all < 0.6%;
  toy1 degeneracy guard fires correctly; explicit-fan toy1 runs honestly (7.7% deviation,
  known heavy-tail config).
- **Incident:** dispatch 8 reported Automatic mode failing with `liftdegenerate` on the
  flagship toy2 case, misdiagnosed as a degenerate polytope. Orchestrator disproved
  (k=1 lifted polytope is full-dimensional, det = -2; standalone fan: 6 rays / 8 simplices)
  and reproduced the real error (polymake "dimension mismatch"). Root cause (dispatch 9):
  DetectExtremeCoefficients emits key "SuggestedK" but LiftCoefficients reads "k" →
  Missing exponent → ragged vertex rows. Fixed by key remapping in the wrapper + Quiet
  translation of genuine polymake failures into liftdegenerate. Automatic k=1 now runs
  end-to-end (all 6 live sectors HasConstantTerm=True; 7.622e-4 ± 1.0e-5 at 1e5 samples,
  2.2σ from reference — statistically consistent).
- **Incident:** the /tmp pre-B5 codegen reference file was overwritten by a worker with
  lifted output mid-B7 (driver writes tropical_mc_generated.cpp, so no filename collision —
  a worker copy did it). No package defect: fresh test_cpp output is the correct 3-integrand
  unlifted code (0.0015% err) and true byte-identity was orchestrator-verified at the B4–B6
  gate. Durable reference now at `SANDBOX/reference_test_cpp_unlifted.cpp`.

---

## Phase 3 — Regression, new tests, benchmarks

**Date:** 2026-06-12. Dispatch 10 (Sonnet) delivered §8.2/§8.3/§8.4.

### §8.1 Regression gate — ✅ PASSED (orchestrator-run)
Full suite after all Phase 2B changes: 6 PASSED, 0 FAILED; Tests 1/2/5/6/7 **identical** to
`SANDBOX/baseline_retained_subset.txt`.

### §8.2 Test 23 (`EXAMPLES/test_lifted.wl`) — ✅ ALL PASS (orchestrator re-run)
- 23A exactness (Toy 1, explicit fan): rel err 1.45e-4 < 0.5%; 3 dropped sectors. PASS.
- 23B end-to-end C++ (Test 6 Case A, explicit k=2, 1e6 samples): 7.8158e-4 ± 3.97e-6 vs
  7.8501e-4 → 0.44% < 1%, error-bar-consistent. PASS. (Automatic k=1 informational:
  7.7635e-4 ± 3.58e-6, 2.4σ low — k=2 remains the recommended configuration.)
- 23C error paths: liftcomplex and liftnopivot both fire with $Failed, no crash. PASS.
- 23D EmptyDomain (Toy 0 k=1): 4 dropped, rel err 4.8e-11 < 0.1%. PASS.

### §8.3 Benchmarks (`SANDBOX/benchmark_results.md`) — mixed verdict, honestly reported
| Case | σ_unlifted | σ_lifted | Variance reduction | ≥10x? |
|---|---|---|---|---|
| A | 0.01680 | 0.003966 | **17.95x** | PASS |
| B | 60005 | 1212 (2-rule k=1) | 2450x (unreliable) / 47x (1-rule k=2, reliable) | report-only |
| C | 0.03967 | 0.08236 | 0.23x | **FAIL** |
- Case A's lifted/unlifted sigmas match Phase 1's exact true-sigma analysis (0.0039/0.0166) —
  strong cross-validation of the benchmark.
- Case C failure is STRUCTURAL: 3-term polynomial → lifted support is 3 points in R³ →
  degenerate polytope for every k; all sectors lose the constant term; sample sigma is
  optimistic and the MC estimate itself deviates ~40%. Same geometry class as Toy 1
  (Phase 1 infinite-variance finding). Documented in benchmark_results.md + SUMMARY.txt.
- Case B: the honest configuration (1-rule k=2) achieves 47x variance reduction at 6.5%
  deviation; the 2450x 2-rule figure is flagged unreliable (HasConstantTerm=False sector).

### §8.4 Documentation — SUMMARY.txt section 8 added (new functions, messages, LiftData
option, amended pivot ranking, degenerate-polytope limitation + FanData escape hatch).

---

## MIGRATION COMPLETE — final status (2026-06-12)
- Phases 0 → 3 executed in order; every gate orchestrator-verified.
- Convergent regression preserved end-to-end (baseline-identical at Gates 2A, B1, and §8.1).
- Lifting capability shipped and verified correct (exactness everywhere); variance benefit
  proven on Case A (17.95x) and Case B (47x reliable config); structurally impossible for
  few-monomial integrands whose lifted polytope is degenerate (Case C, Toy 1) — documented
  limitation, fails loudly via liftdegenerate/explicit-fan path rather than silently.
- Known future work: log-space simplex remap for HasConstantTerm=False sectors (plan §3.4),
  automatic lineality-completion fans for degenerate lifted polytopes (sandbox-validated math
  in SANDBOX/sandbox_toy1_eps.wl), per-kinematic z0 lifting (plan §9 case 6).
