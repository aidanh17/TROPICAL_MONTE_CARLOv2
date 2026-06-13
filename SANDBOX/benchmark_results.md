# TROPICAL_MONTE_CARLO v2: Lift-Variance Benchmark

**Date:** 2026-06-12  
**Samples per run:** 1,000,000  
**Caveat:** Wall-clock includes C++ compilation and optional NIntegrate; MC-loop-only time not separately isolated by the driver.  Timings are total AbsoluteTiming of the driver call.

---

## Case A: `P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2`, B={-2}

Lift rule: explicit `{PolyIndex=1, ExponentVector={2,0}, k=2}`.  Automatic fan (Polymake, full-dimensional 3D Newton polytope).

| Mode | MC Value | StdErr | Per-sample sigma | NIntegrate | Rel.Dev. | Sectors | Dropped | Wall(s) |
|------|----------|--------|-----------------|------------|----------|---------|---------|---------|
| Unlifted | 7.897e-4 | 1.680e-5 | 0.01680 | 7.850e-4 | 0.00597 | 4 | 0 | 4.10 |
| Lifted k=2 | 7.816e-4 | 3.966e-6 | 0.003966 | 7.850e-4 | 0.00437 | 6 | 2 | 5.59 |

**Variance-reduction factor = 17.95**  -->  **PASS** (>= 10x)

CFM (20 samples/sector): Unlifted min/max = 1.51e-13 / 3.89e-5; Lifted k=2 min/max = 0 / 5.86e-3.

**HasConstantTerm=False** in at least one lifted sector (cone 8; see Warning above).  The sample sigma for the lifted run may be OPTIMISTIC for that sector.

---

## Case B: `P = 10^-4 + 10^4 x1^2 + 10^-4 x2^2 + 10^4 x1 x2^2 + x1^2 x2`, B={-2}

Two lift strategies tried: **(i)** two-rule k=1 (primary {2,0} + secondary {1,2} with residual); **(ii)** single-rule k=2 on primary {2,0}.  Case B is report-only (no >= 10x gate).

| Mode | MC Value | Per-sample sigma | NIntegrate | Rel.Dev. | Sectors | Dropped | Var.Red. | Wall(s) |
|------|----------|-----------------|------------|----------|---------|---------|----------|---------|
| Unlifted | 175.7 | 60005 | 144.06 | 0.2195 | 5 | 0 | — | 1.80 |
| Lifted 2-rule k=1 | 52.33 | 1212 | 144.06 | 0.6367 | 9 | 3 | 2450 | 5.70 |
| Lifted 1-rule k=2 | 153.4 | 8723 | 144.06 | 0.0652 | 6 | 2 | 47.3 | 50.1 |

Notes: The 2-rule k=1 run achieves 2450x sample-variance reduction but the MC estimate (52.3) deviates 64% from NIntegrate (144.1) — one sector has HasConstantTerm=False, making the sample sigma optimistic.  The 1-rule k=2 run gives 47x reduction with a more reliable 6.5% deviation.  HasConstantTerm=True in all 1-rule-k=2 sectors.

---

## Case C: `P = 1 + 10^8 x1^3*x2 + x2^3`, B={-3}

The lifted Newton polytope (3 monomials in 3D) is lower-dimensional for both k=2 and k=4.  Automatic fan detection fires `liftdegenerate`.  **Explicit fans supplied** (5 rays: {1,0,0},{0,1,0},{-1,-3,0} plus the lineality pair {2,0,-3}/{-2,0,3} for k=2 and {4,0,-3}/{-4,0,3} for k=4; 6 simplices shared).

| Mode | MC Value | Per-sample sigma | NIntegrate | Rel.Dev. | Sectors | Dropped | Var.Red. | Wall(s) |
|------|----------|-----------------|------------|----------|---------|---------|----------|---------|
| Unlifted | 1.687e-3 | 0.03967 | 1.684e-3 | 0.00153 | 3 | 0 | — | 2.49 |
| Lifted k=2 | 1.015e-3 | 0.08236 | 1.684e-3 | 0.3975 | 3 | 3 | 0.232 | 1.72 |
| Lifted k=4 | 1.015e-3 | 0.08236 | 1.684e-3 | 0.3975 | 3 | 3 | 0.232 | 2.14 |

**Best variance-reduction = 0.232**  -->  **FAIL** (< 10x)

**Structural root cause:** All 6 lifted sectors (both k=2 and k=4) have HasConstantTerm=False.  After delta-elimination via the z-coordinate pivot, the re-cleared polynomial in (y1,y2) lacks a constant term.  Consequently:
1. The sample-based sigma (0.0824) underestimates the true sigma (the MC cannot see the heavy tail near the boundary of the domain).
2. The MC estimate (1.015e-3) misses ~40% of the integral mass relative to NIntegrate (1.684e-3).
3. The lifted sigma is LARGER than the unlifted sigma (0.232 < 1 = variance goes up, not down).

This is a **known limitation** of the current lifting implementation when the lifted Newton polytope is lower-dimensional and all pivot choices produce HasConstantTerm=False.  It is documented as a degenerate-lifted-polytope limitation in SUMMARY.txt.

---

## Summary Table

| Case | sigma_unlifted | sigma_lifted (best) | Var.reduction | >= 10x? |
|------|---------------|---------------------|---------------|---------|
| A    | 0.01680 | 0.003966 | 17.95 | **PASS** |
| B    | 60005 | 1212 | 2450 | report-only |
| C    | 0.03967 | 0.08236 | 0.232 | **FAIL** |

---

## Important Caveats on Variance Estimation

**Sample-based sigma can severely underestimate true sigma for heavy-tailed integrands.**
See `SANDBOX/phase1_results.txt` for a full true-sigma analysis: the unlifted Toy-1 MC
missed the spike mass entirely (sample means ~927 vs true sum 5e5); true-sigma analysis
found infinite-variance configurations in BOTH lifted and unlifted Toy-1 sectors.
For Toy 2, one lifted sector had divergent I2 (0.2% contribution); excluding it, the
remaining 5 sectors showed 4.2x true-sigma improvement.

**When HasConstantTerm = False** for a lifted sector (re-clearing loses the constant term),
the sector polynomial Qtilde vanishes at a cube corner, making Qtilde^B non-integrable
in L^2 for B << 0.  The sample-based sigma (from finite MC) is then an optimistic
(possibly dramatically underestimated) proxy for the true variance.  Per-case status:

- Case A lifted sectors: 1 of 6 sectors has HasConstantTerm=False (cone 8).  Other 5 are True.  Sample sigma likely reliable for the dominant contribution.
- Case B lifted sectors: All HasConstantTerm=True for the 1-rule-k=2 run.  One sector False in the 2-rule-k=1 run.
- Case C lifted sectors: ALL 6 sectors (both k=2 and k=4) have HasConstantTerm=False.  Sample sigma is unreliable; variance comparison invalid.  FAIL is structural, not a sampling artifact.
