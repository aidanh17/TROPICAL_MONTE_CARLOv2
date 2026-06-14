# Lifting Procedure Fix Plan

## Diagnosis

Two independent failures block lifting on CALC2-class integrands (products of
degree-1 atom factors `(u_j + v_j·x^{α_j})^{γ_j}` with complex `γ_j`).

### Failure 1 — `TropicalEval::liftidentity` (warning, non-fatal)

**Location:** `LiftCoefficients`, lines 810–823 of `tropical_eval.wl`

**Root cause:** The fallback identity check is

```mathematica
If[!TrueQ[Simplify[liftedSubbed - original] === 0], ...]
```

`Simplify` cannot reduce floating-point residuals to the symbol `0`.  For a
degree-1 atom `u + v·x^α` with floating-point `u` and `v`:

- `z0 = |v|^{1/k}` is computed as a `Real` (not an exact root)
- `c = Simplify[v / z0^k]` carries machine-epsilon rounding
- `c · z0^k ≠ v` in exact symbolic arithmetic even though the identity holds to
  machine precision

So the message fires as a false positive.  The code does **not abort** (comment
says "Do not abort"), but the spurious message misleads debugging.

### Failure 2 — `TropicalEval::liftcomplex` (fatal, returns `$Failed`)

**Location:** `ProcessSectorLifted`, lines 956–993 of `tropical_eval.wl`

**Root cause:** `atildeRaw` is computed as

```mathematica
atildeRaw[[j]] = [real piece] + Sum_k  B_k · rcMin_{k,j}
```

When `B_k = γ_k = γ_re + i·γ_im` (complex), `atildeRaw` has imaginary parts
`~γ_im · rcMin_{k,j}`.  The pivot acceptance test at line 962 requires

```mathematica
Abs[Im[av]] < 10^-12    (* absolute threshold *)
```

For CALC2 `γ_im ≈ −14` and `rcMin > 0`, the imaginary part is O(14) — the
threshold is never met.  No admissible pivot is found → `liftcomplex` fires →
`$Failed` is returned for every cone → the whole run returns `$Failed`.

The failure is structurally unavoidable with the current algorithm because the
domain indicator `{x : atilde_j > 0 for all j}` requires **real** `atilde`.

**Wasted work path:** `EvaluateTropicalMCLifted` builds the full (n+1)-dimensional
fan before `EvaluateTropicalMC` calls `ProcessSectorLifted`, so all polymake fan
computation happens for nothing.

---

## Fix 1 — Numerical tolerance in `liftidentity` check

**File:** `tropical_eval.wl`, lines 809–823

**Change:** Replace the `Simplify[...] === 0` fallback with a
coefficient-wise relative magnitude check.

```mathematica
(* BEFORE (line 813-815): *)
If[!TrueQ[liftedSubbed === original],
  If[!TrueQ[Simplify[liftedSubbed - original] === 0],
    Message[TropicalEval::liftidentity, j]; ...

(* AFTER: *)
If[!TrueQ[liftedSubbed === original],
  Module[{parsedDiff, scaleOrig, maxResid},
    parsedDiff = ParsePolynomial[Expand[liftedSubbed - original], vars];
    scaleOrig  = Max[1., Max[Abs[N[#[[1]]]] & /@ ParsePolynomial[Expand[original], vars]]];
    maxResid   = Max[Abs[N[#[[1]]]] & /@ parsedDiff];
    If[maxResid > 1*^-8 * scaleOrig,
      Message[TropicalEval::liftidentity, j]; ...
    ]
  ]
```

**Effect:** False positives for floating-point atom coefficients are
eliminated.  A genuine mismatch (> 1e-8 relative) still triggers the message.

**Tests to add (in `test_complex_lifted.wl`):**
- Case 24E: polynomial `3.7 + 2.3e-6 x[1]^2 + x[2]^2` (non-unit constant
  term, real exponent).  Assert `liftidentity` does **not** fire and
  `ValidateLiftedDecomposition` passes.

---

## Fix 2 — Pre-flight complex-exponent detection in `EvaluateTropicalMCLifted`

**File:** `tropical_eval.wl`, function `EvaluateTropicalMCLifted` (~line 2595)

**Change:** Before `LiftCoefficients` is called, inspect
`integrandSpec["PolynomialExponents"]` for complex values and return `$Failed`
immediately with a new, actionable message.

```mathematica
(* Insert after line 2630 (rules are finalized), before line 2633: *)
Module[{bVals},
  bVals = N[integrandSpec["PolynomialExponents"]];
  If[AnyTrue[bVals, (NumericQ[#] && Abs[Im[#]] >= 1*^-12) &],
    Message[TropicalEval::liftcomplexexponents,
            Count[bVals, _?(Abs[Im[N[#]]] >= 1*^-12 &)],
            Length[bVals]];
    Return[$Failed]
  ]
]
```

Add the new message template near line 183:

```mathematica
TropicalEval::liftcomplexexponents =
  "EvaluateTropicalMCLifted: `1` of `2` polynomial exponents (B_k / γ_k) are \
complex.  The current lifting algorithm requires real exponents to construct a \
real-valued domain indicator.  Options: (a) use \"ComplexExponentMode\"->\"SplitRealImag\" \
(Fix 3 in plan.md) to fold imaginary parts into an oscillatory phase weight; \
(b) condition the exponent representation before lifting.";
```

**Effect:** Saves the entire polymake fan computation (~seconds) that currently
runs before the inevitable `$Failed`.  Replaces the confusing `liftcomplex`
message (which fires per-cone, far downstream) with a single upfront diagnostic.

---

## Fix 3 — Complex-exponent splitting mode (core algorithmic extension)

This is the fix that makes lifting **actually work** on CALC2-class integrands.

### Mathematical basis

Decompose each polynomial exponent:

```
B_k = Re(B_k) + i·Im(B_k)

∏_k P_k^{B_k} = ∏_k P_k^{Re(B_k)}  ·  exp(i · Sum_k Im(B_k) · log P_k)
                 ^^^^^^^^^^^^^^^^^^^^    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                  "real weight"              "oscillatory phase"
```

The tropical fan and sector decomposition are computed for the **real weight**
only (exponents `Re(B_k)` are real → standard algorithm works).  The oscillatory
phase is evaluated per MC sample and multiplied into the sample value.

### 3a — New option `"ComplexExponentMode"`

Add to `Options[EvaluateTropicalMCLifted]` and `Options[EvaluateTropicalMC]`:

```mathematica
"ComplexExponentMode" -> "Reject"   (* existing behavior *)
(* new valid value: "SplitRealImag" *)
```

### 3b — Spec pre-processing in `EvaluateTropicalMCLifted`

When `"ComplexExponentMode" -> "SplitRealImag"` and complex exponents are
detected:

1. Build `realSpec`: copy `integrandSpec`, replace
   `"PolynomialExponents"` with `Re /@ integrandSpec["PolynomialExponents"]`.
2. Compute `imagExps = Im /@ integrandSpec["PolynomialExponents"]`.
3. Store `imagExps` in the `LiftData` association under key `"ImagExponents"`.
4. Pass `realSpec` (not `integrandSpec`) to `LiftCoefficients` and the fan
   builder.

```mathematica
(* In EvaluateTropicalMCLifted, after complex-exponent detection: *)
If[complexExpMode === "SplitRealImag" && hasComplexExps,
  imagExps  = Im[N[#]] & /@ integrandSpec["PolynomialExponents"];
  realSpec  = ReplacePart[integrandSpec,
    "PolynomialExponents" -> Re[N[#]] & /@ integrandSpec["PolynomialExponents"]];
  ,
  imagExps = ConstantArray[0., Length[integrandSpec["PolynomialExponents"]]];
  realSpec = integrandSpec
];
(* then pass realSpec to LiftCoefficients and fan builder *)
```

### 3c — Propagate `ImagExponents` through `LiftData`

In `LiftCoefficients`, forward the `"ImagExponents"` key from the input spec
(if present) to the returned `liftData` association.  This requires adding
`imagExps = Lookup[integrandSpec, "ImagExponents", None]` and appending it.

### 3d — C++ code generator changes

**File:** the `GenerateCPPCode` function (or equivalent) in `tropical_eval.wl`

When `liftData["ImagExponents"] =!= None` and any `imagExps` entry is non-zero:

1. The generated integrand function returns `std::complex<double>` instead of
   `double`.
2. Append to the sample value computation:

```cpp
// oscillatory phase: exp(i * Sum_k imB_k * log(P_k(x)))
double logPhaseArg = 0.0;
for (int k = 0; k < nPolys; k++) {
    logPhaseArg += imB[k] * std::log(std::abs(poly_k(x)));  // poly_k already evaluated
}
std::complex<double> phase = std::exp(std::complex<double>(0.0, logPhaseArg));
result *= phase;
```

3. The `imB[]` array is emitted as a `constexpr double imB[] = {...};` literal
   from the Mathematica values of `imagExps`.
4. The result accumulator in the MC driver changes from `double` to
   `std::complex<double>`, and the output JSON emits both `Re` and `Im`.

**Note:** The MC driver already handles complex results (used by test 24A–C).
The main addition is the per-sample phase computation.

### 3e — Variance caveat

The imaginary exponent `γ_im ≈ −14` produces rapid oscillations in
`exp(i · γ_im · log P(x))` when `P(x)` spans many decades.  This may not
reduce variance compared to the unlifted case — it lifts the coefficient
extremes while leaving the oscillatory cancellation.  Document this limitation
in the new message text (Fix 2) and in a code comment near the option.

### Tests to add

- **Case 24F:** A two-variable integral with a single degree-1 atom
  `(u + v·x)^{γ_re + i·γ_im}` where `v ~ 1e-6` (extreme) and `γ_im = −2`
  (moderate, tractable).  Expected behavior with
  `"ComplexExponentMode"->"SplitRealImag"`:
  - Pre-flight check no longer fires `liftcomplexexponents`
  - `ProcessSectorLifted` succeeds (real exponents only)
  - MC result within 5σ of `NIntegrate` reference (complex value)

---

## Fix 4 — Linear-atom normalization (optional, improves conditioning)

**Location:** `LiftCoefficients`, before identity check

For a polynomial of the form `u + v·x^α` (single non-constant term plus a
constant, detected when `ParsePolynomial` returns exactly 2 terms and one has
zero exponent vector), optionally pre-normalize to `1 + (v/u)·x^α` before
lifting.  This is not strictly required after Fix 1 removes the false-positive,
but it improves the z0 anchor: the relative coefficient `v/u` reflects the
actual dimensionless ratio driving numerical difficulty.

Gate this behind a new option `"NormalizeAtoms" -> False` (off by default for
backward compatibility) so it doesn't silently change existing results.

---

## Test Group 25 — Higher-degree Euler integrals with complex small-magnitude coefficients

The existing lifting tests (23, 24) are all 2-variable, degree ≤ 3, uniform
Euler measure `A = {0,0}`, and single-polynomial.  Test group 25 covers the
structural gaps: 3–4 variables, degree 4–6, non-trivial Euler measure, and
multiple-polynomial integrands.

New file: `EXAMPLES/test_complex_lifted_hd.wl`

---

### Test 25A — 3-variable degree-4, single extreme complex coefficient

```
P  = 1 + (2+3I)·1e-5·x1^4 + x2^2 + x3^2 + x1·x2·x3
A  = {0, 0, 0}
B  = {-2}
```

**New territory:** first lifting test in 3 dimensions; first degree-4 monomial.

**Coefficient arithmetic:**
`|C| = √13·1e-5 ≈ 3.61e-5`;
`k* = ⌈|log(3.61e-5)|/log(1000)⌉ = ⌈10.23/6.91⌉ = 2`.
`z0 = |C|^{1/2} ≈ 6.01e-3`, residual `c = (2+3I)·1e-5 / z0^2`, `|c| = 1`.

**Lifted spec:** 4 variables, `P_lifted = 1 + c·z^2·x1^4 + x2^2 + x3^2 + x1·x2·x3`,
newton polytope in ℝ⁴.

**PASS gates:**
1. `DetectExtremeCoefficients` flags exactly 1 entry with `SuggestedK = 2`
   and `|Magnitude - √13·1e-5| < 1e-9`.
2. `LiftCoefficients` identity holds — **no** `liftidentity` message fired
   (validates Fix 1 in a 3-variable floating-point context).
3. `|residual| = 1` to 1e-9.
4. `ValidateLiftedDecomposition` `relErr < 1e-2` (reference via `NIntegrate`
   with `MaxRecursion->25, PrecisionGoal->6`).
5. `EvaluateTropicalMCLifted` (C++, 5×10⁵ samples) within 5σ of reference in
   both Re and Im.

---

### Test 25B — 3-variable degree-5, very small coefficient (k* = 3 regime)

```
P  = 1 + (1+I)·1e-8·x1^5 + x2^2 + x3^2 + x1·x2·x3
A  = {0, 0, 0}
B  = {-3}
```

**New territory:** k*=3 in a 3-variable context (previously only 2-variable in
24B); degree-5 monomial produces a high-degree vertex in the Newton polytope.

**Coefficient arithmetic:**
`|C| = √2·1e-8 ≈ 1.41e-8`;
`k* = ⌈|log(1.41e-8)|/log(1000)⌉ = ⌈18.08/6.91⌉ = 3`.
`z0 = |C|^{1/3} ≈ 2.42e-3`.

The degree-5 lifted vertex `{5,0,0,0}` in ℝ⁴ exercises polymake's 4D fan
machinery with a richer Newton polytope than any previous lifting test.

**PASS gates:**
1. `SuggestedK = 3` from `DetectExtremeCoefficients`.
2. `LiftCoefficients` identity holds to relative tolerance 1e-8 (k=3 means
   `z0^3` roundtrip — floating-point rounding is larger than for k=2).
3. `ValidateLiftedDecomposition` `relErr < 5e-2`.
4. C++ MC result within 5σ of NIntegrate reference.

---

### Test 25C — 3-variable, two polynomials, complex coefficients in each

```
P1 = 1 + (1+2I)·1e-5·x1^3 + x2^2 + x3
P2 = 1 + x1^2 + (3-I)·1e-4·x2^4 + x3^2
A  = {0, 0, 0}
B  = {-2, -1}
```

**New territory:** multi-polynomial lifting with complex extreme coefficients
living in *different* polynomials.  (Existing test 24C is multi-rule but all
rules target the same polynomial.)

**Anchor selection:**
- P1 extreme: `|C1| = √5·1e-5 ≈ 2.24e-5`, `|log|C1|| ≈ 10.7`
- P2 extreme: `|C2| = √10·1e-4 ≈ 3.16e-4`, `|log|C2|| ≈ 8.1`

Primary rule = C1 (largest `|log-mag|`).  `z0 = |C1|^{1/2} ≈ 4.73e-3`.
Residual for P2: `|c2| = |C2|/z0^2 = |C2|/|C1| = √10/√5 · 10 ≈ 14.1` — O(10),
inside the non-extreme band `[1/1000, 1000]`.

**Lift rules:**
```
{ <|"PolyIndex"->1, "ExponentVector"->{3,0,0}, "k"->2|>,
  <|"PolyIndex"->2, "ExponentVector"->{0,4,0}, "k"->2|> }
```

**PASS gates:**
1. Detection finds both extremes; P1 monomial has the larger `|log-mag|` and
   becomes the primary.
2. `z0` matches `|C1|^{1/2}` to 1e-9.
3. P2 residual magnitude in `[1, 20]` (rescaled but not extreme).
4. Identity holds for both `P1` and `P2` separately.
5. `ValidateLiftedDecomposition` `relErr < 5e-2`.

---

### Test 25D — 4-variable, non-trivial Euler measure, degree-4 complex coefficient

```
P  = 1 + (2+I)·1e-6·x1^4 + x2^2 + x3^2 + x4^2 + x1·x2·x3·x4
A  = {1/3, 0, 0, 0}      (* MonomialExponents: x1^(1/3) Euler prefactor *)
B  = {-3}
```

The integral is a **genuine generalized Euler integral**:

```
I = ∫_{[0,∞)^4}  x1^(1/3) · P(x)^{-3}  dx1 dx2 dx3 dx4
```

**New territory:** non-trivial `MonomialExponents` in a lifting test.  All
existing lifting tests use `A = {0,...,0}` (uniform measure).  The non-integer
exponent `1/3` stresses the `rawAVals = (monoExps + 1) . mMatrix` computation
in `ProcessSector` and must survive the lifting spec extension (`newMonoExps =
Append[{1/3,0,0,0}, 0] = {1/3,0,0,0,0}`).

**Coefficient arithmetic:**
`|C| = √5·1e-6 ≈ 2.24e-6`;
`k* = ⌈|log(2.24e-6)|/log(1000)⌉ = ⌈13.01/6.91⌉ = 2`.

**PASS gates:**
1. `SuggestedK = 2`.
2. Lifted spec `MonomialExponents = {1/3, 0, 0, 0, 0}` — non-trivial Euler
   exponent is preserved and appended with 0 for the auxiliary variable.
3. `ValidateLiftedDecomposition` `relErr < 5e-2`; the reference NIntegrate
   must include the `x1^(1/3)` weight.
4. C++ MC result within 5σ of NIntegrate reference.

---

### Test 25E — 3-variable stress test: degree-6, three competing complex extremes

```
P  = 1 + (1+I)·ε·x1^6 + (2+I)·ε·x1^3·x2^3 + (1-3I)·ε·x2^6
     + x1^2 + x2^2 + x3^2 + x1·x2·x3
ε  = 3·1e-5
A  = {0, 0, 0}
B  = {-2}
```

**New territory:** three extreme complex coefficients of *comparable* magnitude
lifted against a single anchor; degree-6 Newton polytope vertices `{6,0,0}`,
`{3,3,0}`, `{0,6,0}` in ℝ³ produce a rich fan; tests primary-rule selection
under a near-flat `|log|C||` spectrum.

**Magnitude spectrum:**
| Monomial | `|C|` | `|log\|C\||` |
|----------|--------|--------------|
| `x1^6` | `√2·ε ≈ 4.24e-5` | 10.07 ← primary |
| `x1^3·x2^3` | `√5·ε ≈ 6.71e-5` | 9.61 |
| `x2^6` | `√10·ε ≈ 9.49e-5` | 9.26 |

`k* = 2` for all three.  `z0 = (√2·ε)^{1/2}`.
Residuals:
- `c1 = (1+I)ε / z0^2 = (1+I)/√2`, `|c1| = 1`
- `c2 = (2+I)ε / z0^2 = (2+I)/√2`, `|c2| = √5/√2 ≈ 1.58`
- `c3 = (1-3I)ε / z0^2 = (1-3I)/√2`, `|c3| = √10/√2 ≈ 2.24`

All residuals are O(1) by construction.

**Lift rules:**
```
{ <|"PolyIndex"->1, "ExponentVector"->{6,0,0}, "k"->2|>,
  <|"PolyIndex"->1, "ExponentVector"->{3,3,0}, "k"->2|>,
  <|"PolyIndex"->1, "ExponentVector"->{0,6,0}, "k"->2|> }
```

**PASS gates:**
1. Detection finds all 3 extreme monomials with `SuggestedK = 2`.
2. Primary rule = `x1^6` monomial (largest `|log|C||`); `z0` matches to 1e-9.
3. All 3 residual magnitudes in `[0.5, 3]`.
4. Identity holds for all 3 monomials simultaneously.
5. `ValidateLiftedDecomposition` `relErr < 1e-1` (3-variable high-degree
   polytope has many sectors; NIntegrate convergence slower).
6. C++ MC result within 5σ of NIntegrate reference.

---

### Summary of coverage gaps closed by Test 25

| Dimension | Degree | Multi-poly | Non-trivial A | Multi-rule | Previously tested |
|-----------|--------|------------|---------------|------------|-------------------|
| 2 | ≤3 | no | no | 24C (yes) | Tests 24A–C |
| **3** | **4** | no | no | no | **25A** |
| **3** | **5** | no | no | no | **25B** |
| **3** | **3/4** | **yes** | no | **yes** | **25C** |
| **4** | **4** | no | **yes** | no | **25D** |
| **3** | **6** | no | no | **yes (3-way)** | **25E** |

---

## Implementation order

| Priority | Fix | Risk | Scope |
|----------|-----|------|-------|
| 1 (quick) | Fix 1 — liftidentity tolerance | Low | ~15 lines, `tropical_eval.wl` |
| 2 (quick) | Fix 2 — pre-flight detection | Low | ~25 lines, `tropical_eval.wl` |
| 3 (core)  | Fix 3a–3c — option + spec split | Medium | ~60 lines, `tropical_eval.wl` |
| 4 (core)  | Fix 3d — C++ generator | Medium | ~40 lines, codegen section |
| 5 (tests) | Test 25A–B — 3D, degree 4–5 | Low | new file, validates Fixes 1–2 |
| 6 (tests) | Test 25C–E — multi-poly, Euler measure, degree-6 | Medium | same file, validates Fix 3 |
| 7 (optional) | Fix 4 — atom normalization | Low | ~20 lines, `tropical_eval.wl` |

Fixes 1 and 2 are independent and can be done in a single commit.  Fix 3
depends on Fix 2 (shared pre-flight logic) but not on Fix 1.  Tests 25A–B can
be written and run before Fix 3 is complete (they use real exponents); 25C–E
also use real exponents and do not require the SplitRealImag mode.

---

## Files touched

- `tropical_eval.wl`: all four fixes
- `EXAMPLES/test_complex_lifted.wl`: cases 24E and 24F
- `EXAMPLES/test_complex_lifted_hd.wl`: cases 25A–E (new file)
