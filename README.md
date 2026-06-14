# Tropical Monte Carlo v2

A numerical integration pipeline for **convergent generalized Euler integrals**

$$I \;=\; \int_{[0,\infty)^n} dx_1\cdots dx_n \; \prod_i x_i^{A_i} \; \prod_j P_j(\mathbf{x})^{B_j}$$

using **tropical geometry** (the Newton-polytope fan) to partition the domain
into simplicial sectors, a symbolic *flattening* that renders each sector
integrand $\mathcal{O}(1)$ on the unit cube $[0,1]^n$, and **uniform Monte Carlo**
(or adaptive **VEGAS**) sampling. A single C++ compilation evaluates thousands of
kinematic points in parallel, because the fan depends only on the exponent
structure, not on the coefficient values.

These integrals arise in quantum field theory (Feynman integrals), cosmological
correlators, and precision physics, where they are high-dimensional and defeat
standard quadrature.

> **Full reference:** [`MANUAL/manual.pdf`](MANUAL/manual.pdf) (source
> [`MANUAL/manual.tex`](MANUAL/manual.tex)) — complete algorithm, data structures,
> options, and worked examples.
> **Developer notes / function reference:** [`SUMMARY.txt`](SUMMARY.txt).

---

## Features

- **Exact tropical decomposition** — fan decomposition + tropical factoring +
  flattening are exact algebraic changes of variables; the Monte Carlo only
  estimates each flat $\mathcal{O}(1)$ sector integral.
- **Complex exponents** — the monomial exponents $A_i$ and polynomial exponents
  $B_j$ may be complex; flattening becomes a contour rotation that also tames the
  oscillation (manual §"Complex Flattening as a Contour Rotation").
- **Complex coefficients** — polynomial coefficients may be complex (carried as
  `std::complex<double>` throughout), including through the lifting module
  (see below).
- **Auxiliary-variable lifting** — for coefficients spanning many orders of
  magnitude (very large *or* very small), lifting moves the extreme magnitude
  into the geometry, cutting Monte Carlo variance by an order of magnitude or
  more. Works for **complex** coefficients too: the magnitude becomes a real
  anchor $z_0=|C|^{1/k}$ and the phase rides along in an $\mathcal{O}(1)$
  residual.
- **Two integrators** — plain uniform `"MC"` (default) and CUBA `"VEGAS"`
  (adaptive importance sampling; 1–3 orders tighter on the flattened integrand).
- **Parallel kinematic scans** — one compile, OpenMP over kinematic points,
  deterministic per-point RNG seeding.
- **Validation built in** — every decomposition is cross-checked against direct
  `NIntegrate`, and the test suite additionally cross-checks against an exact
  Gamma closed form and the independent **CUBA** library.

**Scope (v2):** convergent integrals only. `RegulatorSymbol` must be `None`, and
every sector must have $\mathrm{Re}(a_{\mathrm{eff}})>0$; otherwise the driver
rejects the input (`TropicalEval::noregulator` / `::divergentinput`). For
$\varepsilon$-regulated / divergent integrals (pole extraction by tropical
subtraction or IBP), use the original `TROPICAL_MONTE_CARLO` package.

---

## Dependencies

| Tool | Required? | Used for |
|------|-----------|----------|
| **Wolfram Mathematica 12+** (`wolframscript`) | yes | symbolic engine; runs all `.wl` files |
| **g++** (C++17; OpenMP optional) | for MC runs | compiles/executes the generated integrand |
| **Polymake** | only for fan computation | builds the tropical fan (auto-detected under `/opt/homebrew`, `/usr/local`, `/usr`). Skippable with `$SkipPolymakeLoad = True` + explicit `fanData` |
| **CUBA** | optional | `Integrator -> "VEGAS"` and the independent cross-check examples (`brew install cuba`) |

---

## Quick start

```mathematica
(* from the repository root *)
Get["tropical_eval.wl"];

poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2;
vars = {x[1], x[2]};

spec = <|
  "Polynomials"         -> {poly},
  "MonomialExponents"   -> {0, 0},
  "PolynomialExponents" -> {-2},
  "Variables"           -> vars,
  "KinematicSymbols"    -> {},
  "RegulatorSymbol"     -> None
|>;

(* fan from the Newton polytope (coefficients are irrelevant to the fan) *)
verts   = PolytopeVertices[poly^(-1), vars];
fanData = ComputeDecomposition[verts];

(* exact-decomposition cross-check, then full Monte Carlo *)
ValidateDecomposition[spec, fanData, {}, 3]["RelativeError"]
EvaluateTropicalMC[spec, fanData, {}]
```

**Lifting a small (or large) coefficient** — let the `kStar` anchor rule pick `k`
automatically:

```mathematica
spec = <|"Polynomials" -> {1 + 10^-4 x[1]^2 + x[2]^2 + x[1] x[2]^2},
         "MonomialExponents" -> {0, 0}, "PolynomialExponents" -> {-2},
         "Variables" -> {x[1], x[2]}, "KinematicSymbols" -> {},
         "RegulatorSymbol" -> None|>;
EvaluateTropicalMCLifted[spec, {{}}, "LiftRules" -> Automatic, "NSamples" -> 10^6]
```

The same call works with a **complex** small coefficient, e.g.
`(1 + I) 10^-4 x[1]^2` — the result is correctly complex-valued
(see Examples 24–27).

---

## Examples and tests (`EXAMPLES/`)

Run any of these with `wolframscript -file EXAMPLES/<file>` from the repo root.

| File | Contents |
|------|----------|
| `tropical_eval_examples.wl` | Examples 1–14: core pipeline (2D–4D, kinematic scans, numerators); also *defines* the 6-test suite |
| `run_validation_suite.wl` | Examples 1–14 then `RunAllTests[]` (Tests 1, 2, 3v2, 5, 6, 7); exits non-zero on failure |
| `tropical_eval_examples2.wl` | Examples 15–20: complex exponents, small coefficients, real-coefficient lifting |
| `tropical_eval_examples3.wl` | Examples 21–23: independent cross-checks (exact Gamma closed form + CUBA) |
| `tropical_eval_examples4.wl` | **Examples 24–27: lifting with complex, small-magnitude coefficients** |
| `test_lifted.wl` | Test 23: lifted-integral pipeline (exactness, end-to-end C++, error paths, `EmptyDomain`) |
| `test_complex_lifted.wl` | **Test 24: lifting with complex coefficients (PASS/FAIL; exits non-zero on failure)** |
| `test_complex_lifted_hd.wl` | Test 25: higher-dimensional lifting (3–4 vars) with complex coefficients |
| `test_complex_lifted_vegas_hd.wl` | **Test 26: an 8-D integral with a tiny complex coefficient under VEGAS, with and without lifting, vs an exact reference (requires CUBA)** |
| `test_cuba.wl` | CUBA Cuhre cross-check of every convergent suite integral (requires CUBA) |
| `test_vegas.wl` | MC ↔ VEGAS parity demonstration (requires CUBA) |
| `test_quick.wl`, `test_cpp.wl`, `test_seedbase.wl` | smoke / codegen / RNG-seed tests |

> **Note:** scripts that run the C++ Monte Carlo share `INTERFILES/` with fixed
> filenames — run them **one at a time** in the same working directory.

---

## Repository layout

```
TROPICAL_MONTE_CARLOv2/
├── tropical_eval.wl      Main evaluation package (ProcessSector, lifting,
│                         C++ codegen, EvaluateTropicalMC[ Lifted ])
├── tropical_fan.wl       Tropical fan computation via Polymake (black box)
├── README.md             This file
├── SUMMARY.txt           Developer notes & function reference
├── MANUAL/               Comprehensive reference manual (manual.pdf / .tex)
├── EXAMPLES/             Worked examples and test suites (table above)
├── SANDBOX/              Lifting development experiments & benchmarks
├── TEST/                 Sampler/batching study (MC vs QMC vs CUBA; see TEST/README.md)
├── CC/                   Historical FIESTA cross-checks (original package)
└── INTERFILES/           Generated C++/binaries/results (gitignored)
```
