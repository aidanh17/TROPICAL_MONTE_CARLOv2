# Implementation Plan — VEGAS integrator for the tropical pipeline

**Goal.** Add the option to integrate the *flattened* tropical sectors with **CUBA‑Vegas**
instead of the shipped plain Monte‑Carlo sampler, and add EXAMPLES that run the
**same specs both ways (MC and VEGAS) and demonstrate agreement** (and agreement
with the analytic / NIntegrate reference).

**Why (from `TEST/RESULTS.md` / `TEST/REPORT.md`).** On the post‑decomposition,
post‑flattening integrand, plain MC is the weakest sampler. CUBA‑Vegas wins in
every tested case ($n=2\dots8$), is the only method that scales to $n=7,8$, pairs
the QMC convergence rate ($\sim N^{-1.0\dots1.7}$ vs MC's $\sim N^{-0.5}$) with
adaptive importance sampling that absorbs the $\sqrt{\cdot}$ cube‑edge structure
flattening leaves behind, and parallelizes over kinematic points exactly like the
current code. The decomposition produces the *value*; the sampler is a
second‑order optimization on top — so VEGAS must return the **same answer** as MC,
only more accurately. That equality is exactly what the new examples assert.

This is a **plan only** — do not implement from this document without the
orchestrator dispatching the phases below.

---

## 0. Orientation — files and exact integration points

| File | Role | Key anchors (current line numbers) |
|---|---|---|
| `tropical_eval.wl` | the package | see below |
| `tropical_eval.wl` → `GenerateMonomialSumCpp` | emits per‑poly monomial sum C++ | `1484` |
| `tropical_eval.wl` → `GenerateCppMonteCarlo` | **emits the whole `.cpp`**; `Options` at `1527`, def at `1533`; includes block `1647–1661`; `cx` typedef `1661`; integrand functions `1663–1666`; `integrand_table[]` / `integrand_dim[]` / `N_*` consts `1668–1689`; **`main()` `1692–1840`**; per‑sector MC sampling loop **`1753–1810`**; return assoc `1859–1862` | `1533` |
| `tropical_eval.wl` → `CompileCpp` | `g++` invocation; flags `1874–1881`; OpenMP fallback `1886–1895` | `1869` |
| `tropical_eval.wl` → `EvaluateTropicalMC` | driver; `Options` `1917`; codegen call `2109–2115`; debug compile+run `2123–2154`; release compile `2176`; run `2190–2191`; result parse `2204–2239` | `1928` |
| `tropical_eval.wl` → `EvaluateTropicalMCLifted` | wrapper; `Options` = `Join[..., Options[EvaluateTropicalMC]]` `2327`; pass‑through via `FilterRules` `2352` | `2336` |
| `EXAMPLES/cuba_common.wl` → `cubaFindPrefix[]` | CUBA prefix auto‑detect (reuse) | `33` |
| `TEST/bench.cpp` → `cubaWrap`, `run_cuba`, `cuba_one` | **reference single‑sector Vegas call** (ncomp=2, Re/Im) | wrapper `149`, Vegas call `168` |
| `TEST/bench_batch.cpp` → `cb_ncomp`, `run_cuba_ncomp` | **reference `ncomp`‑batched Vegas** (shared samples across kp), `CUBA_MAXCOMP=512` | `177`, `202`, cap `167` |
| `TEST/gen_kin.wl` | how shared‑basis `KS_*` tables are built (for Phase 4) | `69–119` |

**Verified environment (this machine).**
- CUBA 4.x present: header `/opt/homebrew/include/cuba.h`, static lib `/opt/homebrew/lib/libcuba.a`.
- `Vegas(...)` signature (cuba.h): `void Vegas(int ndim, int ncomp, integrand_t, void* userdata, int nvec, cubareal epsrel, cubareal epsabs, int flags, int seed, int mineval, int maxeval, int nstart, int nincrease, int nbatch, int gridno, const char* statefile, void* spin, int* neval, int* fail, cubareal integral[], cubareal error[], cubareal prob[]);`
- `integrand_t` = `int (*)(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)`. `cubareal` = `double` by default.
- Compiler is **Apple clang** (`/usr/bin/g++`): **no `-fopenmp`**. `CompileCpp` already retries without `-fopenmp` on this error, so the VEGAS kp loop runs **serially** here (which conveniently sidesteps CUBA thread‑safety on this box; see §2.4). On Linux/gcc it parallelizes.

**Invariants that MUST be preserved.**
1. `Integrator -> "MC"` is the **default**. With the default, the generated `.cpp`,
   the compile flags, and all results must be **byte‑identical** to today (regression gate, §5).
2. The 4‑number per‑kp output contract `re im reErr imErr` (written `tropical_eval.wl:1832–1834`,
   parsed `2204–2239`) is **unchanged** — VEGAS writes the same four columns
   (`integral[0] integral[1] error[0] error[1]`, summed in quadrature over sectors).
   So the entire result‑reading / `finalResults` assembly path is untouched.
3. Complex integrands are first‑class (e.g. Example 17 has `B = -(2+I/2)`); the
   VEGAS path uses **`ncomp = 2`** (component 0 = Re, component 1 = Im) — never assume real.

---

## 1. Public API design (decide once, thread everywhere)

Add an `"Integrator"` option to the three public entry points; default `"MC"`.
Accepted values for this plan: `"MC"` (current) and `"VEGAS"`. (Leave room for a
future `"QMC"`; **out of scope** here — note it in usage text only.)

### 1.1 `GenerateCppMonteCarlo` — new options (extend `Options` at `tropical_eval.wl:1527`)
```wl
Options[GenerateCppMonteCarlo] = {
  "NSamples"      -> 1000000,
  "MaxDim"        -> 20,
  "SeedBase"      -> 42,
  "Integrator"    -> "MC",            (* "MC" | "VEGAS" *)
  (* VEGAS tuning — used only when Integrator == "VEGAS" *)
  "VegasEpsRel"   -> Automatic,       (* Automatic => 10^-PrecisionGoal via driver; else a number *)
  "VegasSeed"     -> 0,               (* 0 => Sobol low-discrepancy (recommended by REPORT.md) *)
  "VegasNStart"   -> 1000,
  "VegasNIncrease"-> 500,
  "VegasNBatch"   -> 1000,
  "VegasMinEval"  -> 0
};
```
- Validate `"Integrator"` against `{"MC","VEGAS"}`; on bad value `Message[...]; Return[$Failed]`.
- When `Automatic`, resolve `VegasEpsRel` at the **driver** level from `PrecisionGoal`
  (e.g. `10^-precGoal`, clamped to `[1e-12, 1e-2]`); pass a concrete number into the generator.
- The generator’s return association (`1859–1862`) gains two keys:
  `"Integrator" -> <value>` and `"NeedsCuba" -> (Integrator === "VEGAS")`.
  The driver reads `NeedsCuba` to decide linking.

### 1.2 `EvaluateTropicalMC` — new options (extend `Options` at `tropical_eval.wl:1917`)
Add the same six keys (`"Integrator"`, `"VegasEpsRel"`, `"VegasSeed"`,
`"VegasNStart"`, `"VegasNIncrease"`, `"VegasNBatch"`, `"VegasMinEval"`). Driver responsibilities:
- Read `integrator = OptionValue["Integrator"]`.
- If `integrator === "VEGAS"`: call the new private `findCubaPrefix[]` (§2.1).
  If `$Failed`, `Message[TropicalEval::nocuba]; Return[$Failed]` (do **not** silently
  fall back to MC — the user explicitly asked for VEGAS).
- Resolve `VegasEpsRel` from `Automatic` using `precGoal` here, then forward all
  VEGAS options + `"Integrator"` into the `GenerateCppMonteCarlo[...]` call at `2109`.
- Forward the CUBA prefix into **both** `CompileCpp` calls (debug `2128`, release `2176`).
- The `n_samples` runtime argument (passed at `2139`/`2190`) is **reused unchanged** —
  for VEGAS it is interpreted by the generated `main()` as **`maxeval` per sector** (§2.2).
  No change to the argv contract.

### 1.3 `EvaluateTropicalMCLifted` (`tropical_eval.wl:2327`)
No code change needed for option threading: its `Options` already
`Join[..., Options[EvaluateTropicalMC]]` and it forwards via
`FilterRules[{opts}, Options[EvaluateTropicalMC]]` (`2352`). Adding the keys to
`EvaluateTropicalMC` makes them flow through automatically. **Add a test** that
`EvaluateTropicalMCLifted[..., "Integrator" -> "VEGAS"]` actually reaches the VEGAS
codegen (the lifted sectors are ordinary convergent sectors with a `DomainConstraint`
indicator at `tropical_eval.wl:1586–1612`, which is integrator‑agnostic — it lives
inside `integrand_conv_k`, so VEGAS gets it for free).

### 1.4 New messages (add near `tropical_eval.wl:148–164`)
```wl
TropicalEval::nocuba = "Integrator -> \"VEGAS\" requires the CUBA library (cuba.h + libcuba). \
Not found under /opt/homebrew, /usr/local, or /usr. Install via `brew install cuba` or \
https://feynarts.de/cuba/, or use Integrator -> \"MC\".";
TropicalEval::badintegrator = "Unknown Integrator `1`; expected \"MC\" or \"VEGAS\".";
```

---

## 2. Code generation — the VEGAS branch of `GenerateCppMonteCarlo`

Everything up to and including the integrand functions, `integrand_table[]`,
`integrand_dim[]`, `N_INTEGRANDS`, `N_PARAMS`, `MAX_DIM` (lines `1563–1689`) is
**shared** and unchanged. Only the includes block and `main()` branch on `Integrator`.

### 2.1 Move CUBA detection into the package (private helper)
Port `cubaFindPrefix[]` from `EXAMPLES/cuba_common.wl:33` to a private symbol in
`tropical_eval.wl` (e.g. `findCubaPrefix[]`), returning the prefix dir or `$Failed`.
Leave the copy in `cuba_common.wl` as‑is (that file is standalone for the raw
cross‑check) — just avoid divergence by keeping the logic identical.

### 2.2 Includes + helpers (emit only when `Integrator == "VEGAS"`)
After the existing includes (`1647–1661`) and the `using cx` typedef (`1661`), add:
```cpp
extern "C" {
#include <cuba.h>
}
// VEGAS tuning (compile-time constants, baked from options)
static const double VEGAS_EPSREL    = <VegasEpsRel resolved to a number>;
static const double VEGAS_EPSABS    = 1e-300;
static const int    VEGAS_SEED      = <VegasSeed>;       // 0 => Sobol (recommended)
static const int    VEGAS_NSTART    = <VegasNStart>;
static const int    VEGAS_NINCREASE = <VegasNIncrease>;
static const int    VEGAS_NBATCH    = <VegasNBatch>;
static const int    VEGAS_MINEVAL   = <VegasMinEval>;

struct SectorCtx { IntegrandFunc fn; const double* params; };
static int vegas_wrap(const int* ndim, const cubareal x[], const int* /*ncomp*/,
                      cubareal ff[], void* ud) {
  const SectorCtx* ctx = static_cast<const SectorCtx*>(ud);
  double y[MAX_DIM];
  for (int i = 0; i < *ndim; ++i) y[i] = (double)x[i];
  cx v = ctx->fn(y, ctx->params);
  ff[0] = v.real();
  ff[1] = v.imag();
  return 0;
}
```
Note: `IntegrandFunc` is typedef'd at `1668–1670`; emit the `using IntegrandFunc`
line *before* `vegas_wrap` in the VEGAS branch (currently it is emitted after the
integrand table — move/duplicate so `vegas_wrap` can name the type, or forward‑declare).

### 2.3 `main()` — VEGAS variant
Keep argv parsing, OpenMP setup, kinematic‑file reading, and result writing
**identical** to the MC `main()` (`1692–1741`, `1825–1840`). Replace only the
compute region (the MC body `1743–1823`) with:
```cpp
    const int zero = 0;
    cubacores(&zero, &zero);   // single process: deterministic, macOS-safe

    #pragma omp parallel for schedule(dynamic)
    for (int kp = 0; kp < n_kp; kp++) {
        static double dummy_params[1] = {0.0};
        const double* params = (N_PARAMS > 0) ? kinematic_data[kp].data() : dummy_params;

        double total_re = 0.0, total_im = 0.0;
        double total_var_re = 0.0, total_var_im = 0.0;

        for (int s = 0; s < N_INTEGRANDS; s++) {
            int dim = integrand_dim[s];
            SectorCtx ctx{ integrand_table[s], params };
            int neval = 0, fail = 0;
            cubareal integ[2], err[2], prob[2];
            Vegas(dim, 2, vegas_wrap, &ctx, 1,
                  VEGAS_EPSREL, VEGAS_EPSABS, 0 /*flags*/, VEGAS_SEED,
                  VEGAS_MINEVAL, n_samples /*maxeval per sector*/,
                  VEGAS_NSTART, VEGAS_NINCREASE, VEGAS_NBATCH,
                  0 /*gridno*/, nullptr, nullptr,
                  &neval, &fail, integ, err, prob);
            total_re     += integ[0]; total_im     += integ[1];
            total_var_re += err[0]*err[0]; total_var_im += err[1]*err[1];
#ifdef TROPICAL_MC_DEBUG
            if (kp == 0)
                std::cerr << "Sector " << s << " (Vegas): est=(" << integ[0] << ","
                          << integ[1] << ") err=(" << err[0] << "," << err[1]
                          << ") fail=" << fail << " prob=" << prob[0]
                          << " neval=" << neval << std::endl;
#endif
        }
        results[kp] = {total_re, total_im,
                       std::sqrt(total_var_re), std::sqrt(total_var_im)};
    }
```
- Output identical to MC (same `results[kp]` 4‑tuple → same write loop).
- The MC‑specific debug instrumentation (`nan_count`, `max_mag`, Welford — `1758–1804`)
  is **not** emitted in the VEGAS branch; the small VEGAS debug print above replaces it.
- `seed = 0` makes each `Vegas` call use Sobol (deterministic) ⇒ results are
  **independent of thread count** (each kp is independent). This is what the §5
  thread‑safety check exploits.

### 2.4 Thread‑safety note (must be verified, not assumed)
Concurrent `Vegas` calls from multiple OpenMP threads, with `cubacores(0,0)`,
`gridno=0`, `statefile=nullptr`, are believed reentrant (per‑call state, no stored
grid), but CUBA is not formally documented thread‑safe. **Verification (Phase 5):**
on a gcc/`-fopenmp` build, run the same spec with 1 thread and with N threads and
require bitwise‑equal results (legitimate because seed=0 ⇒ deterministic per kp).
If they differ, the fallback is to wrap the `Vegas(...)` call in
`#pragma omp critical` (serializes only the library call; the integrand work still
parallelizes poorly — acceptable) **or** emit the kp loop without the `parallel for`
when `Integrator==VEGAS`. On this Mac (Apple clang, no OpenMP) the loop is serial,
so the check is a Linux/gcc CI item.

---

## 3. `CompileCpp` — link CUBA when needed

Smallest non‑breaking change: add an **optional 4th positional argument**
`cubaPrefix_: None` to `CompileCpp` (`tropical_eval.wl:1869`). Existing 3‑arg calls
are unaffected.
- When `cubaPrefix` is a string, append to `flags`:
  `"-I" <> FileNameJoin[{cubaPrefix,"include"}]` and to the link line
  `"-L" <> FileNameJoin[{cubaPrefix,"lib"}], "-lcuba"` (before `-lm`).
  (Static `libcuba.a` links fine.)
- Keep the existing `-fopenmp` retry‑without logic (`1886–1895`) — it already handles
  Apple clang. CUBA flags must be preserved across that retry.
- `EvaluateTropicalMC` passes `findCubaPrefix[]` as the 4th arg to both `CompileCpp`
  calls (`2128`, `2176`) **only when** `integrator === "VEGAS"`; otherwise passes
  `None` (or omits) so the MC compile line is byte‑identical to today.

Reference compile line (matches `TEST/run.sh`):
```
g++ -std=c++17 -O3 -I/opt/homebrew/include tropical_mc_generated.cpp \
    -L/opt/homebrew/lib -lcuba -lm -o tropical_mc
```

---

## 4. (OPTIONAL / advanced) `ncomp`‑batched VEGAS for kinematic scans

This is the headline **batch** lever from `TEST/RESULTS.md §8`
(“~30–60× faster to ~1 % at n=8”). It is **optional** — Phases 1–3 + 5 already deliver
“the option to use VEGAS (batched over kinematic points by the existing kp loop)”.
Implement this only if the high‑throughput kinematic‑scan path is wanted.

**Idea (ref `TEST/bench_batch.cpp:177–227`).** One `Vegas` call per *sector* handles a
whole chunk of kp at once via `ncomp`, sharing the sample points (the expensive
transcendentals are coefficient‑independent). For our **complex** integrand use **2
components per kp**: `ff[2*c]=Re`, `ff[2*c+1]=Im` ⇒ `ncomp = 2*chunk`. Cap
`chunk ≤ 256` so `ncomp ≤ 512` (`bench_batch.cpp:167` documents 1024 segfaults Vegas at
n=8; 512 is safe through n=8).

Gate behind `"Integrator" -> "VEGASBATCH"` (a third value) or a sub‑option
`"VegasBatch" -> True`. Generated `main()`:
- For each sector `s`: loop `k0` over kp in steps of `chunk`; set a context giving the
  base pointer `kinematic_data.data()`, `k0`, `cs=min(chunk, n_kp-k0)`; call
  `Vegas(dim, 2*cs, vegas_wrap_ncomp, &ctx, 1, ...)`; scatter
  `integ[2*c],integ[2*c+1]` into `results[k0+c]` and accumulate `err^2`.
- `vegas_wrap_ncomp` computes `y` once, then loops `c=0..cs-1`,
  `cx v = fn(y, base + (k0+c)*N_PARAMS); ff[2*c]=v.real(); ff[2*c+1]=v.imag();`.
- Parallelize over sectors and/or chunks with `#pragma omp parallel for` (same
  thread‑safety caveat as §2.4).
- (Further, optional) the shared coefficient‑independent **monomial basis** reuse
  (`bench_batch.cpp:124–162`, tables built in `gen_kin.wl:69–119`) is a second ~4×
  lever but requires emitting `KS_*` tables from `GenerateCppMonteCarlo`. Defer unless
  needed; note it in docs.

**Recommendation:** ship Phases 1–3 + 5 first (correctness + parity examples), land
Phase 4 as a follow‑up PR.

---

## 5. Validation, regression, and acceptance criteria

**R1 — MC regression (must pass).** For a representative spec, generate the `.cpp`
with `Integrator->"MC"` before and after the change; assert the file is
**byte‑identical** to the pre‑change output (and the compile line unchanged). This
guarantees the default path is untouched.

**R2 — VEGAS smoke.** Generate + compile + run the VEGAS `.cpp` for a small
param‑free spec; assert exit 0 and 4 finite numbers per kp.

**R3 — CUBA‑absent behavior.** Temporarily simulate no CUBA (point `findCubaPrefix`
at empty dirs); assert `EvaluateTropicalMC[..., "Integrator"->"VEGAS"]` emits
`TropicalEval::nocuba` and returns `$Failed` (no silent MC fallback).

**R4 — thread‑safety (gcc/OpenMP only).** §2.4: 1‑thread vs N‑thread bitwise equality.

**R5 — MC≈VEGAS≈truth (the examples, §6).** For every example case, require:
- `|MC − VEGAS|  ≤  K * sqrt(MCerr² + VEGASerr²)` with `K = 5` (combined‑error band), **and**
- both within a budget‑appropriate relative tolerance of the analytic/NIntegrate
  reference (e.g. `≤ 1e-2` at the example budgets; VEGAS will typically be ~1–3 orders
  tighter — that is the point of the demo, not a failure).

**Definition of done.** R1–R3 + R5 green on this Mac (serial VEGAS); R4 green on a
gcc/OpenMP runner (or documented as deferred); docs updated (§7); no new warnings
from `CompileCpp` on the MC or VEGAS source.

---

## 6. EXAMPLES — run each spec BOTH ways (MC and VEGAS)

Create **`EXAMPLES/test_vegas.wl`** (mirror the structure/printing of
`EXAMPLES/test_cuba.wl`). Load the package via the same idiom other examples use
(`Get[FileNameJoin[{ParentDirectory[], "tropical_eval.wl"}]]`), write generated
files under `EXAMPLES/INTERFILES`. Provide a helper:

```wl
compareMCvsVEGAS[label_, spec_, fanData_, kinPoints_, exactList_,
                 mcSamples_, vegasMaxeval_] := Module[{rMC, rVE, ...},
  rMC = EvaluateTropicalMC[spec, fanData, kinPoints,
          "Integrator" -> "MC",    "NSamples" -> mcSamples,
          "RunChecks" -> False, "Verbose" -> False,
          "WorkingDirectory" -> FileNameJoin[{Directory[], "INTERFILES"}]];
  rVE = EvaluateTropicalMC[spec, fanData, kinPoints,
          "Integrator" -> "VEGAS", "NSamples" -> vegasMaxeval,
          "RunChecks" -> False, "Verbose" -> False,
          "WorkingDirectory" -> FileNameJoin[{Directory[], "INTERFILES"}]];
  (* print a table: kp | MC ± err | VEGAS ± err | exact/NI | |MC-VEGAS| | relErr_MC | relErr_VE *)
  (* assert R5 for each kp; collect PASS/FAIL *)
];
```

Cases (each must print a side‑by‑side MC vs VEGAS table and assert R5):

- **V1 — A_simplex, param‑free, n=2 and n=4.** `P=1+∑x_i`, `A_i=0`, `B=-(n+2)`;
  exact `1/(n+1)!` (`TEST/RESULTS.md` table). Fan from `PolytopeVertices[P^-1, vars]`.
  Shows VEGAS reproduces MC and the closed form; n=4 is where MC already struggles.
- **V2 — Af_frac, n=4 (the non‑smooth flattened edge).** `A_i=-1/2`, `B=-(n+1)`,
  exact `π^{n/2} Γ(n/2+1)/n!`. This is the `√(·)` cube‑edge case (`RESULTS.md §4`)
  where MC is weakest and VEGAS shines — strongest demonstration of “same value,
  better accuracy”. (`A_i=-1/2` ⇒ set `MonomialExponents -> {-1/2,...}`.)
- **V3 — complex exponent, kinematic scan (reuse Example 17).**
  `P = 1 + λ x1² + x2² + x1 x2²`, `B = -(2 + I/2)`, `KinematicSymbols={λ}`,
  `kinPoints = {{0.5},{1.},{2.},{4.},{8.}}`; reference via `NIntegrate` (as Example 17,
  `tropical_eval_examples2.wl:216–243`). Exercises **complex** (`ncomp=2`) and the
  **multi‑kp batch loop**; asserts MC≈VEGAS≈NIntegrate for both Re and Im.
- **V4 — coefficient batch scan (the batch regime).** `P = c0+∑c_i x_i`, `n=4`,
  `B=-6`, `KinematicSymbols={c0..c4}`, exact `1/(120 c0² ∏c_i)` (`gen_kin.wl`).
  Use ~24 random coefficient sets (keep runtime modest). Demonstrates MC and VEGAS
  agree across a scan; if Phase 4 lands, add a `"VEGASBATCH"` column here.

Also append a short pointer to `EXAMPLES/` listing (and to the examples index if one
exists) noting `test_vegas.wl` and that it requires CUBA.

---

## 7. Documentation to keep in sync (per repo memory: docs must track behavior)

- **Usage strings** (`tropical_eval.wl:109–134`): document `"Integrator"` and the
  `Vegas*` options on `GenerateCppMonteCarlo`, `EvaluateTropicalMC`,
  `EvaluateTropicalMCLifted`; note CUBA is required for `"VEGAS"`.
- **`MANUAL/manual.tex`** (then rebuild `manual.pdf`): add a subsection “Numerical
  integrator: MC vs VEGAS” summarizing the `TEST/REPORT.md` findings (why VEGAS,
  the `√` edge, dimension scaling, the `ncomp` batch lever) and the new API. Point to
  `TEST/REPORT.md` and `EXAMPLES/test_vegas.wl`.
- **`SUMMARY.txt`**: note the new option and example.
- **`TEST/README.md`**: one line that the study’s recommendation is now implemented
  behind `Integrator -> "VEGAS"`.

---

## 8. Suggested orchestration (phases & rough sizing)

| Phase | Scope | Touches | Depends on |
|---|---|---|---|
| **P1** | API + messages + `findCubaPrefix` + `CompileCpp` CUBA flag | `tropical_eval.wl` (§1, §2.1, §3) | — |
| **P2** | VEGAS codegen branch in `GenerateCppMonteCarlo` + driver wiring | `tropical_eval.wl` (§2.2–2.3, §1.2) | P1 |
| **P3** | `EXAMPLES/test_vegas.wl` V1–V4 + helper | new file (§6) | P2 |
| **P4** | *(optional)* `ncomp`‑batched VEGAS (+ optional shared basis) | `tropical_eval.wl` (§4), example column | P2 |
| **P5** | Validation R1–R5, docs §7 | tests, `manual.tex`, `SUMMARY.txt` | P2 (R4 needs gcc) |

P1–P3 + P5 = MVP that fully satisfies the request. P1 and the V‑case spec authoring
in P3 can start in parallel; P2 must precede running P3. Gate the merge on R1
(byte‑identical MC) and R5 (MC≈VEGAS≈truth).
```
