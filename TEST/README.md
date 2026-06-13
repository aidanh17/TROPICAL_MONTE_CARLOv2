# TEST/ — sampler & batching study on the flattened tropical integrand

How best to integrate the *pre‑processed* (tropically decomposed + flattened)
sectors — **plain MC, quasi‑MC, or CUBA** — and how that changes with dimension
and in batch (kinematic‑scan) mode. Full write‑up: **`REPORT.md`** (and the
typeset **`sum/report.pdf`**); condensed: **`RESULTS.md`**.

**TL;DR**
- **Single integral:** CUBA‑**Vegas** wins (best accuracy/eval, rate stable in
  dimension); Sobol **QMC** is the cheap dependency‑light upgrade over MC; plain
  **MC** is weakest; **Cuhre** only shines when flattening leaves a *smooth*
  integrand (integer effective exponents); **avoid Divonne** (silent bias).
- **Why flattening matters:** it makes the integrand $O(1)$ but *not smooth* — a
  non‑unit effective exponent injects a $\sqrt{\cdot}$ cube‑edge singularity that
  breaks polynomial cubature. The decomposition still pays off even vs Vegas on
  the raw integrand (~340× on a singular case).
- **Batch (7200 coefficient sets):** the big lever is **sample‑sharing + a
  shared coefficient‑independent basis** (~4× MC). Batch CUBA via **chunked
  `ncomp`** (≤512), never per‑kp. Batch winner is dimension‑dependent: Cuhre at
  moderate $n$, **Vegas at $n\gtrsim6$–8**.
- **Wall‑clock (7200 pts, n=8, Apple M2 Pro 12‑core):** ~1 % accuracy in
  **~10–15 s with Vegas vs ~10–15 min with MC** (~30–60× faster).

> **Now implemented.** This study's recommendation ships in the package behind
> `Integrator -> "VEGAS"` (on `GenerateCppMonteCarlo` / `EvaluateTropicalMC` /
> `EvaluateTropicalMCLifted`; default stays `"MC"`, byte‑identical to before).
> MC↔VEGAS parity is demonstrated in `EXAMPLES/test_vegas.wl`. The per‑sector
> Vegas path is in; the `ncomp`‑batched kinematic‑scan lever (§8) remains a
> documented follow‑up. **Lifted sectors** (auxiliary‑variable method) work with
> VEGAS too: for **small coefficients** — the main reason to lift — lifting +
> VEGAS *beats* lifting + MC (e.g. `1e-6`: MC ~13 % with a lying error bar vs
> VEGAS ~3 %); for a **large** coefficient MC edges it (cutoff‑discontinuity
> bias). Both error bars are optimistic in that heavy‑tail regime — cross‑check a
> reference. On the standard continuous sectors VEGAS is unconditionally tighter.
> (cases L1–L3 of `test_vegas.wl`.)

## Files
| file | role |
|---|---|
| `gen_sectors.wl` | runs the pipeline on 15 convergent integrands ($n{=}2$–8); writes flattened sector headers `INTERFILES/sectors_*.hpp` (main() stripped) + exact references + `manifest.csv`. Scale‑invariant fan workaround for thin lattice simplices. |
| `bench.cpp` | integrates the same sectors with MC / Sobol‑QMC / CUBA Vegas, Suave, Divonne, Cuhre; budget sweep; CSV vs exact; QMC canary. |
| `run.sh` | compiles `bench.cpp` per case (Boost + CUBA), runs the sweep → `INTERFILES/results.csv`. |
| `analyze.py` | rates, time‑to‑accuracy, dimension scaling → `INTERFILES/summary.txt`. |
| `plots.py` | `convergence.png`, `dimension_scaling.png`, `canary.png`. |
| `speedup.py` | cost‑to‑target speed‑up factors vs MC. |
| `raw_vs_tropical.wl` | CUBA on the raw vs flattened integrand (is the decomposition worth it?). |
| `gen_kin.wl` | kinematic‑parametrized batch sectors (`KIN_N=n`, default 4): `sectors_kin_n<n>.hpp` + shared‑basis monomial tables. |
| `bench_batch.cpp` | 7200‑coefficient batch test: MC per‑kp / QMC shared / QMC shared‑basis / Vegas per‑kp / Vegas `ncomp` / Cuhre `ncomp`; throughput, accuracy. |
| `REPORT.md`, `RESULTS.md`, `sum/` | full report, condensed report, LaTeX→PDF. |
| `INTERFILES/` | generated headers, binaries, `results.csv`, `batch_results*.txt`, summary, figures (gitignored). |

## Run
```bash
# single-integral study
wolframscript -file gen_sectors.wl && bash run.sh && python3 analyze.py && python3 plots.py
python3 speedup.py ; wolframscript -file raw_vs_tropical.wl
# batch / kinematic-scan study (n=8 example)
KIN_N=8 wolframscript -file gen_kin.wl
g++ -O3 -std=c++17 -I/opt/homebrew/include -DSECTOR_HEADER='"INTERFILES/sectors_kin_n8.hpp"' \
    bench_batch.cpp -L/opt/homebrew/lib -lcuba -lm -o INTERFILES/bench_batch_n8
./INTERFILES/bench_batch_n8 7200 8000
```

## Requirements
g++ (C++17), Boost ≥1.7 (`boost/random/sobol.hpp`), CUBA 4.x (`-lcuba`),
Polymake (fans), Wolfram. All present via Homebrew (`/opt/homebrew`); validated
on Apple M2 Pro (12‑core), CUBA 4.2.2, Boost 1.90, Polymake 4.15.
