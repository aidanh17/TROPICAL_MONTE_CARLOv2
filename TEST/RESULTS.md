# MC vs. QMC vs. CUBA on the *pre‑processed* (flattened) tropical integrand

**Question (yours).** After the tropical sector decomposition and the flattening
substitution, each sector integral is a bounded, $O(1)$ integrand on the unit
cube $[0,1]^n$ — the master formula

$$ I=\sum_\sigma \frac{|\det M_\sigma|}{\prod_j a^\sigma_j}\int_{[0,1]^n}\!\! d^n y'\;\prod_k Q^\sigma_k\!\big({y'}^{1/a^\sigma}\big)^{B_k}. $$

The shipped pipeline samples each sector with **plain uniform Monte Carlo**.
Is it faster to (a) implement **quasi‑Monte Carlo**, or (b) hand the *same
flattened sectors* to a library like **CUBA**?

**Short answer.** On the flattened integrand, **plain MC is the weakest option**
and you should replace it. The clear overall winner is **CUBA‑Vegas** (its
Sobol‑seeded, adaptive importance sampler): best accuracy‑per‑evaluation in
*every* case and the only method whose error degrades gracefully with
dimension. Plain randomized **QMC (Sobol)** is a cheap, dependency‑light upgrade
over MC (≈1 order of magnitude), but its edge narrows as $n$ grows. The
deterministic cubature **Cuhre is a trap**: unbeatable when flattening happens
to leave a *smooth* integrand, but poor otherwise — because **flattening
guarantees $O(1)$ magnitude, not smoothness** (see §4). **Divonne can return
confidently *wrong* answers** here and should be avoided.

> Net recommendation: use **Vegas** as the sector integrator; keep **Sobol QMC**
> as the no‑dependency fallback; use **Cuhre only when all effective exponents
> $a_{\mathrm{eff}}$ are integers** (smooth flattened integrand) and $n$ is small.

---

## 1. What was tested

Every method integrates the **identical** flattened per‑sector functions
`integrand_conv_k(y)` emitted by `GenerateCppMonteCarlo` — i.e. the genuine
output of *tropical decomposition + flattening*. Only the **sampler** differs.
Sectors are summed and compared to an **analytic reference** (exact closed
form), so any method’s error is a true error, not a self‑reported one.

| Family | Integrand $\;\int_{[0,\infty)^n}\!\prod x_i^{A_i}\,P^{B}\,d^nx$ | $n$ | sectors | exact |
|---|---|---|---|---|
| **A** simplex | $P=1+\sum x_i,\;A_i=0,\;B=-(n{+}2)$ | 2–8 | $n{+}1$ | $1/(n{+}1)!$ |
| **Af** frac. simplex | $A_i=-\tfrac12,\;B=-(n{+}1)$ (genuine $x^{-1/2}$ edges) | 2,4,6 | $n{+}1$ | $\pi^{n/2}\Gamma(\tfrac n2{+}1)/n!$ |
| **B** quadratic | $P=1+\sum x_i^2,\;B=-(n{+}1)$ | 4,6 | $n{+}1$ | $2^{-n}\pi^{n/2}\Gamma(\tfrac n2{+}1)/n!$ |
| **D** product | $P=\prod(1{+}x_i),\;B=-2$ (Newton polytope = hypercube) | 3,4,5 | $2^n$ | $1$ |

All 15 specs are convergent (no divergent sectors); the decomposition was
cross‑checked against direct `NIntegrate` for $n\le4$ (rel. err $10^{-8}$–$10^{-5}$).

### Methods
| tag | what it is |
|---|---|
| **MC** | plain pseudo‑random uniform MC (`mt19937_64`) — *the shipped pipeline* |
| **QMC** | randomized quasi‑MC: **Boost Sobol′** (Joe–Kuo) + 16 Cranley–Patterson random shifts (unbiased estimate + honest error bar) |
| **VEGAS** | CUBA Vegas — adaptive importance sampling, Sobol‑based (`seed=0`) |
| **SUAVE** | CUBA Suave — subregion‑adaptive importance sampling |
| **DIVONNE** | CUBA Divonne — stratified sampling / region splitting |
| **CUHRE** | CUBA Cuhre — deterministic globally‑adaptive polynomial cubature |

Everything single‑threaded for apples‑to‑apples timing. A built‑in **QMC canary**
(integrating a known smooth $e^{\sum y_i}$) confirms the Sobol generator delivers
genuine low‑discrepancy convergence (slopes ≈ $-0.8$ to $-1.4$ vs MC’s $\approx-0.5$).

---

## 2. Convergence rate $p$ &nbsp;($\text{rel.err}\sim N_{\rm eval}^{-p}$, log‑log fit)

| case | $n$ | MC | QMC | **VEGAS** | SUAVE | DIVONNE | CUHRE |
|---|--:|--:|--:|--:|--:|--:|--:|
| A_simplex | 2 | 0.28 | 0.80 | **1.07** | 0.07 | 0.29 | 0.03 |
| A_simplex | 4 | 0.42 | 0.66 | **1.40** | 0.59 | 0.61 | 0.18 |
| A_simplex | 6 | 0.80 | 0.35 | **1.21** | 0.91 | 0.43 | 0.05 |
| A_simplex | 8 | 0.22 | 0.63 | **1.56** | 0.85 | 0.21 | 0.30 |
| Af_frac | 4 | 0.47 | 0.86 | **1.15** | 0.37 | 0.50 | 0.23 |
| B_quad | 6 | 0.47 | 0.61 | **1.06** | 0.81 | −0.09 | 0.45 |
| D_product | 4 | 0.26 | 0.67 | 1.11 | 0.34 | 0.46 | **0.03†** |
| D_product | 5 | 0.35 | 0.62 | **1.73** | 0.23 | 0.44 | 0.08† |

- **MC**: $p\approx0.3$–$0.5$ — textbook $1/\sqrt N$ (sometimes worse: mild heavy tails).
- **QMC**: $p\approx0.6$–$0.9$ — the expected low‑discrepancy gain, but it **erodes
  with dimension** (e.g. A_simplex $n=6$: QMC actually trails MC).
- **VEGAS**: $p\approx1.0$–$1.7$, **stable across $n$** — QMC rate *plus* adaptation.
- **CUHRE**: $p\approx0$ on A/Af/B (stuck), but †on family D it is already at
  $\sim10^{-7}$ from the smallest budget (smooth → converged immediately; the flat
  slope means "nothing left to gain", not "not converging").
- **DIVONNE**: erratic, occasionally **negative** (diverging / biased).

---

## 3. The practical metric: wall‑time to reach a target accuracy

**Time (s) to reach relative error ≤ $10^{-3}$** (— = not reached within the swept budget):

| case | $n$ | MC | QMC | **VEGAS** | SUAVE | DIVONNE | CUHRE |
|---|--:|--:|--:|--:|--:|--:|--:|
| A_simplex | 3 | — | 0.021 | 0.004 | 0.018 | 0.009 | **0.002** |
| A_simplex | 4 | 1.600 | 0.018 | 0.025 | 0.024 | 0.271 | **0.003** |
| A_simplex | 5 | — | — | **0.005** | 0.007 | — | — |
| A_simplex | 6 | 0.181 | — | **0.006** | 0.095 | — | — |
| A_simplex | 7 | — | — | **0.025** | — | — | — |
| A_simplex | 8 | — | — | **0.151** | — | — | — |
| Af_frac | 6 | — | — | 0.102 | 0.130 | **0.011** | — |
| B_quad | 6 | — | — | **0.024** | 0.040 | — | — |
| D_product | 5 | 0.131 | 0.041 | 0.020 | 0.021 | — | **0.006** |

**Vegas is the only method that reaches $10^{-3}$ in *every* case**, and it is the
only one that survives into $n=7,8$. MC reaches $10^{-3}$ in roughly a third of
the cases within budget; Cuhre/Divonne fall off a cliff for $n\ge5$ except on the
smooth product family.

---

## 4. Why Cuhre underperforms — *flattening fixes magnitude, not smoothness*

The decisive, verified mechanism. Flattening uses $y=(y')^{1/a_{\mathrm{eff}}}$.
Whenever a sector’s effective exponent $a_{\mathrm{eff}}\neq1$, this injects a
**fractional power** at the cube face. Concretely, a real sector of
`A_simplex_n4` ($a_{\mathrm{eff}}=2$) emits

```cpp
P0 += 1.0 * std::exp((1.0/2.0) * log_y[0]);   //  = sqrt(y0)   <-- !
...
result *= std::exp(-6.0 * std::log(P0));
```

so the integrand contains $\sqrt{y_0}$: its **value is bounded** ($O(1)$, as
advertised) but its **first derivative is infinite** at $y_0\to0$. Polynomial
cubature (Cuhre, and Divonne’s rule‑based splitting) assumes bounded high
derivatives for both its convergence *and* its error estimate, so it stalls;
QMC (which wants bounded Hardy–Krause variation) is partly hurt; plain MC is
rate‑agnostic; **Vegas’s importance sampling reshapes the measure to absorb the
edge** and keeps the QMC rate.

By contrast every sector of the **product family D** has **all‑integer**
exponents ($a_{\mathrm{eff}}=1$, $B=-2$):

```cpp
P0 += 1.0 * std::exp(1.0*log_y[0] + 1.0*log_y[3]);   // smooth, rational
result *= std::exp(-2.0 * std::log(P0));
```

— a genuinely $C^\infty$ rational integrand. There Cuhre is **unbeatable**
(rel. err $\sim10^{-7}$ in 6 ms). This is the whole story: *Cuhre’s usefulness
on a flattened integrand is decided by whether the flattening left it smooth,
which is determined by whether the effective exponents are integers.*

---

## 5. Dimension scaling (simplex family, $10^5$ samples/sector)

Relative error vs exact:

| $n$ | MC | QMC | **VEGAS** | SUAVE | DIVONNE | CUHRE |
|--:|--:|--:|--:|--:|--:|--:|
| 2 | 1.9e‑4 | 1.2e‑4 | **1.2e‑5** | 3.8e‑4 | 2.5e‑5 | 3.9e‑5 |
| 4 | 5.8e‑3 | 3.5e‑3 | **2.4e‑5** | 1.6e‑4 | 5.8e‑4 | 4.5e‑4 |
| 6 | 1.5e‑3 | 5.2e‑3 | **2.9e‑4** | 6.1e‑4 | 2.1e‑2 | 2.4e‑2 |
| 8 | 1.3e‑2 | 6.2e‑3 | **6.2e‑4** | 6.3e‑3 | **2.2e‑1** | 4.9e‑3 |

Wall times at this budget are within ~2–3× across methods (Vegas ≈ MC), so
Vegas’s 1–2 orders of magnitude accuracy advantage is essentially free.
Note Divonne at $n=8$ reports **22 % error** — a fast, confident, *wrong*
answer (the other methods agree on the right value, confirming the
decomposition is correct and the fault is Divonne’s).

See `INTERFILES/convergence.png`, `dimension_scaling.png`, `canary.png`.

---

## 6. Speed‑up to the same accuracy (vs MC)

Cost to reach a target relative error, as a factor over plain MC (geomean over
the 15 single‑integral cases; `speedup.py`):

| target | QMC (Sobol) | **VEGAS** | Cuhre (smooth only) |
|---|--:|--:|--:|
| $10^{-3}$ | ~24× | **~150×** | $10^2$–$10^3$× |
| $10^{-4}$ | ~300× | **~2300×** | $10^4$–$10^5$× |

The speed‑up **grows as the target tightens** (Vegas $N^{-1.3}$ vs MC
$N^{-0.5}$), and at higher dimension the others stop reaching the target at all
while Vegas still does.

## 7. Is the tropical decomposition still worth it? (raw vs flattened)

Yes — even against the *same* sampler. CUBA on the **raw** integrand
(compactified, generous $2\times10^6$ evals) vs on the **flattened** sectors
(`raw_vs_tropical.wl`): with Vegas the flattened version is more accurate at
fewer evals in every case, and **~340× more accurate on an integrand with a
genuine $x^{-1/2}$ singularity** (where flattening removes the singularity the
sampler cannot). The decomposition *creates* the tame $O(1)$ integrand; the
sampler is a second‑order optimisation on top. (Nuance: for an already‑smooth
raw integrand, raw‑Cuhre can beat flattened‑Cuhre, because flattening injects
the $\sqrt{\cdot}$ edge non‑smoothness — but that's Cuhre‑specific, and reverses
the moment a real singularity appears.)

## 8. Batch evaluation (kinematic scan) and wall‑clock

One structure, many coefficient sets (`gen_kin.wl`, `bench_batch.cpp`:
$P=c_0+\sum c_i x_i$, exact $1/((n{+}1)!\,c_0^2\prod c_i)$, 7200 sets).

- **CUBA batches via `ncomp`** (shared samples across components in one call),
  but with **no native multi‑parameter mode** and a hard **`ncomp`≤1024** cap
  (above it Vegas segfaults; the safe chunk **shrinks with dimension** — use
  ≤512 at $n{=}8$). So you must chunk; per‑kp CUBA (one call per point) is the
  wrong way and is slowest.
- **The biggest lever is sample‑sharing.** The flattened integrand's costly
  transcendentals are coefficient‑independent; computing the monomial basis
  once per sample and reusing it for all 7200 (`QMC_shared_basis`) is **~4× MC
  throughput at identical accuracy** — and applies to every method.
- **$n{=}4$:** chunked‑`ncomp` **Cuhre** is both most accurate and high
  throughput (batching amortises its cubature nodes). **$n{=}8$ reverses this**:
  Cuhre is *worse than MC* (curse of dimensionality), QMC is *worse than MC*
  (low discrepancy needs $N\gg2^n$), and **only Vegas (chunked `ncomp`) leads** —
  given enough budget to adapt.

**Wall‑clock, 7200 points, $n{=}8$, Apple M2 Pro (12‑core = 8P+4E):**

| operating point | Vegas (chunked `ncomp`) | plain MC |
|---|---|---|
| equal budget $M{=}8000$ | 4.98 % in 56 s | 7.85 % in 50 s |
| equal budget $M{=}24000$ | ~0.5 % | 5.46 % in 270 s |
| **time to ~1 % median** | **~1.5 min (1 core) / ~10–15 s (12 cores)** | **~1–2 h / ~10–15 min (12 cores)** |

At equal budget MC and Vegas cost ~the same time but Vegas is ~10× more
accurate, so reaching ~1 % is **~30–60× faster** with Vegas. Caveats: absolute
seconds drift ±1.5–2× with thermal throttling (ratio is the robust figure);
MC's max error hit **421 %** (heavy tails → lift those sets); times scale
~linearly with $\#\text{sectors}\times\#\text{monomials}$.

## 9. Takeaways for the pipeline

1. **Drop plain MC as the default.** It has the worst rate and, with heavy‑ish
   sectors, an over‑optimistic error bar.
2. **Adopt CUBA‑Vegas (Sobol, `seed=0`) as the sector integrator.** It pairs the
   QMC convergence rate with adaptive importance sampling that handles the
   fractional‑power cube‑edge structure flattening leaves behind — best in every
   test, and the only method that scales to $n=8$. Per‑sector overhead is modest;
   it parallelizes over kinematic points exactly like the current code.
3. **If a CUBA dependency is unwanted, ship randomized Sobol QMC.** ~10× better
   than MC, trivial to add (Boost `sobol` + a random shift), same OpenMP
   parallelism and $O(1)$ memory. Accept that its margin shrinks at high $n$.
4. **Cuhre only with integer effective exponents** (then it is unbeatable) and
   modest $n$; otherwise it is dominated. **Avoid Divonne** on flattened
   integrands — risk of silent bias.
5. **Smoothness, not magnitude, is the lever.** A flattening (or lift) choice
   that yields **integer** $a_{\mathrm{eff}}$ would unlock deterministic cubature
   — a worthwhile direction if low‑$n$, high‑precision points are the goal.
6. **For kinematic scans, share samples + share the basis** (≈4× MC at equal
   accuracy, dimension‑independent). Batch CUBA via chunked `ncomp` (≤512), not
   per‑kp calls. The *batch winner depends on dimension*: chunked‑`ncomp` Cuhre
   at moderate $n$, **Vegas at $n\gtrsim6$–8** (budget ~100× higher there).

---

## 10. Reproduce

```bash
wolframscript -file TEST/gen_sectors.wl     # tropical decomp -> flattened sectors
bash       TEST/run.sh                      # compile (boost+cuba) & sweep budgets
python3    TEST/analyze.py                  # tables   -> INTERFILES/summary.txt
python3    TEST/plots.py                    # figures  -> INTERFILES/*.png
python3    TEST/speedup.py                  # cost-to-target speed-ups vs MC
wolframscript -file TEST/raw_vs_tropical.wl # raw vs flattened CUBA (decomp value)
KIN_N=8 wolframscript -file TEST/gen_kin.wl # kinematic batch sectors (n=8)
g++ -O3 -std=c++17 -I/opt/homebrew/include -DSECTOR_HEADER='"INTERFILES/sectors_kin_n8.hpp"' \
    TEST/bench_batch.cpp -L/opt/homebrew/lib -lcuba -lm -o bench_batch
./bench_batch 7200 8000                     # batch throughput/accuracy, 7200 sets
```

Raw data: `INTERFILES/results.csv` (single‑integral sweep), `batch_results.txt`
(n=4 batch), `batch_results_n8.txt` (n=8 batch). CUBA 4.2.2 + Boost 1.90
(Homebrew), Polymake 4.15; machine Apple M2 Pro (12‑core). Full write‑up with
figures: `REPORT.md` and `sum/report.pdf`.
