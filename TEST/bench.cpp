// ============================================================================
// TEST/bench.cpp
//
// Benchmark: how best to integrate the *pre-processed* (tropically decomposed
// + flattened) integrand?  Every method below integrates the SAME flattened
// per-sector functions integrand_conv_k(y, params) over [0,1]^dim that the
// tropical pipeline (GenerateCppMonteCarlo) emits.  We then sum sectors and
// compare against the analytic reference baked into the generated header.
//
// Methods compared (single-threaded, for apples-to-apples timing):
//   MC      plain pseudo-random uniform Monte Carlo  (== the shipped pipeline)
//   QMC     randomized quasi-MC: Boost Sobol' (Joe-Kuo) + Cranley-Patterson
//           random shifts (R replicates -> unbiased estimate + honest error)
//   VEGAS   CUBA Vegas    (adaptive importance sampling, Sobol-based, seed 0)
//   SUAVE   CUBA Suave    (subregion-adaptive importance sampling)
//   DIVONNE CUBA Divonne  (stratified sampling / region splitting)
//   CUHRE   CUBA Cuhre    (deterministic globally-adaptive cubature)
//
// The sector header is chosen at compile time via -DSECTOR_HEADER="...".
//
// Build (see run.sh):
//   g++ -O3 -std=c++17 -I/opt/homebrew/include -DSECTOR_HEADER='"..."' \
//       bench.cpp -L/opt/homebrew/lib -lcuba -lm -o bench_xxx
//
// Output: CSV rows on stdout:
//   case,family,dim,nsec,method,budget,neval,wall_s,est_re,est_im,
//   reported_err,abs_err,rel_err
// ============================================================================

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <complex>
#include <random>
#include <vector>
#include <chrono>
#include <string>
#include <cstdlib>

#include <boost/random/sobol.hpp>
extern "C" {
#include <cuba.h>
}

#include SECTOR_HEADER   // defines cx, integrand_table[], integrand_dim[],
                         // N_INTEGRANDS, N_PARAMS, MAX_DIM, CASE_*, REF_RE/IM

static double g_dummy[1] = {0.0};   // params placeholder (cases are param-free)

using clk = std::chrono::steady_clock;
static double secs_since(clk::time_point t0) {
  return std::chrono::duration<double>(clk::now() - t0).count();
}

struct Res { double re, im, err_re, err_im; long long neval; double secs; };

// uint64 -> double in [0,1) using the top 53 bits (exact, one draw per coord)
static inline double u01(uint64_t u) {
  return (double)(u >> 11) * (1.0 / 9007199254740992.0); // 2^53
}

// ---------------------------------------------------------------------------
// 1) Plain pseudo-random Monte Carlo  (faithful to the shipped pipeline:
//    one RNG stream per pass, swept sequentially over sectors, Welford stats)
// ---------------------------------------------------------------------------
static Res run_mc(long long N, uint64_t seed) {
  auto t0 = clk::now();
  double tot_re = 0, tot_im = 0, var_re = 0, var_im = 0;
  long long ne = 0;
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  double y[MAX_DIM];
  for (int s = 0; s < N_INTEGRANDS; ++s) {
    int dim = integrand_dim[s];
    double mean_re = 0, mean_im = 0, M2_re = 0, M2_im = 0;
    for (long long k = 0; k < N; ++k) {
      for (int i = 0; i < dim; ++i) y[i] = U(rng);
      cx v = integrand_table[s](y, g_dummy);
      double dre = v.real() - mean_re; mean_re += dre / (k + 1);
      M2_re += dre * (v.real() - mean_re);
      double dim_ = v.imag() - mean_im; mean_im += dim_ / (k + 1);
      M2_im += dim_ * (v.imag() - mean_im);
    }
    tot_re += mean_re; tot_im += mean_im;
    if (N > 1) {
      var_re += M2_re / ((double)N * (N - 1));
      var_im += M2_im / ((double)N * (N - 1));
    }
    ne += N;
  }
  return {tot_re, tot_im, std::sqrt(var_re), std::sqrt(var_im), ne, secs_since(t0)};
}

// ---------------------------------------------------------------------------
// 2) Randomized QMC: scrambled-by-shift Sobol' (Cranley-Patterson rotation),
//    R independent random shifts -> unbiased estimate with a real error bar.
//    Total evaluations per sector = R * M with R*M ~ N.
// ---------------------------------------------------------------------------
static Res run_qmc(long long N, uint64_t seed) {
  auto t0 = clk::now();
  const int R = 16;
  long long M = N / R; if (M < 1) M = 1;
  double tot_re = 0, tot_im = 0, var_re = 0, var_im = 0;
  long long ne = 0;
  std::mt19937_64 shiftRng(seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  double y[MAX_DIM], shift[MAX_DIM];
  std::vector<double> rep_re(R), rep_im(R);
  for (int s = 0; s < N_INTEGRANDS; ++s) {
    int dim = integrand_dim[s];
    for (int r = 0; r < R; ++r) {
      for (int i = 0; i < dim; ++i) shift[i] = U(shiftRng);
      boost::random::sobol eng(dim);
      double acc_re = 0, acc_im = 0;
      for (long long m = 0; m < M; ++m) {
        for (int i = 0; i < dim; ++i) {
          double x = u01(eng());                 // one Sobol draw per coordinate
          double yy = x + shift[i]; yy -= std::floor(yy);
          y[i] = yy;
        }
        cx v = integrand_table[s](y, g_dummy);
        acc_re += v.real(); acc_im += v.imag();
      }
      rep_re[r] = acc_re / M; rep_im[r] = acc_im / M;
    }
    // mean & sample variance of the R replicate means
    double mr = 0, mi = 0;
    for (int r = 0; r < R; ++r) { mr += rep_re[r]; mi += rep_im[r]; }
    mr /= R; mi /= R;
    double sr = 0, si = 0;
    for (int r = 0; r < R; ++r) {
      sr += (rep_re[r] - mr) * (rep_re[r] - mr);
      si += (rep_im[r] - mi) * (rep_im[r] - mi);
    }
    sr /= (R - 1); si /= (R - 1);          // sample variance of a replicate mean
    tot_re += mr; tot_im += mi;
    var_re += sr / R;                      // variance of the grand mean
    var_im += si / R;
    ne += (long long)R * M;
  }
  return {tot_re, tot_im, std::sqrt(var_re), std::sqrt(var_im), ne, secs_since(t0)};
}

// ---------------------------------------------------------------------------
// 3-6) CUBA routines on the same per-sector flattened integrands
// ---------------------------------------------------------------------------
static IntegrandFunc g_fn = nullptr;
static int g_dim = 0;
static int cubaWrap(const int* ndim, const cubareal xx[], const int* ncomp,
                    cubareal ff[], void* ud) {
  (void)ndim; (void)ncomp; (void)ud;
  double y[MAX_DIM];
  for (int i = 0; i < g_dim; ++i) y[i] = (double)xx[i];
  cx v = g_fn(y, g_dummy);
  ff[0] = v.real(); ff[1] = v.imag();
  return 0;
}

enum CubaKind { K_VEGAS, K_SUAVE, K_DIVONNE, K_CUHRE };

static void cuba_one(CubaKind kind, int ndim, long long maxeval, double epsrel,
                     double integ[2], double err[2], int* neval) {
  int fail = 0, nregions = 0;
  cubareal integral[2], error[2], prob[2];
  const double epsabs = 1e-300;
  switch (kind) {
    case K_VEGAS:
      Vegas(ndim, 2, cubaWrap, nullptr, 1, epsrel, epsabs, 0, 0,
            0, (int)maxeval, 1000, 500, 1000, 0, nullptr, nullptr,
            neval, &fail, integral, error, prob);
      break;
    case K_SUAVE:
      Suave(ndim, 2, cubaWrap, nullptr, 1, epsrel, epsabs, 0, 0,
            0, (int)maxeval, 1000, 2, 25.0, nullptr, nullptr,
            &nregions, neval, &fail, integral, error, prob);
      break;
    case K_DIVONNE:
      Divonne(ndim, 2, cubaWrap, nullptr, 1, epsrel, epsabs, 0, 0,
              0, (int)maxeval, 47, 1, 1, 5, 0.0, 10.0, 0.25,
              0, ndim, nullptr, 0, nullptr, nullptr, nullptr,
              &nregions, neval, &fail, integral, error, prob);
      break;
    case K_CUHRE:
      Cuhre(ndim, 2, cubaWrap, nullptr, 1, epsrel, epsabs, 0,
            0, (int)maxeval, 0, nullptr, nullptr,
            &nregions, neval, &fail, integral, error, prob);
      break;
  }
  integ[0] = integral[0]; integ[1] = integral[1];
  err[0] = error[0]; err[1] = error[1];
}

static Res run_cuba(CubaKind kind, long long maxeval, double epsrel) {
  auto t0 = clk::now();
  double tot_re = 0, tot_im = 0, var_re = 0, var_im = 0;
  long long ne = 0;
  for (int s = 0; s < N_INTEGRANDS; ++s) {
    g_fn = integrand_table[s];
    g_dim = integrand_dim[s];
    double integ[2], err[2]; int nev = 0;
    cuba_one(kind, g_dim, maxeval, epsrel, integ, err, &nev);
    tot_re += integ[0]; tot_im += integ[1];
    var_re += err[0] * err[0]; var_im += err[1] * err[1];
    ne += nev;
  }
  return {tot_re, tot_im, std::sqrt(var_re), std::sqrt(var_im), ne, secs_since(t0)};
}

// ---------------------------------------------------------------------------
// QMC canary: integrate a known smooth function exp(sum y_i) over [0,1]^dim
// (exact = (e-1)^dim) with both MC and QMC, to confirm the Sobol generator
// actually delivers a steeper-than-1/sqrt(N) convergence at this dimension.
// ---------------------------------------------------------------------------
static inline double canary_f(const double* y, int dim) {
  double s = 0; for (int i = 0; i < dim; ++i) s += y[i];
  return std::exp(s);
}
static void run_canary(int dim, const std::vector<long long>& budgets) {
  double ref = std::pow(std::exp(1.0) - 1.0, dim);
  std::mt19937_64 rng(12345), shiftRng(67890);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  double y[MAX_DIM], shift[MAX_DIM];
  for (long long N : budgets) {
    // MC
    {
      auto t0 = clk::now();
      double mean = 0;
      for (long long k = 0; k < N; ++k) {
        for (int i = 0; i < dim; ++i) y[i] = U(rng);
        double v = canary_f(y, dim);
        mean += (v - mean) / (k + 1);
      }
      printf("canary,smooth,%d,1,MC,%lld,%lld,%.4f,%.10g,0,0,%.3e,%.3e\n",
             dim, N, N, secs_since(t0), mean, std::fabs(mean - ref),
             std::fabs(mean - ref) / std::fabs(ref));
    }
    // QMC (single Sobol pass with one random shift)
    {
      auto t0 = clk::now();
      for (int i = 0; i < dim; ++i) shift[i] = U(shiftRng);
      boost::random::sobol eng(dim);
      double acc = 0;
      for (long long m = 0; m < N; ++m) {
        for (int i = 0; i < dim; ++i) {
          double x = u01(eng()); double yy = x + shift[i]; yy -= std::floor(yy);
          y[i] = yy;
        }
        acc += canary_f(y, dim);
      }
      double est = acc / N;
      printf("canary,smooth,%d,1,QMC,%lld,%lld,%.4f,%.10g,0,0,%.3e,%.3e\n",
             dim, N, N, secs_since(t0), est, std::fabs(est - ref),
             std::fabs(est - ref) / std::fabs(ref));
    }
  }
}

// ---------------------------------------------------------------------------
static void emit(const char* method, long long budget, const Res& r) {
  cx ref(REF_RE, REF_IM);
  cx est(r.re, r.im);
  double abs_err = std::abs(est - ref);
  double rel_err = abs_err / std::abs(ref);
  double rep_err = std::sqrt(r.err_re * r.err_re + r.err_im * r.err_im);
  printf("%s,%s,%d,%d,%s,%lld,%lld,%.4f,%.12g,%.12g,%.3e,%.3e,%.3e\n",
         CASE_NAME, CASE_FAMILY, CASE_DIM, CASE_NSEC, method, budget,
         r.neval, r.secs, r.re, r.im, rep_err, abs_err, rel_err);
  fflush(stdout);
}

int main(int argc, char** argv) {
  const int zero = 0;
  cubacores(&zero, &zero);   // single process: deterministic, macOS-safe

  // budgets (per-sector samples / maxeval) from argv, else a default sweep
  std::vector<long long> budgets;
  for (int i = 1; i < argc; ++i) budgets.push_back(atoll(argv[i]));
  if (budgets.empty())
    budgets = {1000, 3162, 10000, 31623, 100000, 316228, 1000000};

  // header (only the first binary's; run.sh dedups)
  printf("case,family,dim,nsec,method,budget,neval,wall_s,"
         "est_re,est_im,reported_err,abs_err,rel_err\n");

  run_canary(CASE_DIM, budgets);

  for (long long b : budgets) {
    emit("MC",      b, run_mc(b, 42));
    emit("QMC",     b, run_qmc(b, 42));
    emit("VEGAS",   b, run_cuba(K_VEGAS,   b, 1e-12));
    emit("SUAVE",   b, run_cuba(K_SUAVE,   b, 1e-12));
    emit("DIVONNE", b, run_cuba(K_DIVONNE, b, 1e-12));
    emit("CUHRE",   b, run_cuba(K_CUHRE,   b, 1e-12));
  }
  return 0;
}
