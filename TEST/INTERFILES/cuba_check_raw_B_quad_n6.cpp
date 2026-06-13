// Auto-generated CUBA cross-check (direct integrand, no tropical
// decomposition).  Compactification x_i = t_i/(1-t_i).
#include <cmath>
#include <complex>
#include <cstdio>
extern "C" {
#include <cuba.h>
}
using cx = std::complex<double>;

static int Integrand(const int *ndim, const cubareal tt[],
                     const int *ncomp, cubareal ff[], void *userdata) {
  (void)ndim; (void)ncomp; (void)userdata;
  double x[6];
  double jac = 1.0;
  for (int i = 0; i < 6; ++i) {
    double t = tt[i];
    if (t < 1e-12) t = 1e-12;
    if (t > 1.0 - 1e-12) t = 1.0 - 1e-12;
    const double u = 1.0 - t;
    x[i] = t / u;
    jac /= (u * u);
  }
  const double P1 = 1. + 1. * std::pow(x[0], 2) + 1. * std::pow(x[1], 2) + 1. * std::pow(x[2], 2) + 1. * std::pow(x[3], 2) + 1. * std::pow(x[4], 2) + 1. * std::pow(x[5], 2);
  const cx logI = cx(-7., 0) * std::log(P1);
  const cx val = std::exp(logI) * jac;
  ff[0] = val.real();
  ff[1] = val.imag();
  return 0;
}

int main() {
  const int zero = 0;
  cubacores(&zero, &zero);   // single process: deterministic, macOS-safe

  const int ndim = 6, ncomp = 2;
  int neval, fail, nregions;
  cubareal integral[2], error[2], prob[2];

  Vegas(ndim, ncomp, Integrand, nullptr, 1,
        1.e-12, 1e-12, 0, 1,
        0, 2000000, 20000, 10000, 1000,
        0, nullptr, nullptr,
        &neval, &fail, integral, error, prob);
  std::printf("VEGAS %.12e %.12e %.3e %.3e %d %d\n",
              integral[0], integral[1], error[0], error[1], neval, fail);

  Cuhre(ndim, ncomp, Integrand, nullptr, 1,
        1.e-12, 1e-12, 0, 0, 2000000,
        0, nullptr, nullptr,
        &nregions, &neval, &fail, integral, error, prob);
  std::printf("CUHRE %.12e %.12e %.3e %.3e %d %d\n",
              integral[0], integral[1], error[0], error[1], neval, fail);

  return 0;
}
