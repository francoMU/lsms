
#include <complex>

namespace lsms {

void sph_bessel_jl(std::complex<double> z, int lmax, std::complex<double> *jl,
                   std::complex<double> *jlp) {
  int l = 0;

  int k;
  std::complex<double> tk, cf1, dcf1, den, c, d, z_1;
  double dummy;

  constexpr double eps = std::numeric_limits<double>::epsilon();
  constexpr double small = std::numeric_limits<double>::min();
  constexpr int max_iter = 2000;

  if (std::fabs(z) < eps) {
    for (l = 0; l <= lmax; l++) {
      jl[l] = 0.0;
      jlp[l] = 0.0;
    }

    jl[0] = 1.0;
    jlp[1] = 1.0;
    return;
  }

  z_1 = 1.0 / z;

  //  // Modified Steed's algorithm
  tk = z_1 * 2.0 * (double) lmax + z_1 * 3.0;
  cf1 = z_1 * (double) lmax;

  if (std::fabs(cf1) < small) {
    cf1 = small;
  }

  den = 1.0;

  c = cf1;
  d = 0.0;

  for (k = 0; k < max_iter; k++) {
    c = tk - 1.0 / c;
    d = tk - d;

    if (std::fabs(c) < small) {
      c = small;
    }
    if (std::fabs(d) < small) {
      d = small;
    }

    d = 1.0 / d;
    dcf1 = d * c;
    cf1 = cf1 * dcf1;

    if (std::real(d) < 0.0) {
      den = -den;
    }

    if (std::fabs(dcf1 - 1.0) < eps) {
      break;
    }

    tk = tk + 2.0 * z_1;
  }

  jl[lmax] = den;
  jlp[lmax] = cf1 * den;

  // Downward recursion
  c = z_1 * (double) lmax;
  for (k = lmax; k > 0; k--) {
    jl[k - 1] = (c + z_1) * jl[k] + jlp[k];
    c = c - z_1;
    jlp[k - 1] = c * jl[k - 1] - jl[k];
  }

  den = jl[0];

  d = std::sin(z) * z_1 / den;
  for (k = 0; k <= lmax; k++) {
    jl[k] = jl[k] * d;
    jlp[k] = jlp[k] * d;
  }

}

}  // namespace lsms
