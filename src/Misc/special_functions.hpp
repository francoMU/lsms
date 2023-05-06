
#include <complex>

namespace lsms {

template<typename T>
void sph_bessel_jl(std::complex<T> z, int lmax, std::complex<T> *jl,
                   std::complex<T> *jlp) {

  int k;
  std::complex<T> fac1, fac2, fac3, den, c, d, z_inv;

  constexpr T eps = std::numeric_limits<T>::epsilon();
  constexpr T small = std::numeric_limits<T>::min();
  constexpr int max_iter = 2000;

  if (std::fabs(z) < eps) {
    for (k = 0; k <= lmax; k++) {
      jl[k] = 0.0;
      jlp[k] = 0.0;
    }

    jl[0] = 1.0;
    jlp[1] = 1.0;
    return;
  }

  z_inv = 1.0 / z;

  fac1 = z_inv * 2.0 * (T) lmax + z_inv * 3.0;
  fac2 = z_inv * (T) lmax;

  if (std::fabs(fac2) < small) {
    fac2 = small;
  }

  den = 1.0;

  c = fac2;
  d = 0.0;

  for (k = 0; k < max_iter; k++) {
    c = fac1 - 1.0 / c;
    d = fac1 - d;

    if (std::fabs(c) < small) {
      c = small;
    }
    if (std::fabs(d) < small) {
      d = small;
    }

    d = 1.0 / d;
    fac3 = d * c;
    fac2 = fac2 * fac3;

    if (std::real(d) < 0.0) {
      den = -den;
    }

    if (std::fabs(fac3 - 1.0) < eps) {
      break;
    }

    fac1 = fac1 + 2.0 * z_inv;
  }

  jl[lmax] = den;
  jlp[lmax] = fac2 * den;

  // Downward recursion
  c = z_inv * (T) lmax;
  for (k = lmax; k > 0; k--) {
    jl[k - 1] = (c + z_inv) * jl[k] + jlp[k];
    c = c - z_inv;
    jlp[k - 1] = c * jl[k - 1] - jl[k];
  }

  den = jl[0];

  d = std::sin(z) * z_inv / den;
  for (k = 0; k <= lmax; k++) {
    jl[k] = jl[k] * d;
    jlp[k] = jlp[k] * d;
  }

}

}  // namespace lsms
