//
// Created by F.Moitzi on 17.12.2021.
//

#ifndef MADELUNG_UTILS_HPP
#define MADELUNG_UTILS_HPP

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

template <typename Iter_T>
auto norm(Iter_T first, Iter_T last) {
  return std::sqrt(std::inner_product(first, last, first, 0.0));
}

template <typename Iter_T>
auto norm_sq(Iter_T first, Iter_T last) {
  return std::inner_product(first, last, first, 0.0);
}

/**
 *
 *  call calGammaFunc to calculate the integral:
 *
 *         inf       2*l
 *  I(l) = int dx * x   * exp(-x**2)
 *          a
 *
 *  and store it in gamma_l(l)
 *
 *  l = 0, 1, ..., lmax.
 *
 *  using: (2*l+1)*I(l) = 2*I(l+1) - c(l) ,
 *
 *
 *         c(l) = a**(2*l+1) * exp(-a**2)
 *
 *         I(0) = sqrt(pi)/2 * erfc(a)
 *
 *  erfc(z) = 1 - erf(z)
 */
template <class T>
std::vector<T> gamma_func(const T &&a, int &lmax) {
  auto r_factor = 0.5 * std::exp(-a * a);

  std::vector<T> gamma(lmax + 1, 0.0);

  // l = 0
  gamma[0] = std::sqrt(M_PI) * 0.5 * std::erfc(a);

  for (int l = 1; l <= lmax; l++) {
    gamma[l] =
        0.5 * (2 * l - 1) * gamma[l - 1] + std::pow(a, (2 * l - 1)) * r_factor;
  }

  return gamma;
}

#endif  // MADELUNG_UTILS_HPP
