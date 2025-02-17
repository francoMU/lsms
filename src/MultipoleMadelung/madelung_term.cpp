//
// Created by F.Moitzi on 07.09.2022.
//

#include "madelung_term.hpp"

#include "Coeficients.hpp"
#include "Indices.hpp"
#include "SphericalHarmonics.hpp"

void lsms::dlsum(std::vector<Real> &aij, matrix<Real> &rslat, int nrslat,
                 int ibegin, matrix<Real> &knlat, int nknlat, double omega,
                 int lmax_mad, double eta, std::vector<Complex> &dlm) {
  auto kmax_mad = (lmax_mad + 1) * (lmax_mad + 1);
  std::vector<Complex> Ylm((lmax_mad + 1) * (lmax_mad + 1), Complex(0.0, 0.0));
  std::vector<Real> Plm((lmax_mad + 1) * (lmax_mad + 2) / 2, 0.0);

  lsms::SphericalHarmonics sph(lmax_mad);

  std::vector<Real> vec(3);

  std::fill(dlm.begin(), dlm.end(), std::complex<double>(0, 0));

  auto aij2 = std::inner_product(aij.begin(), aij.end(), aij.begin(), 0.0);

  for (int i = nrslat - 1; i >= ibegin; i--) {
    for (int j = 0; j < 3; j++) {
      vec[j] = aij[j] + rslat(j, i);
    }

    auto vlen = norm(vec.begin(), vec.end());

    // Ylm
    sph.computeYlm(lmax_mad, vec, Ylm);

    // Gamma
    auto gamma_l = gamma_func(vlen / eta, lmax_mad);

    auto vhalf = 0.5 * vlen;

    for (auto kl = kmax_mad - 1; kl > 0; kl--) {
      auto l = AngularMomentumIndices::lofk[kl];
      dlm[kl] = dlm[kl] + gamma_l[l] * Ylm[kl] / std::pow(vhalf, l + 1);
    }
  }

  auto rfac = 4.0 * std::sqrt(M_PI);
  for (int i = 1; i < kmax_mad; ++i) {
    dlm[i] *= rfac;
  }

  rfac = -eta * eta / 4.0;

  std::vector<Complex> ctmp(kmax_mad, 0.0);

  for (int i = nknlat - 1; i > 0; i--) {
    for (int j = 0; j < 3; j++) {
      vec[j] = knlat(j, i);
    }

    auto vlen = norm(vec.begin(), vec.end());
    auto knlatsq = norm_sq(vec.begin(), vec.end());

    // Ylm
    sph.computeYlm(lmax_mad, vec, Ylm);

    auto expfac = std::exp(rfac * knlatsq) / knlatsq;

    Real tfac, sintfac, costfac;

    if (aij2 > 1e-8) {
      tfac =
          -(knlat(0, i) * aij[0] + knlat(1, i) * aij[1] + knlat(2, i) * aij[2]);
      sintfac = sin(tfac);
      costfac = cos(tfac);
    } else {
      tfac = 0.0;
      sintfac = 0.0;
      costfac = 1.0;
    }

    auto cfac = -vlen;
    for (auto l = 1; l <= lmax_mad; l += 2) {
      for (auto m = -l; m <= l; m++) {
        auto kl = (l + 1) * (l + 1) - l + m - 1;
        ctmp[kl] = ctmp[kl] + expfac * cfac * sintfac * Ylm[kl];
      }
      cfac = -cfac * knlatsq;
    }

    cfac = -knlatsq;
    for (auto l = 2; l <= lmax_mad; l += 2) {
      for (auto m = -l; m <= l; m++) {
        auto kl = (l + 1) * (l + 1) - l + m - 1;
        ctmp[kl] = ctmp[kl] + expfac * cfac * costfac * Ylm[kl];
      }
      cfac = -cfac * knlatsq;
    }
  }

  rfac = 16.0 * M_PI * M_PI / omega;

  for (auto i = 1; i < kmax_mad; i++) {
    dlm[i] = dlm[i] + rfac * ctmp[i];
  }
}
