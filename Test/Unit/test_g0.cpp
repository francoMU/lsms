//
// Created by F.Moitzi on 23.06.2022.
//

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <gtest/gtest.h>

#include <complex>
#include <iostream>
#include <vector>

#include "Faddeeva.hh"

#include "Madelung/Madelung.hpp"
#include "Main/SystemParameters.hpp"
#include "Misc/Coeficients.hpp"
#include "Misc/Indices.hpp"
#include "MultipoleMadelung/MultipoleMadelung.hpp"
#include "MultipoleMadelung/calculateMultipoleMadelung.hpp"
#include "MultipoleMadelung/utils.hpp"
#include "SphericalHarmonics.hpp"
#include "accel_common.hpp"
#include "lattice_utils.hpp"

namespace BareStructureConstants {

using namespace lsms;

inline std::complex<double> Iln(std::complex<double> Ilm1,
                                std::complex<double> Ilm2,
                                int l,
                                std::complex<double> kappa,
                                double nu, double R) {

  auto value = (2.0 * l - 1.0) / 2.0 * Ilm1
      - kappa * kappa / 4.0 * Ilm2
      + 0.5 * std::pow(0.5 * std::sqrt(nu), 2 * l - 1)
          * std::exp(-nu * R * R / 4.0 + kappa * kappa / nu);

  return value / (R * R);

}

TEST(BareStructureConstantsTestSuite, D3) {

  using namespace std::complex_literals;

  /**
   * Only L = {0,0} is non-zero
   */

  // Implementation of 7.7.14 p. 298


  constexpr int n_max = 1000;

  double nu = 0.1;
  std::complex<double> kappa{1.0, 0.1};

  std::complex<double> d3{0.0, 0.0};

  for (int n = 0; n < n_max; n++) {
    d3 += std::pow(kappa * kappa / nu, n) / (std::tgamma(n + 1) * (2 * n - 1));
  }

  d3 *= -std::sqrt(nu) / (2.0 * M_PI);

}

TEST(BareStructureConstantsTestSuite, EwaldIntegral) {

  using namespace std::complex_literals;

  std::
  complex<double> Im1{0.0, 0.0};
  std::
  complex<double> I0{0.0, 0.0};
  std::
  complex<double> I1{0.0, 0.0};
  std::
  complex<double> I2{0.0, 0.0};

  /**
   * Parameters
   */
  double nu = 0.1;
  std::complex<double> kappa{1.0, -0.1};
  double R = 0.5;

  double sqr_nu = std::sqrt(nu);

  /**
   * I_0
   */

  {
    std::complex<double> w1 = Faddeeva::w(0.5 * 1i * sqr_nu * R - kappa / sqr_nu);
    std::complex<double> w2 = Faddeeva::w(0.5 * 1i * sqr_nu * R + kappa / sqr_nu);

    I0
        = std::sqrt(M_PI) / (4.0 * R) * std::exp(-nu * R * R / 4.0 + kappa * kappa / nu) * (w1 + w2);

    EXPECT_NEAR(-31.384265676986605, std::real(I0), 1.0e-10);
    EXPECT_NEAR(-182.98965818096264, std::imag(I0), 1.0e-10);
  }

  /**
   * I_-1
   */

  {
    std::complex<double> w1 = Faddeeva::w(0.5 * 1i * sqr_nu * R - kappa / sqr_nu);
    std::complex<double> w2 = Faddeeva::w(0.5 * 1i * sqr_nu * R + kappa / sqr_nu);

    Im1
        = -std::sqrt(M_PI) / (2.0 * 1i * kappa) * std::exp(-nu * R * R / 4.0 + kappa * kappa / nu) * (w1 - w2);

    EXPECT_NEAR(-1406.359666344661, std::real(Im1), 1.0e-10);
    EXPECT_NEAR(-6400.956919043278, std::imag(Im1), 1.0e-10);
  }

  /**
   * I_1
   */

  I1 = Iln(I0, Im1, 1, kappa, nu, R);

  EXPECT_NEAR(3.2795466323414857, std::real(I1), 1.0e-10);
  EXPECT_NEAR(-5.477836427746981, std::imag(I1), 1.0e-10);

  /**
   * I_2
   */

  I2 = Iln(I1, I0, 2, kappa, nu, R);

  EXPECT_NEAR(22.184650062869274, std::real(I2), 1.0e-10);
  EXPECT_NEAR(-0.3634585189695758, std::imag(I2), 1.0e-10);

}

TEST(BareStructureConstantsTestSuite, Structure1) {

  using namespace std::complex_literals;

  /*
   * Structure constants for a simple lattice
   *
   * D^(1)_L = (k,z)
   *
   */

  double scale = 1;
  double a = 1.0;

  Matrix<double> bravais(3, 3);
  Matrix<double> reciprocal_bravais(3, 3);

  /**
   * Define bravais lattice as row-vectors
   */

  bravais(0, 0) = 1.0 * a;
  bravais(0, 1) = 0.5 * a;
  bravais(0, 2) = 0.25 * a;

  bravais(1, 0) = 0.0;
  bravais(1, 1) = 0.9 * a;
  bravais(1, 2) = 0.1 * a;

  bravais(2, 0) = 0.0;
  bravais(2, 1) = 0.0;
  bravais(2, 2) = 0.8 * a;

  // This we need to fix
  reciprocal_lattice(bravais, reciprocal_bravais);

  // Print lattice
  for (int j = 0; j < 3; j++) {
    fmt::print("{:10.5f} {:10.f} {:10.5f}\n", bravais(j, 0), bravais(j, 1),
               bravais(j, 2));
  }

  // Print reciprocal lattice
  for (int j = 0; j < 3; j++) {
    fmt::print("{:10.5f} {:10.f} {:10.5f}\n", reciprocal_bravais(j, 0),
               reciprocal_bravais(j, 1), reciprocal_bravais(j, 2));
  }

  EXPECT_NEAR(reciprocal_bravais(0, 0), 6.28318531, 1.0e-8);
  EXPECT_NEAR(reciprocal_bravais(1, 1), 6.98131701, 1.0e-8);
  EXPECT_NEAR(reciprocal_bravais(2, 2), 7.85398163, 1.0e-8);

  EXPECT_NEAR(reciprocal_bravais(1, 0), -3.4906585, 1.0e-8);
  EXPECT_NEAR(reciprocal_bravais(2, 0), -1.5271631, 1.0e-8);
  EXPECT_NEAR(reciprocal_bravais(2, 1), -0.87266463, 1.0e-8);

  /**
   * Calculate d1 and d2
   *
   * p. 168
   *
   */

  int lmax = 5;
  int nlm = (lmax + 1) * (lmax + 1);
  std::vector<std::complex<double>> Ylm(nlm);
  std::vector<std::complex<double>> Il(lmax + 1);
  double nu = 0.1;

  double omega = lsms::omega(bravais);

  int rnx = 10;
  int rny = 10;
  int rnz = 10;

  int knx = 10;
  int kny = 10;
  int knz = 10;

  std::vector<double> kn(3);
  std::vector<double> rn(3);

  // k-point
  std::vector<double> k_vec{0.1, 0.1, 0.2};

  // Energy point
  std::complex<double> z{0.1, 2.0};

  std::complex<double> kappa = std::sqrt(z);
  std::complex<double> sqr_nu = std::sqrt(nu);

  std::vector<std::complex<double>> d(nlm);
  std::vector<std::complex<double>> d1(nlm);
  std::vector<std::complex<double>> d2(nlm);
  std::complex<double> d3;

  SphericalHarmonics spherical_harm(lmax);

  /**
   * Reciprocal space summation
   */

  for (int ix = -knx; ix <= knx; ix++) {
    for (int iy = -kny; iy <= kny; iy++) {
      for (int iz = -knz; iz <= knz; iz++) {

        // Create shifted lattice
        kn[0] = reciprocal_bravais(0, 0) * knx + reciprocal_bravais(1, 0) * kny +
            reciprocal_bravais(2, 0) * knz;

        kn[1] = reciprocal_bravais(0, 1) * knx + reciprocal_bravais(1, 1) * kny +
            reciprocal_bravais(2, 1) * knz;

        kn[2] = reciprocal_bravais(0, 2) * knx + reciprocal_bravais(1, 2) * kny +
            reciprocal_bravais(2, 2) * knz;

        // Add k-vector
        kn[0] += k_vec[0];
        kn[1] += k_vec[1];
        kn[2] += k_vec[2];

        double norm_kn = norm(kn.begin(), kn.end());

        spherical_harm.computeYlm(lmax, kn, Ylm);

        for (int l = 0; l <= lmax; l++) {
          std::complex<double> prefac =
              std::pow(norm_kn, l) * std::exp(-norm_kn * norm_kn / nu) / (norm_kn * norm_kn - z);

          for (int m = -l; m <= l; m++) {
            int lm = l * l + l + m;

            d1[lm] += prefac * Ylm[lm];
          }
        }
      }
    }
  }

  for (int l = 0; l <= lmax; l++) {
    for (int m = -l; m <= l; m++) {
      int lm = l * l + l + m;

      d1[lm] *= -(4 * M_PI / omega) * std::pow(1i / kappa, l) * std::exp(z / nu);

    }
  }


  /**
   * Real space summation
   */

  for (int ix = -rnx; ix <= rnx; ix++) {
    for (int iy = -rny; iy <= rny; iy++) {
      for (int iz = -rnz; iz <= rnz; iz++) {

        if (ix == 0 && iy == 0 && iz == 0) {
          continue;
        }

        // Create shifted lattice
        rn[0] = bravais(0, 0) * rnx + bravais(1, 0) * rny +
            bravais(2, 0) * rnz;

        rn[1] = bravais(0, 1) * rnx + bravais(1, 1) * rny +
            bravais(2, 1) * rnz;

        rn[2] = bravais(0, 2) * rnx + bravais(1, 2) * rny +
            bravais(2, 2) * rnz;

        double R = norm(rn.begin(), rn.end());

        spherical_harm.computeYlm(lmax, rn, Ylm);

        std::complex<double> w11 = Faddeeva::w(0.5 * 1i * sqr_nu * R - kappa / sqr_nu);
        std::complex<double> w21 = Faddeeva::w(0.5 * 1i * sqr_nu * R + kappa / sqr_nu);

        std::complex<double>
            I0 = std::sqrt(M_PI) / (4.0 * R) * std::exp(-nu * R * R / 4.0 + kappa * kappa / nu) * (w11 + w21);

        std::complex<double> w12 = Faddeeva::w(0.5 * 1i * sqr_nu * R - kappa / sqr_nu);
        std::complex<double> w22 = Faddeeva::w(0.5 * 1i * sqr_nu * R + kappa / sqr_nu);

        std::complex<double> Im1 =
            -std::sqrt(M_PI) / (2.0 * 1i * kappa) * std::exp(-nu * R * R / 4.0 + kappa * kappa / nu) * (w12 - w22);

        Il[0] = I0;
        Il[1] = Iln(I0, Im1, 1, kappa, nu, R);

        for (int l = 2; l <= lmax; l++) {
          Il[l] = Iln(Il[l - 1], Il[l - 2], l, kappa, nu, R);
        }

        for (int l = 0; l <= lmax; l++) {

          std::complex<double> prefac = std::pow(R, l)
              * std::exp(1i * (k_vec[0] * rn[0] + k_vec[1] * rn[1] + k_vec[2] * rn[2])) * Il[l];

          for (int m = -l; m <= l; m++) {
            int lm = l * l + l + m;

            d2[lm] += prefac * Ylm[lm];

          }

        }

      }
    }
  }

  for (int l = 0; l < lmax; l++) {
    for (int m = -l; m <= l; m++) {
      int lm = l * l + l + m;

      d1[lm] *= -(4 * M_PI / omega) * std::pow(1i / kappa, l) * std::exp(z / nu);

    }
  }

  /**
   *
   */

  constexpr int n_max = 1000;

  for (int n = 0; n < n_max; n++) {
    d3 += std::pow(kappa * kappa / nu, n) / (std::tgamma(n + 1) * (2 * n - 1));
  }

  d3 *= -std::sqrt(nu) / (2.0 * M_PI);


  /**
   * Sum up the structure constants
   */


  for (int l = 0; l <= lmax; l++) {
    for (int m = -l; m <= l; m++) {
      int lm = l * l + l + m;
      d[lm] = d1[lm] + d2[lm];
    }
  }

  d[0] += d3;

}

}  // namespace BareStructureConstants