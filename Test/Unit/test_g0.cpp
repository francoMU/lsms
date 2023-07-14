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

  int nx = 10;
  int ny = 10;
  int nz = 10;

  std::vector<double> kn(3);
  std::vector<double> rn(3);

  std::vector<double> k_vec{0.1, 0.1, 0.2};
  double nu = 0.1;
  std::complex<double> z{0.1, 2.0};

  double omega = lsms::omega(bravais);

  std::vector<std::complex<double>> d1(nlm);
  std::vector<std::complex<double>> d2(nlm);
  std::vector<std::complex<double>> d3(nlm);

  SphericalHarmonics spherical_harm(lmax);

  for (int ix = -nx; ix <= nx; ix++) {
    for (int iy = -ny; iy <= ny; iy++) {
      for (int iz = -nz; iz <= nz; iz++) {
        // Create shifted lattice
        kn[0] = reciprocal_bravais(0, 0) * nx + reciprocal_bravais(1, 0) * ny +
            reciprocal_bravais(2, 0) * nz;

        kn[1] = reciprocal_bravais(0, 1) * nx + reciprocal_bravais(1, 1) * ny +
            reciprocal_bravais(2, 1) * nz;

        kn[2] = reciprocal_bravais(0, 2) * nx + reciprocal_bravais(1, 2) * ny +
            reciprocal_bravais(2, 2) * nz;

        // Add k-vector
        kn[0] += k_vec[0];
        kn[1] += k_vec[1];
        kn[2] += k_vec[2];

        double norm_kn = norm(kn.begin(), kn.end());

        spherical_harm.computeYlm(lmax, kn, Ylm);

        for (int l = 0; l < lmax; l++) {
          double prefac =
              std::pow(norm_kn, l) * std::exp(-norm_kn * norm_kn / nu);

          for (int m = -l; m <= l; m++) {
            int lm = l * l + l + m;

            d1[lm] += prefac * Ylm[lm];
          }
        }
      }
    }
  }

  // Prefactor still need to be added

  /**
   * Real space summation
   */

  for (int ix = -nx; ix <= nx; ix++) {
    for (int iy = -ny; iy <= ny; iy++) {
      for (int iz = -nz; iz <= nz; iz++) {

        // Create shifted lattice
        rn[0] = bravais(0, 0) * nx + bravais(1, 0) * ny +
            bravais(2, 0) * nz;

        rn[1] = bravais(0, 1) * nx + bravais(1, 1) * ny +
            bravais(2, 1) * nz;

        rn[2] = bravais(0, 2) * nx + bravais(1, 2) * ny +
            bravais(2, 2) * nz;

        double norm_rn = norm(rn.begin(), rn.end());

        spherical_harm.computeYlm(lmax, rn, Ylm);

        for (int l = 0; l < lmax; l++) {





//          double prefac = std::pow(norm_rn, l)
//              * std::exp(1i * (k_vec[0] * rn[0] + k_vec[1] * rn[1] + k_vec[2] * rn[2]))
//
//          for (int m = -l; m <= l; m++) {
//            int lm = l * l + l + m;
//
//            d2[lm] += prefac * Ylm[lm];
//
//
//          }

        }

      }
    }
  }

  // Prefactor still need to be added
}

}  // namespace BareStructureConstants