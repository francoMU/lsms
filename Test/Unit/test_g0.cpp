//
// Created by F.Moitzi on 23.06.2022.
//

#include <gtest/gtest.h>

#include <iostream>
#include <vector>
#include <complex>

#include "lattice_utils.hpp"
#include "SphericalHarmonics.hpp"
#include "Madelung/Madelung.hpp"
#include "Main/SystemParameters.hpp"
#include "Misc/Coeficients.hpp"
#include "Misc/Indices.hpp"
#include "MultipoleMadelung/MultipoleMadelung.hpp"
#include "MultipoleMadelung/calculateMultipoleMadelung.hpp"
#include "accel_common.hpp"



namespace BareStructureConstants {

using namespace lsms;

TEST(BareStructureConstantsTestSuite, Structure1) {

  /*
   * Structure constants for a simple lattice
   *
   * D^(1)_L = (k,z)
   *
   */

  /**
   * 1. Calculate recirpocal lattice
   */


  double scale = 1;
  double a = 1.0;

  Matrix<double> bravais(3,3);
  Matrix<double> reciprocal_bravais(3,3);

  bravais(0, 0) = 1.0 * a * scale;
  bravais(0, 1) = 0.0;
  bravais(0, 2) = 0.0;

  bravais(1, 0) = 0.0;
  bravais(1, 1) = 1.0 * a * scale;
  bravais(1, 2) = 0.0;

  bravais(2, 0) = 0.0;
  bravais(2, 1) = 0.0;
  bravais(2, 2) = 1.0 * a * scale;


  // This we need to fix
  reciprocal_lattice(bravais, reciprocal_bravais,scale);


  /**
   * 2. Calculate shells of translation vectors
   *
   * - Is it better to do that in a loop?
   *
   */









//
//  int l = 1;
//
//  int lmax = 8;
//
//  std::vector<double> k;
//
//  double omega = 1.0;
//  double nu = 0.1;
//  std::complex<double> z { .001, 3.0};
//
//  //std:complex<double> D1 = - (4*M_PI/omega) * 1l**l * p**-l * exp(z/nu);
//
//  std::size_t nmax = 0;
//
//  SphericalHarmonics spherical_harm(lmax);
//
//  std::vector<std::complex<double>> Ylm;
//
//  for (auto i = 0; i++; i < nmax) {
//
//    std::vector<double> Kn;
//
//    double k_norm = norm(Kn);
//
//    spherical_harm.computeYlm(lmax, Kn, Ylm);
//
//    std::exp(- k_norm * k_norm / nu) / (k_norm * k_norm - z) * Ylm(idx);
//
//
//  }
//




}

}  // namespace multipole_tests