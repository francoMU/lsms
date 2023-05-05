//
// Created by F.Moitzi on 15.05.2022.
//

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

#include "Misc/special_functions.hpp"
#include "accel_common.hpp"

extern "C" {

void sph_jl(std::complex<double> *z, int *l_max, std::complex<double> *jl, std::complex<double> *jl_p);

}

namespace special_functions {

TEST(SFTests, SpecialFunctions) {

  int lmax = 5;

  std::vector<std::complex<double>> jl(lmax + 1);
  std::vector<std::complex<double>> jlp(lmax + 1);

  std::complex<double> z{10.0, 2.0};

  //sph_jl(&z, &lmax, jl.data(), jlp.data());

  lsms::sph_bessel_jl(z, lmax, jl.data(), jlp.data());

  for (auto j : jl) {
    std::cout << j << std::endl;
  }

  /**
   * l = 0
   * -0.2553223981613775520599791284136009342572587088789620662543539975 - 0.2532550511836289099124826301663632388751174697305568810001608432*I
   *
   * l = 1
   * 0.2361692458827652539151456568849180801266897537544604821814201675 - 0.269868217106676756890224106588811851293717732917153938638684256*I
   *
   * l = 5
   * -0.1883321347952032181125164369378265578983948373608604860560136725 - 0.2100076354765599534043726817839389140508409217077481957003189696*I
   *
   *
   * Derivative:
   *
   * l = 0
   *
   * (-0.23616924588276456+0.2698682171066755j)
   *
   * l = 5
   *
   * (-0.1810419573185277+0.16148534934483577j)
   */

  for (auto j : jlp) {
    std::cout << j << std::endl;
  }
//
}

}  // namespace diff_tests