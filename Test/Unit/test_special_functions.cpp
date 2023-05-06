//
// Created by F.Moitzi on 15.05.2022.
//

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>

#include "Misc/special_functions.hpp"
#include "accel_common.hpp"

extern "C" {

void zsphbesj(int *lmax, std::complex<double> *z,
              std::complex<double> *as,
              std::complex<double> *ac,
              std::complex<double> *x1,
              std::complex<double> *bj,
              std::complex<double> *bscatch);

}

namespace special_functions {

TEST(SpecialFunctionsTests, SphBessel1) {

  int lmax = 5;

  std::vector<std::complex<double>> jl(lmax + 1);
  std::vector<std::complex<double>> jlp(lmax + 1);

  std::complex<double> z{10.0, 2.0};

  lsms::sph_bessel_jl<double>(z, lmax, jl.data(), jlp.data());

  fmt::print("{:50.40f}\n", std::real(jl[0]));

  EXPECT_NEAR(-0.2553223981613775520599791284136009342572587088789620662543539975, std::real(jl[0]), 1.0e-40);
  EXPECT_NEAR(-0.2532550511836289099124826301663632388751174697305568810001608432, std::imag(jl[0]), 1.0e-40);

  EXPECT_NEAR(-0.1883321347952032181125164369378265578983948373608604860560136725, std::real(jl[5]), 3.0e-17);
  EXPECT_NEAR(-0.2100076354765599534043726817839389140508409217077481957003189696, std::imag(jl[5]), 3.0e-19);

  EXPECT_NEAR(0.2361692458827652539151456568849180801266897537544604821814201675, std::real(jl[1]), 3.0e-17);
  EXPECT_NEAR(-0.269868217106676756890224106588811851293717732917153938638684256, std::imag(jl[1]), 3.0e-17);

  EXPECT_NEAR(-0.23616924588276456, std::real(jlp[0]), 1.0e-12);
  EXPECT_NEAR(0.2698682171066755, std::imag(jlp[0]), 1.0e-12);

  std::complex<double> as;
  std::complex<double> ac;
  std::complex<double> x1;

  std::fill(jl.begin(), jl.end(), 0.0);

  zsphbesj(&lmax, &z, &as, &ac, &x1, jl.data(), jlp.data());

  EXPECT_NEAR(-0.2553223981613775520599791284136009342572587088789620662543539975, std::real(jl[0]), 1.0e-40);
  EXPECT_NEAR(-0.2532550511836289099124826301663632388751174697305568810001608432, std::imag(jl[0]), 1.0e-40);

  EXPECT_NEAR(-0.1883321347952032181125164369378265578983948373608604860560136725, std::real(jl[5]), 9.0e-17);
  EXPECT_NEAR(-0.2100076354765599534043726817839389140508409217077481957003189696, std::imag(jl[5]), 9.0e-19);

//
}

TEST(SpecialFunctionsTests, SphBessel2) {

  int lmax = 5;

  std::vector<std::complex<double>> jl(lmax + 1);
  std::vector<std::complex<double>> jlp(lmax + 1);

  std::complex<double> z{1.0e-6, 1.0e-6};

  lsms::sph_bessel_jl<double>(z, lmax, jl.data(), jlp.data());

  fmt::print("{:50.40f}\n", std::real(jl[0]));

  EXPECT_NEAR(-3.848003848004142715523584587186900644963298617096297045342727037e-34, std::real(jl[5]), 3.0e-17);
  EXPECT_NEAR(- 3.848003848003550714931583995438136228863020881055425100055520149e-34, std::imag(jl[5]), 3.0e-19);

  std::complex<double> as;
  std::complex<double> ac;
  std::complex<double> x1;

  std::fill(jl.begin(), jl.end(), 0.0);

  zsphbesj(&lmax, &z, &as, &ac, &x1, jl.data(), jlp.data());

  EXPECT_NEAR(-3.848003848004142715523584587186900644963298617096297045342727037e-34, std::real(jl[5]), 3.0e-17);
  EXPECT_NEAR(- 3.848003848003550714931583995438136228863020881055425100055520149e-34, std::imag(jl[5]), 3.0e-19);

//
}



}  // namespace diff_tests