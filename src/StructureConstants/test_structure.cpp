//
// Created by F.Moitzi on 15.05.2022.
//

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>


TEST(StructureTests, LLLMatrix) {
  constexpr auto N = 1000;
  std::vector<double> fun(N);
  std::vector<double> ref(N);
  std::vector<double> x(N);

  for (auto ir{0}; ir < N; ir++) {
    x[ir] = ir;
    fun[ir] = x[ir] * x[ir] + 1.0;
    ref[ir] = 2.0 * x[ir];
  }

  auto res = lsms::derivative<double>(fun, N);

  for (auto ir{0}; ir < N; ir++) {
    EXPECT_NEAR(ref[ir], res[ir], 10e-6);
  }
}

