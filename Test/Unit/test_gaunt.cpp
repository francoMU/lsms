//
// Created by F.Moitzi on 18.09.2021.
//

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Misc/Indices.hpp"
#include "accel_common.hpp"
#include "common.hpp"

namespace gaunt_tests {

using namespace lsms;

template <typename T>
bool approx_equal(T x, T y, T epsilon) {
  return fabs(x - y) / max(fabs(x), fabs(y)) <= epsilon;
}

template <typename T>
T relative_diff(T ref, T val) {
  return std::fabs(ref - val) / std::fabs(ref);
}

extern "C" {
void cgaunt_(int *lmax, Real *clm, Real *plmg, Real *tg, Real *wg, Real *cgnt,
             int *lofk, int *mofk, int *iprint, char *istop);
}

TEST(GauntFactorTest, Lmax3) {
  int lmax = 3;

  int iprint = 0;
  char istop[256];

  std::vector<int> lofk, mofk, lofj, mofj, kofj;
  int ndlj, ndlm;

  ndlj = (lmax * 2 + 1) * (lmax * 2 + 1);
  ndlm = ((lmax * 2 + 1) * (lmax * 2 + 2)) / 2;

  mofk.resize(ndlj);
  lofk.resize(ndlj);

  mofj.resize(ndlm);
  lofj.resize(ndlm);
  kofj.resize(ndlm);

  int j = 0;
  int k = 0;

  for (int l = 0; l <= lmax * 2; l++) {
    for (int m = 0; m <= l; m++) {
      lofj[j] = l;
      mofj[j] = m;
      kofj[j] = (l + 1) * (l + 1) - l + m - 1;

      j++;
    }

    for (int m = -l; m <= l; m++) {
      lofk[k] = l;
      mofk[k] = m;
      k++;
    }
  }

  std::vector<Real> tg, wg;
  Matrix<Real> plmg;
  Array3d<double> cgnt;
  std::vector<Real> clm(ndlm, 1.0);

  tg.resize(2 * (2 * lmax + 1));
  wg.resize(2 * (2 * lmax + 1));
  cgnt.resize(lmax + 1, (lmax + 1) * (lmax + 1), (lmax + 1) * (lmax + 1));
  plmg.resize(((2 * lmax + 1) * (2 * lmax + 2)) / 2, 2 * lmax + 1);

  cgaunt_(&lmax, &clm[0], &plmg(0, 0), &tg[0], &wg[0], &cgnt(0, 0, 0), &lofk[0],
          &mofk[0], &iprint, &istop[0]);

  EXPECT_NEAR(0.28209479177387842, cgnt(0, 0, 0), 1e-12);
  EXPECT_NEAR(
      0.1802237515728685748747963944546917104709985754562351930533759359,
      cgnt(1, 6, 6), 1e-12);

  /**
   * l1 = 1
   * l2 = 2
   * l3 = 3
   *
   * m1 = 0
   * m2 = 0
   * m3 = 0
   *
   */

  int l1 = 1;
  int l2 = 2;
  int l3 = 3;

  int m1 = 0;
  int m2 = 0;
  int m3 = 0;

  int k1 = l1 * l1 + l1 + m1;
  int k2 = l2 * l2 + l2 + m2;
  int k3 = l3 * l3 + l3 + m3;

  std::cout << cgnt(0, k2, k3) << std::endl;
  std::cout << cgnt(0, k2, k3) << std::endl;

  EXPECT_NEAR(0.2477666950834760618905, cgnt(0, k2, k3), 1e-12);


}

}  // namespace gaunt_tests