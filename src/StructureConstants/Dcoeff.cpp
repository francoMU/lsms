//
// Created by F.Moitzi on 15.07.2023.
//

#include "Dcoeff.h"

void StructConst::Dcoeff::computeDij(std::complex<double> &d1,
                                     std::complex<double> Z,
                                 ) {


  /*
   * Calculate Dij for every poin
   */


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




}
