//
// Created by F.Moitzi on 18.10.2021.
//

#include "MultipoleMoments.hpp"

void lsms::create_multipole_moments(Qlm &qlm,
                                    int lmax_qlm,
                                    int lmax,
                                    std::vector<double> &r_mesh,
                                    double ws_rad
                                    ) {

  /*
   * Variables
   */

  double fac;

  std::vector<double> rl_w;


  // Quick return
  if (lmax_qlm == 0) {
    qlm = 0.0;
  }


  // Simples case means no spin for the moment

  for (int l = 0; l <= lmax_qlm; l++) {


    /*
     * (r / w)^{l}
     *  r == r_{R}
     */
    for (int i = 0; i < r_mesh.size(); i++) {

    }

    rl_w % val = comp % r_mesh % r * *l / ws_rad * *l;

    fac = sqrt(4 * M_PI) / (2 * l + 1)


  } // l loop


  //lmax = comp % lmax
  //nlm = (lmax + 1) * *2
  //nz = size(comp % p_waves, 1)
  //nsp = size(comp % p_waves, 2)

  //sp_fac = spin_factor(nsp)


}



