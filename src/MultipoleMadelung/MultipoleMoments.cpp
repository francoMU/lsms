//
// Created by F.Moitzi on 15.04.2023.
//

#include "MultipoleMoments.hpp"

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>

namespace lsms {

void calculateMultipoleMoments(LSMSCommunication &comm,
                               LSMSSystemParameters &lsms,
                               LocalTypeInfo &local,
                               std::vector<std::vector<NonRelativisticSingleScattererSolution>> &solutionNonRel,
                               std::vector<Matrix<Complex>> &tau00_l,
                               std::vector<Complex> &dele1,
                               int l_mom_max) {

  int maxkkrsz = (lsms.maxlmax + 1) * (lsms.maxlmax + 1);
  int max_mmoms = (l_mom_max + 1) * (l_mom_max + 1);

  for (auto i_atom = 0; i_atom < local.num_local; i_atom++) {

    std::vector<std::complex<Real>> rad(local.atom[i_atom].jmt, 0.0);

    for (auto iie = 0; iie < lsms.energyContour.groupSize(); iie++) {

      std::vector<std::complex<Real>> radial(local.atom[i_atom].jmt, 0.0);

      for (auto isp = 0; isp < lsms.n_spin_pola; isp++) {

        NonRelativisticSingleScattererSolution &SSSol = solutionNonRel[iie][i_atom];

        for (auto l = 0; l <= local.atom[i_atom].lmax; l++) {
          for (auto m = -l; m <= l; m++) {

            std::complex<double> tau00 = 0.0;

            auto kl = l * l + l + m;

            if (isp == 0) {
              tau00 = tau00_l[iie](kl + (maxkkrsz * kl), i_atom);
            } else {
              tau00 = tau00_l[iie](kl + (maxkkrsz * kl), i_atom + local.num_local);
            }

            for (auto ir = 0; ir < local.atom[i_atom].jmt; ir++) {

              radial[ir] += SSSol.zlr(ir, l, isp) * tau00 * SSSol.zlr(ir, l, isp) -
                  SSSol.zlr(ir, l, isp) * SSSol.jlr(ir, l, isp);

            }

          }

        }

      }

      for (auto ir = 0; ir < local.atom[i_atom].jmt; ir++) {

        rad[ir] -= 2.0 * std::imag(radial[ir] * dele1[iie])
            * local.atom[i_atom].r_mesh[ir] * local.atom[i_atom].r_mesh[ir] / M_PI;

      }

    }

    auto res = lsms::radialIntegral(rad, local.atom[i_atom].r_mesh, local.atom[i_atom].jmt);


  }

  for (auto i_atom = 0; i_atom < local.num_local; i_atom++) {

    local.atom[i_atom].multi_moms = std::vector<Complex>(max_mmoms, 0.0);

    for (auto l_mom = 0; l_mom <= l_mom_max; l_mom++) {

      for (auto iie = 0; iie < lsms.energyContour.groupSize(); iie++) {

        for (auto isp = 0; isp < lsms.n_spin_pola; isp++) {

          NonRelativisticSingleScattererSolution &SSSol = solutionNonRel[iie][i_atom];

          Array3d<std::complex<double>> zl_zlp(local.atom[i_atom].lmax + 1,
                                               local.atom[i_atom].lmax + 1,
                                               2);

          Matrix<std::complex<double>> zl_jlp(local.atom[i_atom].lmax + 1, 2);


          // Solve integral Z_l * Z_lp
          for (auto l = 0; l <= local.atom[i_atom].lmax; l++) {

            for (auto lp = 0; lp <= local.atom[i_atom].lmax; lp++) {

              std::vector<std::complex<Real>> radial(local.atom[i_atom].jmt);

              for (auto ir = 0; ir < local.atom[i_atom].jmt; ir++) {

                radial[ir] = std::pow(local.atom[i_atom].r_mesh[ir], l_mom + 2) *
                    SSSol.zlr(ir, l, isp) * SSSol.zlr(ir, lp, isp);

              }

              auto res = lsms::radialIntegral(radial, local.atom[i_atom].r_mesh, local.atom[i_atom].jmt);

              zl_zlp(l, lp, isp) = res;

            }

          }

          // Solve integral Z_l * J_l * delta_l,lp
          for (auto l = 0; l <= local.atom[i_atom].lmax; l++) {

            std::vector<std::complex<Real>> radial(local.atom[i_atom].jmt);

            for (auto ir = 0; ir < local.atom[i_atom].jmt; ir++) {

              radial[ir] = std::pow(local.atom[i_atom].r_mesh[ir], l_mom + 2) *
                  SSSol.zlr(ir, l, isp) * SSSol.jlr(ir, l, isp);

            }

            auto res = lsms::radialIntegral(radial, local.atom[i_atom].r_mesh, local.atom[i_atom].jmt);

            zl_jlp(l, isp) = res;

          }

          for (auto m_mom = -l_mom; m_mom <= l_mom; m_mom++) {
            auto il_mom = l_mom * l_mom + l_mom + m_mom;

            for (auto l = 0; l <= local.atom[i_atom].lmax; l++) {
              for (auto m = -l; m <= l; m++) {

                for (auto lp = 0; lp <= local.atom[i_atom].lmax; lp++) {
                  for (auto mp = -lp; mp <= lp; mp++) {

                    auto kl = l * l + l + m;
                    auto klp = lp * lp + lp + mp;

                    std::complex<double> tau00 = 0.0;

                    if (isp == 0) {
                      tau00 = tau00_l[iie](kl + (maxkkrsz * klp), i_atom);
                    } else {
                      tau00 = tau00_l[iie](kl + (maxkkrsz * klp), i_atom + local.num_local);
                    }

                    auto gf_LLp = zl_zlp(l, lp, isp) * tau00 ;

                    if (kl == klp) {
                      gf_LLp -= zl_jlp(l, isp);
                    }

                    if (((l_mom + l + lp) % 2 == 0) && m - mp == m_mom) {

//                      if(iie == 10 && m_mom == 0 && l_mom == 1) {
//
//                        fmt::print("{:2d} {:2d} {:2d} {:2d}: {:10.6f} {:10.6e} {:10.6e}\n",
//                                   l,
//                                   lp,
//                                   m,
//                                   mp,
//                                   GauntCoeficients::cgnt(l_mom / 2, kl, klp),
//                                   std::real(tau00), std::imag(tau00));
//
//                      }
                      local.atom[i_atom].multi_moms[il_mom] -= 2.0 * GauntCoeficients::cgnt(l_mom / 2, kl, klp)
                          * std::imag(gf_LLp * dele1[iie]) * std::sqrt(4 * M_PI) / ((2.0 * l_mom + 1) * M_PI);


                    }

                  }

                }
              }

            }

          }

        }

      }
    }

  }

}



void printMultipoleMoments(LSMSCommunication &comm,
                           LSMSSystemParameters &lsms,
                           LocalTypeInfo &local,
                           int l_mom_max
) {

  for (auto i_atom = 0; i_atom < local.num_local; i_atom++) {

    auto &moms = local.atom[i_atom].multi_moms;

    for (auto l = 0; l <= l_mom_max; l++) {

      for (auto m = -l; m <= l; m++) {

        int k_mom = l * l + l + m;

        if(std::fabs(std::real(moms[k_mom])) > 1.0e-8) {
          fmt::print(" {:2d} {:2d}: {:20.8f}\n", l, m, std::real(moms[k_mom]));
        }


      }

    }

  }

}

} // lsms
