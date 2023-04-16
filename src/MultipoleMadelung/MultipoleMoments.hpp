//
// Created by F.Moitzi on 15.04.2023.
//

#ifndef LSMS_SRC_MULTIPOLEMADELUNG_MULTIPOLEMOMENTS_HPP_
#define LSMS_SRC_MULTIPOLEMADELUNG_MULTIPOLEMOMENTS_HPP_

#include "SingleSiteScattering.hpp"
#include "LSMSCommunication.hpp"
#include "integrator.hpp"
#include "Coeficients.hpp"

namespace lsms {

void calculateMultipoleMoments(LSMSCommunication &comm,
                               LSMSSystemParameters &lsms,
                               LocalTypeInfo &local,
                               std::vector<Complex> &dele1) {


  // How to access the  wave functions
  int maxkkrsz = (lsms.maxlmax + 1) * (lsms.maxlmax + 1);
  int maxkkrsz_ns = lsms.n_spin_cant * maxkkrsz;
  int l_mom_max = 3;

  int max_mmoms = (l_mom_max + 1) * (l_mom_max + 1);

  // solution{Non}Rel[energy point][local atom idx]
  std::vector<std::vector<NonRelativisticSingleScattererSolution>> solutionNonRel;

  Matrix<Complex> tau00_l(
      maxkkrsz_ns * maxkkrsz_ns,
      local.num_local * lsms.n_spin_pola /
          lsms.n_spin_cant);

  for (auto i_atom = 0; i_atom < local.num_local; i_atom++) {

    std::vector<std::complex<double>> multi_moms(max_mmoms, 0.0);

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

                radial[ir] = std::pow(local.atom[i_atom].r_mesh[ir], l_mom) *
                    SSSol.zlr(ir, l, isp) * std::conj(SSSol.zlr(ir, lp, isp));

              }

              auto res = lsms::radialIntegral(radial, local.atom[i_atom].r_mesh, local.atom[i_atom].jmt);

              zl_zlp(l, lp, isp) = res;

            }

          }

          // Solve integral Z_l * J_l * delta_l,lp
          for (auto l = 0; l <= local.atom[i_atom].lmax; l++) {

            std::vector<std::complex<Real>> radial(local.atom[i_atom].jmt);

            for (auto ir = 0; ir < local.atom[i_atom].jmt; ir++) {

              radial[ir] = std::pow(local.atom[i_atom].r_mesh[ir], l_mom) *
                  SSSol.zlr(ir, l, isp) * std::conj(SSSol.jlr(ir, l, isp));

            }

            auto res = lsms::radialIntegral(radial, local.atom[i_atom].r_mesh, local.atom[i_atom].jmt);

            zl_jlp(l, isp) = res;

          }


          for (auto m_mom = -l_mom; m_mom <= l_mom; m_mom++) {
            auto il_mom = l_mom * l_mom + l_mom + m_mom;

            std::complex<double> factor = 0.0;

            for (auto l = 0; l <= local.atom[i_atom].lmax; l++) {
              for (auto m = -l; m <= l; m++) {

                for (auto lp = 0; lp <= local.atom[i_atom].lmax; lp++) {
                  for (auto mp = -lp; mp <= lp; mp++) {

                    auto kl = l * l + l + m;
                    auto klp = lp * lp + lp + mp;

                    std::complex<double> tau00 = 0.0;

                    if (isp == 0) {
                      tau00 = tau00_l(0, i_atom);
                    } else {
                      tau00 = tau00_l(0, i_atom + local.num_local);
                    }

                    auto gf_LLp = zl_zlp(l, lp, isp) * tau00 - zl_jlp(l, isp);

                    //gaunt = gaunt_real_table(l1, l2, l, m1, m2, m)

                    //std::complex<double> gaunt = GauntCoeficients::cgnt(j3, kl, klp);

                    //multi_moms[il_mom] += gaunt * std::imag(gf_LLp * dele1[iie]);

                  }

                }
              }

            }

          }

        }

      }
    }

  } // Local atoms

} // Energy contour



} // lsms

#endif //LSMS_SRC_MULTIPOLEMADELUNG_MULTIPOLEMOMENTS_HPP_
