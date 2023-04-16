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
                               std::vector<std::vector<NonRelativisticSingleScattererSolution>> &solutionNonRel,
                               std::vector<Matrix<Complex>> &tau00_l,
                               std::vector<Complex> &dele1,
                               int l_mom_max);

void printMultipoleMoments(LSMSCommunication &comm,
                           LSMSSystemParameters &lsms,
                           LocalTypeInfo &local, int l_mom_max);

}

#endif