//
// Created by F.Moitzi on 17.09.2021.
//

#ifndef LSMS_ORDEREDPRINT_HPP
#define LSMS_ORDEREDPRINT_HPP

#include <mpi.h>

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

namespace lsms {

  void gatherEnergies(const LSMSCommunication &comm,
                      const LSMSSystemParameters &lsms,
                      const LocalTypeInfo &local,
                      const CrystalParameters &crystal,
                      std::vector<double> &global_energies,
                      int root);

  void gatherInfo(const LSMSCommunication &comm,
                  const LSMSSystemParameters &lsms,
                  const LocalTypeInfo &local,
                  const CrystalParameters &crystal,
                  std::vector<double> &rmt,
                  std::vector<double> &rws,
                  std::vector<double> &rInscribed,
                  int root);


}

#endif //LSMS_ORDEREDPRINT_HPP
