//
// Created by F.Moitzi on 20.09.2021.
//

#include "io_utils.hpp"
#include "chain_stream.hpp"


std::string lsms::createCoreOutput(LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local) {

  std::stringstream sstream;

  for (int i = 0; i < local.num_local; i++) {

    auto &atom = local.atom[i];

    sstream << "Atom: ["<< i << ":" << local.global_id[i] << "]: " << std::endl;

    for (int is = 0; is < lsms.n_spin_pola; is++) {

      for (int ic = 0; ic < atom.numc; ic++) {

        sstream << std::setw(3) << ic << " ";
        sstream << std::setw(3) << atom.nc(ic, is) << " ";
        sstream << std::setw(3) << atom.lc(ic, is) << " ";
        sstream << std::setw(3) << atom.kc(ic, is) << " ";
        sstream << std::fixed << std::showpoint;
        sstream << std::setw(28) << std::setprecision(18) << atom.ec(ic, is) << " ";
        std::string core_state_type {atom.coreStateType(ic, is)};
        sstream << std::setw(3) << core_state_type;
        sstream << std::endl;

      }

      sstream << std::fixed << std::showpoint;
      sstream << std::setw(28) << std::setprecision(18);
      sstream << atom.ecorv[is] << std::endl;
      sstream << std::fixed << std::showpoint;
      sstream << std::setw(28) << std::setprecision(18);
      sstream << atom.esemv[is] << std::endl;

    }
  }

  sstream << std::endl;

  return sstream.str();
}

void lsms::printCoreOutput(LSMSSystemParameters &lsms,
                           CrystalParameters &crystal,
                           LocalTypeInfo &local,
                           LSMSCommunication &comm) {

  ChainStream cs(std::cout, comm.rank, comm.size, comm.comm);

  auto core_output = createCoreOutput(lsms,
                                      crystal,
                                      local);

  cs << core_output;

}
