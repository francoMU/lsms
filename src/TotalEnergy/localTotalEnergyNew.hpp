//
// Created by F.Moitzi on 10.10.2021.
//

#ifndef LSMS_LOCALTOTALENERGYNEW_HPP
#define LSMS_LOCALTOTALENERGYNEW_HPP

#include "Real.hpp"

#include "Main/SystemParameters.hpp"

namespace lsms {

  void localTotalEnergyNew(LSMSSystemParameters &lsms, AtomData &atom, Real &energy, Real &pressure);

}

#endif //LSMS_LOCALTOTALENERGYNEW_HPP
