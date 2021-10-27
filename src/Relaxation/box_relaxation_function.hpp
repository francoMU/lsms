//
// Created by F.Moitzi on 27.10.2021.
//

#define OPTIM_ENABLE_ARMA_WRAPPERS


#include "armadillo"
#include "optim.hpp"

#include "SystemParameters.hpp"
#include "LSMSCommunication.hpp"
#include "mixing_params.hpp"

#ifndef LSMS_BOX_RELAXATION_FUNCTION_HPP
#define LSMS_BOX_RELAXATION_FUNCTION_HPP

namespace lsms {


  struct BoxRelaxData {
    LSMSSystemParameters *lsms;
    LSMSCommunication *comm;
    CrystalParameters *crystal;
    LocalTypeInfo *local;
    MixingParameters *mix;
    bool iso;
    bool x;
    bool y;
    bool z;
    bool reload_potential;
  };


  /**
   *
   * @param values
   * @param gradient
   * @param data
   * @return
   */
  double total_energy(const arma::vec &values, arma::vec *gradient, void *data);


}


#endif //LSMS_BOX_RELAXATION_FUNCTION_HPP
