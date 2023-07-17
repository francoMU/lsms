//
// Created by F.Moitzi on 15.07.2023.
//

#ifndef LSMS_SRC_STRUCTURECONSTANTS_DCOEFF_H_
#define LSMS_SRC_STRUCTURECONSTANTS_DCOEFF_H_

#include <complex>

namespace StructConst {


/**
 * The only thing that D1 needs to know are the cutoffs and the nu parameter
 *
 * Ideally there should be an algoithm to determine this parameters
 *
 */

class Dcoeff {

 public:

  void computeD1ij(std::complex<double> &d1);





};

}

#endif //LSMS_SRC_STRUCTURECONSTANTS_DCOEFF_H_
