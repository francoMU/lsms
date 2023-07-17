//
// Created by F.Moitzi on 15.07.2023.
//

#ifndef LSMS_SRC_STRUCTURECONSTANTS_STRUCTURE_H_
#define LSMS_SRC_STRUCTURECONSTANTS_STRUCTURE_H_

#include "NDArray.hpp"

namespace SC {


class Structure {

 private:

  NDArray<double, 2> lattice;
  NDArray<double, 2> coords;
  std::vector<int> species;

 public:

  Structure(NDArray<double, 2> &lattice,
      NDArray<double, 2> &coords,
      std::vector<int> &species);



};


}

#endif //LSMS_SRC_STRUCTURECONSTANTS_STRUCTURE_H_
