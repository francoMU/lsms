//
// Created by F.Moitzi on 20.09.2021.
//

#ifndef LSMS_IO_UTILS_HPP
#define LSMS_IO_UTILS_HPP

#include <string>

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"


namespace lsms {

  std::string createCoreOutput(LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local);

  void printCoreOutput(LSMSSystemParameters &lsms,
                       CrystalParameters &crystal,
                       LocalTypeInfo &local,
                       LSMSCommunication &comm);

}



#endif //LSMS_IO_UTILS_HPP
