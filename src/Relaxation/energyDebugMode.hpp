//
// Created by F.Moitzi on 17.09.2021.
//



#ifndef LSMS_ENERGYDEBUGMOD_HPP
#define LSMS_ENERGYDEBUGMOD_HPP

#include <type_traits>

namespace lsms {

   enum class EnergyDebugMode {
      None = 0x00,
      Local = 0x01,
      Global = 0x02
   };

}

#endif //LSMS_ENERGYDEBUGMOD_HPP
