//
// Created by F.Moitzi on 26.08.2021.
//

#ifndef LSMS_POSCARSTRUCTURETYPE_HPP
#define LSMS_POSCARSTRUCTURETYPE_HPP

#include <map>

namespace lsms {

   enum POSCARStructureType {
      Cartesian,
      Direct,
      Fractional
   };

   static std::map <POSCARStructureType, std::string> structureTypeMap = {
         {POSCARStructureType::Cartesian,  "Cartesian"},
         {POSCARStructureType::Direct,     "Direct"},
         {POSCARStructureType::Fractional, "Fractional"},
   };

}

#endif //LSMS_POSCARSTRUCTURETYPE_HPP
