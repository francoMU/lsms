//
// Created by F.Moitzi on 26.08.2021.
//

#ifndef LSMS_POSCARSTRUCTUREIO_HPP
#define LSMS_POSCARSTRUCTUREIO_HPP

#include "StructureIO.hpp"

#include "POSCARStructureType.hpp"


namespace lsms {

   class POSCARStructureIO : public StructureIO<> {

   private:

      POSCARStructureType type;

   public:

      explicit POSCARStructureIO(POSCARStructureType type);

      void writeToStream(std::basic_ostream<char,
                         std::char_traits<char>> &out,
                         const LSMSSystemParameters &lsms,
                         CrystalParameters &crystal) override;

   };


}


#endif //LSMS_POSCARSTRUCTUREIO_HPP
