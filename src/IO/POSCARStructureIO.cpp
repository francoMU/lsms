//
// Created by F.Moitzi on 26.08.2021.
//

#include "POSCARStructureIO.hpp"

void lsms::POSCARStructureIO::writeToStream(std::basic_ostream<char, std::char_traits<char>> &out,
                                            const LSMSSystemParameters &lsms,
                                            CrystalParameters &crystal) {

   constexpr int width = 16;
   constexpr int precision = 12;

   // Write Cubic BN
   out << trim_copy(lsms.systemid) << std::endl;

   out << std::setprecision(2) << std::fixed;
   out << 1.0 << std::endl;

   for (int i = 0; i < 3; i++) {
      out << " ";
      out << std::setprecision(precision) << std::fixed
          << std::setw(width) << crystal.bravais(i, 0) << " ";

      out << crystal.bravais(i, 1) << " ";

      out << crystal.bravais(i, 2) << std::endl;
   }

   /*
    * Atom types
    */

   std::vector<std::string> names;

   for (auto atom_types : crystal.types) {
      names.emplace_back(trim_copy(atom_types.name));
   }


   auto unique_pair = unique(names);

   auto unique_types = unique_pair.first;
   auto unique_counts = unique_pair.second;

   for (auto &unique_type : unique_types) {

      out << std::setw(4);
      out << unique_type << " ";

   }
   out << std::endl;

   for (auto &unique_count : unique_counts) {
      out << std::setw(4);
      out << unique_count << " ";
   }
   out << std::endl;

   /*
    * Cartesian coordinates
    *
    * only supported now
    *
    */

   out << structureTypeMap[type] << std::endl;

   //
   for (auto &unique_type : unique_types) {

      // Find all indices
      auto indices = find_all_indices(names, unique_type);


      for (auto index : indices) {
         out << std::setprecision(precision)
             << std::fixed
             << std::setw(width);
         out << crystal.position(0, index) << " ";
         out << crystal.position(1, index) << " ";
         out << crystal.position(2, index) << std::endl;
      }

   }


}

lsms::POSCARStructureIO::POSCARStructureIO(lsms::POSCARStructureType type)
      : type{type} {

}
