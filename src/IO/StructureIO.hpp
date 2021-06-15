//
// Created by F.Moitzi on 24.08.2021.
//

#ifndef LSMS_STRUCTUREIO_HPP
#define LSMS_STRUCTUREIO_HPP

#include <string>
#include <ostream>
#include <algorithm>
#include <cctype>
#include <locale>

#include "Main/SystemParameters.hpp"


namespace lsms {


   template<typename T>
   static std::vector<int> find_all_indices(const std::vector<T> &vec, T value) {

      std::vector<int> indices;

      auto it = vec.begin();

      while ((it = std::find_if(it,
                                vec.end(),
                                [value](T x) { return x == value; })) != vec.end()) {
         indices.push_back(std::distance(vec.begin(), it));
         it++;
      }

      return indices;

   }

   template<typename T>
   static std::pair<std::vector<T>, std::vector<int>> unique(std::vector<T> vec) {

      auto sorted_vec = vec;

      std::sort(sorted_vec.begin(), sorted_vec.end());
      auto last = unique(sorted_vec.begin(), sorted_vec.end());
      sorted_vec.erase(last, sorted_vec.end());

      std::vector<int> counts;

      for (auto ptr = sorted_vec.begin(); ptr < sorted_vec.end(); ptr++) {
         counts.emplace_back(std::count(vec.begin(), vec.end(), *ptr));
      }

      return std::make_pair(sorted_vec, counts);

   };

   template<class C = char, class T = std::char_traits<C> >
   class StructureIO {

      /*
       * http://www.gotw.ca/gotw/048.htm
       */

      virtual void writeToStream(std::basic_ostream<C, T> &out,
                                 const LSMSSystemParameters &lsms,
                                 CrystalParameters &crystal) = 0;

   };

   // trim from start (in place)
   template<typename C = unsigned char>
   static inline void ltrim(std::string &s) {
      s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](C ch) {
         return !std::isspace(ch);
      }));
   }

   // trim from end (in place)
   template<typename C = unsigned char>
   static inline void rtrim(std::string &s) {
      s.erase(std::find_if(s.rbegin(), s.rend(), [](C ch) {
         return !std::isspace(ch);
      }).base(), s.end());
   }

   // trim from both ends (in place)
   static inline void trim(std::string &s) {
      ltrim(s);
      rtrim(s);
   }

   // trim from start (copying)
   static inline std::string ltrim_copy(std::string s) {
      ltrim(s);
      return s;
   }

   // trim from end (copying)
   static inline std::string rtrim_copy(std::string s) {
      rtrim(s);
      return s;
   }

   // trim from both ends (copying)
   static inline std::string trim_copy(std::string s) {
      trim(s);
      return s;

   }


}


#endif //LSMS_STRUCTUREIO_HPP
