
#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>


namespace TestPOSCAR {

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
   static std::pair<std::vector<T>, std::vector<T>> unique(std::vector<T> vec) {

      auto sorted_vec = vec;

      std::sort(sorted_vec.begin(), sorted_vec.end());
      auto last = unique(sorted_vec.begin(), sorted_vec.end());
      sorted_vec.erase(last, sorted_vec.end());

      std::vector<int> counts = sorted_vec;

      for (int i = 0; i < counts.size(); ++i) {
         counts[i] = std::count(vec.begin(), vec.end(), sorted_vec[i]);
      }

      return std::make_pair(vec, counts);

   };

   bool test_basic_ordering() {

      /*
       *
       */

      std::vector<int> type = {
            2, 1, 3, 4, 3, 3, 3, 5, 1, 2
      };

      auto copy_type = type;

      std::sort(copy_type.begin(), copy_type.end());
      auto last = unique(copy_type.begin(), copy_type.end());
      copy_type.erase(last, copy_type.end());

      for (auto &i : copy_type) {
         std::cout << i << std::endl;
      }

      assert(copy_type == std::vector<int>({1,
                                            2,
                                            3,
                                            4,
                                            5}));

      /*
       *
       */

      std::vector<int> counts = copy_type;

      for (int i = 0; i < counts.size(); ++i) {
         counts[i] = std::count(type.begin(), type.end(), copy_type[i]);
      }


      assert(counts == std::vector<int>({2,
                                         2,
                                         4,
                                         1,
                                         1}));


      /*
       * Indices of counts
       */

      //
      auto index_list = find_all_indices(type, copy_type[0]);
      assert(index_list == std::vector<int>({1, 8}));

      //
      index_list = find_all_indices(type, copy_type[1]);
      assert(index_list == std::vector<int>({0, 9}));


      return true;
   }


}


int main(int argc, char *argv[]) {

   auto test_result =
         TestPOSCAR::test_basic_ordering();

   if (test_result) {
      return EXIT_SUCCESS;
   } else {
      return EXIT_FAILURE;
   }

}