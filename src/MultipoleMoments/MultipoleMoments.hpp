//
// Created by F.Moitzi on 18.10.2021.
//

#ifndef LSMS_MULTIPOLEMOMENTS_HPP
#define LSMS_MULTIPOLEMOMENTS_HPP


#include <vector>

#include "Complex.hpp"

#include "AtomData.hpp"

namespace lsms {


  class Qlm {

  public:

    int lmax;

    int lmax_qlm;
    int nlm;

    std::vector<Complex> value;


    explicit Qlm(int lmax, int lmax_qlm) : lmax{lmax}, lmax_qlm{lmax_qlm} {

      nlm = (lmax_qlm + 1) * (lmax_qlm + 1);

      value = std::vector<Complex>(nlm, 0.0);

    };

    Complex &operator[](int i) {
      return value[i];
    }

    template<class T>
    Qlm &operator=(const T &val) {

      for (auto i = 0; i < nlm; i++) {
        value[i] = val;
      }

    }


  };

  void create_multipole_moments(Qlm &qlm, int lmax_qlm, int lmax, std::vector<double> &r_mesh, double ws_rad);

}


#endif //LSMS_MULTIPOLEMOMENTS_HPP
