//
// Created by F.Moitzi on 24.08.2022.
//

#ifndef LSMS_STATES_HPP
#define LSMS_STATES_HPP

namespace LSMS {

class States {
 public:

  int core_charge;

  std::vector<int> n;
  std::vector<int> l;
  std::vector<int> spin;
  std::vector<int> kappa;
  std::vector<double> occupation;


  static void relativistic_atomic_states(int core_charge, std::vector<int> &n,
                                         std::vector<int> &l,
                                         std::vector<int> &spin,
                                         std::vector<int> &kappa,
                                         std::vector<double> &occupation);

  static void nonrelativistic_atomic_states(int core_charge,
                                            std::vector<int> &n,
                                            std::vector<int> &l,
                                            std::vector<double> &occupation);
};

}  // namespace LSMS

#endif  // LSMS_STATES_HPP
