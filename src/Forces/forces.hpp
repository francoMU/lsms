
#ifndef LSMS_FORCES_HPP
#define LSMS_FORCES_HPP

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"
#include "Madelung/Multipole/higherOrderMadelung.hpp"

namespace lsms {

  /**
   * @brief get the counts and the displacement vectors for the gather operations
   *
   * @param comm
   * @param lsms
   * @param local
   * @param rcounts
   * @param displs
   */
  void gatherStructure(const LSMSCommunication &comm,
                       const LSMSSystemParameters &lsms,
                       const CrystalParameters &crystal,
                       std::vector<int> &rcounts,
                       std::vector<int> &displs);

  /**
   * @brief prints charges on each atom to the screen
   *
   * Here parameter are printed to screen that are relevant for the force calculation
   *
   * @param comm
   * @param lsms
   * @param local
   * @param crystal
   */
  void displayCharges(const LSMSCommunication &comm,
                      const LSMSSystemParameters &lsms,
                      const LocalTypeInfo &local,
                      const CrystalParameters &crystal);


  /**
   * @brief writes the forces to screen
   *
   * @param local
   * @param crystal
   */
  void displayForces(const LSMSCommunication &comm,
                     const LSMSSystemParameters &lsms,
                     const LocalTypeInfo &local,
                     const CrystalParameters &crystal,
                     int rank = 0);


  /**
   * @brief normalize the forces within the unitcell
   *
   * @param comm
   * @param lsms
   * @param local
   * @param crystal
   */
  void normalizeForces(const LSMSCommunication &comm,
                       const LSMSSystemParameters &lsms,
                       LocalTypeInfo &local,
                       const CrystalParameters &crystal);

  /**
   * @brief Gathers all forces form all processes
   *
   * @param comm
   * @param lsms
   * @param local
   * @param crystal
   * @param global_forces
   * @param root
   */
  void gatherForces(const LSMSCommunication &comm,
                    const LSMSSystemParameters &lsms,
                    const LocalTypeInfo &local,
                    const CrystalParameters &crystal,
                    std::vector<std::array<Real, 3>> &global_forces,
                    int root = 0
  );


  /**
   * @brief calculates forces
   *
   * This function calculates the forces on the atoms using the Hellmann-Feynmann
   * Theorem
   *
   * @param comm
   * @param lsms
   * @param local
   * @param crystal
   * @param parameters
   */
  void calculateForces(const LSMSCommunication &comm,
                       const LSMSSystemParameters &lsms,
                       LocalTypeInfo &local,
                       CrystalParameters &crystal,
                       ForceParameters parameters = ForceParameters{});

  /**
   * @brief calculates forces
   *
   * This function calculates the forces
   *
   * @param comm
   * @param lsms
   * @param local
   * @param crystal
   * @param parameters
   */
  void calculateForces(const LSMSCommunication &comm,
                       const LSMSSystemParameters &lsms,
                       LocalTypeInfo &local,
                       CrystalParameters &crystal,
                       HigherOrderMadelung &madelung,
                       ForceParameters parameters = ForceParameters{}
  );


}

#endif //LSMS_FORCES_HPP
