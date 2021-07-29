
#ifndef LSMS_FORCES_HPP
#define LSMS_FORCES_HPP

#include "Communication/LSMSCommunication.hpp"
#include "Main/SystemParameters.hpp"
#include "Madelung/Multipole/higherOrderMadelung.hpp"

namespace lsms {


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
    void displayForces(const LocalTypeInfo &local,
                       const CrystalParameters &crystal
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
