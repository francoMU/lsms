
#ifndef LSMS_RELAXATION_H
#define LSMS_RELAXATION_H

namespace lsms {

    void run_dft_calculation(LSMSSystemParameters &lsms,
                             LSMSCommunication &comm,
                             CrystalParameters &crystal,
                             LocalTypeInfo &local,
                             MixingParameters &mix);


    Real extractTotalEnergy(LSMSSystemParameters &lsms);

    void VectorToCoordinates(LSMSSystemParameters &lsms,
                                LSMSCommunication &comm,
                                CrystalParameters &crystal,
                                LocalTypeInfo &local,
                                const std::vector<Real> &coordinates);

    void CoordinatesToVector(LSMSSystemParameters &lsms,
                             LSMSCommunication &comm,
                             CrystalParameters &crystal,
                             std::vector<Real> &coordinates);

    void ForcesToGradientVector(LSMSSystemParameters &lsms,
                                LSMSCommunication &comm,
                                CrystalParameters &crystal,
                                LocalTypeInfo &local,
                                std::vector<Real> &gradient);

}


#endif //LSMS_RELAXATION_H
