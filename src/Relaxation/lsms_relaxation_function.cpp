
#include "lsms_relaxation_function.hpp"

#include "Relaxation/lsms_relaxation.hpp"

void lsms::LsmsRelaxationFunction::evaluate(const std::vector<Real> &coordinates,
                                            double &result,
                                            std::vector<Real> &gradient) {

    lsms::VectorToCoordinates(lsms,
                              comm,
                              crystal,
                              local,
                              coordinates);

    // Main DFT run
    lsms::run_dft_calculation(lsms, comm, crystal, local, mix);

    // Gives the results
    auto results = lsms::extractTotalEnergy(lsms);
    if (lsms.global.iprint >= 0) {
        std::printf("Results: %lf\n", results);
    }

    // Gives the gradient
    lsms::ForcesToGradientVector(lsms,
                                 comm,
                                 crystal,
                                 local,
                                 gradient);

    // Set the total energy as result

    result = lsms.totalEnergy;


}

lsms::LsmsRelaxationFunction::LsmsRelaxationFunction(LSMSSystemParameters &lsms,
                                                     LSMSCommunication &comm,
                                                     CrystalParameters &crystal,
                                                     LocalTypeInfo &local,
                                                     MixingParameters &mix)
        : lsms(lsms),
          comm(comm),
          crystal(crystal),
          local(local),
          mix(mix) {}
