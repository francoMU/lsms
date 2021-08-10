
#ifndef LSMS_RELAXATION_FUNCTION_HPP
#define LSMS_RELAXATION_FUNCTION_HPP

#include <vector>

#include "Real.hpp"
#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Main/mixing.hpp"

#include "base_function.hpp"

namespace lsms {

    class LsmsRelaxationFunction : public BaseFunction<Real> {

    public:

        LsmsRelaxationFunction(LSMSSystemParameters &lsms,
                               LSMSCommunication &comm,
                               CrystalParameters &crystal,
                               LocalTypeInfo &local,
                               MixingParameters &mix);

        void evaluate(const std::vector<Real> &coordinates_vector,
                      double &result,
                      std::vector<Real> &gradient) override;

        LSMSSystemParameters &lsms;
        LSMSCommunication &comm;
        CrystalParameters &crystal;
        LocalTypeInfo &local;
        MixingParameters &mix;


    };

}


#endif //LSMS_RELAXATION_FUNCTION_HPP
