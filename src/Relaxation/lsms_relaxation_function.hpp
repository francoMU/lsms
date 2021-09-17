
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

   private:

      constexpr static int STARTER_NUMBER_OF_EVALUATIONS = 1;

      int number_of_evaluations{STARTER_NUMBER_OF_EVALUATIONS};
      bool write_to_file{false};

      void writeStructureToFile() const;

      std::string generateFileName() const;

   public:

      LsmsRelaxationFunction(LSMSSystemParameters &lsms,
                             LSMSCommunication &comm,
                             CrystalParameters &crystal,
                             LocalTypeInfo &local,
                             MixingParameters &mix,
                             bool write_to_file = false);

      void evaluate(const std::vector<Real> &coordinates_vector,
                    Real &result,
                    std::vector<Real> &gradient) override;

      LSMSSystemParameters &lsms;
      LSMSCommunication &comm;
      CrystalParameters &crystal;
      LocalTypeInfo &local;
      MixingParameters &mix;

      static std::string generateFileName(int number,
                                          const std::string &suffix = {"POSCAR"},
                                          const std::string &prefix = {"step"});

      int getNumberOfEvaluations() const;



   };

}


#endif //LSMS_RELAXATION_FUNCTION_HPP
