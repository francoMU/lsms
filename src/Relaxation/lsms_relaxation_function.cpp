
#include "lsms_relaxation_function.hpp"

#include <iostream>
#include <fstream>

#include "Relaxation/lsms_relaxation.hpp"
#include "IO/POSCARStructureType.hpp"
#include "IO/POSCARStructureIO.hpp"

void lsms::LsmsRelaxationFunction::evaluate(const std::vector<Real> &coordinates,
                                            double &result,
                                            std::vector<Real> &gradient) {


   lsms::VectorToCoordinates(lsms,
                             comm,
                             crystal,
                             local,
                             coordinates);

   if (write_to_file) {
      writeStructureToFile();
   }

   // Main DFT run
   bool reload_potentials = number_of_evaluations == STARTER_NUMBER_OF_EVALUATIONS;
   lsms::run_dft_calculation(lsms, comm, crystal, local, mix, reload_potentials);
   number_of_evaluations++;

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
                                                     MixingParameters &mix,
                                                     bool write_to_file)
      : write_to_file{write_to_file},
        lsms(lsms),
        comm(comm),
        crystal(crystal),
        local(local),
        mix(mix) {


}

void lsms::LsmsRelaxationFunction::writeStructureToFile() const {

   if (lsms.global.iprint >= 0) {
      std::ofstream file;
      file.open(generateFileName().c_str());
      POSCARStructureIO io(POSCARStructureType::Cartesian);
      io.writeToStream(file, lsms, crystal);
      file.close();
   }

}

std::string lsms::LsmsRelaxationFunction::generateFileName() const {
   return generateFileName(number_of_evaluations);
}

std::string lsms::LsmsRelaxationFunction::generateFileName(int number,
                                                           const std::string &suffix,
                                                           const std::string &prefix) {
   auto sep = std::string(".");
   return prefix + sep + std::to_string(number) + sep + suffix;
}

int lsms::LsmsRelaxationFunction::getNumberOfEvaluations() const {
   return number_of_evaluations;
}



