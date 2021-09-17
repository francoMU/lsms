//
// Created by F.Moitzi on 23.08.2021.
//

#ifndef LSMS_RELAX_PARAMS_HPP
#define LSMS_RELAX_PARAMS_HPP

namespace lsms {


   constexpr int MAX_ITERATIONS = 30;
   constexpr double TOLERANCE = 1e-8;
   constexpr double INITIAL_SIGMA = 0.0001;
   constexpr bool WRITE_TO_FILE = false;
   constexpr int NO_OPTIMIZATION = -1;

   class RelaxParams {


   public:

      explicit RelaxParams(int maxIterations = MAX_ITERATIONS,
                           double tolerance = TOLERANCE,
                           double initialSigma = INITIAL_SIGMA,
                           bool write_to_file = WRITE_TO_FILE);

      int max_iterations;

      double tolerance;

      double initial_sigma;

      bool write_to_file;

      bool isOptimizationRun() const;

   };


}


#endif //LSMS_RELAX_PARAMS_HPP
