//
// Created by F.Moitzi on 23.08.2021.
//

#ifndef LSMS_RELAX_PARAMS_HPP
#define LSMS_RELAX_PARAMS_HPP

namespace lsms {


  constexpr double TOLERANCE = 1e-8;
  constexpr double INITIAL_SIGMA = 0.0001;
  constexpr bool WRITE_TO_FILE = false;
  constexpr int NO_OPTIMIZATION = -1;
  constexpr int MAX_ITERATIONS = NO_OPTIMIZATION;
  constexpr bool BOX_RELAX = true;
  constexpr bool X = true;
  constexpr bool Y = true;
  constexpr bool Z = true;
  constexpr bool ISO = true;

  class RelaxParams {


  public:

    explicit RelaxParams(int maxIterations = MAX_ITERATIONS,
                         double tolerance = TOLERANCE,
                         double initialSigma = INITIAL_SIGMA,
                         bool write_to_file = WRITE_TO_FILE,
                         bool box_relax = true,
                         bool x = X, bool y = Y, bool z = Z,
                         bool iso = ISO);

    int max_iterations;

    double tolerance;

    double initial_sigma;

    bool write_to_file;

    /**
     * Box relaxation parameters
     *
     * if activated the box will be relaxed instead of the coordinates
     *
     */
    bool box_relax;

    /**
     * Box coordinates to relax
     */
    bool x;
    bool y;
    bool z;
    bool iso;

    bool isOptimizationRun() const;

    bool isBoxOptimizationRun() const;

  };


}


#endif //LSMS_RELAX_PARAMS_HPP
