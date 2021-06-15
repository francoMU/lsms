//
// Created by F.Moitzi on 23.08.2021.
//

#include "relax_params.hpp"

lsms::RelaxParams::RelaxParams(int maxIterations,
                               double tolerance,
                               double initialSigma,
                               bool write_to_file) :
      max_iterations(
            maxIterations),
      tolerance(tolerance),
      initial_sigma(initialSigma),
      write_to_file(write_to_file) {}
