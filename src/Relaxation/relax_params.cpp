//
// Created by F.Moitzi on 23.08.2021.
//

#include "relax_params.hpp"

lsms::RelaxParams::RelaxParams(int maxIterations,
                               double tolerance,
                               double initialSigma,
                               bool write_to_file,
                               bool box_relax,
                               bool x,
                               bool y,
                               bool z,
                               bool iso) :
    max_iterations(maxIterations),
    tolerance(tolerance),
    initial_sigma(initialSigma),
    write_to_file(write_to_file),
    box_relax{BOX_RELAX},
    x{X},
    y{Y},
    z{Z},
    iso{ISO} {

}

bool lsms::RelaxParams::isBoxOptimizationRun() const {
  return (max_iterations != NO_OPTIMIZATION) & box_relax;
}


bool lsms::RelaxParams::isOptimizationRun() const {
  return max_iterations != NO_OPTIMIZATION;
}
