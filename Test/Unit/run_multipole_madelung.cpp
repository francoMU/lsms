//
// Created by F.Moitzi on 23.06.2022.
//

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include "Madelung/Madelung.hpp"
#include "Main/SystemParameters.hpp"
#include "Misc/Coeficients.hpp"
#include "Misc/Indices.hpp"
#include "MultipoleMadelung/MultipoleMadelung.hpp"
#include "MultipoleMadelung/calculateMultipoleMadelung.hpp"
#include "accel_common.hpp"

using namespace lsms;

SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;
GauntCoeficients gauntCoeficients;
IFactors iFactors;

int main(int argc, char *argv[]) {


  /**
   * Create structure with 20000 atoms and 10 atoms locally
   */

  int num_atoms = 3;
  int lmax = 3;

  LSMSSystemParameters lsms;
  lsms.global.iprint = -1;

  AngularMomentumIndices::init(2 * lmax);
  SphericalHarmonicsCoeficients::init(2 * lmax);
  GauntCoeficients::init(lsms);
  IFactors::init(lsms, lmax);

  LocalTypeInfo local;
  local.setNumLocal(2);
  local.global_id = {0, 2};
  //
  LocalTypeInfo local_modern = local;
  //
  CrystalParameters crystal;
  crystal.resize(num_atoms);
  crystal.num_atoms = num_atoms;
  crystal.num_types = num_atoms;

  crystal.bravais(0, 0) = 1.1;
  crystal.bravais(1, 0) = 0.2;
  crystal.bravais(2, 0) = 0.5;

  crystal.bravais(0, 1) = -0.1;
  crystal.bravais(1, 1) = 1.1;
  crystal.bravais(2, 1) = 0.1;

  crystal.bravais(0, 2) = 0.0;
  crystal.bravais(1, 2) = 0.0;
  crystal.bravais(2, 2) = 1.6;

  crystal.position(0, 0) = 0.0;
  crystal.position(1, 0) = 0.0;
  crystal.position(2, 0) = 0.1;

  crystal.position(0, 1) = 0.5;
  crystal.position(1, 1) = 0.5;
  crystal.position(2, 1) = 0.4;

  crystal.position(0, 2) = 0.75;
  crystal.position(1, 2) = 0.65;
  crystal.position(2, 2) = 0.4;

  lsms::MultipoleMadelung madelung(lsms, crystal, local_modern, lmax);

  return EXIT_SUCCESS;
}

