//
// Created by F.Moitzi on 27.10.2021.
//

#include "box_relaxation_function.hpp"

#include "lsms_relaxation.hpp"


static void rescale_lattice(CrystalParameters &crystal, const arma::vec &values) {

  for (int i = 0; i < 3; i++) {
    crystal.bravais(0, i) *= values[0];
    crystal.bravais(1, i) *= values[1];
    crystal.bravais(2, i) *= values[2];
  }

}

static void rescale_atoms(CrystalParameters &crystal, const arma::vec &values) {

  for (int i = 0; i < crystal.num_atoms; i++) {
    crystal.position(0, i) *= values[0];
    crystal.position(1, i) *= values[1];
    crystal.position(2, i) *= values[2];
  }

}

double lsms::total_energy(const arma::vec &values, arma::vec *gradient, void *data) {

  // Convert pointer
  auto *data_converter = (BoxRelaxData *) data;

  arma::vec scaling = arma::ones(3, 1);

  // Adapt crystal structure
  CrystalParameters rescaled_crystal = *(data_converter->crystal);


  if (data_converter->iso) {
    scaling *= values[0];
  } else {

    int i = 0;

    if (data_converter->x) {
      scaling *= values[i];
      i++;
    }
    if (data_converter->y) {
      scaling *= values[i];
      i++;
    }
    if (data_converter->z) {
      scaling *= values[i];
      i++;
    }
  }

  rescale_lattice(rescaled_crystal, scaling);
  rescale_atoms(rescaled_crystal, scaling);


  // Perform calculation
  lsms::run_dft_calculation(
      *(data_converter->lsms),
      *(data_converter->comm),
      rescaled_crystal,
      *(data_converter->local),
      *(data_converter->mix),
      data_converter->reload_potential);

  data_converter->reload_potential = false;

  // Evaluate total energy
  auto total_energy = extractTotalEnergy(*(data_converter->lsms));


  return total_energy;
}
