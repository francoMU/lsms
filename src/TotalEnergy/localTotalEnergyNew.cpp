//
// Created by F.Moitzi on 10.10.2021.
//

#include "localTotalEnergyNew.hpp"

#include "integrator.hpp"
#include "poisson.hpp"
#include "PhysicalConstants.hpp"


void lsms::localTotalEnergyNew(LSMSSystemParameters &lsms, AtomData &atom, Real &energy, Real &pressure) {


  auto &energyStruct = atom.energyStruct;

  Real eigenvalueSum = 0.0;
  Real kineticEnergy = 0.0;
  Real coulombEnergy = 0.0;
  Real kohn_sham = 0.0;
  Real xcEnergy = 0.0;
  Real ezpt = 0.0;
  Real tpzpt = 0.0;

  std::vector<Real> integrand(atom.r_mesh.size());
  std::vector<Real> vhartree(atom.r_mesh.size());
  std::vector<Real> vhartreederiv(atom.r_mesh.size());

  energy = 0.0;
  pressure = 0.0;

  Real rSphere;
  switch (lsms.mtasa) {
    case 1:
      rSphere = atom.rws;
      break;
    case 2:
      rSphere = atom.rws;
      break;
    default:
      rSphere = atom.rInscribed;
  }


  energyStruct.zero_point_energy = ezpt;


  /*
   * Kinetic energy
   */


  if (lsms.n_spin_pola == 1) {
    eigenvalueSum = atom.evalsum[0] + atom.esemv[0];

    kineticEnergy = atom.ecorv[0] + atom.esemv[0]; // (1)
    kineticEnergy += atom.evalsum[0]; // (2)

    energyStruct.one_electron = kineticEnergy;


    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = (atom.rhoNew(i, 0) * atom.vr(i, 0)) / (atom.r_mesh[i]);
    }
  } else {
    eigenvalueSum = atom.evalsum[0] + atom.evalsum[1] + atom.esemv[0] + atom.esemv[1];

    kineticEnergy = atom.ecorv[0] + atom.ecorv[1] + atom.esemv[0] + atom.esemv[1]; // (1)
    kineticEnergy += atom.evalsum[0] + atom.evalsum[1]; // (2)

    energyStruct.one_electron = kineticEnergy;

    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] =
          (atom.rhoNew(i, 0) * atom.vr(i, 0) + atom.rhoNew(i, 1) * atom.vr(i, 1)) / (atom.r_mesh[i]);
    }
  }


  kohn_sham = radialIntegral(integrand, atom.r_mesh, rSphere); // (3)
  energyStruct.kohn_sham = kohn_sham;
  kineticEnergy -= kohn_sham;


  if (lsms.global.iprint > 0) {
    printf("evssum                      = %35.25lf Ry\n", eigenvalueSum);
    printf("kinetic Energy              = %35.25lf Ry\n", kineticEnergy);
  }


  /*
   * Coloumbic energy
   */


  if (lsms.n_spin_pola == 1) {
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = - 2.0 * (atom.rhoNew(i, 0)) * atom.ztotss / (atom.r_mesh[i]);
    }
  } else { //spin polarized
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = - 2.0 * (atom.rhoNew(i, 0) + atom.rhoNew(i, 1)) * atom.ztotss / (atom.r_mesh[i]);
    }
  }
  auto ezrho = lsms::radialIntegral(integrand, atom.r_mesh, rSphere);
  energyStruct.core_interaction = ezrho;

  if (lsms.global.iprint > 0) {
    printf("ezrho                       = %35.25lf Ry\n", ezrho);
  }

  if (lsms.n_spin_pola == 1) {
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = atom.rhoNew(i, 0);
    }
  } else { //spin polarized
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = (atom.rhoNew(i, 0) + atom.rhoNew(i, 1));
    }
  }

  radial_poisson(vhartree, vhartreederiv, atom.r_mesh, atom.getRmeshDerivative(), integrand, atom.jmt);
  for (auto i = 0; i < atom.jmt + 1; i++) {
    integrand[i] = vhartree[i] * integrand[i];
  }
  auto erho = radialIntegral(integrand, atom.r_mesh, rSphere);
  energyStruct.hartree = erho;


  if (lsms.global.iprint > 0) {
    printf("erho                        = %35.25lf Ry\n", erho);
  }

  coulombEnergy = erho + ezrho;

  if (lsms.global.iprint > 0) {
    printf("Coulomb Energy              = %35.25lf Ry\n", coulombEnergy);
  }

  /*
   * Exchange-correlation part
   */

  if (lsms.xcFunctional[0] == 0)  // for the built in xc functionals:
  {
    if (lsms.n_spin_pola == 1) {
      for (int i = 0; i < atom.r_mesh.size(); i++) {
        integrand[i] = (atom.rhoNew(i, 0) * atom.exchangeCorrelationEnergy(i, 0));
      }
    } else { // spin polarized
      for (int i = 0; i < atom.r_mesh.size(); i++) {
        integrand[i] = (atom.rhoNew(i, 0) * atom.exchangeCorrelationEnergy(i, 0)
                        + atom.rhoNew(i, 1) * atom.exchangeCorrelationEnergy(i, 1));
      }
    }

    xcEnergy = lsms::radialIntegral(integrand, atom.r_mesh, rSphere);
  } else if (lsms.xcFunctional[0] == 1) { // for libxc functionals
    if (lsms.n_spin_pola == 1) {
      for (int i = 0; i < atom.r_mesh.size(); i++) {
        integrand[i] = atom.exchangeCorrelationEnergy(i, 0) * (atom.rhoNew(i, 0));
      }
    } else { // spin polarized
      for (int i = 0; i < atom.r_mesh.size(); i++) {
        integrand[i] = atom.exchangeCorrelationEnergy(i, 0) * (atom.rhoNew(i, 0) + atom.rhoNew(i, 1));
      }
    }

    xcEnergy = lsms::radialIntegral(integrand, atom.r_mesh, rSphere); // (7)
  } else {  // unknown functional (we never should arrive here!
    printf("Unknown xc function in localTotalEnergy!\n");
    exit(1);
  }

  if (lsms.global.iprint > 0) {
    printf("Exchange-Correlation Energy = %35.25lf Ry\n", xcEnergy);
    printf("ezpt                        = %35.25lf Ry\n\n", ezpt);
  }

  energyStruct.exchange_correlation = xcEnergy;
  energy += kineticEnergy + coulombEnergy + xcEnergy; // + ezpt;

  /*
   * Longitudinal spin fluctuations
   */
  if (lsms.n_spin_pola == 2) {
    auto mag_mom = atom.mvalws;
    energy += -convertKtoRydberg * lsms.temperature * atom.lsf_functional.entropy(mag_mom);
  }


  /*
  * Pressure is not includequ
   */
  pressure = 0.0;

}
