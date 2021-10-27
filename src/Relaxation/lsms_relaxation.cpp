
#include <mpi.h>

#include <Accelerator.hpp>
#include <buildLIZandCommLists.hpp>
#include <distributeAtoms.hpp>
#include "IO/io_utils.hpp"
#include "setupVorpol.hpp"
#include "Madelung.hpp"
#include "interpolatePotential.hpp"
#include "EnergyContourIntegration.hpp"
#include "calculateEvec.hpp"
#include "Coeficients.hpp"
#include "checkAntiFerromagneticStatus.hpp"
#include "calculateDensities.hpp"
#include "calculateChargesPotential.hpp"
#include "calculateTotalEnergy.hpp"
#include "Core/calculateCoreStates.hpp"
#include "calculateChemPot.hpp"
#include "higherOrderMadelung.hpp"
#include "calcHigherOrderMadelung.hpp"
#include "Forces/forces.hpp"
#include "lsms_relaxation.hpp"
#include "utils.hpp"

#include "SystemParameters.hpp"
#include "PotentialIO.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "mixing.hpp"
#include "energyDebugMode.hpp"
#include "orderedPrint.hpp"


static void prepareSpinOnAtoms(const LSMSSystemParameters &lsms,
                               LocalTypeInfo &local) {
  if (lsms.n_spin_cant > 1) {
    for (int i = 0; i < local.num_local; i++) {
      local.atom[i].get_b_basis();
      local.atom[i].spinFlipped = false;
    }
  } else {
    for (int i = 0; i < local.num_local; i++) {
      local.atom[i].reset_b_basis();
      local.atom[i].spinFlipped = false;
    }
  }
}

static void printVolumes(LSMSSystemParameters &lsms,
                         LSMSCommunication &comm,
                         CrystalParameters &crystal,
                         LocalTypeInfo &local) {

  std::vector<double> rmt;
  std::vector<double> rws;
  std::vector<double> rInscribed;
  lsms::gatherInfo(comm,
                   lsms,
                   local,
                   crystal,
                   rmt,
                   rws,
                   rInscribed,
                   0);

  if (comm.rank == 0) {
    for (int i = 0; i < crystal.num_atoms; i++) {
      std::printf("%5d %s [%d.%d]: %20.16f %20.16f %20.16f\n",
                  i + 1,
                  crystal.types[i].name,
                  crystal.types[i].node,
                  crystal.types[i].local_id,
                  rmt[i],
                  rws[i],
                  rInscribed[i]);
    }
  }

}

void lsms::getInterCoreEnergy(LSMSSystemParameters &lsms,
                              LSMSCommunication &comm,
                              CrystalParameters &crystal,
                              LocalTypeInfo &local,
                              Real &total_core_energy) {

  Real local_sum = 0.0;

  for (int i = 0; i < local.num_local; i++) {

    auto &atom = local.atom[i];

    if (lsms.n_spin_pola == 1) {
      local_sum += atom.ecorv[0] + atom.esemv[0];
    } else { //spin polarized
      local_sum += atom.evalsum[0] + atom.evalsum[1] + atom.esemv[0] + atom.esemv[1];
    }
    local_sum += atom.interstitialEnergy;
  }

  globalSum(comm, local_sum);
  total_core_energy = local_sum;

}


void lsms::getSiteLocalEnergy(LSMSSystemParameters &lsms, LSMSCommunication &comm, CrystalParameters &crystal,
                              LocalTypeInfo &local, Real &total_site_energy) {

  Real local_sum = 0.0;

  for (int i = 0; i < local.num_local; i++) {
    auto &atom = local.atom[i];
    local_sum += atom.siteLocalEnergy;
  }

  globalSum(comm, local_sum);
  total_site_energy = local_sum;
}


void lsms::getMagneticMoments(LSMSSystemParameters &lsms,
                              LSMSCommunication &comm,
                              CrystalParameters &crystal,
                              LocalTypeInfo &local,
                              Real &sum_mag) {

  Real local_sum_mag = 0.0;

  for (int i = 0; i < local.num_local; i++) {
    local_sum_mag += local.atom[i].mvalws;
  }

  local_sum_mag /= local.num_local;

  globalSum(comm, local_sum_mag);

  sum_mag = local_sum_mag;

}

void lsms::ForcesToGradientVector(LSMSSystemParameters &lsms,
                                  LSMSCommunication &comm,
                                  CrystalParameters &crystal,
                                  LocalTypeInfo &local,
                                  std::vector<Real> &gradient) {

  std::vector<std::array<Real, 3>> global_forces(crystal.num_atoms);

  gatherForces(comm, lsms, local, crystal, global_forces, 0);

  if (comm.rank == 0) {
    for (int i = 0; i < crystal.num_atoms; i++) {
      gradient[i * 3] = global_forces[i][0];
      gradient[i * 3 + 1] = global_forces[i][1];
      gradient[i * 3 + 2] = global_forces[i][2];
    }
  }

  MPI_Bcast(gradient.data(),
            gradient.size(),
            MPI_DOUBLE,
            0,
            comm.comm);


}

void lsms::CoordinatesToVector(LSMSSystemParameters &lsms,
                               LSMSCommunication &comm,
                               CrystalParameters &crystal,
                               std::vector<Real> &coordinates
) {

  for (int i = 0; i < crystal.num_atoms; i++) {
    coordinates[i * 3] = crystal.position(0, i);
    coordinates[i * 3 + 1] = crystal.position(1, i);
    coordinates[i * 3 + 2] = crystal.position(2, i);
  }

}


void lsms::VectorToCoordinates(LSMSSystemParameters &lsms,
                               LSMSCommunication &comm,
                               CrystalParameters &crystal,
                               LocalTypeInfo &local,
                               const std::vector<Real> &coordinates) {

  for (int i = 0; i < crystal.num_atoms; i++) {
    crystal.position(0, i) = coordinates[i * 3];
    crystal.position(1, i) = coordinates[i * 3 + 1];
    crystal.position(2, i) = coordinates[i * 3 + 2];
  }

}

Real lsms::extractTotalEnergy(LSMSSystemParameters &lsms) {
  return lsms.totalEnergy;
}


void lsms::run_dft_calculation(LSMSSystemParameters &lsms,
                               LSMSCommunication &comm,
                               CrystalParameters &crystal,
                               LocalTypeInfo &local,
                               MixingParameters &mix,
                               bool reload_potentials,
                               EnergyDebugMode mode
) {

  // Ensure that the store id is set to default value
  for (int i = 0; i < crystal.num_atoms; i++) {
    crystal.types[i].store_id = DEFAULT_STORE_ID;
  }


  //communicateParameters(comm, lsms, crystal, mix);

  if (comm.rank != lsms.global.print_node) {
    lsms.global.iprint = lsms.global.default_iprint;
  }

  /*
   * Setup of DFT calculation
   */

  lsms.angularMomentumIndices.init(2 * crystal.maxlmax);


  if (lsms.global.iprint >= 0) {
    std::printf("\nStart DFT calculation:\n");
  }


  auto num_local = distributeTypes(crystal, comm);
  local.setNumLocal(num_local);
  local.setGlobalId(comm.rank, crystal);


  if (reload_potentials) {
    local.setMaxPts(lsms.global.iprpts);
    local.setMaxCore(lsms.global.ipcore);
  }


  if (lsms.global.iprint >= 0) {
    printCompressedCrystalParameters(stdout, crystal);
  }

#if defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  deviceAtoms.resize(local.num_local);
#endif

  // Setup XC functional
  if (lsms.xcFunctional[0] == 1) {
    lsms.libxcFunctional.init(lsms.n_spin_pola, lsms.xcFunctional);
  }

  if (lsms.global.iprint >= 0) {
    if (lsms.xcFunctional[0] == 1) {
      std::printf("The functional [1]:'%s'\n", lsms.libxcFunctional.functional[0].info->name);
      std::printf("The functional [2]:'%s'\n", lsms.libxcFunctional.functional[1].info->name);

      if (lsms.libxcFunctional.needGradients) {
        std::printf("XC: GGA\n");
      } else {
        std::printf("XC: LDA\n");
      }

    }
  }

#if defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  deviceConstants.allocate(lsms.angularMomentumIndices, gauntCoeficients, iFactors);
#endif

  // Build LIZ
  double timeBuildLIZandCommList = MPI_Wtime();
  buildLIZandCommLists(comm, lsms, crystal, local);
  timeBuildLIZandCommList = MPI_Wtime() - timeBuildLIZandCommList;

  if (lsms.global.iprint >= 0) {
    printf("time for buildLIZandCommLists [num_local=%d]: %lf sec\n",
           local.num_local, timeBuildLIZandCommList);
    fflush(stdout);
  }

  // initialize the potential accelerators (GPU)
  // we need to know the max. size of the kkr matrix to invert: lsms.n_spin_cant*local.maxNrmat()
  // which is only available after building the LIZ
  acceleratorInitialize(lsms.n_spin_cant * local.maxNrmat(), lsms.global.GPUThreads);
  local.tmatStore.pinMemory();
#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
  // deviceStorage = allocateDStore();
  deviceStorage = new DeviceStorage;
#endif
#ifdef BUILDKKRMATRIX_GPU
  deviceConstants.resize(local.num_local);
  // for(int i=0; i<local.num_local; i++) deviceConstants[i]=allocateDConst();
#endif

  for (int i = 0; i < local.num_local; i++) {
    local.atom[i].pmat_m.resize(lsms.energyContour.groupSize());
  }


// set maximal number of radial grid points and core states if reading from bigcell file
  if (lsms.global.iprint >= 1) {
    for (int i = 0; i < local.num_local; i++) {
      fprintf(stdout, "LIZ for atom %d.%d on %d\n",
              i, local.global_id[i], comm.rank);
      printLIZInfo(stdout, local.atom[i]);
      fflush(stdout);
    }

  }

  if (lsms.global.iprint >= 1) {
    printLSMSGlobals(stdout, lsms);
    printLSMSSystemParameters(stdout, lsms);
    fprintf(stdout, "LIZ for atom 0 on this node\n");
    printLIZInfo(stdout, local.atom[0]);
    if (local.atom[0].forceZeroMoment)
      fprintf(stdout, "\nMagnetic moment of atom 0 forced to be zero!\n\n");
    printCommunicationInfo(stdout, comm);
  }
  fflush(stdout);

  /*
   * 1. Load the potentials
   */
  if (lsms.global.iprint >= 1) {
    std::cout << "Load potentials...." << std::endl;
  }

  if (reload_potentials) {
    loadPotentials(comm, lsms, crystal, local);
  } else {
    readjustPotential(comm, lsms, crystal, local);
  }

  /*
   * 2. Setup Voronoi tesselation
   */
  if (lsms.global.iprint >= 1) {
    std::cout << "Setup voronoi polyhedras...." << std::endl;
  }

  //if (reload_potentials) {
  setupVorpol(lsms, crystal,
              local, sphericalHarmonicsCoeficients);
  //}

  /*
   * 3. Interpolate the potentials
   */
  for (int i = 0; i < local.num_local; i++) {
    if (local.atom[i].generateNewMesh)
      interpolatePotential(lsms, local.atom[i]);
  }

  /*
   * 4. Calculate volumes
   */
  if (lsms.global.iprint >= 1) {
    std::cout << "Calculate volumes...." << std::endl;
  }
  calculateVolumes(comm, lsms, crystal, local);

  /*
   * 5. Initialize mixing
   */
  Mixing *mixing;
  setupMixing(mix, mixing, lsms.global.iprint);

  /*
   * 6. Initialize madelung
   */
  calculateMadelungMatrices(lsms, crystal, local);

  if (lsms.global.iprint >= 1) {
    printLocalTypeInfo(stdout, local);
  }


  /*
   * 7. Calculate core states
   */
  if (lsms.global.iprint >= 1) {
    std::cout << "Calculate core states...." << std::endl;
  }
  calculateCoreStates(comm, lsms, local);

#ifdef LSMS_CORE_DEBUG
  // 8. Print the core states
  printCoreOutput(lsms, crystal, local, comm);
#endif // LSMS_CORE_DEBUG

  prepareSpinOnAtoms(lsms, local);

  // 10. Prepare mixing
  mixing->prepare(comm, lsms, local.atom);

  printVolumes(lsms, comm, crystal, local);

  /*
   * SCF Cycle
   */

  Real eband{0.0};


  double timeScfLoop = MPI_Wtime();
  double timeCalcPotentialsAndMixing = 0.0;
  double timeCalcChemPot = 0.0;

  bool converged = false;
  int iteration;

  //
  if (comm.rank == 0) {
    if (lsms.n_spin_pola == 1) {
      std::printf("\n");
      std::printf("%4s %12s %16s %16s %16s %16s %10s\n",
                  "#.",
                  "Fermi",
                  "Band",
                  "Core+Inter",
                  "Site local",
                  "Total Energy",
                  "RMS");
    } else if (lsms.n_spin_pola == 2) {
      std::printf("\n");
      std::printf("%4s %12s %16s %16s %16s %16s %10s %10s\n",
                  "#.",
                  "Fermi",
                  "Band",
                  "Core+Inter",
                  "Site local",
                  "Total Energy",
                  "RMS",
                  "Mag.");
    }

  }


  // Reset
  local.tmatStore = 0.0;
  for (int i = 0; i < local.num_local; i++) {
    local.atom[i].resetLocalDensities();
  }
  //


  for (iteration = 0; iteration < lsms.nscf && !converged; iteration++) {


    // Calculate band energy
    energyContourIntegration(comm, lsms, local);

    double dTimeCCP = MPI_Wtime();

    // Calculate chemical potential
    calculateChemPot(comm, lsms, local, eband);
    dTimeCCP = MPI_Wtime() - dTimeCCP;
    timeCalcChemPot += dTimeCCP;


    // Calculate magnetic moments for each site and check if spin has flipped
    calculateEvec(lsms, local);
    mixEvec(lsms, local, 0.0);
    for (int i = 0; i < local.num_local; i++) {
      local.atom[i].newConstraint();
      checkIfSpinHasFlipped(lsms, local.atom[i]);
    }

    double dTimePM = MPI_Wtime();

    // Calculate charge densities, potentials, and total energy
    calculateAllLocalChargeDensities(lsms, local);
    calculateChargesPotential(comm, lsms, local, crystal, 0);
    checkAllLocalCharges(lsms, local);

    // Ensure that everything is calculated on the same potential
    calculateTotalEnergy(comm, lsms, local, crystal);

    // Calculate charge density rms
    calculateLocalQrms(lsms, local);

    // Mix charge density
    mixing->updateChargeDensity(comm, lsms, local.atom);
    dTimePM = MPI_Wtime() - dTimePM;
    timeCalcPotentialsAndMixing += dTimePM;


    // Recalculate core states
    // - swap core state energies for different spin channels first if spin has flipped
    //   (from LSMS 1: lsms_main.f:2101-2116)
    for (int i = 0; i < local.num_local; i++) {
      if (local.atom[i].spinFlipped) {
        checkIfSpinHasFlipped(lsms, local.atom[i]);
        if (!local.atom[i].spinFlipped)
          swapCoreStateEnergies(local.atom[i]);
      }
    }
    calculateCoreStates(comm, lsms, local);


    dTimePM = MPI_Wtime();
    // If charge is mixed, recalculate potential and mix (need a flag for this from input)
    calculateChargesPotential(comm, lsms, local, crystal, 1);
    mixing->updatePotential(comm, lsms, local.atom);
    dTimePM = MPI_Wtime() - dTimePM;
    timeCalcPotentialsAndMixing += dTimePM;

    Real rms = 0.0;
    for (int i = 0; i < local.num_local; i++) {
      Real averaged_rms = 0.0;
      for (int is = 0; is < lsms.n_spin_pola; is++) {
        averaged_rms += local.atom[i].qrms[is];
      }
      averaged_rms /= lsms.n_spin_pola;
      rms = std::max(rms, averaged_rms);
    }
    globalMax(comm, rms);

    if (iteration > MIN_ITERATIONS) {
      converged = rms < lsms.rmsTolerance;
    }

    auto sum_mag = 0.0;
    if (lsms.n_spin_pola == 2) {
      getMagneticMoments(lsms, comm, crystal, local, sum_mag);
    }
    auto e_core_inter = 0.0;
    getInterCoreEnergy(lsms, comm, crystal, local, e_core_inter);

    auto e_site = 0.0;
    getSiteLocalEnergy(lsms, comm, crystal, local, e_site);

    if (comm.rank == 0) {
      if (lsms.n_spin_pola == 1) {

        std::printf("%4d %12.8f %16.8f %16.8f %16.8f %16.8f %10.4e\n",
                    iteration,
                    lsms.chempot,
                    eband / lsms.num_atoms,
                    e_core_inter / lsms.num_atoms,
                    e_site / lsms.num_atoms,
                    lsms.totalEnergy / lsms.num_atoms,
                    rms);

      } else if (lsms.n_spin_pola == 2) {

        std::printf("%4d %12.8f %16.8f %16.8f %16.8f %16.8f %10.4e %12.6f\n",
                    iteration,
                    lsms.chempot,
                    eband / lsms.num_atoms,
                    e_core_inter / lsms.num_atoms,
                    e_site / lsms.num_atoms,
                    lsms.totalEnergy / lsms.num_atoms,
                    rms,
                    sum_mag);
      }
    }

    // Recalculate core states for new potential if we are performing scf calculations
    calculateCoreStates(comm, lsms, local);

  }

  timeScfLoop = MPI_Wtime() - timeScfLoop;

  if (lsms.global.iprint >= 0) {

    if (converged) {
      std::cout << std::endl;
      std::cout << "SCF is converged" << std::endl;
      std::cout << std::endl;
    } else {
      std::cout << std::endl;
      std::cout << "SCF is not converged" << std::endl;
      std::cout << std::endl;
    }

    std::printf("%s: %lf", "Time in SCF Loop", timeScfLoop);
  }


  delete mixing;

  /*
   * Print local energies
   */

  if (mode == EnergyDebugMode::Local) {

    /*
     *  Total energy per atom
     */

    std::vector<double> global_energies;

    lsms::gatherEnergies(comm, lsms, local, crystal, global_energies, 0);

    if (lsms.global.iprint >= 0) {

      std::cout << std::endl;

      for (int i = 0; i < global_energies.size(); i++) {
        std::printf("[%2d:%2d]     : %30.18lf Ry\n",
                    crystal.types[i].node,
                    crystal.types[i].local_id,
                    global_energies[i]);
      }

    }


    auto total_interstitial = 0.0;
    auto total_band_sum = 0.0;
    auto core_energies = 0.0;
    auto xc_energies = 0.0;
    auto coloumb_energies = 0.0;
    auto hartree_energies = 0.0;

    for (int i = 0; i < local.num_local; i++) {

      auto &atom = local.atom[i];

      total_interstitial += atom.interstitialEnergy;
      xc_energies += atom.energyStruct.exchange_correlation;
      hartree_energies += atom.energyStruct.hartree;
      coloumb_energies += atom.energyStruct.core_interaction;

      if (lsms.n_spin_pola == 1) {
        total_band_sum += atom.evalsum[0];
        core_energies += atom.ecorv[0] + atom.esemv[0];
      } else {
        total_band_sum += atom.evalsum[0] + atom.evalsum[1];
        core_energies += atom.ecorv[0] + atom.ecorv[1] + atom.esemv[0] + atom.esemv[1]; // (1)
      }
    }

    globalSum(comm, total_interstitial);
    globalSum(comm, total_band_sum);
    globalSum(comm, core_energies);
    globalSum(comm, xc_energies);
    globalSum(comm, coloumb_energies);
    globalSum(comm, hartree_energies);

    if (lsms.global.iprint >= 0) {

      std::printf(" Coloumb energies : %30.18lf Ry\n", coloumb_energies);
      std::printf(" Hartree energies : %30.18lf Ry\n", hartree_energies);
      std::printf(" XC energies      : %30.18lf Ry\n", xc_energies);
      std::printf(" Band sum         : %30.18lf Ry\n", total_band_sum);
      std::printf(" Core energies    : %30.18lf Ry\n", core_energies);
      std::printf(" Interstitial     : %30.18lf Ry\n", total_interstitial);
      std::printf(" Madelung         : %30.18lf Ry\n", lsms.u0);
      std::printf(" Total Energy     : %30.18lf Ry\n", lsms.totalEnergy);
      std::printf(" Total Energy/Atom: %30.18lf Ry\n", lsms.totalEnergy / lsms.num_atoms);

    }


  }


/*
 * Forces
 */

  lsms::HigherOrderMadelung higherOrderMadelung(local.num_local, crystal.num_atoms);

  lsms::calculateHigherOrderMadelung(lsms,
                                     crystal,
                                     local,
                                     higherOrderMadelung
  );

  lsms::calculateForces(comm,
                        lsms,
                        local,
                        crystal,
                        higherOrderMadelung,
                        lsms
                            .forceParams);

  if (lsms.global.iprint > 0) {
    lsms::displayForces(comm, lsms, local, crystal,
                        0);
  }


  lsms::normalizeForces(comm,
                        lsms,
                        local,
                        crystal
  );


  lsms::displayForces(comm, lsms, local, crystal,
                      0);

  if (comm.rank == 0) {
    std::printf("\nEnd DFT calculation:\n\n");
  }

  local.tmatStore.

      unpinMemory();

}




