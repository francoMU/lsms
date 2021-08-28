
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>


#ifdef _OPENMP

#include <omp.h>

#endif

#ifdef USE_PAPI
#include <papi.h>
#endif

#ifdef USE_GPTL
#include "gptl.h"
#endif

#include <hdf5.h>
#include <POSCARStructureIO.hpp>
#include <fstream>

#include "lua.hpp"

#include "SystemParameters.hpp"
#include "PotentialIO.hpp"
#include "Communication/distributeAtoms.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Core/calculateCoreStates.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "Madelung/Madelung.hpp"
#include "EnergyContourIntegration.hpp"
#include "Accelerator/Accelerator.hpp"
#include "calculateChemPot.hpp"
#include "calculateDensities.hpp"
#include "mixing.hpp"
#include "calculateEvec.hpp"
#include "Potential/calculateChargesPotential.hpp"
#include "Potential/interpolatePotential.hpp"
#include "Potential/PotentialShifter.hpp"
#include "TotalEnergy/calculateTotalEnergy.hpp"
#include "SingleSite/checkAntiFerromagneticStatus.hpp"
#include "Forces/forces.hpp"
#include "Madelung/Multipole/higherOrderMadelung.hpp"
#include "Madelung/Multipole/calcHigherOrderMadelung.hpp"
#include "Main/read_input.h"
#include "LuaInterface/LuaInterface.hpp"
#include "Misc/readLastLine.hpp"
#include "VORPOL/setupVorpol.hpp"
#include "VORPOL/VORPOL.hpp"
#include "buildLIZandCommLists.hpp"


#include "writeInfoEvec.cpp"
#include "write_restart.hpp"
#include "mixing_params.hpp"

#include "lsms_relaxation.hpp"
#include "lsms_relaxation_function.hpp"
#include "nl_optimization.hpp"

SphericalHarmonicsCoeficients sphericalHarmonicsCoeficients;
GauntCoeficients gauntCoeficients;
IFactors iFactors;

#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
#include "Accelerator/DeviceStorage.hpp"
// void * deviceStorage;
DeviceStorage *deviceStorage;
DeviceConstants deviceConstants;
#endif
#ifdef BUILDKKRMATRIX_GPU
#include "Accelerator/buildKKRMatrix_gpu.hpp"

std::vector<DeviceConstants> deviceConstants;
// void *allocateDConst(void);
// void freeDConst(void *);
#endif


#define LSMS_RELAX_DEBUG

int main(int argc, char *argv[]) {

   LSMSSystemParameters lsms;
   LSMSCommunication comm;
   CrystalParameters crystal;
   LocalTypeInfo local;
   MixingParameters mix;
   PotentialShifter potentialShifter;
   AlloyMixingDesc alloyDesc;
   AlloyAtomBank alloyBank;

   char inputFileName[128];

   Real eband;

   lua_State *L = luaL_newstate();
   luaL_openlibs(L);
   initLSMSLuaInterface(L);

#ifdef USE_GPTL
   GPTLinitialize();
#endif
   initializeCommunication(comm);
   H5open();

   // set input file name (default 'i_lsms')
   strncpy(inputFileName, "i_lsms", 10);
   if (argc > 1)
      strncpy(inputFileName, argv[1], 120);

   lsms.global.iprpts = 1051;
   lsms.global.ipcore = 30;
   lsms.global.setIstop("main");
   lsms.global.iprint = 0;
   lsms.global.default_iprint = -1;
   lsms.global.print_node = 0;
   lsms.ngaussr = 10;
   lsms.ngaussq = 40;
   lsms.vSpinShiftFlag = 0;
#ifdef _OPENMP
   lsms.global.GPUThreads = std::min(12, omp_get_max_threads());
#else
   lsms.global.GPUThreads = 1;
#endif
   if (comm.rank == 0) {
      lsms.global.iprint = 0;
   }

   if (comm.rank == 0) {
      printf("LSMS_3: Program started\n");
      printf("Using %d MPI processes\n", comm.size);
#ifdef _OPENMP
      printf("Using %d OpenMP threads\n", omp_get_max_threads());
#endif
      acceleratorPrint();
#ifdef BUILDKKRMATRIX_GPU
      printf("Using GPU to build KKR matrix.\n");
#endif
#ifdef LSMS_NO_COLLECTIVES
      printf("\nWARNING!!!\nCOLLECTIVE COMMUNICATION (ALLREDUCE etc.) ARE SKIPPED!\n");
      printf("THIS IS FOR TESTING ONLY!\nRESULTS WILL BE WRONG!!!\n\n");
#endif
      printf("Reading input file '%s'\n", inputFileName);
      fflush(stdout);

      if (luaL_loadfile(L, inputFileName) || lua_pcall(L, 0, 0, 0)) {
         fprintf(stderr, "!! Cannot run input file!!\n");
         exit(1);
      }

      printf("Loaded input file!\n");
      fflush(stdout);

      if (readInput(L, lsms, crystal, mix, potentialShifter, alloyDesc)) {
         fprintf(stderr, "!! Something wrong in input file!!\n");
         exit(1);
      }

      printf("System information:\n");
      printf("===================\n");
      printf("Number of atoms        : %10d\n", crystal.num_atoms);
      printf("Number of atomic types : %10d\n", crystal.num_types);
      switch (lsms.mtasa) {
         case 1:
            printf("Performing Atomic Sphere Approximation (ASA) calculation\n");
            break;
         case 2:
            printf("Performing Atomic Sphere Approximation + Muffin-Tin (ASA-MT) calculation\n");
            break;
         default:
            printf("Performing Muffin-Tin (MT) calculation\n");
      }
      fflush(stdout);
   }


   // Parameters are communicated here
   communicateParameters(comm, lsms, crystal, mix, alloyDesc);

   if (comm.rank != lsms.global.print_node) {
      lsms.global.iprint = lsms.global.default_iprint;
   }

   // Setup up expansion coeffiencts
   lsms.angularMomentumIndices.init(2 * crystal.maxlmax);
   sphericalHarmonicsCoeficients.init(2 * crystal.maxlmax);
   gauntCoeficients.init(lsms, lsms.angularMomentumIndices, sphericalHarmonicsCoeficients);
   iFactors.init(lsms, crystal.maxlmax);


   lsms::LsmsRelaxationFunction relax_function(lsms, comm, crystal, local,
                                               mix, lsms.relaxParams.write_to_file);

   /*
    */

   std::vector<double> coordinates(crystal.num_atoms * 3);
   std::vector<double> gradient(crystal.num_atoms * 3);

   lsms::CoordinatesToVector(lsms,
                             comm,
                             crystal,
                             coordinates);

   auto starting_coordinates = coordinates;


   /*
    *
    */

   if (lsms.global.iprint >= 0) {
      std::ofstream start_file;
      start_file.open("start.POSCAR");
      lsms::POSCARStructureIO io(lsms::POSCARStructureType::Cartesian);
      io.writeToStream(start_file, lsms, crystal);
      start_file.close();
   }

   /*
    *
    */

   lsms::NLOptimization relaxation(relax_function,
                                   starting_coordinates,
                                   lsms.relaxParams.max_iterations,
                                   lsms.relaxParams.tolerance,
                                   lsms.relaxParams.initial_sigma);


   auto x_0 = coordinates;
   auto x_1 = coordinates;
   auto grad_0 = gradient;
   auto grad_1 = gradient;

   int iters;

   bool converged = false;

   if (lsms.global.iprint >= 0) {
      std::cout << " Start of iterations: " << std::endl;
   }

   relaxation.start(x_0, x_1, grad_0, grad_1);

   // Update all internal convergence parameters
   relaxation.update_step(x_0, x_1);

   // Update the gradient to save one evaluation per step
   grad_0 = grad_1;

   for (iters = 0; iters <= lsms.relaxParams.max_iterations; iters++) {

      if (lsms.global.iprint >= 0) {

         std::cout << std::endl;
         std::cout << std::endl;
         std::cout << " ------------------------------ " << std::endl;
         std::cout << "           Iterations: " << iters << std::endl;
         std::cout << " ------------------------------ " << std::endl;
         std::cout << std::endl;

      }

      relaxation.iteration(x_0, x_1, grad_0, grad_1);

      if (relaxation.check_convergence(grad_1)) {
         converged = true;
         break;
      }

      grad_0 = grad_1;
   }

   if (lsms.global.iprint >= 0) {
      if (converged) {
         std::cout << " Optimization is converged " << std::endl;
      } else {
         std::cout << " Optimization is not converged " << std::endl;
      }
   }

   // Update last step to coordinates in crystal
   lsms::VectorToCoordinates(lsms,
                             comm,
                             crystal,
                             local,
                             x_1);


   if (lsms.global.iprint >= 0) {

      std::ofstream end_file;
      end_file.open("end.POSCAR");
      lsms::POSCARStructureIO io(lsms::POSCARStructureType::Cartesian);
      io.writeToStream(end_file, lsms, crystal);
      end_file.close();

      printCompressedCrystalParameters(stdout, crystal);
      relaxation.printRelaxationHistory();
   }

   if (lsms.pot_out_type >= 0) {
      if (comm.rank == 0) std::cout << "Writing new potentials.\n";
      writePotentials(comm, lsms, crystal, local);
      if (comm.rank == 0) {
         std::cout << "Writing restart file.\n";
         writeRestart("i_lsms.restart", lsms, crystal, mix, potentialShifter, alloyDesc);
      }
   }

   if (lsms.pot_out_type >= 0) {
      if (comm.rank == 0) {
         std::cout << "Simulation is done!!!!" << std::endl;
      }
   }

#ifdef BUILDKKRMATRIX_GPU
   // for(int i=0; i<local.num_local; i++) freeDConst(deviceConstants[i]);
#endif

#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C) || defined(ACCELERATOR_HIP)
   // freeDStore(deviceStorage);
   delete deviceStorage;
#endif
#ifdef BUILDKKRMATRIX_GPU
   deviceConstants.clear();
#endif

   acceleratorFinalize();

#ifdef USE_GPTL
   GPTLpr(comm.rank);
#endif

   H5close();
   finalizeCommunication();
   lua_close(L);
   return 0;
}
