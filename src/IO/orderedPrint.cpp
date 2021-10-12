//
// Created by F.Moitzi on 18.09.2021.
//

#include "orderedPrint.hpp"

#include "Forces/forces.hpp"


void lsms::gatherInfo(const LSMSCommunication &comm,
                      const LSMSSystemParameters &lsms,
                      const LocalTypeInfo &local,
                      const CrystalParameters &crystal,
                      std::vector<double> &rmt,
                      std::vector<double> &rws,
                      std::vector<double> &rInscribed,
                      int root) {
  //  Number of local atoms
  auto num_local_atoms = local.num_local;
  // Gives the total number of atoms
  auto num_atoms = lsms.num_atoms;

  MPI_Datatype mpi_double_type;
  auto status = MPI_Type_contiguous(1, MPI_DOUBLE, &mpi_double_type);
  MPI_Type_commit(&mpi_double_type);

  std::vector<double> local_rmt(num_local_atoms);
  std::vector<double> local_rws(num_local_atoms);
  std::vector<double> local_rInscribed(num_local_atoms);

  rmt.resize(num_atoms);
  rws.resize(num_atoms);
  rInscribed.resize(num_atoms);

  for (int i_local = 0; i_local < num_local_atoms; i_local++) {
    local_rmt[i_local] = local.atom[i_local].rmt;
    local_rws[i_local] = local.atom[i_local].rws;
    local_rInscribed[i_local] = local.atom[i_local].rInscribed;
  }

  // number of elements in the message to receive from each process
  std::vector<int> rcounts(comm.size);
  std::vector<int> displs(comm.size, 0);

  gatherStructure(comm, lsms, crystal,
                  rcounts, displs);

  status = MPI_Gatherv(local_rmt.data(),
                       num_local_atoms,
                       mpi_double_type,
                       rmt.data(),
                       rcounts.data(),
                       displs.data(),
                       mpi_double_type,
                       root,
                       comm.comm);

  status = MPI_Gatherv(local_rws.data(),
                       num_local_atoms,
                       mpi_double_type,
                       rws.data(),
                       rcounts.data(),
                       displs.data(),
                       mpi_double_type,
                       root,
                       comm.comm);

  status = MPI_Gatherv(local_rInscribed.data(),
                       num_local_atoms,
                       mpi_double_type,
                       rInscribed.data(),
                       rcounts.data(),
                       displs.data(),
                       mpi_double_type,
                       root,
                       comm.comm);

  status = MPI_Type_free(&mpi_double_type);

}

void lsms::gatherEnergies(const LSMSCommunication &comm,
                          const LSMSSystemParameters &lsms,
                          const LocalTypeInfo &local,
                          const CrystalParameters &crystal,
                          std::vector<double> &global_energies,
                          int root) {

  //  Number of local atoms
  auto num_local_atoms = local.num_local;
  // Gives the total number of atoms
  auto num_atoms = lsms.num_atoms;

  MPI_Datatype mpi_double_type;
  auto status = MPI_Type_contiguous(1, MPI_DOUBLE, &mpi_double_type);
  MPI_Type_commit(&mpi_double_type);

  std::vector<double> local_energies(num_local_atoms);
  global_energies.resize(num_atoms);

  for (int i_local = 0; i_local < num_local_atoms; i_local++) {
    local_energies[i_local] = local.atom[i_local].localEnergy;
  }

  // number of elements in the message to receive from each process
  std::vector<int> rcounts(comm.size);
  std::vector<int> displs(comm.size, 0);

  gatherStructure(comm, lsms, crystal,
                  rcounts, displs);

  status = MPI_Gatherv(local_energies.data(),
                       num_local_atoms,
                       mpi_double_type,
                       global_energies.data(),
                       rcounts.data(),
                       displs.data(),
                       mpi_double_type,
                       root,
                       comm.comm);


  status = MPI_Type_free(&mpi_double_type);

}