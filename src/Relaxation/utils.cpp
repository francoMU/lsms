
#include "utils.hpp"

#include <icecream.hpp>

void lsms::printLocalTypeInfo(const LSMSCommunication &comm,
                              const LocalTypeInfo &local,
                              int rank) {

   std::printf("RANK: %d\n", rank);

   if (comm.rank == rank) {

      IC(local.num_local);

      IC(local.global_id);

      IC(local.n_per_type);

      IC(local.lDimTmatStore);
      IC(local.blkSizeTmatStore);

      IC(local.tmatStoreGlobalIdx);

      IC(local.tmatStore.n_row());
      IC(local.tmatStore.n_col());

      auto maxNrmat = local.maxNrmat();
      IC(maxNrmat);

      auto maxjws = local.maxjws();
      IC(maxjws);

   }

   std::printf("----------------\n");

}

void lsms::printKKRMatrix(Matrix<Complex> &kkr_matrix,
                          std::string name,
                          int rank) {

   int n_col = kkr_matrix.n_col();
   int n_row = kkr_matrix.n_row();

   MPI_Comm comm = MPI_COMM_WORLD;

   FILE *stream;
   stream = fopen(name.c_str(), "w");

   for (int i = 0; i < n_row; i++) {

      for (int j = 0; j < n_col; j++) {

         std::fprintf(stream, "  %18.10f %18.10f  ",
                      std::real(kkr_matrix(i, j)),
                      std::imag(kkr_matrix(i, j)));

      }

      std::fprintf(stream, "\n");

   }

   fclose(stream);

}
