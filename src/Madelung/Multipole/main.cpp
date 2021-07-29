

#include <vector>
#include <iostream>

#include "mpi.h"

int main(int argc, char *argv[]) {
    int rank;
    int np;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    /*
     *
     */

    std::vector<int> recv_buffer(np * 2, 0.0);
    std::vector<int> send_buffer(2, rank + 1);

    MPI_Allgather(
            send_buffer.data(),
            2,
            MPI_INTEGER,
            recv_buffer.data(),
            2,
            MPI_INTEGER,
            MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    for (auto &a : recv_buffer) {
        std::cout << "R: " << rank << ": " << a << std::endl;
    }

    /*
     *
     */

    MPI_Finalize();

    return 0;
}
