//
// Created by F.Moitzi on 20.09.2021.
//

#ifndef LSMS_MPI_CHAIN_HPP
#define LSMS_MPI_CHAIN_HPP

#include <mpi.h>

#include "Real.hpp"

namespace lsms {

  constexpr int MPI_CHAIN_TAG = 233423;

  class MPIChain {

    /**
     * Uses a chained MPI message (T) to coordinate
     * serial execution of code (the content of the message is irrelevant).
     */

  private:
    int message_out; // The messages aren't really used here
    int message_in {0};
    int size;
    int rank;
    MPI_Comm comm;

  public:

    MPIChain(int message_init,
             int c_rank,
             int c_size,
             MPI_Comm comm) :
        message_out(message_init),
        size(c_size),
        rank(c_rank),
        comm(comm) {}

    void next() {
      // Send message to next core (if there is one)
      if (rank + 1 < size) {
        // MPI_Send - Performs a standard-mode blocking send.
        MPI_Send(&message_out, 1,
                 MPI_INT,
                 rank + 1, MPI_CHAIN_TAG, comm);
      }
    }

    void wait(int &msg_count) {
      // Waits for message to arrive. Message is well-formed if msg_count = 1
      MPI_Status status;

      // MPI_Probe - Blocking test for a message.
      MPI_Probe(MPI_ANY_SOURCE, MPI_CHAIN_TAG, comm, &status);
      // MPI_Get_count - Gets the number of top level elements.
      MPI_Get_count(&status, MPI_INT, &msg_count);

      if (msg_count == 1) {
        // MPI_Recv - Performs a standard-mode blocking receive.
        MPI_Recv(&message_in, msg_count, MPI_INT, MPI_ANY_SOURCE, MPI_CHAIN_TAG, comm, &status);
      }
    }


    int get_rank() const { return rank; }

    int get_size() const { return size; }

    MPI_Comm get_comm() const { return comm; }

  };
}

#endif //LSMS_MPI_CHAIN_HPP
