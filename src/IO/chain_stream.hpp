//
// Created by F.Moitzi on 20.09.2021.
//

#ifndef LSMS_CHAIN_STREAM_HPP
#define LSMS_CHAIN_STREAM_HPP

#include <mpi.h>

#include "mpi_chain.hpp"

namespace lsms {

  class ChainStream : public MPIChain {
    /*
     * Uses the MPIChain class to implement a ostream with a serial operator<< implementation.
     */
  private:
    std::ostream &s_out;

  public:

    ChainStream(std::ostream &os,
                int c_rank,
                int c_size,
                MPI_Comm comm = MPI_COMM_WORLD)
    : MPIChain(0, c_rank, c_size, comm), s_out(os) {};

    ChainStream &operator<<(const std::string &os);;
  };

}

#endif //LSMS_CHAIN_STREAM_HPP
