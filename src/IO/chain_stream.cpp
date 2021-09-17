//
// Created by F.Moitzi on 20.09.2021.
//


#include "chain_stream.hpp"

namespace lsms {

  ChainStream &ChainStream::operator<<(const std::string &os) {

    MPI_Barrier(this->get_comm());

    if (this->get_rank() == 0) {
      this->s_out << os << std::flush;
      // Initiate chain of MPI messages
      this->next();
    } else {

      int msg_count = 0;

      // Wait until a message arrives (MPIChain::wait uses a blocking test)
      this->wait(msg_count);
      if (msg_count == 1) {
        // If the message is well-formed (i.e. only one message is received): output string
        this->s_out << os << std::flush;
        // Pass onto the next member of the chain (if there is one)
        this->next();
      }

    }

    // Ensure that the chain is resolved before returning the stream
    MPI_Barrier(this->get_comm());

    // Don't output the ostream! That would break the serial-in-time exeuction.
    return *this;
  }
}

