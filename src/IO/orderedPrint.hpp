//
// Created by F.Moitzi on 17.09.2021.
//

#ifndef LSMS_ORDEREDPRINT_HPP
#define LSMS_ORDEREDPRINT_HPP

#include <mpi.h>

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

template <std::size_t Size = 128>
void orderedPrint(const std::string & message_string,
                  LSMSCommunication &comm) {

   char recv_buffer[Size];

   if (comm.rank == 0) {
      std::cout << "RANK: " << comm.rank << " " << message_string << std::endl;

      // Send message to next process
      MPI_Send(message_string.c_str(), Size, MPI_CHAR, 1, 0, comm.comm);

   } else {

      int counter;
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, 0, comm.comm, &status);
      MPI_Get_count(&status, MPI_CHAR, &counter);

      auto rank = comm.rank;

      if (counter >= 0) {

         std::cout << "RANK: " << comm.rank << " " << message_string << std::endl;

         MPI_Recv(recv_buffer, counter, MPI_CHAR, MPI_ANY_SOURCE, 0, comm.comm, &status);

         if (rank + 1 != comm.size) {
            MPI_Send(message_string.c_str(), Size, MPI_CHAR, ++rank, 0, comm.comm);
         }

      };
   };


}


#endif //LSMS_ORDEREDPRINT_HPP
