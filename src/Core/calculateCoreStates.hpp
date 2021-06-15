
#ifndef MUST_CALCULATECORESTATES_HPP
#define MUST_CALCULATECORESTATES_HPP

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

void calculateCoreStates(LSMSCommunication &comm, LSMSSystemParameters &lsms, LocalTypeInfo &local);

void calculateCoreStates(LSMSCommunication &comm, LSMSSystemParameters &lsms, AlloyAtomBank &alloyBank);

#endif //MUST_CALCULATECORESTATES_HPP
