
#ifndef LSMS_BUILDLIZANDCOMMLISTS_HPP
#define LSMS_BUILDLIZANDCOMMLISTS_HPP

#include "SystemParameters.hpp"
#include "AtomData.hpp"
#include "Communication/LSMSCommunication.hpp"


void buildLIZandCommLists(LSMSCommunication &comm,
                          LSMSSystemParameters &lsms,
                          CrystalParameters &crystal,
                          LocalTypeInfo &local);


#endif
