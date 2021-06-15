
#ifndef MUST_SETUPVORPOL_HPP
#define MUST_SETUPVORPOL_HPP

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"
#include "Misc/Coeficients.hpp"
#include "VORPOL.hpp"

void setupVorpol(LSMSSystemParameters &lsms,
                 CrystalParameters &crystal,
                 LocalTypeInfo &local,
                 SphericalHarmonicsCoeficients &shc);

void calculateVolumes(LSMSCommunication &comm,
                      LSMSSystemParameters &lsms,
                      CrystalParameters &crystal,
                      LocalTypeInfo &local);


#endif //MUST_SETUPVORPOL_HPP
