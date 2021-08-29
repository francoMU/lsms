
#ifndef MUST_CALCHIGHERORDERMADELUNG_HPP
#define MUST_CALCHIGHERORDERMADELUNG_HPP

#include "higherOrderMadelung.hpp"
#include "Main/SystemParameters.hpp"


namespace lsms {

    void calculateHigherOrderMadelung(
            LSMSSystemParameters &lsms,
            CrystalParameters &crystal,
            LocalTypeInfo &local,
            HigherOrderMadelung &madelung);

}


#endif //MUST_CALCHIGHERORDERMADELUNG_HPP
