
#ifndef LSMS_UTILS_HPP
#define LSMS_UTILS_HPP


#include "SystemParameters.hpp"
#include "LSMSCommunication.hpp"

namespace lsms {

    void printLocalTypeInfo(const LSMSCommunication &comm,
                            const LocalTypeInfo &local,
                            int rank = 0);

    void printKKRMatrix(Matrix<Complex> &kkr_matrix,
                        std::string name,
                        int rank = 0);

}


#endif //LSMS_UTILS_HPP
