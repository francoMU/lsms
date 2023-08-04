//
// Created by F.Moitzi on 30.07.2023.
//

#ifndef LSMS_SPGLIBWRAPPER_H
#define LSMS_SPGLIBWRAPPER_H

#include "spglib.h"
#include "Eigen/Dense"

// https://arne-mertz.de/2018/10/calling-cpp-code-from-c-with-extern-c/
#ifdef __cplusplus
extern "C" {
#include "spglib.h"
}
#endif

namespace SC {

    constexpr double SYMPREC_DEFAULT = 1.0e-5;

    /**
     * Class that wraps the sglib dataset from c
     */
    class SpglibWrapper {

    private:

        double symprec;

        SpglibDataset *dataset;

    public:

        SpglibWrapper(const Eigen::Matrix<double, 3, 3> &lattice,
                      const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordinates,
                      const Eigen::Vector<int, Eigen::Dynamic> &specie,
                      double symprec = SYMPREC_DEFAULT);

        SpglibWrapper(const SpglibWrapper &other);

        SpglibWrapper &operator=(const SpglibWrapper &rhs);

        SpglibWrapper(SpglibWrapper &&other) noexcept;

        SpglibWrapper &operator=(SpglibWrapper &&rhs);

        ~SpglibWrapper();


    };

} // SC

#endif //LSMS_SPGLIBWRAPPER_H
