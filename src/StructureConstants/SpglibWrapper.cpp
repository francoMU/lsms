//
// Created by F.Moitzi on 30.07.2023.
//

#include "SpglibWrapper.h"

namespace SC {


    SpglibWrapper::SpglibWrapper(const Eigen::Matrix<double, 3, 3> &lattice,
                                 const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordinates,
                                 const Eigen::Vector<int, Eigen::Dynamic> &specie,
                                 double symrec) {

        this->dataset = spg_get_dataset(lat, position, types, num_atom, symprec);


    }

    SpglibWrapper::SpglibWrapper(const SpglibWrapper& other) {

    }

    SpglibWrapper& SpglibWrapper::operator=(const SpglibWrapper& rhs) {

    }

    SpglibWrapper::SpglibWrapper(SpglibWrapper&& other) noexcept {

    }

    SpglibWrapper& SpglibWrapper::operator=(SpglibWrapper&& rhs) {

    }



    SpglibWrapper::~SpglibWrapper() {
        spg_free_dataset(dataset);
    }


} // SC