//
// Created by F.Moitzi on 15.07.2023.
//

#ifndef LSMS_SRC_STRUCTURECONSTANTS_STRUCTURE_H_
#define LSMS_SRC_STRUCTURECONSTANTS_STRUCTURE_H_

#include <Eigen/Dense>
#include <map>
#include<tuple>

namespace SC {


    constexpr double DEFAULT_DELTA = 0.75;

    class Structure {


    private:


        // Lattice
        Eigen::Matrix<double, 3, 3> lattice;

        // Coordinates
        Eigen::Matrix<double, Eigen::Dynamic, 3> coordinates;

        // Species
        Eigen::Vector<int, Eigen::Dynamic> species;

        // Map for lll-reduction
        std::map<double, std::tuple<Eigen::Matrix<double, 3, 3>, Eigen::Matrix<double, 3, 3>>> map;

    public:

        Structure(const Eigen::Matrix<double, 3, 3> &lattice,
                  const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordinates,
                  const Eigen::Vector<int, Eigen::Dynamic> &species);


        void calculateLLLreduction(
                Eigen::Matrix<double, 3, 3> &lll,
                Eigen::Matrix<double, 3, 3> &lll_mapping,
                double delta = DEFAULT_DELTA);


        const Eigen::Matrix<double, 3, 3> &lll_matrix();

        const Eigen::Matrix<double, 3, 3> &lll_mapping();

        void get_distances(const Eigen::Vector<double, 3> &fcoords1, const Eigen::Vector<double, 3> &fcoords2);

    };


}

#endif //LSMS_SRC_STRUCTURECONSTANTS_STRUCTURE_H_
