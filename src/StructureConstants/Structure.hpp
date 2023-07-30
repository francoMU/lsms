/**
 * @file
 * @brief This file was authored by F. Moitzi.
 */

#ifndef LSMS_SRC_STRUCTURECONSTANTS_STRUCTURE_H_
#define LSMS_SRC_STRUCTURECONSTANTS_STRUCTURE_H_

#include <Eigen/Dense>
#include <map>
#include<tuple>

namespace SC {

    constexpr double DEFAULT_DELTA = 0.75;

    /**
    * @class Structure
    * @brief Class representing a crystal structure.
    */
    class Structure {


    private:

        Eigen::Matrix<double, 3, 3> lattice;      /**< The lattice matrix. */
        Eigen::Matrix<double, 3, 3> inv_lattice;  /**< The inverse of the lattice matrix. */

        Eigen::Matrix<double, Eigen::Dynamic, 3> coordinates;  /**< The fractional coordinates of atoms. */

        Eigen::Vector<int, Eigen::Dynamic> species;  /**< The species of each atom in the structure. */

        std::map<double, std::tuple<Eigen::Matrix<double, 3, 3>, Eigen::Matrix<double, 3, 3>>> map; /**< A map for LLL-reduction. */

    public:

        /**
         * @brief Constructor for the Structure class.
         * @param lattice The lattice matrix.
         * @param coordinates The fractional coordinates of atoms.
         * @param species The species of each atom in the structure.
         */
        Structure(const Eigen::Matrix<double, 3, 3> &lattice,
                  const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordinates,
                  const Eigen::Vector<int, Eigen::Dynamic> &species);

        /**
         * @brief Calculate the LLL-reduction of the lattice matrix.
         * @param lll The LLL-reduced lattice matrix (output parameter).
         * @param lll_mapping The transformation matrix for LLL-reduction (output parameter).
         * @param delta The tolerance for LLL-reduction (optional, default is DEFAULT_DELTA).
         */
        void calculateLLLreduction(
                Eigen::Matrix<double, 3, 3> &lll,
                Eigen::Matrix<double, 3, 3> &lll_mapping,
                double delta = DEFAULT_DELTA);

        /**
         * @brief Get the LLL-reduced lattice matrix.
         * @return The LLL-reduced lattice matrix.
         */
        const Eigen::Matrix<double, 3, 3> &lll_matrix();

        /**
         * @brief Get the transformation matrix for LLL-reduction.
         * @return The LLL-reduction transformation matrix.
         */
        const Eigen::Matrix<double, 3, 3> &lll_mapping();

        /**
         * @brief Calculate the distances between two fractional coordinates.
         * @param fcoords1 The first fractional coordinates.
         * @param fcoords2 The second fractional coordinates.
         * @param dist_vec The vector representing the distance between the coordinates (output parameter).
         * @param dist The distance between the coordinates (output parameter).
         */
        void get_distances(const Eigen::Vector<double, 3> &fcoords1,
                           const Eigen::Vector<double, 3> &fcoords2,
                           Eigen::Vector<double, 3> &dist_vec,
                           double &dist);

        /**
         * @brief Calculate the distances between two Cartesian coordinates.
         * @param coords1 The first Cartesian coordinates.
         * @param coords2 The second Cartesian coordinates.
         * @param dist_vec The vector representing the distance between the coordinates (output parameter).
         * @param dist The distance between the coordinates (output parameter).
         */
        void get_distances_cart(const Eigen::Vector<double, 3> &coords1,
                                const Eigen::Vector<double, 3> &coords2,
                                Eigen::Vector<double, 3> &dist_vec,
                                double &dist);

        /**
         * @brief Convert fractional coordinates to Cartesian coordinates.
         * @param frac_coords The fractional coordinates to convert.
         * @param cart_coords The resulting Cartesian coordinates (output parameter).
         */
        void get_fractional_coords(
                const Eigen::Matrix<double, Eigen::Dynamic, 3> &frac_coords,
                Eigen::Matrix<double, Eigen::Dynamic, 3> &cart_coords
        );


        /**
         * @brief Convert Cartesian coordinates to fractional coordinates.
         * @param frac_coords The Cartesian coordinates to convert.
         * @param cart_coords The resulting fractional coordinates (output parameter).
         */
        void get_cartesian_coords(
                const Eigen::Matrix<double, Eigen::Dynamic, 3> &frac_coords,
                Eigen::Matrix<double, Eigen::Dynamic, 3> &cart_coords
        );


    };


}

#endif //LSMS_SRC_STRUCTURECONSTANTS_STRUCTURE_H_
