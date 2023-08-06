//
// Created by F.Moitzi on 30.07.2023.
//

#ifndef LSMS_SPGLIBWRAPPER_H
#define LSMS_SPGLIBWRAPPER_H

#include "spglib.h"
#include "Eigen/Dense"

#include <boost/algorithm/string.hpp>

// https://arne-mertz.de/2018/10/calling-cpp-code-from-c-with-extern-c/
#ifdef __cplusplus
extern "C" {
#include "spglib.h"
}
#endif

namespace SC {

    namespace ba = boost::algorithm;

    constexpr double SYMPREC_DEFAULT = 1.0e-5;

   /**
    * @class SpglibWrapper
    * @brief A class that wraps functionality from the Spglib library for crystal symmetry analysis.
    */
    class SpglibWrapper {

    private:

        double symprec{};

        SpglibDataset *dataset;

        /**
            * @brief Resizes a SpglibDataset.
            * @param rhs_dataset The source dataset to be resized from.
            * @param n_ops Number of symmetry operations.
            * @param n_atoms Number of atoms.
            * @param n_std_atoms Number of standard atoms.
            */
        static void resize(SpglibDataset *rhs_dataset, int n_ops, int n_atoms, int n_std_atoms);

        /**
         * @brief Copies a SpglibDataset.
         * @param rhs_dataset The source dataset to be copied from.
         * @param copy_dataset The destination dataset to copy to.
         */
        static void copy(const SpglibDataset *rhs_dataset, SpglibDataset *copy_dataset);


    public:

        /**
         * @brief Constructs a SpglibWrapper instance using lattice parameters, coordinates, and species information.
         * @param lattice The lattice parameters as a 3x3 matrix.
         * @param coordinates The atomic coordinates as a matrix with rows containing (x, y, z) coordinates.
         * @param specie The species information as a vector of integers.
         * @param symprec The symmetry precision. Default is SYMPREC_DEFAULT.
         */
        SpglibWrapper(const Eigen::Matrix<double, 3, 3> &lattice,
                      const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordinates,
                      const Eigen::Vector<int, Eigen::Dynamic> &specie,
                      double symprec = SYMPREC_DEFAULT);

        /**
         * @brief Copy constructor for SpglibWrapper.
         * @param other The source object to be copied from.
         */
        SpglibWrapper(const SpglibWrapper &other);

        /**
         * @brief Copy assignment operator for SpglibWrapper.
         * @param rhs The source object to be copied from.
         * @return Reference to the modified object after the copy.
         */
        SpglibWrapper &operator=(const SpglibWrapper &rhs);

        /**
         * @brief Move constructor for SpglibWrapper.
         * @param other The source object to be moved from.
         */
        SpglibWrapper(SpglibWrapper &&other) noexcept;

        /**
         * @brief Move assignment operator for SpglibWrapper.
         * @param rhs The source object to be moved from.
         * @return Reference to the modified object after the move.
         *
         * `noexcept` serves compiler optimizations similarly to const,
         * though rarely. It primarily helps you detect compile-time exceptions.
         * It doesn't dictate optimization, but guides usage. move_if_noexcept
         * template returns const& instead of && if not noexcept, suggesting safe moves.
         *
         */
        SpglibWrapper &operator=(SpglibWrapper &&rhs) noexcept;

        /**
         * @brief Returns the international symbol representing the crystal symmetry.
         * @return The international symbol as a string.
         */
        std::string international_symbol();

        /**
         * @brief Destructor for SpglibWrapper.
         */
        ~SpglibWrapper();


    };

} // SC

#endif //LSMS_SPGLIBWRAPPER_H
