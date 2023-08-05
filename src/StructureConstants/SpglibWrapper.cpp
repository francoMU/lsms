//
// Created by F.Moitzi on 30.07.2023.
//

#include "SpglibWrapper.h"

namespace SC {


    SpglibWrapper::SpglibWrapper(const Eigen::Matrix<double, 3, 3> &lattice,
                                 const Eigen::Matrix<double, Eigen::Dynamic, 3> &frac_coordinates,
                                 const Eigen::Vector<int, Eigen::Dynamic> &species,
                                 double symprec) : symprec{symprec} {


        auto num_atom = frac_coordinates.rows();


        /**
         * 1.
         *  float *p[2] defines an array of float*, with size 2.
         *  float (*p)[2] defines a pointer to an array of float, with size 2:
         *
         * 2. The entry of the array will always decay
         */

        this->dataset = spg_get_dataset(reinterpret_cast<const double (*)[3]>(lattice.data()),
                                        reinterpret_cast<const double (*)[3]>(frac_coordinates.data()),
                                        reinterpret_cast<const int *>(species.data()),
                                        num_atom, symprec);


    }

    SpglibWrapper::SpglibWrapper(const SpglibWrapper &other) {


        this->dataset = new SpglibDataset();
        this->symprec = other.symprec;

        deepCopy(this->dataset, other.dataset);


    }

//    SpglibWrapper& SpglibWrapper::operator=(const SpglibWrapper& rhs) {
//
//    }
//
//    SpglibWrapper::SpglibWrapper(SpglibWrapper&& other) noexcept {
//
//    }
//
//    SpglibWrapper& SpglibWrapper::operator=(SpglibWrapper&& rhs) {
//
//    }
//


    SpglibWrapper::~SpglibWrapper() {
        spg_free_dataset(dataset);
    }

    void SpglibWrapper::deepCopy(const SpglibDataset *rhs_dataset, SpglibDataset *copy_dataset) {

        copy_dataset->spacegroup_number = rhs_dataset->spacegroup_number;
        copy_dataset->hall_number = rhs_dataset->hall_number;

        std::copy(std::begin(rhs_dataset->international_symbol),
                  std::end(rhs_dataset->international_symbol),
                  std::begin(copy_dataset->international_symbol));
        std::copy(std::begin(rhs_dataset->hall_symbol),
                  std::end(rhs_dataset->hall_symbol),
                  std::begin(copy_dataset->hall_symbol));
        std::copy(std::begin(rhs_dataset->choice),
                  std::end(rhs_dataset->choice),
                  std::begin(copy_dataset->choice));

        std::copy_n(&rhs_dataset->transformation_matrix[0][0], 9, &copy_dataset->transformation_matrix[0][0]);

        std::copy(std::begin(rhs_dataset->origin_shift),
                  std::end(rhs_dataset->origin_shift),
                  std::begin(copy_dataset->origin_shift));

        /* Operations */
        auto n_ops = rhs_dataset->n_operations;
        copy_dataset->n_operations = n_ops;

        copy_dataset->rotations = new int[n_ops][3][3];
        copy_dataset->translations = new double[n_ops][3];

        std::copy_n(&rhs_dataset->rotations[0][0][0], n_ops * 9, &copy_dataset->rotations[0][0][0]);
        std::copy_n(&rhs_dataset->translations[0][0], n_ops * 3, &copy_dataset->translations[0][0]);


        /* Atom related quantity */
        auto n_atoms = rhs_dataset->n_atoms;
        copy_dataset->n_atoms = n_atoms;

        copy_dataset->wyckoffs = new int[n_atoms];
        std::copy_n(rhs_dataset->wyckoffs,
                    n_atoms,
                    copy_dataset->wyckoffs);

        copy_dataset->equivalent_atoms = new int[n_atoms];
        std::copy_n(rhs_dataset->equivalent_atoms,
                    n_atoms,
                    copy_dataset->equivalent_atoms);

        copy_dataset->crystallographic_orbits = new int[n_atoms];
        std::copy_n(rhs_dataset->crystallographic_orbits,
                    n_atoms,
                    copy_dataset->crystallographic_orbits);

        std::copy_n(&rhs_dataset->primitive_lattice[0][0], 9, &copy_dataset->primitive_lattice[0][0]);

        copy_dataset->mapping_to_primitive = new int[n_atoms];
        std::copy_n(rhs_dataset->mapping_to_primitive,
                    n_atoms,
                    copy_dataset->mapping_to_primitive);

        /* Standardized lattice operation */
        auto n_std_atoms = rhs_dataset->n_std_atoms;

        copy_dataset->n_std_atoms = n_std_atoms;
        copy_dataset->std_types = new int[n_std_atoms];
        copy_dataset->std_mapping_to_primitive = new int[n_std_atoms];

        std::copy_n(&rhs_dataset->std_lattice[0][0], 9, &copy_dataset->std_lattice[0][0]);
        std::copy_n(rhs_dataset->std_types, n_std_atoms, copy_dataset->std_types);
        std::copy_n(&rhs_dataset->std_positions[0][0], n_std_atoms * 3, &copy_dataset->std_positions[0][0]);
        std::copy_n(&rhs_dataset->std_rotation_matrix[0][0], 9, &copy_dataset->std_rotation_matrix[0][0]);
        std::copy_n(rhs_dataset->std_mapping_to_primitive, n_std_atoms, copy_dataset->std_mapping_to_primitive);

        /* Point group symbol */
        std::copy(std::begin(rhs_dataset->pointgroup_symbol),
                  std::end(rhs_dataset->pointgroup_symbol),
                  std::begin(copy_dataset->pointgroup_symbol));


    }


} // SC