//
// Created by F.Moitzi on 30.07.2023.
//

#include "SpglibWrapper.h"

namespace SC {


    SpglibWrapper::SpglibWrapper(const Eigen::Matrix<double, 3, 3> &lattice,
                                 const Eigen::Matrix<double, Eigen::Dynamic, 3> &frac_coordinates,
                                 const Eigen::Vector<int, Eigen::Dynamic> &species,
                                 double symprec) : symprec{symprec} {


        int num_atom = static_cast<int>(frac_coordinates.rows());

        this->dataset = spg_get_dataset(reinterpret_cast<const double (*)[3]>(lattice.data()),
                                        reinterpret_cast<const double (*)[3]>(frac_coordinates.data()),
                                        reinterpret_cast<const int *>(species.data()),
                                        num_atom, symprec);


    }

    SpglibWrapper::SpglibWrapper(const SpglibWrapper &other)
            : symprec(other.symprec), dataset(new SpglibDataset()) {
        copy(other.dataset, dataset);
    }


    SpglibWrapper &SpglibWrapper::operator=(const SpglibWrapper &rhs) {

        if (this != &rhs) // not a self-assignment
        {
            resize(dataset, rhs.dataset->n_operations, rhs.dataset->n_atoms, rhs.dataset->n_std_atoms);
            copy(rhs.dataset, this->dataset);
        }

        return *this;

    }

    SpglibWrapper::SpglibWrapper(SpglibWrapper &&other) noexcept
            : dataset(nullptr) {
        *this = std::move(other);
    }


    SpglibWrapper &SpglibWrapper::operator=(SpglibWrapper &&other) noexcept {

        if (this != &other) {

            if (dataset != nullptr) {
                spg_free_dataset(dataset);
            }

            dataset = other.dataset;
            symprec = other.symprec;

            other.dataset = nullptr;
        }
        return *this;

    }


    void SpglibWrapper::copy(const SpglibDataset *rhs_dataset, SpglibDataset *copy_dataset) {

        /* Space group type */
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

        /* Transformation matrix and origin shift */
        std::copy_n(&rhs_dataset->transformation_matrix[0][0], 9, &copy_dataset->transformation_matrix[0][0]);

        std::copy(std::begin(rhs_dataset->origin_shift),
                  std::end(rhs_dataset->origin_shift),
                  std::begin(copy_dataset->origin_shift));

        /* Symmetry operations */
        auto n_ops = rhs_dataset->n_operations;
        copy_dataset->n_operations = n_ops;

        std::copy_n(&rhs_dataset->rotations[0][0][0], n_ops * 9, &copy_dataset->rotations[0][0][0]);
        std::copy_n(&rhs_dataset->translations[0][0], n_ops * 3, &copy_dataset->translations[0][0]);

        /* Wyckoff positions and symmetrically equivalent atoms */
        auto n_atoms = rhs_dataset->n_atoms;
        copy_dataset->n_atoms = n_atoms;

        std::copy_n(rhs_dataset->wyckoffs,
                    n_atoms,
                    copy_dataset->wyckoffs);

        std::copy_n(rhs_dataset->equivalent_atoms,
                    n_atoms,
                    copy_dataset->equivalent_atoms);

        std::copy_n(rhs_dataset->crystallographic_orbits,
                    n_atoms,
                    copy_dataset->crystallographic_orbits);

        /* Intermediate data in symmetry search */

        std::copy_n(&rhs_dataset->primitive_lattice[0][0], 9, &copy_dataset->primitive_lattice[0][0]);
        std::copy_n(rhs_dataset->mapping_to_primitive,
                    n_atoms,
                    copy_dataset->mapping_to_primitive);

        /* Standardized lattice operation */
        auto n_std_atoms = rhs_dataset->n_std_atoms;

        copy_dataset->n_std_atoms = n_std_atoms;

        std::copy_n(rhs_dataset->std_types, n_std_atoms, copy_dataset->std_types);
        std::copy_n(&rhs_dataset->std_lattice[0][0], 9, &copy_dataset->std_lattice[0][0]);
        std::copy_n(&rhs_dataset->std_positions[0][0], n_std_atoms * 3, &copy_dataset->std_positions[0][0]);
        std::copy_n(&rhs_dataset->std_rotation_matrix[0][0], 9, &copy_dataset->std_rotation_matrix[0][0]);
        std::copy_n(rhs_dataset->std_mapping_to_primitive, n_std_atoms, copy_dataset->std_mapping_to_primitive);

        /* Crystallographic point group */
        std::copy(std::begin(rhs_dataset->pointgroup_symbol),
                  std::end(rhs_dataset->pointgroup_symbol),
                  std::begin(copy_dataset->pointgroup_symbol));


    }

    SpglibWrapper::~SpglibWrapper() {
        spg_free_dataset(dataset);
    }

    void SpglibWrapper::resize(SpglibDataset *rhs_dataset, int n_ops, int n_atoms, int n_std_atoms) {

        if (n_ops /= rhs_dataset->n_operations) {

            if (rhs_dataset->n_operations != 0) {
                delete[] rhs_dataset->rotations;
                delete[] rhs_dataset->translations;
            }

            rhs_dataset->rotations = new int[n_ops][3][3];
            rhs_dataset->translations = new double[n_ops][3];
        }

        if (n_atoms /= rhs_dataset->n_atoms) {

            if (rhs_dataset->n_atoms != 0) {
                delete[] rhs_dataset->crystallographic_orbits;
                delete[] rhs_dataset->wyckoffs;
                delete[] rhs_dataset->equivalent_atoms;
                delete[] rhs_dataset->mapping_to_primitive;
            }

            rhs_dataset->crystallographic_orbits = new int[n_atoms];
            rhs_dataset->wyckoffs = new int[n_atoms];
            rhs_dataset->equivalent_atoms = new int[n_atoms];
            rhs_dataset->mapping_to_primitive = new int[n_atoms];
        }

        if (n_std_atoms != rhs_dataset->n_std_atoms) {

            if (rhs_dataset->n_std_atoms != 0) {
                delete[] rhs_dataset->std_types;
                delete[] rhs_dataset->std_mapping_to_primitive;
            }

            rhs_dataset->std_types = new int[n_std_atoms];
            rhs_dataset->std_mapping_to_primitive = new int[n_std_atoms];
        }

    }

    std::string SpglibWrapper::international_symbol() {
        std::string value{dataset->international_symbol};
        ba::trim(value);
        return value;
    }


} // SC