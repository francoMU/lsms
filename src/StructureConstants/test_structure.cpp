//
// Created by F.Moitzi on 15.05.2022.
//

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

#include <Eigen/Dense>

#include "Structure.hpp"
#include "SpglibWrapper.h"

bool CheckMatrixEquality(const Eigen::Matrix<double, Eigen::Dynamic, 3> &lhs,
                         const Eigen::Matrix<double, Eigen::Dynamic, 3> &rhs) {
    return lhs.isApprox(rhs, 1e-4);
}

TEST(StructureTests, LLLMatrix1) {

    using namespace SC;
    using namespace Eigen;

    Matrix<double, 3, 3> lattice;
    Matrix<double, 3, 3> lll_matrix_ref;

    lattice
            << 8.0, 0.0, 0.0,
            0.0, 3.0, 0.0,
            -0.1, 10.0, 1.0;

    lll_matrix_ref
            << -0.1, 1, 1,
            0.1, 2, -1,
            8, 0, 0;


    Matrix<double, Dynamic, 3> coordinates(2, 3);
    Vector<int, Dynamic> species(2);

    coordinates << 0.5, 0.5, 0.5, 1.5, 0.5, 0.5;
    species << 1, 0;

    SC::Structure structure(lattice, coordinates, species);

    auto lll_matrix = structure.lll_matrix();
    ASSERT_PRED2(CheckMatrixEquality, lll_matrix_ref, lll_matrix);

}

TEST(StructureTests, LLLMatrix2) {


    using namespace SC;
    using namespace Eigen;

    Matrix<double, 3, 3> lattice;
    Matrix<double, 3, 3> lll_matrix_ref;

    lattice
            << 1.0, 2.0, 3.0,
            8.0, 7.0, -6.0,
            -4.0, 3.0, 2.0;

    lll_matrix_ref
            << 1.0, 2.0, 3.0,
            -5.0, 1.0, -1.0,
            3.0, 8.0, -7.0;


    Matrix<double, Dynamic, 3> coordinates(2, 3);
    Vector<int, Dynamic> species(2);

    coordinates << 0.5, 0.5, 0.5, 1.5, 0.5, 0.5;
    species << 1, 0;

    SC::Structure structure(lattice, coordinates, species);

    auto lll_matrix = structure.lll_matrix();
    ASSERT_PRED2(CheckMatrixEquality, lll_matrix_ref, lll_matrix);

    Matrix<double, Dynamic, 3> frac_coords(2, 3);
    Matrix<double, Dynamic, 3> cart_coords(2, 3);


    frac_coords << 0.5, 0.5, 0.5, 0, 0, 0.2;

    Matrix<double, Dynamic, 3> cart_coords_ref(2, 3);

    cart_coords_ref << 2.5, 6., -0.5, -0.8, 0.6, 0.4;


    structure.get_cartesian_coords(frac_coords, cart_coords);


    ASSERT_PRED2(CheckMatrixEquality, cart_coords_ref, cart_coords);

}


TEST(StructureTests, Distance1) {


    using namespace SC;
    using namespace Eigen;

    Matrix<double, 3, 3> lattice;
    Matrix<double, 3, 3> lll_matrix_ref;

    lattice
            << 2.0, 0.0, 0.0,
            0.0, 2.0, 0.0,
            0.0, 0.0, 2.0;

    Matrix<double, Dynamic, 3> coordinates(2, 3);
    Vector<int, Dynamic> species(2);

    coordinates << 0.5, 0.5, 0.5, 1.5, 0.5, 0.5;
    species << 1, 0;

    SC::Structure structure(lattice, coordinates, species);

    Vector<double, Dynamic> fcoords1(3);
    Vector<double, Dynamic> fcoords2(3);
    Vector<double, 3> dist_vector;
    double dist;

    fcoords1 << 0, 0, 0;
    fcoords2 << 0.9, 0.9, 0.9;

    structure.get_distances(fcoords1, fcoords2, dist_vector, dist);

    ASSERT_NEAR(0.34641016, dist, 1.0e-8);


}

TEST(StructureTests, Symmetry1) {


    using namespace SC;
    using namespace Eigen;

    Matrix<double, 3, 3> lattice;
    Matrix<double, 3, 3> lll_matrix_ref;

    lattice
            << 4.0, 0.0, 0.0,
            0.0, 4.0, 0.0,
            0.0, 0.0, 4.0;

    Matrix<double, Dynamic, 3> coordinates(2, 3);
    Vector<int, Dynamic> species(2);

    coordinates << 0.5, 0.5, 0.5, 0.0, 0.0, 0.0;
    species << 1, 0;

    SpglibWrapper wrapper(lattice, coordinates, species);

    ASSERT_EQ(wrapper.international_symbol(), "Pm-3m");

}

TEST(StructureTests, MappingOfStorage) {

    using namespace Eigen;

    int (*array_storage)[2][2] = new int[2][2][2]();

    int *array;

    array = new int[8]{1, 2, 3, 4, 5, 6, 7, 8};

    std::cout << "Row-major using stride:\n" <<
              Map<Matrix<int, 2, 4>, Unaligned, Stride<1, 4> >(array) << std::endl;

    int o = 0;

    for (auto i = 0; i < 2; i++) {
        for (auto j = 0; j < 2; j++) {
            for (auto k = 0; k < 2; k++) {
                o++;
                array_storage[i][j][k] = o;
            }
        }
    }

    std::cout << "Row-major using stride:\n" <<
              Map<Matrix<int, 2, 4>, Unaligned, Stride<1, 4> >(reinterpret_cast<int *>(array_storage)) << std::endl;
}


