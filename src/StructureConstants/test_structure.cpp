//
// Created by F.Moitzi on 15.05.2022.
//

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

#include <Eigen/Dense>

#include "Structure.hpp"

bool CheckMatrixEquality(const Eigen::Matrix<double, 3, 3> &lhs, const Eigen::Matrix<double, 3, 3> &rhs) {
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
