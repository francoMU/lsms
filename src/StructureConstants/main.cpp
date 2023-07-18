

#include <iostream>

#include "Structure.hpp"
#include <Eigen/Dense>

int main(int argc, char *argv[]) {

    using namespace SC;
    using namespace Eigen;

    int k;
    int q;

    Matrix<double, 3, 3> lattice;
    Matrix<double, 3, 3> lll;
    Matrix<double, 3, 3> lll_mapping;
    Matrix<double, 3, 3> lll_inverse;

    lattice << 2.0, 0.0, 0.0,
            0.1, 1.8, 0.0,
            0.1, 0.2, 0.9;

    Matrix<double, 3, 3> a;
    Matrix<double, 3, 3> b;
    Matrix<double, 3, 3> u;
    Vector<double, 3> m;
    double mu;
    Matrix<double, 3, 3> mapping;
    Vector<double, 10> uu;
    Vector<double, 3> v;
    Vector<double, 3> v_m;
    Matrix<double, 3, 3> qq;
    Vector<double, 3> pp;

    double delta = 0.75;

    // Basis vectors need to be in column representation
    a = lattice.transpose();

    b.setZero(); // # Vectors after the Gram-Schmidt process
    u.setZero(); // Gram-Schmidt coefficients
    m.setZero(); //

    b.col(0) = a.col(0);
    m(0) = b.col(0).dot(b.col(0));

    for (auto i = 1; i <= 2; i++) {
        u(i, seq(0, i - 1)) = a(indexing::all, i).transpose() * b(indexing::all, seq(0, i - 1));
        u(i, seq(0, i - 1)) = u(i, seq(0, i - 1)).array() / m(seq(0, i - 1)).transpose().array();
        b(indexing::all, i) = a(indexing::all, i) -
                              b(indexing::all, seq(0, i - 1)) * u(i, seq(0, i - 1)).transpose();
        m(i) = b(indexing::all, i).transpose() * b(indexing::all, i);
    }



    k = 2;

    mapping.setIdentity();

    while (k <= 3) {


        for (int i = k - 1; i > 0; i--) {

            mu = u(k - 1, i - 1);

            if (std::fabs(mu) > 0.5) {

                q = std::rint(u(k - 1, i - 1));
                a(indexing::all, k - 1) = a(indexing::all, k - 1) - q * a(indexing::all, i - 1);
                mapping(indexing::all, k - 1) = mapping(indexing::all, k - 1) - q * mapping(indexing::all, i - 1);

                uu.setZero();

                if (i >= 2) {
                    uu(seq(0, (i - 2))) = u(i - 1, seq(0, (i - 2)));
                }

                uu(i - 1) = 1.0;

                u(k - 1, seq(0, i - 1)) = u(k - 1, seq(0, i - 1)) - q * uu(seq(0, i - 1));

            }
        }


        double norm = b(indexing::all, k - 1).squaredNorm();
        std::cout << norm << std::endl;

//        if ()


    }


    //       do i = 1, 2
    //       u(i, 0:i - 1) = matmul(a(:, i), b(:, 0:i - 1)) / m(0:i - 1)
    //       b(:, i) = a(:, i) - matmul(b(:, 0:i - 1), u(i, 0:i - 1))
    //       m(i) = dot_product(b(:, i), b(:, i))
    //    end do

    // b[:, 0] = a[:, 0]


}