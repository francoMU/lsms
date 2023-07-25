

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


        /*
         * Size reduction
         */
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

                u(k - 1, seq(0, i - 1)).array() = u(k - 1, seq(0, i - 1)).array()
                                                  - q * uu(seq(0, i - 1)).transpose().array();

            }
        }


        double norm1 = b(indexing::all, k - 1).squaredNorm();
        double norm2 = (delta - std::pow(std::fabs(u(k - 1, k - 2)), 2)) * b(indexing::all, k - 2).squaredNorm();

        if (norm1 >= norm2) {
            k += 1;
        } else {

            //v = a(indexing::all, k - 1);
            //a(indexing::all, k - 1) = a(indexing::all, k - 2);
            //a(indexing::all, k - 2) = v;

            a.col(k - 1).swap(a.col(k - 2));
            mapping.col(k - 1).swap(mapping.col(k - 2));

            // H.col(i).swap(H.col(j))

            //v_m = mapping(indexing::all, k - 1);
            //mapping(indexing::all, k - 1) = mapping(indexing::all, k - 2);
            //mapping(indexing::all, k - 2) = v_m;


            for (int s = k - 1; s < k + 1; s++) {
                u(s - 1, seq(0, s - 2)) = a(indexing::all, s - 1).transpose() * b(indexing::all, seq(0, s - 2));
                u(s - 1, seq(0, s - 2)) = u(s - 1, seq(0, s - 2)).array() / m(seq(0, s - 2)).transpose().array();

                b(indexing::all, s - 1) = a(indexing::all, s - 1) -
                                          b(indexing::all, seq(0, s - 2)) * u(s - 1, seq(0, s - 2)).transpose();
                m(s - 1) = b(indexing::all, s - 1).transpose() * b(indexing::all, s - 1);

            }

            if (k > 2) {
                k -= 1;
            }


        }


    }

    std::cout << a.transpose() << std::endl;
    std::cout << mapping.transpose() << std::endl;

    /**
     * 0.1 0.2 0.9
  2   0   0
0.1 1.8   0
0 0 1
1 0 0
0 1 0

     */

}