
#include <iostream>

#include "spglib.h"

#include "Structure.hpp"

SC::Structure::Structure(const Eigen::Matrix<double, 3, 3> &lattice,
                         const Eigen::Matrix<double, Eigen::Dynamic, 3> &coordinates,
                         const Eigen::Vector<int, Eigen::Dynamic> &species) :

        lattice{lattice}, coordinates{coordinates}, species{species} {


    if (coordinates.rows() != species.rows()) {
        throw std::runtime_error("Length of species and coordinates is not the same!");
    }

    inv_lattice = lattice.inverse();

}

/**
 * @brief Perform the Lenstra-Lenstra-Lovász (LLL) lattice reduction algorithm.
 *
 * The LLL reduction algorithm is used to transform a given lattice matrix into a
 * nearly orthogonal basis with short, nearly integral vectors. It is particularly
 * useful for lattice problems and can simplify many computations involving lattices.
 *
 * The LLL algorithm is based on the Gram-Schmidt orthogonalization process with
 * additional optimizations to improve the orthogonality of the reduced basis.
 *
 * @param[out] lll The LLL-reduced lattice matrix will be stored in this parameter.
 * @param[out] lll_mapping The transformation matrix for LLL reduction will be stored
 *                         in this parameter. It represents the change of basis from
 *                         the input lattice to the reduced lattice.
 * @param[in] delta A parameter controlling the reduction quality. It should be a value
 *                  in the range (0, 1) where smaller values produce more orthogonal
 *                  and shorter vectors at the cost of increased computation time.
 *                  The default value is recommended for most cases.
 *
 * @note The input lattice matrix should be provided in the parameter `lattice` of the
 *       `SC::Structure` class. After calling this function, the output lattice matrix
 *       will be available in the `lll` parameter, and the transformation matrix will
 *       be stored in the `lll_mapping` parameter.
 *
 * @see SC::Structure::lattice, SC::Structure::lll_matrix, SC::Structure::lll_mapping
 */
void SC::Structure::calculateLLLreduction(
        Eigen::Matrix<double, 3, 3> &lll,
        Eigen::Matrix<double, 3, 3> &lll_mapping,
        double delta) {


    using namespace Eigen;

    int k;
    int q;

    Matrix<double, 3, 3> a;
    Matrix<double, 3, 3> b;
    Matrix<double, 3, 3> u;
    Vector<double, 3> m;
    double mu;
    Matrix<double, 3, 3> mapping;
    Vector<double, 10> uu;
    Vector<double, 3> v;
    Matrix<double, 3, 3> qq;
    Vector<double, 3> pp;

    // Basis vectors need to be in column representation
    a = lattice.transpose();

    b.setZero(); // Vectors after the Gram-Schmidt process
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

            a.col(k - 1).swap(a.col(k - 2));
            mapping.col(k - 1).swap(mapping.col(k - 2));

            for (int s = k - 1; s < k + 1; s++) {
                u(s - 1, seq(0, s - 2)) = a(indexing::all, s - 1).transpose() * b(indexing::all, seq(0, s - 2));
                u(s - 1, seq(0, s - 2)) = u(s - 1, seq(0, s - 2)).array() / m(seq(0, s - 2)).transpose().array();

                b(indexing::all, s - 1) = a(indexing::all, s - 1) -
                                          b(indexing::all, seq(0, s - 2)) * u(s - 1, seq(0, s - 2)).transpose();
                m(s - 1) = b(indexing::all, s - 1).transpose() * b(indexing::all, s - 1);

            }

            if (k > 2) {
                k -= 1;
            } else {


                Matrix<double, Dynamic, 1> p = (a(indexing::all, seq(k, 2)).transpose() *
                                                b(indexing::all, seq(k - 2, k - 1))).transpose();
                Matrix<double, Dynamic, Dynamic> q = m(seq(k - 2, k - 1)).asDiagonal();


                Matrix<double, Dynamic, Dynamic> result = q.colPivHouseholderQr().solve(p).transpose();

                u(seq(k, 2), seq(k - 2, k - 1)) = result;


            }


        }


    }


    lll = a.transpose();
    lll_mapping = mapping.transpose();

}

const Eigen::Matrix<double, 3, 3> &SC::Structure::lll_matrix() {

    if (map.count(DEFAULT_DELTA) == 0) {

        Eigen::Matrix<double, 3, 3> lll;
        Eigen::Matrix<double, 3, 3> lll_mapping;

        calculateLLLreduction(lll, lll_mapping, DEFAULT_DELTA);

        //std::tuple<Eigen::Matrix<double, 3, 3>, Eigen::Matrix<double, 3, 3>> entry;

        auto entry = std::make_tuple(lll, lll_mapping);

        map.insert({DEFAULT_DELTA, entry});

    }

    return std::get<0>(map[DEFAULT_DELTA]);
}

const Eigen::Matrix<double, 3, 3> &SC::Structure::lll_mapping() {
    if (map.count(DEFAULT_DELTA) == 0) {

        Eigen::Matrix<double, 3, 3> lll;
        Eigen::Matrix<double, 3, 3> lll_mapping;

        calculateLLLreduction(lll, lll_mapping, DEFAULT_DELTA);

        auto entry = std::make_tuple(lll, lll_mapping);

        map.insert({DEFAULT_DELTA, entry});

    }

    return std::get<1>(map[DEFAULT_DELTA]);
}

/**
 * @brief Calculate the shortest distance between two points in fractional coordinates.
 *
 * This function calculates the shortest distance between two points represented by their
 * fractional coordinates in the crystal lattice. The fractional coordinates are expressed
 * as a three-dimensional vector, representing the position of an atom within the unit cell
 * relative to the lattice vectors.
 *
 * @param[in] fcoords1 The fractional coordinates of the first point.
 * @param[in] fcoords2 The fractional coordinates of the second point.
 * @param[out] dist_vec A three-dimensional vector storing the shortest distance vector
 *                      between the two points. The vector points from the first point
 *                      (fcoords1) to the second point (fcoords2).
 * @param[out] dist The shortest distance between the two points in the crystal lattice.
 *
 * @note The input fractional coordinates `fcoords1` and `fcoords2` should be provided
 *       as three-dimensional vectors representing the fractional positions of the two
 *       points within the unit cell. After calling this function, the shortest distance
 *       vector between the points will be available in the `dist_vec` parameter, and the
 *       shortest distance itself will be stored in the `dist` parameter.
 *
 * @see SC::Structure::calculateLLLreduction, SC::Structure::get_distances_cart
 */
void SC::Structure::get_distances(const Eigen::Vector<double, 3> &fcoords1,
                                  const Eigen::Vector<double, 3> &fcoords2,
                                  Eigen::Vector<double, 3> &dist_vec,
                                  double &dist) {

    using namespace Eigen;

    Matrix<double, 3, 3> lll = lll_matrix();
    Matrix<double, 3, 3> lll_inverse = lll_mapping().inverse();

    Vector<double, 3> lll_fcoords1;
    Vector<double, 3> lll_fcoords2;

    Vector<double, 3> pre_image;
    Vector<double, 3> cart1;
    Vector<double, 3> cart2;
    Vector<double, 3> dv;

    Matrix<double, 3, 27> images;
    Matrix<double, 3, 27> cart_images;

    images << -1, -1, -1,
            -1, -1, 0,
            -1, -1, 1,
            -1, 0, -1,
            -1, 0, 0,
            -1, 0, 1,
            -1, 1, -1,
            -1, 1, 0,
            -1, 1, 1,

            0, -1, -1,
            0, -1, 0,
            0, -1, 1,
            0, 0, -1,
            0, 0, 0,
            0, 0, 1,
            0, 1, -1,
            0, 1, 0,
            0, 1, 1,

            1, -1, -1,
            1, -1, 0,
            1, -1, 1,
            1, 0, -1,
            1, 0, 0,
            1, 0, 1,
            1, 1, -1,
            1, 1, 0,
            1, 1, 1;


    lll_fcoords1 = fcoords1.transpose() * lll_inverse;
    lll_fcoords2 = fcoords2.transpose() * lll_inverse;

    lll_fcoords1.unaryExpr([](const int x) { return x % 2; });
    lll_fcoords2.unaryExpr([](const int x) { return x % 2; });


    cart1 = lll_fcoords1.transpose() * lll;
    cart2 = lll_fcoords2.transpose() * lll;

    cart_images = lll.transpose() * images;

    pre_image = cart2 - cart1;

    double best = 1.0e20;
    int bestk = 0;

    for (auto i = 0; i < 27; i++) {
        dv = pre_image + cart_images(indexing::all, i);
        double d = dv.squaredNorm();

        if (d < best) {
            best = d;
            bestk = i;
        }

    }


    dist_vec = pre_image + cart_images(indexing::all, bestk);
    dist = std::sqrt(best);

}


void SC::Structure::get_distances_cart(const Eigen::Vector<double, 3> &coords1,
                                       const Eigen::Vector<double, 3> &coords2,
                                       Eigen::Vector<double, 3> &dist_vec,
                                       double &dist) {


    get_distances(coords1.transpose() * inv_lattice, coords2.transpose() * inv_lattice, dist_vec, dist);

}

void
SC::Structure::get_fractional_coords(
        const Eigen::Matrix<double, Eigen::Dynamic, 3> &cart_coords,
        Eigen::Matrix<double, Eigen::Dynamic, 3> &frac_coords
) {
    frac_coords = cart_coords * inv_lattice;
}

void
SC::Structure::get_cartesian_coords(
        const Eigen::Matrix<double, Eigen::Dynamic, 3> &frac_coords,
        Eigen::Matrix<double, Eigen::Dynamic, 3> &cart_coords
) {
    cart_coords = frac_coords * lattice;
}


