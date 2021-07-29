

#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <complex>
#include <cmath>


#include "multipole_madelung.h"

namespace TestMultipoleMadelung {

    bool initialization() {

        int num_local_atoms = 2;
        int num_atoms = 2;
        int gindex[] = {1, 2};
        int lmax_rho = 3;
        int lmax_pot = 3;

        double bravais[] = {
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
        };

        double position[] = {
                0.0, 0.0, 0.0,
                0.5, 0.5, 0.5
        };

        initMadelung(num_local_atoms,
                     num_atoms,
                     gindex,
                     lmax_rho,
                     lmax_pot,
                     bravais,
                     position,
                     1);

        endMadelung();

        return true;
    }

    bool madelungMatrix() {

        int num_local_atoms = 2;
        int num_atoms = 2;
        int gindex[] = {1, 2};
        int lmax_rho = 2;
        int lmax_pot = 2;

        double bravais[] = {
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
        };

        double position[] = {
                0.0, 0.0, 0.0,
                0.5, 0.5, 0.5
        };

        initMadelung(num_local_atoms,
                     num_atoms,
                     gindex,
                     lmax_rho,
                     lmax_pot,
                     bravais,
                     position,
                     1);


        double *madelung_matrix = nullptr;
        int madelung_matrix_size;

        getMadelungMatrix(1,
                          &madelung_matrix,
                          &madelung_matrix_size);

        assert (madelung_matrix_size == 2);

        auto ref_value = -2.837;
        auto comp_value = std::round(madelung_matrix[0] * 1000.0) / 1000.0;

        assert (ref_value == comp_value);

        endMadelung();

        return true;
    }

    bool DLFactor() {

        int num_local_atoms = 2;
        int num_atoms = 2;
        int gindex[] = {1, 2};
        int lmax_rho = 2;
        int lmax_pot = 2;

        double bravais[] = {
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
        };

        double position[] = {
                0.0, 0.0, 0.0,
                0.5, 0.5, 0.5
        };

        initMadelung(num_local_atoms,
                     num_atoms,
                     gindex,
                     lmax_rho,
                     lmax_pot,
                     bravais,
                     position,
                     1);


        double *matrix;
        int size;

        getDLFactor(1,
                    &matrix,
                    &size);

        assert (size == 25);

        endMadelung();

        return true;
    }

    bool MadlungMatrixComponents() {

        int num_local_atoms = 2;
        int num_atoms = 4;
        int gindex[] = {1, 2, 3, 4};
        int lmax_rho = 1;
        int lmax_pot = 1;

        double bravais[] = {
                3.0, 0.0, 0.0,
                0.0, 3.0, 0.0,
                0.0, 0.0, 3.0,
        };

        double position[] = {
                0.0, 0.0, 0.0,
                0.5, 0.5, 0.0,
                0.0, 0.5, 0.5,
                0.5, 0.0, 0.5
        };

        initMadelung(num_local_atoms,
                     num_atoms,
                     gindex,
                     lmax_rho,
                     lmax_pot,
                     bravais,
                     position,
                     1);


        int atom_size;
        int k_size;
        double factor;


        //double *matrix_2;
        //testMatrix(&matrix_2);
        //std::cout << matrix_2[0] << std::endl;
        //std::cout << matrix_2[1] << std::endl;
        //std::cout << matrix_2[2] << std::endl;
        //std::cout << matrix_2[(2-1)+(3-1)*2+(4-1)*3*2] << std::endl;

        /*
         * Compare indexing
         */

        int id = 0;
        int jd = 0;
        int k = 0;
        std::complex<double> *matrix;

        getDLMatrix(id,
                    &matrix,
                    &atom_size,
                    &k_size,
                    &factor);

        std::complex<double> M00_1 = DL_matrix_ptr[0];
        std::cout << "M00_11" << " " << M00_1 << std::endl;

        id = 1;
        jd = 0;
        k = 0;
        auto M00_12 = DL_matrix_ptr[jd + num_atoms * k + num_atoms * 9 * id];
        std::cout << "M00_12" << " " << M00_12 << std::endl;

        id = 1;
        jd = 3;
        k = 3;
        auto M11_23 = DL_matrix_ptr[jd + num_atoms * k + num_atoms * 9 * id];
        std::cout << "M11_23" << " " << M11_23 << std::endl;

        /*
         * Get force matrix
         */


        std::complex<double> qlm_00[] = {15, 13, 14, 12};
        double z[] = {12, 13, 14, 15};


        auto alat = *alat_ptr;

        id = 0;

        auto force_mm_x = 0.0;
        auto force_mm_y = 0.0;
        auto force_mm_z = 0.0;

        for (int ig = 0; ig < num_atoms; ig++) {

            // Prefactor L,L_prime = {1,0},{0,0}
            auto prefactor_10 = DL_factor_ptr[9 * 2];
            std::cout << prefactor_10 << std::endl;


            // Prefactor L,L_prime = {1,1},{0,0}
            auto prefactor_11 = DL_factor_ptr[9 * 3];
            std::cout << prefactor_11 << std::endl;

            // Reduced madelung matrix {l+l_prime},{m-m_prime}
            // L = {1,0} (2)
            auto M_10 = DL_matrix_ptr[ig + num_atoms * 2 + num_atoms * 9 * id];
            std::cout << M_10 << std::endl;

            // L = {1,1} (3)
            auto M_11 = DL_matrix_ptr[ig + num_atoms * 3 + num_atoms * 9 * id];
            std::cout << M_11 << std::endl;

            force_mm_x += -std::sqrt(3.0 / (2.0 * M_PI))
                          * std::real(M_11 * prefactor_11 * (qlm_00[ig] - z[ig]) / alat);

            force_mm_y += std::sqrt(3.0 / (2.0 * M_PI))
                          * std::imag(M_11 * prefactor_11 * (qlm_00[ig] - z[ig]) / alat);

            force_mm_z += std::sqrt(3.0 / (4.0 * M_PI))
                          * std::real(M_10 * prefactor_10 * (qlm_00[ig] - z[ig]) / alat);

        }

        std::cout << force_mm_x << std::endl;
        std::cout << force_mm_y << std::endl;
        std::cout << force_mm_z << std::endl;

        /*
         * x direction
         */




        endMadelung();

        return true;
    }

    bool ForceComponents() {

        int num_local_atoms = 2;
        int num_atoms = 2;
        int gindex[] = {1, 2};
        int lmax_rho = 1;
        int lmax_pot = 1;

        double bravais[] = {
                3.0, 0.0, 0.0,
                0.0, 3.0, 0.0,
                0.0, 0.0, 3.0,
        };

        double position[] = {
                0.0, 0.0, 0.0,
                0.55, 0.55, 0.55
        };

        initMadelung(num_local_atoms,
                     num_atoms,
                     gindex,
                     lmax_rho,
                     lmax_pot,
                     bravais,
                     position,
                     1);


        std::complex<double> *matrix;

        double qlm_00[] = {3.0, 13.5};
        double qis[] = {0.25, 0.25};
        double z[] = {4.0, 13.0};

        auto alat = *alat_ptr;

        constexpr auto kmax = 4;
        constexpr auto jmax = 3;

        /*
         * Prefactor is real quantity
         */

        // Prefactor L,L_prime = {1,0},{0,0}
        auto prefactor_10 = DL_factor_ptr[kmax * 1];
        std::cout << prefactor_10 << " " << 9.4031597257959509E-002 << std::endl;

        // Prefactor L,L_prime = {1,1},{0,0}
        auto prefactor_11 = DL_factor_ptr[kmax * 2];
        std::cout << prefactor_11 << " " << -9.4031597257959509E-002 << std::endl;


        for (int id = 0; id < num_local_atoms; id++) {

            std::cout << "Main index: " << id << std::endl;

            auto force_mm_x = 0.0;
            auto force_mm_y = 0.0;
            auto force_mm_z = 0.0;

            for (int ig = 0; ig < num_atoms; ig++) {

                std::cout << "Subindex: " << ig << std::endl;

                /*
                 * Reduced Madelung matrix or Green's function is complex number
                 */

                // Reduced madelung matrix {l+l_prime},{m-m_prime}
                // L = {1,0} (2)
                auto G_10 = DL_matrix_ptr[ig + num_atoms * 2 + num_atoms * kmax * id];
                std::cout << G_10 << std::endl;

                // L = {1,1} (1)
                // Note: This seems to be then defined
                auto G_11 = DL_matrix_ptr[ig + num_atoms * 1 + num_atoms * kmax * id];
                std::cout << G_11 << std::endl;

                // L = {1,-1} (3)
                auto G_1m1 = DL_matrix_ptr[ig + num_atoms * 3 + num_atoms * kmax * id];
                std::cout << G_1m1 << std::endl;

                force_mm_x += -std::sqrt(3.0 / (2.0 * M_PI))
                              * std::real(G_11 * prefactor_11) * (qlm_00[ig] - z[ig] - qis[ig]) / alat
                              * (z[id] - qlm_00[id] + qis[id]);

                force_mm_y += std::sqrt(3.0 / (2.0 * M_PI))
                              * std::imag(G_11 * prefactor_11) * (qlm_00[ig] - z[ig] - qis[ig]) / alat
                              * (z[id] - qlm_00[id] + qis[id]);

                force_mm_z += std::sqrt(3.0 / (4.0 * M_PI))
                              * std::real(G_10 * prefactor_10) * (qlm_00[ig] - z[ig] - qis[ig]) / alat
                              * (z[id] - qlm_00[id] + qis[id]);

            }

            std::cout << "Force x: " << force_mm_x << std::endl;
            std::cout << "Force y: " << force_mm_y << std::endl;
            std::cout << "Force z: " << force_mm_z << std::endl;

        }

        endMadelung();

        return true;
    }

}


int main(int argc, char *argv[]) {

    auto test_result =
            TestMultipoleMadelung::ForceComponents();

    /*auto test_result =
            TestMultipoleMadelung::initialization() &&
            TestMultipoleMadelung::madelungMatrix() &&
            TestMultipoleMadelung::DLMatrix() &&
            TestMultipoleMadelung::DLFactor() &&
            TestMultipoleMadelung::MadlungMatrixComponents() &&
            true;
    */

    if (test_result) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }

}