

#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <complex>
#include <cmath>


#include "multipole_madelung.h"

namespace TestImprovedMadelung {

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
                0.0, 0.5, 0.3,
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

        /*
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
        std::cout << "M00_11: " << " " << M00_1 << std::endl;

        id = 1;
        jd = 0;
        k = 0;
        auto M00_12 = DL_matrix_ptr[jd + num_atoms * k + num_atoms * 9 * id];
        std::cout << "M00_12: " << " " << M00_12 << std::endl;

        id = 1;
        jd = 3;
        k = 3;
        auto M11_23 = DL_matrix_ptr[jd + num_atoms * k + num_atoms * 9 * id];
        std::cout << "M11_23: " << " " << M11_23 << std::endl;

        // Prefactor L,L_prime = {1,0},{0,0}
        auto prefactor_10 = DL_factor_ptr[9 * 2];
        std::cout << "P 1000: " << prefactor_10 << std::endl;

        // Prefactor L,L_prime = {1,1},{0,0}
        auto prefactor_11 = DL_factor_ptr[9 * 3];
        std::cout << "P 1100: " << prefactor_11 << std::endl;

        endMadelung();

        return true;
    }

}


int main(int argc, char *argv[]) {

    auto test_result =
            TestImprovedMadelung::MadlungMatrixComponents();

    if (test_result) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }

}