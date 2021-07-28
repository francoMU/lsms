

#include <cstdlib>

#include "Madelung/Multipole/multipole_madelung.h"

namespace TestMultipoleMadelung {

    bool initalization() {

        int num_local_atoms = 1;
        int num_atoms = 1;
        int gindex[] = {0};
        int lmax_rho = 1;
        int lmax_pot = 1;

        double bravais[] = {
                1.0, 0.0, 0.0,
                0.0, 1.0, 0.0,
                0.0, 0.0, 1.0,
        };

        double position[] = {
                0.0, 0.0, 0.0
        };

        initMadelung(num_local_atoms, num_atoms, gindex, lmax_rho, lmax_pot,
                     bravais, position, 1);

        return true;
    }

}


int main(int argc, char *argv[]) {

    auto test_result =
            TestMultipoleMadelung::initalization() && true;

    if (test_result) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }

}