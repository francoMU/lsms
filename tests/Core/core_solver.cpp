/**
 *
 * RadialDFTLib
 *
 * core_solver.cpp
 *
 * Created by Franco Moitzi on 6/26/21.
 *
 * Copyright (c) 2021 University of Leoben. All rights reserved.
 *
 */

#include <cstdlib>

namespace TestCoreSolvers {

    bool basic() {

        return true;
    }

}


int main (int argc, char *argv[]) {

    auto test_result = TestCoreSolvers::basic();

    if (test_result) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }

}