


#include "higherOrderMadelung.hpp"

lsms::HigherOrderMadelung::HigherOrderMadelung(int num_local_atoms, int num_atoms) {

    G_11 = Matrix<Complex>(num_atoms, num_local_atoms);

    G_10 = Matrix<Complex>(num_atoms, num_local_atoms);

    G_1m1 = Matrix<Complex>(num_atoms, num_local_atoms);

}
