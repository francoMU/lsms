
#include <algorithm>
#include <numeric>

#include "calcHigherOrderMadelung.hpp"

#include "Main/SystemParameters.hpp"
#include "multipole_madelung.h"

void lsms::calculateHigherOrderMadelung(LSMSSystemParameters &lsms,
                                        CrystalParameters &crystal,
                                        LocalTypeInfo &local,
                                        lsms::HigherOrderMadelung &madelung) {

    constexpr auto i_print = 0;
    constexpr auto kmax = 4;

    // Number of atoms
    auto num_atoms = lsms.num_atoms;
    auto num_local_atoms = local.num_local;

    // Global index vector has to be converted to id = id + 1
    std::vector<int> global_index = local.global_id;
    std::transform(global_index.begin(),
                   global_index.end(),
                   global_index.begin(),
                   [&](auto &value) {
                       return value + 1;
                   });

    initMadelung(
            num_atoms,
            num_local_atoms,
            global_index.data(),
            0,
            1,
            &crystal.bravais[0],
            &crystal.position[0],
            i_print);

    madelung.pre_factor_10_00 = DL_factor_ptr[kmax * 1];
    madelung.pre_factor_11_00 = DL_factor_ptr[kmax * 2];

    madelung.lattice_scale = *alat_ptr;

    for (int i = 0; i < num_atoms; i++) {
        for (int j = 0; j < num_local_atoms; j++) {
            madelung.G_11(i, j) = DL_matrix_ptr[i + num_atoms * 1 + num_atoms * kmax * j];
            madelung.G_10(i, j) = DL_matrix_ptr[i + num_atoms * 2 + num_atoms * kmax * j];
            madelung.G_1m1(i, j) = DL_matrix_ptr[i + num_atoms * 3 + num_atoms * kmax * j];
        }
    }

    endMadelung();

}
