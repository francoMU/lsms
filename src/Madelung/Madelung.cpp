#include "Madelung.hpp"
#include "Main/SystemParameters.hpp"

void calculateMadelungMatrices(LSMSSystemParameters &lsms, CrystalParameters &crystal, LocalTypeInfo &local) {
    std::vector<Real> atom_position_1(crystal.num_atoms);
    std::vector<Real> atom_position_2(crystal.num_atoms);
    std::vector<Real> atom_position_3(crystal.num_atoms);
    for (int i = 0; i < crystal.num_atoms; i++) {
        atom_position_1[i] = crystal.position(0, i);
        atom_position_2[i] = crystal.position(1, i);
        atom_position_3[i] = crystal.position(2, i);
    }

    printf("Madelung matrix is calculated\n");

    // int lmax=1;
    // int ndlmadv=((lmax+1)*(lmax+2))/2;
    int num_atoms = crystal.num_atoms;
    for (int i = 0; i < local.num_local; i++) {
        local.atom[i].madelungMatrix.resize(num_atoms);
        // local.atom[i].madelungMatrixJ.resize(ndlmadv,num_atoms);
    }

    printf("Number of atoms: %d \n", num_atoms);

// cal_madelung_matrix appears not to be thread save.
    for (int i = 0; i < local.num_local; i++) {
        int mynod = local.global_id[i];
        cal_madelung_matrix_(&mynod, &num_atoms,
                             &crystal.bravais(0, 0),
                             atom_position_1.data(),
                             atom_position_2.data(),
                             atom_position_3.data(),
                             &local.atom[i].madelungMatrix[0],
                             &lsms.global.iprint, lsms.global.istop, 32);

        for (int j = 0; j < crystal.num_types; j++) {
            printf("Madelung matrix elements: %f\n", local.atom[i].madelungMatrix[j]);
        }
    }


}

