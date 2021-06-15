
#include "PotentialShifter.hpp"

#include "Real.hpp"
#include "Main/SystemParameters.hpp"

void PotentialShifter::applyShifts(LocalTypeInfo &local) {


    for (int i = 0; i < local.num_local; i++) {
        for (int ir = 0; ir < local.atom[i].r_mesh.size(); ir++) {
//        Real deltaPotential=local.atom[i].vSpinShift*local.atom[i].r_mesh[ir]*
//                            (local.atom[i].rhotot(ir,0)-local.atom[i].rhotot(ir,1));
            Real deltaPotential = 0.5 * local.atom[i].vSpinShift * local.atom[i].r_mesh[ir];
            local.atom[i].vr(ir, 0) = vr0[i](ir, 0) - deltaPotential;
            local.atom[i].vr(ir, 1) = vr0[i](ir, 1) + deltaPotential;
        }
    }
}
