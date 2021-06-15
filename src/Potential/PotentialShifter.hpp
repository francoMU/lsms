#ifndef LSMS_SHIFT_POTENTIALS
#define LSMS_SHIFT_POTENTIALS

#include <vector>

#include "Real.hpp"
#include "Matrix.hpp"

#include "Main/SystemParameters.hpp"

class PotentialShifter {

public:

    bool vSpinShiftFlag{false};
    double minShift{0.0};
    double maxShift{0.0};

    void resetPotentials(LocalTypeInfo &local) {
        for (int i = 0; i < local.num_local; i++)
            vr0[i] = local.atom[i].vr;
    }

    void resize(int n) { vr0.resize(n); }

    void applyShifts(LocalTypeInfo &local);

    void restorePotentials(LocalTypeInfo &local) {
        for (int i = 0; i < local.num_local; i++)
            local.atom[i].vr = vr0[i];
    }

    std::vector<Matrix<Real> > vr0{}; // the unshifted potential

};

#endif
